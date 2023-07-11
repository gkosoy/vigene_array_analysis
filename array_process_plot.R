library(tidyverse)
library(ggpmisc)
library(readxl)

#raw data input
###begin####
#order of experiments 
order_5_16 <- readxl::read_xlsx("5-30-23/5-16-23_chip_order.xlsx") 
#layout on all experimental chips
chip_map_5_16 <- readxl::read_xlsx("5-30-23/5-16-23_vir_chipmap.xlsx")

df_ex_15 <- read.delim(file="5-30-23/5_16_23_vic_ser_15ms.spots_2.txt", sep = "\t", 
                                header = TRUE, skip = 15) %>% 
  mutate(exposure = 15) %>%
  right_join(chip_map_5_16) %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  filter(!type == 11)
df_ex_20 <- read.delim(file="5-30-23/5_16_23_vic_ser_20ms.spots_2.txt", sep = "\t", 
                                header = TRUE, skip = 15) %>% 
  mutate(exposure = 20) %>%
  right_join(chip_map_5_16) %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  filter(!type == 11)
df_ex_40 <- read.delim(file="5-30-23/5_16_23_vic_ser_40ms.spots_2.txt", sep = "\t", 
                                header = TRUE, skip = 15) %>% 
  mutate(exposure = 40) %>%
  right_join(chip_map_5_16) %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  filter(!type == 11)

full_ser <- rbind(df_5_16_vicser_15, df_5_16_vicser_20, 
                      df_5_16_vicser_40)
ser_vir_5_16 <- order_5_16  %>% right_join(full_ser) %>%
  filter(!exposure %in% c(75, 200))
ser_vir_5_16_cut <- ser_vir_5_16 %>%
  select(thickness_median_total, protein, uL, concentration,
         Main_Column)

####end####

#average spots across arrays, sub-neg ctrl chip, sub-neg ctrl on chip, propagate error
vir_sum_func <- function(df) {
  df_info <- df %>% 
    group_by(serum, dilution, Main_Column, protein, concentration) %>%
    summarize(mean_thickness = mean(thickness_median_total),
              std_thickness = sd(thickness_median_total),
              count = n(),
              sem_thickness = std_thickness/(sqrt(count)),
              .groups = 'drop') %>%
    group_by(protein, concentration) %>%
    mutate(ctrl_change = mean_thickness - mean(mean_thickness[serum == "con" & Main_Column == 3]),
           ctrl_sd = sqrt((std_thickness^2)+(std_thickness[serum == "con" & Main_Column == 3]^2)))%>% 
    group_by(serum, Main_Column, dilution) %>%
    mutate(FITC_change = ctrl_change - (ctrl_change %>% 
                                          keep(protein == "FITC" & concentration == 350)),
           FITC_sd = sqrt((ctrl_sd^2)+(ctrl_sd %>% 
                                         keep(protein == "FITC" & concentration == 350))^2)) 
  return(df_info)
}
vir_sum_df <- vir_sum_func(ser_vir_5_16)
vir_sum_df$dilution <- as.numeric(vir_sum_df$dilution)

#plot binding vs serum dilution for probe/protein of interest with binding model
plot <- ggplot(vir_sum_df %>% filter(protein == "anti-IgG"), 
               aes(x=dilution, y=FITC_change, 
                   color=as.factor(concentration))) + 
  geom_point(size = 8)+
  facet_wrap(~as.factor(concentration))+
  geom_pointrange(aes(ymin=(FITC_change)-FITC_sd, 
                      ymax=(FITC_change)+FITC_sd)) +
  theme_bw()+
  xlab("dilution")+ylab("thickness")+labs(col="concentration")+
  theme(text=element_text(size = 35),
        legend.title=element_text(size=5))+
  geom_smooth(method="nls", formula=y~(Bmax*x/(Kd + x)), 
              method.args = list(start = list(Bmax=25,Kd=0.005)), se = FALSE, show.legend = TRUE)+
  stat_fit_tidy(method = "nls", 
                method.args = list(formula = my_eq,
                                   start = list(Bmax=25,KD=0.005)),
                label.x = "right",
                label.y = "bottom",
                aes(label = paste("B[max]~`=`~", signif(after_stat(Bmax_estimate), digits = 2),
                                  #  "%+-%", signif(after_stat(Bmax_se), digits = 2),
                                  "~~~~K[D]~`=`~", signif(after_stat(KD_estimate), digits = 2),
                                  #"%+-%", signif(after_stat(KD_se), digits = 2),
                                  sep = "")),
                parse = TRUE, size = 10)
tiff("5-30-23/7_10_vicser_IgG.tiff", units="in", width=18, height=11, res=300) #add file path
print(plot)
dev.off()
