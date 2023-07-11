# vigene_array_analysis
processing, analyzing, and plotting array data output from microvigene

1)Input vigene file with array information, join to array layout file, eliminate flagged spots. Join all data with experimental information file. 

2)Average spots on each array (group as needed, protein and concentration used here). Subtract spots from same protein/concentration on negative control chip and subtract from internal control (anti-FITC here). Propagate error across substractions. 

3)Plot binding vs dilution for protein of interest. With one-site binding model. 
