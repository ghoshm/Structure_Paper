# Structure_Paper

Code for Ghosh & Rihel (2019). 

For a detailed description please see: 
https://www.biorxiv.org/content/10.1101/694471v1 

Data is avaiable at: 
https://zenodo.org/record/3344770#.XYYowShKiUm

For questions & queries please contact: 
marcus.ghosh.11@ucl.ac.uk 

Flow diagram depicting the steps of our analysis framework: 

![Analysis_Framework](https://user-images.githubusercontent.com/26411096/65374080-9a581f80-dc7d-11e9-9230-55e011ccab18.png)

Data is output from our behavioural set-up (ViewPoint) in the form of a .xls file. perl_batch_192.m organises this data to a .txt format. Experiment metadata (e.g. animal genotypes) is supplied in the form of a .txt file. The 1min bin method uses sleep_analysis2.m to produce figures and statistics from these two .txt files. The 25Hz method exports .raw data from ViewPoint to produce .xls files. Vp_Extract.m reorganises these, using .txt data, to a .mat file which can be input to either Vp_Analyse.m or Bout_Clustering.m. Vp_Analyse.m produces figures and statistics. Bout_Clustering.m uses the clustering function gmm_sample_ea.m to assign data to modules, produce figures and calculate statistics, Bout_Clustering.m’s output can be input to Bout_Transitions.m, which compresses full modular sequences by calling Batch_Compress.m and Batch_Grammar_Freq.m. The motifs identified from this approach can be input to Batch_Transitions_Hours.m which compresses 500 module chunks and uses Batch_Grammar_Freq.m to count motif occurrences per hour. With the exception of the 1min bin method (sleep_analysis2.m), two example figures are shown for each figure producing step. All code can be run locally, though for speed several steps (indicated in green) are best run on a cluster computer.
