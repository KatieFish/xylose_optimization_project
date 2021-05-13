Contents of scripts folder:

Markdowns: most recent markdowns. Mostly deal with CO pipeline - NOT phylogenetic or bioinformatic steps. 

Xylose_optim_project-KJF-spring2021.Rmd

Xylose_optim_project-KJF-spring2021.html

Scripts: 

Processing_xylose_MSA_tsvs-KJF-spring2021.R - R script that reads in a multiple sequence file (fasta file) and converts to a dataframe in R in which rows are seqeunces and columns are base positions of sequene. Gaps are converted to NAs. 

Processing_xylose_MSA_tsvs.R - older version of above script. Difference is in how sequence names are processed. Use newer version. 

Pull_out_subgroup_cds_from_orffidner.sh - script I wrote to pull out the corresponding coding sequences from orfFinder annotations.

RLN_raxml.submit - script that Rishitha can modify and run on scarcity to produce ML trees. Set for PROTGAMMAAUTO and 100 bootstraps. 
*Rishitha - remember to update the script fully by reading my comments and update the arguments I specify. 
Then to run, simply run condor_submit RLN_raxML.submit (or whatver you name it). You'll know it's done by either looking at the output file 
or by waiting until all tree files are generated*

The scripts below shouldn't be needed at all: 
Quality_control_homologs_in_subgroups-V2.R - Initial attempts at QC - we abandoned after moving to orfFinder. 
Quality_control_homologs_in_subgroups.R
Round1_treeprocessing_for_iTOL.R
Xylose_optim_project.Rmd
Xylose_optim_project.html
convert_MSA_to_matrices.R
