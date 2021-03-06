


FILES present: 
________________________________________

hmmer_searches directory: 

alnmnts_for_hmms directory: 

.hmm files used for hmmer queries

hmmer_raw_results directory:

orf_finder_xyl*_hits.fasta - all pep annotations identified as significant in hmmer search

orf_finder_xyl*_hits-aln-cleaned1.fasta	- all significant hmmer hits aligned with mafft and with sites with less than 5% occpancy removed. Alignment used to build Fasttree. 

all_hmmer_result_fastTrees directory:

alignments and trees built using fastTree. 
-------------------------------------------
xyl_gene_subgroups directory: 

xyl*_subgroup.txt - file of ORFfinder sequence IDs of pep sequences of monophyletic group of hmmer hits that recieved xyl* kegg annotations

xyl*_subgroup.fasta - fasta file of pep sequences of monophyletic group of hmmer hits that recieved xyl* kegg annotations 

xyl*_subgroup.aln - mafft aligned xyl* pep sequences

xyl*_subgroup-cds.fasta - corresponding cds sequences of 
xyl*_subgroup.fasta file

xyl*_subgroup-codon_aligned_cds.fasta - pal2nal codon aligned cds file

xyl*_subgroup-codon_aligned_cds-cleaned.fasta - codon aligned cds file with codons <80% occupancy removed 

xyl*_subgroup-codon_aligned_cds-cleaned.txt - tab deliminated text file of cds sequences that can be imported into R as dataframe. Gaps are NAs. 

-------------------------

data directory: 

*******Note - the figshare assemblies have naming differences from the 332 tree and various files generated by Labella et al and D. Opulente data. I have tried to rename all files according to our data (the figshare assembly names). I know these are old species names and we can fix them all at the end. 

**Note - wi values, genome-wide stAI values, and s values from Labella et al 2019 are located in /fall_2020/data/labella_et_al/. I always have to unzip the genome_wide_tAI_all_spp.txt.gz file before use and rezip it before git pushing. 


files either used for or generated by CO pipeline

- Processing_xylose_MSA_tsvs-KJF-spring2021.R - processes the cds sequence files and growth data files 
- Xylose_optim_project-KJF-spring2021.Rmd - runs analysis pipeline and generates figures 

XYLpthwy_gene_presence_absence_matrix.txt - binary matrix for gene presence after eliminating identical sequences of the same species


spp_by_gene_stAI_vals.txt - stAI vals for each gene for each species

spp_by_gene_estAI_vals.txt - estAI vals for each gene for each species

spp_by_gene_maximum_paralog_stAI_vals.txt - max stAI val for each species

spp_by_gene_maximum_paralog_estAI_vals.txt - max estAI val for each species

spp_in_90th_percentile_codon_opt.txt - file of all species names appearing in >=90th percentile of estAI vals for all 3 genes

growth_data_taxon_name_switches.txt	- file that contains nomenclature disparities between Dana's growth data and our data

growth_rate_master_df.txt - tab deliminated file of all species, xylose growth rate, and max estAI vals for all genes.

tree_to_seq_data_nameswaps.txt - file of nomenclature disparities between our data and the 332 tree. 


iTOL_files

332_tree.nwk - 332 tree with names changed to match our data (figshare assembly names)

Spring2021_growthrate_orthopresence_iTOL_tree.pdf 

XYLpthwy_gene_presence_absence_matrix-2.txt - iTOL annotation file of gene presence/absence 

Xylose_growth_rate_matrix-nonzero_growthrates.txt - iTOl annotation file of growth rates from D. Opulente. 