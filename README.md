# xylose_optimization_project

Rishitha summer '21 - 

All of the data you should need (sequences pulled out by HMMER, alignments and fast trees for all sequences, cds sequences of XYL gene groups, amino acid alignments of XYL gene proteins, & codon aware alignments of XYL gene cds sequences are in the Rishitha_summer_21 directory. The readme file explains what all the files are in case you forget how we named things. 

In the scripts folder you'll find all the scripts and the most recent markdown that we used to push the data through our R pipeline. 

I think the outstanding work to be done first, however, is curating our set of XYL homologs and spending some time understanding their patterns of evolution. Abbe will certainly have much better thoughts on this than I do. 

I think the weakness of our ORF-finder pipeline was that it pulled out more paralogs per genome than I think are really there. This might be because it pulls out fragments of genes that are there because of assembly errors. I think you should spend the first few weeks building good maximum liklihood XYL trees and using these to identify sequences that don't belong. Maybe Abbe has some thoughts on how to best do this. 

*A note of caution here - the SOR1 and SOR2 genes of S. cerevisiae cluster within the XYL2 clade on our larger tree. I think we need to exercise some caution with removing genes from XYL2 until we understand the evolutionary relationships between these genes.*

Once you have a list of sequences that you are confident in - I think the next step would be again to build ML trees with these sequences and look at tree topology and compare it to the whole genome 332 topology in the Shen et al. paper. Can we test the gene trees against the species tree and find anything interesting? What is going on with XYL gene evolution in the Saccharomycetaceae?

I put a modifiable RAxML script in the /scripts directory. Rishitha should be able to run this script on scarcity, download her tree files (remember the scp command!) and visualize them on iTOL. 
