#!/bin/bash

sed -e 's/ /_/g' orf_finder_xyl1_hits.fasta > fixed1.fasta

sed -e 's/-/_/g' fixed1.fasta > fixed2.fasta

sed -e 's/|/_/g' fixed2.fasta > fixed3.fasta

sed -e 's/:/_/g' fixed3.fasta > fixed4.fasta

sed -e 's/__/_/g' fixed4.fasta > xyl1.hits.fixedIDs.fasta

rm fixed*

/opt/bifxapps/samtools/samtools faidx xyl1.hits.fixedIDs.fasta


for ID in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16

do

sed -e 's/ /_/g' clades/clade$ID\.txt > fixed1.txt

sed -e 's/-/_/g' fixed1.txt > fixed2.txt

sed -e 's/|/_/g' fixed2.txt > fixed3.txt

sed -e 's/:/_/g' fixed3.txt > fixed4.txt

sed -e 's/__/_/g' fixed4.txt > clade$ID\-fixedIDs.txt

rm fixed* 

	cat clade$ID\-fixedIDs.txt | while read line 
	do
  	
 	seq="$(grep "$line" xyl1.hits.fixedIDs.fasta)"
 	seq1="$(echo "$seq" | sed 's/>//g')"

 	/opt/bifxapps/samtools/samtools faidx xyl1.hits.fixedIDs.fasta $seq1 > $seq1\.fna

	done
	
	cat *.fna > xyl1_clade_$ID\.fasta

	rm *.fna

done
