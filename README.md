
## HORIZONA is a pipeline coded essentially in python3.5 through the Horizon project (all python codes are made by hand and calling other packages already written). Its aim is to find horizontal gene transfers between two distantly related species.

### Prerequisites

What softwares you need to install:

* BUSCO V2 or v3 (BUSCO genes prediction): [BUSCO software](http://busco.ezlab.org/)      
* AUGUSTUS (genes prediction): [AUGUSTUS](http://bioinf.uni-greifswald.de/augustus/downloads/) 
* PYTHON 3 : [PYTHON 3](https://www.python.org/downloads/)
(use conda to install all dependency such as numpy, pandas, biopython etc.)
* Gffread (deal with augustus output): [Gffread](https://github.com/gpertea/gffread)
* Diamond (make a quick blastp to find homologous genes): [Diamond](https://github.com/bbuchfink/diamond)
* Silix (creat cluster within genes): [SiliX](http://lbbe.univ-lyon1.fr/Download.html) ask the authors for the version "silix-1.2.10-p1" 
* Seaview (build a tree): [SeaView](http://doua.prabi.fr/software/seaview)
* Muscle (alignment) : [Muscle](https://www.drive5.com/muscle/downloads.htm)

Databases: 
The nr database available here: (ftp://ftp.ncbi.nlm.nih.gov/blast/db/) (huge file).
For the taxid informations, see Peter Thorpe's instuctions here: [Instructions](https://github.com/peterthorpe5/public_scripts/tree/master/Diamond_BLAST_add_taxonomic_info)

### Run example

###For the process to be more accessible and understandable by everyone, we will name in this example the species 1: ```0035``` and the species 2: ```0042```

Let's run Busco to find conserved unique genes inside our 2 genomes: (please change informations inside these files to fit with your own data)
```
bash Busco_run_sp1.sh
bash Busco_run_sp2.sh
```
Outfiles:
* run_sp1_BUSCO_v2 
* run_sp1_BUSCO_v2

With inside these folders, the Busco sequences in AA and DNA formats "single_copy_busco_sequences", a short summary of Busco found and the retraining parameters for Augustus "/augustus_output/retraining_parameters" (we will use these informations later in the process). 

Now that we got the Busco output files, we will keep only completes and conserved sequences into both species.
```
python3 Order_paired_seq.py  -s1 0035 -s2 0042 -d1 /Buscou.out/single_copy_busco_sequences_0035/ -d2 Buscou.out/single_copy_busco_sequences_0042/ -t1 full_table_ACG-0035_BUSCO_v2.tsv -t2 full_table_ACG-0042_BUSCO_v2.tsv
```

Outfiles:
* sp1.faa (amino acide format)
* sp2.faa
* sp1.fna (nucleotide format)
* sp2.fna
* concatened_sp1.fst (contains the aa and dna sequences of sp1 in order)
* concatened_sp2.fst (contains the aa and dna sequences of sp2 in order)

Then, you will get 4 fasta files ready to be analyzed with the divergence.py program: 
```
python3 divergence.py -f1 0035.fna -f2 0042.fna -f3 0035.faa -f4 0042.fna -m ML -a /Users/etudiant/Downloads/muscle3.8.31_i86darwin64 -o dn_ds_Busco.out 
```
Outfiles:
* dn_ds_Busco.out  
#### Table format

| seq.id | dN | dS | Dist_third_pos | Dist_brute | Length_seq_1 |Length_seq2 | GC_content_seq1 | GC_content_seq2 | GC | Mean_length |
| ------ |:--:| :-:| :-------------:| :---------:| :-----------:| :---------:| :--------------:| :--------------:| :-:| -----------:|


Where:
* seq.id : Busco seq name 
* dN : Non-synoymous distance 
* dS : Synonymous distance 
* Dist_third_pos : Distance at the third codon position 
* Dist_brute : Distance at all sites 
* Length_seq : Length of the sequence 
* GC_content_seq : GC content of the sequence 
* GC : GC mean between the two sequences 
* Mean_length : The smallest length 


#Here you got the distances of Busco paired sequences with a Maximum Likelihood (ML) method.


Now, we will have to try to find all possible genes present in our 2 genomes with Augustus: (please change information inside these files to fit with your own data)
```
bash augustus_run_sp1_training_sp1.sh #outfile : run_augustus_sp1_training_sp1.out (gff format)
bash augustus_run_sp1_training_sp2.sh #outfile : run_augustus_sp1_training_sp2.out (gff format)
bash augustus_run_sp2_training_sp2.sh #outfile : run_augustus_sp2_training_sp2.out (gff format)
bash augustus_run_sp2_training_sp1.sh #outfile : run_augustus_sp2_training_sp1.out (gff format)
```

Now that you got the Augustus ouptufiles in gff format, you will have to get the fasta format of the CDS sequences (takes the phase into account as expected -x option and so does the -y option for protein).
* -g ACG-0035-scaffold.fa  is the genome file of the sp1 for exemple.

```
/Users/etudiant/Downloads/gffread-0.9.9-0/bin/gffread /Users/etudiant/Desktop/Horizon_project/Augustus/Augustus_out/run_augustus_sp1_training_sp1.out -g /Users/etudiant/Desktop/Horizon_project/Buscou.out/ACG-0035-scaffold.fa -x run_augustus_sp1_training_sp1.out.fna -y run_augustus_sp1_training_sp1.out.faa

/Users/etudiant/Downloads/gffread-0.9.9-0/bin/gffread /Users/etudiant/Desktop/Horizon_project/Augustus/Augustus_out/run_augustus_sp1_training_sp2.out -g /Users/etudiant/Desktop/Horizon_project/Buscou.out/ACG-0035-scaffold.fa -x run_augustus_sp1_training_sp2.out.fna -y run_augustus_sp1_training_sp2.out.faa

/Users/etudiant/Downloads/gffread-0.9.9-0/bin/gffread /Users/etudiant/Desktop/Horizon_project/Augustus/Augustus_out/run_augustus_sp2_training_sp2.out -g /Users/etudiant/Desktop/Horizon_project/Buscou.out/ACG-0042-scaffold.fa -x run_augustus_sp2_training_sp2.out.fna -y run_augustus_sp2_training_sp2.out.faa

/Users/etudiant/Downloads/gffread-0.9.9-0/bin/gffread /Users/etudiant/Desktop/Horizon_project/Augustus/Augustus_out/run_augustus_sp2_training_sp1.out -g /Users/etudiant/Desktop/Horizon_project/Buscou.out/ACG-0042-scaffold.fa -x run_augustus_sp2_training_sp1.out.fna -y run_augustus_sp2_training_sp1.out.faa
```

Outfiles are the -x and -y parameters

Now we will have to find homologous sequences into all these predicted sequences, to do so, first we will have to concatenate all sequences into one AA file. 
#But first, please name and tag your sequences with their respective number to recognize them at the end of the process. To do so:
#Add the specie number after all ".t1" patterns, please replace your own number in the code below.

```
sed 's/\(.t1\)/\1_0035_0035/' run_augustus_sp1_training_sp1.out.faa > augustus_sp1_training_sp1.out.faa 
sed 's/\(.t1\)/\1_0035_0042/' run_augustus_sp1_training_sp2.out.faa > augustus_sp1_training_sp2.out.faa
sed 's/\(.t1\)/\1_0042_0042/' run_augustus_sp2_training_sp2.out.faa > augustus_sp2_training_sp2.out.faa
sed 's/\(.t1\)/\1_0042_0035/' run_augustus_sp2_training_sp1.out.faa > augustus_sp2_training_sp1.out.faa

sed 's/\(.t1\)/\1_0035_0035/' run_augustus_sp1_training_sp1.out.fna > augustus_sp1_training_sp1.out.fna
sed 's/\(.t1\)/\1_0035_0042/' run_augustus_sp1_training_sp2.out.fna > augustus_sp1_training_sp2.out.fna
sed 's/\(.t1\)/\1_0042_0042/' run_augustus_sp2_training_sp2.out.fna > augustus_sp2_training_sp2.out.fna
sed 's/\(.t1\)/\1_0042_0035/' run_augustus_sp2_training_sp1.out.fna > augustus_sp2_training_sp1.out.fna
```
Concatenate these 4 files into one:
```
cat augustus_sp1_training_sp1.out.faa augustus_sp1_training_sp2.out.faa augustus_sp2_training_sp2.out.faa augustus_sp2_training_sp1.out.faa > concatenate_0035_0042.faa 

cat augustus_sp1_training_sp1.out.fna augustus_sp1_training_sp2.out.fna augustus_sp2_training_sp2.out.fna augustus_sp2_training_sp1.out.fna > concatenate_0035_0042.fna
```
Outfiles:
* concatenate_0035_0042.faa 
* concatenate_0035_0042.fna

#Then, we will perform a Diamond "all against all" and a SiLix run on these amino acide sequences to find homologous sequences and then, to only keep the best ones within each cluster :

#Let's make a Blast database of or sequences:

DIAMOND
```
Diamond=/Users/etudiant/Downloads/diamond-master/bin/diamond
```
Database (-d is the output db file in .dmnd extention)
```protein_fasta=concatenate_0035_0042.faa 
$Diamond makedb --in $protein_fasta -d Augustus_diamond_0035_0042
```

Let's now make a Blastp all against all with Diamond:

#DIAMOND
```Diamond=/Users/etudiant/Downloads/diamond-master/bin/diamond
nr=Augustus_diamond_0035_0042.dmnd
out=matches_Augustus_0035_0042.m8
```

#Database
```protein_fasta=concatenate_0035_0042.faa
$Diamond blastp -d $nr -q $protein_fasta  -o $out
```
Finnaly, we'll have to get cluster of these sequences with Silix:
#SOFTWARE
```
SILIX=/Users/etudiant/Downloads/silix-1.2.10-p1/src/silix
```

#DATASET USED
```BLAST=matches_Augustus_0035_0042.m8 
FASTA=concatenate_0035_0042.faa

$SILIX  $FASTA $BLAST -f cluster_ -n > cluster_Augustus_0035_0042.fnodes
```

Now we only want to get the 2 paired sequences into cluster with the best percentage of identity (pident) of the mean High-scoring Segment Pair (HSP): 
```
python3 cluster_silix_merging.py -b matches_Augustus_0035_0042.m8 -c cluster_Augustus_0035_0042.fnodes -d matches_Augustus_0035_0042.net -f1 concatenate_0035_0042.faa -f2 concatenate_0035_0042.fna 
```
Outfiles: These files are the files were are only present homologous sequences whith the best pident within each cluster. In order. 
* clusters1_aa.fasta
* clusters2_aa.fasta
* clusters1_dna.fasta
* clusters2_dna.fasta

Then, we will calculate the distances values of all these predicted paired genes: (the fact to add -names will re-organise the dN_dS output tab)
```
python3 divergence.py -f1 clusters1_dna.fasta -f2 clusters2_dna.fasta -f3 clusters1_aa.fasta -f4 clusters2_aa.fasta -m ML -a /Users/etudiant/Downloads/muscle3.8.31_i86darwin64 -n1 0035 -n2 0042 -o dn_ds_Augustus.out
```
Outfile:
* dn_ds_Augustus.out 
#### Table format

| seq1.id | seq2_id | dN | dS | Dist_third_pos | Dist_brute | Length_seq_1 |Length_seq2 | GC_content_seq1 | GC_content_seq2 | GC | Mean_length |
| ------- |:-------:|:--:|:-:| :-------------:| :---------:| :-----------:| :---------:| :--------------:| :--------------:| :-:| -----------:|

Where seqn_id is the name of the sequence n.

We will only keep within all these paired sequences the ones wich have passed through the dS and the Cov thresholds:
```
python3 gff_cov.py -d1 dn_ds_Busco.out -d2 dn_ds_Augustus.out  -s1 0035 -s2 0042 -c1 cov_GC_0035.tab -c2 cov_GC_0042.tab -g1 run_augustus_0035.out -g2 run_augustus_0035_training_0042.out -g3 run_augustus_0042.out -g4 run_augustus_0042_training_0035.out -g5 gff_Busco_0035 -g6 gff_Busco_0042
```

Outfiles:
* candidates_genes_pvalue_cov_0035_0042.tab (Final output with candidates genes predicted by Augustus after passing through the filter p-value and coverage)
* Augustus_clustering.tab (Output table of Augustus predicted genes after passing through the clustering filter (max pident within each cluster)
* candidates_genes_pvalue_0035_0042.tab (Final output with candidates genes predicted by Augustus after passing through the filter p-value:)

Final fasta output with candidates genes predicted by Augustus after passing through the filter dS (this file will be usefull for the blast against the nr db step).
* candidates_aa_pvalue_0035.fasta
* candidates_aa_pvalue_0042.fasta
* candidates_dna_pvalue_0035.fasta
* candidates_dna_pvalue_0042.fasta

#The following steps will have to be done with a program developped by Peter Thorpe, all needed files are available here : (https://github.com/peterthorpe5/public_scripts/tree/master/Diamond_BLAST_add_taxonomic_info)

Now that we got the candidates genes, to characterize them, we will perform a blastp against the nr database: (please change information inside these files to fit with your own data)
Note that this will take into account all genes passed through the dS filter only (if you only want those with dS + cov, please change the input fasta file)
```
bash Diamond_blastp_sp1_candidates.sh
bash Diamond_blastp_sp2_candidates.sh
```

Outfiles:
* matches_sp1_candidates.m8 (blast output table of the sp1)
* matches_sp2_candidates.m8 (blast output table of the sp2)

Finnaly to be more readable, we can get taxid informations about these HGT candidates genes.(cf Peter Thorpe program above).
```
bash tax_name_sp1.sh
bash tax_name_sp2.sh
```

Output files for sp1:
* outfile_0035.tab (Blast output table with taxid informations)
* outfile_0035.tab_top_blast_hits.out (Best blast output table with taxid informations)
* outfile_0035.tab_top_blast_hits.out_based_on_order_tax_king.tab
* outfile_0035.tab_top_blast_hits.out_histogram.png (Summary plot of the blast result)

You can use this script to also include sequence with no blast hit and get a final table with all information: .
```
python3 get_dataframe_candidates.py -c candidates_genes_pvalue_0035_0042.tab -b1 outfile_0035.tab_top_blast_hits.csv -b2 outfile_0042.tab_top_blast_hits.csv -s1 0035 -s2 0042
```
Output files: 
* All_candidates_cov_all_inf.tab (Table with all information with all candidats passed through the dS and covergae threshold).
* Summarize_All_candidates_with_all_inf.tab (Table with more important information with all candidats passed through the dS threshold).
* All_candidates_with_all_inf.tab (Table with all information with all candidats passed through the dS threshold).

### Table content of the summary file:
| seq1.id | seq2.id | dN | dS | Dist_third_pos | Dist_brute | Mean_length| Length_seq_1 |Length_seq2 | GC_content_seq1 | GC_content_seq2 | caf_name_seq1 | GC_scaff_seq1scaf_name_seq1 | cov_depth_seq1 |scaf_name_seq2 | GC_scaff_seq2 | cov_depth_seq2 | sseqid_sp1 | evalue_sp1 |
| salltitles_sp1 | scientific_name_sp1 | Order_sp1	sseqid_sp2 | sseqid_sp2 | evalue_sp2 | salltitles_sp2 | scientific_name_sp2 | Order_sp2 |



Now that you got all the candidate genes, you will have to do multiples thing to check if your candidates are real HGT and integrated into the other specie's genome. 

## Here are the steps:

1 - Construct a phylogenetic tree and see if there is an incongruence (It does mean that you need a resolved phylogenetic tree of your species built with ubiquitous genes (ex: Busco genes) for instance):
1.1 - In order to creat a phylogeny of your gene, you will need to have several sequences coming from several species closely related to your two species. To do so take your target sequence in amino acide format and make a research of homologous sequences against the nr databe and pick up around 15-20 sequences for both species.
1.2 - Align all these sequences with Seaview (Gouy M et al.,2010 available here : http://doua.prabi.fr/software/seaview) for exemple, align all these sequences and use Gblocks program to select only blocks of evolutionarily conserved sites.
1.3 - Construct your tree with a Maximum likelihood method (or a baysian method, not availabe on seaview but you can use MrBaye) and add boostrap analysis (or a LRT if you want to be faster).
1.4 - If there is an incongruence and a good bootstrap support, then your sequence is probably an HGT since the two species are diverging from a long time. 

2 - Check if this sequence is a contamination, not integrated one or if it is incorporated into the host genome:
2.1- If this sequence is found in a scaffold where you have found other predicted genes or even BUSCO ones, it is a good clue for an integrated HGT sequence.
2.2- Check if the 2 sequences have some differences at the nucleotid level, if they do not and are 100% identical, it can be a contamination.


## Contributors 
-GUINET Benjamin
-VARALDI Julien
-CHARLAT Sylvain	

## Acknowledgments
This work was performed using the computing facilities of the CC LBBE/PRABI. 	
Heartfelt thanks to all LBBE contributors for their help during the realization of this pipeline.	

