#PBS -S /bin/bash
#PBS -l nodes=1:ppn=8
#PBS -e /pandata/bguinet/LEPIWASP/blast_database/Diamond_blastp_sp2_candidate.error
#PBS -o /pandata/bguinet/LEPIWASP/blast_database/Diamond_blastp_sp2_candidate.out
#PBS -q q1day
#PBS -N diamond_blastp


diamond=/panhome/bguinet/miniconda3/bin/diamond 
nr=/panbanques/diamond/nr/nr.dmnd
candidates_aa_sp2=candidates_aa_pvalue_sp2.fasta

$diamond blastp -d $nr -q $candidates_aa_sp2 -o /pandata/bguinet/LEPIWASP/blast_database/matches_sp2_candidates.m8

