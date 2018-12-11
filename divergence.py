#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio import codonalign
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import AlignIO
from Bio import pairwise2
from Bio.codonalign.codonseq import _get_codon_list, CodonSeq, cal_dn_ds
import scipy
from Bio.Align.Applications import MuscleCommandline
from Bio.Align import MultipleSeqAlignment
from scipy.linalg import expm
from Bio import SeqIO
import sys
import numpy as np
from Bio.SeqUtils import GC
import traceback
from docopt import docopt
#import pandas as pd




"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                %
%           Programme de calcul de Divergence entre paires de séquences     					 % 
%                    																			 %
%             ##INPUT###                                                                         %
%             Prend en entrée 4 fichiers fasta ordonnés par ordre alphabétique:                  %
%			  (Chaque fichier doit contenir la séquence homologue de l'autre fichier)            %
%             2 fichiers nucléotidiques.fa                                                       %
%             2 fichiers protéiques.fa               											 %
%     												 											 %
%			 ##	Ce que fait le programme ##													     %
%			-Crée un alignement de 2 séquences protéiques (alogrithme MUSCLE)    				 %
%			-Crée un alignement codon de 2 séquences nucléotidiques				      		     %
%																							     %
%            Analyse de l'alignement de codon:                                                   %
%                                                                                                %
%         	_Calcul la distance brute   														 % 
%.                                                                                               %
%			_Calcul la distance à la 3ieme position du codon   									 %
%                                  																 %
%			_Calcul la divergence synonyme (dN) et non-synonyme (dS):                            %
%                                                                                                %
%             4 méthodes :                           											 %
%																								 %
%			  - NG86  - `Nei and Gojobori (1986)`_ (PMID 3444411). 								 %
%	          - LWL85 - `Li et al. (1985)`_ (PMID 3916709). 									 %
%     		  - ML    - `Goldman and Yang (1994)`_ (PMID 7968486). 								 %
%   	      - YN00  - `Yang and Nielsen (2000)`_ (PMID 10666704). 							 %
% 																								 %
%        	  _`Nei and Gojobori (1986)`: http://www.ncbi.nlm.nih.gov/pubmed/3444411 			 %
%    	      _`Li et al. (1985)`: http://www.ncbi.nlm.nih.gov/pubmed/3916709 					 %
%     	      _`Goldman and Yang (1994)`: http://mbe.oxfordjournals.org/content/11/5/725 		 %
%             _`Yang and Nielsen (2000)`: https://doi.org/10.1093/oxfordjournals.molbev.a026236	 %                                               %
%                                                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


"""

if len(sys.argv) < 7:
    sys.exit("Usage: python3 divergence.py <sp1.dna file> <sp2.dna file> <sp1.aa file> <sp2.aa file> <outfile name> <method> <path to muscle.exe> ")

"""

Exemple de run:

python3 divergence.py concatenate_0035_fna_renamed.fst concatenate_0042_fna_renamed.fst concatenate_0035_faa_renamed.fst concatenate_0042_faa_renamed.fst dn_ds.out ML /Users/etudiant/Downloads/muscle3.8.31_i86darwin64

python3 divergence.py seq1_dna_test.fa seq2_dna_test.fa seq1_aa_test.fa seq2_aa_test.fa dn_ds.out ML /Users/etudiant/Downloads/muscle3.8.31_i86darwin64

help= Divergence

Usage:    dico.py <DNA_FILE_specie_1> 
		  dico.py <DNA_FILE_specie_2> 
          dico.py <AA_FILE_specie_1> 
          dico.py <AA_FILE_specie_2> 
          dico.py <OUPUT_FILE_dN_dS> 
          dico.py <METHOD> 

 
Options:
 <DNA_FILE_specie_1>	Fichier de séquence adn de l'espèce 1 format fasta
 <DNA_FILE_specie_2> 	Fichier de séquence adn de l'espèce 2 format fasta
 <AA_FILE_specie_1> 	Fichier de séquence aa de l'espèce 1 format fasta
 <AA_FILE_specie_2> 	Fichier de séquence aa de l'espèce 2 format fasta
 <OUPUT_FILE_dN_dS> 	Fichier de sortie des valeurs de dN er dS
 <METHOD> 				Méthode utilisée 


if __name__ == '__main__':
    arguments = docopt(help)
    print(arguments)

"""
def divergence():

	########################
	## Arguments d'entrée ##
	########################
	fic1dna=sys.argv[1] #fichier des séquences adn de l'espèce 1
	fic2dna=sys.argv[2] #fichier des séquences adn de l'espèce 2
	fic1prot=sys.argv[3] #fichier des séquences protéiques de l'espèce 1
	fic2prot=sys.argv[4] #fichier des séquences protéiques de l'espèce 2


	

	#outfile_unaligned="outfile_unaligned.fa"
	#outfile_unaligned=open(outfile_unaligned,"w",encoding='utf-8')
	outfile_dn_ds=sys.argv[5] #fichier de sortie format tableau, sep = ";"
	outfile_dn_ds=open(outfile_dn_ds,"w",encoding='utf-8')
	method=sys.argv[6] #Methode utilisée
	muscle_exe=sys.argv[7] #Chemin vers le fichier executable de MUSCLE

	#Transformation des séquences en format SeqIO
	seq1dna = list(SeqIO.parse(fic1dna, "fasta",alphabet=IUPAC.IUPACUnambiguousDNA()))
	seq2dna = list(SeqIO.parse(fic2dna, "fasta",alphabet=IUPAC.IUPACUnambiguousDNA()))
	seq1prot = list(SeqIO.parse(fic1prot, "fasta",alphabet=IUPAC.protein))
	seq2prot= list(SeqIO.parse(fic2prot, "fasta",alphabet=IUPAC.protein))
	


	#Première ligne du tableau "titres"
	"""print("seq.id",";","dN",";","dS",";","Dist_third_pos",";","Dist_brute",";","Length_seq_1",";","Length_seq2",
		";","GC_content_seq1",";","GC_content_seq2",";","GC",";","Mean_length",file=outfile_dn_ds)"""

	print("Nombre de paires de sequences a analyser: ",len(seq1dna))

	print("seq.id",";","dN",";","dS",";","Dist_third_pos",";","Dist_brute",";","Length_seq_1",";","Length_seq2",
		";","GC_content_seq1",";","GC_content_seq2",";","GC",";","Mean_length")



	"""df2 = pd.DataFrame(columns=("seq.id","dN","dS","Dist_third_pos","Dist_brute","Length_seq_1","Length_seq2",
		"GC_content_seq1","GC_content_seq2","GC","Mean_length"))"""



	#Boucle sur chaque paire de séquence
	u=0
	while u < (len(seq1dna)): 


		try:

			###########################################################
			#.    Alignement entre chaque paire de séquence           #
			###########################################################

		

			nuc1=str(seq1dna[u].seq)  #Récupère la séquence u et la transforme en string
			nuc2=str(seq2dna[u].seq)
			prot1=str(seq1prot[u].seq)
			prot2=str(seq2prot[u].seq)


			protein2 = SeqRecord(Seq(prot2,alphabet=IUPAC.protein), id='protein2') #Transformation de la séquence protéique en format SeqRecord
			protein1 = SeqRecord(Seq(prot1,alphabet=IUPAC.protein), id='protein1')

		

			with open("outfile_unaligned.fa", "w",encoding='utf-8') as output_handle: #Permet de créer un fichier de deux séquences non-alignées (format fasta)
				SeqIO.write(protein1, output_handle, "fasta")
				SeqIO.write(protein2, output_handle, "fasta")


			muscle_cline = MuscleCommandline(muscle_exe, input="outfile_unaligned.fa", out="outfile_aligned.aln") #Prend en entrée le fichier de séquences non-alignées et sort un fichier de séquences alignées
			stdout, stderr = muscle_cline()
			alns = AlignIO.read("outfile_aligned.aln", "fasta") #Lecture du fichier de séquences alignées 

			
			prot1=str(alns[0].seq) #Récupère la séquence protéique 1 alignée
			prot2=str(alns[1].seq) #Récup§re la séquence protéique 2 alignée

			nuc2 = SeqRecord(Seq(nuc2,alphabet=IUPAC.IUPACUnambiguousDNA()), id='nuc2') #Transformation de la séquence nucléique en format SeqRecord
			nuc1 = SeqRecord(Seq(nuc1,alphabet=IUPAC.IUPACUnambiguousDNA()), id='nuc1')


			prot1 = SeqRecord(Seq(prot1, alphabet=IUPAC.protein),id='pro1') #Transformation de la séquence protéique en format SeqRecord
			prot2 = SeqRecord(Seq(prot2, alphabet=IUPAC.protein),id='pro2')

			aln = MultipleSeqAlignment([prot1, prot2]) #Créer format alignement des 2 séquences protéiques préalablement alignées
			

			codon_aln = codonalign.build(aln, [nuc1, nuc2]) #Créer un alignement de codon 



			#Fichier d'alignement 
			#AlignIO.write(codon_aln,"outfile_aligned", 'fasta')

			lengthseq1=len(nuc1.seq)
			lengthseq2=len(nuc2.seq)
			GCcontentseq1=GC(nuc1.seq)
			GCcontentseq2=GC(nuc2.seq)

			GC_mean=((GCcontentseq1 + GCcontentseq2)/2)

			if lengthseq1 >= lengthseq2:
				Min_length=lengthseq2
			if lengthseq1 < lengthseq2:
				Min_length=lengthseq1
		
			##########################################################
			#           CALCULS DES INDICES DE DIVERGENCE            #      
			##########################################################

			#Calcul de divergence synonyme et non-synonyme 

			#Supression des gaps
			seq1=""
			seq2=""
			for x, z in zip(codon_aln[0],codon_aln[1]):
				if z =="-":
					continue
				if x =="-":
					continue
				else:
					seq1+=x
					seq2+=z

			#################################################################
			#.	        Comptage du nombre de site polymorhe brute          #
			#################################################################

			#Compteur de différences par site
			compteur0=0
			for i, e in zip(seq1,seq2):
				if i!=e:
					compteur0+=1

			distance_brute=round(float((compteur0)/len(seq1)),3)


			seq1_third_pos=""
			seq2_third_pos=""

			compteur1=0
			for i in seq1[2::3]: 
			    if  i.isalpha():
			    	seq1_third_pos+=i
			    	compteur1+=1

			compteur2=0
			for i in seq2[2::3]:
			    if  i.isalpha():
			    	seq2_third_pos+=i
			    	compteur2+=1

			####################################################################
			#	Comptage du nombre de site polymorphe en troisième position    #
			####################################################################

			#Compteur de différences par site (3ieme position)
			compteur3=0
			for i, e in zip(seq1_third_pos,seq2_third_pos):
				if i!=e:
					compteur3+=1

			distance_third_pos=round(float((compteur3)/compteur2),3)

			####################################################################
			#			Calcul dN et dS selon la méthode utilisée 			   #
			####################################################################

			try:
	        
				dN, dS = cal_dn_ds(codon_aln[0], codon_aln[1], method=method)

				"""print(seq1dna[u].id,";",dN,";",dS,";",distance_third_pos,";",distance_brute,";",lengthseq1,
					";",lengthseq2,";",GCcontentseq1,";",GCcontentseq2,";",GC_mean,";",Min_length,file=outfile_dn_ds)"""
				print(seq1dna[u].id,";",dN,";",dS,";",distance_third_pos,";",distance_brute,";",lengthseq1,
					";",lengthseq2,";",GCcontentseq1,";",GCcontentseq2,";",GC_mean,";",Min_length)

				"""df2=df2.append({"seq.id":seq1dna[u].id,"dN":dN,"dS":dS,"Dist_third_pos":distance_third_pos,"Dist_brute":distance_brute,"Length_seq_1":lengthseq1,
		"Length_seq2":lengthseq2,"GC_content_seq1":GCcontentseq1,"GC_content_seq2":GCcontentseq2,"GC":GC_mean,"Mean_length":Min_length}, ignore_index=True)"""

	

			except ValueError:
				result = 9.999 #Saturation trop importante pour calculer les indices.

				"""print(seq1dna[u].id,";",result,";",result,";",distance_third_pos,";",distance_brute,";",lengthseq1,
					";",lengthseq2,";",GCcontentseq1,";",GCcontentseq2,";",GC_mean,";",Min_length,file=outfile_dn_ds)"""
				print(seq1dna[u].id,";",result,";",result,";",distance_third_pos,";",distance_brute,";",lengthseq1,
					";",lengthseq2,";",GCcontentseq1,";",GCcontentseq2,";",GC_mean,";",Min_length)

				"""df2=df2.append({"seq.id":seq1dna[u].id,"dN":result,"dS":result,"Dist_third_pos":distance_third_pos,"Dist_brute":distance_brute,"Length_seq_1":lengthseq1,
		"Length_seq2":lengthseq2,"GC_content_seq1":GCcontentseq1,"GC_content_seq2":GCcontentseq2,"GC":GC_mean,"Mean_length":Min_length}, ignore_index=True)"""



			except ZeroDivisionError:
				result = 9.999
				
				"""print(seq1dna[u].id,";",result,";",result,";",distance_third_pos,";",distance_brute,";",lengthseq1,
					";",lengthseq2,";",GCcontentseq1,";",GCcontentseq2,";",GC_mean,";",Min_length,file=outfile_dn_ds)"""
				print(seq1dna[u].id,";",result,";",result,";",distance_third_pos,";",distance_brute,";",lengthseq1,
					";",lengthseq2,";",GCcontentseq1,";",GCcontentseq2,";",GC_mean,";",Min_length)

				"""df2=df2.append({"seq.id":seq1dna[u].id,"dN":result,"dS":result,"Dist_third_pos":distance_third_pos,"Dist_brute":distance_brute,"Length_seq_1":lengthseq1,
		"Length_seq2":lengthseq2,"GC_content_seq1":GCcontentseq1,"GC_content_seq2":GCcontentseq2,"GC":GC_mean,"Mean_length":Min_length}, ignore_index=True)"""



			except KeyError:
				result = 9.999
				
				"""print(seq1dna[u].id,";",result,";",result,";",distance_third_pos,";",distance_brute,";",lengthseq1,
					";",lengthseq2,";",GCcontentseq1,";",GCcontentseq2,";",GC_mean,";",Min_length,file=outfile_dn_ds)"""
				print(seq1dna[u].id,";",result,";",result,";",distance_third_pos,";",distance_brute,";",lengthseq1,
					";",lengthseq2,";",GCcontentseq1,";",GCcontentseq2,";",GC_mean,";",Min_length)

				"""df2=df2.append({"seq.id":seq1dna[u].id,"dN":result,"dS":result,"Dist_third_pos":distance_third_pos,"Dist_brute":distance_brute,"Length_seq_1":lengthseq1,
		"Length_seq2":lengthseq2,"GC_content_seq1":GCcontentseq1,"GC_content_seq2":GCcontentseq2,"GC":GC_mean,"Mean_length":Min_length}, ignore_index=True)"""



			u+=1

		except:
			traceback.print_exc()
			print("Une erreur est survenue pour la sequence: ", seq1dna[u].id, "vs", seq2dna[u].id)
			"""df2=df2.append({"seq.id":seq1dna[u].id,"dN":"NA","dS":"NA","Dist_third_pos":"NA","Dist_brute":"NA","Length_seq_1":"NA",
		"Length_seq2":"NA","GC_content_seq1":"NA","GC_content_seq2":"NA","GC":"NA","Mean_length":"NA"}, ignore_index=True)"""

			u+=1

	#df2.to_csv(outfile_dn_ds, sep='\t')
	outfile_dn_ds.close() #Fermeture du fichier ouvert 

if __name__ == '__main__':
	print(divergence())





