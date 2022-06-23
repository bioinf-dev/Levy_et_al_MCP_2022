import sys
import pandas as pd
import numpy as np

dir = sys.argv[1]

#inputs
pepfile=dir + "/peptides.txt"
msfile=dir + "/msms.txt"

#############################################################################################################
# This program takes MaxQuant peptides and msms txt outputs, and prepares input files for other programs    #
# (for example for MS2 %fragmentation calculation, peptide ambiguity calculation,                           #
#  retention time outlier identification and etc.)                                                          #
#############################################################################################################


#outputs
out_pep_list=dir + "/peptide.list.full.txt"
out_fasta=dir + "/peptide.fasta"
out_pep_intens=dir + "/pep.identifications.intensities.txt"
out_ambiguity=dir + "/best.msms.for.mimcry.txt"
out_fragmentation=dir + "/best.msms.for.fragmentation.txt"
out_rt=dir + "/best.msms.for.RT.txt"
out_il_fasta=dir + "/il.fasta"
out_msms_pep_and_score=dir + "/msms_pep_and_score.txt"
out_prot_list=dir + "/pep.prot.txt"

#get peptide list and save to file
df0 = pd.read_csv(pepfile,sep='\t',dtype=str,header=(0))


cond1=df0['Reverse'].isnull()
cond2=df0['Potential contaminant'].isnull()
df=df0[ cond1 & cond2 ]

seq_prots = df[["Sequence","Proteins"]]
seq_prots.to_csv(out_prot_list, sep='\t', index=False, header=False)

df1 = df[["Sequence"]]
df1.to_csv(out_pep_list, sep='\t', index=False, header=False)

#convert peptide list to fasta
n = 0
with open(out_pep_list, 'r') as f: #read
	with open(out_fasta, 'w') as out: #write
		for line in f:
			out.write('>' + line.strip() + '\n' + line.strip() + '\n')




#convert peptide list to fasta I/L permutations
# dict with replacements
il_dict = {'I': 'L', 'L': 'I'}
with open(out_pep_list, 'r') as f: #read
	with open(out_il_fasta, 'w') as out: #write
		for line in f:
			pep=line.strip()
			all_perm = []  # list to hold all current peptide permutations
			all_perm.append(pep)  # add original peptide
			aa = np.array(list(pep))  # break string to list of amino acids
			il_index = np.isin(aa, ['I', 'L']).nonzero()[0]  # find index of I\L
			for i in il_index:  # for each I or L aa
				i_perm = []  # list of this i permutations
				for perm in all_perm:  # iterate over existing permutations
					# insert the current replacment of i aa, prefix of permutated pep and suffix of original pep
					i_perm.append(perm[:i] + il_dict[aa[i]] + pep[i + 1:])
				all_perm = all_perm + i_perm  # append all i permutation to peptide permutations

			all_perm = pd.Series(list(set(all_perm)))  # to series
			non_wt_perm=all_perm[all_perm!=pep]
			for nwt in non_wt_perm:
				out.write('>' + pep + '_' + nwt + '\n' + nwt + '\n')




pep_intens=df.filter(regex=r'^(Sequence|Intensity |Identification )', axis=1)
pep_intens.to_csv(out_pep_intens, sep='\t', index=False, header=True)

pep = pd.DataFrame()
pep["Sequence.id"] = df['Sequence'].map(str) + '.' +df['Best MS/MS'].map(str)

ms0 = pd.read_csv(msfile,sep='\t',dtype=str,header=(0))
ms0["RID"]= ms0['Raw file'].map(str) + '@' +ms0['id'].map(str)
cond=ms0['Reverse'].isnull()
ms=ms0[cond]

msms_pep_and_score=ms[["RID","Score","PEP"]] #select needed columns
msms_pep_and_score.to_csv(out_msms_pep_and_score, sep='\t', index=False, header=True)


bestms=ms #pd.merge(pep,ms,on='Sequence.id')

rt=bestms[["Raw file", "RID","Sequence","Retention time"]] #select needed columns
rt.to_csv(out_rt, sep='\t', index=False, header=True)

amb=bestms[["RID", "Sequence","Score", "All scores", "All sequences"]]
amb.to_csv(out_ambiguity, sep='\t', index=False, header=False)

frag=bestms[["RID", "Sequence", "Matches"]]
frag.to_csv(out_fragmentation, sep='\t', index=False, header=False)

