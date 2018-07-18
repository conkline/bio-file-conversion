import numpy as np

filein = open("Nuclear_Haplotypes_FINAL.ima", "r")
fileout = open("Nuclear_Haplotypes_FINAL.nex", "w")
pop_dict = {}
tmp_dict = {}
lenSeq = 0

for line in filein.readlines()[5:]:

	if line.startswith("dDocent"):

		if len(tmp_dict.keys()) > 0:

			sample_cons = {}
			alis = []
			samples = sorted(tmp_dict)
			for sample in samples:

				sample_cons[sample] = ""
				for pos in range(len(tmp_dict[sample][0])):
					#compare haplotypes at each base - see if homozygous or heterozygous
					hapA_base = tmp_dict[sample][0][pos]
					hapB_base = tmp_dict[sample][1][pos]

					#if homozygous, keep nucleotide as is for now
					#if heterozygous, change to 1
					if hapA_base != hapB_base:
						sample_cons[sample] += "1"
					else:
						sample_cons[sample] += hapA_base
				alis.append(list(sample_cons[sample]))

			#join consensus sequences together
			alignment = np.array(alis)

			for ali_pos in range(0,alignment.shape[1]):
				col = alignment[:,ali_pos]
				unique_bp = np.unique(col)

				#see if we have a variable site
				if len(unique_bp) > 1:
					#if the only variations are the sites that have heteros
					#set all other values to 1
					no_hets = [x for x in unique_bp if x != "1"] #no hets...a perfect world
					if len(no_hets) == 1:
						alignment[:,ali_pos] = np.array(["1" if x=="1" else "0" for x in list(col)])
					#otherwise, encode different homos as 0 or 2
					else:
						encoding = ['0','2']
						for i in range(0,2):
							alignment[:,ali_pos][col==no_hets[i]] = encoding[i]

				#if we don't have a variable site, convert to NaN
				else:
					alignment[:,ali_pos] = np.array([np.nan for x in list(col)])

			#remove non-variable sites
			ix = ~np.isin(alignment, "n")
			rows, cols = np.where(ix)
			alignment = alignment[:,list(np.unique(cols))]

			#add SNPs to master dictionary
			for sam in range(len(samples)):
				if samples[sam] in pop_dict.keys():
					pop_dict[samples[sam]] += "".join(list(alignment[sam,:]))
				else:
					pop_dict[samples[sam]] = ("?" * lenSeq) + "".join(list(alignment[sam,:]))

			#for samples that were not in this contig, add in missing data
			for miss in [x for x in pop_dict.keys() if x not in samples]:
				pop_dict[miss] += ("?" * len("".join(list(alignment[0,:]))))

			lenSeq += len("".join(list(alignment[0,:])))

		#set up new contig
		tmp_dict = {}

	else:
		parts = line.split()
		sample = parts[0][0:-1]
		sample_hap = parts[0]
		sample_seq = parts[1].strip("\n")

		if sample_hap[-1] == "A":
			tmp_dict[sample] = [sample_seq]
		elif sample_hap[-1] == "B":
			tmp_dict[sample].append(sample_seq)

#get last alignment
sample_cons = {}
alis = []
samples = sorted(tmp_dict)
for sample in samples:

	sample_cons[sample] = ""
	for pos in range(len(tmp_dict[sample][0])):
		#compare haplotypes at each base - see if homozygous or heterozygous
		hapA_base = tmp_dict[sample][0][pos]
		hapB_base = tmp_dict[sample][1][pos]

		#if homozygous, keep nucleotide as is for now
		#if heterozygous, change to 1
		if hapA_base != hapB_base:
			sample_cons[sample] += "1"
		else:
			sample_cons[sample] += hapA_base
	alis.append(list(sample_cons[sample]))

#join consensus sequences together
alignment = np.array(alis)

for ali_pos in range(0,alignment.shape[1]):
	col = alignment[:,ali_pos]
	unique_bp = np.unique(col)

	#see if we have a variable site
	if len(unique_bp) > 1:
		#if the only variations are the sites that have heteros
		#set all other values to 1
		no_hets = [x for x in unique_bp if x != "1"] #no hets...a perfect world
		if len(no_hets) == 1:
			alignment[:,ali_pos] = np.array(["1" if x=="1" else "0" for x in list(col)])
		#otherwise, encode different homos as 0 or 2
		else:
			encoding = ['0','2']
			for i in range(0,2):
				alignment[:,ali_pos][col==no_hets[i]] = encoding[i]

	#if we don't have a variable site, convert to NaN
	else:
		alignment[:,ali_pos] = np.array([np.nan for x in list(col)])

#remove non-variable sites
ix = ~np.isin(alignment, "n")
rows, cols = np.where(ix)
alignment = alignment[:,list(np.unique(cols))]

#add SNPs to master dictionary
for sam in range(len(samples)):
	if samples[sam] in pop_dict.keys():
		pop_dict[samples[sam]] += "".join(list(alignment[sam,:]))
	else:
		pop_dict[samples[sam]] = ("?" * lenSeq) + "".join(list(alignment[sam,:]))

#for samples that were not in this contig, add in missing data
for miss in [x for x in pop_dict.keys() if x not in samples]:
	pop_dict[miss] += ("?" * len("".join(list(alignment[0,:]))))

lenSeq += len("".join(list(alignment[0,:])))

#write data to file
fileout.write("#NEXUS\nBEGIN TAXA;\n  DIMENSIONS NTAX="+str(len(pop_dict.keys()))+";\n")
fileout.write("  TAXLABELS\n")
for sample in sorted(pop_dict):
	fileout.write("    "+sample+"\n")
fileout.write("    ;\nEND;\nBEGIN CHARACTERS;\n")
fileout.write("  DIMENSIONS NCHAR="+str(lenSeq)+";\n  FORMAT\n    DATATYPE=SNP\n    MISSING=?;\n  MATRIX\n")

for sample in sorted(pop_dict):
	fileout.write("    "+sample+"\t"+pop_dict[sample]+"\n")

fileout.write("   ;\nEND;\nBEGIN SETS;\n")
fileout.write("  TaxSet Pop_1 = 1-"+str(len(pop_dict.keys()))+";\n")
fileout.write("  TaxPartition populations=\n    Pop_1:1-"+str(len(pop_dict.keys()))+";\nEND;\n")

fileout.close(); filein.close()