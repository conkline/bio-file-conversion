#Emily Conklin 2018
import random

def imaToPhylip(filein, fileout):
	filein = open(filein, "r"); fileout = open(fileout, "w")
	lines = filein.readlines()
	pops = lines[2].split() #extracts population names from the top of file
	phyDict = dict()
	runnSeqLen = 0
	seq = ""
	pops_this_contig = []
	count = 0

	for line in lines[5::]:
		if not line.startswith(tuple(pops)):
			contig = line

			for key in phyDict.keys():
				if key not in pops_this_contig:
					phyDict[key] += "N" * len(seq)

			runnSeqLen += len(seq)
			pops_this_contig = []

		else:
			full_popname = line.split()[0]
			temp_popname = full_popname[:-1]

			if full_popname[-1] == "A":
				seqA = line.split()[1]
			elif full_popname[-1] == "B":
				pops_this_contig.append(temp_popname)
				seqB = line.split()[1]
				seq = random.choice([seqA, seqB])

				if temp_popname in phyDict.keys():
					phyDict[temp_popname] += seq
				else:
					phyDict[temp_popname] = ("N" * runnSeqLen) + seq

	#make sure we get last contig
	for key in phyDict.keys():
		if key not in pops_this_contig:
			phyDict[key] += "N" * len(seq)

	runnSeqLen += len(seq)

	fileout.write(str(len(phyDict.keys())) + " " + str(runnSeqLen) + "\n")
	for key in sorted(phyDict.keys()):
		fileout.write(key + "\t" + phyDict[key] + "\n")


imaToPhylip("example_IMA2_format_fixed.fasta", "example_phylip_fixed.phy")