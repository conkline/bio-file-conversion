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

	#iterates through rest of file
	for line in lines[5::]:
		if not line.startswith(tuple(pops)):

			#if any population wasn't in the last contig, add N-string of the proper length
			for key in phyDict.keys():
				if key not in pops_this_contig:
					phyDict[key] += "N" * len(seq)

			#update total sequence length and count of pops per contig
			runnSeqLen += len(seq)
			pops_this_contig = []

			#randomly choose between hapA and hapB for this contig
			hapChoose = random.choice(["A", "B"])

		else:
			full_popname = line.split()[0]
			temp_popname = full_popname[:-1]

			#if we're using HapA
			if full_popname[-1] == "A" and hapChoose == "A":
				seq = line.split()[1]

				if temp_popname in phyDict.keys():
					phyDict[temp_popname] += seq #if it's in the dictionary, add the new sequence
				else:
					phyDict[temp_popname] = ("N" * runnSeqLen) + seq # if it's not, add N-string if needed before sequence

			#if we're using HapB
			elif full_popname[-1] == "B" and hapChoose == "B":
				pops_this_contig.append(temp_popname)
				seq = line.split()[1]

				if temp_popname in phyDict.keys():
					phyDict[temp_popname] += seq #if it's in the dictionary, add the new sequence
				else:
					phyDict[temp_popname] = ("N" * runnSeqLen) + seq # if it's not, add N-string if needed before sequence

	#make sure we get last contig
	for key in phyDict.keys():
		if key not in pops_this_contig:
			phyDict[key] += "N" * len(seq)

	runnSeqLen += len(seq)

	fileout.write(str(len(phyDict.keys())) + " " + str(runnSeqLen) + "\n")
	for key in sorted(phyDict.keys()):
		fileout.write(key + "\t" + phyDict[key] + "\n")
		print(phyDict[key])


imaToPhylip("example_IMA2_format_fixed.fasta", "example_phylip_fixed.phy")

