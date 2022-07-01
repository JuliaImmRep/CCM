import pandas as pd
import numpy as np
import csv, os, sys
from collections import Counter
from statistics import mode



def numbering_seq(listedSeq):
	mapping = {"A":0, "T":1, "C":2, "G":3, "N":4}
	numed = list(map(lambda base:mapping[base], listedSeq))
	return numed

def Cal(seq):
	base = np.array(["A", "T", "C", "G", "N"])
	baseArr = np.array([0, 1,2,3,4])
	array = np.array(list(map(lambda x:numbering_seq(x), seq)))
	Arr  = (array == baseArr[:, np.newaxis, np.newaxis]).sum(1)
	maxvalue = np.max(Arr, axis = 0)
	issingle = (Arr == maxvalue).sum(0)
	result = (issingle > 1).any()
	if result:
		consensus = np.nan
	else:
		seqArr = Arr.argmax(0)
		if (seqArr == 4).any():
			consensus = np.nan
		else:
			consensus = seqArr.choose(base)
			consensus = "".join(consensus)
	return consensus

def GetAC(group):
	Seqcount = np.array(Counter(group["Variable"]).most_common())
	SeqArr, NumArr = Seqcount[:,0], Seqcount[:,1].astype(np.int16)
	if NumArr.shape[0] == 1:
		A, AS = SeqArr[0], NumArr[0]
	else:
		if NumArr[0] == NumArr[1]:
			A, AS = np.nan, np.nan
		else:
			A, AS = SeqArr[0], NumArr[0]
	C = Cal(group["Variable"].values)
	CS = group.shape[0]
	return A, AS, C, CS
	
	

def main(infile, threfile, output):
	df = pd.read_csv(infile, sep = "\t", usecols = ["UMI5", "UMI3", "Sequence"])
	df.rename(columns = {"CDR3":"CDR3nt", "Sequence":"Variable"}, inplace = True)
	df.dropna(how = 'any', inplace = True)
	df["Target"] = df["UMI5"] + "|" + df["UMI3"]
	result = Counter(df["Target"])
	df["UMInum"] = list(map(lambda x:result[x], df["Target"]))
	del df["Target"]
	df = df[df["UMInum"] >= threfile]
	df["Length"] = df["Variable"].apply(len)
	groups = df.groupby(["UMI5", "UMI3"])
	out = csv.writer(open(output, 'w'), delimiter = "\t")
	out.writerow(["UMI5", "UMI3", "Totalsize", "Top1", "Top1size", "Consensus", "Consensussize", "Flag"])
	for (umi5, umi3), group in groups:
		TotalSize = group.shape[0]
		VRlengthinfo = np.bincount(group["Length"].values)
		Mostlength = np.max(VRlengthinfo)
		Ismost = VRlengthinfo == Mostlength
		if Ismost.sum() >= 2:
			A, AS = np.nan, np.nan
			C, CS = np.nan, np.nan
		else:
			mostlength = np.argmax(VRlengthinfo)
			usegroup = group[group["Length"] == mostlength]
			A, AS, C, CS= GetAC(usegroup)
		flag = "Y" if A == C else "N"
		out.writerow([umi5, umi3, TotalSize, A, AS, C, CS, flag])

if __name__ == "__main__":
	infile, threfile, output = sys.argv[1], int(sys.argv[2]), sys.argv[3]
	main(infile, threfile, output)