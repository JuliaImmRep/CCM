import pandas as pd
import numpy as np
import csv, os, sys
import re
from Bio import SeqIO
import Levenshtein
import glob, argparse
from itertools import chain

def TreatN(bardic):
	adddic = {}
	for i in bardic:
		for j in range(4):
			newid = list(i)
			newid[j] = "N"
			addid = "".join(newid)
		adddic.setdefault(addid, bardic[i])
	finaldic = dict(**bardic, **adddic)
	return finaldic

def Readbarcode(barConf):
	df = pd.read_csv(barConf, sep = ",")
	df["Ori"] = df["BarName"].str.split('_', expand = True)[0]
	bar5df = df[df['Ori'] == 'B5']
	bar3df = df[df['Ori'] == 'B3']
	bar5dic = dict(zip(bar5df["BarSeq"], bar5df["BarName"]))
	bar3dic = dict(zip(bar3df["BarSeq"], bar3df["BarName"]))
	return bar5dic, bar3dic

def Splitprimer(df):
	final = {}
	groups = df.groupby("PrimersetID")
	for primerset, group in groups:
		result = dict(zip(group['PrimerSeq'], group['PrimerID']))
		final.setdefault(primerset, result)
	return final

def ReadPrimer(PrimConf, primer5setID, primer3setID):
	df = pd.read_csv(PrimConf, sep = ",")
	df.dropna(how = 'any', inplace = True)
	Primer5df = df[df["PrimersetID"].isin(primer5setID)]
	Primer3df = df[df["PrimersetID"].isin(primer3setID)]
	Primer5dic, Primer3dic = Splitprimer(Primer5df), Splitprimer(Primer3df)
	return Primer5dic, Primer3dic

def Readconf(libraryConf, library):
	df = pd.read_csv(libraryConf, sep = ",")
	df.dropna(how = 'any', inplace = True)
	usedf = df[df["Library"] == library]
	usedf["BarcodePair"] = usedf["Bar5"] + "-" + usedf["Bar3"]
	sampledic = dict(zip(usedf["BarcodePair"], usedf["SampleName"]))
	Use5setdic = dict(zip(usedf["SampleName"], usedf["Primer5"]))
	Use3setdic = dict(zip(usedf["SampleName"], usedf["Primer3"]))
	Bar5Lendic = dict(zip(usedf["SampleName"], usedf["Barcode5len"].astype(int)))
	Bar3Lendic = dict(zip(usedf["SampleName"], usedf["Barcode3len"].astype(int)))
	return sampledic, Use5setdic, Use3setdic, Bar5Lendic, Bar3Lendic

def find_barcode(barcode_dict, sequence, start, end):
	label, Tbar, mismatch, best_pos = "-", "NNNNNNNNNN", 100, "-"
	for barcode in barcode_dict.keys():
		barcode_len = len(barcode)
		candidate_bar, candidate_mis, candidate_pos = "NNNNNNNNNN", 100, "-"
		for i in range(start, end):
			misnum = Levenshtein.distance(barcode, sequence[i:i+barcode_len])
			if misnum < candidate_mis:
				candidate_bar, candidate_mis, candidate_pos = sequence[i: i + barcode_len], misnum, i
			else:
				continue
			if candidate_mis <= 5:
				ide1, seq1, pos1 = candidate_mis, candidate_bar, candidate_pos
			else:
				ide1, seq1, pos1 = 100, "NNNNNNNNNN", "-"
			if ide1 < mismatch:
				label, Tbar, mismatch, best_pos = barcode_dict[barcode], seq1, ide1, pos1
			else:
				continue
	return [mismatch, label, Tbar, best_pos]

def RemoveN(result):
	N_num = result[2].count("N")
	if N_num >= 0:
		result[0] = result[0] - N_num
	return result

def main():
	'''
	In library-free data, the oriention of the sequences is always the same, only 3'-end is necessary to disginguish isotype
	'''
	Sampleinfo, primer5s, primer3s, Barcodelen5, Barcodelen3 = Readconf(libraryConf, LibName)
	Use5s, Use3s = set(list(primer5s.values())), set(list(primer3s.values()))
	pri5dic, pri3dic = ReadPrimer(Primerfile, Use5s, Use3s)
	bar5dic, bar3dic = Readbarcode(barfile)
	print(Sampleinfo)
	fq1, fq2 = SeqIO.parse(fastq1, 'fastq'), SeqIO.parse(fastq2, 'fastq')
	os.system("mkdir -p %s"%outdir)
	out = csv.writer(open("%s/%s"%(outdir, outfile), 'w'), delimiter = "\t", lineterminator="\n")
	out.writerow(["SeqId", "Sample_id", "Out5ID", "Out3ID", "Prim5Mis", "Prim5ID", "Prim3Mis", "Prim3ID", "UMI5", "UMI3", "CDR3"])
	try:
		while True:
			rec1, rec2 = next(fq1), next(fq2)
			R1seq, R2seq = str(rec1.seq), str(rec2.seq)
			out5 = R1seq[:4]
			out3 = R2seq[:4]
			out5id, out3id = bar5dic.get(out5, "No"), bar3dic.get(out3, "No")
			foundpair = "-".join([out5id, out3id])
			sampleid = Sampleinfo.get(foundpair, "No")
			P5set, P3set = primer5s.get(sampleid, False), primer3s.get(sampleid, False)
			barcodelen5 = Barcodelen5.get(sampleid, 0)
			barcodelen3 = Barcodelen3.get(sampleid, 0)
			target5len = barcodelen5 + UMI5len
			target3len = barcodelen3 + UMI3len
			if P3set:
				P3info = find_barcode(pri3dic[P3set], R2seq, target3len - 5, target3len + 5)
			else:
				P3info = [100, '-', 'NNNN', '-1']
			if P5set:
				P5info = find_barcode(pri5dic[P5set], R1seq, target5len - 5, target5len + 5)
			else:
				P5info = [100, '-', 'NNNN', '-1']
			if P5info[0] <= 2:
				UMI5 = R1seq[barcodelen5:P5info[3]]
			else:
				UMI5 = "NNNN"
			if P3info[0] == 0:
				UMI3 = R2seq[barcodelen3:P3info[3]]
			else:
				UMI3 = "NNNN"
			if P3info[0] != 100:
				inferedCDR3 = R2seq[P3info[3] + len(P3info[2]) + 50 : 150]
			else:
				inferedCDR3 = 'NNNN'
			out.writerow([rec1.id] + [sampleid] + [out5, out3] + P5info[:2] + P3info[:2] + [UMI5, UMI3] + [inferedCDR3])
	except StopIteration:
		print("Done")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog='python find_UMI_primer_barcode.py')
	parser.add_argument('-LibName', '--LibraryName')
	parser.add_argument("-f1", '--fastq1', help = "Input fastq1")
	parser.add_argument("-f2", '--fastq2', help = "Input fastq2")
	parser.add_argument('-Libconf', '--LibraryConfFile')
	parser.add_argument('-PFile', '--Primerfile')
	parser.add_argument('-BFile', '--Barcodefile')
	parser.add_argument('-UMI5Len', '--UMI5length', type = int)
	parser.add_argument('-UMI3Len', '--UMI3length', type = int)
	parser.add_argument('-d', '--outdir')
	parser.add_argument('-o', '--output')
	args = parser.parse_args()
	fastq1, fastq2 = args.fastq1, args.fastq2
	LibName = args.LibraryName
	libraryConf, Primerfile, barfile = args.LibraryConfFile, args.Primerfile, args.Barcodefile
	UMI5len, UMI3len = args.UMI5length, args.UMI3length
	outdir, outfile = args.outdir, args.output
	main()
