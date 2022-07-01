# Installation

## Software

    In-house scripts were written in Python (v3.7.4) using modules including NumPy (v1.16.5), Biopython (v1.73), Levenshtein (v0.12.0), SciPy (v1.3.0), Collections, scikit-learn (v1.0.2) and Pandas (v0.25.1). The algorithm of pairwise2.align.globalms from Biopython (v1.73) was applied for the global alignment. Modules including seaborn (v0.10.1) and Matplotlib (v3.3.2) were used for visualization. The statistical methods used for each analysis were labeled in the figure legends accordingly.

# Data preprocess

## primers and barcode sequences

<style>
</style>

***python FindSequenceInformation.py -LibName LibID -f1 test_r1.fastq -f2 test_r2.fastq -Libconf library_conf.csv -PFile primer_conf.csv -BFile barcode_conf.csv -UMI5Len 12 -UMI3Len 12 -d output directory -o seq.result.txt***



    The script, “FindSequenceInformation.py”, implements the sequence information of UMI, barcode and primer. The detailed parameters were described as follows:

-LibName: the library information of samples, such as: “Lib6” (The details in library_conf.csv).

-f1: forward fastq file, such as: “data/test_r1.fastq” (test_r1.fastq.gz).

-f2: reversed fastq file, such as: “data/test_r2.fastq” (test_r2.fastq.gz).

-Libconf: the information of library: “data/library_conf.csv”.

-PFile: the information of primer: “data/primer_conf.csv”.

-BFile: the information of barcode: “data/barcode_conf.csv”.

-UMI5Len: the length of 5’ UMI sequences: 12

-UMI3Len: the length of 3’ UMI sequences: 12

-d: The directory of output result, such as: “./data”.

-o: The filename of the result, such as: “seq.result.txt”.

**Output:**

-seq.result.txt:
The information of sequences.
