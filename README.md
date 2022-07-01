# README

## Overview

This repository contains all scripts and preliminary statistics related to the key analyses in each section in the Result of the manuscript (MS) " Multi-tissue architecture of the adaptive immune receptor repertoire in the cynomolgus macaque". These analyses can be classified into five categories, including data preprocess, novel alleles, gene usage, clone network and SHM, and class switch, And the scripts and test data are stored in separate folders. More details can be found in the following sections.

## Software

In-house scripts were written in Python (v3.7.4) using modules including NumPy (v1.16.5), Biopython (v1.73), Levenshtein (v0.12.0), SciPy (v1.3.0), Collections,
scikit-learn (v1.0.2) and Pandas (v0.25.1). The algorithm of pairwise2.align.globalms from Biopython (v1.73) was applied for the global alignment. Modules including seaborn (v0.10.1) and Matplotlib (v3.3.2) were used for visualization. The statistical methods used for each analysis were labeled in the figure legends accordingly.

# Data preprocess

**Obtain primers and barcode sequences**

***python FindSequenceInformation.py -LibName LibID -f1 test_r1.fastq -f2 test_r2.fastq -Libconf library_conf.csv -PFile primer_conf.csv -BFile barcode_conf.csv -UMI5Len 12 -UMI3Len 12 -d output directory -o seq.result.txt***

The script, “FindSequenceInformation.py”, implements the sequence information of UMI, barcode, and primer. The detailed parameters were described as follows:

-LibName: the library information of samples, such as “Lib6” (The details in
library_conf.csv).

-f1: forward fastq file, such as: “data/test_r1.fastq” (test_r1.fastq.gz).

-f2: reversed fastq file, such as: “data/test_r2.fastq” (test_r2.fastq.gz).

-Libconf: the information of library: “data/library_conf.csv”.

-PFile: the information of primer: “data/primer_conf.csv”.

-BFile: the information of barcode: “data/barcode_conf.csv”.

-UMI5Len: the length of 5’ UMI sequences: 12

-UMI3Len: the length of 3’ UMI sequences: 12

-d: The directory of the output result, such as: “./data”.

-o: The filename of the result, such as: “seq.result.txt”.

**Output:**

-seq.result.txt:
The information of sequences.

**Assemble consensus sequences**

***python GetConsensus.py isotype-specific seq.result.txt threshold output_name***

The script, “GetConsensus.py”, was used to obtain consensus sequences. The detailed parameters were described as follows:

-isotype-specific.seq.result.txt: The information of specific-isotype sequences from each sample. Such as: “data/IgD.txt”, due to the limitation of file size, the file was compressed as “IgD.txt.gz”.

-threshold: The threshold value of UMIs’ size: 5

-output_name: The file name of the output result: “Consensus.txt”.

**Output:**

-Consensus.txt: The assembled consensus sequences by dual UMIs from sequences.

**Percentage of multi-isotype UMI groups**

***python CalUMIRatio.py seq.result.txt output_directory***

The script, “CalUMIRatio.py”, was used to calculate the percentage of multi-isotype UMI groups within each sample. The detailed parameters were described as follows:

-seq.result.txt: “data/seq.result.txt”, due to the limitation of file size, only selected 2500 lines of the original file as example data.

-output_directory: The directory of the output result, such as: “./data”.

**Output:**

-UMI.isotype.txt: The output recorded the percentage of multi-isotype UMI groups.

# Novel alleles

**Inferring novel alleles**

***python GetRealNovelAlleles.py KIMDB IMGTDB OutputFile***

The script, "GetRealNovelAlleles.py", implements the real V novel alleles of BCR from the results of the IgDiscover. The detailed parameters were described as follows:

-KIMDB: A file contains 615 V allele sequences from the KIMDB database: “data/Macaca-fascicularis_Ig_Heavy_V_615.fasta”.

-IMGTDB: A file contains V allele sequences from the IMGT database: “data/fascicularis.imgt.ighv.dedup.fasta”.

-OutputFile: A file contains the real V novel alleles: “data/Real-novel-V-BCR.txt”.

It takes some files that contain the IgDiscover results of V alleles in different donors, such as “data/IgDiscover-CCM1_V.fasta”. It also takes some input files that are the genotyping results in different donors, such as “data/Database-CCM1-V.txt”.

It also outputs one file that contains the merge sequence data of the KIMDB and IMGT (“data/IMGT_615_DF.txt”).

**Distribution of positional SNPs and SAPs and calculation of normalized SNP(s) or SAP(s)**

***python TCRFindClosestAllele.py output_directory germline_TCR_file novel_TCR_file***

The script, “TCRFindClosestAllele.py”, was used to find the closest allele from the germline database for each novel allele. The detailed parameters were described as follows:

-output_directory:
The directory of the output result, such as: “./data”.

-germline_TCR_file:
The information of TRBV germline sequences: “data/complete.TRBV.fasta”

-novel_TCR_file:
The information of novel TRBV sequences: “data/Novel_TRBV_seq.txt”

**Output:**

-FigS2B.TCR.closest.txt:
The result of the closest allele for each novel allele.

***python TCRCalSNPAndSAP.py output_directory novel_TCR_file nogap_file gap_file
closest_file***

The script, “TCRCalSNPAndSAP.py”, was used to calculate SNP and SAP between the closest allele and novel allele. The detailed parameters were described as follows:

-output_directory: The directory of the output result, such as: “./data”.

-novel_TCR_file: The information of novel TRBV sequences: “data/Novel_TRBV_seq.txt”.

-nogap_file: TCR germline sequences without IMGT-gap: “data/TRBV_without_gap.fasta”.

-gap_file: TCR germline sequences with IMGT-gap: “data/TRBV_with_gap.fasta”.

-closest_file: FigS2B.TCR.closest.txt: “data/FigS2B.TCR.closest.txt”.

**Output:**

-FigS2B.TCR.SNPs.txt: The SNP numbers from FR1 to FR3 for each novel allele.

-Summary.TCR.SNP.SAP.result.txt: The summary of detailed SNP and SAP of novel alleles.

-FigS3C.TCR.PlotSNP.txt: The distribution of SNPs.

***python TCRSNPsAndSAPsNormalization.py output_directory summay_file***

The script, “TCRSNPsAndSAPsNormalization.py”, was used to calculate the normalized SNPs and SAPs based on the length of the region. The detailed parameters were described as follows:

-output_directory: The directory of the output result, such as: “./data”.

-summay_file: “data/Summary.TCR.SNP.SAP.result.txt”.

**Output:**

-FigS2E.SNP-SAP.TCR.normalize.length.txt: The result of the normalized number of SNPs and SAPs.

**Genomic Validation of IGHV novel alleles**

***python NovelValidation.py OutputFile***

The script, NovelValidation.py, implements the validation of V novel alleles of BCR from
the real V novel alleles above. The detailed parameters were described as follows:

-OutputFile: Name of the result file.

It takes four input files including the real V novel alleles (data/Real-novel-V-BCR.txt), the database reference sequences (data/IMGT_615_DF.txt), the genotyping results (such as data/Database-CCM1-V.txt), and the Sanger sequencing results (such as data/CM1-Sanger-Cross-Main.fasta).

It outputs a file (data/ValidationAlleles.txt) containing the message of alleles ID and relevant sequences. The column named “Validation_in_donor” is the validation results of the alleles and NA is representing these alleles that cannot be validated in our data.

# Gene usage

## Random selection of clones in IGH and TRB samples

***python FetchCloneAndGeneUsage.py isotype threshold output_directory clone_files***

The script, FetchCloneAndGeneUsage.py, implements the V gene usage in subsample and sample. The detailed parameters were described as follows:

-isotype: The isotype of the file (BCR: IgM, IgA, IgG, TCR: TRB), such as “IgM”.

-threshold: The number of clones from the sample, such as 100.

-output_directory: The directory of the output result, such as: “./data”.

-clone_files: The clone files of samples. It could input more than 1 clone file of the sample. Such as: “data/Clone/CCM1-*-IgM.ChangeO.Clone.txt”. Due to the limitation of
file size, only 1000 records from every 3 original samples were provided.

**Output:**

-IgM.Random.100.origin-clone.txt: The original files of randomly selected clones at the input number.

-IgM.Random.100.origin-geneusage.txt: The V gene usages of each subsample.

-IgM.All.geneusage.txt: The V gene usage of the sample.

**Clustering samples by core V gene usage**

***python GeneClusterBCR.py output_directory gene_usage_file gene_name***

The script, “GeneClusterBCR.py”, was used to cluster samples by V gene usage. The detailed parameters were described as follows:

-output_directory: The directory of the output result, such as: “./data”.

-gene_usage_file: The core IGHV gene usages from samples: “data/core.BCR.V.geneusage.txt”

-gene_name: The name of subjected gene: “V”

**Output:**

-Fig3A.Dis.V.txt: The result of Euclidean distance between samples.

-Fig3A.Dis.V.heatmap.pdf: The output figures of Fig 3A.

**Detection of differentially used core V gene among CMs**

***python GeneUsageOneWayANOVA.py*** ***output_directory gene_usage_file***

The script, “GeneUsageOneWayANOVA.py”, was used to detect differentially used core V gene. The detailed parameters were described as follows:

-output_directory: The directory of the output result, such as: “./data”.

-gene_usage_file: The core IGHV gene usages from samples: “data/core.BCR.V.geneusage.txt”

**Output:**

-Fig3C.BCR.V.ANOVA.table.txt: The result of differentially used core V gene.

-Fig3C.BCR.ANOVA.V.Gene_differ.heatmap.pdf: The output figures of Fig 3C.

**Principal component analysis of gene usage**

***python BCR.PCA.GeneUsage.py*** ***output_directory gene_usage_file***

The script, “GeneUsageOneWayANOVA.py”, was used to do the PCA analysis by core V gene usage. The detailed parameters were described as follows:

-output_directory: The directory of the output result, such as: “./data”.

-gene_usage_file: The core IGHV gene usages from samples: “data/core.BCR.V.geneusage.txt”

**Output:**

The output figures of Figs S3B-D:

-FigS3B-D.PCA.BCR.Isotype.V.usage.png

-FigS3B-D.PCA.BCR.Tissue.V.usage.png

-FigS3B-D.PCA.BCR.Tissue_Group.V.usage.png

**Calculation of gene family usage**

***python GeneFamilyUsage.py*** ***output_directory gene_usage_file***

The script, “GeneFamilyUsage.py”, was used to compare the IGHV family usage among different groups. The detailed parameters were described as follows:

-output_directory: The directory of the output result, such as: “./data”.

-gene_usage_file: The core IGHV gene usages from samples: “data/core.BCR.V.geneusage.txt”

**Output:**

-Fig3E.Isotype.t-test.txt: The t-test result of Fig 3E.

-Fig3F.Tissue-group.t-test.txt: The t-test result of Fig 3F.

-Fig3E.Family.usage.V.isotype.png: The output figures of Fig 3E.

-Fig3F.Family.usage.V.tissue-group.png: The output figures of Fig 3F.

# Clone network and SHM

## Overlapping CDR3s between CM and RAPID (Fig. 4B)

***python OverlapCDR3sCMRAPID.py RAPIDAnnoFile AllOrPublic***

The script, OverlapCDR3sCMRAPID.py, implements the overlapping CDR3s between CM and human. The detailed parameters were described as follows:

-RAPIDAnnoFile: The file contains the annotation data from RAPID and the example file only contains the top 10000 lines of the original file because the original file is too large: “example-Abs.Dedup.txt”

-AllOrPublic: The datatype of CM’s CDR3s needs to be analyzed. If it is set as “all (or public)”, the script will compare the overlapping between CM’s all (or public) CDR3s and RAPID’s CDR3s.

In addition, it needs to put the clone files of all samples in the working directory, and we will supply partial clone files for the test because the number of all clone files is too much. The example clone files are in this directory (BCRExampleClone/BCR-CCM*-BM-IgM-ChangeO-GenoNovel.txt), due to the limitation of file size, the file named “BCRExampleClone/BCR-CCM2-BM-IgM-ChangeO-GenoNovel.txt” was compressed as “BCRExampleClone/BCR-CCM2-BM-IgM-ChangeO-GenoNovel.txt.gz”.

It outputs four files including the meta CDR3s (example-meta-aa.txt), and the number of overlapping CDR3s (such as example_Anno_venn_data_in_all.txt), and the message of annotation sequences in Known (example_Known_CDR3_aa_in_all.txt) and Therapeutic Abs (example_Thera_CDR3_aa_in_all.txt).

## Length distribution of CDR3s (Fig. 4D)

***python CM-RAPID-len.py MetaAAFile RAPIDAnnoFile OutputFile***

The script, CM-RAPID-len.py, implements the length distribution of CM and RAPID.

The detailed parameters were described as follows:

-MetaAAFile: File generated by the script named “OverlapCDR3sCMRAPID.py” above, which contains the meta CDR3s data of CM: “example-meta-aa.txt”.

-RAPIDAnnoFile: The file contains the annotation data from RAPID and the example file only contains the top 10000 lines of the original file because the original file is too large: “example-Abs.Dedup.txt”

-OutputFile: File contains the length distribution of CM and RAPID dataset: “example-len-distribution-CM-RAPID.txt”.

***python Human-len.py MetaAAFile HuamanAllClone HuamanPublicClone***

The script, Human-len.py, implements the length distribution of human. The detailed parameters were described as follows:

-MetaAAFile: File generated by the script named “OverlapCDR3sCMRAPID.py” above, which contains the meta CDR3s data of CM: “example-meta-aa.txt”.

-HuamanAllClone: File contains the CDR3s of all clones of a human dataset and the example file only contains the top 10000 lines of the original file because the original file is too large: “example-Totaltranslate.txt”.

-HuamanPublicClone: File contains the CDR3s of public clones of a human dataset and the example file only contains the top 10000 lines of the original file because the original file is too large: “example-human-public.txt”.

The output files include length distribution of human named “example-len-distribution-human.txt” and number of overlapping CDR3s between CM and human named “example-CDR3aa-overlap-CM-human.txt”.

## Length distribution of clones in isotype groups (Fig. S4A,B)

***python LenDistributionIso.py***

The script, LenDistributionIso.py, implements the length distribution of clones of different isotype groups. It needs to put the clone files of all samples in the working directory, and we will supply partial clone files for the test because the number of all clone files is too much. The example clone files are in this directory (BCRExampleClone/BCR-CCM*-BM-IgM-ChangeO-GenoNovel.txt).

The output files include the meta clone data of isotype groups (example-Meta-clonedata-isotype.txt) and the length distribution of different isotypes in all and public clones (example-Len-isotype-All.txt and example-Len-isotype-Public.txt).

## Tissue group composition of inter-tissue clones (Fig. 5A)

***python TissueGroupComposition.py BT***

The script, TissueGroupComposition.py, implements the tissue group composition of
inter-tissue clones. The detailed parameters were described as follows:

-BT: The input file is BCR data (B) or TCR data (T).

It needs to put the two files (“example-BCR.AA.0.06.Overlapped.cloneNum.txt” and “example-TCR.AA.0.Overlapped.cloneNum.txt”) in the working directory, which contain partial clone data (200000) for the test as the original file is too large.

It contains two output files including the shared clone data (“example-BCR-Fig5A-sharedData.txt” or “example-TCR-Fig5A-sharedData.txt”) and the composition data (“example-BCR-Fig5A-percentage.txt” or “example-TCR-Fig5A-percentage.txt”).

## Similarity score in the network diagram (Fig. 5C-E, S5)

***python SimilarityNetworkAll.py BT***

The script, SimilarityNetworkAll.py, implements the similarity score between different
tissues (Fig. 5C-E). The detailed parameters were described as follows:

-BT: The input file is BCR data (B) or TCR data (T).

It needs to put the two files (“example-BCR.AA.0.06.Overlapped.cloneNum.txt” and “example-TCR.AA.0.Overlapped.cloneNum.txt”) in the working directory, which contain partial clone data (200000) for the test as the original file is too large.

The output file named “example-Network-figure-BCR.txt” (or “example-Network-figure-TCR.txt”) contains the clone number and similarity scores of different tissues of BCR (or TCR) data.

***python SimilarityNetworkDonor.py BT CM***

The script, SimilarityNetworkDonor.py, implements the similarity score between different tissues in different donors (Fig. S5). The detailed parameters were described as follows:

-BT: The input file is BCR data (B) or TCR data (T).

-CM: Donor (1, 2 or 3) It needs to put the two files (“example-BCR.AA.0.06.Overlapped.cloneNum.txt” and “example TCR.AA.0.Overlapped.cloneNum.txt”) in the working directory, which contain partial clone data (200000) for the test as the original file is too large.

The output files contain the clone number and similarity scores of different tissues in different donors, such as “example-Network-figure-BCR-CCM1.txt”.

## SHM rates of clones (Fig. 6A-B, 6H, S6A-C)

***python Fig6SHMMeta.py step datatype***

The script, Fig6SHMMeta.py, implements the SHM comparison between different groups. The detailed parameters were described as follows:

-step: Number of step (1, 2, 3, 4), different steps represent different analyses.

-datatype: All clones (all) or public clones (public).

Step 1: *python Fig6SHMMeta.py 1 all*

When the first parameter (step) was set to 1, this script will be run for getting an SHM data file, which is the input file of the scripts for the next SHM analyses. It needs to put the clone files of all samples in the working directory, and we will supply partial clone files for the test because the number of all clone files is too much. The example clone
files are in this directory (BCRExampleClone/BCR-CCM*-BM-IgM-ChangeO-GenoNovel.txt). The output file contains the SHM meta data (Partial-BCR-SHM-meta-data.txt), which is the input file for the next analyses. But because the example output file contains SHM data of a little of clones, we got SHM data of partial clone data (300000) from the original SHM meta data file for the next analyses, which named “example-BCR-SHM-meta-data.txt”.

Step 2: *python Fig6SHMMeta.py 2 all*

When the first parameter (step) was set to 2, this script will calculate the SHM of clones in different isotypes between different tissues (Fig. S6C). This command needs the input file (example-BCR-SHM-meta-data.txt) to be put in the working directory.

The output file named “example-BCR-SHM-tissue-isotype-allclone.txt” (or “example-BCR-SHM-tissue-isotype-publicclone.txt”), which contains the message of cloneid, SHM, tissue, and Isotype.

Step 3: *python Fig6SHMMeta.py 3 all*

When the first parameter (step) was set to 3, this script will calculate the SHM of all clones (all) or public clones (public) in different isotypes (Fig. 6A-B). This command needs the input file (example-BCR-SHM-meta-data.txt) to be put in the working directory.

The output file is named “example-BCR-SHM-only-isotype-allclone.txt” (or “example-BCR-SHM-only-isotype-publicclone.txt”). The T-test results can be obtained from the script named “Fig6ACohenEffectsize.py”.

Step 4: *python Fig6SHMMeta.py 4 all*

When the first parameter (step) was set to 4, this script will calculate the SHM of all clones (all) or public clones (public) of different isotypes in different tissue groups (Fig. 6C-D, S6C-D). This command needs the input file (example-BCR-SHM-meta-data.txt) to be put in the working directory.

The output files named “example-BCR-SHM-tissuegroup-isotype-allclone.txt” or “example-BCR-SHM-tissuegroup-isotype-publicclone.txt”.

***python Fig6SHMTissueShared.py***

The script, Fig6SHMTissueShared.py, implements the SHM of tissue-specific and inter-tissues shared clones of different tissues (Fig. 6H). This command needs the input file (example-BCR-SHM-meta-data-6H.txt) to be put in the working directory.

The output file (example-SHM-tissue-shared.txt) contains the message of the clone id, SHM, tissue, and if the clone is the tissue-specific clone.

# Class switch

## Class switch analysis

***python ClassSwitch.py output_directory Clone.summary.txt***

The script, “ClassSwitch.py”, was used to calculate the class switch of BCR.

-output_directory: The directory of the output result, such as: “./data”.

-Clone.summary.txt: The input file is the statistical result of overlapped clones among different isotype-specific BCR samples: “data/BCR.AA.0.06.Overlapped.cloneNum.txt”, due to the limitation of file size, the file was compressed as: “BCR.AA.0.06.Overlapped.cloneNum.txt.gz”.

**Output:**

-Fig6E.result.txt: The percentage of class switch in public clones and intra-individual clones.

The output figures of Figs 6F-G and Fig S6D:

-FigS6D.All-IgMD.Individual.Class-Switch.png

-Fig6F.All-IgMD.Public.Class-Switch.png

-Fig6G.GasIntes-IgMD.Public.Class-Switch.png
