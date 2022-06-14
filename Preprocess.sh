for k in $1;do cd ${k};python Rename.py;cd ../;done


for k in $1;do cd ${k};for i in {IgM,IgD,IgG,IgA,IgE};do 
if [ -f $i.txt ];then
bsub -n 16 -q max_34_cpu -e err -o out python GetConsensus_V2.py ${i}.txt 5 Consensus_${i}.txt;
fi;done;cd ../;done


for k in $1;do cd ${k};for i in {IgM,IgD,IgG,IgA,IgE};do 
if [ -f "Consensus_${i}.txt" ];then
mv Consensus_${i}.txt ${i}/Consensus.txt;
fi;
done;cd ../;done


for k in $1;do cd ${k};for i in {IgM,IgD,IgG,IgA,IgE};do 
if [ -f "${i}/Consensus.txt" ];then 
cd ${i};bsub -n 4 -q cpu -e err -o out python Csv2fastaofConseusns.py Consensus.txt Consensus.fasta;
else
 continue
fi;cd ../;done;cd ../;done


for k in $1;do cd ${k};for i in {IgM,IgD,IgG,IgA,IgE};do cd ${i};mkdir -p subfile_fasta;bsub -n 4 -q cpu -e err -o out python split_to_subfiles.py Consensus.fasta subfile_fasta;cd ../;done;cd ../;done


for k in $1;do cd ${k};for i in `echo Ig*/subfile_fasta`;do cd ${i};for j in `echo merge*.fasta`;do bsub -n 8 -q cpu -e err -o out sh IgBLAST4CMBCR_KIMDB.sh ${j};done;cd ../../;done;cd ../;done


for k in $1;do cd ${k};for i in `echo Ig*/subfile_fasta`;do cd ${i};for j in `echo merge*.fasta`;do bsub -n 2 -q cpu -e err_${j} -o out_${j} python new_ParseIgBLAST20210810.py -f ${j} -i GenoNovel-20220124-IgBLAST.${j}.m7.txt -o GenoNovel20220124-Seq${j}.txt;done;cd ../../;done;cd ../;done


for k in $1;do cd ${k};for i in `echo Ig*/subfile_fasta`;do cd ${i};bsub -n 2 -q cpu -e err -o out python ConsensusIgblast2.py GenoNovel20220124-Consensusoutput.txt;cd ../../;done;cd ../;done