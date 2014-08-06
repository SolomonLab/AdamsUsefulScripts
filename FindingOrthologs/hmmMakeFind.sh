#!/bin/bash
#A bash pipeline for building HMMs from many multi-fastas containing amino-acid seqs for gene families, and querying those
#HMMs against a reference amino-acid sequence set to find best matches.
#Dependancies: TrimAL, hmmbuild, hmmsearch, MUSCLE, awk, sed

echo "Usage hmmMakeFind.sh path/to/multi.fasta path/TargetDB.fasta output.txt"
# i.e. $ hmmMakeFind.sh clusters PsnTransCDS.fa output.csv

rm $3

rm -rf ./TRIMEDFASTA
mkdir ./TRIMEDFASTA

rm -rf ./MUSCLE
mkdir ./MUSCLE
mkdir ./MUSCLE/TMP

rm -rf ./FUNY_HMM
mkdir ./FUNY_HMM

rm -rf ./FOUND_HMM
mkdir ./FOUND_HMM

rm -rf TMP
mkdir TMP

for i in $(find $1 -name '*.fasta');
do
file=$(basename $i)
name="${file%.*}"
echo $name': Begin alignment.'
#Copy source fasta, and trim seq names to first column.
awk '{print $1}' $1/$file > TRIMEDFASTA/$name.fa
#Align new multifasta
muscle -in TRIMEDFASTA/$name.fa -out ./MUSCLE/TMP/$name'_first'.fa -quiet
echo $name': First align complete!'
# Discard spurious seqs from alignment.
#-resoverlap: Minimum overlap of a positions with other positions in the column to be considered a "good position". Range: [0 - 1]
#-seqoverlap: Minimum percentage of "good positions" that a sequence must have in order to be conserved. Range: [0 - 100]
trimal -in ./MUSCLE/TMP/$name'_first'.fa -out ./MUSCLE/TMP/$name'_noSpur'.aln -resoverlap 0.25 -seqoverlap 60 
echo $name': Spurious seq removal complete!'
#Realign after discarding spurious sequences.
muscle -in ./MUSCLE/TMP/$name'_noSpur'.aln -out ./MUSCLE/TMP/$name'_second'.aln -quiet
echo $name': Realign complete!'
#Remove columns with gaps in more than 80% of positions.
#-gt = fraction of seqs with gap allowed i.e. out of 20 seqs 4 must have residues in column.
trimal -in ./MUSCLE/TMP/$name'_second'.aln -out ./MUSCLE/$name.aln -gt 0.2
echo $name': Trimmed poorly conserved blocks from alignment.'
done

#Delete intermediate alignments.
rm -rf ./MUSCLE/TMP

for i in $(find $1 -name '*.fasta');
do
file=$(basename $i)
name="${file%.*}"
#Check how many sequences made it into the final alignment.
seqcount=$(awk '/>/ {print $0}' ./MUSCLE/$name.aln | wc -l)
echo $seqcount
#Build the HMMs from the quality checked alignments.
hmmname=$name'_'$seqcount
hmmname=$(echo $hmmname | awk '{gsub(/[[:space:]]/,"")} {print $0}')
echo $hmmname
hmmbuild -n $hmmname -o ./FUNY_HMM/$name.log --amino ./FUNY_HMM/$name.hmm ./MUSCLE/$name.aln
echo $name': HMM created!'
done

for i in $(find ./FUNY_HMM -name '*.hmm');
do
file=$(basename $i)
name="${file%.*}"
#Search for HMM matches in a multi-fasta list.
hmmsearch --noali --tblout ./FOUND_HMM/$name.tab ./FUNY_HMM/$name.hmm $2
echo $name': Queried against DB '$2
#Extract hits from search report, add to list.
awk '!/#/ {print $0}' ./FOUND_HMM/$name.tab >> TMP/$3.out
echo $name': HMM query results written to output!'
done

#Reformat found HMM matches as csv.
echo ',,,,full sequence,,,best 1 domain,,,domain number estimation,,,,,,,,' > $3
echo 'target_name,accession,query_name,accession,E-value,score,bias,E-value,score,bias,exp,reg,clu,ov,env,dom,rep,inc,description of target' >> $3
sed 's/[[:space:]]\{1,\}/,/gp' TMP/$3.out >> $3
echo 'All HMMSearch results saved as formatted CSV file!'

#Dump tmp dir
echo 'Cleaning up temp files.'
rm -rf TMP
echo 'Done!'
 



