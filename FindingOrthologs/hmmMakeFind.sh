#!/bin/bash
#A bash pipeline for building HMMs from many multi-fastas containing amino-acid seqs for gene families, and querying those
#HMMs against a reference amino-acid sequence set to find best matches.
#Dependancies: TrimAL, hmmbuild, hmmsearch, MUSCLE, awk, sed

echo "Usage hmmMakeFind.sh data/AA_Family.fasta TargetDB.fasta output.csv hmmScoreTreshold"
# i.e. $ hmmMakeFind.sh clusters myTransCDS.fa output.csv

rm $3

rm -rf ./TRIMEDFASTA
mkdir ./TRIMEDFASTA

rm -rf ./MUSCLE
mkdir ./MUSCLE
mkdir ./MUSCLE/TMP

rm -rf ./CREATED_HMM
mkdir ./CREATED_HMM

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
#Use this instead of last Trimal step if do not want to remove highly gapped blocks from alignment
#cp ./MUSCLE/TMP/$name'_second'.aln ./MUSCLE/$name.aln
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
hmmbuild -n $hmmname -o ./CREATED_HMM/$name.log --amino ./CREATED_HMM/$name.hmm ./MUSCLE/$name.aln
echo $name': HMM created!'
done

for i in $(find ./CREATED_HMM -name '*.hmm');
do
file=$(basename $i)
name="${file%.*}"
#Search for HMM matches in a multi-fasta list.
hmmsearch --incT $4 --noali --tblout ./FOUND_HMM/$name.tab ./CREATED_HMM/$name.hmm $2
echo $name': Queried against DB '$2
#Extract hits from search report, add to list.
awk '!/#/ {print $0}' ./FOUND_HMM/$name.tab >> TMP/$3.raw
echo $name': HMM query results written to output!'
done

#Convert hmmsearch results to CSV, sort by query name then bitscore (highest to lowest)
sed 's/[[:space:]]\{1,\}/,/gp' TMP/$3.raw > TMP/$3.cols
#Return only rows with a score (col6) above threshold from a comma delimited file. 
awk -F '[,]' -v x=$4 '$6 >= x {print $0}' TMP/$3.cols > TMP/$3.threshold
#Sort comma delimited file by col3 then col6 (as number, reverse order)
sort --field-separator=',' -k3,3 -k6,6nr TMP/$3.threshold > TMP/$3.sorted

#Reformat found HMM matches as csv.
echo ',,,,full sequence,,,best 1 domain,,,domain number estimation,,,,,,,,' > $3
echo 'target_name,accession,query_name,accession,E-value,score,bias,E-value,score,bias,exp,reg,clu,ov,env,dom,rep,inc,description of target' >> $3
awk '{print $0}' TMP/$3.sorted >> $3
echo 'All HMMSearch results saved as formatted CSV file!'

#Dump tmp dir
echo 'Cleaning up temp files.'
rm -rf TMP
echo 'Done!'
 



