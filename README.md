AdamsUsefulScripts
==================

A collection of handy scripts and config files for common lab tasks.

===Please Add a Description and Use Case for Each New Script===

==hmmMakeFind.sh==
A bash pipeline for building HMMs from many multi-fastas containing amino-acid seqs for gene families, and querying those
HMMs against a reference amino-acid sequence set to find best matches. Returns csv summary of all hmmsearch hits for each query hmm in the target database.

Dependancies: TrimAL, hmmbuild, hmmsearch, MUSCLE, awk, sed

Usage: hmmMakeFind.sh <path/to/multi.fasta> <path/TargetDB.fasta> <output.csv>
