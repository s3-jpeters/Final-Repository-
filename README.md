# Tutorial for evolutionary gene history analysis of STT3B in BIO312 
  ### Contents  
  1. [Introduction](#1-introduction)
  2. [Lab03](#2-lab03)
  4. [Lab04](#4-lab04) 
  5. Lab05
  6. Lab06
  7. Lab07
  8. Lab08  
# 1. Introduction
This repository contains the analysis pipeline and results for studying the evolutionary history of the STT3B gene family. It includes all commands, scripts, and output files necessary to replicate the study, along with detailed documentation. The primary focus of this analysis is to explore the inheritance and conservation of STT3B gene copies across various vertebrate lineages, using multiple computational tools and methods such as BLASTP, sequence alignment, phylogenetic tree construction, and reconciliation analysis. By following the steps and information provided, users can reproduce the study’s findings and gain insights into the evolutionary dynamics of the STT3B gene family.

Before each lab clone that labs repository 

```
git clone https://github.com/Bio312/lab(insert lab)-$MYGIT
```

Next move to this new repository 

```
cd lab(insert lab)-$MYGIT
```
This is the same first step for all labs 3-8

At anytime to check what directory you are in use the pwd command 

```
pwd
```
# 2. Lab03:  Finding homologs with BLAST KEY

## Create the BLAST database 

First go to the lab 3 directory 

```
cd ~/lab03-$MYGIT
```
Next in order to uncompress the proteomes. Run the following command:

```
gunzip proteomes/*.gz
```
In order to put all the protein sequences into a single file use the cat command

```
cat  proteomes/*.faa > allprotein.fas
```
Make the BLAST database using this command: 

```
makeblastdb -in allprotein.fas -dbtype prot
```
To organize the BLAST searches, we will first create a directory for the STT3B search and navigate into it:

We first created a folder for the STT3B BLAST search using mkdir command 

```
mkdir ~/lab03-$MYGIT/STT3B
```
Check you are in the right directory

## Download the query protein

```
ncbi-acc-download -F fasta -m protein "NP_849193.1"
```
## Perform the BLAST search 

```
blastp -db ../allprotein.fas -query NP_849193.1.fa -outfmt 0 -max_hsps 1 -out STT3B.blastp.typical.out
```
## Perform a BLAST search, and request tabular output

We will create an easier to analyse output. 
This command will allow us to better interpret our BLAST search 

```
blastp -db ../allprotein.fas -query NP_849193.1.fa  -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle"  -max_hsps 1 -out STT3B.blastp.detail.out
```
Filtering the BLAST output for high scoring homologs 

We only want high scoring matches so we used a e-value of 1e-30 for our gene STT3B

```
awk '{if ($6< 1e-30)print $1 }' Stt3B.blastp.detail.out > STT3B.blastp.detail.filtered.out
```
In order to actually count the usedul homologs in the BLAST we used this command 

```
grep -o -E "^[A-Z]\.[a-z]+" STT3B.blastp.detail.filtered.out  | sort | uniq -c
```
The output for the above command should give you a table of 11 species and the number of paralogs found in each species. 

# 3. Lab04- Gene Family Sequence Alighnment 












