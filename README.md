# Tutorial for evolutionary gene history analysis of STT3B in BIO312 
  ### Contents  
  1. [Introduction](#1-introduction)
  2. [BLAST](#2-BLAST)
  4. [Allighnment](#4-Allighnment) 
  5. [IQ-Tree](#5-IQ-Tree) 
  6. Lab06
  7. Lab07
  8. Lab08
  9. [Conclusion](#9conclusion) 
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
# 2. BLAST

## Lab03:  Finding homologs with BLAST KEY

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

# 3. Allighnment 
## Lab04: Gene family sequence alignment

Clone lab 04 and use the cd command to move into this directory 

Make sure the VS code extention used to visualize pdf's is downloaded 

## Alighn STT3B gene family members 

create a folder for the STT3B gene family in lab 04 using this command 

```
mkdir ~/lab04-$MYGIT/STT3B
```
go into this directory using the cd command 

```
cd ~/lab04-$MYGIT/STT3B
```
Obtain the sequence that were in the BLAST output file from lab03 using this code

```
seqkit grep --pattern-file ~/lab03-$MYGIT/STT3B/STT3B.blastp.detail.filtered.out ~/lab03-$MYGIT/allprotein.fas | seqkit grep -v -p "carpio" > ~/lab04-$MYGIT/STT3B/STT3B.homologs.fas
```
## Perform a global multiple sequence alignment in muscle:

Use the following command to make a multiple sequence allighnment using muscle:

```
muscle -align ~/lab04-$MYGIT/STT3B/STT3B.homologs.fas -output ~/lab04-$MYGIT/STT3B/STT3B.homologs.al.fas
```
Now we are going to view this allighnment then convert it to a pdf so it is easier to analyze

The following is the command used to view the allighnment: 

```
alv -kli  ~/lab04-$MYGIT/STT3B/STT3B.homologs.al.fas | less -RS
```
Use the majority option for easier viewing 

```
alv -kli --majority ~/lab04-$MYGIT/globins/globins.homologs.al.fas | less -RS
```

This command shows the alignment file with a majority consensus sequence, highlighting the most common residues at each position. The output is formatted for easy viewing and scrolling in the terminal.

Next we will run an R-script code for printing the allighnment:

```
 Rscript --vanilla ~/lab04-$MYGIT/plotMSA.R  ~/lab04-$MYGIT/STT3B/STT3B.homologs.al.fas
```

## Get information about this allighnment 

First we will caluclate the width of this allighnment using this command: 

```
alignbuddy  -al  ~/lab04-$MYGIT/STT3B/STT3B.homologs.al.fas
```
Now we will calculate the allghnment removing any of the columns that have gaps in them:

```
alignbuddy -trm all  ~/lab04-$MYGIT/STT3B/STT3B.homologs.al.fas | alignbuddy  -al
```
Next we will caluclate the length of the allighnment after removing completly conserved positions:

```
alignbuddy -dinv 'ambig' ~/lab04-$MYGIT/STT3B/STT3B.homologs.al.fas | alignbuddy  -al
```
## Caluclate average percent identity 
First we are gonna calulate the average percent identity using t-coffee: 

```
t_coffee -other_pg seq_reformat -in ~/lab04-$MYGIT/STT3B/STT3B.homologs.al.fas -output sim
```
Here is how to calculate the percent identity using alighnbuddy: 

```
 alignbuddy -pi ~/lab04-$MYGIT/STT3B/STT3B.homologs.al.fas | awk ' (NR>2)  { for (i=2;i<=NF  ;i++){ sum+=$i;num++} }
     END{ print(100*sum/num) } '
```
# 4. IQ-Tree




















