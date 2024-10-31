# Tutorial for evolutionary gene history analysis of STT3B in BIO312 
  ### Contents  
  1. [Introduction](#1-introduction)
  2. [Lab3](#2-lab3)
  4. Lab04 
  5. Lab05
  6. Lab06
  7. Lab07
  8. Lab08  
# 1. Introduction
This repository contains the analysis pipeline and results for studying the evolutionary history of the STT3B gene family. It includes all commands, scripts, and output files necessary to replicate the study, along with detailed documentation. The primary focus of this analysis is to explore the inheritance and conservation of STT3B gene copies across various vertebrate lineages, using multiple computational tools and methods such as BLASTP, sequence alignment, phylogenetic tree construction, and reconciliation analysis. By following the steps and information provided, users can reproduce the study’s findings and gain insights into the evolutionary dynamics of the STT3B gene family.

Before each lab we would clone that labs repository 

```
git clone https://github.com/Bio312/lab(insert lab)-$MYGIT
```

Next we will need to move to this new repository 

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
Check you are in the write directory

## Download the query protein

```


NP_849193.1









