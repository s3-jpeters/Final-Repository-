# Tutorial for Evolutionary Gene History Analysis of STT3B in BIO312 
  ### Contents  
  1. [Introduction](#1-introduction)
  2. [BLAST](#2-BLAST)
  3. [Alignment](#3-Alignment) 
  4. [IQ-Tree](#4-IQ-Tree) 
  5. [Reconciling](#5-Reconciling)
  6. [Protein-Domain](#6-Protein-Domain)
  7. [Results](#7-Results)
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

The following steps provide a structured approach for preparing and organizing proteomic data for a BLAST search focused on the STT3B gene family. We will navigate to the lab directory, decompress the protein files, consolidate them into a single file, create a BLAST database, and set up a dedicated directory for the STT3B analysis.

```
# Step 1: Navigate to the Lab 3 Directory
# This command changes the working directory to the Lab 3 folder.
cd ~/lab03-$MYGIT

# Step 2: Decompress the Proteome Files
# The gunzip command is used to uncompress all the .gz proteome files within the 'proteomes' folder.
gunzip proteomes/*.gz

# Step 3: Consolidate Protein Sequences into a Single File
# The cat command merges all protein sequence files (.faa) into one file named 'allprotein.fas'.
cat proteomes/*.faa > allprotein.fas

# Step 4: Create a BLAST Database
# This command generates a BLAST database from the consolidated protein sequences file.
makeblastdb -in allprotein.fas -dbtype prot

# Step 5: Set Up a Directory for STT3B BLAST Search
# Here, we create a new directory for organizing the STT3B BLAST search results.
mkdir ~/lab03-$MYGIT/STT3B

# Step 6: Confirm You Are in the Correct Directory
# Ensure that the current working directory is set for the STT3B analysis.
pwd
```
This sequence of commands sets up the environment for analyzing the STT3B gene family using BLAST. We start by navigating to the relevant directory and decompressing proteomic data. The individual protein sequences are consolidated, and a BLAST database is created to enable efficient searching. Finally, a dedicated directory is prepared to organize the subsequent analysis of the STT3B protein family.

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

# 3. Alignment 
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

The following commands will create a multiple sequence alignment (MSA) for the STT3B sequences using MUSCLE, view the alignment in the terminal, and generate a PDF for easier analysis.

```
# Step 1: Create a Multiple Sequence Alignment (MSA) using MUSCLE
# This command aligns the unaligned STT3B protein sequences to identify conserved regions.
muscle -align ~/lab04-$MYGIT/STT3B/STT3B.homologs.fas -output ~/lab04-$MYGIT/STT3B/STT3B.homologs.al.fas

# Step 2: View the MSA in the terminal
# This command displays the alignment file in the terminal using 'alv' for easy scrolling.
alv -kli ~/lab04-$MYGIT/STT3B/STT3B.homologs.al.fas | less -RS

# Step 3: View the MSA with Majority Consensus Sequence
# This command highlights the most common residues at each position for easier interpretation.
alv -kli --majority ~/lab04-$MYGIT/STT3B/STT3B.homologs.al.fas | less -RS

# Step 4: Generate a PDF of the MSA using an R script
# This command uses an R script to create a visual PDF of the alignment for detailed analysis.
Rscript --vanilla ~/lab04-$MYGIT/plotMSA.R ~/lab04-$MYGIT/STT3B/STT3B.homologs.al.fas
```
The MUSCLE command aligns the STT3B sequences to identify conserved regions. The alv tool allows us to view the alignment in the terminal, with the --majority option showing the most common residues for easier interpretation. Finally, the R script generates a PDF of the alignment for more detailed analysis.

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

The following steps outline the procedure for calculating the average percent identity of the aligned STT3B protein sequences. We will first utilize T-Coffee to determine sequence similarity, followed by AlignBuddy to compute the percent identity. These commands provide a quantitative measure of the similarity between the aligned sequences, offering insights into the conservation of the STT3B gene family.

The following steps outline how to calculate the average percent identity using T-Coffee and AlignBuddy with explanations:

```
# Step 1: Calculate Sequence Similarity using T-Coffee
# This command uses T-Coffee to reformat the aligned sequences and compute similarity scores.
t_coffee -other_pg seq_reformat -in ~/lab04-$MYGIT/STT3B/STT3B.homologs.al.fas -output sim

# Step 2: Calculate Average Percent Identity using AlignBuddy
# This command uses AlignBuddy to compute the percent identity of the aligned sequences.
# It then uses an awk script to calculate the average percent identity across all pairwise comparisons.
alignbuddy -pi ~/lab04-$MYGIT/STT3B/STT3B.homologs.al.fas | awk '(NR>2) { for (i=2; i<=NF; i++) { sum+=$i; num++ } } END { print(100*sum/num) }'
```
This workflow calculates the average percent identity of the aligned STT3B sequences, providing a metric to evaluate sequence conservation. T-Coffee generates initial similarity scores, while AlignBuddy and the accompanying awk script determine the overall percent identity. These results help assess the level of evolutionary conservation within the STT3B gene family.

# 4. IQ-Tree
## Lab 5: Gene Family Phylogeny using IQ-TREE
### Constructing a Phylogenetic Tree for STT3B Homologs 

Use the mkdir command to create a directory in lab 05 for STT3B go to that new directory using the cd command

The commands for lab05 work off the alighnment we created in lab 04 

This command standardizes labels and removes duplicate sequences from the STT3B alignment, saving the cleaned file to your lab05 directory:

```
sed 's/ /_/g'  ~/lab04-$MYGIT/STT3B/STT3B.homologs.al.fas | seqkit grep -v -r -p "dupelabel" >  ~/lab05-$MYGIT/STT3B/STT3B.homologsf.al.fas
```
This command processes the STT3B homologs file by first replacing any spaces with underscores in the sequence labels and then removing any sequences with the tag "dupelabel". The cleaned output is saved in your lab05 directory as STT3B.homologsf.al.fas.

The following command runs IQ-TREE to estimate the maximum likelihood tree for the STT3B alignment, including bootstrap support:

```
iqtree -s ~/lab05-$MYGIT/STT3B/STT3B.homologsf.al.fas -bb 1000 -nt 2
```
The input file is specified with -s. IQ-TREE determines the best substitution model, builds the tree, and estimates branch lengths. The -bb 1000 flag performs 1,000 bootstrap replicates for branch support, while -nt 2 uses 2 CPU threads to speed up the process.

## The Unrooted Tree

These commands display the phylogenetic tree from the STT3B alignment, first using a basic text graphic and then with a cleaner, unrooted graphical representation:

```
# Display the STT3B phylogenetic tree as an ASCII (text-based) graphic using Newick formatted tree file
nw_display ~/lab05-$MYGIT/STT3B/STT3B.homologsf.al.fas.treefile

# Use R script to generate an unrooted graphical version of the STT3B phylogenetic tree
# --vanilla: Runs R without loading or saving the workspace (clean session)
# Input: Newick formatted tree file (STT3B.homologsf.al.fas.treefile)
# Output: PDF file of the tree (STT3B.homologsf.al.fas.treefile.pdf)
# 0.4: Sets the text label size
# 15: Limits label length to fit the visualization
Rscript --vanilla ~/lab05-$MYGIT/plotUnrooted.R ~/lab05-$MYGIT/STT3B/STT3B.homologsf.al.fas.treefile ~/lab05-$MYGIT/STT3B/STT3B.homologsf.al.fas.treefile.pdf 0.4 15
```
The first command shows a quick text-based tree with nw_display, but it may look rooted incorrectly. The second command uses an R script to create a clearer, unrooted graphical tree, adjusting label size and length for readability, and saves it as a PDF.

## Root the optimal phylogeny

### These commands root the phylogenetic tree using the midpoint method, display it graphically, and then provide an alternative cladogram view for easier interpretation.

```
# Midpoint rooting the tree
gotree reroot midpoint -i ~/lab05-$MYGIT/STT3B/STT3B.homologsf.al.fas.treefile -o ~/lab05-$MYGIT/STT3B/STT3B.homologsf.al.mid.treefile

# Display the rooted tree in a basic text format
nw_order -c n ~/lab05-$MYGIT/STT3B/STT3B.homologsf.al.mid.treefile | nw_display -

# Create a graphical SVG image of the rooted tree
nw_order -c n ~/lab05-$MYGIT/STT3B/STT3B.homologsf.al.mid.treefile | nw_display -w 1000 -b 'opacity:0' -s > ~/lab05-$MYGIT/STT3B/STT3B.homologsf.al.mid.treefile.svg

# Convert the SVG image to a PDF
convert ~/lab05-$MYGIT/STT3B/STT3B.homologsf.al.mid.treefile.svg ~/lab05-$MYGIT/STT3B/STT3B.homologsf.al.mid.treefile.pdf

# Create a cladogram view of the tree
nw_order -c n ~/lab05-$MYGIT/STT3B/STT3B.homologsf.al.mid.treefile | nw_topology - | nw_display -s -w 1000 > ~/lab05-$MYGIT/STT3B/STT3B.homologsf.al.midCl.treefile.svg

# Convert the cladogram SVG to a PDF
convert ~/lab05-$MYGIT/STT3B/STT3B.homologsf.al.midCl.treefile.svg ~/lab05-$MYGIT/STT3B/STT3B.homologsf.al.midCl.treefile.pdf
```
  The first command roots the tree at the midpoint of the longest branch using gotree, producing a more meaningful view of divergence events. The rooted tree is initially displayed in text format using nw_display, but this is replaced with a cleaner SVG graphic, which is then converted to a PDF for easier viewing. The cladogram version of the tree, generated by nw_topology, shows all branches with equal lengths, making it simpler to see the relationships between clades without being influenced by branch lengths. The output files (SVG and PDF) are ready for analysis and comparison.

  For the STT3B gene family, we first visualized the tree and rooted it using midpoint rooting, estimating the oldest divergence based on the longest branch. Then, we examined both a phylogram (branch lengths show substitution rates) and a cladogram (all branches are equal), which made the relationships clearer. Next, we refined the tree by rooting it with an outgroup, using a known distant relative for better accuracy. Finally, we checked bootstrap support values, with anything above 80% indicating strong confidence in the branching pattern. This approach helps provide a clear and reliable view of the evolutionary history of the STT3B gene family.

# 5. Reconciling
## Lab 6. Reconciling a Gene and Species Tree

The commands set up a Python 2.7 environment with ETE3 for phylogenetic analysis:

```
mamba create -n my_python27_env python=2.7
conda activate my_python27_env
mamba install -y ete3
```
The first command creates a new environment named my_python27_env with Python 2.7 using Mamba, a faster alternative to Conda. Then, we activate this environment with conda activate. Finally, ETE3, a toolkit for tree analysis and visualization, is installed. This setup only needs to be done once and prepares the environment for future phylogenetic work.

## Preparation and Setup for STT3B Gene Family Reconciliation

In this lab, we are setting up the analysis for the STT3B gene family by copying over the midpoint-rooted gene tree from Lab 5 and organizing our files in the Lab 6 directory. This step ensures we are using the correct gene tree for subsequent reconciliation and analysis tasks.

The following code consolidates the entire process, including directory creation, file copying, pruning (if needed), and visualization:

```
# Step 1: Create a new directory for the STT3B gene family in Lab 6
mkdir ~/lab06-$MYGIT/STT3B

# Step 2: Copy the midpoint-rooted gene tree from Lab 5 to the new Lab 6 directory
cp ~/lab05-$MYGIT/STT3B/STT3B.homologs.al.mid.treefile ~/lab06-$MYGIT/STT3B/STT3B.homologs.al.mid.treefile

# Step 3: Verify that the gene tree file has been copied correctly
ls ~/lab06-$MYGIT/STT3B/STT3B.homologs.al.mid.treefile

# Step 4: (Optional) Prune unwanted sequences from the gene tree
gotree prune -i ~/lab06-$MYGIT/STT3B/STT3B.homologs.al.mid.treefile \
-o ~/lab06-$MYGIT/STT3B/STT3B.homologs.pruned.treefile \
<sequence1> <sequence2> <sequence3>

# Step 5: Display the pruned gene tree for verification (if pruning was performed)
nw_display ~/lab06-$MYGIT/STT3B/STT3B.homologs.pruned.treefile
```
This setup process for the STT3B gene family involves creating a new directory in the Lab 6 folder and copying the midpoint-rooted gene tree from Lab 5 to maintain consistency in our analysis. After copying, we verify that the gene tree file was transferred correctly. If necessary, we prune unwanted sequences using the gotree prune command to focus on the ingroup sequences of interest. Finally, we use nw_display to visualize the pruned tree, confirming that the structure is correct before proceeding with reconciliation. This streamlined preparation ensures that our analysis is based on a clean and accurate gene tree, ready for evolutionary investigation. 

## Reconcile the gene and species tree using Notung

In this section, we perform reconciliation of the STT3B gene tree using Notung. This analysis allows us to compare the gene tree with the species tree, helping us detect key evolutionary events such as duplications and losses. Below is the full workflow, including all necessary commands.

```
# Step 1: Create a directory for STT3B in the Lab 6 folder
mkdir ~/lab06-$MYGIT/STT3B

# Step 2: Copy the midpoint-rooted gene tree from Lab 5 to Lab 6
cp ~/lab05-$MYGIT/STT3B/STT3B.homologs.al.mid.treefile ~/lab06-$MYGIT/STT3B/STT3B.homologs.al.mid.treefile

# Step 3: Verify the copied file
ls ~/lab06-$MYGIT/STT3B/STT3B.homologs.al.mid.treefile

# Step 4: (Optional) Prune unwanted sequences from the gene tree
gotree prune -i ~/lab06-$MYGIT/STT3B/STT3B.homologs.al.mid.treefile \
-o ~/lab06-$MYGIT/STT3B/STT3B.homologs.pruned.treefile \
<sequence1> <sequence2> <sequence3>

# Step 5: Display the pruned gene tree (if pruning was performed)
nw_display ~/lab06-$MYGIT/STT3B/STT3B.homologs.pruned.treefile

# Step 6: Perform reconciliation using Notung
java -jar ~/tools/Notung-3.0_24-beta/Notung-3.0_24-beta.jar \
-s ~/lab05-$MYGIT/species.tre \
-g ~/lab06-$MYGIT/STT3B/STT3B.homologs.pruned.treefile \
--reconcile --speciestag prefix --savepng --events \
--outputdir ~/lab06-$MYGIT/STT3B/

# Step 7: Examine the reconciliation results for duplications and losses
less ~/lab06-$MYGIT/STT3B/STT3B.homologs.pruned.treefile.rec.events.txt

# Step 8: Display internal node names assigned by Notung (if needed)
grep NOTUNG-SPECIES-TREE ~/lab06-$MYGIT/STT3B/STT3B.homologs.pruned.treefile.rec.ntg | \
sed -e "s/^\[&&NOTUNG-SPECIES-TREE//" -e "s/\]/;/" | nw_display -

# Step 9: Generate RecPhyloXML object for detailed visualization
python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py \
-g ~/lab06-$MYGIT/STT3B/STT3B.homologs.pruned.treefile.rec.ntg \
--include.species

# Step 10: Create a reconciliation graphic using thirdkind
thirdkind -Iie -D 40 -f ~/lab06-$MYGIT/STT3B/STT3B.homologs.pruned.treefile.rec.ntg.xml \
-o ~/lab06-$MYGIT/STT3B/STT3B.homologs.pruned.treefile.rec.svg

# Step 11: Convert the SVG file to a PDF for easy viewing
convert -density 150 ~/lab06-$MYGIT/STT3B/STT3B.homologs.pruned.treefile.rec.svg \
~/lab06-$MYGIT/STT3B/STT3B.homologs.pruned.treefile.rec.pdf
```
This workflow involves setting up the analysis for the STT3B gene family by creating a new directory, copying over the midpoint-rooted gene tree from Lab 5, and optionally pruning any unwanted sequences. We then perform reconciliation using Notung, comparing the gene tree against the species tree to identify evolutionary events such as duplications and losses. The process also includes verifying the results and checking internal node names. Finally, we generate a RecPhyloXML object and create a detailed visualization using thirdkind, followed by converting the graphic to PDF format for easy viewing. This streamlined approach provides a clear and thorough analysis of the STT3B gene family, allowing us to understand its evolutionary changes across different species.

# 6. Protein-Domain
## Lab 8: Protein Domain Prediction

Clone this lab and move into this new directory. 

### Domain Identification 

In this lab, we will use RPS-BLAST to identify Pfam domains within unaligned protein sequences from the STT3B gene family. Pfam is a database for protein family classification, allowing us to detect key functional domains. This analysis will help us understand the structural features of STT3B.

The following commands will set up the necessary files, remove stop codons, and prepare the STT3B sequences for RPS-BLAST analysis using the Pfam database:

```
# Setting up the globin sequences directory
mkdir ~/lab08-$MYGIT/globins && cd ~/lab08-$MYGIT/globins

# Copy the raw, unaligned globin sequence and remove any stop codons (asterisks) using sed
sed 's/*//' ~/lab04-$MYGIT/globins/globins.homologs.fas > ~/lab08-$MYGIT/globins/globins.homologs.fas

# Note: The Pfam database required for RPS-BLAST analysis is pre-downloaded at ~/data/Pfam/
# The Pfam models will be used to classify protein families and predict functional domains.
```
We processed the unaligned STT3B protein sequences by removing stop codons. Using RPS-BLAST and the Pfam database, we identified functional domains in STT3B. The findings will provide insights into the evolutionary and functional significance of these domains.

### Running RPS-BLAST

The following command runs RPS-BLAST on the unaligned STT3B sequences to identify Pfam domains, generating an output file for further analysis.

```
rpsblast -query ~/lab08-$MYGIT/STT3B/STT3B.homologs.fas -db ~/data/Pfam/Pfam -out ~/lab08-$MYGIT/STT3B/STT3B.rps-blast.out -outfmt "6 qseqid qlen qstart qend evalue stitle" -evalue 0.0000000001
```
Review the RPS-BLAST output file. It contains results from the bioinformatics analysis, identifying conserved domains in the STT3B gene family across various species based on sequence homology.

#### E-value Discussion:
The e-value threshold influences the sensitivity of domain detection. A stricter e-value (e.g., 0.0000000001) ensures higher confidence in detected domains, while a higher e-value (e.g., 0.1) allows weaker domain hits, broadening the range of detected proteins. Use the strict e-value for precise questions but adjust as needed for broader exploration of your specific gene family.

# 7. Results

## Lab03: Finding homologs with BLAST KEY

### BLAST Search Summary

The BLAST analysis was conducted with an e-value filter set at 1e-30, targeting high-confidence homologs for the STT3B gene family. Below is a summary of the identified homologs across various species.

### Homolog Counts per Species:

The following table lists the number of STT3B homologs found for each species after filtering:

| Species      |   Count   |
|--------------|-----------|
| C.carcharias |    2       |
| C.mydas      |     2      |
| D.rerio      |     2       |
| E.caballus   |      2     |
| F.catus      |      2     |
| G.aculeatus  |      3     |
| G.gallus     |       2    |
| H.sapiens    |      2     |
| S.salar      |      6    |
| S.townsendi  |       2    |
| X.laevis     |      4     |

### Zero Homologs Check:

None of the analyzed species were found to have zero STT3B homologs in the filtered BLAST output.

### Notes for Analysis:

For a comprehensive study, aim to work with a dataset containing between 20 and 85 total homologs, with each species contributing 2 to 15 copies. If the homolog counts fall outside these ranges, you may need to adjust the e-value threshold to capture an appropriate number of matches. Please consult your TA before making any adjustments to ensure consistency in your analysis.

## Lab04: Gene family sequence alignment

The following table provides a summary of key metrics from the multiple sequence alignment and analysis of the STT3B gene family, highlighting sequence conservation and alignment characteristics.

#### Results Summary for STT3B Alignment

| Metric                                         | Value                                   |
|------------------------------------------------|-----------------------------------------|
| **Unaligned sequence file name and location**  | `STT3B.homologs.fas` at `~/lab04-$MYGIT/STT3B` |
| **Aligned sequence file name and location**    | `STT3B.homologs.al.fas` at `~/lab04-$MYGIT/STT3B` |
| **Length of the alignment**                    | 927 columns                             |
| **Number of columns with gaps**                | 695 columns                             |
| **Number of invariant columns**                | 539 columns                             |
| **Average percent identity (T-Coffee)**        | 76.69%                                  |
| **Average percent identity (AlignBuddy)**      | 71.24%                                  |

### Summary
The analysis shows high sequence conservation within the STT3B gene family, with substantial alignment length and percent identity, indicating strong evolutionary relationships among the homologs.

## Lab05 Gene Family Phylogeny using IQ-TREE

The STT3 gene family analysis focused on constructing a midpoint-rooted phylogenetic tree to understand the evolutionary history and divergence of the STT3A and STT3B subunits across various species. Using the best-fit model (Q.plant+I+R3) determined by BIC, we examined sequence alignments and evaluated the phylogenetic relationships with bootstrap support values. This analysis helped identify key divergence patterns and areas of low support, highlighting potential issues in certain taxa placements.

| **Parameter**                  | **Value/Description**                                         |
|--------------------------------|---------------------------------------------------------------|
| Best-fit Model                 | Q.plant+I+R3 (chosen based on BIC)                            |
| Rate Multiplier                | R3 (gamma categories)                                         |
| Frequencies                    | Empirical + invariant sites (I)                               |
| Low Bootstrap Nodes            | *Equus caballus* STT3A (0.000002), *Felis catus* STT3A (0.000002), *Homo sapiens* STT3B and *Equus caballus* STT3B (0.00499), *Salmo salar* LOC123727364 (0.000002) |
| Outgroup Selection             | No external evidence for outgroup selection                   |
| General Phylogeny Assessment   | Clear STT3A/STT3B separation with some unresolved nodes       |

The phylogenetic tree supports an ancient duplication event separating STT3A and STT3B subunits. High bootstrap values affirm the relationships in major clades, though low support in specific nodes suggests potential alignment or sequence data issues. Future work could refine these areas by incorporating more comprehensive data or alternative alignment methods.

## Lab06: Reconciling a Gene and Species Tree

The reconciliation of the STT3B gene family tree was performed using Notung, aligning the gene tree with the species tree based on the midpoint-rooted tree. The analysis aimed to identify duplication and loss events, providing insights into the evolutionary history of the STT3B gene family. Necessary pruning was applied to focus on relevant homologs.

### Reconciliation Analysis Results for the STT3B Gene Family

| **Parameter**                    | **Value/Description**                         |
|----------------------------------|-----------------------------------------------|
| Cost of Reconciled Tree          | 83.5                                          |
| Common Ancestor Copies (Gnathostomes) | 7 copies                                    |
| Events Leading to Tetrapods      | 2 duplications, 5 losses                      |
| Visualization Files              | `.svg` and `.pdf` files included in the repository |

The reconciliation analysis highlights significant duplication and loss events in the evolution of the STT3B gene family. These findings suggest complex evolutionary changes, especially during the transition from the common ancestor of Gnathostomes to Tetrapods.

## Lab08: Protein Domain Prediction









































