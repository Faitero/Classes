
---
title: "NGS Practical 3: German 2011 E. coli Outbreak"
author: "Igor Ruiz de los Mozos & Irilenia Nobeli"
date: '2021'
output:
  html_document: default
  pdf_document: default
---


# An introduction to NGS practical 3 (Genome Comparison)  

This part of the document was written to help students understand the context of the third NGS practical and should be read ideally *before* starting to work on that practical. It contains the background of the problem and an outline of the steps we aim to follow to solve the problem.

## Background  

In 2011, a very infectious strain of *E. coli* spread through Germany. Nearly 4000 people were affected and over 50 people died as a result. In such cases, the focus is often on identifying potential drug targets among the genes of the new strain. In this practical, we carry out a genome annotation of this strain and examine the conservation of antibiotic-resistant genes in this and other strains of *E. coli*.



## Obtaining a consensus sequence for the new strain  

We already have a file containing a number of long reads mapped to the *E. coli* reference genome (*Long.bam* from the previous practical. We can use the program `freebayes` to identify variants in the mapped sequences, i.e. positions in the reads that are different to the reference genome. The result is a VCF file containing a list of these variants. We can then use the program `vcf2fasta` to produce a consensus sequence (in fasta format), keeping the most likely nucleotide at each site where a variation is observed. Although variant calling in human sequences is a fairly complex process, the approach described here for bacteria is very easy to run for simple SNP calling. Note that this approach does not assume any ploidy (i.e. does not take into account different alleles arising from chromosomes of different origin in sexually reproducing organisms).

   

## Annotating the consensus sequence

We now have a FASTA file of the most likely genomic sequence of our new strain **but we have no information on where the genes are in this sequence**. Hence, the next step is to annotate the consensus sequence using the server DFAST at: https://dfast.nig.ac.jp/dfc/

This server allows us to annotate prokaryotic genomic sequences using homology searches.
The result isa  a number of files containing the predicted annotation of genes in our genomic sequence (genbank, gff etc).


   
## Producing a list of sequences for antibiotic-resistant genes   


The next step is to identify antibiotic-resistant genes and examine whether they are present in our strain and other strains of *E. coli*. This part is slightly more complicated due to the various file tranformations required to get to the final list of fasta sequences.

We will gather antibiotic-resistant genes for our list from two sources:   
      
 
- A list of relevant gene sequences has been obtained for you from an online resource (https://megares.meglab.org) and you simply need to clean it up and select a few terms as suggested in the practical guidelines.    
   
- A second list is created by searching the annotation of the strain we carried out in the previous step for terms referring to antibiotic resistance. This requires a bit more work as the original *gff* information must be turned into a *bed* file and finally the *fasta* file must be obtained from that, so there are several steps involved.   
 
  
Finally, the two lists must be concatenated into one *fasta* file which will be used as input for the comparison between genomes carried out in the next step.



## Comparing antibiotic-resistance genes across genomes using BRIG  
   
In this part of the practical the BLAST Ring Image Generator (BRIG) is used to visualise the results of comparing antibiotic-resistant genes in our strain and other *E. coli* strains. The input will be the concatenated list of sequences we created in the section above.
Most of the work here is done to ensure that the final result will offer a nice visualisation of the comparison of genes across different *E. coli* strains.

   
## BLAST Tutorial 
    
We expect that all students who have taken the Sequence Analysis and Genomics module were introduced to BLAST. For anyone not familiar with BLAST or anyone wanting to refresh their memory, we have put together a short tutorial that shows you how to use a single hypothetical gene sequence (resulting from an ORF prediction) to predicting its function and finding its orthologues (and paralogues) in NCBI databases using BLAST.   
  

# Practical 2 guidelines

## Set up your directory structure and working environment   

As in previous practicals, we will first login to *pandora* and then *thoth* (using `ssh -X`). We will set a shortcut for the path where we keep the practical work for this set of sessions. We will then create a new directory for this practical and work in that directory.

```
# I assume you have already connected to the server pandora.
# In any commands shown below, ex001 should be replaced with your own login name.

# connect to thoth from pandora using ssh
ssh -X thoth

# Change directory to the folder allocated to you for these practicals
cd /d/projects/u/ex001

# Set the path to your "omics" directory created in the previous practical.
setenv st_path /d/projects/u/user001/omics/

# Check that you have the following files created from the previous practical:
ls $st_path/results_NGS1/Long.bam
ls $st_path/course_materials/genomes/AFPN02.1/AFPN02.1_merge.fasta

# If you do not have these files, then copy them from my directories (see COPY FROM IRILENIA'S DIRECTORIES)

# Create a new directory called results_GC (for Genome Comparison) for the current practical.
mkdir $st_path/results_GC

# In that directory, create a directory called alignment and copy there the files we need
mkdir $st_path/results_GC/alignment
cp $st_path/results_NGS1/Long.bam $st_path/results_GC/alignment/
cp $st_path/course_materials/genomes/AFPN02.1/AFPN02.1_merge.fasta $st_path/results_GC/alignment/

# If you didn't have the files COPY FROM IRILENIA's DIRECTORY
cp /d/in4/u/ubcg71a/teaching/omics/results_NGS1/Long.bam $st_path/results_GC/alignment/
cp /d/in4/u/ubcg71a/teaching/omics/course_materials/genomes/AFPN02.1/AFPN02.1_merge.fasta $st_path/results_GC/alignment/

# Then change to the new alignment directory.
cd $st_path/results_GC/alignment/

```   



## Obtain a consensus sequence for the new strain  

We have a file containing a number of long reads mapped to the *E. coli* reference genome (*Long.bam* ) from the previous practical. We can use the program `freebayes` to identify variants in the mapped sequences with reference to the genome. The result is a VCF file containing a list of these variants. We can then use the program `vcf2fasta` to produce a consensus sequence (in fasta format), keeping the most likely nucleotide at each site where a variation is observed. 

```

# Load the relevant modules for freebayes and vcflib
module load freebayes
module load vcflib

# Get the consensus sequence with freebayes
# Create a variant call for all aligned nucleotides. We will create Single Nucleotide Polymorphisms and frecuencies for each nucleotide aligned.
freebayes -f AFPN02.1_merge.fasta -p 1 Long.bam > Long.vcffile

# Examine the VCF file (to exit the 'less' viewer, just press 'q').
less Long.vcffile

# Get a consensus fasta file from the VCF file 
# Merge all the possible variants to a single sequence and keep the most probable nucleotide for each site. 
vcf2fasta -f AFPN02.1_merge.fasta -P 1 Long.vcffile -p AFPN02.1_consensus

# Remame consensus file
mv AFPN02.1_consensusunknown_AFPN02.1_merge:0.fa AFPN02.1_merge_consensus.fa


# Examine the resulting consensus fasta file (To exit less press 'q').
less AFPN02.1_merge_consensus.fa

```   

** Note on the *fasta* format **

Fasta records (for gene, transcript or protein) start with a ">" symbol followed by the header of the sequence.
The next line contains the sequence itself, for example:  

```
>Sequence header
GTCACCAGCTGTCGCCCAGCAAAACGATGATAATGAGATCATAGTGTCTGCCAGCCGCAGCAATCGAACTGTAGCGGAGA
TGGCGCAAACCACCTGGGTTATCGAAAATGCCGAACTGGAGCAGCAGATTCAGGGCGGTAAAGAGCTGAAAGACGCACTG
GCTCAGTTAATCCCCGGCCTTGATGTCAGCAGCCAGAGTCGAACCAACTACGGTATGAACATGCGTGGCCGTCCACTGGT
TGTCCTGATTGACGGTGTGCGCCTCAACTCTTCACGTTCCGACAGCCGACAACTGGACTCTGTCGATCCTTTTAATATCG
ACCATATTGAAGTGATCTCCGGCGCGACGGCCCTGTACGGTGGCGGGAGTACCGGAGGGTTGATCAACATCGTGACCAAA
AAAGGTCAGCCGGAAACCATGATGGAGTTTGAGGCTGGCACAAAAAGTGGCTTTAACAGCAGTAAAGATCACGATGAGCG
CATTGCCGGTGCTGTCTCCGGCGGAAATGACCATATCTCCGGACGTCTTTCCGTGGCATATCAGAAATTTGGCGGCTGGT...
>Next record
GCGCAAACCACCTGGGTTATCGAAAATGCCGAACTGGAGCAGCAGATTCAGGGCGGTAAAGAGCTGAAAGAACTGTACGT...

```

The consensus fasta file we created should look like this:
![ fasta format ](figures/consensus_fasta.png)


## Annotate the consensus sequence

We now have a FASTA file of the most likely genomic sequence of our strain but we have no information on where the genes are in this sequence. Hence, the next step is to annotate the consensus sequence using the server DFAST at: https://dfast.nig.ac.jp/dfc/

This produces a number of files containing the predicted annotation of genes in our genomic sequence (genbank, gff etc).

*You can try to annotate your consensus sequence with DFast. But for time efficiency we provide you with a precompiled annotation (see code further down)*

**THIS PART NEEDS UPDATING AS DFAST HAS CHANGED ITS INTERFACE SLIGHTLY**

If you want to repeat the annotation, go to the server, browse to your consensus sequence, provide a job title and your email address. Select "E. coli" additional DB, "100" minimum sequence length, "Enable HMM scan against TIGRFAM", "Enable RPSBLAST against COG" and "Rotate/flip the chromosome so that the dnaA gene comes first".


![ DFast](figures/dfast.png)
![ DFast2](figures/dfast2.png)

When the job is finished, download the "annotation.zip" file and save it in the directory results_NGS1.
However, we advise you to try this another time and use instead the annotation provided by us to save time.

```
# Change directory
cd "$st_path"/results_GC/

# Download the annotated file from here:
wget "https://www.dropbox.com/s/yh19885owchp94h/annotation.zip" -O annotation.zip

# Unzip the annotation folder
unzip annotation.zip

```

### Annotation formats
We can take a quick look at the various files that DFast produced.

```
# Change directory to the annotation results
cd "${st_path}"/results_GC/annotation/

## Explore the file formats present on the annotation folder

head annotation.gbk
head annotation.gff
head cds.fna
head features.tsv
head protein.faa
head rna.fna
head statistics.txt

```
A quick explanation of the files produced is given below.

**GenBank**

Designed to store annotation, genomic position, features, metadata, taxonomy, sequence and reference authors.
Flat file format with unique identifiers for proteins and genes.  
  
https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html.


![ GenBank ](figures/GeneBank.png)


**GFF annotation format**   

Designed to store genomic position, annotation and metadata.
  
https://www.ensembl.org/info/website/upload/gff.html

![ gff ](figures/GFF.png)

**CDS.fna**
  
Multifasta format with protein sequences.  
  
**features.tsc**  
   
Table with genomic coordinates and minimal annotation.
  


## Produce a list of sequences for antibiotic-resistant genes   

Now we have annotated the sequence of our strain, the next step is to identify antibiotic-resistant genes and examine whether they are present in our strain and other strains of *E. coli*. This part is slightly more complicated due to the various file tranformations required to get to the final list of fasta sequences.

We will gather antibiotic-resistant genes for our list from two sources:   
       
1) A list of relevant gene sequences has been obtained for you from an online resource (https://megares.meglab.org) and you simply need to clean it up and select a few terms as suggested in the practical guidelines.    
   
2) A second list is created by searching the annotation of the strain we carried out in the previous step for terms referring to antibiotic resistance. This requires a bit more work as the original *gff* information must be turned into a *bed* file and finally the *fasta* file must be obtained from that, so there are several steps involved.   
  
Finally, **the two lists must be concatenated into one *fasta* file** which will be used as input for the comparison between genomes carried out in the next step.

You will need to use the `cat` command to concatenate files and direct the output to the screen or to a new file, so below is a brief reminder of how this works.
#### Examples of cat  
  
The command `cat` will concatenate files and redirect to standard output. We will be using this command below. 

https://www.gnu.org/software/coreutils/manual/html_node/cat-invocation.html  

  
![ cat GNU ](figures/cat.png)   


We will now proceed to build the list of antibiotic-resistant genes.

### 1) Antibiotic resistance genes from a database
  
Antibiotic resistance genes have been downloaded for you from https://megares.meglab.org .  

We are going to extract genes from the downloaded file. We will __grep__ antibiotic terms and then clean the symbols "--" and the fasta header to obtain the final file: **selected\_antibiotic\_resistance\_genes\_meglab.fasta**  

To give you an overview of how we will do this, below is a schematic of the steps we will take, followed by the actual code: 
   
![ Antibiotic resistance genes from database](figures/selected.png) 

```
# Start off by creating an "antibiotics" directory in the results_GC directory and then
# copy the file we obtained from the database to your own directories.

mkdir "${st_path}"/results_GC/antibiotics
cd "${st_path}"/results_GC/antibiotics
cp /d/in4/u/ubcg71a/teaching/omics/results_GC/antibiotics/antibiotic_resistance_genes_meglab.fasta .

# Read file obtained from meglab.org containing antibiotic resistance genes in fasta format.
less antibiotic_resistance_genes_meglab.fasta

# Parse antibiotic genes
# Below we use cat in combination with other commands to progressively get what we want out of that file.

# Concatenate file (continue beyond end of line "\n")
# cat manual can be found here:
#    https://www.gnu.org/software/coreutils/manual/html_node/cat-invocation.html
cat antibiotic_resistance_genes_meglab.fasta | head -5

# The file contains data from many bacteria. We will extract the lines that contain
# the word "Escherichia" using grep
# grep manual https://www.gnu.org/software/grep/manual/grep.html
# Will output all the lines that contain 'Escherichia' string 
# Note the use of | to pipe between commands
cat antibiotic_resistance_genes_meglab.fasta | grep 'Escherichia'

# Grep a single word and count number of lines
# wc manual https://www.gnu.org/software/coreutils/manual/html_node/wc-invocation.html
cat antibiotic_resistance_genes_meglab.fasta | grep 'Escherichia' | wc -l

# Grep multiple words (-E) and count number of lines (total count 40 genes)
cat antibiotic_resistance_genes_meglab.fasta | grep -E 'Escherichia|carbapene|CEPH|NDM|QnrB9|Metronidazole' | wc -l

# Grep multiple words and output also next line -nucleotide sequence (-A1)
cat antibiotic_resistance_genes_meglab.fasta | grep -A1 -E 'Escherichia|carbapene|CEPH|NDM|QnrB9|Metronidazole' | wc -l

# We have more than double the number of results
cat antibiotic_resistance_genes_meglab.fasta | grep -A1 -E 'Escherichia|carbapene|CEPH|NDM|QnrB9|Metronidazole'

# Inspect output!! What is there? We note that the line containing "--" is included in the output.
# We will need to remove it and will do so below using "sed".
```

#### Manipulating text with sed  
  
```

# Let's see quickly how sed substitution works

# Print a string
echo "line to play with sed substitute"

# substitute "play" string
# substitute in sed have the following patern 's/patern_to_search/substitute_to/g'. Note s (substitute) and g at the end (global)
# Substitute play with learn
echo "line to play with sed substitute" | sed 's/play/learn/g'
# Substitute space with _
echo "line to play with sed substitute" | sed 's/ /_/g'
# Substitute space with nothing
echo "line to play with sed substitute" | sed 's/ //g'

```   
   

#### Continuing the editing of file with antibiotic resistance genes from database 
  
```
# Make sure we are in the right directory.
cd "${st_path}"/results_GC/antibiotics

# To remove the "--" separator, we use the sed command:
# sed manual https://www.gnu.org/software/sed/manual/sed.html

cat antibiotic_resistance_genes_meglab.fasta | grep -A1 -E 'Escherichia|carbapene|CEPH|NDM|QnrB9|Metronidazole' | sed '/^--$/d' | wc -l

# Let's see what we just did...
# We are using regular expressions (abbreviated to regex) to search 
# (/whatever in between forward slashes/) for '--' at the start of the 
# line (denoted by '^') and at the end of the line (denoted by '$').
# Regex are quite complex even for computer scientists and are out of the scope of this tutorial, 
# but if you want to learn more:
# manuals https://en.wikipedia.org/wiki/Regular_expression
# manual and try it https://medium.com/factory-mind/regex-tutorial-a-simple-cheatsheet-by-examples-649dc1c3f285
# test your regex https://regexr.com/

# This time we have 80 lines so we can redirect output (>) to a new file
cat antibiotic_resistance_genes_meglab.fasta | grep -A1 -E 'Escherichia|carbapene|CEPH|NDM|QnrB9|Metronidazole' | sed '/^--$/d' > selected_antibiotic_resistance_genes_meglab.fasta


# Let's now inspect the headers of this fasta file:
head -n 1 selected_antibiotic_resistance_genes_meglab.fasta

# output is: ">Met|nimB_1_X71443|Metronidazole|nim_nitroimidazole_reductase|NIMB"
# This fasta header contains more information that we need plus it will be 
# cumbersome to use later on.
# We will need to edit it.  
# We need to keep only the ">" fasta symbol, the first gene acronym "Met" 
# and the antibiotic family "NIMB"

# We will use awk specifiying "|" as a field/column separator "FS" and then print the 
# first element $1 and the last element $NF
# Do you notice the difference?
cat selected_antibiotic_resistance_genes_meglab.fasta | awk 'BEGIN {FS="|"}  {print $1 "_" $NF }' | head -n 1

# We output the result to temp.fasta using (>)
cat selected_antibiotic_resistance_genes_meglab.fasta | awk 'BEGIN {FS="|"}  {print $1 "_" $NF }' > temp.fasta
# and move "temp.fasta" file to "selected_antibiotic_resistance_genes.fasta" 
mv temp.fasta selected_antibiotic_resistance_genes_meglab.fasta
# Answer "y" to "mv: overwrite ‘selected_antibiotic_resistance_genes.fasta’?"

# Inspect final selected antibiotic resistance genes downloaded from meglab.org database
less selected_antibiotic_resistance_genes_meglab.fasta

```   

###  2) Finding antibiotic genes in the predicted annotation    
   
Starting from the annotation of our genome AFPN02.1 we will **grep** antibiotic terms. Then we will transform the gff file to a bed file. We will use the bed file to retrieve the nucleotide sequence helped with BedTools. Finally we will clean the fasta header to obtain the file `present_in_AFPN02_antibiotic_resistance_genes.fasta`. The process is summarised in the diagram below:

  
![ Antibiotic resistance genes from database](figures/present.png)  

```
# Go to the annotation directory
cd "$st_path"/course_materials/results_GC/annotation

# Awk is stand-alone scripting language. You can do incredible tasks with it in a single line. 
# Replicating the same with python will take several lines of code > perl > R
# If you want to refresh your memory of awk from the first practical, see here:
# awk manual https://en.wikipedia.org/wiki/AWK
# gawk manual https://www.gnu.org/software/gawk/manual/gawk.html

# The easiest awk call is to retrieve the 1st column. It will automatically detect field/column separator in a tabular text file. In this case it separates columns by "blank space" (default)
awk '{print $1}' annotation.gff | head -n 2
# Compare it with the original. Note that with head we are only taking the first 2 lines (-n 2)
head -n 2 annotation.gff


# awk call to retrieve 1st and 2nd columns.
# On the print output statement we include " " space as separator.
awk '{print $1 " " $2}' annotation.gff | head -n 2
# Compare it with the original. 
head -n 2 annotation.gff


# awk call to retrieve 1st and 2nd columns but change the field separator 
# to a tab via the command FS="\t" - sometimes this will be necessary to 
# parse the text. Note that on the print statement we also separate fields with tab "\t".
awk 'BEGIN {FS="\t"} {print $1 "\t" $2}' annotation.gff | head -n 2
# Compare it with the original.
head annotation.gff | head -n 2

```
Next, we want to tranform the gff file to the bed format so let's take a quick look at bed below.

#### A note on the BED format

BED is a flat file format separated with tab, used to store, retieve and interact with genomic positions. 
It usually includes some annotation.
Check out the manual here:  https://www.ensembl.org/info/website/upload/bed.html   

```
Required fields:

The first three fields in each feature line are required:

    chrom - name of the chromosome or scaffold. Any valid seq_region_name can be used, and chromosome names can be given with or without the 'chr' prefix.
    chromStart - Start position of the feature in standard chromosomal coordinates (i.e. first base is 0).
    chromEnd - End position of the feature in standard chromosomal coordinates

chr1  213941196  213942363
chr1  213942363  213943530
chr1  213943530  213944697
chr2  158364697  158365864
chr2  158365864  158367031
chr3  127477031  127478198
chr3  127478198  127479365
chr3  127479365  127480532
chr3  127480532  127481699

Optional fields

Nine additional fields are optional. Note that columns cannot be empty - lower-numbered fields must always be populated if higher-numbered ones are used.

    name - Label to be displayed under the feature, if turned on in "Configure this page".
    score - A score between 0 and 1000. See track lines, below, for ways to configure the display style of scored data.
    strand - defined as + (forward) or - (reverse).
    thickStart - coordinate at which to start drawing the feature as a solid rectangle
    thickEnd - coordinate at which to stop drawing the feature as a solid rectangle
    itemRgb - an RGB colour value (e.g. 0,0,255). Only used if there is a track line with the value of itemRgb set to "on" (case-insensitive).
    blockCount - the number of sub-elements (e.g. exons) within the feature
    blockSizes - the size of these sub-elements
    blockStarts - the start coordinate of each sub-element

chr7  127471196  127472363  Pos1  0  +  127471196  127472363  255,0,0
chr7  127472363  127473530  Pos2  0  +  127472363  127473530  255,0,0
chr7  127473530  127474697  Pos3  0  +  127473530  127474697  255,0,0
chr7  127474697  127475864  Pos4  0  +  127474697  127475864  255,0,0
chr7  127475864  127477031  Neg1  0  -  127475864  127477031  0,0,255
chr7  127477031  127478198  Neg2  0  -  127477031  127478198  0,0,255
chr7  127478198  127479365  Neg3  0  -  127478198  127479365  0,0,255
chr7  127479365  127480532  Pos5  0  +  127479365  127480532  255,0,0
chr7  127480532  127481699  Neg4  0  -  127480532  127481699  0,0,255


```   

#### Transform GFF format to bed  
   
```
# Go to the annotation directory
cd "${st_path}"/results_GC/annotation

# Grep lines with the word  "antibiotic"
cat annotation.gff | grep -E 'antibiotic'


# Transform gff file to bed 6 columns
# Grep lines with "antibiotic" string and transform to bed format
# on the print statement we can include any string between quotes like I'm doing with "sequence1"
cat annotation.gff | grep -E 'antibiotic' | awk 'BEGIN {FS="\t"};  {print "sequence1" "\t" $4 "\t" $5 "\t" $3 "\t" $6 "\t" $7}'


# Grep several antibiotics' names and transform to bed with awk.

# Note how we use here the split command. We are splitting field $9 (9th) into everything 
# that start with "="" and end with ";".
# We are trying to capture the gene name that has the format ";gene=gene_name;".
# Then we assign captured array to "captured" variable and require that the split 
# has more than 10 elements ">=10".
# Then we select the 10th element "captured[10]", which is gene acronym.
# We also retrieve the "captured[4]" element that will be the protein product.
cat annotation.gff | grep -E 'antibio|penici|lactama|macrolide|tetracycli|streptomycin|sulfonamide|ampheni|tetracycline|Tellurite|nalidixic' \
| awk 'BEGIN {FS="\t"}  split($9, captured, /[(=);]/) >=10  {print "sequence1" "\t" $4 "\t" $5 "\t" captured[10] "\t" captured[4] "\t" $7}'


# Send output to a bed file "present_antiotic_resistance_genes.bed"
cat annotation.gff | grep -E 'antibio|penici|lactama|macrolide|tetracycli|streptomycin|sulfonamide|ampheni|tetracycline|Tellurite|nalidixic' \
| awk 'BEGIN {FS="\t"}  split($9, captured, /[(=);]/) >=10  {print "sequence1" "\t" $4 "\t" $5 "\t" captured[10] "\t" captured[4] "\t" $7}' \
> present_in_AFPN02_antibiotic_resistance_genes.bed

# BedTools
# Once we have the annotation and genomic position we need to include the fasta sequence; 
# for this task we will use bedtools.
# Load the module bedtools so that all relevant executables can be found by the shell.
module load bedtools


# bedtools manual https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html
# bedtools getfasta manual https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html
bedtools getfasta -name -s -fi "$st_path"/results_GC/annotation/genome.fna -bed present_in_AFPN02_antibiotic_resistance_genes.bed -fo present_in_AFPN02_antibiotic_resistance_genes.fasta


# Note the addition of nucleotide sequence in the fasta file
head -n 1 present_in_AFPN02_antibiotic_resistance_genes.bed
head -n 2 present_in_AFPN02_antibiotic_resistance_genes.fasta


# Realize that we have ">tetA(-)" in the header of the fasta.
head -n 1 present_in_AFPN02_antibiotic_resistance_genes.fasta

# We need to modify the header to end up with ">tetA" only. 
# (i.e. just the fasta beginning ">" and gene symbol "tetA")
# This can be done with another sed one-liner job.
# Find "(" -  parentheses. Then match anything (.) zero or more times (*) and substitute with nothing.
cat present_in_AFPN02_antibiotic_resistance_genes.fasta | sed 's/(.*//g' > temp.fasta

# Move the temporal file to present_in_AFPN02_antibiotic_resistance_genes.fasta
mv temp.fasta present_in_AFPN02_antibiotic_resistance_genes.fasta

# Move all the antibiotic files to antibiotics folder
mv present_in_AFPN02_antibiotic_resistance_genes*  "${st_path}"/results_GC/antibiotics/

``` 

### Merge the two fasta files of antibiotics

Finally we will merge **selected** genes from the meglab database and those that we have queried from the annotation (**present**) in a final file: **"final\_comparison\_antibiotics.fasta "**.

![ Antibiotic resistance genes from database](figures/merge.png) 


```
cd "${st_path}"/results_GC/antibiotics/

# Merge the two antibiotics-resistance files we have obtained.

cat present_in_AFPN02_antibiotic_resistance_genes.fasta selected_antibiotic_resistance_genes_meglab.fasta > final_comparison_antibiotics.fasta

# Copy all the new files to wholeGenomeExamples folder
# First, create this new folder

mkdir "${st_path}"/results_GC/wholeGenomeExamples
cp present_in_AFPN02_antibiotic_resistance_genes.fasta selected_antibiotic_resistance_genes_meglab.fasta final_comparison_antibiotics.fasta "${st_path}"/results_GC/wholeGenomeExamples

```
