Day 2 Afternoon
===============
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

High-throughput BLAST and pan-genome analysis
---------------------------------------------

This morning we learned how to perform basic genome annotation and comparison using Prokka and ACT. Now we will up the ante and do some more sophisticated comparative genomics analyses! 
First, we will create custom BLAST databases to identify specific antibiotic resistance genes of interest in a set of genomes. 
Second, we will use the tool [ARIBA](https://github.com/sanger-pathogens/ariba/wiki) to identify the complete antibiotic resistome in our genomes. 
Third, we will move beyond antibiotic resistance, and look at the complete set of protein coding genes in our input genomes. 
Finally, we will go back to ACT to understand the sorts of genomic rearrangements underlying observed variation in gene content.

For BLAST and ARIBA, we will be looking at 8 *Klebsiella pneumoniae* genomes from human and environmental sources. Six of these genomes are from [this paper](https://www.pnas.org/content/112/27/E3574), and the other two are sequences from our lab. We are interested in learning more about potential differences in the resistomes of human and environmental isolates.

For the pan-genome analysis, we will be looking at four closely related *Acinetobacter baumannii* strains. However, despite being closely related, these genomes have major differences in gene content, as *A. baumannii* has a notoriously flexible genome! In fact, in large part due to its genomic flexibility, *A. baumannii* has transitioned from a harmless environmental contaminant to a pan-resistant super-bug in a matter of a few decades. If you are interested in learning more, check out this nature [review](http://www.nature.com/nrmicro/journal/v5/n12/abs/nrmicro1789.html) or [this](http://www.pnas.org/content/108/33/13758.abstract) paper I published a few years back analyzing the very same genomes you will be working with.

Execute the following command to copy files for this afternoon’s exercises to your scratch directory:

```  

cd /scratch/micro612w19_fluxod/username

or

wd

cp -r /scratch/micro612w19_fluxod/shared/data/day2_after/ ./

```

Determine which genomes contain KPC genes using [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
----------------------------------------------------
[[back to top]](day2_afternoon.html)
[[HOME]](index.html)

Before comparing full genomic content, lets start by looking for the presence of particular genes of interest. Some *K. pneumoniae* harbor a KPC gene that confers resistance to carbapenems, a class of antibiotics of last resort (more information [here](https://www.sciencedirect.com/science/article/pii/S1473309913701907?via%3Dihub) and [here](https://academic.oup.com/jid/article/215/suppl_1/S28/3092084)). We will see if any of our samples have a KPC gene, by comparing the genes in our genomes to KPC genes extracted from the antibiotic resistance database ([ARDB](http://ardb.cbcb.umd.edu/)). These extracted genes can be found in the file `ardb_KPC_genes.pfasta`, which we will use to generate a BLAST database.

> ***i. Run makeblastdb on the file of KPC genes to create a BLAST database.***

makeblastdb takes as input: 

1) an input fasta file of protein or nucleotide sequences (`ardb_KPC_genes.pfasta`) and 

2) a flag indicating whether to construct a protein or nucleotide database (in this case protein: `-dbtype prot`).

```
#change directory to day2_after
d2a


makeblastdb -in ardb_KPC_genes.pfasta -dbtype prot

```

> ***ii. BLAST K. pneumoniae protein sequences against our custom KPC database.***

Run BLAST! 

The input parameters are: 

1) query sequences (`-query kpneumo_all.pfasta`), 

2) the database to search against (`-db ardb_KPC_genes.pfasta`), 

3) the name of a file to store your results (`-out KPC_blastp_results`), 

4) output format (`-outfmt 6`), 

5) e-value cutoff (`-evalue 1e-100`), 

6) number of database sequences to return (`-max_target_seqs 1`) (Note that when using large databases, this might not give you the best hit. See [here](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty833/5106166) for more details.)


```
blastp -query kpneumo_all.pfasta -db ardb_KPC_genes.pfasta -out KPC_blastp_results.tsv -outfmt 6 -evalue 1e-100 -max_target_seqs 1
```

Use `less` to look at `KPC_blastp_results.tsv`.

```
less KPC_blastp_results.tsv
```

- **Exercise:** Experiment with the `–outfmt` parameter, which controls different output formats that BLAST can produce. 

- **Exercise:** In this exercise you will try a different type of blasting – blastx. Blastx compares a nucleotide sequence to a protein database by translating the nucleotide sequence in all six frames and running blastp. Your task is to determine which Enterococcus  genomes are vancomycin resistant by blasting against a database of van genes. The required files are located in VRE_van_blast folder under day2_after directory.

Your steps should be:

1) Concatenate .fasta files (VRE/VSE genomes) into a single file (your blast query file)
2) Create a blastp database from ardb_van.pfasta
3) Run blastx
4) Verify that only the VRE genomes hit the database
5) For extra credit, determine which van genes were hit by using grep to search for the hit gene ID in ardb_van.pfasta



<details>
  <summary>Solution</summary>
  
```
cat VRE_*.fasta VSE_ERR37*.fasta > VRE_VSE_genomes.fasta

makeblastdb -in ardb_van.pfasta -dbtype prot

blastx -query VRE_VSE_genomes.fasta -db ardb_van.pfasta -out van_blastp_results.tsv -outfmt 6 -evalue 1e-100 -max_target_seqs 1

```
</details>

Identify antibiotic resistance genes with [ARIBA](https://github.com/sanger-pathogens/ariba) directly from paired end reads
----------------------------------------------------------
[[back to top]](day2_afternoon.html)
[[HOME]](index.html)

Now let's look at the full spectrum of antibiotic resistance genes in our *Klebsiella* genomes!

ARIBA (Antimicrobial Resistance Identification By Assembly) is a tool that identifies antibiotic resistance genes by running local assemblies. The input is a FASTA file of reference sequences (can be a mix of genes and noncoding sequences) and paired sequencing reads. ARIBA reports which of the reference sequences were found, plus detailed information on the quality of the assemblies and any variants between the sequencing reads and the reference sequences.

ARIBA is compatible with various databases and also contains a utility to download different databases such as: argannot, card, megares, plasmidfinder, resfinder, srst2_argannot, vfdb_core. Today, we will be working with the [card](https://card.mcmaster.ca/) database, which has been downloaded and placed in the `/scratch/micro612w19_fluxod/shared/bin/ariba/database/CARD/` directory.

<!---
Note: There is an issue with downloading the database. They are in a process to fix the broken CARD database link issue. For now, I am using my own downloaded database.
>i. Run this command to download CARD database
```
/nfs/esnitkin/bin_group/anaconda3/bin/python /nfs/esnitkin/bin_group/ariba/scripts/ariba getref card out.card
```
>ii. Prepare this downloaded card database for ARIBA
```
/nfs/esnitkin/bin_group/anaconda3/bin/python /nfs/esnitkin/bin_group/ariba/scripts/ariba prepareref -f out.card.fa -m out.card.tsv out.card.prepareref
```
-->

> ***i. Run ARIBA on input paired-end fastq reads for resistance gene identification.***

The fastq reads are placed in the `kpneumo_fastq` directory. Since ARIBA is a memory intensive, we will enter an interactive flux session to run this exercise. Start the interactive session and change directories to `day2_after` and run the short for loop that will  commands below to start ARIBA jobs in the background.

<!---
module load python-anaconda3/latest
-->

```
iflux

cd /scratch/micro612w19_fluxod/username/day2_after

#or 

d2a

#Load dependencies

module load python-anaconda3/latest-3.6   bowtie2/2.1.0   cd-hit/4.6.4   mummer/3.23  ariba/2.13.3


#ARIBA commands

# List forward end fastq files in the directory and save the filenames into the variable samples. 
samples=$(ls kpneumo_fastq/*1.fastq.gz) #forward reads

# Set ARIBA dabase directory to the CARD database that we downloaded in the below folder
db_dir=/scratch/micro612w19_fluxod/shared/bin/ariba/database/CARD/out.card.proteus.prepareref/ #reference database

# Run for loop, where it generates ARIBA command for each of the forward end files.
for samp in $samples; do   
samp2=${samp//1.fastq/2.fastq} #reverse reads   
outdir=$(echo ${samp//.fastq.gz/} | cut -d/ -f2) #output directory 
ariba run --force $db_dir $samp $samp2 $outdir & #ariba command 
done


```

The "&" in the above commands(at the end) is a little unix trick to run commands in background. 

You can run multiple commands in background and make full use of parallel processing. You can check the status of these background jobs by typing:

```
jobs
```

> ***ii. Run ARIBA summary function to generate a summary report.***

ARIBA has a summary function that summarises the results from one or more sample runs of ARIBA and generates an output report with various level of information determined by the `-preset` parameter. The parameter `-preset minimal` will generate a minimal report showing only the presence/absence of resistance genes whereas `-preset all` will output all the extra information related to each database hit such as reads and reference sequence coverage, variants and their associated annotations (if the variant confers resistance to an antibiotic) etc.

```

ariba summary --preset minimal kpneumo_ariba_minimal_results */report.tsv

ariba summary --preset all kpneumo_ariba_all_results */report.tsv

```

The ARIBA summary generates three output:

1. kpneumo_ariba*.csv file that can be viewed in your favourite spreadsheet program.
2. kpneumo_ariba*.phandango.{csv,tre} that allow you to view the results in [Phandango](http://jameshadfield.github.io/phandango/#/). You can drag-and-drop these files straight into Phandango.

Lets copy these  files, along with a metadata file, to the local system using cyberduck or scp.

```
scp username\@flux-xfer.arc-ts.umich.edu:/scratch/micro612w19_fluxod/username/day2_after/kpneumo_ariba* ~/Desktop/
scp username\@flux-xfer.arc-ts.umich.edu:/scratch/micro612w19_fluxod/username/day2_after/kpneumo_source.tsv ~/Desktop/
```

Drag and drop these two files onto the [Phandango](http://jameshadfield.github.io/phandango/#/) website. What types of resistance genes do you see in these *Klebsiella* genomes? This [review]() may help interpret.

> ***iii. Explore full ARIBA matrix in R***

- Now, fire up RStudio and read in the ARIBA full report "kpneumo_ariba_all_results.csv"

```
ariba_full  = read.csv(file = '~/Desktop/kpneumo_ariba_all_results.csv', row.names = 1)
rownames(ariba_full) = gsub('_1|_R1|/report.tsv','',rownames(ariba_full))
```

- Subset to get description for each gene

```
ariba_full_match = ariba_full[, grep('match',colnames(ariba_full))]
```

- Make binary for plotting purposes

```
ariba_full_match[,] = as.numeric(ariba_full_match != 'no')
```

- Make a heatmap!

```
# install pheatmap
install.packages('pheatmap')

# load pheatmap
library(pheatmap)

# load metadata about sample source
annots = read.table('~/Desktop/kpneumo_source.tsv',row.names=1)
colnames(annots) = 'Source'

# plot heatmap
pheatmap(ariba_full_match,annotation_row = annots)
```

- **Exercise:** Bacteria of the same species can be classified into different sequence types (STs) based on the sequence identity of certain housekeeping genes using a technique called [multilocus sequence typing (MLST)](https://en.wikipedia.org/wiki/Multilocus_sequence_typing). The different combination of these house keeping sequences present within a bacterial species are assigned as distinct alleles and, for each isolate, the alleles at each of the seven genes define the allelic profile or sequence type (ST). Sometimes, different sequence types are associated with different environments or different antibiotic resistance genes. We want to know what sequence type(s) our genomes come from, and if there are certain ones that are associated with certain sources or certain antibiotic resistance genes. 

Using the [ARIBA MLST manual](https://github.com/sanger-pathogens/ariba/wiki/MLST-calling-with-ARIBA), write and run a script (similar to the one above) to perform MLST calling with ARIBA on all 8 of our *K. pneumonia* genomes. Then, use this information to add a second annotation column to the heatmap we created above to visualize the results. Running ARIBA mlst requires a MLST species database, so dont forget to download "Klebsiella pneumoniae" databsae and give the path to these database while running ARIBA MLST detection.

Did you find anything interesting?


<details>
  <summary>Solution</summary>
  
```
# Load modules
module load python-anaconda3/latest-3.6   bowtie2/2.1.0   cd-hit/4.6.4   mummer/3.23  ariba/2.13.3

# Check if you have a mlst database for your species of interest
ariba pubmlstspecies

# Download your species mlst database
ariba pubmlstget "Klebsiella pneumoniae" get_mlst

# Set ARIBA database directory to the get_mlst database that we just downloaded.
db_dir=./get_mlst/ref_db/

# Run ariba mlst with this database
samples=$(ls kpneumo_fastq/*1.fastq.gz) #forward reads

# Run for loop, where it generates ARIBA command for each of the forward end files.
for samp in $samples; do   
samp2=${samp//1.fastq/2.fastq} #reverse reads   
outdir=$(echo ${samp//.fastq.gz/} | cut -d/ -f2) #output directory 
ariba run --force $db_dir $samp $samp2 $outdir & #ariba command 
done

# Once the run completes, run summarize_mlst.sh script to print out mlst reports that are generated in current directory
bash summarize_mlst.sh .

```
</details>


Perform pan-genome analysis with [Roary](https://sanger-pathogens.github.io/Roary/)
----------------------------------------

Roary is a pan genome pipeline, which takes annotated assemblies in GFF3 format and calculates the pan genome. The pan-genome is just a fancy term for the full complement of genes in a set of genomes. 

The way Roary does this is by: 
1) Roary gets all the coding sequences from GFF files, convert them into protein, and create pre-clusters of all the genes, 
2) Then, using BLASTP and MCL, Roary will create gene clusters, and check for paralogs. and 
3) Finally, Roary will take every isolate and order them by presence/absence of genes.

> ***i. Generate pan-genome matrix using Roary and GFF files***

Make sure you are on an interactive node, as this will be even more computationally intensive!

```
iflux
```

Change your directory to day2_after

```

> Make sure to change username with your uniqname

cd /scratch/micro612w19_fluxod/username/day2_after/

or 

d2a

```

Load all the required dependencies and run roary on GFF files placed in Abau_genomes_gff folder.

```
module load samtools
module load bedtools2
module load cd-hit
module load ncbi-blast
module load mcl
module load parallel
module load mafft
module load fasttree
module load perl-modules
module load R/3.3.0
module load roary

#Run roary
roary -p 4 -f Abau_genomes_roary_output -r -n -e -v Abau_genomes_gff/*.gff 
```

The above roary command will run pan-genome pipeline on gff files placed in Abau_genomes_gff(-v) using 4 threads(-p), save the results in an output directory Abau_genomes_roary_output(-f), generate R plots using .Rtab output files and align core genes(-n)

Change directory to Abau_genomes_roary_output to explore the results.

```
cd Abau_genomes_roary_output

ls
```

Output files:

1. summary_statistics.txt: This file is an overview of your pan genome analysis showing the number of core genes(present in all isolates) and accessory genes(genes absent from one or more isolates or unique to a given isolate). 

2. gene_presence_absence.csv: This file contain detailed information about each gene including their annotations which can be opened in any spreadsheet software to manually explore the results. It contains plethora of information such as gene name and their functional annotation, whether a gene is present in a genome or not, minimum/maximum/Average sequence length etc.

3. gene_presence_absence.Rtab: This file is similar to the gene_presence_absence.csv file, however it just contains a simple tab delimited binary matrix with the presence and absence of each gene in each sample. It can be easily loaded into R using the read.table function for further analysis and plotting. The first row is the header containing the name of each sample, and the first column contains the gene name. A 1 indicates the gene is present in the sample, a 0 indicates it is absent.

4. core_gene_alignment.aln: a multi-FASTA alignment of all of the core genes that can be used to generate a phylogenetic tree.

<!--
#Plots are not very useful. Seems like a waste of time.
>ii. Generate a phylogenetic tree and plot pan-genome matrix.
we will use core_gene_alignment.aln multi-fasta core gene alignment as an input to generate a phylogenetic tree using FastTree tool. 
This tree along with the pan-genome matrix can then be used to plot some nice plots. Roary comes with a python script called roary_plots.py that can be used for visualizing pan-genome analysis results. 
```
module load fasttree
FastTree core_gene_alignment.aln > core_gene_alignment.tree
module load python-anaconda3/latest
python /scratch/micro612w19_fluxod/shared/bin/roary/contrib/roary_plots/roary_plots.py core_gene_alignment.tree gene_presence_absence.csv
```
-->

> ***ii. Explore pan-genome matrix gene_presence_absence.csv and gene_presence_absence.Rtab using R***

<!---
Note:plots generated by roary_plots.py doesn't seem to be very useful and is completely a waste of time. query_pan_genome script provided by roary doesn't work and generates an empty result which seems like a bug and the link to the issue raised for this bug can be found [here](https://github.com/sanger-pathogens/Roary/issues/298) 
-->

**Modify gene_presence_absence.Rtab file to include annotations**

- Get column names from gene_presence_absence.csv file

```
head -n1 gene_presence_absence.csv | tr ',' '\n' | cat --number
```
- Pull columns of interest

```
cut -d "," -f 3 gene_presence_absence.csv | tr '"' '_' > gene_presence_absence_annot.csv
```

- Paste it into pan-genome matrix

```
paste -d "" gene_presence_absence_annot.csv gene_presence_absence.Rtab > gene_presence_absence_wannot.Rtab
```

- Check gene_presence_absence_wannot.Rtab file

```
less gene_presence_absence_wannot.Rtab
```

**Read matrix into R, generate exploratory plots and query pan-genome**

Use scp or cyberduck to get gene_presence_absence_wannot.Rtab onto your laptop desktop folder.

> ***i. Prepare and clean data***

- Fire up RStudio and read gene_presence_absence_wannot.Rtab into matrix.

```
pg_matrix = read.table('~/Desktop/gene_presence_absence_wannot.Rtab', sep = "\t", quote = "", row.names = 1, skip = 1)
```

- Add column names back

```
colnames(pg_matrix) = c('ACICU', 'AbauA', 'AbauB', 'AbauC')
```

- Use head, str, dim, etc. to explore the matrix.

> ***ii. Generate exploratory heatmaps.***

- Make a heatmap for the full matrix

```
heatmap(as.matrix(pg_matrix), , scale = "none", distfun = function(x){dist(x, method = "manhattan")}, margin = c(10,10), cexCol = 0.85, cexRow = 0.5, col= c('black', 'red'))
```

- Make a heatmap for variable genes (present in at least one, but not all of the genomes)

```

pg_matrix_subset = pg_matrix[rowSums(pg_matrix > 0) > 0 & rowSums(pg_matrix > 0) < 4 ,] 
heatmap(as.matrix(pg_matrix_subset), , scale = "none", distfun = function(x){dist(x, method = "manhattan")}, margin = c(10,10), cexCol = 0.85, cexRow = 0.5, col= c('black', 'red'))

```

> ***iii. Query pan-genome***

-  Which genomes are most closely related based upon shared gene content?

We will use the outer function to determine the number of genes shared by each pair of genomes. 

<!--
Here we are arbitrarily deciding that a gene is present if the BSR is greater than 0.4. 
-->

Look at the help page for outer to gain additional insight into how this is working.

```
help(outer)
```

```
outer(1:4,1:4, FUN = Vectorize(function(x,y){sum(pg_matrix_subset[,x] > 0 & pg_matrix_subset[,y] > 0)}))
```

- What is the size of the core genome?

Lets first get an overview of how many genes are present in different numbers of genomes (0, 1, 2, 3 or 4) by plotting a histogram. Here, we combine hist with rowSums to accomplish this.

```
hist(rowSums(pg_matrix > 0), col="red")
```

Next, lets figure out how big the core genome is (e.g. how many genes are common to all of our genomes)?

```
sum(rowSums(pg_matrix > 0) == 4)
```

- What is the size of the accessory genome?

Lets use a similar approach to determine the size of the accessory genome (e.g. those genes present in only a subset of our genomes).

```
sum(rowSums(pg_matrix > 0) < 4 & rowSums(pg_matrix > 0) > 0)
```

- What types of genes are unique to a given genome?

So far we have quantified the core and accessory genome, now lets see if we can get an idea of what types of genes are core vs. accessory. Lets start by looking at those genes present in only a single genome. 

```
row.names(pg_matrix[rowSums(pg_matrix > 0) == 1,])
```

What do you notice about these genes?

- What is the number of hypothetical genes in core vs. accessory genome?

Looking at unique genes we see that many are annotated as “hypothetical”, indicating that the sequence looks like a gene, but has no detectable homology with a functionally characterized gene. 

Determine the fraction of “hypothetical” genes in unique vs. core. 

```
sum(grepl("hypothetical" , row.names(pg_matrix[rowSums(pg_matrix > 0) == 1,]))) / sum(rowSums(pg_matrix > 0) == 1)
sum(grepl("hypothetical" , row.names(pg_matrix[rowSums(pg_matrix > 0) == 4,]))) / sum(rowSums(pg_matrix > 0) == 4)
```

Why does this make sense?

Perform genome comparisons with [ACT](http://www.sanger.ac.uk/science/tools/artemis-comparison-tool-act)
-------------------------------------
[[back to top]](https://github.com/alipirani88/Comparative_Genomics/blob/master/day2_afternoon/README.md)
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

In the previous exercises we were focusing on gene content, but losing the context of the structural variation underlying gene content variation (e.g. large insertions and deletions). 
Here we will use ACT to compare two of our genomes (note that you can use ACT to compare more than two genomes if desired). 

> ***i. Create ACT alignment file with BLAST***

As we saw this morning, to compare genomes in ACT we need to use BLAST to create the alignments. We will do this on flux.

```

cd scratch/micro612w19_fluxod/username/day2_after
blastall -p blastn -i ./Abau_genomes/AbauA_genome.fasta -d ./Abau_BLAST_DB/ACICU_genome.fasta -m 8 -e 1e-20 -o AbauA_vs_ACICU.blast

```

> ***ii. Read in genomes, alignments and annotation files***

Use scp or cyberduck to transfer Abau_ACT_files folder onto your laptop


1. Abau_genomes/AbauA_genome.fasta 
2. Abau_genomes/ACICU_genome.fasta 
3. AbauA_vs_ACICU.blast 
4. Abau_ACT_files/AbauA_genome_gene.gff 
5. Abau_ACT_files/ACICU_genome_gene.gff


> ***iii. Explore genome comparison and features of ACT***

Read in genomes and alignment into ACT

```

Go to File -> open 
Sequence file 1  = ACICU_genome.fasta 
Comparison file 1  = AbauA_vs_ACICU.blast
Sequence file 2  = AbauA_genome.fasta

```

Before we use annotation present in genbank files. Here we will use ACT specific annotation files so we get some prettier display (resistance genes = red, transposable elements = bright green)  

```

Go to File -> ACICU_genome.fasta -> Read an entry file = ACICU_genome_gene.gff

Go to File -> AbauA_genome.fasta -> Read an entry file = AbauA_genome_gene.gff

```

Play around in ACT to gain some insight into the sorts of genes present in large insertion/deletion regions. 
See if you can find:

1) differences in phage content, 
2) membrane biosynthetic gene cluster variation and 
3) antibiotic resistance island variation.

