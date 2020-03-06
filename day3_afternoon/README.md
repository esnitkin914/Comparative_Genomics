Day 3 Afternoon
===============
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

Klebsiella pneumoniae comparative genomic analysis 
--------------------------------------------------

To finish up the workshop we are going to go through the process of working up a complete dataset, from start to finish.  This set of genomes originated from a regional outbreak of bla-KPC carrying Klebsiella pneumoniae – one of the most concerning healthcare associated pathogens. 
The goal is to follow up on a previously [published](http://cid.oxfordjournals.org/content/53/6/532.abstract) epidemiologic analysis, and see if genomics supports prior epidemiologic conclusions and can provide additional insights. 

The results of this genomics analysis were published in [this](http://stm.sciencemag.org/content/9/417/eaan0093) paper.

We have our genomes, and we know in which regional facility each isolate originated. 

The goal of this exercise is to:

1) process our genomes (QC, variant calling), 

2) perform a phylogenetic analysis and 

3) overlay our meta-data. 

To make this more difficult, the instructions will be much more vague than in previous sessions, and you will be challenged to use what you have learned, both in the past three days and in the prior workshop, to complete this analysis. 

Hopefully we’ve prepared you to take on the challenge, but remember this is an open book test! 

Feel free to lean on materials from the workshops, manuals of tools and Google (and of course instructors and neighbors). 

Execute the following command to copy files for this afternoon’s exercises to your scratch directory:

```

cd /scratch/micro612w20_class_root/micro612w20_class/username

or 

wd

cp -r /scratch/micro612w20_class_root/micro612w20_class/shared/data/day3pm ./

```

Perform QC on fastq files
-------------------------
[[back to top]](https://github.com/alipirani88/Comparative_Genomics/blob/master/day3pmnoon/README.md)
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

On the first morning you ran FastQC to evaluate the quality of a single genome. However, a typical project will include many genomes and you will want to check the quality of all of your samples. From the bash workshop, I hope you can appreciate that you do not want to process 100 genomes by typing 100 commands – rather you want to write a short shell script to do the work for you!


> ***i. Edit the shell script fastqc.sh located in /scratch/micro612w20_class_root/micro612w20_class/your username/day3pm to run FastQC on all fastq files.***

**Important info about this shell script** 
- The shell script includes a for loop that loops over all of the genomes in the target directory
- The tricky part of this shell script is that each fastq command contains two files (forward and reverse reads). So, you need to take advantage of the fact that the forward and reverse read files both have the same prefix, and you can loop over these prefixes. 
- You should be able to get prefixes by piping the following unix commands: ls, cut, sort, uniq
- The prefix should be a part of both forward and reverse reads. For example, the file_prefix for samples Rush_KPC_264_1_sequence.fastq.gz and Rush_KPC_264_2_sequence.fastq.gz should be Rush_KPC_264
- when you are testing your shell script, comment out (using #) the lines below echo so you can see that if the script is 'echo'-ing the correct commands.
- Try running multiqc inside the script by adding the multiqc command with appropriate out directory 
- Don't run multiqc inside for loop and should be run only after the for loop ends.


The fastq files are located in:

```
/scratch/micro612w20_class_root/micro612w20_class/shared/data/day3pm_fastq/
```

Rather than copying these to your directory, analyze the files directly in that directory, so everyone doesn’t have to copy 25G to their home directories. 

Copy and paste commands to run fastqc.sh as slurm script, into a slurm script and submit this slurm script as a job to great lakes.

Your slurm script wil contain the following command after the slurm preamble stuff(Make sure your $SLURM_SUBMIT_DIR is set inside the slurm script):

```bash fastqc.sh /scratch/micro612w20_class_root/micro612w20_class/shared/data/day3pm_fastq/ ```


> ***ii. Examine output of FastQC to verify that all samples are OK***

Check the multiqc report of your fastq files.

Run ARIBA CARD and MLST typing on day3pm_fastq samples
------------------------------------------------------
[[back to top]](https://github.com/alipirani88/Comparative_Genomics/blob/master/day3pmnoon/README.md)
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

On Day 2 afternoon, you ran ARIBA to perform sequence typing of a single genome using ARIBA tool. However, a typical public epidimiology project will include many genomes sampled from different sites and you will want to sequence type each of these samples to study their genetic diversity and determine if these samples contain antibiotic resistance genes. 

Since both MLST typing and CARD resistance gene detection takes a while to run on these samples, We already processed them and have placed the results in day3pm folder.

> i. Explore MLST typing results. 

Try exploring mlst.sh and mlst.sbat scripts that were used to generate these results and try to understand the implementation of for loop in running them sequentially.

You can find the processed MLST results for day3pm_fastq samples in this folder:

/scratch/micro612w20_class_root/micro612w20_class/shared/data/day3pm/MLST_results

Run summarize_mlst.sh and check the dominant sequence type that these samples belong to.

> ii. Summarize ARIBA CARD results

Try exploring card.sh and card.sbat scripts that were used to generate these results and try to understand the implementation of for loop in running them sequentially.

Use Ariba's summary function and summarize ARIBA CARD results that are placed in:

/scratch/micro612w20_class_root/micro612w20_class/shared/data/day3pm/CARD_results

> iii. Drag and drop Phandango associated files onto the Phandango website and check what type of actibiotic resistance genes are present in these genomes?

> iv. Explore full ARIBA matrix in R and plot a heatmap to visualize the presence/absence of various resistancve genes.

how many of these samples contain carbapenam resistant gene - KPC?


Examine results of [SPANDx](http://www.ncbi.nlm.nih.gov/pubmed/25201145) pipeline
---------------------------
[[back to top]](https://github.com/alipirani88/Comparative_Genomics/blob/master/day3pmnoon/README.md)
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

On the afternoon of day 1 we saw how many steps are involved in calling variants relative to a reference genome. However, the same steps are applied to every sample, which makes this very pipeline friendly!  So, you could write your own shell script to string together these commands, or take advantage of one of several published pipelines. Here, we will use the output of the SPANDx pipeline, which takes as input a directory of fastq files and produces core variant and indel calls.

More information on SPANDx pipeline can be obtained from [this](https://sourceforge.net/projects/spandx/files/SPANDx%20Manual_v3.1.pdf/download) manual.

A snapshot of the pipeline is shown below:

![alt tag](https://github.com/alipirani88/Comparative_Genomics/blob/master/_img/spandx.jpg)

Because it takes a while to run, we have pre-run it for you. Your task will be to sort through the outputs of SPANDx. The detailed information about how to interpret the output is in SPANDx manual(section INTERPRETING THE OUTPUTS). 

> ***i. Look at overall statistics for variant calling in excel***

SPANDx produces an overall summary file of its run that includes:

1) numbers of SNPs/indels,

2) numbers of filtered SNPs/indels and 

3) average coverage across the reference genome. 

This summary file is in:  Outputs/Single_sample_summary.txt

Use less to look at this file and then apply unix commands to extract and sort individual columns 

**HINTS**
The following unix commands can be used to get sorted lists of coverage and numbers of SNPs/indels: tail, cut, sort

> ***ii. Look at filtered variants produced by SPANDx in excel***

SPANDx also produces a summary file of the variants/indels it identified in the core genome. 

This summary file is: 
```/scratch/micro612w20_class_root/micro612w20_class/username/day3pm/SPANDx_output/Outputs/All_SNPs_annotated.txt ```

Use cyberduck/scp to download this file and view in excel

- View SPANDx manual for interpretation of different columns which can be found [here](https://sourceforge.net/projects/spandx/files/SPANDx%20Manual_v3.1.pdf/download)
- Back on great lakes, use grep to pull SNPs that have HIGH impact
- What types of mutations are predicted to have “HIGH” impact?
- How many genomes do these HIGH impact mutations tend to be present in? How do you interpret this?

Recombination detection and tree generation
-------------------------------------------
[[back to top]](https://github.com/alipirani88/Comparative_Genomics/blob/master/day3pmnoon/README.md)
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

> ***i. Plot the distribution of variants across the genome in R***

The positions of variants are embedded in the first column of Outputs/Comparative/All_SNPs_annotated.txt, but you have to do some work to isolate them! 

**HINTS**  

- You will need to pipe together two “cut” commands: the first command will use tab as a delimiter and the second will use _. 
- Note that for cut you can specify tab as the delimiter as follows: cut –d$’\t’ and _ as: cut -d ‘_’
- You should redirect the output of your cut commands (a list of SNP positions) to a file called ‘snp_positions.txt’. For example, the first line of your snp_positions.txt should be:
```
12695
```
- Finally, download this file, read it into R using ‘read.table’ and use ‘hist’ to plot a histogram of the positions
- Do you observe clustering of variants that would be indicative of recombination?

> ***ii.  Create fasta file of variants from nexus file***

SPANDx creates a file of core SNPs in a slightly odd format (transposed nexus). 
This file is called: 
```/scratch/micro612w20_class_root/micro612w20_class/username/day3pm/SPANDx_output/Outputs/Comparative/Ortho_SNP_matrix.nex ```

For convenience, apply the custom perl script located in the same directory to convert it to fasta format

```
perl transpose_nex_to_fasta.pl Ortho_SNP_matrix.nex
```

This file Outputs/Comparative/Ortho_SNP_matrix.fasta should now exist

> ***iii. Create a neighboring-joining tree in R***

```

Download Ortho_SNP_matrix.fasta to your home computer
Read the fasta file into R, create a distance matrix, and make a neighbor joining tree. Make sure you load in the necessary packages! 
Save tree for later analysis using write.tree()

```

Phylogenetic tree annotation and visualization
----------------------------------------------
[[back to top]](https://github.com/alipirani88/Comparative_Genomics/blob/master/day3pmnoon/README.md)
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)
<!---

> ***i. Load the neighbor-joining tree into iTOL***

Instead of annotating your tree in R, this time, let's use iTOL to make a publication quality tree. iTOL is a web-based tool to visualize trees and metadata of the taxa. 

Note that because the out-group is so distantly related it is difficult to make out the structure of the rest of the tree. 

**To remedy this:**

- Click on the KPNIH1 leaf, go to the “tree structure” menu and “delete leaf” 
- Click on the extended branch leading to where KPNIH1 was, go to the “tree structure” menu and click “collapse branch”

> ***ii. Load the annotation file ‘Rush_KPC_facility_codes_iTOL.txt’ to view the facility of isolation, play with tree visualization properties to understand how isolates group by facility, Circular vs. normal tree layout, Bootstrap values, Ignoring branch lengths***
-->

```

Which facilities appear to have a lot of intra-facility transmission based on grouping of isolates from the same facility? 
Which patient’s infections might have originated from the blue facility?

```

<!---
Assessment of genomic deletions
-------------------------------
[[back to top]](https://github.com/alipirani88/Comparative_Genomics/blob/master/day3pmnoon/README.md)
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

> ***i. Download genome coverage bed file and load into R***

This file is located in: Outputs/Comparative/Bedcov_merge.txt
This file contains information regarding locations in the reference genome that each sequenced genome does and does not map to.

The first 3 columns of the file are:

1) the name of the reference, 

2) the start coordinate of the window and 

3) the end coordinate of the window

The remaining columns are your analyzed genomes, with the values indicating the fraction of the window covered by reads in that genome.

In essence, this file contains information on parts of the reference genome that might have been deleted in one of our sequenced genomes.

After you download this file, read it into R

**HINTS**
- Use the read.table function with the relevant parameters being: header and sep

> ***ii. Plot heatmap of genome coverage bed file***

**HINTS**

- The first 3 columns of the bed file specify the name of the chromosome and the genome coordinates – therefore you want to subset your matrix to not include these columns 
- Use the heatmap3 function to make your heatmap with the following parameters: scale = “none” (keeps original values), Rowv = NA (suppress clustering by rows – why might we not want to cluster by rows for this analysis?)

- Note a large genomic deletion among a subset of isolates. Does this deletion fit with the phylogeny from above?

iii. Explore genomic deletion in more detail with ACT

- Use abacus to orient contigs from Rush_KPC_298 to KPNIH 
- Load KPNIH.gb, Rush_KPC_298_ordered and the .crunch alignment into ACT

```

What genes appear to have been lost?

```
-->
