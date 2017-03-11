# Day 2 Morning
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

On day 1 we worked through a pipeline to map short-read data to a pre-existing assembly and identify single-nucleotide variants (SNVs) and small insertions/deletions. However, what this sort of analysis misses is the existence of sequence that is not present in your reference. Today we will tackle this issue by assembling our short reads into larger sequences, which we will then analyze to characterize the functions unique to our sequenced genome.   

Execute the following command to copy files for this morning’s exercises to your scratch directory: 

```
> Note: Make sure you change 'username' in the commands below to your 'uniqname'. 

cd /scratch/micro612w17_fluxod/username

> Note: Check if you are in your home directory(/scratch/micro612w17_fluxod/username) by executing 'pwd' in terminal. 'pwd' stands for present working directory and it will display the directory you are in.

pwd

> Note: Copy files for this morning's exercise in your home directory.

cp -r /scratch/micro612w17_fluxod/shared/data/day2_morn ./
```

## Genome Assembly using [Spades](http://bioinf.spbau.ru/spades) Pipeline
[[back to top]](https://github.com/alipirani88/Comparative_Genomics/blob/master/day2_morning/README.md)
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

![alt tag](https://github.com/alipirani88/Comparative_Genomics/blob/master/_img/day2_morning/intro.png)

There are a wide range of tools available for assembly of microbial genomes. These assemblers fall in to two general algorithmic categories, which you can learn more about [here](?). In the end, most assemblers will perform well on microbial genomes, unless there is unusually high GC-content or an over-abundance of repetitive sequences, both of which make accurate assembly difficult. 

Here we will use the Spades assembler with default parameters. Because genome assembly is a computationally intensive process, we will submit our assembly jobs to the cluster, and move ahead with some pre-assembled genomes, while your assemblies are running. 

>i. Create directory to hold your assembly output.

Create a new directory for the spades output in your day2_morn folder

```
> Note: Make sure you change 'username' in the below command with your 'uniqname'. 

cd /scratch/micro612w17_fluxod/username/day2_morn

> We will create a new directory in day2_morn to save genome assembly results:

mkdir Rush_KPC_266_assembly_result

```

Now, we will use a genome assembly tool called Spades for assembling the reads.

>ii. Test out Spades to make sure its in your path

To make sure that your paths are set up correctly, try running Spades with the –h (help) flag, which should produce usage instruction.

```
> Run the below commands to load a python module and check if spades is working. 

module load python
spades.py -h     

```

>iii. Submit a cluster job to assemble 

Since it takes huge amount of memory and time to assemble genomes using spades, we will run a pbs script on cluster for this step.

Now, Open the spades.pbs file residing in day2_morning folder with nano and add the following spades command to the bottom of the file. 

```
> Open spades.pbs file using nano:

nano spades.pbs

> Now replace the EMAIL_ADDRESS in spades.pbs file with your actual email-address. This will make sure that whenever the job starts, aborts or ends, you will get an email notification.

> Copy and paste the below command to the bottom of spades.pbs file.

spades.py --pe1-1 forward_paired.fq.gz --pe1-2 reverse_paired.fq.gz --pe1-s forward_unpaired.fq.gz --pe1-s reverse_unpaired.fq.gz -o Rush_KPC_266_assembly_result/ --careful

```

>iv. Submit your job to the cluster with qsub

```
qsub -V spades.pbs
```

>v. Verify that your job is in the queue with the qstat command

```
qstat –u username 
```

## Assembly evaluation using [QUAST](http://bioinf.spbau.ru/quast)
[[back to top]](https://github.com/alipirani88/Comparative_Genomics/blob/master/day2_morning/README.md)
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

The output of an assembler is a set of contigs (contiguous sequences), that are composed of the short reads that we fed in. Once we have an assembly we want to evaluate how good it is. This is somewhat qualitative, but there are some standard metrics that people use to quantify the quality of their assembly. Useful metrics include: i) number of contigs (the fewer the better), ii) N50 (the minimum contig size that at least 50% of your assembly belongs, the bigger the better). In general you want your assembly to be less than 200 contigs and have an N50 greater than 50 Kb, although these numbers of highly dependent on the properties of the assembled genome. 

To evaluate some example assemblies we will use the tool quast. Quast produces a series of metrics describing the quality of your genome assemblies. 

>i. Run quast on a set of previously generated assemblies

Now to check the example assemblies residing in your day2_morn folder, run the below quast command. Make sure you are in day2_morn folder in your home directory using 'pwd'

```
quast.py -o quast sample_264_contigs.fasta sample_266_contigs.fasta
```

The command above will generate a report file in /scratch/micro612w17_fluxod/username/day2_morn/quast

>ii. Explore quast output

QUAST creates output in different formats such as html, pdf and text. Now lets check the report.txt file residing in quast folder for assembly statistics. Open report.txt using nano.

```
nano quast/report.txt
```
Check the difference between each assembly statistics. Also check different type of report formats it generated.

## Compare assembly to reference genome and post-assembly genome improvement
[[back to top]](https://github.com/alipirani88/Comparative_Genomics/blob/master/day2_morning/README.md)
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

Now that we feel confident in our assembly, lets compare it to our reference to see if we can identify any large insertions/deletions using a graphical user interface called Artemis Comparison Tool (ACT) for visualization. 

To do this we need to first align our genome assembly to our reference. We will accomplish this using command-line BLAST.

>i. Align unordered contigs to reference

Create a BLAST database from your reference genome using the makeblastdb command.

```

> Make sure you are in /scratch/micro612w17_fluxod/username/day2_morn directory
cd /scratch/micro612w17_fluxod/username/day2_morn

makeblastdb -in KPNIH1.fasta -dbtype nucl -out KPNIH1.fasta

```

>ii. Stitch together your contigs into a single sequence

```
echo ">sample_266_contigs_concat" > sample_266_contigs_concat.fasta 
grep -v ">" sample_266_contigs.fasta >> sample_266_contigs_concat.fasta 
```

BLAST your stitched together contigs against your reference. 
The input parameters are: 
1) query sequences (-query sample_266_contigs_concat.fasta), 
2) the database to search against (-db KPNIH1.fasta), 
3) the name of a file to store your results (-out blastn_results), 
4) output format (-outfmt 6), 
5) e-value cutoff (-evalue 1e-20)

```
blastn -outfmt 6 -evalue 1e-20 -db KPNIH1.fasta -query sample_266_contigs_concat.fasta -out concat_comp.blast
```

>ii. Use ACT(Installed in your local system) to compare stitched together contigs to reference.

For these, first we will create a seperate directory called ACT_contig_comparison in day2_morn folder and copy all the necessary ACT input to this directory.

```

mkdir ACT_contig_comparison 
cp KPNIH.gb KPNIH1.fasta concat_comp.blast sample_266_contigs_concat.fasta ACT_contig_comparison/

```

Use scp to get sequences and BLAST alignments onto your laptop 

```

> Note: Make sure you change 'username' in the below command with your 'uniqname'.

scp -r username@flux-xfer.arc-ts.umich.edu:/scratch/micro612w17_fluxod/username/day2_morn/ACT_contig_comparison/ /path-to-local-directory/

```

>iii. Read these Input files in ACT_contig_comparison folder into ACT

```

Start ACT and set your working directory to ACT_contig_comparison(wherever it is stored on your local system)
Go to File on top left corner of ACT window -> open 

Sequence file 1 = KPNIH.gb
Comparison file 1  = concat_comp_blast 
Sequence file 2  = sample_266_contigs_concat.fasta

Click Apply button
```

> Notice that it is a complete mess!!!! The reason is that the contigs are in random order, so it is very difficult to visually compare to the reference. 

![alt tag](https://github.com/alipirani88/Comparative_Genomics/blob/master/_img/day2_morning/mess.png)

iv. Run abacas to orient contigs to reference

To orient our contigs relative to the reference we will use a tool called abacas. [ABACAS](http://www.sanger.ac.uk/science/tools/pagit) aligns contigs to a reference genome and then stitches them together to form a “pseudo-chromosome”. 

Go back to flux and into the directory where the assembly is located.

```
cd /scratch/micro612w17_fluxod/username/day2_morn/
```

Before running Abacas, add the following path to your ~/.bashrc file:

```
> Open bashrc file using nano

nano ~/.bashrc

> Add this path to the bottom of bashrc file

export PATH=$PATH:/scratch/micro612w17_fluxod/shared/bin/MUMmer3.23/

>  Exit and save this file. Source your bashrc file.

source ~/.bashrc

```

Now, we will run abacas using these input parameters: 

1) your reference sequence (-r KPNIH.fasta), 
2) your contig file (-q sample_266_contigs.fasta), 
3) the program to use to align contigs to reference (-p nucmer), 
4) append unmapped contigs to end of file (-b), 
5) use default nucmer parameters (-d), 
6) append contigs into pseudo-chromosome (-a), 
7) the prefix for your output files (–o sample_266_contigs_ordered) 

```
perl abacas.1.3.1.pl -r KPNIH1.fasta -q sample_266_contigs.fasta -p nucmer -b -d -a -o sample_266_contigs_ordered
```

v. Use ACT to view contig alignment to reference genome

> Use scp to get ordered fasta sequence and .cruch file onto your laptop 

```
> Dont forget to change username and /path-to-local-ACT_contig_comparison-directory/ in the below command

scp username@flux-xfer.arc-ts.umich.edu:/scratch/micro612w17_fluxod/username/day2_morn/sample_266_contigs_ordered* /path-to-local-ACT_contig_comparison-directory/

```

> Read files into ACT

```
Go to File on top left corner of ACT window -> open 
Sequence file 1 = KPNIH.gb 
Comparison file 1  = sample_266_contigs_ordered.crunch 
Sequence file 2  = sample_266_contigs_ordered.fasta

Click Apply button

> Dont close the ACT window
```

> Notice that the alignment is totally beautiful now!!! Scan through the alignment and play with ACT features to look at genes present in reference but not in assembly. Keep the ACT window open for further visualizations.

![alt tag](https://github.com/alipirani88/Comparative_Genomics/blob/master/_img/day2_morning/beautiful.png)

## Map reads to the final ordered assembly
[[back to top]](https://github.com/alipirani88/Comparative_Genomics/blob/master/day2_morning/README.md)
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

You already know the drill/steps involved in reads mapping. Here, we will map the reads to the final ordered assembly genome instead of KPNIH1.fasta.

First create bwa index of ordered fasta file.

```
> Only proceed further if everything worked uptil now. Make sure you are in day2_morn directory.

cd /scratch/micro612w17_fluxod/username/day2_morn/
bwa index sample_266_contigs_ordered.fasta
samtools faidx sample_266_contigs_ordered.fasta

```

Align the trimmed reads which we used for genome assembly to this ordered assembly using BWA mem. Convert SAM to BAM. Sort and index it.

```

bwa mem -M -R "@RG\tID:96\tSM:Rush_KPC_266_1_combine.fastq.gz\tLB:1\tPL:Illumina" -t 8 sample_266_contigs_ordered.fasta forward_paired.fq.gz reverse_paired.fq.gz > sample_266_contigs_ordered.sam

samtools view -Sb sample_266_contigs_ordered.sam > sample_266_contigs_ordered.bam

samtools sort sample_266_contigs_ordered.bam sample_266_contigs_ordered_sort

samtools index sample_266_contigs_ordered_sort.bam

```

Lets visualize the alignments against our ordered assembly.
Copy this sorted and indexed BAM files to local ACT_contig_comparison directory.

```
> Dont forget to change username and /path-to-local-ACT_contig_comparison-directory/ in the below command

scp username@flux-xfer.arc-ts.umich.edu:/scratch/micro612w17_fluxod/username/day2_morn/sample_266_contigs_ordered_sort* /path-to-local-ACT_contig_comparison-directory/

```

```
Go back to ACT where your ordered contigs are still open in the window.

Select File -> sample_266_contigs_ordered.fasta -> Read BAM/VCF > select sorted bam file(sample_266_contigs_ordered_sort.bam) you just copied from flux.
```

![alt tag](https://github.com/alipirani88/Comparative_Genomics/blob/master/_img/day2_morning/aligned_reads_deletion.png)


## Genome Annotation
[[back to top]](https://github.com/alipirani88/Comparative_Genomics/blob/master/day2_morning/README.md)
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

**Identify protein-coding genes with [Prokka](http://www.vicbioinformatics.com/software.prokka.shtml)**

From our ACT comparison of our assembly and the reference we can clearly see that there is unique sequence in our assembly. However, we still don’t know what that sequence encodes! To try to get some insight into the sorts of genes unique to our assembly we will run a genome annotation pipeline called Prokka. Prokka works by first running denovo gene prediction algorithms to identify protein coding genes and tRNA genes. Next, for protein coding genes Prokka runs a series of comparisons against databases of annotated genes to generate putative annotations for your genome. 

>i. Run Prokka on assembly

Load modules required for Prokka

```
module load med perl-modules prokka 
prokka –setupdb
```

Execute Prokka on your ordered assembly 

```
> Make sure you are in day2_morn directory.

cd /scratch/micro612w17_fluxod/username/day2_morn/
mkdir sample_266_prokka 

> Dont forget to change username in the below command

prokka -kingdom Bacteria -outdir /scratch/micro612w17_fluxod/username/day2_morn/sample_266_prokka -force -prefix sample_266 sample_266_contigs_ordered.fasta

> Use scp to get Prokka annotated genome on your laptop. 

scp -r username@flux-xfer.arc-ts.umich.edu:/scratch/micro612w17_fluxod/username/day2_morn/sample_266_prokka/ /path-to-local-ACT_contig_comparison-directory/

```

>ii. Reload comparison into ACT now that we’ve annotated the un-annotated!

Read files into ACT

```
Go to File on top left corner of ACT window -> open 
Sequence file 1  = KPNIH.gb 
Comparison file 1  = sample_266_contigs_ordered.crunch 
Sequence file 2  = sample_266_contigs_ordered.gbf
```

>Play around with ACT to see what types of genes are unique to sample 266!!! 
