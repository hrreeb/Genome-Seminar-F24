---
title: "Gouldian Finch Genome Seminar"
author: "Hannah Reeb"
date: "2024-12-09"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
### Project Background 
Gouldian finches are brightly colored birds, with some morphs having yellow plumage on both their belly and head feathers. In those morphs, the yellow color does not come from the same carotenoid pigments in both patches. In the belly feathers, lutein is the source of the yellow, and in the mask, 3'-dehydrolutein is the source. In House Finches, DHRS13 is involved in converting lutein to dehydrolutein, and we have som initial evidence that it is differentially expressed in head versus belly feathers of Gouldian finches. Here, our goal is to assess whether DHRS13 is upregulated in feathers where dehydrolutein is present, implicating its importance in carotenoid conversion in Gouldian finches. We are also interested in finding other differentially expressed genes between the two feather types, and whether any are genes known to be associated with carotenoid metabolism. 

![Gerhard Hofmann and Dr. Claudia Mettke-Hofmann.](/Users/loreeb/Downloads/gouldian finch work/to push to GIt/gouldianmorphs.png)







We predict that DHRS13 will be upregulated in head feathers. To quantify this, we take RNA sequences from Gouldian finch head and belly feathers, and map them to a reference genome to find upregulated genes related to carotenoid expression. 

### Data Sources
This analysis uses 8 samples of single-end, Illumina reads gathered from Gouldian Finch (GOFI) feathers. 4 samples from the belly, 4 from the head. This data was collected in 2018 by the Corbo Lab at WashU. 

Throughout the project I number the samples as below:

#0 = belly

#1 = head

#2 = belly

#3 = head

#4 = belly

#5 = belly

#6 = head

#9 = head

### Methods Overview
1. Fast QC
2. Trim reads
3. Align and count with kallisto
4. Utilize sleuth in R to find DEGs
5. Visualize DEGs using volcano plots and heatmaps
6. Align reads to genome using Hisat2
7. Get BAM file for manipulation
8. Load BAM file to IGV to visualize splicing

## Analysis and Results

## Step 1: FastQC
We will start with an initial survey of read quality using fastqc. I used an array job for this step. Overall, the reads look pretty good-- adapter content was low before trimming, and quality scores across bases was high. I included one link to the html output from fastQC as an example, and a snapshot of quality score-- all samples were similar in quality. 

[fastqc.sbatch](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/fastq_array.sbatch)
[fastqc.args](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/fastqc_array.args)
[fastqc.sh](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/fastqc_array.sh)

[example of output](/Users/loreeb/Downloads/gouldian finch work/FastQC output/run_2609_s_1_withindex_sequence.txt_CCGATTA_fastqc.html)

![](/Users/loreeb/Downloads/gouldian finch work/to push to GIt/fastQCexample.png)



## Step 2: Trimming 
We use trimgalore to trim reads. I initially tried an array job for this step, but there were some issues with opening files, which I only discovered after running each job individually. Since that worked best, those are the files I link here. I trim reads to be 100 bp or longer. 


[trim0.sbatch](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/trimming_gould.sbatch)
[trim0.sh](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/trimming_gould.sh)


[trim1.sbatch](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/trimming_gould1.sbatch)
[trim1.sh](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/trimming_gould1.sh)


[trim2.sbatch](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/trimming_gould2.sbatch)
[trim2.sh](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/trimming_gould2.sh)


[trim3.sbatch](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/trimming_gould3.sbatch)
[trim3.sh](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/trimming_gould3.sh)


[trim4.sbatch](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/trimming_gould4.sbatch)
[trim4.sh](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/trimming_gould4.sh)


[trim5.sbatch](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/trimming_gould5.sbatch)
[trim5.sh](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/trimming_gould5.sh)



[trim6.sbatch](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/trimming_gould6.sbatch)
[trim6.sh](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/trimming_gould6.sh)


[trim9.sbatch](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/trimming_gould9.sbatch)
[trim9.sh](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/trimming_gould9.sh)


After trimming, adapter content is decreased. We can now have more confidence that our alignments will not be contaminated with adapter sequences!


![](/Users/loreeb/Downloads/gouldian finch work/to push to GIt/trimBa.png)


## Step 3: Pseudoalignment with Kallisto
Now that I have trimmed reads, we want to align them to known reference to identify genes of interest. The first step is to create a k-mer indexed reference. I used the gouldian finch cDNA from ensembl.

[index_kallisto.sbatch](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/push2/index_kallisto.sbatch)
[index_kallisto.sh](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/push2/index_kallisto.sh)

The next step is to actually align our reads to the indexed reference. This is pseudoalignment. I do this task as an array job. Since we have single end reads, kallisto will not estimate fragment lenght and standard deviation for us. Since we have Illumina reads, I follow the Pachter Lab's recommendation on settings. I use average fragment length = 200, standard deviation of fragment length = 20. Kallisto gives us abundance.tsv and abundance.h5 files as output, which will let us understand up or down regulation of genes later.
 

[Pachter lab kallisto documentation](https://pachterlab.github.io/kallisto/manual)


[kallisto_quant.args](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/push2/kallisto_quant.args)
[kallisto_quant.sbatch](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/push2/kallisto_quant.sbatch)
[kallisto_quant.sh](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/push2/kallisto_quant.sh)


## Step 4: Sleuthing for DEGs
The sleuth package in R is useful for understanding kallisto output. First, though, we need to pull useful information from the cDNA file we use. I use a regular expression to pull out the Gene ID and Transcript ID, then manipulate the output in excel until I am left with columsn of Gene ID and Transcript ID to use as a header.

regex for headers
```{}
ln -s /scratch/biol726311/gouldian/Erythrura_gouldiae.GouldianFinch.cdna.all.fa Erythrura_gouldiae.GouldianFinch.cdna.all.fa

grep  '^>' Erythrura_gouldiae.GouldianFinch.cdna.all.fa | sed -E 's/>([^ ]+) .* \(([^)]+)\).*/\1 \2/' > GOULD_headers.txt
```


This lists the study sample in 2 experimental groups
[ExpTable_GOULD.txt](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/push2/ExpTable_GOULD.txt)


Now that we have organized information and our experimental groups, we can use sleuth to find differentially expressed genes. I used the code Dr. Toomey gave us in class to find them. I try it using data aggregated by gene name, as well as unaggregated. In the unaggregated, every transcript ID is counted by itself. In the aggregated, transcipt IDs with the same gene name are grouped, to give us higher statistical power. 

The volcano plot generated by these analyses uses unaggregated data, as does the heatmap, linked below. 

![](/Users/loreeb/Downloads/gouldian finch work/to push to GIt/unannotated_Volcano.png)
![](/Users/loreeb/Downloads/gouldian finch work/to push to GIt/unannotated_heatmap.png)

Oh no! the volcano plot and heat map show differential gene expression, but the transcripts are not annotated. That is not of great use to us. 

## Step 5: Uniprot search and diamond blastx 

In order to avoid redoing the last 2 steps of the analysis, I download all of the sequenced proteins for Zebra finch (uniprot), make a database, then query with Gouldian finch DNA using a diamond blastx. 


[diamond_mkdb.sh](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/push2/diamond_mkdb.sh)
[diamond_mkdb.sbatch](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/push2/diamond_mkdb.sbatch)


[diamond_blastx.sbatch](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/push2/diamond_blastx.sbatch)
[diamond_blastx.sh](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/push2/diamond_blastx.sh)

Now, I use a regex to pull headers from the output. Now transcript ID is linked to gene name, not “ENSEGOG...”

```{}
awk -F'\t' '{split($6, arr, " "); for (i in arr) if (arr[i] ~ /^GN=/) {split(arr[i], gene, "="); print $1 "\t" gene[2]}}' ZeFiGoFi_blastx.tsv > ENSEGOT_GN.tsv
```

![](/Users/loreeb/Downloads/gouldian finch work/to push to GIt/awkregexOut.png)



## Step 6: Redo sleuthing for DEGs
Now, I will repeat Step 4 and retry the sleuth work in R, this time with the zebra finch genome, which is annotated. This will give us information about genes that are up or downregulated across tissues in a way that is useful for us to search with candidate genes. 

Now, we have plots with gene names labeled! Much better!

![](/Users/loreeb/Downloads/gouldian finch work/to push to GIt/annotated_volcano.png)

![](/Users/loreeb/Downloads/gouldian finch work/to push to GIt/annotated_heatmap.png)

Here are the top 6 most differentially expressed genes!

![](/Users/loreeb/Downloads/gouldian finch work/to push to GIt/top6DEGs.png)

Most of them are genes related to keratin, which makes sense given we sampled different feather types, with different microstructures. 

## Step 7: Look at Gene Ontogeny with ShinyGO
Now that we have gene names, I can also look at gene ontogeny with ShinyGO. Here is a quick overview of some of the general functions of upregulated genes in the head feathers. 

![](/Users/loreeb/Downloads/gouldian finch work/to push to GIt/shinyGO.png)




## Step 8: Search for differentially expressed genes of interest
To start, I make sure that the gene name is what I expect it to be. I do a blastx search with the DHRS13 sequence querying Zebra finch, and it appears the name is the same in Zebra finch. Great!

![](/Users/loreeb/Downloads/gouldian finch work/to push to GIt/dhrs13namesearch.png)


From the sleuth analysis, I have a file that lists all significantly upregulated genes in the head feathers, compared to the belly. I can then search by gene names of interest! For this step, transcripts are aggregated by gene name to increase statistical power. 

I enter the following queries:

DHRS

BDH

SCARB

CYP

There is significant upregulation for the following: 

DHRS3: dehydrogenase/reductase 3

SCARB2: scavenger receptor class B member 2

CYP36B1: cytochrome P450 family 26 subfamily B member 1

CYP27A1: cytochrome P450 family 27 subfamily A member 1

CYP27C1: cytochrome P450 family 27 subfamily C member 1

CYP24C1: cytochrome P450 family 24 subfamily A member 1


![](/Users/loreeb/Downloads/gouldian finch work/to push to GIt/upregulatedNums.png)

Here is a link to the R code that supports the stats and graphs in steps 4 and 8:

[DEGs_Rcode](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/DEGSandplot.R)

## Step 9: Dig deeper on DHRS13
We have found carotenoid candidate genes that are upregulated in the head (see above), but let's explore further on our gene of highest interest, DHRS13. 

Searching through all genes present, not just those that were significantly differentially expressed, we find the following results for DHRS13. 

![](/Users/loreeb/Downloads/gouldian finch work/to push to GIt/dhrs13specific.png)


Differential expression is not significant for DHRS13 (p=0.11)
But, it is technically upregulated in the head (β=0.45). It has a positive coefficient, meaning that it is more highly expressed in the head, just not to a statistically significant degree. 





## Step 10: Align trimmed reads to GOFI genome using hisat2

Hisat2 stands for Hierarchical Indexing for Spliced Alignment of Transcripts 2. It allows us to align RNA-Seq reads to a reference genome, which is great for our purposes here (we have transcript reads for Gouldian finch and a reference genome to map to) 

The following code for establishing the hisat2 environment I got from Caleb and Salil. These steps allow us to install hisat2 within an environment using the OSCER login. Adding these channels lets us utilize the necessary packages later in our analysis. 
```{}
$ conda create -n hisat2_env
$ conda activate hisat2_env
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
$ conda config --add channels defaults
$ conda config --show channels
$ conda install hisat2

```

Now, we will use the function "hisat2-build" to create an indexed Gouldian Finch genome that is compatible with Hisat. Hisat is great for this, because unlike some other alignment tools, it takes the fact that introns exist into account during mapping. This works well for samples of RNA-Seq data being aligned to a reference genome, because transcripts will be 'missing' the introns present in the reference genome. Using hisat2 lets us eventually visualize splicing zones and introns.


[index_hisat.sh](/Users/loreeb/Downloads/gouldian finch work/to push to GIt/index_hisat.sh)
[index_hisat.sbatch](/Users/loreeb/Downloads/gouldian finch work/to push to GIt/hisat_index.sbatch)

Once we have our indexed genome, the next task is to align using hisat2! I did this as an array job. For single end reads, we align each read sample to the indexed genome individually (paired ends get submitted at once as -1 and -2.) Here, we need to specify that we are using single end reads with -U and then the file name of our trimmed reads. 

[hisat2_full_map.args](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/push2/hisat2_full_mapping.args)
[hisat2_full_map.sbatch](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/push2/hisat2_full_mapping.sbatch)
[hisat2_full_map.sh](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/push2/hisat2_full_mapping.sh)

This step was pretty quick! It took about 9 minutes to run. 
Exciting! Now we have 8 .sam files. 

### Step 11: .sam to .bam
Our next step is to convert each of these .sam files to .bam files. We will switch back to the mamba class environment for this step, and use samtools view to convert the files. Again, we'll do this step as an array job, since we have many (8) files. 

[sam_to_bam.args](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/push3/sam_to_bam.args)
[sam_to_bam.sbatch](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/push3/sam_to_bam.sbatch)
[sam_to_bam.sh](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/push3/sam_to_bam.sh)


### Step 12: Sort the .bam files
Now, we need to sort the .bam files by read name. This is what we'll need for the next couple of our processing steps

[sort_bam.sh](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/sort_bam.sh)
[sort_bam.sbatch](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/sort_bam.sbatch)
[sort_bam.args](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/sort_bam.args)

### Step 13: fixmate the name-sorted .bam files
We will tell samtools to mark the name sorted bam files with tags.

[fixmate.sh](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/fixmate.sh)
[fixmate.sbatch](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/fixmate.sbatch)
[fixmate.args](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/fixmate.args)

### Step 14: Re-sort the fixmate output
Re-sort by chromosomal position instead of name. This is necessary for next steps. 

[resort_fixmate.sh](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/resort_fixmate.sh)
[resort_fixmate.sbatch](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/resort_fixmate.sbatch)
[resort_fixmate.args](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/resort_fixmate.args)


### Step 15: Mark and remove duplicates
Tag and mark for removal any potential PCR duplicates

[markdup.args](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/markdup.args)
[markdup.sbatch](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/markdup.sbatch)
[markdup.sh](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/markdup.sh)


### Step 16: Create index with duplicate removed .bam file
Array job to add an index to our bam file. Now, we're ready for stats!

[indexed_markdup.args](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/indexed_markdup.args)
[indexed_markdup.sbatch](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/indexed_markdup.sbatch)
[indexed_markdup.sh](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/indexed_markdup.sh)


### Step 17: Mapping statistics using flagstats

[map_stats.args](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/map_stats.args)
[map_stats.sbatch](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/map_stats.sbatch)
[map_stats.sh](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/map_stats.sh)

[example of mapping stats output](https://github.com/hrreeb/Genome-Seminar-F24/blob/main/to%20push%20to%20GIt/0mapstats.txt)

```{}
8547656 + 0 in total (QC-passed reads + QC-failed reads)
752786 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
7263450 + 0 mapped (84.98% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```
For the reads from sample 1, we can see that we had 752786 total reads, and all of them passed QC! About 85% of reads mapped to the reference genome. 

## Step 18 (last one!): Use IGV to visualize splicing
In IGV, we load in our genome (Gouldian finch) and the .bam file we made in the steps above, being sure that the index is in teh same directory. Below is one example of some output from the program-- it is zoomed in on reads from the 1st sample, and shows information about scaffold 1. 

![](/Users/loreeb/Downloads/gouldian finch work/to push to GIt/scaffold1read1.png)



## Conclusions, Future Directions, and Reflection

DHRS13 is not a significantly differentially expressed carotenoid candidate gene (p=0.11), but it is present in head plumage, and has higher levels of expression in the head feathers than the belly. The head is where dehydrolutein is expressed, and since we think DHRS13 may be involved in the conversion of lutein to dehydrolutein, this would follow that prediction. 

There are other candidate genes with differential expression for head and belly feathers: DHRS3, SCARB2, CYP36B1, CYP27A1, CYP27C1, CYP24C1 are all upregulated in the head. This is interesting, and if further explorations indicate that the pathway of carotenoid conversion may involve multiple enzymes, or different enzymes than we currently analyze, these may be of note.

Going forward, I will look into ways to redownload the samples originally collected, and see if that doesn't resolve the corruption issues seen at the beginning of the project. For the purposes of the class project, the 8-sample study should be okay, but to make research conclusions it will be important to utilize the full data set. 

Through this process and coursework there were some difficulties I encountered, and lessons I learned. For one, it is difficult to understand the broader picture in detail when you are just learning the methods. But, from practice, I have a bit better grasp on how genomic and transcriptomic analyses work, and what checkpoints are useful to set. For instance, with the lack of an annotated genome for Gouldian finches-- I am sure there is some way to check that ahead of time, and if I had to do similar work in the future, I would know to double-check the resources I am using before getting started. 







