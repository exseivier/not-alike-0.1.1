# Not-Alike pipeline

### Overview

This pipeline identifies genomic regions which are the most dissimilar compared at least to one genome of a database composed with several genomes sequences. It uses bioinformatic tools to acomplish this task such as: blastn, hisat2, samtools, stringtie and gffread. This pipeline has three steps: **(1)** split query genome, **(2)** iterative searching and **(3)** mapping & assembling. The result is a GTF file generated by stringtie which assembles the dissimilar fragments not-alike found. Using gffread you can get the sequences of those assembled fragments.

### Introduction

This pipeline finds dissimilar regions of the genome of interest by mean of the interative comparision between all posible fragments of size ws at each step size ss of the genome. It is written in python language and uses external bioinformatic tools such as: blastn, hisat2, samtools, stringtie and gffread. Also it uses the following python packages: click, os, random and subprocess that handles the communication with python and the external tools. As we pointed out above, it has three steps.

**Split genome step**. In this step the genome of interest is split into fragments of defined size (ws) with a sliding window method that cuts fragments at each step of determined size (ss). In the **iterative searching step**, fragments are used as query to search the genomes database, genome by genome, and fragments that hit a subsequence of the searched genome are elimineated form fragments file. This procedure is performed iterativelly until the end of genomes database. And the final step, the **mapping & assembly**. In this step the remained fragments are subject to be mapped to reference genome, which is the same genome where fragments come from, with hisat2, samtools and stringtie. The resulting GTF file contains the coordinates of the assembled fragments. Those assembled fragments are the dissimilar regions.


### Usage

To ask for help.

```bash
not-alike --help
```

To ask for help for the **search** command.

```bash
not-alike search --help
```

Task execution.

```bash
not-alike search -g test/seqs/Cgla.genome.fasta -ws 1000 -ss 100 -db test/db/noCgla2.txt -e 10 -i 50 -q 50 -t megablast -c 'Leave a comment encolsed by single quotes'
```

This pipeline requires two input files, one is the fasta file of the genome of interest (or query genome), and the another is a file with the name of the BLAST DB formated files of the genomes that are wanted to be in the database (in this example noCgla2.txt). This file must be inside the folder that holds all BLAST DB formatted files. The first file is set with the -g flag and the second with the -db flag. E-value (-e), identity percentage (-i) and query HSP coverage pewrcentage (-q) cutoff values are also set with the corresponding flags.

### Commands

## database-makeblastdb

This function...

## database-makefiledb

## search

## show-database

## show-exp

### Download & install


Download the zipped file at [exseivier's github](https://www.github.com/exseivier/not-alike)

unzip not-alike-main.zip

cd not-alike-main

sudo pip3 install dist/Not-Alike-0.0.1.tar.gz
