# find_te_ins
find_te_ins is designed to find Transposon Element (TE) insertions using long reads (nanopore), by alignment directly. (minimap2)

## Install

```bash
$ git clone https://github.com/bakerwm/find_te_ins.git
$ cd find_te_ins
```

Change the following variables upon your condition: `genome_fa` and `te_fa` in line-10 and line-11;

```
$ bash run_pipe.sh
run_pipe.sh <fa/q_dir> <out_dir>
```

## Prerequisite

+ [minimap2](https://github.com/lh3/minimap2) - 2.17-r974-dirty, align long reads to reference genome 
+ [featureCounts]() - v2.0.0, quantification
+ [samtools](https://github.com/samtools/samtools) - v1.12, working with BAM files
+ python 3.8+
+ [pysam](https://github.com/pysam-developers/pysam) 0.16.0.1, python module, working with BAM files

## Getting Started

### 1 Prepare input files 

+ **genome_fa** - reference genome in fasta format, in script `run_pipe.sh`, line-10 
+ **te_fa** - TE consensus sequence in fasta format, in script `run_pipe.sh`, line-11 
+ **long reads** - Long reads from `NanoPore` or `Pacbio`, in fasta or fastq format

### 2 Run pipe

```bash
$ cd ~/work/te_ins
# specify the path of long reads data: <path-to_long-reads>/
$ git clone https://github.com/bakerwm/find_te_ins.git 
$ bash find_te_ins/run_pipe.sh <path-to-long-reads>/ results

[1/9] align to reference genome
[2/9] extract raw insertions from BAM, by CIGAR
[3/9] convert raw insertions to fasta format
[4/9] align raw_insertion to transposon
[5/9] extract transposon name for insertions
[6/9] merge raw_insertions by window=100
[7/9] count reads for each insertion
[8/9] save final insertions to file
[9/9] Done!
```

### 3 Output

The following files listed below are the output of the pipeline, the TE insertions saved in file `*.te_ins.final.bed`

```
$ tree -L 2 results/ONT_sample-1
.
├── ONT_sample-1
│   ├── ONT_sample-1.bam
│   ├── ONT_sample-1.bam.bai
│   ├── ONT_sample-1.raw_ins.bed
│   ├── ONT_sample-1.raw_ins.fa
│   ├── ONT_sample-1.raw_ins.fa.bam
│   ├── ONT_sample-1.raw_ins.fa.bam.bai
│   ├── ONT_sample-1.te_ins.bed
│   ├── ONT_sample-1.te_ins.final.bed
│   ├── ONT_sample-1.te_ins.final.bed6
│   ├── ONT_sample-1.te_ins.gtf
│   ├── ONT_sample-1.te_ins.quant.stderr
│   ├── ONT_sample-1.te_ins.quant.stdout
│   ├── ONT_sample-1.te_ins.quant.txt
│   ├── ONT_sample-1.te_ins.quant.txt.summary
│   ├── ONT_sample-1.te_ins.raw.txt
│   ├── run_minimap2.dm6.stderr
│   └── run_minimap2.dm6_transposon.stderr
...
```

**{sample_name}.te_ins.final.bed** 

```
column 1. chr name of reference 
column 2. start pos of Insertion 
column 3. end pos of Insertion 
column 4. insertion name 
column 5. a fixed integer [255]  
column 6. strand # in current version, not consider the dirction of TE insertions !!!
column 7. name of TE consensus 
column 8. length of TE consensus  
column 9. proportion of the TE consensus identified  
column 10. number of supported reads for the insertion 
column 11. number of all reads cover the insertion 
column 12. proportion TE supported reads 
column 13. type of the TE insertions [full, p3, p5]
```

**{sample_name}.te_ins.raw.txt** 

`column 16` (last column), is the type of TE insertions: [full, p3, p5] 

+ full, more then cutoff [60%] of the TE consensus were detected  
+ p3, only the 3' end of the TE consensus were detected  
+ p5, only the 5' end of the TE consensus were detected 

In the `.final.bed` file, ONLY **full** TE insertions were saved for further analysis


### Change criteria

TE types were defined in `run_pipe.sh` by `anno_te.py`, the criteria `-c 0.6` could be changed to [0-1] float number based on your condition.
see line-100 in file `run_pipe.sh` 

```
# line-100 of run_pipe.sh
[[ ! -f ${te_ins_txt} ]] && python ${src_dir}/anno_te.py -x ${te_fa_fai} ${te_bam} | sort -k4,4 -k5,5n > ${te_ins_txt}

# change criteria to 0.7
[[ ! -f ${te_ins_txt} ]] && python ${src_dir}/anno_te.py -x ${te_fa_fai} -c 0.7 ${te_bam} | sort -k4,4 -k5,5n > ${te_ins_txt}

# remove te_ins files, and run the command again
$ rm results/ONT_sample-1.te_ins*
$ bash find_te_ins/run_pipe.sh <path-to-long-reads>/ results
```


## How it works? 

1. extract INSERTIONS 




















