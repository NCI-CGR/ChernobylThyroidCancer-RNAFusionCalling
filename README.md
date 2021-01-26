# Consensus RNA-seq Fusion Calling

## I. Description
To improve precision and recall when detecting gene fusions, we used a rank-based scoring method based on the approach by [Liu et al 2015](https://doi.org/10.1093/nar/gkv1234). The four callers used in this study were selected based on their performance in Liu et al 2015 and [Kumar et al 2016](https://doi.org/10.1038/srep21597):
* Chimerascan
* EricScript
* Mapsplice2
* STAR-Fusion

Major steps in the workflow include:
1. Each fusion calling tool has a Snakemake subworkflow that runs the tool and calculates the scores for each fusion based on the number of supporting reads. For each sample, the fusion with the most supporting reads will have a score of 1.0, and the remaining fusions will be scored relative to the total number of fusions down to zero. 
    * For instance, if there are 4 fusions detected by STAR-Fusion in Sample01.fastq.gz, the fusions will be scored 1.0, 0.75, 0.5, and 0.25 from most reads to fewest.
    
2. The main Snakefile takes the scores for each tool and sums them for the total score of each fusion in a sample.
    * For instance, a fusion in Sample01.fastq.gz that scored 1.0, 1.0, 0.75, and 0.5 across the four tools will have a total score of 3.25.
    * If a tool did not detect a fusion that was found by another tool, it receives no score. A fusion only detected by three tools might have a score of 1.0, 1.0, 0.5, and 0 for a total score of 2.5. 

3. The final output table is a compilation of all fusions and scores for all samples. 
    


## II. Dependencies
  - [Snakemake](https://snakemake.readthedocs.io/en/stable/)
  - [pandas](https://pandas.pydata.org/docs/)
  - [Samtools](http://www.htslib.org/)
  - [Chimerascan](https://code.google.com/archive/p/chimerascan/)
  - [Ericscript](https://sites.google.com/site/bioericscript/)
  - [MapSplice2](http://www.netlab.uky.edu/p/bioinfo/MapSplice2)
  - [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion)
  
## III. Input
The workflow inputs zipped fastq files as designated in the fusion_config.yaml file. Some tools require unzipped fastqs, which is handled by the Snakefile. 
 
## IV. Output
Each tool has its own output file directory that generally include bam files of the supporting reads and a file with the fusion calls. A table of fusions and scores is created for each tool individually, and there is one summary table with the final summed scores for all samples and fusions. 

