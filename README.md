# code-examples

Code examples from my research at Syracuse University.

run_cufflinks.pl contains a pipeline to run Tophat, Cufflinks, Cuffmerge, Cuffquant, Cuffdiff, and Cuffnorm.

fastq_processing.pl performs adapter trimming, low quality read trimming, and bowtie allignment for micro-RNA data.

clustering_algorithm_1.pl performs clustering of genes across the genome based on similarities in expression, where clusters are defined by a minimum number of genes and minimum percentage of genes in the cluster that share a specific expression characteristic.

chromatin_reads_2_3.pl maps histone ChIP-seq data and bins read count information by chromosome position for downstream analysis.
