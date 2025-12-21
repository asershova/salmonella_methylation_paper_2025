# Scripts for processing of nanopore data 

### Basecalling and demultiplexing of Nanopore data.

m6A, m5C, and m4C models are used

dorado_pod5_ref_bac_meth_14-12-23.sh


### Basecalled and splitted by barcodes reads mapped to the _S.enterica_ ST4/74 reference genome
dorado_aligner_march_2025_sal_all_mod.sh


### modBam files are sorted and indexed for modkit
samtools_sort_coverage_modkit_march_2025.sh


### modBam files are transformed into bedMethyl files
modkit_pileup_2025.sh


### Find motif sequences
modkit_search_motif_4.sh


### Methylation statistics
modkit_summary.sh


### Extract modification information for individual reads
modkit_full_table.sh

