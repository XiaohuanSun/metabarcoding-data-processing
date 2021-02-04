# metabarcoding-data-processing
A full procedure of metabarcoding data including needed python packages are provided here for analyzing sequences generated from DNBSEQ-G400 with SE400 and PE200 modes, and Illumina MiSeq with PE300 mode.
For instance:
snakefile_amplicon_SE400 is the snakemake script with the corresponding config_SE400.yaml for SE400 sequencing analysis; snakemake pipline_PE200.smk and config_PE200.yaml for PE200 analysis; snakemake_amplicon_illumina_pe and config_illumina_pe.yaml for MiSeq PE300 analysis

Reads from the same sample can be analysed by simply following the snakemake file based on sequencing modes. OTU table can be generated after that, the rarefraction step could be taken after this analysis based on the rarefraction curve.

The following analyzing steps are written in short amplicon sequencing analysis.Rmd.
