#############################################
# Amplicon sequencing processing pipeline   #
# Used for Single End 400 reads             #
# Burst for closed-reference clustering     #
# Winner take all for combined two targets  #
#############################################

# This is a Snakemake pipeline for processing amplicon sequencing data based on 
# BGISEQ single-end method. The insert library is build without controling the sequencing
# direction. So Read1 can start from either the forward primer or reverse primer, 
# randomly and in equal number.

# If you are not familar with amplicon sequencing, or you have only done amplicon
# sequencing on illumina's platform, I would suggest you to read the following
# instructions carefully (at least once). Most of the illumina based protocol controls
# the direction of amplicons, so Read1 always occurs from the forward primers used.
# This is NOT the case for our BGISEQ based protocol, and will result in some differences
# in data processing (and some intersting aftermath). For single end mode, we will have
# only one read for each sample.

# Two types of amplicons appears in the SE data:
	# TYPE 1: Read1-forward primer-insert-reverse primer
	# TYPE 2: Read1-reverse primer-insert-forward primer
# For SE data, these two types of amplicon may cover quite different regions.
# For the commands in this pipeline, refer the forward and reverse primers as forward and reverse.
# Most of these labels appear as suffixes.

# You will need the sequence of the forward and reverse primers (or the conservative region
# you like to trim off). These primers are different based on the amplified target.
# Change the primer sequences in config.yaml to the one you used (and ask around if you don't know).
# Also check if you have the right Burst database (bacteria or fungi?).
# Check cluster.yaml if you use qsub.

# Get the sample label from data/samples
# The sample file should be in the format [SampleName].[ReadDirection].fq (For SE data, only r1)
# Put all sample files in data/samples
# A sample should have one corresponding files

configfile: 'config.yaml'

import fnmatch
import os
SAMPLES = {}
for file in os.listdir('data/samples/'):
	if fnmatch.fnmatch(file, '*.fq'):
		string = file.split('.')
		sampleName = string[0]
		SAMPLES[sampleName] = 'data/samples/' + file

# Target of this pipeline: combined.biom is the final output
rule target:	
	input:
		'data/combined.biom',
		'data/combined.taxa.biom',
		'data/combined.taxa.tsv'

# Find the two directions in the SE data		
rule direction:
	input:
		lambda wildcards: SAMPLES[wildcards.sample]		
	output:
		fw_cut1=temp('data/split_reads/{sample}.forward.cut1.fq'),
		fw_cut2=temp('data/split_reads/{sample}.forward.cut2.fq'),
		rv_cut1=temp('data/split_reads/{sample}.reverse.cut1.fq'),
		rv_cut2=temp('data/split_reads/{sample}.reverse.cut2.fq')
	params:
		fw=config['primers']['fw'],
		rv=config['primers']['rv']
	log:
		'logs/split_reads/{sample}.log'
	threads: 2
	run:
		base = {'A':'T','T':'A','C':'G','G':'C','R':'Y','Y':'R','K':'M','M':'K','S':'S','W':'W','B':'V','V':'B','D':'H','H':'D','N':'N'}
		fwrc = ''.join([base[i] for i in list({params.fw})[0][::-1]])
		rvrc = ''.join([base[i] for i in list({params.rv})[0][::-1]])
		shell('cutadapt {input} -g {params.fw} -o {output.fw_cut1} --discard-untrimmed -e 0.1 --overlap 60 -j {threads} > {log}')
		shell('cutadapt {output.fw_cut1} -a rvrc -o {output.fw_cut2} -e 0.1 -j {threads} >> {log}'.replace('rvrc', rvrc))
		shell('cutadapt {input} -g {params.rv} -o {output.rv_cut1} --discard-untrimmed -e 0.1 --overlap 60 -j {threads} >> {log}')
		shell('cutadapt {output.rv_cut1} -a fwrc -o {output.rv_cut2} -e 0.1 -j {threads} >> {log}'.replace('fwrc', fwrc))

# Remove low quality reads		
rule quality:
	input:
		fw_cut2='data/split_reads/{sample}.forward.cut2.fq',
		rv_cut2='data/split_reads/{sample}.reverse.cut2.fq'
	output:
		fw_qc='data/qc_reads/{sample}.forward.qc.fq',
		rv_temp=temp('data/qc_reads/{sample}.reverse.temp.fq'),
		rv_qc='data/qc_reads/{sample}.reverse.qc.fq'
	params:
		truncLen=config['qc']['truncLen']
	log:
		'logs/qc_reads/{sample}.log'
	threads: 2
	run:
		shell('vsearch --fastq_filter {input.fw_cut2} --fastq_trunclen {params.truncLen} --fastq_truncee 1 --fastq_minlen 100 --fastq_maxns 0 --fastqout {output.fw_qc} --fasta_width 0 --threads {threads} 2> {log}')
		shell('vsearch --fastq_filter {input.rv_cut2} --fastq_trunclen {params.truncLen} --fastq_truncee 1 --fastq_minlen 100 --fastq_maxns 0 --fastqout {output.rv_temp} --fasta_width 0 --threads {threads} 2>> {log}')
		shell('seqtk seq -r {output.rv_temp} > {output.rv_qc}')

# Trim off sequences not in database and convert to FASTA
rule trim_tail:
	input:
		fw_qc='data/qc_reads/{sample}.forward.qc.fq',
		rv_qc='data/qc_reads/{sample}.reverse.qc.fq'
	output:
		fw_qc=temp('data/qc_reads/{sample}.forward.qc.fa'),
		rv_qc=temp('data/qc_reads/{sample}.reverse.qc.fa')
	params:
		fw_trim = config['trim']['fw'],
		rv_trim = config['trim']['rv']
	run:
		shell('seqtk trimfq -b {params.fw_trim} {input.fw_qc} | seqtk seq -A -L 50 - > {output.fw_qc}')
		shell('seqtk trimfq -e {params.rv_trim} {input.rv_qc} | seqtk seq -A -L 50 - > {output.rv_qc}')
		
# Align reads to the database
rule align:
	input:
		fw_qc='data/qc_reads/{sample}.forward.qc.fa',
		rv_qc='data/qc_reads/{sample}.reverse.qc.fa'
	output:
		fw_aln='data/alignments/{sample}.forward.b6',
		rv_aln='data/alignments/{sample}.reverse.b6'
	params:
		acx=config['database']['acx'],
		edx=config['database']['edx'],
		mode=config['align']['mode'],
		id=config['align']['id']
	log:
		'logs/alignments/{sample}.log'
	threads: 4
	run:
		shell('burst12 -q {input.fw_qc} -a {params.acx} -r {params.edx} -o {output.fw_aln} -i {params.id} -m {params.mode} -t {threads} > {log}')
		shell('burst12 -q {input.rv_qc} -a {params.acx} -r {params.edx} -o {output.rv_aln} -i {params.id} -m {params.mode} -t {threads} >> {log}')

# Combine alignments using winner take all method
rule winnerTakeAll:
	input:
		fw_aln='data/alignments/{sample}.forward.b6',
		rv_aln='data/alignments/{sample}.reverse.b6'
	output:
		tsv='data/profiles/{sample}.tsv',
		biom='data/profiles/{sample}.biom',
		tsv_fw='data/profiles/{sample}.fw.tsv',
		tsv_rv='data/profiles/{sample}.rv.tsv',
		biom_fw='data/profiles/{sample}.fw.biom',
		biom_rv='data/profiles/{sample}.rv.biom',
	params:
		sampleName='{sample}'
	log:
		'logs/winnerTakeAll/{sample}.log'
	run:
		shell('amplicon_winnerTakeAll.py -i {input.fw_aln} {input.rv_aln} -sn {params.sampleName} -t {output.tsv} -b {output.biom} -l > {log}')
		shell('amplicon_winnerTakeAll.py -i {input.fw_aln} -sn {params.sampleName} -t {output.tsv_fw} -b {output.biom_fw} -g >> {log}')
		shell('amplicon_winnerTakeAll.py -i {input.rv_aln} -sn {params.sampleName} -t {output.tsv_rv} -b {output.biom_rv} -g >> {log}')

# Concatenate all profiles into one
rule concat:
	input:
		both=expand('data/profiles/{sample}.biom', sample=SAMPLES),
		fw=expand('data/profiles/{sample}.fw.biom', sample=SAMPLES),
		rv=expand('data/profiles/{sample}.rv.biom', sample=SAMPLES)
	output:
		combined = 'data/combined.biom',
		biom_taxa = 'data/combined.taxa.biom',
		tsv_taxa = 'data/combined.taxa.tsv',
		combined_fw = 'data/combined.fw.biom',
		biom_taxa_fw = 'data/combined.taxa.fw.biom',
		tsv_taxa_fw = 'data/combined.taxa.fw.tsv',
		combined_rv = 'data/combined.rv.biom',
		biom_taxa_rv = 'data/combined.taxa.rv.biom',
		tsv_taxa_rv = 'data/combined.taxa.rv.tsv'
	params:
		taxa=config['database']['tax']
	log:
		'logs/concat/concate.log'
	run:
		shell('amplicon_concat.py -i {input.both} -biom_out {output.combined} > {log}')
		shell('biom add-metadata -i {output.combined} -o {output.biom_taxa} --observation-metadata-fp {params.taxa} --observation-header OTUID,taxonomy --output-as-json --sc-separated taxonomy')
		shell('biom convert -i {output.biom_taxa} -o {output.tsv_taxa} --to-tsv --header-key taxonomy')
		shell('amplicon_concat.py -i {input.fw} -biom_out {output.combined_fw} > {log}')
		shell('biom add-metadata -i {output.combined_fw} -o {output.biom_taxa_fw} --observation-metadata-fp {params.taxa} --observation-header OTUID,taxonomy --output-as-json --sc-separated taxonomy')
		shell('biom convert -i {output.biom_taxa_fw} -o {output.tsv_taxa_fw} --to-tsv --header-key taxonomy')
		shell('amplicon_concat.py -i {input.rv} -biom_out {output.combined_rv} > {log}')
		shell('biom add-metadata -i {output.combined_rv} -o {output.biom_taxa_rv} --observation-metadata-fp {params.taxa} --observation-header OTUID,taxonomy --output-as-json --sc-separated taxonomy')
		shell('biom convert -i {output.biom_taxa_rv} -o {output.tsv_taxa_rv} --to-tsv --header-key taxonomy')