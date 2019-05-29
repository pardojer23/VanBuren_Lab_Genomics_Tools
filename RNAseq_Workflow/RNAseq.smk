"""SnakeMake workflow to process RNAseq reads.

Takes input of transcriptome fasta file and RNAseq fastq files. Outputs TPM and differential expression
"""

## Created by: Jeremy Pardo 04/01/2019
from os.path import join
import glob
import re
import subprocess
## set global variables

with open('RNAseq_snakeMake_config.txt','r') as config_file:
    for line in config_file:
        if line.startswith("##"):
            pass
        elif line.startswith("script_dir="):
            sd = line.strip().split('=')
            script_dir = sd[1]
        elif line.startswith("transcriptome="):
            t = line.strip().split('=')
            transcriptome = t[1]
        elif line.startswith('fastq_dir='):
            f = line.strip().split('=')
            fastq_dir = f[1]
        elif line.startswith('output_dir='):
            o = line.strip().split("=")
            output_dir = o[1]
        elif line.startswith('prefix='):
            p = line.strip().split('=')
            prefix = p[1]
        elif line.startswith('sample_pattern'):
            sp = line.strip().split('=')
            sample_pattern = sp[1]
        elif line.startswith("read1="):
            r1 = line.strip().split("=")
            read1 = r1[1]
        elif line.startswith("read2="):
            r2 = line.strip().split("=")
            read2 = r2[1]
        elif line.startswith("transcript_suffix="):
            ts = line.strip().split('=')
            transcript_suffix = ts[1]
        elif line.startswith("compare_list="):
            cl = line.strip().split("=")
            compare_list = cl[1]
samples, = glob_wildcards(join(fastq_dir,'{sample,'+prefix+sample_pattern+'}'+read1))
Pattern_R1 = '{sample}'+read1
Pattern_R2 = '{sample}'+read2
print(compare_list)

def make_tx2gene(transcript_file):
    ''' Make tx2gene file describing mapping between transcripts and genes needed for tximport.

    Output is a 2 column text file where the left column is the transcript_ID and right column is gene_ID.

    Parameters:
        transcript_file: path to transcriptome in fasta format.'''
    with open(output_dir+"tx2gene.txt",'w+') as tx2gene:
        tx2gene.write("transcript_ID\tgene_ID\n")
        print(transcript_file)
        with open(transcript_file,'r') as infile:
            for line in infile:
                if line.startswith(">"):
                    line.strip()
                    transcript = line.split()[0]
                    transcript = transcript.replace(">","")
                    gene = re.sub(transcript_suffix,"",transcript)
                    tx2gene.write(transcript+'\t'+gene+'\n')

def make_sample_table(outputDir):
    '''make sample table containing metadata for all samples and paths to the quant.sf files output from salmon
    make_sample_table function makes an empty sampleTable with column names Sample, Condition, TimePoint, Replicate, FilePath.

    Parameters:
        outputDir: Path to output directory where sample table will be created
    '''
    with open(outputDir+"SampleTable.txt",'w+') as SampleTable:
        SampleTable.write("Sample\tCondition\tTimePoint\tReplicate\tFilePath\n")
    sample_files = glob.glob(output_dir+'salmon_aln_*')
    for i in range(len(sample_files)):
        sample = sample_files[i].replace(outputDir,"")
        print(sample)
        add_to_sample_table(sample)

def add_to_sample_table(sample):
    '''add_to_sample_table function adds sample information based on file name.

    Function expects samples to be named with format [prefix]_[Condition]_[Timepoint]_[Replicate]_[suffix]

    Parameters:
        sample: name of sample following [prefix]_[Condition]_[Timepoint]_[Replicate] pattern.
    '''
    RNA_sample = sample.strip().split(sep="_")
    with open(output_dir+"SampleTable.txt",'a') as SampleTable:
        SampleTable.write(RNA_sample[3]+RNA_sample[4]+RNA_sample[5]+'\t'+
                          RNA_sample[3]+'\t'+RNA_sample[4]+'\t'+RNA_sample[5]+'\t'+
                          output_dir+sample+'/quant.sf'+'\n')



###
## snakemake rules to execute workflow
###
## rule all has the final output as the dependency insures that the whole pipeline is run
rule all:
    input:
        output_dir+"DESeq2.RData"

## create salmon transcriptome index
rule salmon_index:
    input:
         transcriptome

    output:
          directory(output_dir+prefix+"_transcripts_index")
    threads:
            24
    shell:
         "salmon index -t {input} -i {output} -k 31"
## map RNAseq files to transcriptome using salmon in pseudoalignment mode
rule salmon_map:
    input:
         transcript_index= output_dir+prefix+"_transcripts_index",
         fwd = fastq_dir+Pattern_R1,
         rev = fastq_dir+Pattern_R2
    threads:
            48
    output:
          directory(output_dir+'salmon_aln_{sample}')

    shell:
         "salmon quant --threads 48 \
         -i {input.transcript_index} \
         -l A --seqBias --gcBias --validateMappings\
         -1 {input.fwd} -2 {input.rev} \
         -o {output}"
## create tx2gene file
rule tx2gene:
    input:
         output_dir+prefix+"_transcripts_index"
    output:
          output_dir+'tx2gene.txt'
    threads:
           1
    run:
        print(transcriptome)
        make_tx2gene(transcriptome)
## create sampleTable
rule sampleTable:
    input:
         expand(output_dir+"salmon_aln_{sample}",sample=samples)
    output:
          output_dir+"SampleTable.txt"
    threads:
            1
    run:
        make_sample_table(output_dir)
## run tximport to combine salmon outputs and calculate TPM at the gene level
rule run_tximport:
    input:
         sampleTable = output_dir+"SampleTable.txt",
         tx2gene = output_dir+'tx2gene.txt',
         outputDir = output_dir,
         scriptDir = script_dir

    output:
          output_dir+'txi.RData'
    threads:
            1
    shell: "Rscript {input.scriptDir}/tximport.r {input.outputDir}"
## run DESeq2 to calculate differential expression using the formula y~Condition
rule run_DESeq2:
    input:
         txi = output_dir+'txi.RData',
         sampleTable = output_dir+"SampleTable.txt",
         outputDir = output_dir,
         scriptDir = script_dir


    output:
          output_dir+"DESeq2.RData"
    threads:
           24
    shell: "Rscript {input.scriptDir}/DESeq2.r {input.outputDir} "+compare_list





