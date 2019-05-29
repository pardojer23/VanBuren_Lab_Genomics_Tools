import datetime
import subprocess
import argparse
import os
import shutil

class SnakeMakeWrapper:
    '''Class defines a wrapper to pass command line arguments to snakeMake via a config file.

    Arguments:
        -t, --transcriptome',Path to transcript fasta file
        -f, --fastq_dir', Path to directory containing trimmed RNAseq fastq files.
                        Samples should be named with the pattern [Prefix]_[Condition]_[Timepoint]_[Replicate]
        -o, --output_dir', Path to output directory
        -s, --sample_pattern, Regular expression defining the sample name,default: _.*_S.*
        -p, --prefix, Species prefix, must match sample filenames. Example: -p Enin for samples named with pattern Enin_.*_S.*
        -r1, --read1', Fastq file suffix for read1 files,default: _R1_001_trimmed.fastq.gz
        -r2, --read2',Fastq file suffix for read2 files,default:_R2_001_trimmed.fastq.gz
        -tx, --transcript_suffix, Regular expression matching transcript suffix appended to geneIDs to make them transcript_IDs, default:-RA
        -sd, --script_dir,Path to directory containing executables for pipeline
        -cl, --condition_list, comma seperated list of conditions for differential expression. The first listed condition will be treated as the reference'''
    def __init__(self):
        '''__init__ method gets command line arguments'''
        parser = argparse.ArgumentParser()
        parser.add_argument('-t', '--transcriptome', help="Path to transcript fasta file")
        parser.add_argument('-f', '--fastq_dir', help="Path to directory containing trimmed RNAseq fastq files. "
                                                      "Samples should be named with the pattern [Prefix]_[Condition]_[Timepoint]_[Replicate]")
        parser.add_argument('-o', '--output_dir', help="Path to output directory")
        parser.add_argument('-s', '--sample_pattern', help="Regular expression defining the sample name",
                            default=".*_")
        parser.add_argument('-p', '--prefix',
                            help="Species prefix, must match sample filenames. Example: -p Enin for samples named with pattern Enin_.*_S.*")
        parser.add_argument('-r1', '--read1', help="fastq file suffix for read1 files",
                            default="_R1_001_trimmed.fastq.gz")
        parser.add_argument('-r2', '--read2', help="fastq file suffix for read2 files",
                            default="_R2_001_trimmed.fastq.gz")
        parser.add_argument('-tx', '--transcript_suffix',help="regular expression matching transcript suffix appended to geneIDs to make them transcript_IDs",default="-RA")
        parser.add_argument('-sd','--script_dir',help="Path to directory containing executables for pipeline")
        parser.add_argument('-cl','--compare_list',help='Comma seperated list of comparision for DE. Each comparision should include two dash seperated values the 1st is the reference level')
        args = parser.parse_args()
        self.transcriptome = args.transcriptome
        self.fastq_dir = args.fastq_dir
        self.output_dir = args.output_dir
        self.sample_pattern = args.sample_pattern
        self.prefix = args.prefix
        self.read1 = args.read1
        self.read2 = args.read2
        self.transcript_suffix = args.transcript_suffix
        self.script_dir = args.script_dir
        self.compare_list = args.compare_list
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)

    def get_config_file(self):
        '''Create config file containing the command line arguments to pass to snakeMake'''
        with open("RNAseq_snakeMake_config.txt",'w+') as configFile:
            configFile.write('## Path to script directory \n'+
                             'script_dir='+self.script_dir+'\n'+
                             '## Path to transcriptome fasta file \n'+
                             'transcriptome='+self.transcriptome+'\n'+
                             '## Path to fastq directory containing trimmed RNAseq reads \n'+
                             '## Note: samples should be named with the pattern [Prefix]_[Condition]_[Timepoint]_[Replicate] \n'+
                             'fastq_dir='+self.fastq_dir+'\n'+
                             '## Path to output directory \n'+
                             'output_dir='+self.output_dir+'\n'+
                             '## Regular expression defining sample name without prefix. Default =  _.*_S.*'+'\n'+
                             'sample_pattern='+self.sample_pattern+'\n'+
                             '## species prefix for sample name. Example: Enin_[Condition]_[Timepoint]_[Replicate]_S[]_R1_001_trimmed.fastq.gz would have the prefix Enin \n'+
                             'prefix='+self.prefix+'\n'+
                             '## suffix for read1 (forward) samples. Example: [Prefix]_[Condition]_[Timepoint]_[Replicate]_S[]_R1_001_trimmed.fastq.gz would have the suffix _R1_001_trimmed.fastq.gz \n'+
                             'read1='+self.read1+'\n'+
                             '## suffix for read1 (forward) samples. Example: [Prefix]_[Condition]_[Timepoint]_[Replicate]_S[]_R2_001_trimmed.fastq.gz would have the suffix _R2_001_trimmed.fastq.gz \n'+
                             'read2='+self.read2+'\n'+
                             '## regular expression defining suffix for transcript ID that differentiates it from gene ID \n'+
                             'transcript_suffix='+self.transcript_suffix+'\n'+
                             'compare_list='+self.compare_list)


    def get_launch_script(self):
        '''Create shell script to launch snakeMake'''
        with open('run_snakeMake.sh','w+') as launchScript:
            launchScript.write('#!/usr/bin/sh\n'+
                                'snakemake -s '+self.script_dir+'/RNAseq.smk -j 500 '
                                '--cluster-config '+self.script_dir+'/cluster.json '
                                '--latency-wait 30 '
                                '--cluster \'sbatch --mem={cluster.mem} -t {cluster.time} -c {threads}\'')



def main():
    '''Main method creates config file and launches snakeMake.'''
    start = datetime.datetime.now()
    print("Starting the pipeline at: {0}".format(start))
    my_SnakeMakeWrapper = SnakeMakeWrapper()
    my_SnakeMakeWrapper.get_config_file()
    my_SnakeMakeWrapper.get_launch_script()
    os.chmod('./run_snakeMake.sh',0o755)
    shutil.copy(my_SnakeMakeWrapper.script_dir+'/cluster.json','./')
    subprocess.run("./run_snakeMake.sh")
    print("Pipeline finished in: {0}".format(datetime.datetime.now()-start))
if __name__ == '__main__':
    main()
