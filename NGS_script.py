# WGS analysis for Candida
# Workload manager: Slurm
# sbatch is used to submit a job script for later execution.
# The script will typically contain one or more srun commands to launch parallel tasks.

# Data preparation:
# 2 fastq files: forward and reverse
# a reference fasta file for alignment

# The pipeline contains following software: fastqc, trimmomatic, spades, quast, bwa, samtools，GATK and snpEff
# The outcome are stored in 4 directories:
# fastqc results -> C.para_fastqc
# quast results -> C.para_quast
# assembled fasta files -> C.para_fasta
# The rest of results -> C.para_result/{}.format(sample_id)

import subprocess
import os


def sbatch(job_name, command, time=4, mem=60, tasks=20, dep=''):
    if dep != '':
        dep = '--dependency=afterok:{} --kill-on-invalid-dep=yes '.format(dep)

    sbatch_command = "sbatch -J {} -o {}.out -e {}.err --mail-user=$USER@ntu.edu.sg --mail-type=FAIL -t {}:00:00 --mem={}000 --ntasks-per-node={} --wrap='{}' {}".format(
        job_name, job_name, job_name, time, mem, tasks, command, dep)
    sbatch_response = subprocess.getoutput(sbatch_command)
    print(sbatch_response)
    job_id = sbatch_response.split(' ')[-1].strip()
    return job_id


def fastqc():
    command = "fastqc {}_R1.fastq -o /home/zzhang082/C.para_fastqc".format(sample_id)
    job_id = sbatch('fastqc', command)
    return job_id


def trimmomatic():
    command = "java -jar trimmomatic-0.39.jar PE {}_R1.fastq.gz {}_R2.fastq.gz \
    {}_R1_paired.fastq.gz {}_R1_unpaired.fq.gz {}_R2_paired.fastq.gz {}_R2_unpaired.fq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36".format(sample_id, sample_id, sample_id, sample_id, sample_id, sample_id)
    job_id = sbatch('trimmomatic', command)
    return job_id


def spades(dep=''):
    command = "spades.py -t 8 -1 {}_R1_paired.fastq.gz -2 {}_R2_paired.fastq.gz -o /home/zzhang082/C.para_fasta/{}".format(sample_id, sample_id, sample_id)
    job_id = sbatch('spades', command, dep=dep)
    return job_id


def quast(dep=''):
    command = "quast.py --fungus /home/zzhang082/C.para_fasta/{}/contigs.fasta -o /home/zzhang082/C.para_quast/{}".format(sample_id, sample_id)
    job_id = sbatch('quast', command, dep=dep)
    return job_id


def alignment(dep=''):
    rg = '@RG\tID:group_n\tLB:library_n\tPL:illumina\tPU:unit1\tSM:sample_n'
    command = "bwa mem -R {} -t 8 {} {}_R1_paired.fastq.gz {}_R2_paired.fastq.gz > {}.sam".format(rg, ref, sample_id, sample_id, sample_id)
    job_id = sbatch('align', command, time=8, mem=120, dep=dep)
    return job_id


def convert(dep=''):
    command = "samtools view -S -B {}.sam > {}.bam".format(sample_id,sample_id)
    job_id = sbatch('convert', command, dep=dep)
    return job_id


def sort(dep=''):
    command = "samtools sort {}.bam -o {}.sorted.bam".format(sample_id, sample_id)
    job_id = sbatch('sort', command, dep=dep)
    return job_id


def mark_duplicates(dep=''):
    command = "gatk MarkDuplicates -I {}.sorted.bam -O {}.sorted_marked.bam -M metrics.txt".format(sample_id, sample_id)
    job_id = sbatch('mark_duplicates', command, dep=dep)
    return job_id


def index_bam(dep=''):
    command = "gatk BuildBamIndex -I {}.sorted_marked.bam".format(sample_id)
    job_id - sbatch('index_bam', command, dep=dep)
    return job_id


def variant(dep=''):
    command = "gatk HaplotypeCaller -ploidy 2 -R {} -I {}.sorted_marked.bam -o {}.vcf".format(ref, sample_id, sample_id)
    job_id = sbatch('variant', command, dep=dep)
    return job_id


def Annotation(dep=''):
    command = "java -Xmx8g -jar snpEff.jar Candida_parapsilosis_cdc317 {}.vcf > {}.ann.vcf".format(sample_id, sample_id)
    job_id = sbatch('Annotation', command, dep=dep)
    return job_id


ref = '/home/zzhang082/C.para_ref/*.fasta'

os.system("module load bwa")
os.system("module load samtools")
os.system("module load anaconda3")
os.system("export PATH=$PATH:$PWD/SPAdes-3.15.5-Linux/bin/")

# index fasta file
os.system("bwa index {}".format(ref))
os.system("samtools faidx {}".format(ref))


files = os.listdir("/home/zzhang082/C.para")
for file in files:
    # extract sample_id
    prefix = os.path.splitext(file)[0]
    sample_id = prefix[:-9]

    # create directory for each sample
    os.system("mkdir /home/zzhang082/C.para_result/{}".format(sample_id))

    # copy fastq file to work directory
    os.chdir("/home/zzhang082/C.para")
    os.system("cp {}_R1.fastq.gz /home/zzhang082/C.para_result/{}".format(sample_id, sample_id))
    os.system("cp {}_R2.fastq.gz /home/zzhang082/C.para_result/{}".format(sample_id, sample_id))

    # change work directory
    os.chdir("/home/zzhang082/C.para_result/{}".format(sample_id))

    # run the pipeline
    fastqc_jobid = fastqc()
    trimmomatic_jobid = trimmomatic()

    # create directory to store assembled fasta file from spades
    os.system("mkdir /home/zzhang082/C.para_fasta/{}".format(sample_id))

    spades_jobid = spades(trimmomatic_jobid)

    # create directory to store quast results
    os.system("mkdir /home/zzhang082/C.para_quast/{}".format(sample_id))
    quast_jobid = quast(spades_jobid)

    alignment_jobid = alignment(trimmomatic_jobid)
    convert_jobid = convert(alignment_jobid)
    sort_jobid = sort(convert_jobid)
    mark_duplicates_jobid = mark_duplicates(sort_jobid)
    index_bam_jobid = index_bam(mark_duplicates_jobid)
    variant_jobid = variant(index_bam_jobid)
    SNP_jobid = SNP(variant_jobid)
    Annotation_jobid = Annotation(SNP_jobid)




