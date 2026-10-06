### FastQC

sudo apt update
sudo apt install fastqc default-jre

fastqc --version
FastQC v0.12.1

java -version
openjdk version "21.0.12.1" 2026-08-18

### fastp

sudo apt install fastp

fastp --version
fastp 0.3.7
/home/flavio/miniforge3/envs/genome_env/bin/fastqc


conda create -n genome_env \
    --override-channels \
    -c conda-forge -c bioconda \
    --strict-channel-priority \
    fastqc fastp samtools bcftools bwa-mem2 gatk4 plink


conda activate genome_env

(genome_env) flavio@flavio:~$ conda list -n genome_env '^(fastqc|fastp|samtools|bcftools|bwa-mem2|gatk4.*|plink2?)$'
# packages in environment at /home/flavio/miniforge3/envs/genome_env:
#
# Name                     Version          Build            Channel
bcftools                   1.24             h118bc1c_2       bioconda
bwa-mem2                   2.3              he70b90d_0       bioconda
fastp                      1.3.7            h43da1c4_0       bioconda
fastqc                     0.12.1           hdfd78af_0       bioconda
gatk4                      4.6.2.0          py310hdfd78af_1  bioconda
plink                      1.90b7.7         h18e278d_1       bioconda
samtools                   1.24             h9dcdb79_1       bioconda



### BWA-MEM2, GATK, PLINK

conda create -n genome_env \
    --override-channels \
    -c conda-forge -c bioconda \
    --strict-channel-priority \
    fastqc fastp samtools bcftools bwa-mem2 gatk4 plink


conda activate genome_env

command -v fastqc fastp samtools bcftools bwa-mem2 gatk plink
