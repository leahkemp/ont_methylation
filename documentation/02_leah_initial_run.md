# First methylation adaptive sampling run

*Preliminary run*

Author: Leah Kemp
Created: 2022/07/28 17:19:25
Last modified: 2022/09/05 13:25:44

- **Aim:** This document documents/describes a preliminary analysis of sample AB1052B
- **Prerequisite software:** [conda](https://docs.conda.io/en/latest/) v4.10.3, [mamba](https://mamba.readthedocs.io/en/latest/index.html) v0.17.0, [samtools](https://www.htslib.org/) v1.13, [Docker](https://www.docker.com/) v20.10.11, [bgzip](https://www.htslib.org/doc/bgzip.html) v1.13+ds, [GNU tar](https://www.gnu.org/software/tar/) v1.34
- **OS:** Leviathan (Debian GNU/Linux bookworm/sid) (ESR research network)

## Table of contents

- [First methylation adaptive sampling run](#first-methylation-adaptive-sampling-run)
  - [Table of contents](#table-of-contents)
  - [Preface](#preface)
  - [Basecalling and initial alignment](#basecalling-and-initial-alignment)
  - [Variant calling](#variant-calling)
  - [Phasing with whatshap](#phasing-with-whatshap)
  - [modbam2bed](#modbam2bed)

## Preface

This is a preliminary run to run another sample at `/data/ont_methylation/data/AB1052B/run6/fast5` through the documentation/analysis that Miles did [here](./01_miles_initial_run.md). The only modification is that we *won't* subset out on-target reads.

## Basecalling and initial alignment

Run bonito with SUP (super accuracy calling) and Remora for methylation calling. I didn't have permissions to access the bonito Miles pre-installed on Leviathan, so I installed from scratch in a conda environment (I also couldn't create a python virtual environments because of missing dependencies on Leviathan).

```bash
# get bonito
cd /data/ont_methylation/software/
git clone https://github.com/nanoporetech/bonito.git
cd bonito/

# install bontio in conda env
mamba create -n bonito python=3.7 pip
conda activate bonito
pip install --upgrade pip
pip install -r requirements.txt
python setup.py develop

# create screen to run in
screen -S bonito

# reactivate conda env I created with bontio installed
conda activate bonito

# process the fast5 data
bonito basecaller dna_r9.4.1_e8_sup@v3.3 \
/data/ont_methylation/data/AB1052B/run6/fast5 \
--modified-bases 5mC \
--reference /data/publicData/genomes/human/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.mmi \
--recursive \
--alignment-threads 16 | \
samtools view -u | samtools sort -@ 8 > /data/ont_methylation/data/AB1052B/run6/bam/adipose_1052_sup_modbases_GRCh38_20220729.bam

# index
cd /data/ont_methylation/data/AB1052B/run6/
samtools index ./bam/adipose_1052_sup_modbases_GRCh38_20220729.bam
```

## Variant calling

Get Clair3 models (http://www.bio8.cs.hku.hk/clair3/clair3_models/clair3_models.tar.gz)

```bash
cd /data/ont_methylation/data/AB1052B/run6/
wget http://www.bio8.cs.hku.hk/clair3/clair3_models/clair3_models.tar.gz
tar -xvf clair3_models.tar.gz --directory clair3_models
rm clair3_models.tar.gz
```

Get nextflow

```bash
# install nextflow in conda env
mamba create -n nextflow nextflow=22.04.0
conda activate nextflow
```

Modify the resources to blast Levithan

In line 78 of the file at `/home/lkemp/.nextflow/assets/epi2me-labs/wf-human-snp/nextflow.config`, I changed this...

```txt
executor {
    $local {
        cpus = 4
        memory = "8 GB"
    }
}
```

...to this...

```txt
executor {
    $local {
        cpus = 32
        memory = "120 GB"
    }
}
```

Run pipeline

```bash
cd /data/ont_methylation/data/AB1052B/run6/
nextflow run epi2me-labs/wf-human-snp \
-resume \
-w clair3_adipose_1052 \
-profile standard \
--model ./clair3_models/ont_guppy5 \
--bam /data/ont_methylation/data/AB1052B/run6/bam/adipose_1052_sup_modbases_GRCh38_20220729.bam \
--bed /data/ont_methylation/data/AB1052A/run1/bed/illumina450K_hg38_cpgsites_generegions_2Kpad_controls_noY.collapsed.sorted.bed \
--ref /data/publicData/genomes/human/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
--out_dir vcf_adipose_1052 \
--phase_vcf
```

## Phasing with whatshap

Get software

```bash
mamba create -n whatshap -c bioconda -c conda-forge whatshap=1.4 nomkl=3.0
conda activate whatshap
```

Run

```bash
cd /data/ont_methylation/data/AB1052B/run6/
# run the phasing
whatshap haplotag \
--ignore-read-groups \
--output ./bam/adipose_1052_sup_modbases_GRCh38_20220729.hp.bam \
--reference /data/publicData/genomes/human/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
./vcf_adipose_1052/all_contigs.vcf.gz \
./bam/adipose_1052_sup_modbases_GRCh38_20220729.bam
# index
samtools index ./bam/adipose_1052_sup_modbases_GRCh38_20220729.hp.bam
```

## modbam2bed

Aggregate modified base counts stored in a modified-base BAM file to a bedMethyl file

Get software

```bash
# add epi2melabs conda channel in order to get software
conda config --add channels epi2melabs
# create conda env with modbam2bed installed
mamba create -n modbam2bed -c bioconda -c conda-forge -c epi2melabs modbam2bed=0.6.2
conda activate modbam2bed
```

Run

```bash
for HP in 1 2; do
    modbam2bed \
        -e \
        -m 5mC \
        --cpg \
        -t 10 \
        --haplotype ${HP} \
        /data/publicData/genomes/human/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
        ./bam/adipose_1052_sup_modbases_GRCh38_20220729.hp.bam \
        | bgzip -c > ./bed/adipose_1052_bonito.hp${HP}.cpg.bed.gz
done;
```
