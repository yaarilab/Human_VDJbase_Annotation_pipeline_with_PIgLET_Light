# AIRR-seq repertoire annotation and downstream analysis


Installation and usage guide are adapted from [Nextflow FLAIRR pipeline based on pipeAIRR](https://github.com/williamdlees/flairr_dsl2) written by William Lees

## Installation

The pipeline will run on Linux, or Windows Subsystem for Linux (WSL).

You will need to install nextflow as well as docker.

- [Nextflow](https://www.nextflow.io/)
- [Docker client](https://www.docker.com/) or [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html)


## Usage

- Clone this repo
- Review the configuration in nextflow.config and make any changes to number of threads, etc.
- Pull the relevant dockers
  1. peresay/suite:latest
  2. williamlees/ogrdbstats:latest
- Run the pipeline with the following line
 
```bash
sample="sample.fasta"
outDir=${pwd}
vRef="orgdb_germlineset/IGHV_asc.fasta" # ogrdb V reference set
dRef="orgdb_germlineset/IGHD.fasta" # ogrdb D reference set
jRef="orgdb_germlineset/IGHJ.fasta" # ogrdb J reference set
pigletThreshold="piglet_files/allele_threshold_table_ogrdb.tsv" # piglet allele based genotype threshold table

nextflow run main.nf --airr_seq ${sample} --v_germline ${vRef} --d_germline ${dRef} --j_germline ${jRef} --allele_threshold_table ${pigletThreshold} --outdir ${outDir}

# optional, a seperate work directory can be defined using -w flag
nextflow run main.nf --airr_seq ${sample} --v_germline ${vRef} --d_germline ${dRef} --j_germline ${jRef} --allele_threshold_table ${pigletThreshold} --outdir ${outDir} -w ${work}
```

## Pipeline Product

The final product of the pipeline will results in these folders:

- initial_annotation - contains the initial annotated repertoire to the input reference.
- pre_genotype - contains the sequences used to infer the genotype.
- genotype_report - contains the PIgLET genotype reports for V, D, and J inference.
- final_annotation - contains the annotated repertoire to the inferred personal reference.
- ogrdb_report - contains the OGRDB statistics reports on the final annotated repertoire.


## Annotation pipeline overview:

The pipeline performs anotation and downstrean analysis of AIRR-seq.

The pipeline can be devided into seven main componenet:

**1. Initial repertoire alignment and annotation based reference set**

> In this section, the repertoire sequences are annotated using IgBlast and MakeDb (presto) against the supplied reference set and collapse tham.

**2. Undocumented allele inference**

> In this section, inference for undocumented V allele (not found in the initial reference set) is performed using TIgGER.

**3. Second repertoire alignment and annotation with the discovered undocumented alleles.**

> In this section, in case undocumented alleles were inferred the repertoire is re-aligned and annotated with the additional alleles.

**4. Clonal inference and selection of colonal representative**

> In this section, clones are infered for the annotated repertoire, and a single representative for each clone whith the least number of mutation is chosen.

**5. Genotype inference**

> In this section, genotypes are inferred for each of the IG calls: V, D, and J using PIgLET inferernce tool, and a personal reference set is created.

**6. Third repertoire alignment and annotation with the personal reference set.**

> In this section, the repertoire sequences are annotated using IgBlast and MakeDb (presto) against the personal reference set from step 5.

**7. OGRDB statistics.**

> the reperotire statistics are deduced using ogrdbstats.


### Input files:

1. An AIRR-seq in fasta format.
2. heavy_chain (yes/no)
4. Reference set files for IGHV, IGHD, and IGHJ alleles in fasta format for the right chain.
5. auxiliary_data and custom_internal_data for the right chain.

### Output:

1. The aligned repertoire with the personal reference set.
2. The aligned repertoire with the init reference set.
3. A genotype reports
4. A log files
5. An ogrdb statistics report for the init aligned repertoire
6. An ogrdb statistics report for the final aligned repertoire
7. pipeline_statistic.csv - table of pass and fail reads for some of the steps.

### Docker images: 

The pipeline uses two docker images:

1. peresay/suite
2. williamlees/ogrdbstats
