# Illumina SARS-CoV-2 data processing using the [SIGNAL pipeline](https://github.com/jaleezyy/covid-19-signal)

This repository provides a WDL wrapper for running the [SIGNAL pipeline](https://github.com/jaleezyy/covid-19-signal) to process Illumina paired-end SARS-CoV-2 sequencing data.


## Workflow inputs

An input template file with some defaults pre-defined can be found [here](./workflows/inputs.json).

| Input | Description |
|:-|:-|
| `accession` | Sample ID |
| `fastq_R1s`, `fastq_R2s` | Array of paired FASTQ file locations; paired files should be at the same index in each array |
| `scheme_bed` | The BED-format primer scheme used to prepare the library |
| `viral_reference_genome` | [The SARS-CoV-2 reference genome](https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3) |
| `viral_reference_feature_coords` | [Feature coordinates for the SARS-CoV-2 reference genome](https://storage.googleapis.com/dnastack-data-ingestion-storage/resources/MN908947.3.gff3) |
| `viral_reference_contig_name` | [`MN908947.3`] |
| `primer_pairs_tsv` | Primer pair TSV file; used for iVar's amplicon filter. This file is a headerless TSV containing one row per primer pair, with the LEFT primer names in column 1 and the RIGHT column names in column 2. |
| `amplicon_bed` | BED-formatted amplicon locations |
| `container_registry` | Registry that hosts workflow containers. All containers are hosted in [DNAstack's Dockerhub](https://hub.docker.com/u/dnastack) [`dnastack`] |


### Primer schemes

Primer schemes will differ based on the protocols used by the sequencing lab. Some common schemes can be downloaded from the official [artic-network github](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019). Additional schemes can be found in the [SIGNAL repository](https://github.com/jaleezyy/covid-19-signal/tree/master/resources/primer_schemes). Primers from these locations map to the inputs as follows:

- `scheme_bed`: ends in `.primer.bed`
- `amplicon_bed`: ends in `.scheme.bed`
- `primer_pairs_tsv`: this file is not provided directly, but can be generated from the `amplicon_bed` file

Example command to generate `primer_pairs_tsv` using the [ARTIC V3 scheme bed](https://github.com/artic-network/artic-ncov2019/blob/master/primer_schemes/nCoV-2019/V3/nCoV-2019.scheme.bed):

```bash
paste \
	<(cut -f 4 nCoV-2019.scheme.bed | sort -t _ -k 2 -g | grep LEFT) \
	<(cut -f 4 nCoV-2019.scheme.bed | sort -t _ -k 2 -g | grep RIGHT) \
> nCoV-2019.primer_pairs.tsv
```


## Workflow outputs

| Output | Description |
|:-|:-|
| `ivar_vcf`, `ivar_vcf_index` | Variants and index output by iVar |
| `ivar_assembly` | Genome assembly generated by iVar |
| `freebayes_vcf`, `freebayes_vcf_index` | Variants and index output by Freebayes |
| `freebayes_assembly` | Genome assembly generated by Freebayes |
| `summary` | Pipeline metrics |
| `lineage_metadata` | Pangolin lineage assignment metadata |
| `bam` | Reads aligned to the SARS-CoV-2 reference genome |


## Containers

Docker image definitions can be found in [dockerfiles](./dockerfiles).

The pipeline will always be pegged to a specific SIGNAL commit hash to avoid breaking due to updates to SIGNAL. The pipeline is periodically updated to the most recent version of the SIGNAL pipeline.

All containers are publicly hosted in [DNAstack's container registry](https://hub.docker.com/u/dnastack).

N.B. that the SIGNAL Docker container is ~10 GB to allow it to be used at scale in AWS, where EBS auto-scaling can sometimes not expand rapidly enough to accomodate running hundreds of samples in parallel. Including reference data in the Docker container seems to solve this issue, but does make it somewhat unwieldly.