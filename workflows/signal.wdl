version 1.0

workflow run_signal {
	input {
		String accession
		Array[File] fastq_R1s
		Array[File] fastq_R2s
		File scheme_bed
		File viral_reference_genome
		File viral_reference_feature_coords
		String viral_reference_contig_name
		File primer_pairs_tsv
		File amplicon_bed

		String container_registry
	}

	call signal {
		input:
			accession = accession,
			fastq_R1s = fastq_R1s,
			fastq_R2s = fastq_R2s,
			scheme_bed = scheme_bed,
			viral_reference_genome = viral_reference_genome,
			viral_reference_feature_coords = viral_reference_feature_coords,
			viral_reference_contig_name = viral_reference_contig_name,
			primer_pairs_tsv = primer_pairs_tsv,
			amplicon_bed = amplicon_bed,
			container_registry = container_registry
	}

	output {
		File ivar_vcf = signal.ivar_vcf
		File ivar_vcf_index = signal.ivar_vcf_index
		File ivar_assembly = signal.ivar_assembly
		File freebayes_vcf = signal.freebayes_vcf
		File freebayes_vcf_index = signal.freebayes_vcf_index
		File freebayes_assembly = signal.freebayes_assembly
		File summary = signal.summary
		File lineage_metadata = signal.lineage_metadata
		File bam = signal.bam
	}

	meta {
		author: "Heather Ward"
		email: "heather@dnastack.com"
	}
}

task signal {
	input {
		String accession
		Array[File] fastq_R1s
		Array[File] fastq_R2s
		File scheme_bed
		File viral_reference_genome
		File viral_reference_feature_coords
		String viral_reference_contig_name
		File primer_pairs_tsv
		File amplicon_bed

		String container_registry
	}

	Int threads = 8
	Int num_paired_fastqs = length(fastq_R1s)

	command <<<
		set -exEo pipefail

		ln -s /covid-19-signal/scripts/ "$(pwd)"

		viral_reference_contig_name='~{viral_reference_contig_name}'
		export viral_reference_contig_name

		# Generate config
		echo \
			"samples: 'sample_table.csv'
			result_dir: $(pwd)/output/
			min_qual: 20
			min_len: 20
			scheme_bed: '~{scheme_bed}'
			composite_reference: '$COMPOSITE_REFERENCE'
			viral_reference_contig_name: '~{viral_reference_contig_name}'
			viral_reference_genome: '~{viral_reference_genome}'
			viral_reference_feature_coords: '~{viral_reference_feature_coords}'
			run_breseq: False
			breseq_reference: ''
			run_freebayes: True
			kraken2_db: '$KRAKEN2_DB'
			primer_pairs_tsv: '-f ~{primer_pairs_tsv}'
			mpileup_depth: 100000
			var_freq_threshold: 0.75
			var_min_coverage_depth: 10
			var_min_freq_threshold: 0.25
			var_min_variant_quality: 20
			amplicon_loc_bed: '~{amplicon_bed}'
			phylo_include_seqs: ''
			negative_control_prefix: []
			pangolin_fast: False
			pangolin: '4.0.6'
			constellations: '0.1.9'
			scorpio: '0.3.17'
			pangolearn: ''
			pango-designation: ''
			pangolin-data: '1.8'
			nextclade-data: ''
			nextclade-include-recomb: True" \
		| tr -d '\t' > config.yaml

		count=0
		while [[ "$count" -lt ~{num_paired_fastqs} ]]; do
			echo "~{accession}" >> samplename.tmp
			set +e
			((count++))
			set -e
		done
		echo -e "~{sep='\n' fastq_R1s}" >> fastq_R1s.tmp
		echo -e "~{sep='\n' fastq_R2s}" >> fastq_R2s.tmp

		echo sample,r1_path,r2_path > sample_table.csv
		paste -d , samplename.tmp fastq_R1s.tmp fastq_R2s.tmp >> sample_table.csv && rm samplename.tmp fastq_R1s.tmp fastq_R2s.tmp

		conda run \
			-n snakemake \
			snakemake \
			--verbose \
			--use-conda \
			--conda-prefix "$HOME/.snakemake" \
			--cores ~{threads} \
			-s /covid-19-signal/Snakefile \
			all

		conda run \
			-n snakemake \
			snakemake \
			--verbose \
			--use-conda \
			--conda-prefix "$HOME/.snakemake" \
			--cores ~{threads} \
			-s /covid-19-signal/Snakefile \
			postprocess

		# iVar
		## Variants and index
		cp output/~{accession}/core/~{accession}_ivar_variants.tsv ./~{accession}.tsv
		ivar_variants_to_vcf.py ~{accession}.tsv ~{accession}.ivar.vcf
		bgzip ~{accession}.ivar.vcf && tabix ~{accession}.ivar.vcf.gz

		# having an issue with some corrupted VCFs causing variant transforms to hang
		# try reindexing until the index shows the proper block size, or exit
		set +e
		count=0
		while [ "$count" -lt 3 ]; do
			validate_index.sh \
				-v "~{accession}.ivar.vcf.gz" \
				-i "~{accession}.ivar.vcf.gz.tbi"
			retc=$?
			if [ "$retc" -ne 0 ]; then
				tabix -f "~{accession}.ivar.vcf.gz"
			else
				break
			fi
			count=$((count + 1))
		done
		if [ "$retc" -ne 0 ]; then
			echo "Failed to validate index"
			echo "Index dump:"
			bgzip -d < "~{accession}.ivar.vcf.gz.tbi" | xxd
			echo; echo "VCF dump:"
			xxd "~{accession}.ivar.vcf.gz"
			exit 1
		fi
		set -e

		## Assembly
		sed '1s/>\(.*\)/>~{accession} \1/' output/~{accession}/core/~{accession}.consensus.fa > ~{accession}.ivar.fa


		# Freebayes
		## Variants and index
		echo -e "unknown\t~{accession}" > samples.txt
		bcftools reheader \
			--samples samples.txt \
			output/~{accession}/freebayes/~{accession}.variants.norm.vcf \
		| bgzip > ~{accession}.freebayes.vcf.gz
		tabix ~{accession}.freebayes.vcf.gz

		set +e
		count=0
		while [ "$count" -lt 3 ]; do
			validate_index.sh \
				-v "~{accession}.freebayes.vcf.gz" \
				-i "~{accession}.freebayes.vcf.gz.tbi"
			retc=$?
			if [ "$retc" -ne 0 ]; then
				tabix -f "~{accession}.freebayes.vcf.gz"
			else
				break
			fi
			count=$((count + 1))
		done
		if [ "$retc" -ne 0 ]; then
			echo "Failed to validate index"
			echo "Index dump:"
			bgzip -d < "~{accession}.freebayes.vcf.gz.tbi" | xxd
			echo; echo "VCF dump:"
			xxd "~{accession}.freebayes.vcf.gz"
			exit 1
		fi
		set -e


		## Assembly
		cp output/~{accession}/freebayes/~{accession}.consensus.fasta ~{accession}.freebayes.fa


		# Lineage assignment
		# Only looking at the assignment for the iVar consensus sequence
		cut -f 2- output/lineage_assignments.tsv | sed 's/\t\([^\t]*,[^\t]*\)/\t"\1"/g' | tr '\t' ',' > ~{accession}.lineage_assignments.csv

		# Summary info
		mkdir ~{accession}_summary
		unzip output/summary.zip -d ~{accession}_summary/
		mv ~{accession}_summary/summary.html ~{accession}_summary/~{accession}/
		cp output/*.png ~{accession}_summary/
		cp output/lineage_assignments.tsv ~{accession}_summary/
		cp output/freebayes_lineage_assignments.tsv ~{accession}_summary/
		zip -r ~{accession}_summary.zip ~{accession}_summary/
	>>>

	output {
		File ivar_vcf = "~{accession}.ivar.vcf.gz"
		File ivar_vcf_index = "~{accession}.ivar.vcf.gz.tbi"
		File ivar_assembly = "~{accession}.ivar.fa"
		File freebayes_vcf = "~{accession}.freebayes.vcf.gz"
		File freebayes_vcf_index = "~{accession}.freebayes.vcf.gz.tbi"
		File freebayes_assembly = "~{accession}.freebayes.fa"
		File summary = "~{accession}_summary.zip"
		File lineage_metadata = "~{accession}.lineage_assignments.csv"
		File bam = "output/~{accession}/core/~{accession}_viral_reference.mapping.primertrimmed.sorted.bam"
	}

	runtime {
		docker: "~{container_registry}/signal:e6cae1e"
		cpu: threads
		memory: "32 GB"
		disks: "local-disk 500 HDD"
		bootDiskSizeGb: 30
	}
}
