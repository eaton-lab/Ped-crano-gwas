#!/usr/bin/env bash
#
# trimmed for low quality bases and adapters using default settings in fastp v0.24.0
# (Chen, 2023), and mapped to the P. cranolopha "Liuliu" reference genome assembly
# (Xu et al., 2022) with BWA-MEM2 (Vasimuddin et al., 2019) while adding read groups
# for downstream consideration. Bam files were sorted with samtools v1.21 (Danecek et al., 2021

set -euo pipefail

# get sample name and group
#DATA="/media/hippo5/All_Pedicularis_WGS/2024-06-21_WGS_Pedicularis-106-10G/usftp21.novogene.com/01.RawData/S41_49/*.gz"
DATA=$1
REF=$2
THREADS=$3
OUTDIR=$4

# get read1 and read2 files
shopt -s nullglob
# Temporary associative array
declare -A r1_files
declare -A r2_files

# Loop over all gzipped fastq files
for file in $DATA; do
	fname=$(basename "$file")
	echo "Checking: $fname"
	# Match *_1.fq.gz or *_2.fq.gz or *_1.fastq.gz or *_2.fastq.gz
	if [[ "$fname" =~ ^(.+)_([12])\.f(ast)?q\.gz$ ]]; then
		prefix="${BASH_REMATCH[1]}"
		readnum="${BASH_REMATCH[2]}"
		ext=".fq.gz"
		[[ "${BASH_REMATCH[3]}" == "ast" ]] && ext=".fastq.gz"
		key="${prefix}${ext}"

		echo "  Match! key=$key readnum=$readnum"

		if [[ "$readnum" == "1" ]]; then
			r1_files["$key"]="$file"
		elif [[ "$readnum" == "2" ]]; then
			r2_files["$key"]="$file"
		fi
	else
		echo " No match"
	fi
done

# ensure ref is indexed
bwa-mem2 index "$REF"

# get unique keys
all_keys=("${!r1_files[@]}" "${!r2_files[@]}")
unique_keys=($(printf "%s\n" "${all_keys[@]}" | sort -u))
pair_id=0
declare -A bams
for key in "${unique_keys[@]}"; do
        ((pair_id++))
	r1="${r1_files[$key]}"
	r2="${r2_files[$key]}"
	NAME=$(basename "$(dirname "$r1")")
	GROUP=$(echo "$r1" | tr '/' '\n' | grep '^20' | head -n1)
	ID="${NAME}_${pair_id}"
	RG="@RG\tID:${ID}\tSM:${NAME}\tLB:${GROUP}\tPL:ILLUMINA"

	if [[ -n "$r1" && -n "$r2" ]]; then
        	echo -e "${RG}"
	elif [[ -n "$r1" ]]; then
        	>&2 echo "❌ Error: Found R1 file with no matching R2: $r1"
	elif [[ -n "$r2" ]]; then
        	>&2 echo "❌ Error: Found R2 file with no matching R1: $r2"
	fi

	# call methods
	fastp -i "${r1}" -I "${r2}" --stdout -w "${THREADS}" -j "${OUTDIR}/${ID}.fastp.json" --detect_adapter_for_pe \
	  | bwa-mem2 mem -p -t "${THREADS}" -R "${RG}" "$REF" - | \
	  | samtools sort -o "${OUTDIR}/${ID}.sorted.tmp.bam"

	# append bam to the array
	bams+=("${OUTDIR}/${ID}.sorted.tmp.bam")
done

# merge samples
samtools merge -@ "${THREADS}" -o "${OUTDIR}/${NAME}.merged.tmp.bam" "${bams[@]}"
samtools sort -@ "${THREADS}" -o "${OUTDIR}/${NAME}.sorted.bam" "${OUTDIR}/${NAME}.merged.tmp.bam"
samtools index "${OUTDIR}/${NAME}.sorted.bam"

# remove tmp bam files
rm -f "${OUTDIR}/*.tmp.bam"

# echo statement that final bam file is done.
echo "Final BAM file ready: ${NAME}.sorted.bam"
