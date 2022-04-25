#!/bin/bash

set -e

function usage() {
cat << EOF

  Check that the block end position specified by the index file is equal to the block size specified by the VCF file
  Usage: $0 -v vcf_file -i vcf_index_file

  OPTIONS
    -h    Display this message and exit
    -v    VCF file
    -i    VCF index file (TBI)

EOF
}

while getopts "hv:i:" OPTION; do
  case $OPTION in
    h) usage; exit;;
    v) VCF="$OPTARG";;
    i) VCF_INDEX="$OPTARG";;
    \?) usage; exit 1;;
  esac
done

if [ ! -s "${VCF}" ]; then
  usage
  echo "[ERROR] Must provide VCF file"
  exit 1
fi

if [ ! -s "${VCF_INDEX}" ]; then
  usage
  echo "[ERROR] Must provide VCF index file"
  exit 1
fi

vcf_bgzipped_size=$(stat --format="%s" "${VCF}")
vcf_bgzipped_size_no_eof_decimal=$(($vcf_bgzipped_size - 28))
vcf_bgzipped_size_no_eof_hex="$(printf '%012x\n' "${vcf_bgzipped_size_no_eof_decimal}" | sed 's/.\{2\}/& /g;s/ $//' | awk '{print $6$5$4$3$2$1}')"

total_bgzipped_size_in_index=$(cat "${VCF_INDEX}" | bgzip -d | xxd -p | tr -d '\n' | grep "${vcf_bgzipped_size_no_eof_hex}" || [[ $? == 1 ]])

if [ -z "${total_bgzipped_size_in_index}" ]; then
  echo "[ERROR] Expected bgzipped VCF file length [${vcf_bgzipped_size_no_eof_decimal} (${vcf_bgzipped_size_no_eof_hex})] not found in index file"
  exit 1
else
  echo "Index valid"
  exit 0
fi
