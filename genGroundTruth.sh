#!/usr/bin/env bash

output_directed=false

print_usage () {
    echo "Usage: genGroundTruth.sh [-d] <min_ovlp> <sample-basename>"
    echo "Example: genGroundTruth.sh 15 SRX533603-sample"
    exit 1
}

if [[ $# -eq 2 ]]
then
    min_ovlp="$1"
    dataset_name="$2"
elif [[ $# -eq 3 ]]
then
    if [[ $1 = "-d" ]]
    then
        output_directed=true
        min_ovlp="$2"
        dataset_name="$3"
    else
	print_usage
    fi
else
    print_usage
fi

sample_dir="sample-reads"
suffix=".efa"
sample_file="$sample_dir/$dataset_name$suffix"

if [[ -z "$sample_file" ]]
then
	echo "Cannot find sample file in $sample_dir/"
	exit 1
fi

awk '($1 ~ ">"){print $0}' "$sample_file" | sort -t$'\t' -k6,6 -k7,7n -k8,8n -o "$sample_file"-sorted
if [[ "$output_directed" = true ]]
then
    k8 ./genGroundTruth.js -d -l "$min_ovlp" "$sample_file"-sorted
    sort -k1,1n -k2,2n -o "$sample_dir/$dataset_name".truepairs-directed "$sample_file"-sorted.pairs
else
    k8 ./genGroundTruth.js -l "$min_ovlp" "$sample_file"-sorted
    sort -k1,1n -k2,2n -o "$sample_dir/$dataset_name".truepairs "$sample_file"-sorted.pairs
fi
rm "$sample_file"-sorted "$sample_file"-sorted.pairs


