#!/usr/bin/env bash

min_ovlp=15

if [[ $# -eq 1 ]]
then
    dataset_name="$1"
elif [[ $# -eq 2 ]]
then
    min_ovlp="$1"
    dataset_name="$2"
else
    echo "Usage: genGroundTruth.sh [min_ovlp] <sample-basename>"
    echo "Example: genGroundTruth.sh 15 SRX533603-sample"
    exit 1
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
k8 ./genGroundTruth.js -l "$min_ovlp" "$sample_file"-sorted
sort -k1,1n -k2,2n -o "$sample_dir/$dataset_name".truepairs "$sample_file"-sorted.pairs
rm "$sample_file"-sorted "$sample_file"-sorted.pairs


