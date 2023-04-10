#!/usr/bin/env bash

dir="sample-reads"

if [[ $# -ne 2 && $# -ne 3 ]]
then
    echo "Usage: drawStringGraph.sh <sample-basename> <overlap-basename> [overlap-extension]"
    echo "Example: drawStringGraph.sh sim-e-coli-pb-le20k-nosub sim-e-coli-pb-le20k-nosub-n60-k42-t0/unionPos 10"
    echo "Will use files $dir/sim-e-coli-pb-le10k-nosub.{header-sorted|irr-edges} and $dir/sim-e-coli-pb-le20k-nosub-n60-k42-t0/unionPos.{all-pair|correct|irr-edges}10"
    exit 1
fi

#check all needed files exist
header="$dir/$1.header-sorted"
if [[ ! -f $header ]]
then
    echo "cannot read node file: $header" >&2
    exit 1
fi

all_irr="$dir/$1.irr-edges"
if [[ ! -f $all_irr ]]
then
    echo "cannot read all irreducible file: $all_irr" >&2
    exit 1
fi

found_prefix="$dir/$2"
found_irr="$found_prefix.irr-edges"
found_correct="$found_prefix.correct"
found_all="$found_prefix.all-pair"
if [[ $# -eq 3 ]]
then
    found_irr=$found_irr$3
    found_correct=$found_correct$3
    found_all=$found_all$3
fi

if [[ ! -f $found_irr ]]
then
    echo "cannot read found irreducible file: $found_irr" >&2
    exit 1
fi
if [[ ! -f $found_correct ]]
then
    echo "cannot read found correct file: $found_correct" >&2
    exit 1
fi
if [[ ! -f $found_all ]]
then
    echo "cannot read found all file: $found_all" >&2
    exit 1
fi

#generate dot file, from dot to svg
#dot -Knop2 -Tsvg xxx.dot > xxx.svg
echo 'graph{'
echo 'graph [splines=line];'
echo 'node [width=1.2, height=.1];'

#store y-coordinate of all the reads for later use
loc=()
#plot nodes from the sorted header file
cur_y=0
step=-80
while read name remain
do
    name=${name#>}
    echo "r$name [pos=\"80,$cur_y\", label=\"$name\"];"
    loc[$name]=$cur_y
    ((cur_y+=step))
done <$header

echo 'edge [penwidth=2.0];'
#plot all irreducible edges -- found:green, missed:black
found=false;
while read -u3 read1 read2 len
do
    if [[ "$found" = false ]]
    then
	if read -u4 fread1 fread2 len weight
	then
	    found=true
	fi
    fi

    if [[ "$found" = true && $read1 -eq $fread1 && $read2 -eq $fread2 ]]
    then
	echo "r$read1 -- r$read2 [color=\"lime\", fontcolor=\"limegreen\", label=\"$weight\"];"
	found=false
    else
	echo "r$read1 -- r$read2 [color=\"black\"];"
    fi
done 3<$all_irr 4<$found_irr

#plot all transitive edges -- i.e., correct \setminus irr-edges
found=false;
while read -u3 read1 read2 len weight
do
    if [[ "$found" = false ]]
    then
	if read -u4 fread1 fread2 remain
	then
	    found=true
	fi
    fi

    if [[ "$found" = true && $read1 -eq $fread1 && $read2 -eq $fread2 ]]
    then
	#is an irreducible edge, has been plotted
	found=false
    else
	((cur_y=(loc[$read1]+loc[$read2])/2))
	if [[ ${loc[$read1]} -gt ${loc[$read2]} ]]
	then
	    ((cur_x=140-loc[$read1]+loc[$read2]))
	else
	    ((cur_x=140+loc[$read1]-loc[$read2]))
	fi
	echo "r$read1""midr$read2 [pos=\"$cur_x,$cur_y\", shape=\"plaintext\", fontcolor=\"blue\", label=\"$weight\"];"
	echo "r$read1:w -- r$read1""midr$read2:c [color=\"cyan\"];"
	echo "r$read1""midr$read2:c -- r$read2:w [color=\"cyan\"];"
    fi
done 3<$found_correct 4<$found_irr

#plot all wrong edges
found=false;
while read -u3 read1 read2 weight
do
    if [[ "$found" = false ]]
    then
	if read -u4 fread1 fread2 remain
	then
	    found=true
	fi
    fi

    if [[ "$found" = true && $read1 -eq $fread1 && $read2 -eq $fread2 ]]
    then
	#is a correct edge, has been plotted
	found=false
    else
	((cur_y=(loc[$read1]+loc[$read2])/2))
	if [[ ${loc[$read1]} -gt ${loc[$read2]} ]]
	then
	    ((cur_x=80+(loc[$read1]-loc[$read2])/40))
	else
	    ((cur_x=80-(loc[$read1]-loc[$read2])/40))
	fi
	echo "r$read1""midr$read2 [pos=\"$cur_x,$cur_y\", shape=none, fontcolor=\"red\", label=\"$weight\"];"
	echo "r$read1:e -- r$read1""midr$read2:c [color=\"red\"];"
	echo "r$read1""midr$read2:c -- r$read2:e [color=\"red\"];"
    fi
done 3<$found_all 4<$found_correct

echo '}' #end of graph
