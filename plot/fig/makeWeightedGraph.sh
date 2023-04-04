#!/usr/bin/env bash

echo "graph{"
echo "node [shape=\"point\"];"

ct=10

while read l
do
    readarray -td$'\t' r <<< $l
    name=${r[0]#>}
    ypos[$name]=$ct
    
    x=${r[1]:0:-1}.${r[1]: -1}
    echo r"$name"st" "[pos=\"$x,-$ct\"]\;

    y=${r[2]:0:-2}.${r[2]: -2:1}
    y=$(echo "scale=1; ${r[2]%$'\n'}/10" | bc)
    echo r"$name"ed" "[pos=\"$y,-$ct\"]\;

    x=$(echo "scale=1; ($x+$y)/2" | bc)
    echo r"$name"mid" "[pos=\"$x,-$ct\", style=\"invis\"]\;
    echo r"$name"st--r"$name"ed" "[label=\"$name\"]\;
    ((ct+=30))
done < ../SRX533603-pri-sample10k.sorted

echo "edge [color=\"red\", style=\"dashed\"];"
while read l
do
    readarray -td' ' r <<< $l
    echo r"${r[0]}"mid--r"${r[1]}"mid" "[label=\""${r[2]%$'\n'}"\"]\;
done < ../raw/overlap-n10000.all-pair

echo "}"
