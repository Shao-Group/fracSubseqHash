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
    x=${r[2]:0:-2}.${r[2]: -2:1}
    x=$(echo "scale=1; ${r[2]%$'\n'}/10" | bc)

    echo r"$name"ed" "[pos=\"$x,-$ct\"]\;
    echo r"$name"st--r"$name"ed" "[label=\"$name\"]\;
    ((ct+=30))
done < ../SRX533603-pri-sample10k.sorted

echo "node [shape=\"rnastab\", label=\"\", color=\"red\"];"
echo "edge [color=\"red\", style=\"dashed\"];"
ct=1
while read l
do
    readarray -td' ' r <<< $l
    y=${ypos[${r[0]}]}
    x=${r[4]:0:-1}.${r[4]: -1}
    echo seed"$ct"m1" "[pos=\"$x,-$y\"]\;
    y=${ypos[${r[1]}]}
    x=${r[5]:0:-2}.${r[5]: -2:1}
    echo seed"$ct"m2" "[pos=\"$x,-$y\"]\;

    echo seed"$ct"m1--seed"$ct"m2\;
    ((ct+=1))
done < ../overlap-pairs-seedloc.txt

echo "}"
