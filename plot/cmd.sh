#!/usr/bin/bash

while read l; do
    readarray -td' ' r <<< $l;

    x=$(grep -m1 "^>${r[0]} " SRX533603-pri-10k-loc.txt);
    readarray -td' ' a <<< $x;
    x=$(grep -m1 "^>${r[1]} " SRX533603-pri-10k-loc.txt);
    readarray -td' ' b <<< $x;

    if [ ${a[3]} -le ${b[3]} ]; then
	if [ ${a[4] -le ${b[4]} ]; then
	       

    break;
done < overlap-pairs-seedloc.txt
