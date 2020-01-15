#!/bin/bash
tr="1";
for strc in `cat eterna.csv | cut -f 4 -d ','` ; 

do 
	if [ "$strc" != "$tr" ]; 
	then 
	 `python rnaevol.py --target="$strc" --job=1 -g 500 -n 100 >> eterna_net.txt` ;
	fi
done
