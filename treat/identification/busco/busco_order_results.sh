#!/bin/bash

# Usage: order results according to busco scores 

for i in *odb10.out # BUSCO output directories
do
	cd $i
	grep -r "C:" short_summary*.txt | cut -c 4-7 >> ../buscoscore.txt
	cd ..
done

for i in *odb10.out
do
	echo "$i" >> name.txt
done

cat buscoscore.txt
cat name.txt
paste buscoscore.txt name.txt > buscoworm.txt

rm name.txt
rm buscoscore.txt

sort -k 2 -n -r buscoworm.txt > buscoscore_wormbase.txt

rm buscoworm.txt
