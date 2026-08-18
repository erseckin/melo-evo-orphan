#!/bin/bash

DB=/path/to/NR_diamond_taxo.dmnd
QUERY=query.fa
OUT=output.tsv

diamond blastp --more-sensitive -e 0.0001 --id 25 -p16 -d $DB -q $QUERY -o $OUT -f 102 --top 100
