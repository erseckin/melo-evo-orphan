#!/bin/bash

LIN=/path/to/metazoa_odb10 


for n in *fa
do
       busco -i ${n} -o ${n}_busco_metazoa_odb10.out -l $LIN -m prot
done



