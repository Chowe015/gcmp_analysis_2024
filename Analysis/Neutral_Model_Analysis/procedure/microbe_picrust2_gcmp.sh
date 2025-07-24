#!/bin/bash

#import data into fastp
echo "Picrust2 pipeline for GCMP Neutral Psuedo table Analysis..."
echo "How many cpu would you like to use?"
read speed

for i in {mucus,tissue,skeleton}; do
        echo $i
        
	 R1=$(ls ${i}/${i}_psudo_feature_table_hdf5.biom);

### Generate metagenome predictions
       picrust2_pipeline.py -s ./merged_sequence.fasta -i $R1 -o ${i}_picrust2_out_pipeline -p $speed


done

echo "finished... at $date..."
