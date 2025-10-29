#!/bin/bash

echo "Converting Psudo tables into biom files"
biom convert -i ./mucus_picrust2_out_pipeline/mucus_KO_pred_metagenome_unstrat.tsv -o ./mucus_picrust2_out_pipeline/mucus_KO_table.biom --table-type="OTU table" --to-hdf5
biom convert -i ./tissue_picrust2_out_pipeline/tissue_KO_pred_metagenome_unstrat.tsv -o ./tissue_picrust2_out_pipeline/tissue_KO_table.biom --table-type="OTU table" --to-hdf5
biom convert -i ./skeleton_picrust2_out_pipeline/skeleton_KO_pred_metagenome_unstrat.tsv -o ./skeleton_picrust2_out_pipeline/skeleton_KO_table.biom --table-type="OTU table" --to-hdf5

echo "Converting psudo tables into biom files"
qiime tools import   --input-path ./tissue_picrust2_out_pipeline/tissue_KO_table.biom   --type 'FeatureTable[Frequency]'   --input-format BIOMV210Format  --output-path ./tissue_picrust2_out_pipeline/tissue_picrust_KO_table.qza
qiime tools import   --input-path ./skeleton_picrust2_out_pipeline/skeleton_KO_table.biom --type 'FeatureTable[Frequency]'   --input-format BIOMV210Format  --output-path ./skeleton_picrust2_out_pipeline/skeleton_picrust_KO_table.qza
qiime tools import   --input-path ./mucus_picrust2_out_pipeline/mucus_KO_table.biom   --type 'FeatureTable[Frequency]' --input-format BIOMV210Format  --output-path ./mucus_picrust2_out_pipeline/mucus_picrust_KO_table.qza

echo "finished"

