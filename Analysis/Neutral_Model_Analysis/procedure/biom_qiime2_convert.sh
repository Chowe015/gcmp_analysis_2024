#!/bin/bash

echo "Converting Psudo tables into biom files"
biom convert -i ./mucus/mucus_combined_psudo_table.tsv -o ./mucus/mucus_psudo_feature_table_hdf5.biom --table-type="OTU table" --to-hdf5
biom convert -i ./tissue/tissue_combined_psudo_table.tsv -o ./tissue/tissue_psudo_feature_table_hdf5.biom --table-type="OTU table" --to-hdf5
biom convert -i ./skeleton/skeleton_combined_psudo_table.tsv -o ./skeleton/skeleton_psudo_feature_table_hdf5.biom --table-type="OTU table" --to-hdf5

echo "Converting psudo tables into biom files"
qiime tools import   --input-path ./tissue/tissue_psudo_feature_table_hdf5.biom   --type 'FeatureTable[Frequency]'   --input-format BIOMV210Format  --output-path ./tissue/tissue_feature_table.qza
qiime tools import   --input-path ./skeleton/skeleton_psudo_feature_table_hdf5.biom   --type 'FeatureTable[Frequency]'   --input-format BIOMV210Format  --output-path ./skeleton/skeleton_feature_table.qza
qiime tools import   --input-path ./mucus/mucus_psudo_feature_table_hdf5.biom   --type 'FeatureTable[Frequency]' --input-format BIOMV210Format  --output-path ./mucus/mucus_feature_table.qza

echo "finished"

