#!/bin/bash

echo "Converting Psudo tables into biom files"
biom convert -i ./pseudo/mucus/mucus_combined_pseudo_table.tsv -o ./pseudo/mucus/mucus_pseudo_feature_table_hdf5.biom --table-type="OTU table" --to-hdf5
biom convert -i ./pseudo/tissue/tissue_combined_pseudo_table.tsv -o ./pseudo/tissue/tissue_pseudo_feature_table_hdf5.biom --table-type="OTU table" --to-hdf5
biom convert -i ./pseudo/skeleton/skeleton_combined_pseudo_table.tsv -o ./pseudo/skeleton/skeleton_pseudo_feature_table_hdf5.biom --table-type="OTU table" --to-hdf5

echo "Converting psudo tables into biom files"
qiime tools import   --input-path ./pseudo/tissue/tissue_pseudo_feature_table_hdf5.biom   --type 'FeatureTable[Frequency]'   --input-format BIOMV210Format  --output-path ./pseudo/tissue/tissue_feature_table.qza
qiime tools import   --input-path ./pseudo/skeleton/skeleton_pseudo_feature_table_hdf5.biom   --type 'FeatureTable[Frequency]'   --input-format BIOMV210Format  --output-path ./pseudo/skeleton/skeleton_feature_table.qza
qiime tools import   --input-path ./pseudo/mucus/mucus_pseudo_feature_table_hdf5.biom   --type 'FeatureTable[Frequency]' --input-format BIOMV210Format  --output-path ./pseudo/mucus/mucus_feature_table.qza

echo "finished"

