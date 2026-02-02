#!/bin/bash
## important to remember combined_psudo_tables must contain #OTU ID in the first column to locate ASV with dna-sequences.fasta when running picrust2.
echo "Converting Psudo tables into biom files"
biom convert -i ./pseudo/mucus/M_T_combined_pseudo_table.tsv -o ./pseudo/mucus/t_disp_m_pseudo_feature_table_hdf5.biom --table-type="OTU table" --to-hdf5
biom convert -i ./pseudo/tissue/T_S_combined_pseudo_table.tsv -o ./pseudo/tissue/s_disp_t_pseudo_feature_table_hdf5.biom --table-type="OTU table" --to-hdf5
biom convert -i ./pseudo/mucus/M_S_combined_pseudo_table.tsv -o ./pseudo/mucus/s_disp_m_pseudo_feature_table_hdf5.biom --table-type="OTU table" --to-hdf5
#echo "Converting psudo tables into biom files"
qiime tools import --input-path ./pseudo/mucus/t_disp_m_pseudo_feature_table_hdf5.biom --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path ./pseudo/mucus/t_disp_m_pseudo_table
qiime tools import --input-path ./pseudo/tissue/s_disp_t_pseudo_feature_table_hdf5.biom --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path ./pseudo/tissue/s_disp_t_pseudo_table
qiime tools import --input-path ./pseudo/mucus/s_disp_m_pseudo_feature_table_hdf5.biom --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path ./pseudo/mucus/s_disp_m_pseudo_table
echo "finished"
