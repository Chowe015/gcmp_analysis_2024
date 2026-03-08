#!/bin/bash
echo "Converting Pathway metagenome predictions"
biom convert -i ./final_picrust2/s_disp_m_picrust2_out_pipeline/pathway_out/path_abun_unstrat.tsv -o ./final_picrust2/s_disp_m_picrust2_out_pipeline/pathway_out/s_disp_m_pathway_table.biom --table-type="OTU table" --to-hdf5
biom convert -i ./final_picrust2/s_disp_t_picrust2_out_pipeline/pathway_out/path_abun_unstrat.tsv -o ./final_picrust2/s_disp_t_picrust2_out_pipeline/pathway_out/s_disp_t_pathway_table.biom --table-type="OTU table" --to-hdf5
biom convert -i ./final_picrust2//t_disp_m_picrust2_out_pipeline/pathway_out/path_abun_unstrat.tsv -o ./final_picrust2/t_disp_m_picrust2_out_pipeline/pathway_out/t_disp_m_pathway_table.biom --table-type="OTU table" --to-hdf5

echo "Converting biom files into phyloseq object for downstream analysis"
qiime tools import --input-path ./final_picrust2/s_disp_m_picrust2_out_pipeline/pathway_out/s_disp_m_pathway_table.biom \
                   --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path ./final_picrust2/s_disp_m_picrust2_out_pipeline/pathway_out/s_disp_m_picrust_pathway_table.qza
qiime tools import --input-path ./final_picrust2/s_disp_t_picrust2_out_pipeline/pathway_out/s_disp_t_pathway_table.biom \
                   --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path ./final_picrust2/s_disp_t_picrust2_out_pipeline/pathway_out/s_disp_t_picrust_pathway_table.qza
qiime tools import --input-path ./final_picrust2//t_disp_m_picrust2_out_pipeline/pathway_out/t_disp_m_pathway_table.biom \
                   --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path ./final_picrust2/t_disp_m_picrust2_out_pipeline/pathway_out/t_disp_m_picrust_pathway_table.qza

echo "Converting KO metagenome predictions into biom files"
biom convert -i ./final_picrust2/s_disp_m_picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv -o ./final_picrust2/s_disp_m_picrust2_out_pipeline/KO_metagenome_out/s_disp_m_KO_table.biom --table-type="OTU table" --to-hdf5
biom convert -i ./final_picrust2/s_disp_t_picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv -o ./final_picrust2/s_disp_t_picrust2_out_pipeline/KO_metagenome_out/s_disp_t_KO_table.biom --table-type="OTU table" --to-hdf5
biom convert -i ./final_picrust2/t_disp_m_picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv -o ./final_picrust2/t_disp_m_picrust2_out_pipeline/KO_metagenome_out/t_disp_m_KO_table.biom --table-type="OTU table" --to-hdf5


echo "Converting biom files into phyloseq object for downstream analysis"
qiime tools import --input-path ./final_picrust2/s_disp_m_picrust2_out_pipeline/KO_metagenome_out/s_disp_m_KO_table.biom \
                   --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path ./final_picrust2/s_disp_m_picrust2_out_pipeline/KO_metagenome_out/s_disp_m_KO_table.qza

qiime tools import --input-path ./final_picrust2/s_disp_t_picrust2_out_pipeline/KO_metagenome_out/s_disp_t_KO_table.biom \
                   --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path ./final_picrust2/s_disp_t_picrust2_out_pipeline/KO_metagenome_out/s_disp_t_KO_table.qza

qiime tools import --input-path ./final_picrust2/t_disp_m_picrust2_out_pipeline/KO_metagenome_out/t_disp_m_KO_table.biom \
                   --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path ./final_picrust2/t_disp_m_picrust2_out_pipeline/KO_metagenome_out/t_disp_m_KO_table.qza

echo "finished"


