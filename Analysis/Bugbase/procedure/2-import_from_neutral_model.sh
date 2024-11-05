# Convert output table from Neutral Model for reclassification with Greengenes and then for use with Bugbase

# Request an interactive allocation.
salloc -N 1 -n 20 --mem-per-cpu=250 -t 12:00:00

# Load R module
module load r

# Run split tables using R script
cd /storage/home/yvl6147/scratch/GCMP/scripts

Rscript 1-split_neutral_model_outputs.R \
    /storage/home/yvl6147/scratch/GCMP/data/neutral_model_output/All/Mucus_Neutralmodel.csv \
    /storage/home/yvl6147/scratch/GCMP/data/neutral_model_output/All/Mucus_Neutralmodel_below.csv \
    /storage/home/yvl6147/scratch/GCMP/data/neutral_model_output/All/Mucus_Neutralmodel_neutral.csv \
    /storage/home/yvl6147/scratch/GCMP/data/neutral_model_output/All/Mucus_Neutralmodel_above.csv

Rscript 1-split_neutral_model_outputs.R \
    /storage/home/yvl6147/scratch/GCMP/data/neutral_model_output/All/Tissue_Neutralmodel.csv \
    /storage/home/yvl6147/scratch/GCMP/data/neutral_model_output/All/Tissue_Neutralmodel_below.csv \
    /storage/home/yvl6147/scratch/GCMP/data/neutral_model_output/All/Tissue_Neutralmodel_neutral.csv \
    /storage/home/yvl6147/scratch/GCMP/data/neutral_model_output/All/Tissue_Neutralmodel_above.csv

Rscript 1-split_neutral_model_outputs.R \
    /storage/home/yvl6147/scratch/GCMP/data/neutral_model_output/All/Skeleton_Neutralmodel.csv \
    /storage/home/yvl6147/scratch/GCMP/data/neutral_model_output/All/Skeleton_Neutralmodel_below.csv \
    /storage/home/yvl6147/scratch/GCMP/data/neutral_model_output/All/Skeleton_Neutralmodel_neutral.csv \
    /storage/home/yvl6147/scratch/GCMP/data/neutral_model_output/All/Skeleton_Neutralmodel_above.csv

Rscript 1-split_neutral_model_outputs.R \
    /storage/home/yvl6147/scratch/GCMP/data/neutral_model_output/All/all_compart_Neutralmodel.csv \
    /storage/home/yvl6147/scratch/GCMP/data/neutral_model_output/All/all_compart_Neutralmodel_below.csv \
    /storage/home/yvl6147/scratch/GCMP/data/neutral_model_output/All/all_compart_Neutralmodel_neutral.csv \
    /storage/home/yvl6147/scratch/GCMP/data/neutral_model_output/All/all_compart_Neutralmodel_above.csv

# Activate qiime2 environment
conda activate qiime2

# Extract all seq IDs from neutral model output tables. Then add a header row that says "feature-id".
cd /storage/home/yvl6147/scratch/GCMP/data/neutral_model_output/All
for i in *.csv; 
    do 
        awk -F "," '{print $1}' ${i} > ${i%.csv}_seqIDs.txt;
        sed -i '1d' ${i%.csv}_seqIDs.txt;
        sed -i 1i"feature-id" ${i%.csv}_seqIDs.txt;
        sed -i 's/"//g' ${i%.csv}_seqIDs.txt
    done

# Import split MST feature tables as Qiime2 artefacts.

### Convert csv to tsv
sed 's/,/\t/g' M_feature_table.csv > M_feature_table.tsv
sed 's/,/\t/g' T_feature_table.csv > T_feature_table.tsv
sed 's/,/\t/g' S_feature_table.csv > S_feature_table.tsv
### Convert tsv to biom
biom convert -i M_feature_table.tsv -o M_feature_table.biom --table-type="OTU table" --to-hdf5
biom convert -i T_feature_table.tsv -o T_feature_table.biom --table-type="OTU table" --to-hdf5
biom convert -i S_feature_table.tsv -o S_feature_table.biom --table-type="OTU table" --to-hdf5
### Import biom to qza
qiime tools import \
  --input-path M_feature_table.biom \
  --output-path M_feature_table.qza \
  --type "FeatureTable[Frequency]"
qiime tools import \
  --input-path T_feature_table.biom \
  --output-path T_feature_table.qza \
  --type "FeatureTable[Frequency]"
qiime tools import \
  --input-path S_feature_table.biom \
  --output-path S_feature_table.qza \
  --type "FeatureTable[Frequency]"

# Extract the respective lists of seq IDs (above, at and below neutral) from rep-seqs. Extract the same lists of seq IDs from feature table.
cd /storage/home/yvl6147/scratch/GCMP/data
mkdir Bugbase_input
mkdir Bugbase_input/rep-seqs
mkdir Bugbase_input/rep-seqs/Mucus
mkdir Bugbase_input/rep-seqs/Tissue
mkdir Bugbase_input/rep-seqs/Skeleton
mkdir Bugbase_input/feature-table
mkdir Bugbase_input/feature-table/Mucus
mkdir Bugbase_input/feature-table/Tissue
mkdir Bugbase_input/feature-table/Skeleton
cd /storage/home/yvl6147/scratch/GCMP/data/neutral_model_output/All

# Filter feature table and rep seqs for just mucus samples
qiime feature-table filter-samples \
    --i-table /storage/home/yvl6147/scratch/GCMP/data/nomito/merged_table.qza \
    --m-metadata-file /storage/home/yvl6147/scratch/GCMP/data/gcmp_complete_mapping2024_v1.txt \
    --p-where '[tissue_compartment]="M"' \
    --o-filtered-table /storage/home/yvl6147/scratch/GCMP/data/nomito/M_feature_table_bugbase.qza

# Filter feature table for just tissue samples
qiime feature-table filter-samples \
    --i-table /storage/home/yvl6147/scratch/GCMP/data/nomito/merged_table.qza \
    --m-metadata-file /storage/home/yvl6147/scratch/GCMP/data/gcmp_complete_mapping2024_v1.txt \
    --p-where '[tissue_compartment]="T"' \
    --o-filtered-table /storage/home/yvl6147/scratch/GCMP/data/nomito/T_feature_table_bugbase.qza

# Filter feature table for just skeleton samples
qiime feature-table filter-samples \
    --i-table /storage/home/yvl6147/scratch/GCMP/data/nomito/merged_table.qza \
    --m-metadata-file /storage/home/yvl6147/scratch/GCMP/data/gcmp_complete_mapping2024_v1.txt \
    --p-where '[tissue_compartment]="S"' \
    --o-filtered-table /storage/home/yvl6147/scratch/GCMP/data/nomito/S_feature_table_bugbase.qza

# Filter feature table and rep-seqs by above, below and neutral.
for i in Mucus_*_seqIDs.txt;
    do
        qiime feature-table filter-seqs \
            --i-data /storage/home/yvl6147/scratch/GCMP/data/nomito/merged_rep-seqs.qza \
            --o-filtered-data /storage/home/yvl6147/scratch/GCMP/data/Bugbase_input/rep-seqs/Mucus/${i%_seqIDs.txt}.qza \
            --m-metadata-file $i;
        qiime feature-table filter-features \
            --i-table /storage/home/yvl6147/scratch/GCMP/data/nomito/M_feature_table_bugbase.qza \
            --o-filtered-table /storage/home/yvl6147/scratch/GCMP/data/Bugbase_input/feature-table/Mucus/${i%_seqIDs.txt}.qza \
            --m-metadata-file $i;
    done

for i in Tissue_*_seqIDs.txt;
    do
        qiime feature-table filter-seqs \
            --i-data /storage/home/yvl6147/scratch/GCMP/data/nomito/merged_rep-seqs.qza \
            --o-filtered-data /storage/home/yvl6147/scratch/GCMP/data/Bugbase_input/rep-seqs/Tissue/${i%_seqIDs.txt}.qza \
            --m-metadata-file $i;
        qiime feature-table filter-features \
            --i-table /storage/home/yvl6147/scratch/GCMP/data/nomito/T_feature_table_bugbase.qza \
            --o-filtered-table /storage/home/yvl6147/scratch/GCMP/data/Bugbase_input/feature-table/Tissue/${i%_seqIDs.txt}.qza \
            --m-metadata-file $i;
    done

for i in Skeleton_*_seqIDs.txt;
    do
        qiime feature-table filter-seqs \
            --i-data /storage/home/yvl6147/scratch/GCMP/data/nomito/merged_rep-seqs.qza \
            --o-filtered-data /storage/home/yvl6147/scratch/GCMP/data/Bugbase_input/rep-seqs/Skeleton/${i%_seqIDs.txt}.qza \
            --m-metadata-file $i;
        qiime feature-table filter-features \
            --i-table /storage/home/yvl6147/scratch/GCMP/data/nomito/S_feature_table_bugbase.qza \
            --o-filtered-table /storage/home/yvl6147/scratch/GCMP/data/Bugbase_input/feature-table/Skeleton/${i%_seqIDs.txt}.qza \
            --m-metadata-file $i;
    done

# Do the same for the all_compart portion.

cd /storage/home/yvl6147/scratch/GCMP/data
mkdir Bugbase_input/rep-seqs/all_compart
mkdir Bugbase_input/feature-table/all_compart
cd /storage/home/yvl6147/scratch/GCMP/data/neutral_model_output/All
for i in all_compart_*_seqIDs.txt;
    do
        qiime feature-table filter-seqs \
            --i-data /storage/home/yvl6147/scratch/GCMP/data/nomito/merged_rep-seqs.qza \
            --o-filtered-data /storage/home/yvl6147/scratch/GCMP/data/Bugbase_input/rep-seqs/all_compart/${i%_seqIDs.txt}.qza \
            --m-metadata-file $i;
        qiime feature-table filter-features \
            --i-table /storage/home/yvl6147/scratch/GCMP/data/nomito/merged_table.qza \
            --o-filtered-table /storage/home/yvl6147/scratch/GCMP/data/Bugbase_input/feature-table/all_compart/${i%_seqIDs.txt}.qza \
            --m-metadata-file $i;
    done