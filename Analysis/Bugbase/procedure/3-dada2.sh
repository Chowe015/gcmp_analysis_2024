### step 5. import ref database

qiime tools import \
--input-path 97_otus.fasta \
--output-path gg_97_otus.qza \
--type 'FeatureData[Sequence]'

### step 6. assign taxonomy
cd /storage/home/yvl6147/scratch/GCMP/data/Bugbase_input/rep-seqs/Mucus
for i in *.qza;
    do
        mkdir /storage/home/yvl6147/scratch/GCMP/otu_picked_output/${i%.qza};
        qiime vsearch cluster-features-closed-reference \
            --i-table /storage/home/yvl6147/scratch/GCMP/data/Bugbase_input/feature-table/Mucus/$i \
            --i-sequences $i \
            --i-reference-sequences /storage/home/yvl6147/scratch/GCMP/Greengenes_ref/gg_97_otus.qza \
            --p-perc-identity 0.87 \
            --o-clustered-table /storage/home/yvl6147/scratch/GCMP/otu_picked_output/${i%.qza}/table-97.qza \
            --o-unmatched-sequences /storage/home/yvl6147/scratch/GCMP/otu_picked_output/${i%.qza}/unmatched.qza \
            --o-clustered-sequences /storage/home/yvl6147/scratch/GCMP/otu_picked_output/${i%.qza}/clustered-features.qza;
    done

cd /storage/home/yvl6147/scratch/GCMP/data/Bugbase_input/rep-seqs/Tissue
for i in *.qza;
    do
        mkdir /storage/home/yvl6147/scratch/GCMP/otu_picked_output/${i%.qza};
        qiime vsearch cluster-features-closed-reference \
            --i-table /storage/home/yvl6147/scratch/GCMP/data/Bugbase_input/feature-table/Tissue/$i \
            --i-sequences $i \
            --i-reference-sequences /storage/home/yvl6147/scratch/GCMP/Greengenes_ref/gg_97_otus.qza \
            --p-perc-identity 0.87 \
            --o-clustered-table /storage/home/yvl6147/scratch/GCMP/otu_picked_output/${i%.qza}/table-97.qza \
            --o-unmatched-sequences /storage/home/yvl6147/scratch/GCMP/otu_picked_output/${i%.qza}/unmatched.qza \
            --o-clustered-sequences /storage/home/yvl6147/scratch/GCMP/otu_picked_output/${i%.qza}/clustered-features.qza;
    done

cd /storage/home/yvl6147/scratch/GCMP/data/Bugbase_input/rep-seqs/Skeleton
for i in *.qza;
    do
        mkdir /storage/home/yvl6147/scratch/GCMP/otu_picked_output/${i%.qza};
        qiime vsearch cluster-features-closed-reference \
            --i-table /storage/home/yvl6147/scratch/GCMP/data/Bugbase_input/feature-table/Skeleton/$i \
            --i-sequences $i \
            --i-reference-sequences /storage/home/yvl6147/scratch/GCMP/Greengenes_ref/gg_97_otus.qza \
            --p-perc-identity 0.87 \
            --o-clustered-table /storage/home/yvl6147/scratch/GCMP/otu_picked_output/${i%.qza}/table-97.qza \
            --o-unmatched-sequences /storage/home/yvl6147/scratch/GCMP/otu_picked_output/${i%.qza}/unmatched.qza \
            --o-clustered-sequences /storage/home/yvl6147/scratch/GCMP/otu_picked_output/${i%.qza}/clustered-features.qza;
    done

cd /storage/home/yvl6147/scratch/GCMP/data/Bugbase_input/rep-seqs/all_compart
for i in *.qza;
    do
        mkdir /storage/home/yvl6147/scratch/GCMP/otu_picked_output/${i%.qza};
        qiime vsearch cluster-features-closed-reference \
            --i-table /storage/home/yvl6147/scratch/GCMP/data/Bugbase_input/feature-table/all_compart/$i \
            --i-sequences $i \
            --i-reference-sequences /storage/home/yvl6147/scratch/GCMP/Greengenes/gg_97_otus.qza \
            --p-perc-identity 0.97 \
            --o-clustered-table /storage/home/yvl6147/scratch/GCMP/otu_picked_output/${i%.qza}/table-97.qza \
            --o-unmatched-sequences /storage/home/yvl6147/scratch/GCMP/otu_picked_output/${i%.qza}/unmatched.qza \
            --o-clustered-sequences /storage/home/yvl6147/scratch/GCMP/otu_picked_output/${i%.qza}/clustered-features.qza;
    done

### step 7. export biom
cd /storage/home/yvl6147/scratch/GCMP/otu_picked_output/
for i in */;
    do
        qiime tools export \
            --input-path ${i}table-97.qza \
            --output-path ${i%/}_table-97_biom;
    done

### step 8. add taxonomy to biom

##### Observation metadata file must contain headers:
##### #OTUID    taxonomy
cd /storage/home/yvl6147/scratch/GCMP/otu_picked_output
for i in *_biom/;
    do
        biom add-metadata \
            -i ${i}feature-table.biom \
            -o ${i}final-table-with-taxonomy.biom \
            --observation-metadata-fp ../Greengenes_ref/97_otu_taxonomy.txt \
            --sc-separated taxonomy;
    done

### step 9. convert biom to tsv

for i in *_biom/;
    do
        biom convert -i ${i}final-table-with-taxonomy.biom -o ${i}final-table-with-taxonomy.txt --to-tsv;
    done

### step 10. transfer files to local computer to run Bugbase.

