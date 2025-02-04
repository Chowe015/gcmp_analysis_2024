##### STEP 5: BUGBASE #####

# Download Bugbase
wget https://github.com/knights-lab/BugBase/archive/refs/tags/v.0.1.0.tar.gz
# Uncompress
tar -xzvf v.0.1.0.tar.gz

# Modify ~/.bash_profile or ~/.bashrc with the following paths to your BugBase directory.

nano ~/.bash_profile
nano ~/.bashrc

# Add the following:
# export BUGBASE_PATH=/Path/to/my/BugBase
# export PATH=$PATH:$BUGBASE_PATH/bin

# Check install and install dependencies.
run.bugbase.r -h 
# Some R packages that run into errors here may be installed manually in RStudio into the same BugBase-v.0.1.0/R_lib directory.
# The biom package must be re-installed. To do this, first delete the existing biom folder in the R_lib directory, then run in RStudio:
download.file(url="https://cran.r-project.org/src/contrib/Archive/biom/biom_0.3.12.tar.gz", destfile="~/Documents/PhD/Projects/GCMP/BugBase/BugBase-v.0.1.0/R_lib/biom_0.3.12.tar.gz")
install.packages("~/Documents/PhD/Projects/GCMP/BugBase/BugBase-v.0.1.0/R_lib/biom_0.3.12.tar.gz", type="source", repos=NULL, lib="~/Documents/PhD/Projects/GCMP/BugBase/BugBase-v.0.1.0/R_lib")

# Convert feature table from .biom to .txt
biom convert -i otu_picked_output/indo_carib_table-with-taxonomy.biom -o otu_picked_output/indo_carib_table-with-taxonomy.txt --to-tsv

# Run BugBase:
cd /Users/yifanli/Documents/PhD/Projects/GCMP/data/neutral_model_GG/

## First, find the universal threshold value for the entire dataset (Mucus + Skeleton + Tissue).
### default traits
for i in all_compart_*/;
    do
        /Users/yifanli/Documents/PhD/Projects/GCMP/BugBase/BugBase-v.0.1.0/bin/run.bugbase.r -i ${i}final-table-with-taxonomy.txt -t 5 -o /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/${i%_biom/}_defaults -a;
    done
### custom traits from KEGG ID map
for i in all_compart_*/;
    do
        /Users/yifanli/Documents/PhD/Projects/GCMP/BugBase/BugBase-v.0.1.0/bin/run.bugbase.r -i ${i}final-table-with-taxonomy.txt -t 5 -o /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/${i%_biom/}_custom_traits -a -k -p M00175_Nitrogen_fixation_nitrogen_ammonia_,M00176_Sulfur_reduction_sulfate_H2S_,M00506_CheA_CheYBV_chemotaxis_two_component_regulatory_system_,M00453_QseC_QseB_quorum_sensing_two_component_regulatory_system_,M00513_LuxQN_CqsS_LuxU_LuxO_quorum_sensing_two_component_regulatory_system_;
    done

## Next, run through the individual above, neutral, and below tables using the threshold value found previously. Threshold values will be different for each trait (phenotype) predicted, so this code has to be run for each trait of interest individually.
TRAIT=Aerobic
mkdir /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/$TRAIT
for i in *Neutralmodel_*_biom/;
    do
        /Users/yifanli/Documents/PhD/Projects/GCMP/BugBase/BugBase-v.0.1.0/bin/run.bugbase.r -i ${i}final-table-with-taxonomy.txt -t 5 -o /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/$TRAIT/${i%_biom/}_$TRAIT -a -p $TRAIT -T 0.001;
    done

TRAIT=Anaerobic
mkdir /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/$TRAIT
for i in *Neutralmodel_*_biom/;
    do
        /Users/yifanli/Documents/PhD/Projects/GCMP/BugBase/BugBase-v.0.1.0/bin/run.bugbase.r -i ${i}final-table-with-taxonomy.txt -t 5 -o /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/$TRAIT/${i%_biom/}_$TRAIT -a -p $TRAIT -T 0.001;
    done

TRAIT=Contains_Mobile_Elements
mkdir /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/$TRAIT
for i in *Neutralmodel_*_biom/;
    do
        /Users/yifanli/Documents/PhD/Projects/GCMP/BugBase/BugBase-v.0.1.0/bin/run.bugbase.r -i ${i}final-table-with-taxonomy.txt -t 5 -o /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/$TRAIT/${i%_biom/}_$TRAIT -a -p $TRAIT -T 0.2;
    done

TRAIT=Facultatively_Anaerobic
mkdir /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/$TRAIT
for i in *Neutralmodel_*_biom/;
    do
        /Users/yifanli/Documents/PhD/Projects/GCMP/BugBase/BugBase-v.0.1.0/bin/run.bugbase.r -i ${i}final-table-with-taxonomy.txt -t 5 -o /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/$TRAIT/${i%_biom/}_$TRAIT -a -p $TRAIT -T 0.001;
    done

TRAIT=Forms_Biofilms
mkdir /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/$TRAIT
for i in *Neutralmodel_*_biom/;
    do
        /Users/yifanli/Documents/PhD/Projects/GCMP/BugBase/BugBase-v.0.1.0/bin/run.bugbase.r -i ${i}final-table-with-taxonomy.txt -t 5 -o /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/$TRAIT/${i%_biom/}_$TRAIT -a -p $TRAIT -T 0.2;
    done

TRAIT=Gram_Negative
mkdir /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/$TRAIT
for i in *Neutralmodel_*_biom/;
    do
        /Users/yifanli/Documents/PhD/Projects/GCMP/BugBase/BugBase-v.0.1.0/bin/run.bugbase.r -i ${i}final-table-with-taxonomy.txt -t 5 -o /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/$TRAIT/${i%_biom/}_$TRAIT -a -p $TRAIT -T 0.001;
    done

TRAIT=Gram_Positive
mkdir /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/$TRAIT
for i in *Neutralmodel_*_biom/;
    do
        /Users/yifanli/Documents/PhD/Projects/GCMP/BugBase/BugBase-v.0.1.0/bin/run.bugbase.r -i ${i}final-table-with-taxonomy.txt -t 5 -o /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/$TRAIT/${i%_biom/}_$TRAIT -a -p $TRAIT -T 0.001;
    done

TRAIT=Potentially_Pathogenic
mkdir /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/$TRAIT
for i in *Neutralmodel_*_biom/;
    do
        /Users/yifanli/Documents/PhD/Projects/GCMP/BugBase/BugBase-v.0.1.0/bin/run.bugbase.r -i ${i}final-table-with-taxonomy.txt -t 5 -o /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/$TRAIT/${i%_biom/}_$TRAIT -a -p $TRAIT -T 0.2;
    done

TRAIT=Stress_Tolerant
mkdir /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/$TRAIT
for i in *Neutralmodel_*_biom/;
    do
        /Users/yifanli/Documents/PhD/Projects/GCMP/BugBase/BugBase-v.0.1.0/bin/run.bugbase.r -i ${i}final-table-with-taxonomy.txt -t 5 -o /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/$TRAIT/${i%_biom/}_$TRAIT -a -p $TRAIT -T 0.6;
    done

TRAIT=M00175_Nitrogen_fixation_nitrogen_ammonia_
mkdir /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/$TRAIT
for i in *Neutralmodel_*_biom/;
    do
        /Users/yifanli/Documents/PhD/Projects/GCMP/BugBase/BugBase-v.0.1.0/bin/run.bugbase.r -i ${i}final-table-with-taxonomy.txt -t 5 -o /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/$TRAIT/${i%_biom/}_$TRAIT -a -k -p $TRAIT -T 0.9;
    done

TRAIT=M00176_Sulfur_reduction_sulfate_H2S_
mkdir /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/${TRAIT}
for i in *Neutralmodel_*_biom/;
    do
        /Users/yifanli/Documents/PhD/Projects/GCMP/BugBase/BugBase-v.0.1.0/bin/run.bugbase.r -i ${i}final-table-with-taxonomy.txt -t 5 -o /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/${TRAIT}/${i%_biom/}_$TRAIT -a -k -p $TRAIT -T 0.7;
    done

TRAIT=M00453_QseC_QseB_quorum_sensing_two_component_regulatory_system_
mkdir /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/$TRAIT
for i in *Neutralmodel_*_biom/;
    do
        /Users/yifanli/Documents/PhD/Projects/GCMP/BugBase/BugBase-v.0.1.0/bin/run.bugbase.r -i ${i}final-table-with-taxonomy.txt -t 5 -o /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/$TRAIT/${i%_biom/}_$TRAIT -a -k -p $TRAIT -T 0.6;
    done

TRAIT=M00506_CheA_CheYBV_chemotaxis_two_component_regulatory_system_
mkdir /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/$TRAIT
for i in *Neutralmodel_*_biom/;
    do
        /Users/yifanli/Documents/PhD/Projects/GCMP/BugBase/BugBase-v.0.1.0/bin/run.bugbase.r -i ${i}final-table-with-taxonomy.txt -t 5 -o /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/$TRAIT/${i%_biom/}_$TRAIT -a -k -p $TRAIT -T 0.9;
    done

TRAIT=M00513_LuxQN_CqsS_LuxU_LuxO_quorum_sensing_two_component_regulatory_system_
mkdir /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/$TRAIT
for i in *Neutralmodel_*_biom/;
    do
        /Users/yifanli/Documents/PhD/Projects/GCMP/BugBase/BugBase-v.0.1.0/bin/run.bugbase.r -i ${i}final-table-with-taxonomy.txt -t 5 -o /Users/yifanli/Documents/PhD/Projects/GCMP/Bugbase_output/$TRAIT/${i%_biom/}_$TRAIT -a -k -p $TRAIT -T 0.8;
    done


# -i is the input OTU table.
# -m is the mapping file. Not required. Use -a flag if not using a mapping file.
# -c is the treatment groups column header. Only required if mapping file is provided.
# -t tells it to classify taxonomy to the Family (5) level.
# -o is the output directory.
# -l is CLR transformed.
# -p are the specific traits (phenotypes) to predict, separated by commas, no spaces [default NULL].

## I then tried running Bugbase without CLR but then I get this error during the plotting threshold step:
# Error in trait_thresholds[, 1] <- as.numeric(thresholds[which.threshold]) : 
#   replacement has length zero
# Calls: single.cell.predictions

