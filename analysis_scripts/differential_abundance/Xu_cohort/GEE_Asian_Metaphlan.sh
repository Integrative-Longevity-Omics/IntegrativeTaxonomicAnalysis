#!/bin/bash -l

FLPATH="/restricted/projectnb/uh2-sebas/data/metagenomics/external_cohorts/Han_Chinese_Centenarian/processed_data/data_library/"
OUTDIR="/restricted/projectnb/uh2-sebas/analysis/metagenomics/ye_analyses/Output/Asian_Metaphlan/"
BIOMDIR="metaphlan4_renorm_chinese_phyloseq.10.23.2024.rds"
PHENODIR="/restricted/projectnb/uh2-sebas/data/metagenomics/external_cohorts/Han_Chinese_Centenarian/metadata/DifferentialAbundance_Supplement.xlsx"
OUTNAME="AsianMetaphlan_11062024"
METHOD="GEE"

# Define the parameters for the jobs
declare -a arg1_values=($FLPATH $FLPATH $FLPATH $FLPATH $FLPATH $FLPATH $FLPATH)
declare -a arg2_values=($OUTDIR $OUTDIR $OUTDIR $OUTDIR $OUTDIR $OUTDIR $OUTDIR)
declare -a arg3_values=($BIOMDIR $BIOMDIR $BIOMDIR $BIOMDIR $BIOMDIR $BIOMDIR $BIOMDIR)
declare -a arg4_values=("kingdom" "phylum" "class" "order" "family" "genus" "species")
declare -a arg5_values=($OUTNAME $OUTNAME $OUTNAME $OUTNAME $OUTNAME $OUTNAME $OUTNAME)
declare -a arg6_values=($METHOD $METHOD $METHOD $METHOD $METHOD $METHOD $METHOD)
declare -a arg7_values=($PHENODIR $PHENODIR $PHENODIR $PHENODIR $PHENODIR $PHENODIR $PHENODIR)

# Loop over the arrays and submit jobs
for i in "${!arg4_values[@]}"; do
    arg1=${arg1_values[$i]}
    arg2=${arg2_values[$i]}
    arg3=${arg3_values[$i]}
    arg4=${arg4_values[$i]}
    arg5=${arg5_values[$i]}
    arg6=${arg6_values[$i]}
    arg7=${arg7_values[$i]}

    # Create a job script for each job
    job_script="job_$i.sh"
    cat <<EOT > $job_script
#!/bin/bash -l

#$ -P uh2-sebas
#$ -l h_rt=12:00:00
#$ -N GEE_Asian_Metaphlan
#$ -j y

module load R/4.3.1

Rscript /restricted/projectnb/uh2-sebas/analysis/metagenomics/ye_analyses/Scripts/Differential_Abundance/Asian_Metaphlan/GEE_Asian_Metaphlan.R -f "$arg1" -o "$arg2" -b "$arg3" -t "$arg4" -n "$arg5" -m "$arg6" -p "$arg7"

EOT

    # Submit the job
    qsub $job_script
done
