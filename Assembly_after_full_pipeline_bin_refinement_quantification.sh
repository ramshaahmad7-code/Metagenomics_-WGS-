#!/bin/bash
#SBATCH -t 160:00:00
#SBATCH -J assembly
#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --mail-user=gm.grimaud@gmail.com
#SBATCH --mail-type=END,FAIL
#SBATCH -A lu2024-2-85



### ----------------------------------------------------------------------------------------
### General configurations
##Activate conda environment
module load Anaconda3/2024.02-1
# conda activate /home/ggrimaud/amr1

dir2="/home/ggrimaud/project"
util="/home/ggrimaud/util"
dir3="/lunarc/nobackup/projects/lu2024-12-47/shotgun1_HC"
metaW="/home/ggrimaud/util/metaWRAP/bin"
util2="/lunarc/nobackup/projects/lu2024-12-47/util"
dir="/lunarc/nobackup/projects/lu2024-12-47/Assembly/test"
KDdir="/lunarc/nobackup/projects/lu2024-12-47/Assembly/test/KNEADDATA/unpaired"
MAGSdir="/lunarc/nobackup/projects/lu2024-12-47/Assembly/test/KNEADDATA/unpaired/MAGS"


cd $dir/data

for sra_id in $(cat sra_ids_leukemia_control_leukemia_v1_save); do
    
    conda activate /home/ggrimaud/amr1
 conda activate /home/ggrimaud/amr1

    ###Downloading data######################
    cd $dir/data
    echo "Processing: $sra_id"
    #Download SRA file
    $util/sratoolkit.3.1.1-ubuntu64/bin/prefetch "$sra_id"
    echo "Downloaded sample $sra_id"

    #Uncompress SRA file using fastq-dump
    $util/sratoolkit.3.1.1-ubuntu64/bin/fasterq-dump "$sra_id"
    echo "uncompress SRA from $sra_id"
    rm -r "$sra_id" #remove the folder that is created, keep only the fastq file

    fname=${sra_id}
    R1="_1.fastq"
    R2="_2.fastq"

    # Create output folder for the sample
    mkdir -p "$dir/KNEADDATA/$fname"
    cd "$dir/KNEADDATA/$fname"
    kneaddata -i "$dir/data/$fname$R1" -i "$dir/data/$fname$R2" \
    -db "$util/KD_HOMO/" --run-fastqc-start --run-fastqc-end \
    --output "$dir/KNEADDATA/$fname" \
    --bowtie2-options="--very-fast -p 12" -t 2 -p 6 \
    --trimmomatic "$util/Trimmomatic-0.39"

    cd $dir/KNEADDATA/$fname
    find . -iname *paired_1.fastq* -exec mv {} $dir/CLEAN/ \;
    find . -iname *paired_2.fastq* -exec mv {} $dir/CLEAN/ \;

    cd $dir/KNEADDATA/$fname
    find . -iname *paired_1.fastq* -exec mv {} $dir/CLEAN/ \;
    find . -iname *paired_2.fastq* -exec mv {} $dir/CLEAN/ \;

    cd $dir/CLEAN
    fname=${sra_id}
    p1="${fname}_1_kneaddata_paired_1.fastq"
    p2="${fname}_1_kneaddata_paired_2.fastq"
    p3="${fname}.fastq"

    mv $p1 $p2 $dir/KNEADDATA/unpaired/

    # Removed unwanted files 
    cd $dir/data
    rm *.fastq
    cd $dir/KNEADDATA
    rm -r $fname
    cd $dir/CLEAN
    rm *.fastq

    ###fix the reads######################
    cd $KDdir
    for g in $KDdir/*_paired_1.fastq
    do
        a=$(basename "${g}") # Extracts the filename that is a = "SRR_paired_1.fastq"
        b=${g%_paired_1.fastq}_paired_2.fastq # removes _paired_1.fastq and adds _paired_2.fastq
        c=$(echo ${a} | cut -d_ -f1) # Extracts the sample name from the filename using cut.It splits the filename at underscores (_) and takes the first field.
        echo repair.sh in1=${g} in2=${b} out=$MAGSdir/sortedFQ/${c}_sorted_1.fastq out2=$MAGSdir/sortedFQ/${c}_sorted_2.fastq 
    done > fastq_sorted.sh
    sh fastq_sorted.sh


    cd $MAGSdir/metawrap/
    SFQ=$MAGSdir/sortedFQ
    ######## Binning###########
    f="${MAGSdir}/metawrap/trimmed_scaffolds/${sra_id}_scaffolds_trimmed.fasta"
    # a=$(basename "${f}")
    # b=$(dirname "${f}")
    # c=$(echo ${a} | cut -d_ -f1)
    # metawrap binning -t 13 --metabat2 --maxbin2 -a ${f} -o $MAGSdir/metawrap/bins/${c}_bins \
    # $SFQ/${c}_sorted_1.fastq $SFQ/${c}_sorted_2.fastq > ${sra_id}_metawrap_binning_log.txt

    # ### ########Bin refinement##################
    module load checkm/1.0.18
    cd $MAGSdir/metawrap/
    # g="${MAGSdir}/metawrap/bins/${sra_id}_bins"
    # a=$(basename "${g}")
    # metawrap bin_refinement -c 50 -x 10 -t 13 -o $MAGSdir/metawrap/refinement/${a}_refined -A ${g}/metabat2_bins -B ${g}/maxbin2_bins > ${sra_id}_metawrap_refinement_log.txt


    ## ########Bin quantification##############
    conda deactivate
    conda activate $util2/metawrap_env

    e="$MAGSdir/metawrap/trimmed_scaffolds"
    g="$MAGSdir/sortedFQ"
    f="${MAGSdir}/metawrap/refinement/${sra_id}_bins_refined/metawrap_50_10_bins"
    b=$(dirname "${f}")
    c=$(echo ${b} | cut -d / -f13)
    d=$(echo ${c} | cut -d_ -f1)
    metawrap quant_bins -t 20 -b ${f} -o ${b}/quant_bins -a ${e}/${d}_scaffolds_trimmed.fasta ${g}/${d}_sorted_1.fastq ${g}/${d}_sorted_2.fastq > ${sra_id}_metawrap_quan_log.txt

    conda deactivate
  conda deactivate
    
    cd $KDdir
    rm *.fastq
    cd $MAGSdir/sortedFQ
    rm *.fastq
done


# ###change names of bins
# for f in $MAGSdir/metawrap/refinement/*refined/metawrap_50_10_bins/*fa
# do
#   a=$(basename "${f}")
#   b=$(dirname "${f}")
#   c=$(echo ${b} | cut -d / -f13)
#   d=$(echo ${c} | cut -d_ -f1)
#   mv "${b}/${a}" "${b}/${d}_${a}"
# done


# ###cp MAGs to folder for GTDB
# cd $dir
# mkdir GTDB
# cd GTDB
# mkdir MAGS
# cd $MAGSdir/metawrap/refinement/
# find . -iname *.fa -exec cp {} $dir/GTDB/MAGS/ \;


##TO DO
# ####change name of stats and quant_bins
# for f in $MAGSdir/metawrap/refinement/*_refined/metawrap_50_10_bins.stats
# do
#   a=$(dirname "${f}")
#   b=$(basename "${f}")
#   c=$(echo ${a} | cut -d / -f13)
#   d=$(echo ${c} | cut -d_ -f1)
#   mv "${a}/${b}" "${a}/${d}_bins.stats"
#   mv "${a}/quant_bins/bin_abundance_table.tab" "${a}/quant_bins/${d}_bin_abundance_table.tab"
# done


####change name of stats and quant_bins
for f in $MAGSdir/metawrap/refinement/*_refined/metawrap_50_10_bins.stats
do
  a=$(dirname "${f}")
  b=$(basename "${f}")
  c=$(echo ${a} | cut -d / -f13)
  d=$(echo ${c} | cut -d_ -f1)
  mv "${a}/${b}" "${a}/${d}_bins.stats"
  mv "${a}/quant_bins/bin_abundance_table.tab" "${a}/quant_bins/${d}_bin_abundance_table.tab"
done
