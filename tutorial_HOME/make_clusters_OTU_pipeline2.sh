
####    PIPELINE 2 - HOME

####    CLUSTERING SEQUENCES  AND  MAKE ALIGNMENTS


## You first need to install PYTHON3, VSEARCH, CUTADAPT and SWARM
## see https://github.com/frederic-mahe/swarm/wiki/Fred's-metabarcoding-pipeline for details


#####################################################################################################
######################  STEP0: PREPRE THE DATABASE FOR THE TAXONOMIC ASSIGNATION USING SILVA ########
#####################################################################################################

## to see all the details:  https://github.com/frederic-mahe/stampa

path="YOUR_WORKING_DIRECTORY"

cd $path

RELEASE=138
URL="https://ftp.arb-silva.de/release_${RELEASE}/Exports"
INPUT="SILVA_${RELEASE}_SSURef_NR99_tax_silva.fasta.gz"

# Download and check
wget -c ${URL}/${INPUT}{,.md5} && md5 -r ${INPUT}.md5

# Define variables and output files
OUTPUT="${INPUT/.fasta.gz/_515F_926R.fasta}"
LOG="${INPUT/.fasta.gz/_515F_926R.log}"


# Precise the sequences of the primers you used for the metabarcoding (here 16S V1-3 regions)
PRIMER_F="AGAGTTTGATCCTGGCTCAG"
PRIMER_R="TGCTGCCTCCCGTAGGAGT"
ANTI_PRIMER_R="ACTCCTACGGGAGGCAGCA"
MIN_LENGTH=32
MIN_F=$(( ${#PRIMER_F} * 2 / 3 ))
MIN_R=$(( ${#PRIMER_R} * 2 / 3 ))
CUTADAPT="cutadapt --discard-untrimmed --minimum-length ${MIN_LENGTH}"

# Trim forward & reverse primers, format
zcat < "${INPUT}" | sed '/^>/ ! s/U/T/g' | \
${CUTADAPT} -g "${PRIMER_F}" -O "${MIN_F}" - 2> "${LOG}" | \
${CUTADAPT} -a "${ANTI_PRIMER_R}" -O "${MIN_R}" - 2>> "${LOG}" | \
sed '/^>/ s/;/|/g ; /^>/ s/ /_/g ; /^>/ s/_/ /1' > "${OUTPUT}"



###############################################################################################
###############################   STEP1:  SWARM CLUSTERING         ############################
###############################################################################################

## to see all the details: https://github.com/frederic-mahe/swarm/wiki/Fred's-metabarcoding-pipeline

path="YOUR_WORKING_DIRECTORY"

cd $path

# Prepare your working directory

# add the FASTA file "reads_16s.fa" (single file containing all the metabarcoding reads) in your working directory.

# add also FASTA files "reads_Sk_16s.fa" (for every sample "Sk", you must have a separate file containing all the metabarcoding reads from this sample).



swarm=/Users/appli/swarm/swarm-3.0.0-macos-x86_64/bin/swarm #PATH SWARM v3


# Step 1-A: Dereplicate of fasta file of every  individual sample  using VSEARCH

declare -a list_samples=(S1 S2 S3 S4) # declare the list of the samples names "Sk"

for k in "${list_samples[@]}"
do
    vsearch \
        --derep_fulllength reads_S"$k"_16s.fa \
        --sizein \
        --sizeout \
        --relabel_sha1 \
        --fasta_width 0 \
        --output reads_S"$k"_16s_derep.fa
done
# sha1 (encoding system) is giving the same names to the identical amplicons across samples


##  Step 1-B: Dereplication of fasta file of all the reads using VSEARCH

vsearch \
    --derep_fulllength reads_16s.fa \
    --sizein \
    --sizeout \
    --relabel_sha1 \
    --fasta_width 0 \
    --output reads_16s_derep.fa
# sha1 (encoding system) is giving the same names to the identical amplicons across samples



##   Step 1-C: SWARM clustering

FINAL_FASTA="swarm_otu.fas"
THREADS=2  #number of cores

$swarm \
    -d 1 -f -t ${THREADS} -z \
    -i ${FINAL_FASTA/.fas/_1f.struct} \
    -s "stats_swarm_1f.txt" \
    -w "reads_16S_swarm_1f.fa" \
    -o "reads_16S_mapped_swarm_1f.txt" < reads_16S_derep.fa


# Sort representatives
vsearch --fasta_width 0 \
    --sortbysize reads_16S_swarm_1f.fa \
    --output reads_16S_swarm_1f_final.fa


##  Step 1-D: Chimera checking

vsearch --uchime_denovo reads_16S_swarm_1f_final.fa --uchimeout reads_16S_swarm_1f.uchime  --nonchimeras reads_16S_swarm_1f_nonchimeras.fa


##  Step 1-E:  Assign taxonomy to OTUs using known reference sequences (obtained from SILVA)

vsearch --usearch_global "reads_16S_swarm_1f_nonchimeras.fa" \
    --threads 4 \
    --dbmask none \
    --qmask none \
    --rowlen 0 \
    --notrunclabels \
    --userfields query+id1+target \
    --maxaccepts 0 \
    --maxrejects 32 \
    --top_hits_only \
    --output_no_hits \
    --db SILVA_138_SSURef_NR99_tax_silva_515F_926R.fasta \
    --id 0.5 \
    --iddef 1 \
    --userout "taxonomy_otu_16S_swarm_1f.txt"



##  Step 1-F:  Make the OTU table

STATS="stats_swarm_1f.txt"
SWARMS="reads_16S_mapped_swarm_1f.txt"
REPRESENTATIVES="reads_16S_swarm_1f_final.fa"
UCHIME="reads_16S_swarm_1f.uchime"
ASSIGNMENTS="taxonomy_otu_16S_swarm_1f.txt"
OTU_TABLE="OTU_table_16S_swarm_1f.txt"


SCRIPT="OTU_contingency_table.py"  # this script is located https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/OTU_contingency_table.py

python3 \
    "${SCRIPT}" \
    "${REPRESENTATIVES}" \
    "${STATS}" \
    "${SWARMS}" \
    "${UCHIME}" \
    "${ASSIGNMENTS}" \
    reads_S[0-9]*_16s_derep.fa > "${OTU_TABLE}"


# Filter per OTU size or spread, quality and chimeric status:

TABLE="OTU_table_16S_swarm_1f.txt"
FILTERED="${TABLE/.txt/_filtered.txt}"

head -n 1 "${TABLE}" > "${FILTERED}"
cat "${TABLE}" | awk '$5 == "N" && $4 >= 150 && $2 >= 5 && $6 >= 2' >> "${FILTERED}"
# remove chimera, remove OTU shorter than 150 bp, only keep OTUs represented by at least 5 sequences, and spread in at least 2 samples



# RUN THE R SCRIPT TO CLEAN THE OTU TABLE (REMOVE CONTAMINANTS) AND SELECT THE CORE OTUs
# the script can be found here: https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/select_core_OTUs.R




###############################################################################################
###############################   STEP2:  MAKE OTU ALIGNMENTS         #########################
###############################################################################################


# This step rely on bash script, VSEARCH, and python script (available in https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/fasta_extract.py)
# it also relies on MAFFT to align the sequences (but any other software can be used)

# it requires to first run the R script select_core_OTUs. R (available in https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/select_core_OTUs.R) to generate the files  list_core_OTU.txt and list_sample_OTU_$OTU.txt


MAPPED_OTUS="reads_16S_mapped_swarm_1f.txt"

while read OTU; do
    echo $OTU
    touch alignment_$OTU.fas
    rm alignment_$OTU.fas
    touch alignment_$OTU.fas
    while read sample; do
        echo $sample

        list_reads=$(grep "$OTU" $MAPPED_OTUS )

        # Set space as the delimiter
        IFS=' '

        #Read the split words into an array based on space delimiter
        read -a strarr <<< "$list_reads"

        # Print each value of the array by using the loop
        touch "list_sequences.txt"
        rm "list_sequences.txt"
        touch "list_sequences.txt"
        for val in "${strarr[@]}";
        do
            sequence=$(echo $val | sed 's/;size=*.*[0-9]//')
            echo $sequence >> "list_sequences.txt"
        done
        
        python3 fasta_extract.py -f reads_S"$sample"_16s_derep.fa -k "list_sequences.txt" > "sequences_by_"$sample"_by_"$OTU".fas"
        # This script  fasta_extract.py is available here: https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/fasta_extract.py
        
        ###########   get the most abundant sequence per sample
        vsearch -sortbysize "sequences_by_"$sample"_by_"$OTU".fas" -output "sequences_by_"$sample"_by_"$OTU"_sorted.fas" -minsize 1 -quiet
        
        vsearch --fasta_width 0 --quiet \
            --sortbysize "sequences_by_"$sample"_by_"$OTU"_sorted.fas" \
            --output "sequences_by_"$sample"_by_"$OTU".fas"
        
        sed "s/>/>"$sample"_/" "sequences_by_"$sample"_by_"$OTU".fas" > "sequences_by_"$sample"_by_"$OTU"_label.fas"
        
        head -n 2 "sequences_by_"$sample"_by_"$OTU"_label.fas" >> alignment_$OTU.fas
        
        rm "sequences_by_"$sample"_by_"*
        
    done <list_sample_OTU_$OTU.txt


    ###########   Align the sequences (use MAFFT or any other software)
    
    mafft --quiet alignment_$OTU.fas > alignment_mafft_$OTU.fas

    rm alignment_$OTU.fas

done <list_core_OTU.txt


