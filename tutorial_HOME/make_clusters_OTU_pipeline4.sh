
####    PIPELINE 4 - PERFORM ASV CLUSTERING WITH VSEARCH (e.g. to assess phylosymbiosis - but not to run HOME !)

####    CLUSTERING SEQUENCES


## You first need to install PYTHON3, VSEARCH, and CUTADAPT
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
###############################   STEP1:  OTU CLUSTERING         ##############################
###############################################################################################

## to see all the details: https://github.com/frederic-mahe/swarm/wiki/Fred's-metabarcoding-pipeline

path="YOUR_WORKING_DIRECTORY"

cd $path

# Prepare your working directory

# add the FASTA file "reads_16s.fa" (single file containing all the metabarcoding reads) in your working directory.

# add also FASTA files "reads_Sk_16s.fa" (for every sample "Sk", you must have a separate file containing all the metabarcoding reads from this sample).


##  Step 1-A:  Dereplicate of fasta file of every  individual sample  using VSEARCH

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


##  Step 1-B:  Dereplication of fasta file of all the reads using VSEARCH

vsearch \
    --derep_fulllength reads_16s.fa \
    --sizein \
    --sizeout \
    --relabel_sha1 \
    --fasta_width 0 \
    --output reads_16s_derep.fa
# sha1 (encoding system) is giving the same names to the identical amplicons across samples


# Sort by size
vsearch -sortbysize reads_16s_derep.fa -output reads_16s_sorted.fa -minsize 1


##  Step 1-C:  Denoising into amplicons single variants (ASV)

vsearch --cluster_unoise  reads_16S_sorted.fa --minsize 8 --centroids reads_16S_ASV.fa --uc clusters_16S_ASV.uc  --sizein --sizeout


# Make file that mapped reads to OTUs
python3 map2qiime.py clusters_16S_ASV.uc > reads_16S_mapped_ASV.txt # python script map2qiime.py available in https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/map2qiime.py


# One line per OTU sequence
vsearch --fasta_width 0 \
    --sortbysize reads_16S_ASV.fa \
    --output reads_16S_ASV_final.fa

# Make stats file
python3 make_stats.py reads_16S_ASV_final.fa > stats_16S_ASV.txt # python script make_stats.py available in https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/make_stats.py


##  Step 1-D:  Chimera filtering (use uchime3_denovo for denoised OTUs)

vsearch --uchime3_denovo reads_16S_ASV.fa --uchimeout reads_16S_ASV.uchime


##  Step 1-E:  Assign taxonomy

vsearch --usearch_global reads_16S_ASV_final.fa \
    --threads 2 \
    --dbmask none \
    --qmask none \
    --rowlen 0 \
    --notrunclabels \
    --userfields query+id1+target \
    --maxaccepts 0 \
    --maxrejects 32 \
    --top_hits_only \
    --output_no_hits \
    --db SILVA_138_SSURef_NR99_tax_silva_515F_806R.fasta \
    --id 0.5 \
    --iddef 1 \
    --userout "taxonomy_16S_ASV.txt"


##  Step 1-F:  Make OTU table

STATS="stats_16S_ASV.txt"
MAPPED_OTUS="reads_16S_mapped_ASV.txt"
REPRESENTATIVES="reads_16S_ASV_final.fa"
UCHIME="reads_16S_ASV.uchime"
ASSIGNMENTS="taxonomy_16S_ASV.txt"
OTU_TABLE="OTU_table_16S_ASV.txt"

SCRIPT="OTU_contingency_table.py" # this script is located https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/OTU_contingency_table.py


python3 \
    "${SCRIPT}" \
    "${REPRESENTATIVES}" \
    "${STATS}" \
    "${MAPPED_OTUS}" \
    "${UCHIME}" \
    "${ASSIGNMENTS}" \
    reads_S[0-9]*_16s_derep.fa > "${OTU_TABLE}"

#NB: check that all the individual fasta files from each sample is present by typing "ls reads_S[0-9]*_16s.fa"


# Filter per OTU size or spread, quality and chimeric status:

# Filter per OTU size or spread, quality and chimeric status:
TABLE="OTU_table_16S_ASV.txt"
FILTERED="${TABLE/.txt/_filtered.txt}"

head -n 1 "${TABLE}" > "${FILTERED}"
cat "${TABLE}" | awk '$5 == "N" && $4 >= 150 && $2 >= 5 && $6 >= 2' >> "${FILTERED}"
# remove chimera, remove OTU shorter than 150 bp, only keep OTUs represented by at least 5 sequences, and spread in at least 2 samples



# RUN THE R SCRIPT TO CLEAN THE OTU TABLE (REMOVE CONTAMINANTS) AND SELECT THE CORE OTUs
# the script can be found here: https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/select_core_OTUs.R


