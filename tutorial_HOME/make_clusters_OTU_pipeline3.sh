
####    PIPELINE 3 - HOME

####    CLUSTERING SEQUENCES  AND  MAKE ALIGNMENTS


## You first need to install PYTHON, VSEARCH, and CUTADAPT
## see https://github.com/frederic-mahe/swarm/wiki/Fred's-metabarcoding-pipeline for details


#####################################################################################################
######################  STEP0: PREPRE THE DATABASE FOR THE TAXONOMIC ASSIGNATION USING SILVA ########
#####################################################################################################

## to see all the details:  https://github.com/frederic-mahe/stampa

path="YOUR_WORKING_DIRECTORY"

cd $path

RELEASE=132
URL="https://www.arb-silva.de/fileadmin/silva_databases/release_${RELEASE}/Exports"
INPUT="SILVA_${RELEASE}_SSURef_Nr99_tax_silva.fasta.gz"

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
zcat "${INPUT}" | sed '/^>/ ! s/U/T/g' | \
${CUTADAPT} -g "${PRIMER_F}" -O "${MIN_F}" - 2> "${LOG}" | \
${CUTADAPT} -a "${ANTI_PRIMER_R}" -O "${MIN_R}" - 2>> "${LOG}" | \
sed '/^>/ s/;/|/g ; /^>/ s/ /_/g ; /^>/ s/_/ /1' > "${OUTPUT}"

gunzip -c  "${INPUT}" | sed '/^>/ ! s/U/T/g' | \
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



## Dereplication using VSEARCH

vsearch \
    --derep_fulllength reads_16s.fa \
    --sizein \
    --sizeout \
    --relabel_sha1 \
    --fasta_width 0 \
    --output reads_16s_derep.fa

# Sort by size

vsearch -sortbysize reads_16s_derep.fa -output reads_16s_sorted.fa -minsize 1

## OTU clustering (at 97%)

vsearch -cluster_fast  reads_16s_sorted.fa --id 0.97 --centroids reads_16s_OTU_97.fa --uc clusters_97.uc


# Relabel
sed 's/>/>OTU/g' reads_16s_OTU_97.fa > reads_16s_OTU_97_2.fa
cat reads_16s_OTU_97_2.fa > reads_16s_OTU_97.fa
rm reads_16s_OTU_97_2.fa

sed 's/;size=[0-9]*//g' reads_16s_OTU_97.fa > reads_16s_OTU_97_clean.fa
vsearch --fasta_width 0 \
--sortbysize reads_16s_OTU_97_clean.fa \
--output reads_16s_OTU_97_final.fa

# Assign taxonomy
vsearch --usearch_global reads_16s_OTU_97_final.fa \
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
--db ../../alignments/taxonomy/SILVA_132_SSURef_Nr99_tax_silva_515F_926R.fasta \
--id 0.5 \
--iddef 1 \
--userout "taxonomy_16s_OTU_97.txt"

# Chimera filtering

vsearch --uchime_denovo reads_16s_OTU_97.fa --uchimeout reads_16s_OTU_97.uchime
sed 's/;size=[0-9]*//g' reads_16s_OTU_97.uchime > reads_16s_OTU_97_clean.uchime

# Map the raw reads back to OTU (with a 97% threshold)
vsearch -usearch_global reads_16s.fa -db reads_16s_OTU_97_final.fa --threads 4 -strand plus -id 0.97 -uc reads_16s_mapped_otu.uc
cut -f 9- reads_16s_mapped_otu.uc > reads_16s_mapped_otu.txt

python bmp-map2qiime.py reads_16s_mapped_otu.uc > reads_16s_mapped_otu_temp.txt # python script present in https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/bmp-map2qiime.py ; from the BMP projects https://www.brmicrobiome.org)

grep ">" reads_16s_OTU_97_clean.fa > reads_16s_mapped_otu_clean2.txt
sed 's/>//g' reads_16s_mapped_otu_clean2.txt > reads_16s_mapped_otu_clean3.txt
sed '1,$s/$/    20/' reads_16s_mapped_otu_clean3.txt > reads_16s_mapped_otu_clean2.txt


# Make OTU table
STATS="reads_16s_mapped_otu_clean2.txt"
SWARMS="reads_16s_mapped_otu_temp.txt"
REPRESENTATIVES="reads_16s_OTU_97_final.fa"
UCHIME="reads_16s_OTU_97_clean.uchime"
ASSIGNMENTS="taxonomy_16s_OTU_97.txt"
OTU_TABLE="OTU_table_16S_97.txt"


SCRIPT="../../alignments/script/OTU_contingency_table_3.py" # this script is located https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/OTU_contingency_table_3.py


python3 \
    "${SCRIPT}" \
    "${REPRESENTATIVES}" \
    "${STATS}" \
    "${SWARMS}" \
    "${UCHIME}" \
    "${ASSIGNMENTS}" \
    reads_S[0-9]*_16s.fa > "${OTU_TABLE}"

# Filter per OTU size or spread, quality and chimeric status:

TABLE="OTU_table_16S_97.txt"
FILTERED="${TABLE/.txt/_filtered.txt}"

head -n 1 "${TABLE}" > "${FILTERED}"
cat "${TABLE}" | awk '$7 == "N" && $5 >= 250 && $8 >= 2' >> "${FILTERED}"


#  Dereplicate the fasta file of every sample

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



# Map the raw reads back to OTU

rm list_sample.txt
for k in "${list_samples[@]}"
do
echo $k >> list_sample.txt
done

while read sample; do
echo $sample
vsearch -usearch_global reads_S"$sample"_16s.fa -db reads_16s_OTU_97_final.fa --threads 4 -strand plus -id 0.97 -uc reads_S"$sample"_16s_mapped_otu.uc
cut -f 9- reads_S"$sample"_16s_mapped_otu.uc > reads_S"$sample"_16s_mapped_otu.txt
done <list_sample.txt



# RUN THE R SCRIPT TO CLEAN THE OTU TABLE (REMOVE CONTAMINANTS) AND SELECT THE CORE OTUs
# the script can be found here: https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/select_core_OTUs.R




###############################################################################################
###############################   STEP2:  MAKE OTU ALIGNMENTS         #########################
###############################################################################################


# This step rely on bash script, VSEARCH, python script (https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/fasta_extract.py) and the FASTX-Toolkit (http://hannonlab.cshl.edu/fastx_toolkit/commandline.html)
# it also relies on QIIME1 to align the sequences (but any other software can be used)


while read OTU; do
    echo $OTU
    touch alignment_$OTU.fas
    rm alignment_$OTU.fas
    touch alignment_$OTU.fas
    while read sample; do
        echo $sample
        grep -E "$OTU" reads_S"$sample"_16s_mapped_otu.txt | cut -f -1 > "name_seq_OTU_"$OTU".txt"
        python fasta_extract.py -f reads_S"$sample"_16s.fa -k "name_seq_OTU_"$OTU".txt" > sequences_by_sample_by_OTU.fas

        ###########   get the most abundant sequence
        vsearch -derep_fulllength sequences_by_sample_by_OTU.fas -output sequences_by_"$sample"_by_"$OTU"_derep_1.fas -sizeout -quiet -relabel $sample"_"
        vsearch -sortbysize sequences_by_"$sample"_by_"$OTU"_derep_1.fas -output sequences_by_"$sample"_by_"$OTU"_derep.fas -minsize 2 -quiet
        /Users/appli/fastx/fasta_formatter -i sequences_by_"$sample"_by_"$OTU"_derep.fas -o sequences_by_"$sample"_by_"$OTU".fas  # FASTX-Toolkit
        rm sequences_by_"$sample"_by_"$OTU"_derep*
        head -n 2 sequences_by_"$sample"_by_"$OTU".fas | sed "1s/.*/>$sample/" >> alignment_$OTU.fas
        rm sequences*
    done <list_sample_OTU_$OTU.txt

    ###########   Align the sequences (use QIIME1 or any other software)
    source activate qiime1
    align_seqs.py -i alignment_$OTU.fas -m muscle -o output/
    source deactivate
    cp output/alignment_"$OTU"_aligned.fasta alignment_"$name"_"$OTU".fas


    rm alignment*
    rm name_seq*
    rm -r output
done <list_core_OTU.txt








