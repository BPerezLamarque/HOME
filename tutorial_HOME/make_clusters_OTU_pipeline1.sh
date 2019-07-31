
####    PIPELINE 1 - HOME

####    CLUSTERING SEQUENCES  AND  MAKE ALIGNMENTS


path="YOUR_WORKING_DIRECTORY"
path_usearch=/Users/applications/uparse/usearch  #PATH USEARCH
path_pipelines=/Users/appli/pipelines/scripts/ #PATH BMP PIPELINES  (downloaded on http://www.brmicrobiome.org/16sillumina)
path_db=/Users/database/ #PATH GREENGENES AND RDP DATABASES

source /macqiime/configs/bash_profile.txt  #LOAD QIIME (see http://qiime.org)


#### GREAT APES DATASET (https://doi.org/10.5061/dryad.023s6/3 (Sanders, et al., 2014)).


seqs_file=$path/ochman_original_data/seqs.fna #FASTA READS (FASTA file of all reads with sample name specified in the header; see example https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/seqs.fna)

map_file=$path/metadata/samples.txt #METADATA FILE (list of sample; see example https://github.com/BPerezLamarque/HOME/blob/master/tutorial_HOME/samples.txt)


while read sample; do
    grep -A1 $sample $seqs_file  > $path/ochman_original_data/seqs_$sample.fas
done <$map_file




###############################################################################################
###############################   STEP1:  OTU CLUSTERING    (E.G. 97%)   ######################
###############################################################################################



###########################
#  OTU CLUSTERING
###########################

thre=97  #THRESHOLD VALUE
mkdir $path/analysis/
mkdir $path/analysis/97/
cd $path/analysis/97/

echo "Dereplication"
$path_usearch -derep_fulllength $seqs_file -fastaout reads_16s_derep.fa -sizeout

echo "Abundance sort and discard singletons"
$path_usearch -sortbysize reads_16s_derep.fa -fastaout reads_16s_sorted.fa -minsize 2

echo "OTU Clustering de novo:"
$path_usearch -cluster_otus reads_16s_sorted.fa -otu_radius_pct 3 -otus cluster_reads_OTU_97.fa # RADIUS  OF OTU = 100 - THRESHOLD

echo "Chimera filtering:"
$path_usearch -uchime_ref cluster_reads_OTU_97.fa -db $path_db/rdp_gold.fa -strand plus -nonchimeras cluster_reads_OTU_97_clean.fa

echo "Fasta Formatter:"
/Users/appli/fastx/fasta_formatter -i cluster_reads_OTU_97_clean.fa  -o cluster_reads_OTU_97_format.fa

echo "Renamer:"
perl $path_pipelines/bmp-otuName.pl -i cluster_reads_OTU_97_format.fa -o cluster_OTU_97.fa

echo "Taxonomy:"
assign_taxonomy.py -i cluster_OTU_97.fa -o $path/analysis/97/ -r $path_db/97_otus.fasta -t $path_db/97_otu_taxonomy.txt

echo "Map reads back to OTU database:"
$path_usearch -usearch_global $seqs_file -db cluster_OTU_97.fa -strand plus -id 0.97 -uc cluster_OTU_97_reads.uc
cut -f 9- cluster_OTU_97_reads.uc > cluster_OTU_97_reads.txt

echo "Convert OTU Table:"
python $path_pipelines/bmp-map2qiime.py cluster_OTU_97_reads.uc > OTU_table_97_temp.txt
make_otu_table.py -i OTU_table_97_temp.txt -t cluster_OTU_97_tax_assignments.txt -o OTU_table_97.biom
biom summarize-table -i  OTU_table_97.biom -o OTU_table_summary.txt
biom convert -i  OTU_table_97.biom  -o OTU_table_97.txt --to-tsv



###############################################################################################
###############################   STEP2:  MAKE OTU ALIGNMENTS         #########################
###############################################################################################


#  SELECT THE 75% CORE OTUs

echo "Compute core OTUs:"
compute_core_microbiome.py -i OTU_table_97.biom -o cores/

core=75
mkdir $path/alignments/
mkdir $path/alignments/$thre
cp cores/core_otus_75.txt $path/alignments/$thre

grep -o "OTU[0-9][0-9]*" cores/core_otus_75.txt > $path/alignments/$thre/list_core_OTU_97_75.txt #LIST OF CORES OTUs

# MAKE ALIGNMENTS

thre=97
cd $path/alignments/$thre

mkdir $path/alignments/alignments_$thre

while read OTU; do
    echo $OTU
    grep -E "$OTU$" $path/analysis/97/cluster_OTU_97_reads.txt | cut -f -1 >  $path/alignments/$thre/"name_seq_"$thre"_OTU_"$OTU".txt"
    grep -o "^[A-Z]*[0-9]*[A-Z]*" $path/alignments/$thre/"name_seq_"$thre"_OTU_"$OTU".txt" | sort -u | sed '1d' > $path/alignments/$thre/samples_$OTU.txt
    touch $path/alignments/$thre/alignment_97_$OTU.fas
    rm $path/alignments/$thre/alignment_97_$OTU.fas
    touch $path/alignments/$thre/alignment_97_$OTU.fas
        while read sample; do
        echo $sample
        grep $sample $path/alignments/$thre/"name_seq_"$thre"_OTU_"$OTU".txt" > $path/alignments/$thre/list_reads_OTU_samples.txt
        touch $path/alignments/$thre/sequences_by_sample_by_OTU.fas
        rm $path/alignments/$thre/sequences_by_sample_by_OTU.fas
        touch $path/alignments/$thre/sequences_by_sample_by_OTU.fas
        name=$(cat $path/alignments/$thre/list_reads_OTU_samples.txt | paste -sd "|" -)
        awk "/$name/{print;getline;print}" $path/ochman_original_data/seqs_$sample.fas >> $path/alignments/$thre/sequences_by_sample_by_OTU.fas
        $path_usearch -derep_fulllength $path/alignments/$thre/sequences_by_sample_by_OTU.fas -fastaout $path/alignments/$thre/sequences_by_"$sample"_by_"$OTU"_derep_1.fas -sizeout -quiet -relabel $sample"_"
        $path_usearch -sortbysize $path/alignments/$thre/sequences_by_"$sample"_by_"$OTU"_derep_1.fas -fastaout $path/alignments/$thre/sequences_by_"$sample"_by_"$OTU"_derep.fas -minsize 2 -quiet
        /Users/appli/fastx/fasta_formatter -i $path/alignments/$thre/sequences_by_"$sample"_by_"$OTU"_derep.fas -o $path/alignments/$thre/sequences_by_"$sample"_by_"$OTU".fas
        rm $path/alignments/$thre/sequences_by_"$sample"_by_"$OTU"_derep*
        ###########   MAKE THE MAIN ALIGNMENT (regroup sequences)
        head -n 2 $path/alignments/$thre/sequences_by_"$sample"_by_"$OTU".fas | sed "1s/.*/>$sample/" >> $path/alignments/$thre/alignment_97_$OTU.fas
        rm $path/alignments/$thre/sequences*
    done <$path/alignments/$thre/samples_$OTU.txt

    ###########   MAKE THE MAIN ALIGNMENT (align sequences)
    align_seqs.py -i $path/alignments/$thre/alignment_97_$OTU.fas -m muscle -o $path/alignments/$thre/output/
    filter_alignment.py -i "$path/alignments/$thre/output/alignment_97_"$OTU"_aligned.fasta" -g 0.75 --suppress_lane_mask_filter
    cp "$path/alignments/"$thre"/alignment_97_"$OTU"_aligned_pfiltered.fasta" $path"/alignments/alignments_"$thre"/alignment_97_"$OTU".fas"


    rm $path/alignments/$thre/alignment*
    rm $path/alignments/$thre/list_reads*
    rm $path/alignments/$thre/name_seq*
    rm $path/alignments/$thre/samples_$OTU*
    rm -r output

done <$path/alignments/$thre/list_core_OTU_97_75.txt



