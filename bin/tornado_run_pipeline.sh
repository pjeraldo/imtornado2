#!/usr/bin/env bash
#author = "Patricio Jeraldo"
#copyright = "Copyright 2014, Mayo Foundation for Medical Education and Research "
#credits = "Patricio Jeraldo,  Nicholas Chia"
#license = "GPL"
#version = "2.0.3"
#maintainer = "Patricio Jeraldo"
#email = "jeraldo.patricio@mayo.edu"
#status = "Production"

#helper functions
set -e

#die function
warn () {
    echo "$0:" "$@" >&2
}
die () {
    rc=$1
    shift
    warn "$@"
    exit $rc
}

isitthere () {
  [ -s $1 ] || die 1 "File $1 not found or zero sized."
}

#get the params for this run
source tornado-params.sh

#make results folder
[ -d $RESULTS ] || mkdir $RESULTS

#make results folder
[ -d $WORKSPACE ] || mkdir $WORKSPACE

isitthere $MAPPING

#get the first column of the mapping file, to get the IDs of the files involved
cut -f 1 $MAPPING | sed '/^#/d' > $WORKSPACE/ids.txt
#SANITY CHECKS
#count the number of samples
NSAMPLES=`cat $WORKSPACE/ids.txt| wc -l`
#count the number of fastq files in the directory. Assume paired end. COunt only R1.
NREAD1=`ls *R1*fastq|wc -l`
#are they equal?
[ $NSAMPLES == $NREAD1 ] || die 1 "Number of samples in mapping does not match number of fastq files."

#test the gnu parallel naming expression...
#if the sample names do not have a reasonable format, this will fail
NFOUND=0
for f in $(cat $WORKSPACE/ids.txt)
do
  N=$(ls ${f}${SPACER}*R1*fastq|wc -l)
  NFOUND=$[NFOUND+N]
done
#same as number of samples?
[ $NSAMPLES == $NFOUND ] || die 1 "Number of samples found with unambiguous ids do not mach number of fastq files found. Check the naming of your samples."


#Filter with Trimmomatic
#The TrimmomaticPE command

for name in `cat $WORKSPACE/ids.txt`
do
java -jar $TRIMMOMATIC SE -threads 4 -phred33 ${name}${SPACER}*R1*fastq $WORKSPACE/${name}_R1.cfastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$MINIMUM_LENGTH
java -jar $TRIMMOMATIC SE -threads 4 -phred33 ${name}${SPACER}*R2*fastq $WORKSPACE/${name}_R2.cfastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$MINIMUM_LENGTH
done

#We won't need to be somewhere else until the very end
cd $WORKSPACE
#Convert FASTQ to FASTA
echo "Converting FASTQ to FASTA"
for name in `cat ids.txt`
do
tornado_fastq2fasta.py ${name}_R1.cfastq ${name}_R1.fasta
tornado_fastq2fasta.py ${name}_R2.cfastq ${name}_R2.fasta
done

#Quality filtering (according to EMP protocol). Q3 trimming, 75% of length as cutoff, no ambigs, nothing on homopolymers
echo "Remove ambigs..."
#no check for polys... let OTUing tacke care of this
#there has to be a faster way to do this.
for name in `cat ids.txt`
do
tornado_remove_ambigs.sh ${name}_R1.fasta ${name}_R1.trim.fasta
tornado_remove_ambigs.sh ${name}_R2.fasta ${name}_R2.trim.fasta
done

#Rename read IDs accoring to sample ID (QIIME convention).
#Trim read lengths to specified numbers
#do not discard orphans
#output the common ids to PREFIX.common.accnos
[ -e ${PREFIX}.groups ] && rm ${PREFIX}.groups
[ -e ${PREFIX}_R1.fasta ] && rm ${PREFIX}_R1.fasta
[ -e ${PREFIX}_R2.fasta ] && rm ${PREFIX}_R2.fasta
[ -e ${PREFIX}.common.accnos ] && rm ${PREFIX}.common.accnos
echo "Rename and trim"
tornado_rename_ids.py --trim_r1 $R1_TRIM --trim_r2 $R2_TRIM $PREFIX *R1.trim.fasta
#remove empty lines from accnos file
sed -i'' '/^$/d' ${PREFIX}.common.accnos
#check newly created files
isitthere ${PREFIX}.groups
isitthere ${PREFIX}_R1.fasta
isitthere ${PREFIX}_R2.fasta
isitthere ${PREFIX}.common.accnos

#Concatenate paired reads. This is only used for R1+R2 OTUing.
#First, pick out the common ids
echo "Pick common ids..."

tornado_read_picker.py ${PREFIX}.common.accnos ${PREFIX}_R1.fasta ${PREFIX}_R1.common.fasta
tornado_read_picker.py ${PREFIX}.common.accnos ${PREFIX}_R2.fasta ${PREFIX}_R2.common.fasta
#flatten out the reads
echo "Flatten read files..."
tornado_flatten_fasta.py -i ${PREFIX}_R1.common.fasta -o ${PREFIX}_R1.common.flat.fasta
tornado_flatten_fasta.py -i ${PREFIX}_R2.common.fasta -o ${PREFIX}_R2.common.flat.fasta

#and concatenate
echo "Concatenate read pairs"
#concatenate padded

paste -d 'N' ${PREFIX}_R1.common.flat.fasta ${PREFIX}_R2.common.flat.fasta | awk 'NR%2 ==0; NR%2 ==1 {sub(/N.*/,""); print}' > ${PREFIX}.padded.fasta
awk 'NR%2 == 1; NR%2 == 0 {sub(/N/,""); print}' ${PREFIX}.padded.fasta > ${PREFIX}_paired.fasta

#calculate overall taxonomy if consensus is set
if [[ $CONSENSUS_TAXONOMY == "1" ]]
  then
  #unique samples
  mothur "#unique.seqs(fasta=${PREFIX}.padded.fasta)"
  mothur "#classify.seqs(fasta=${PREFIX}_R1.unique.fasta-${PREFIX}_R2.unique.fasta-${PREFIX}.padded.unique.fasta, name=${PREFIX}_R1.names-${PREFIX}_R2.names-${PREFIX}.padded.names, taxonomy=${DATA}/${TAXONOMY}.taxonomy, template=${DATA}/${TAXONOMY}.fna, probs=false, processors=${NPROC})"
fi

#OTU R1 and R2, separately
#dereplicate

#is VSEARCH is set
if [[ -n $VSEARCH ]]
then

  $VSEARCH -derep_fulllength ${PREFIX}_R1.fasta -output ${PREFIX}_R1.derep.fasta -sizeout
  $VSEARCH -derep_fulllength ${PREFIX}_R2.fasta -output ${PREFIX}_R2.derep.fasta -sizeout
  $VSEARCH -derep_fulllength ${PREFIX}_paired.fasta -output ${PREFIX}_paired.derep.fasta -sizeout

else
  #use mothur for this for now
  mothur "#unique.seqs(fasta=${PREFIX}_R1.fasta)"
  mothur "#unique.seqs(fasta=${PREFIX}_R2.fasta)"
  mothur "#unique.seqs(fasta=${PREFIX}_paired.fasta)"

  #get the "counts" .. returns PREFIX_R?.seq.count
  mothur "#count.seqs(name=${PREFIX}_R1.names)"
  mothur "#count.seqs(name=${PREFIX}_R2.names)"
  mothur "#count.seqs(name=${PREFIX}_paired.names)"

  #annotate the uniques with the sizes
  #and remove reads shorter than the specified trims
  echo "Annotating unique sizes..."
  echo "R1..."

  tornado_annotate_read_sizes.py $R1_TRIM ${PREFIX}_R1.unique.fasta ${PREFIX}_R1.count_table ${PREFIX}_R1.derep.fasta
  echo "R2..."

  tornado_annotate_read_sizes.py $R2_TRIM ${PREFIX}_R2.unique.fasta ${PREFIX}_R2.count_table ${PREFIX}_R2.derep.fasta
  echo "Paired..."

  tornado_annotate_read_sizes.py $[R1_TRIM + R2_TRIM] ${PREFIX}_paired.unique.fasta ${PREFIX}_paired.count_table ${PREFIX}_paired.derep.fasta

fi

#remove singletons
$USEARCH7 -sortbysize ${PREFIX}_R1.derep.fasta -output ${PREFIX}_R1.derep2.fasta -minsize 2
$USEARCH7 -sortbysize ${PREFIX}_R2.derep.fasta -output ${PREFIX}_R2.derep2.fasta -minsize 2
$USEARCH7 -sortbysize ${PREFIX}_paired.derep.fasta -output ${PREFIX}_paired.derep2.fasta -minsize 2

#OTU the reads!
$USEARCH7 -cluster_otus ${PREFIX}_R1.derep2.fasta -otus ${PREFIX}_R1.otus.fasta
$USEARCH7 -cluster_otus ${PREFIX}_R2.derep2.fasta -otus ${PREFIX}_R2.otus.fasta
$USEARCH7 -cluster_otus ${PREFIX}_paired.derep2.fasta -otus ${PREFIX}_paired.otus.fasta

echo "Rename otu ids..."
#rename the OTU ids...

tornado_otu_renamer.py ${PREFIX}_R1.otus.fasta ${PREFIX}_R1.otus2.fasta
tornado_otu_renamer.py ${PREFIX}_R2.otus.fasta ${PREFIX}_R2.otus2.fasta
tornado_otu_renamer.py ${PREFIX}_paired.otus.fasta ${PREFIX}_paired.otus2.fasta
#ok, this is weird... for now extract the name pairings (OTU and original ID)
#this is only for paired... we will use it for taxonomy and alignment
grep ">" ${PREFIX}_paired.otus2.fasta | tr -d '>'|  sed 's/;.*// ; s/ .*=/ /' > ${PREFIX}_paired.otu_id.map

#get only the original ids
cut -f 2 -d ' ' ${PREFIX}_paired.otu_id.map > ${PREFIX}_paired.original.accnos
#and pick it from the gapped read set

tornado_read_picker.py ${PREFIX}_paired.original.accnos ${PREFIX}.padded.fasta ${PREFIX}_paired.otus.gapped.fasta
#Rename it again? will it work? I'm assuming the order is kept.
tornado_otu_renamer.py ${PREFIX}_paired.otus.gapped.fasta ${PREFIX}_paired.tax.fasta
#it works!!

#get the taxonomy of the representatives
mothur "#classify.seqs(fasta=${PREFIX}_R1.otus2.fasta-${PREFIX}_R2.otus2.fasta-${PREFIX}_paired.tax.fasta, taxonomy=${DATA}/${TAXONOMY}.taxonomy, template=${DATA}/${TAXONOMY}.fna, processors=${NPROC}, iters=1000)"

#copy the taxonomy files to the results
cp ${PREFIX}_R1.otus2.${TAXONOMY}.wang.taxonomy ../$RESULTS/${PREFIX}_R1.probs.taxonomy
cp ${PREFIX}_R2.otus2.${TAXONOMY}.wang.taxonomy ../$RESULTS/${PREFIX}_R2.probs.taxonomy
cp ${PREFIX}_paired.tax.${TAXONOMY}.wang.taxonomy ../$RESULTS//${PREFIX}_paired.probs.taxonomy

#do this again, with no support values
mothur "#classify.seqs(fasta=${PREFIX}_R1.otus2.fasta-${PREFIX}_R2.otus2.fasta-${PREFIX}_paired.tax.fasta, taxonomy=${DATA}/${TAXONOMY}.taxonomy, template=${DATA}/${TAXONOMY}.fna, probs=false, processors=${NPROC})"

echo "Filter bad reads..."
#find the ones that are obviously not 16S
grep "unknown" ${PREFIX}_R1.otus2.${TAXONOMY}.wang.taxonomy| cut -f 1 > bad_taxonomy_R1.accnos
grep "unknown" ${PREFIX}_R2.otus2.${TAXONOMY}.wang.taxonomy| cut -f 1 > bad_taxonomy_R2.accnos
grep "unknown" ${PREFIX}_paired.tax.${TAXONOMY}.wang.taxonomy| cut -f 1 > bad_taxonomy_paired.accnos

#And remove those reads from the otu file
tornado_read_remover.py bad_taxonomy_R1.accnos ${PREFIX}_R1.otus2.fasta ${PREFIX}_R1.otus3.fasta
tornado_read_remover.py bad_taxonomy_R2.accnos ${PREFIX}_R2.otus2.fasta ${PREFIX}_R2.otus3.fasta
tornado_read_remover.py bad_taxonomy_paired.accnos ${PREFIX}_paired.otus2.fasta ${PREFIX}_paired.otus3.fasta
#In the paired set, prepare to align separately
grep ">" ${PREFIX}_paired.otus3.fasta | tr -d '>'| sed 's/;.*// ; s/ .*=/ /' | awk '{print $2"\t"$1}' > ${PREFIX}_paired.align.map
#get the ids to look for in the files...
cut -f 1 ${PREFIX}_paired.align.map > ${PREFIX}.to_align.accnos
#get those IDs from the original files
tornado_read_picker.py ${PREFIX}.to_align.accnos ${PREFIX}_R1.fasta ${PREFIX}_paired.R1.fasta
tornado_read_picker.py ${PREFIX}.to_align.accnos ${PREFIX}_R2.fasta ${PREFIX}_paired.R2.fasta
#rename the ids to the OTU ids
tornado_custom_read_renamer.py ${PREFIX}_paired.align.map ${PREFIX}_paired.R1.fasta ${PREFIX}_paired.R1.to_align.fasta
tornado_custom_read_renamer.py ${PREFIX}_paired.align.map ${PREFIX}_paired.R2.fasta ${PREFIX}_paired.R2.to_align.fasta
#align the reads... there may be still some bad reads in there
#probably this will need to be done using STK instead of fasta
echo "Align R1"
cmalign --cpu $NPROC -g --notrunc --sub --dnaout --noprob --sfile ${PREFIX}_R1.scores -o ${PREFIX}_R1.otus3.aligned.stk $DATA/$MODEL ${PREFIX}_R1.otus3.fasta
echo "Align R2"
cmalign --cpu $NPROC -g --notrunc --sub --dnaout --noprob --sfile ${PREFIX}_R2.scores -o ${PREFIX}_R2.otus3.aligned.stk $DATA/$MODEL ${PREFIX}_R2.otus3.fasta
echo "Align R1+R2: R1..."
cmalign --cpu $NPROC -g --notrunc --sub --dnaout --noprob --sfile ${PREFIX}_paired.R1.scores -o ${PREFIX}_paired.R1.aligned.stk $DATA/$MODEL ${PREFIX}_paired.R1.to_align.fasta
echo "Align R1+R2: R2..."
cmalign --cpu $NPROC -g --notrunc --sub --dnaout --noprob --sfile ${PREFIX}_paired.R2.scores -o ${PREFIX}_paired.R2.aligned.stk $DATA/$MODEL ${PREFIX}_paired.R2.to_align.fasta

#convert stk to fasta
echo "Convert to fasta"

tornado_stk2fasta.py ${PREFIX}_R1.otus3.aligned.stk ${PREFIX}_R1.otus3.aligned.fasta
tornado_stk2fasta.py ${PREFIX}_R2.otus3.aligned.stk ${PREFIX}_R2.otus3.aligned.fasta
tornado_stk2fasta.py ${PREFIX}_paired.R1.aligned.stk ${PREFIX}_paired.R1.aligned.fasta
tornado_stk2fasta.py ${PREFIX}_paired.R2.aligned.stk ${PREFIX}_paired.R2.aligned.fasta

echo "Filter bad reads"
#look for the badly aligned reads
sed '/#/d' ${PREFIX}_R1.scores | awk '{if ($7 < 0) print $2}' > bad_align_R1.accnos
sed '/#/d' ${PREFIX}_R2.scores | awk '{if ($7 < 0) print $2}' > bad_align_R2.accnos
sed '/#/d' ${PREFIX}_paired.R1.scores | awk '{if ($7 < 0) print $2}' > bad_align_paired.R1.accnos
sed '/#/d' ${PREFIX}_paired.R2.scores | awk '{if ($7 < 0) print $2}' > bad_align_paired.R2.accnos
cat bad_align_paired.R1.accnos bad_align_paired.R2.accnos |sort -u > bad_align_paired.accnos

#remove them from the OTU file 
tornado_read_remover.py bad_align_R1.accnos ${PREFIX}_R1.otus3.fasta ${PREFIX}_R1.otus.final.fasta
tornado_read_remover.py bad_align_R2.accnos ${PREFIX}_R2.otus3.fasta ${PREFIX}_R2.otus.final.fasta
tornado_read_remover.py bad_align_paired.accnos ${PREFIX}_paired.otus3.fasta ${PREFIX}_paired.otus.final.fasta

#copy the final OTUs into the results file
#because they are actually a results
cp ${PREFIX}_R1.otus.final.fasta ../$RESULTS
cp ${PREFIX}_R2.otus.final.fasta ../$RESULTS
cp ${PREFIX}_paired.otus.final.fasta ../$RESULTS

#and from the alignment
tornado_read_remover.py bad_align_R1.accnos ${PREFIX}_R1.otus3.aligned.fasta ${PREFIX}_R1.otus3.aligned.clean.fasta
tornado_read_remover.py bad_align_R2.accnos ${PREFIX}_R2.otus3.aligned.fasta ${PREFIX}_R2.otus3.aligned.clean.fasta
tornado_read_remover.py bad_align_paired.accnos ${PREFIX}_paired.R1.aligned.fasta ${PREFIX}_paired.R1.aligned.clean.fasta
tornado_read_remover.py bad_align_paired.accnos ${PREFIX}_paired.R2.aligned.fasta ${PREFIX}_paired.R2.aligned.clean.fasta
echo "Remove gapped columns"
#remove gapped columns from the alignment
tornado_remove_gaps.py ${PREFIX}_R1.otus3.aligned.clean.fasta ${PREFIX}_R1.otus3.aligned.clean2.fasta
tornado_remove_gaps.py ${PREFIX}_R2.otus3.aligned.clean.fasta ${PREFIX}_R2.otus3.aligned.clean2.fasta
tornado_remove_gaps.py ${PREFIX}_paired.R1.aligned.clean.fasta ${PREFIX}_paired.R1.aligned.clean2.fasta
tornado_remove_gaps.py ${PREFIX}_paired.R2.aligned.clean.fasta ${PREFIX}_paired.R2.aligned.clean2.fasta

#change periods to dashes in alignment, as well as change letters to uppercase
cat ${PREFIX}_R1.otus3.aligned.clean2.fasta| tr '.' '-' | tr 'a-z' 'A-Z' > ${PREFIX}_R1.aligned.fasta
cat ${PREFIX}_R2.otus3.aligned.clean2.fasta| tr '.' '-' | tr 'a-z' 'A-Z' > ${PREFIX}_R2.aligned.fasta
cat ${PREFIX}_paired.R1.aligned.clean2.fasta| tr '.' '-' | tr 'a-z' 'A-Z' > ${PREFIX}_paired.R1.aligned.fasta
cat ${PREFIX}_paired.R2.aligned.clean2.fasta| tr '.' '-' | tr 'a-z' 'A-Z' > ${PREFIX}_paired.R2.aligned.fasta

#copy the alignments too
cp ${PREFIX}_R1.aligned.fasta ../$RESULTS
cp ${PREFIX}_R2.aligned.fasta ../$RESULTS
cp ${PREFIX}_paired.R1.aligned.fasta ../$RESULTS
cp ${PREFIX}_paired.R2.aligned.fasta ../$RESULTS

#remove gaps from the paired alignments in the results, to create the independent OTUs
cat ${PREFIX}_paired.R1.aligned.fasta | tr -d '-' > ../$RESULTS/${PREFIX}_paired.R1.otus.final.fasta
cat ${PREFIX}_paired.R2.aligned.fasta | tr -d '-' > ../$RESULTS/${PREFIX}_paired.R2.otus.final.fasta

#Join the paired alignments
#find common reads... just in case
tornado_common_reads.py ${PREFIX}_paired.R1.aligned.fasta ${PREFIX}_paired.R2.aligned.fasta ${PREFIX}.aligned.common.accnos
tornado_read_picker.py ${PREFIX}.aligned.common.accnos ${PREFIX}_paired.R1.aligned.fasta ${PREFIX}_paired.R1.aligned.common.fasta
tornado_read_picker.py ${PREFIX}.aligned.common.accnos ${PREFIX}_paired.R2.aligned.fasta ${PREFIX}_paired.R2.aligned.common.fasta

#flatten reads
#cat ${PREFIX}_paired.R1.aligned.common.fasta | sed '/^>/s/$/xXx/ ; /^>/s/^/xXx/'| tr -d '\n'| sed 's/xXx/\n/g' | sed '/^$/d' > ${PREFIX}_paired.R1.aligned.flat.fasta
#cat ${PREFIX}_paired.R1.aligned.common.fasta | perl -pe 's/$/xXx/ if /^>/; s/^/xXx/ if /^>/' | tr -d '\n' | perl -pe 's/xXx/\n/g' | sed '/^$/d' > ${PREFIX}_paired.R1.aligned.flat.fasta
tornado_flatten_fasta.py -i ${PREFIX}_paired.R1.aligned.common.fasta -o ${PREFIX}_paired.R1.aligned.flat.fasta
#cat ${PREFIX}_paired.R2.aligned.common.fasta | sed '/^>/s/$/xXx/ ; /^>/s/^/xXx/'| tr -d '\n'| sed 's/xXx/\n/g' | sed '/^$/d' > ${PREFIX}_paired.R2.aligned.flat.fasta
#cat ${PREFIX}_paired.R2.aligned.common.fasta | perl -pe 's/$/xXx/ if /^>/; s/^/xXx/ if /^>/' | tr -d '\n' | perl -pe 's/xXx/\n/g' | sed '/^$/d' > ${PREFIX}_paired.R2.aligned.flat.fasta
tornado_flatten_fasta.py -i ${PREFIX}_paired.R2.aligned.common.fasta -o ${PREFIX}_paired.R2.aligned.flat.fasta

#and concatenate
echo "Concatenate read pairs"
#paste -d 'N' ${PREFIX}_paired.R1.aligned.flat.fasta ${PREFIX}_paired.R2.aligned.flat.fasta| sed '1~2 s/N>.*//'| sed '2~2 s/N//' > ${PREFIX}_paired.aligned.fasta
paste -d 'N' ${PREFIX}_paired.R1.aligned.flat.fasta ${PREFIX}_paired.R2.aligned.flat.fasta| awk 'NR%2 == 1 {sub(/N>.*/,""); print}; NR%2 == 0 {sub(/N/,""); print}' > ${PREFIX}_paired.aligned.fasta
echo "Make trees"
#now we can tree
OMP_NUM_THREADS=$NPROC $FASTTREE -nt -gtr -gamma -spr 4 -mlacc 2 -slownni -out ${PREFIX}_R1.tree ${PREFIX}_R1.aligned.fasta
OMP_NUM_THREADS=$NPROC $FASTTREE -nt -gtr -gamma -spr 4 -mlacc 2 -slownni -out ${PREFIX}_R2.tree ${PREFIX}_R2.aligned.fasta
OMP_NUM_THREADS=$NPROC $FASTTREE -nt -gtr -gamma -spr 4 -mlacc 2 -slownni -out ${PREFIX}_paired.tree ${PREFIX}_paired.aligned.fasta

#copy the trees to the results directory
cp *.tree ../$RESULTS

isitthere ../$RESULTS/${PREFIX}_R1.tree
isitthere ../$RESULTS/${PREFIX}_R2.tree
isitthere ../$RESULTS/${PREFIX}_paired.tree

if [[ -n $VSEARCH ]]
then
  echo "Map OTUs R1"
  $VSEARCH --threads $NPROC --usearch_global ${PREFIX}_R1.fasta --db ${PREFIX}_R1.otus.final.fasta --strand plus --id 0.97 --uc ${PREFIX}_R1.uc $VSEARCH_OPTS
  echo "Map OTUs R2"
  $VSEARCH --threads $NPROC --usearch_global ${PREFIX}_R2.fasta --db ${PREFIX}_R2.otus.final.fasta --strand plus --id 0.97 --uc ${PREFIX}_R2.uc $VSEARCH_OPTS
  echo "Map OTUs R1+R2"
  $VSEARCH --threads $NPROC --usearch_global ${PREFIX}_paired.fasta --db ${PREFIX}_paired.otus.final.fasta --strand plus --id 0.97 --uc ${PREFIX}_paired.uc $VSEARCH_OPTS
else
  #now map the reads into the OTU buckets
  echo "Map OTUs R1"
  #split the fasta into 1 GB chunks, 
  gt splitfasta -targetsize 1000 ${PREFIX}_R1.fasta
  #how many we have?
  COUNT=`ls ${PREFIX}_R1.fasta.*|wc -l`
  for i in $(seq $COUNT)
  do
  $USEARCH7 -threads $NPROC -usearch_global ${PREFIX}_R1.fasta.${i} -db ${PREFIX}_R1.otus.final.fasta -strand plus -id 0.97 -uc ${PREFIX}_R1.uc.${i}
  done
  #and merge the outputs
  cat ${PREFIX}_R1.uc.* > ${PREFIX}_R1.uc
  echo "Map OTUs R2"
  #split the fasta into 1 GB chunks, 
  gt splitfasta -targetsize 1000 ${PREFIX}_R2.fasta
  #how many we have?
  COUNT=`ls ${PREFIX}_R2.fasta.*|wc -l`
  for i in $(seq $COUNT)
  do
  $USEARCH7 -threads $NPROC -usearch_global ${PREFIX}_R2.fasta.${i} -db ${PREFIX}_R2.otus.final.fasta -strand plus -id 0.97 -uc ${PREFIX}_R2.uc.${i}
  done
  #and merge the outputs
  cat ${PREFIX}_R2.uc.* > ${PREFIX}_R2.uc
  echo "Map OTUs R1+R2"
  #split the fasta into 1 GB chunks, 
  gt splitfasta -targetsize 1000 ${PREFIX}_paired.fasta
  #how many we have?
  COUNT=`ls ${PREFIX}_paired.fasta.*|wc -l`
  for i in $(seq $COUNT)
  do
  $USEARCH7 -threads $NPROC -usearch_global ${PREFIX}_paired.fasta.${i} -db ${PREFIX}_paired.otus.final.fasta -strand plus -id 0.97 -uc ${PREFIX}_paired.uc.${i}
  done
  #and merge the outputs
  cat ${PREFIX}_paired.uc.* > ${PREFIX}_paired.uc
fi
mv ${PREFIX}_paired.tax.${TAXONOMY}.wang.taxonomy ${PREFIX}_paired.otus2.${TAXONOMY}.wang.taxonomy
#parse the OTU clusters
echo "Parsing OTU clusters..."

tornado_parse_otu_clusters.py ${PREFIX}_R1.otus.final.fasta ${PREFIX}_R1.uc ${PREFIX}_R1.otus.txt ${PREFIX}_R1.failures.txt
tornado_parse_otu_clusters.py ${PREFIX}_R2.otus.final.fasta ${PREFIX}_R2.uc ${PREFIX}_R2.otus.txt ${PREFIX}_R2.failures.txt
tornado_parse_otu_clusters.py ${PREFIX}_paired.otus.final.fasta ${PREFIX}_paired.uc ${PREFIX}_paired.otus.txt ${PREFIX}_paired.failures.txt

#if consensus taxonomy is set
if [ $CONSENSUS_TAXONOMY = "1" ]
then
#convert the otu txt file into a mothur list file
echo "Calculating consensus taxonomy"
#first, count the OTUs
OTU_COUNT=$(wc -l ${PREFIX}_R1.otus.txt|cut -f 1 -d ' ')
#and now convert the OTU files
cut -f 2- ${PREFIX}_R1.otus.txt | tr '\t' ','|tr '\n' ' ' | sed "s/^/consensus   ${OTU_COUNT} /" > ${PREFIX}_R1.otus.list
#do the same for the other two
OTU_COUNT=$(wc -l ${PREFIX}_R2.otus.txt|cut -f 1 -d ' ')
cut -f 2- ${PREFIX}_R2.otus.txt | tr '\t' ','|tr '\n' ' ' | sed "s/^/consensus   ${OTU_COUNT} /" > ${PREFIX}_R2.otus.list
OTU_COUNT=$(wc -l ${PREFIX}_paired.otus.txt|cut -f 1 -d ' ')
cut -f 2- ${PREFIX}_paired.otus.txt | tr '\t' ','|tr '\n' ' ' | sed "s/^/consensus   ${OTU_COUNT} /" > ${PREFIX}_paired.otus.list

#now calculate the consensus taxonomy
mothur "#classify.otu(list=${PREFIX}_R1.otus.list, taxonomy=${PREFIX}_R1.unique.${TAXONOMY}.wang.taxonomy, name=${PREFIX}_R1.names, reftaxonomy=${DATA}/${TAXONOMY}.taxonomy, basis=sequence)"
mothur "#classify.otu(list=${PREFIX}_R2.otus.list, taxonomy=${PREFIX}_R2.unique.${TAXONOMY}.wang.taxonomy, name=${PREFIX}_R2.names, reftaxonomy=${DATA}/${TAXONOMY}.taxonomy, basis=sequence)"
mothur "#classify.otu(list=${PREFIX}_paired.otus.list, taxonomy=${PREFIX}.padded.unique.${TAXONOMY}.wang.taxonomy, name=${PREFIX}.padded.names, reftaxonomy=${DATA}/${TAXONOMY}.taxonomy, basis=sequence)"

#construct taxonomy tables
#R1
#Nice tables for human-readable results
sed -n '2,$p' ${PREFIX}_R1.otus.consensus.cons.taxonomy| cut -f 2- > R1.tax.tmp
cut -f 1 ${PREFIX}_R1.otus.txt > R1.otu.tmp
echo -e 'OTU\tSize\tTaxonomy' > ${PREFIX}_R1.cons.probs.taxonomy
paste R1.otu.tmp R1.tax.tmp >> ${PREFIX}_R1.cons.probs.taxonomy
#Table for BIOM file
sed 's/([\.0-9]\{1,4\});/;/g' R1.tax.tmp | cut -f 2- > R1.tax2.tmp
#overwrite OTU taxonomy
paste R1.otu.tmp R1.tax2.tmp > ${PREFIX}_R1.otus2.${TAXONOMY}.wang.taxonomy
#R2
sed -n '2,$p' ${PREFIX}_R2.otus.consensus.cons.taxonomy| cut -f 2- > R2.tax.tmp
cut -f 1 ${PREFIX}_R2.otus.txt > R2.otu.tmp
echo -e 'OTU\tSize\tTaxonomy' > ${PREFIX}_R2.cons.probs.taxonomy
paste R2.otu.tmp R2.tax.tmp >> ${PREFIX}_R2.cons.probs.taxonomy
#Table for BIOM file
sed 's/([\.0-9]\{1,4\});/;/g' R2.tax.tmp | cut -f 2- > R2.tax2.tmp
#overwrite OTU taxonomy
paste R2.otu.tmp R2.tax2.tmp > ${PREFIX}_R2.otus2.${TAXONOMY}.wang.taxonomy
#paired
sed -n '2,$p' ${PREFIX}_paired.otus.consensus.cons.taxonomy| cut -f 2- > paired.tax.tmp
cut -f 1 ${PREFIX}_paired.otus.txt > paired.otu.tmp
echo -e 'OTU\tSize\tTaxonomy' > ${PREFIX}_paired.cons.probs.taxonomy
paste paired.otu.tmp paired.tax.tmp >> ${PREFIX}_paired.cons.probs.taxonomy
#Table for BIOM file
sed 's/([\.0-9]\{1,4\});/;/g' paired.tax.tmp |cut -f 2- > paired.tax2.tmp
#overwrite OTU taxonomy
paste paired.otu.tmp paired.tax2.tmp > ${PREFIX}_paired.otus2.${TAXONOMY}.wang.taxonomy
#copy human readable tables to results, overwrite previous
cp ${PREFIX}_R1.cons.probs.taxonomy ../$RESULTS/${PREFIX}_R1.probs.taxonomy
cp ${PREFIX}_R2.cons.probs.taxonomy ../$RESULTS/${PREFIX}_R2.probs.taxonomy
cp ${PREFIX}_paired.cons.probs.taxonomy ../$RESULTS/${PREFIX}_paired.probs.taxonomy

fi

#finally, make the biom tables
#Clean up taxonomy
sed 's/unclassified;//g' ${PREFIX}_R1.otus2.${TAXONOMY}.wang.taxonomy > tmp_R1.taxonomy
sed 's/unclassified;//g' ${PREFIX}_R2.otus2.${TAXONOMY}.wang.taxonomy > tmp_R2.taxonomy
sed 's/unclassified;//g' ${PREFIX}_paired.otus2.${TAXONOMY}.wang.taxonomy > tmp_paired.taxonomy
mv tmp_R1.taxonomy ${PREFIX}_R1.otus2.${TAXONOMY}.wang.taxonomy
mv tmp_R2.taxonomy ${PREFIX}_R2.otus2.${TAXONOMY}.wang.taxonomy
mv tmp_paired.taxonomy ${PREFIX}_paired.otus2.${TAXONOMY}.wang.taxonomy

echo "Making OTU tables..."
tornado_make_biom_table.py ${PREFIX}_R1.otus.txt ../$MAPPING ${PREFIX}_R1.otus2.${TAXONOMY}.wang.taxonomy ${PREFIX}_R1.biom
tornado_make_biom_table.py ${PREFIX}_R2.otus.txt ../$MAPPING ${PREFIX}_R2.otus2.${TAXONOMY}.wang.taxonomy ${PREFIX}_R2.biom
tornado_make_biom_table.py ${PREFIX}_paired.otus.txt ../$MAPPING ${PREFIX}_paired.otus2.${TAXONOMY}.wang.taxonomy ${PREFIX}_paired.biom

#copy the failures to the results
cp ${PREFIX}*failures.txt ../$RESULTS
#copy the OTUs to the results too 
cp ${PREFIX}_*otus.txt ../$RESULTS

echo "Copy results"
#and copy the bioms into the results file
cp ${PREFIX}*.biom ../$RESULTS/
#copy the mapping file there too
cp ../$MAPPING ../$RESULTS/
echo "Fix taxonomy"
tornado_clean_taxonomy.py ../$RESULTS/${PREFIX}_R1.otus.final.fasta ../$RESULTS/${PREFIX}_R1.probs.taxonomy ../$RESULTS/tmp_R1.taxonomy
tornado_clean_taxonomy.py ../$RESULTS/${PREFIX}_R2.otus.final.fasta ../$RESULTS/${PREFIX}_R2.probs.taxonomy ../$RESULTS/tmp_R2.taxonomy
tornado_clean_taxonomy.py ../$RESULTS/${PREFIX}_paired.otus.final.fasta ../$RESULTS/${PREFIX}_paired.probs.taxonomy ../$RESULTS/tmp_paired.taxonomy

mv ../$RESULTS/tmp_paired.taxonomy ../$RESULTS/${PREFIX}_paired.probs.taxonomy
mv ../$RESULTS/tmp_R1.taxonomy ../$RESULTS/${PREFIX}_R1.probs.taxonomy
mv ../$RESULTS/tmp_R2.taxonomy ../$RESULTS/${PREFIX}_R2.probs.taxonomy

#Compress the large fasta files
$GZIP ${PREFIX}_R1.fasta
$GZIP ${PREFIX}_R2.fasta
$GZIP ${PREFIX}_paired.fasta

mv *.gz ../$RESULTS/

echo "Cleaning up..."
if [ $CLEAN = 'all' ]
then
cd ..
rm -rf $WORKSPACE
elif [ $CLEAN = 'normal' ]
then
mv ${PREFIX}_R1.fasta ${PREFIX}_R1.fasta.save
mv ${PREFIX}_R2.fasta ${PREFIX}_R2.fasta.save
mv ${PREFIX}_paired.fasta ${PREFIX}_paired.fasta.save
rm -f *.fasta *.fasta.[0-9]* *.uc *.uc.* *.cfastq *accnos *biom *groups *logfile *map *scores *stk *summary *taxonomy *tree *txt *names *count_table
mv ${PREFIX}_R1.fasta.save ${PREFIX}_R1.fasta
mv ${PREFIX}_R2.fasta.save ${PREFIX}_R2.fasta
mv ${PREFIX}_paired.fasta.save ${PREFIX}_paired.fasta
elif [ $CLEAN = 'no' ]
then
echo "No cleanup performed."
fi
echo "Done"
