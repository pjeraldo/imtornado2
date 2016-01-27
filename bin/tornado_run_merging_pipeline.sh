#!/usr/bin/env bash
#author = "Patricio Jeraldo"
#copyright = "Copyright 2014-2016, Mayo Foundation for Medical Education and Research "
#credits = "Patricio Jeraldo,  Nicholas Chia"
#license = "GPL"
#version = "2.0.3.3"
#maintainer = "Patricio Jeraldo"
#email = "jeraldo.patricio@mayo.edu"
#status = "Alpha"

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
NSAMPLES=$(cat $WORKSPACE/ids.txt| wc -l)
#count the number of fastq files in the directory. Assume paired end. COunt only R1.
NREAD1=$(ls *${SPACER}*R1*fastq|wc -l)
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


#Merge reads with USEARCH7

for name in $(cat $WORKSPACE/ids.txt)
do
$USEARCH7 -fastq_mergepairs ${name}${SPACER}*R1*fastq -reverse ${name}${SPACER}*R2*fastq -fastq_truncqual 3 -fastq_minovlen $MINIMUM_OVERLAP -minhsp $MINIMUM_OVERLAP -fastqout $WORKSPACE/${name}_M.fastq
done

#We won't need to be somewhere else until the very end
cd $WORKSPACE
#Convert FASTQ to FASTA
echo "Converting FASTQ to FASTA, and rename IDs according to sample name"
for name in $(cat ids.txt)
do
tornado_convert_and_rename_ids.py ${name}_M.fastq ${name}_M.fasta
done

#Quality filtering (according to EMP protocol). Q3 trimming, 75% of length as cutoff, no ambigs, nothing on homopolymers
echo "Remove ambigs..."
#no check for polys... let OTUing tacke care of this
#there has to be a faster way to do this.
for name in $(cat ids.txt)
do
tornado_remove_ambigs.sh ${name}_M.fasta ${name}_M.trim.fasta
done

#We do not need to trim on this merging pipeline by definition
#only merging the different reads
#No groups file
cat *M.trim.fasta > ${PREFIX}_M.fasta
isitthere ${PREFIX}_M.fasta

#calculate overall taxonomy if consensus is set
if [ $CONSENSUS_TAXONOMY = "1" ]
then
#unique samples
mothur "#classify.seqs(fasta=${PREFIX}_M.unique.fasta name=${PREFIX}_M.names, taxonomy=${DATA}/${TAXONOMY}.taxonomy, template=${DATA}/${TAXONOMY}.fna, probs=false, processors=${NPROC}, iters=100, cutoff=${TAXCUTOFF})"
fi

if [[ -n $VSEARCH ]]
then
  $VSEARCH -derep_fulllength ${PREFIX}_M.fasta -output ${PREFIX}_M.derep.fasta -sizeout
else
  #OTU R1 and R2, separately
  #dereplicate
  #use mothur for this for now
  mothur "#unique.seqs(fasta=${PREFIX}_M.fasta)"

  #get the "counts" .. returns PREFIX_R?.seq.count
  mothur "#count.seqs(name=${PREFIX}_M.names)"

  #annotate the uniques with the sizes
  #and remove reads shorter than the specified trims
  echo "Annotating unique sizes..."
  #changing minimum length to 1 to not exclude any reads, for now... can leave up to user
  tornado_annotate_read_sizes.py 1 ${PREFIX}_M.unique.fasta ${PREFIX}_M.count_table ${PREFIX}_M.derep.fasta
fi

#remove singletons
$USEARCH7 -sortbysize ${PREFIX}_M.derep.fasta -output ${PREFIX}_M.derep2.fasta -minsize 2

#OTU the reads!
$USEARCH7 -cluster_otus ${PREFIX}_M.derep2.fasta -otus ${PREFIX}_M.otus.fasta

echo "Rename otu ids..."
#rename the OTU ids...

tornado_otu_renamer.py ${PREFIX}_M.otus.fasta ${PREFIX}_M.otus2.fasta

#get the taxonomy of the representatives
mothur "#classify.seqs(fasta=${PREFIX}_M.otus2.fasta, taxonomy=${DATA}/${TAXONOMY}.taxonomy, template=${DATA}/${TAXONOMY}.fna, processors=${NPROC}, iters=100, cutoff=${TAXCUTOFF})"

#copy the taxonomy files to the results
cp ${PREFIX}_M.otus2.${TAXONOMY}.wang.taxonomy ../$RESULTS/${PREFIX}_M.probs.taxonomy

#do this again, with no support values
mothur "#classify.seqs(fasta=${PREFIX}_M.otus2.fasta, taxonomy=${DATA}/${TAXONOMY}.taxonomy, template=${DATA}/${TAXONOMY}.fna, probs=false, processors=${NPROC}, iters=100, cutoff=${TAXCUTOFF})"

echo "Filter bad reads..."
#find the ones that are obviously not 16S
grep "unknown" ${PREFIX}_M.otus2.${TAXONOMY}.wang.taxonomy| cut -f 1 > bad_taxonomy_M.accnos

#And remove those reads from the otu file
tornado_read_remover.py bad_taxonomy_M.accnos ${PREFIX}_M.otus2.fasta ${PREFIX}_M.otus3.fasta

#align the reads... there may be still some bad reads in there
#probably this will need to be done using STK instead of fasta
echo "Align reads"
cmalign --cpu $NPROC -g --notrunc --sub --dnaout --noprob --sfile ${PREFIX}_M.scores -o ${PREFIX}_M.otus3.aligned.stk $DATA/$MODEL ${PREFIX}_M.otus3.fasta


#convert stk to fasta
echo "Convert to fasta"

tornado_stk2fasta.py ${PREFIX}_M.otus3.aligned.stk ${PREFIX}_M.otus3.aligned.fasta


echo "Filter bad reads"
#look for the badly aligned reads
#for the merging pipeline, this can be improved by looking at alignment start and ends
tornado_overlap_score_filter.py ${PREFIX}_M.scores bad_align_M.accnos

#convert all lowercase bases to uppercase
awk '{print (substr($0,0,1) == ">")? $0 : toupper($0)}' ${PREFIX}_M.otus3.fasta > ${PREFIX}_M.otus4.fasta

#remove them from the OTU file 
tornado_read_remover.py bad_align_M.accnos ${PREFIX}_M.otus4.fasta ${PREFIX}_M.otus.final.fasta

#copy the final OTUs into the results file
#because they are actually a results
cp ${PREFIX}_M.otus.final.fasta ../$RESULTS

#and from the alignment
tornado_read_remover.py bad_align_M.accnos ${PREFIX}_M.otus3.aligned.fasta ${PREFIX}_M.otus3.aligned.clean.fasta

echo "Remove gapped columns"
#remove gapped columns from the alignment
tornado_remove_gaps.py ${PREFIX}_M.otus3.aligned.clean.fasta ${PREFIX}_M.otus3.aligned.clean2.fasta

#change periods to dashes in alignment, as well as change letters to uppercase
cat ${PREFIX}_M.otus3.aligned.clean2.fasta| tr '.' '-' | tr 'a-z' 'A-Z' > ${PREFIX}_M.aligned.fasta

#copy the alignments too
cp ${PREFIX}_M.aligned.fasta ../$RESULTS

echo "Make trees"
#now we can tree
OMP_NUM_THREADS=$NPROC $FASTTREE -nt -gtr -gamma -spr 4 -mlacc 2 -slownni -out ${PREFIX}_M.tree ${PREFIX}_M.aligned.fasta

#copy the trees to the results directory
cp *.tree ../$RESULTS

isitthere ../$RESULTS/${PREFIX}_M.tree

#now map the reads into the OTU buckets
echo "Map OTUs"

if [[ -n $VSEARCH ]]
then
  $VSEARCH --threads $NPROC --usearch_global ${PREFIX}_M.fasta --db ${PREFIX}_M.otus.final.fasta --strand plus --id 0.97 --uc ${PREFIX}_M.uc $VSEARCH_OPTS
else
  #split the fasta into 1 GB chunks, 
  gt splitfasta -targetsize 1000 ${PREFIX}_M.fasta
  #how many we have?
  COUNT=$(ls ${PREFIX}_M.fasta.*|wc -l)
  for i in $(seq $COUNT)
  do
  $USEARCH7 -threads $NPROC -usearch_global ${PREFIX}_M.fasta.${i} -db ${PREFIX}_M.otus.final.fasta -strand plus -id 0.97 -uc ${PREFIX}_M.uc.${i}
  done
  #and merge the outputs
  cat ${PREFIX}_M.uc.* > ${PREFIX}_M.uc
fi
#parse the OTU clusters
echo "Parsing OTU clusters..."

tornado_parse_otu_clusters.py ${PREFIX}_M.otus.final.fasta ${PREFIX}_M.uc ${PREFIX}_M.otus.txt ${PREFIX}_M.failures.txt

#if consensus taxonomy is set
if [ $CONSENSUS_TAXONOMY = "1" ]
then
#convert the otu txt file into a mothur list file
echo "Calculating consensus taxonomy"
#first, count the OTUs
OTU_COUNT=$(wc -l ${PREFIX}_M.otus.txt|cut -f 1 -d ' ')
#and now convert the OTU files
cut -f 2- ${PREFIX}_M.otus.txt | tr '\t' ','|tr '\n' ' ' | sed "s/^/consensus   ${OTU_COUNT} /" > ${PREFIX}_M.otus.list


#now calculate the consensus taxonomy
mothur "#classify.otu(list=${PREFIX}_M.otus.list, taxonomy=${PREFIX}_M.unique.${TAXONOMY}.wang.taxonomy, name=${PREFIX}_M.names, reftaxonomy=${DATA}/${TAXONOMY}.taxonomy, basis=sequence)"

#construct taxonomy tables
#R1
#Nice tables for human-readable results
sed -n '2,$p' ${PREFIX}_M.otus.consensus.cons.taxonomy| cut -f 2- > M.tax.tmp
cut -f 1 ${PREFIX}_M.otus.txt > M.otu.tmp
echo -e 'OTU\tSize\tTaxonomy' > ${PREFIX}_M.cons.probs.taxonomy
paste M.otu.tmp M.tax.tmp >> ${PREFIX}_M.cons.probs.taxonomy
#Table for BIOM file
sed 's/([\.0-9]\{1,4\});/;/g' M.tax.tmp | cut -f 2- > M.tax2.tmp
#overwrite OTU taxonomy
paste M.otu.tmp M.tax2.tmp > ${PREFIX}_M.otus2.${TAXONOMY}.wang.taxonomy

#copy human readable tables to results
cp ${PREFIX}_M.cons.probs.taxonomy ../$RESULTS/${PREFIX}_M.probs.consensus.taxonomy

fi

echo "Making OTU tables..."
tornado_make_biom_table.py ${PREFIX}_M.otus.txt ../$MAPPING ${PREFIX}_M.otus2.${TAXONOMY}.wang.taxonomy ${PREFIX}_M.biom

#copy the failures to the results
cp ${PREFIX}*failures.txt ../$RESULTS
#copy OTUs to results
cp ${PREFIX}*otus.txt ../$RESULTS

#compress clean reads
$GZIP ${PREFIX}_M.fasta
#and move them to the results area
mv ${PREFIX}_M.fasta ../$RESULTS

echo "Copy results"
#and copy the bioms into the results file
cp ${PREFIX}*.biom ../$RESULTS/
#copy the mapping file there too
cp ../$MAPPING ../$RESULTS/
echo "Fix taxonomy"
tornado_clean_taxonomy.py ../$RESULTS/${PREFIX}_M.otus.final.fasta ../$RESULTS/${PREFIX}_M.probs.taxonomy ../$RESULTS/tmp_M.taxonomy

mv ../$RESULTS/tmp_M.taxonomy ../$RESULTS/${PREFIX}_M.probs.taxonomy

echo "Cleaning up..."
if [ $CLEAN = 'all' ]
then
cd ..
rm -rf $WORKSPACE
elif [ $CLEAN = 'normal' ]
then
mv ${PREFIX}_M.fasta ${PREFIX}_M.fasta.save
rm -f *.fasta *.fasta.[0-9]* *.uc *.uc.* *.cfastq *accnos *biom *groups *logfile *map *scores *stk *summary *taxonomy *tree *txt *names *count_table
mv ${PREFIX}_M.fasta.save ${PREFIX}_M.fasta
elif [ $CLEAN = 'no' ]
then
echo "No cleanup performed."
fi
echo "Done"
