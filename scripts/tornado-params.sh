## IM-TORNADO parameters. REVIEW BEFORE RUNNING THE PIPELINE.

#Prefix for the output files.
PREFIX=test

#Name of the mapping file
MAPPING=mapping.file.txt

#Model for alignment of reads. Allowed values are 'bacteria.cm' and 'archaea.cm'
#Custom models allowed if placed in the DATA directory (see end of configuration file)
#An example custom model would be a model for ITS
MODEL=bacteria.cm

#Minimum length of reads to keep after quality filtering.
#Suggested: 75% of sequence length (112 for 150bp reads, 187 for 250 bp reads)
#Not used for merging pipeline
MINIMUM_LENGTH=187

#Taxonomy database to use.
#Currently we provide 'rdp9' (RDP database, training set 9)
#Assignments are not reliable below genus level.
#For other taxonomies, download them from the mothur website.
#See the INSTALL.md document for instructions.
#NAme files as TAXONOMY.fna and TAXONOMY.taxonomy,
#where TAXONOMY is defined below
TAXONOMY=rdp9

#Trimming length for OTU steps... YMMV here. These are suggestions
#for 500 cycle Illumina kits.
#Not used for merging pipeline
R1_TRIM=250
R2_TRIM=200

#Set this to 1 if you wish to run consensus taxonomy 
#instead of taxonomy of the OTU representative.
#Not tested extensively.
#This will take significantly longer, and will use considerable more resources.
#Otherwise, leave as 0.
CONSENSUS_TAXONOMY=0

#Merging pipeline parameters
#keep this one larger than 8
#make this larger if you know your construct allows for a larger overlap
MINIMUM_OVERLAP=16
#Directory for files created during the pipeline run
WORKSPACE=workspace/

#Directory to place resulting files.
RESULTS=results/

#Maximum number of processors to be requested. biocluster has a maximum of 8.
NPROC=4

#Spacer for the sample name in the sample file.
#Common spacers are '.' and '_'
#Example: SAMPLE1.FC000_R1_L001.fastq ('.' spacer)
#Example: SAMPLE1_FC000_R1_L001.fastq ('_' spacer)
SPACER='.'

#Clean up the workspace?
#'all' => deletes everything in the workspace after completion
#'no' => do not delete anything, useful for debugging
#'normal' => delete everything but the grouped fasta files after cleanup (useful for, e.g., picrust)
CLEAN=normal

#Name of the usearch version 7 binary. Could be the actual name,
#or some other name. It needs to be in the PATH
USEARCH7=usearch7

#Name of the Fast Tree executable
FASTTREE=FastTreeMP

#PATH to Trimmomatic jar file
TRIMMOMATIC=/path/to/trimmomatic-0.30.jar

#Advanced settings. Do not modify unless you know what you are doing.
TORNADO2=/path/to/IMTORNADO-2.0.0
DATA=$TORNADO2/data
