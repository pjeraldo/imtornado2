# IM-TORNADO: A pipeline for 16S reads from paired-end libraries

This is the documentation for the IM-TORNADO pipeline, currently at version 2.0.3

Documentation will be available in [Read the Docs](http://readthedocs.org) soon.

## Changelog

### Version 2.0.3.1

* Bigfix - Minor problem when cleaning after the pipeline corrected now.
* Bugfix - Correct README now included in the tarball

### Version 2.0.3

* Bugfix - Fixed issue with large datasets not being processed properly, hopefully preserving POSIX compliance
* Bugfix - Fixed issue with empty files when processing files with low quality R2 reads (leaving no paired reads but some R1 and R2 reads)
* Pipeline now outputs the clean grouped reads (for use in other tools, e.g. picrust or for closed-reference OTU)
* Added support for [VSEARCH](https://github.com/torognes/vsearch) when appropriate. If using 64-bit USEARCH or VSEARCH, processing should be faster for very large datasets. See rest of documentation for more information

### Version 2.0.2

* **Pipeline runs in Apple OSX now** (tested only in Mavericks). Pipeline is now POSIX-compliant-enough to properly run in OSX. Instructions to follow soon
* Pipeline supports BIOM format APIs 1.0, 2.0 and 2.1, so it can run regardless of what version of biom-format is installed (as long as the API doesn't change again). The BIOM output is still JSON only
* Minor fixes in the dependency check scripts
* Other minor fixes

### Version 2.0.1

* Fix hardcoded parameters
* Minor bugfixes

## Requirements

IM_TORNADO was designed to run on Linux systems. It can run in OSX too (tested only in Mavericks). Windows is not supported (and it may be difficult or impossible to support).

IM-TORNADO depends on the following programs to function:

* [USEARCH](http://drive5.com/usearch "USEARCH website") (version 7.0 series; version 8.0 not tested)
* [mothur](http://mothur.org "mothur website") (version 1.28 and above)
* [FastTree](http://www.microbesonline.org/fasttree "FastTree website")
* [genometools](http://genometools.org "genometools website")
* [Infernal](http://infernal.janelia.edu "Infernal website") (version 1.1)
* [Trimmomatic](http://usadellab.org/cms/?page=trimmomatic "Trimmomatic website") (version 0.28 and above)
* python (version 2.7, version 3 and higher not supported)
* java (for Trimmomatic to work)
* perl, sed and awk
* [genometools](http://genometools.org "genometools website") (not required if using 64-bit USEARCH or VSEARCH in version 2.0.3 of the pipeline)
* [VSEARCH](https://github.com/torognes/vsearch) (optional for 2.0.3 and above, not really needed if using 64-bit USEARCH)

And the following python libraries are required:

* [Biopython](http://biopython.org "Biopython website")
* [BIOM-format](http://biom-format.org "BIOM-format website") (version 1.0 and above)
* [bitarray](https://pypi.python.org/pypi/bitarray "bitarray package website")

Before installing, please make sure all software is installed and their commands in the PATH, and the respective python libraries available to the interpreter.

Be careful of software you may have installed that bundles their own python interpreter, as they will most likely conflict with IM-TORNADO's operation due to missing libraries. A known example is QIIME. They can coexist, but attention is needed.

## Installation

These are the minimal install instructions. 

* Unpack the package into a directory of your choice, for example /home/yourusername. The directory IM-TORNADO-<version> will be created.
* IMPORTANT: under the IM-TORNADO-<version> directory, go to the scripts directory and edit the file "tornado-params.sh"
* In the "tornado-params.sh" file, go to the USEARCH7 variable, and enter the name of the USEARCH executable (for example, usearch7.0.1090_i86linux32). Do the same with the FASTREE variable.
* For the variable TRIMMOMATIC, indicate the location of the JAR file of Trimmomatic. For example "/opt/Trimmomatic-0.32/trimmomatic-0.32.jar"
* Finally, for the TORNADO2 variable, indicate the location of the main directory of the pipeline. In this case, /home/<yourusername>/IM-TORNADO-<version>
* Run the "check_deps.sh" script as "./check_deps.sh", to check for missing dependencies. If all goes well, all dependencies should be recognized.
* Add the bin directory to the PATH, so the commands can be recognized: "PATH=/home/<yourusername>/IM-TORNADO-<version>:$PATH"
* If using 64-bit USEARCH or if using VSEARCH, in the params file enter the name of the corresponding executable. For VSEARCH only, enter any neccesary extra options in the VSEARCH_OPTS variable.

Of course, you can choose to install it in a different directory. Just make sure you have write permissions if you want to do so.

If you install it in a directory that most users don't have write permissions, **run the pipeline with the example data at least once**, so the taxonomy database index is created by mothur. This has to be done for other new or custom taxonomies, and possibly for every new version of mothur.

## Running the pipeline

_These instructions are being expanded_

After installing you can run the pipeline:

* Create a directory where your data is going to be located. Copy all the fastq files into that directory.
* Copy your metadata/mapping file into your work directory. Only the fastq files referenced in the metadata/mapping file must be in the directory. If there are more (or fewer) files than what is declared in the metadata file, then the pipeline will throw an error.
* Copy the "tornado-params.sh" from the scripts directory in the pipeline.
* Edit the "tornado-params.sh" file, modify the run parameters accordingly, such as the prefix names, read lengths, taxonomy to use, maximum number of processors to use, spacer character, and output directories.
* To run the pipeline, execute "tornado_run_pipeline.sh"
* If all goes well, after much text scrolling, the pipeline will finish with no errors.

### What if there is an error?

Usually, errors are due to missing files or misnamed samples.

* Check your file names and sample names declared in the metadata file, and make sure they match.
* Make sure the sample names have no reserved characters. In particular, **the underscore "_" is not allowed.** It will cause problems in the pipeline, and most likely in downstream analysis with other software such as QIIME.

## Pipeline output

The pipeline will produce a set of files, all stored in the results directory:

* .biom files, BIOM-formatted OTU tables. You can use these directly many specialized packages, or easily converted into tables.
* .tree files. Newick-formatted tree files for the OTU reprentatives. Useful for calculating metrics such as UniFrac.
* .taxonomy files. Text based taxonomies with confidence scores.
* .final.fasta files. Fasta files of the OTU representatives.
* .aligned.fasta files. Fasta formatted multiple sequence alignments of the OTU representatives.
* .failures.txt files. List of reads that failed to map to the OTU representatives. These most likely are low quality reads, potentially chimeric sequences or singletons far from any known OTU.
* Your metadata/mapping file.

The pipeline deletes most of its intermediate files, except for the merged fasta files with the clean reads from all samples. These are stored in the workspace directory. These files can be used as input to software such as _picrust_ (for this case, the file must be pre-processed through _QIIME_ first. See the _picrust_ tutorial).

### What happens next?

Now you can use your BIOM, tree and metadata files for proper statistical analysis, testing and discovery. _QIIME_'s "core\_diversity\_analysis.py" script can give you a nice overview of your data. You can also try _phyloseq_ if you are familiar with R, or use the _metagenassist.ca_ web interface for analysis.

## A note on VSEARCH

VSEARCH is a relatively new package designed to be a drop-in replacement for some core function of the USEARCH program. Most importantly, it is not hindered by the 32-bit limitations of the free version of USEARCH. It can be used in the pipeline for two steps: dereplication and mapping. It definitely improves running times for large datasets, and the code that manages it is easier to maintain. If you have a license for the 64-bit version, you can now take full advantage of this speedup: just enter the name of the USEARCH binary in the params file in the VSEARCH variable

Although the dereplication step give identical results, the mapping (usearch_global) step gives slightly different results, I presume due to different heuristics being used to find matches. In my short experience with VSEARCH, it seems that for each read being searched, it generally chose the same or slightly better (higher similarity) OTU bin to add the read to. For example:

	< H     643     207     97.1    +       0       0       207M43I ICF-10_10262    680
	< H     791     249     97.6    +       0       0       96MD152M2I      ICF-10_10263    832
	< H     645     250     99.2    +       0       0       249MD   ICF-10_10264    682
	< H     761     250     97.2    +       0       0       250M    ICF-10_10265    801
	---
	> H     149     207     97.6    +       0       0       207M43I ICF-10_10262    150
	> H     28      249     98.8    +       0       0       249MI   ICF-10_10263    29
	> H     641     250     99.2    +       0       0       249MD   ICF-10_10264    677
	> H     803     250     98.0    +       0       0       250M    ICF-10_10265    843

In this example, the top four lines are the choices made by USEARCH and the bottom lines are the choices made by VSEARCH. In this case, for three of the reads (numbers 1, 3 and 4) it chose higher similarity (4th column) OTU bins (OTU IDs are the last column). I've seen a only few cases where it chose more poorly. So for now it seems to be performing somewhat better. Feel free to experiment. YMMV.

## How to cite

If you use IM-TORNADO for your project, please cite the following manuscript:

Jeraldo P, Kalari K, Chen X, Bhavsar J, Mangalam A, White B, *et al*. *IM-TORNADO: A Tool for Comparison of 16S Reads from Paired-End Libraries.* PLOS ONE 9 (12):e114804. Available from: [http://dx.plos.org/10.1371/journal.pone.0114804](http://dx.plos.org/10.1371/journal.pone.0114804)
