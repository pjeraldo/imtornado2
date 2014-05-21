## Installation

These are the minimal install instructions. 

* Unpack the package into a directory of your choice, for example /home/yourusername. The directory IM-TORNADO-<version> will be created.
* IMPORTANT: under the IM-TORNADO-<version> directory, go to the scripts directory and edit the file "tornado-params.sh"
* In the "tornado-params.sh" file, go to the USEARCH7 variable, and enter the name of the USEARCH executable (for example, usearch7.0.1090_i86linux32). Do the same with the FASTREE variable.
* For the variable TRIMMOMATIC, indicate the location of the JAR file of Trimmomatic. For example "/opt/Trimmomatic-0.32/trimmomatic-0.32.jar"
* Finally, for the TORNADO2 variable, indicate the location of the main directory of the pipeline. In this case, /home/<yourusername>/IM-TORNADO-<version>
* Run the "check_deps.sh" script as "./check_deps.sh", to check for missing dependencies. If all goes well, all dependencies should be recognized.
* Add the bin directory to the PATH, so the commands can be recognized: "PATH=/home/<yourusername>/IM-TORNADO-<version>:$PATH"

Of course, you can choose to install it in a different directory. Just make sure you have write permissions if you want to do so.

If you install it in a directory that most users don't have write permissions, **run the pipeline with the example data at least once**, so the taxonomy database index is created by mothur. This has to be done for other new or custom taxonomies, and possibly for every new version of mothur.