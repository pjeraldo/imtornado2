#!/usr/bin/env bash

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

source tornado-params.sh

#CHECK FOR DEPENDENCIES
echo -e "Dependency checks"
function deps {
#dependency checking function
#Credit: http://www.snabelb.net

        DEPENDENCIES="mothur cmalign gt python2.7 $USEARCH7 java sed awk $FASTTREE"
 
        deps_ok=YES
        for dep in $DEPENDENCIES
        do
                if ! which $dep &>/dev/null;  then
                        echo -e "This script requires $dep to run but it is not installed"
                        deps_ok=NO
		else
		    echo -e "Dependency $dep found!"
                fi
        done
        if [[ "$deps_ok" == "NO" ]]; then
                echo -e "Unmet dependencies ^"
                echo -e "Aborting!"
                exit 1
        else
                return 0
        fi
}

#check some deps deps
deps

#check for python deps
$TORNADO2/bin/tornado_check_python_deps.py

#check for Trimmomatic's existence
echo -e "Checking for Trimmomatic"
[ -s $TRIMMOMATIC ] || die 1 "This script requires Trimmomatic to run but it is not installed"

echo -e "Done with dependency checks."
