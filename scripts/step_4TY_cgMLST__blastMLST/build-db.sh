#!/bin/bash
set -euo pipefail

if [ -z "$1" ] ;
    then
        echo -e "\n usage: get-db <dbname>\n"
        exit 1
fi

cd /opt/conda/db/pubmlst &&\
tar xf "$1" &&\
/opt/conda/scripts/mlst-make_blast_db