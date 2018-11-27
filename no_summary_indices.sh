#!/bin/bash
if [ ! -d $1 ]; then
    echo "Error: fast5 directory not found"
    exit 1
fi

	PRG="$BASH_SOURCE"

    while [ -h "$PRG" ] ; do
       ls=`ls -ld "$PRG"`
       link=`expr "$ls" : '.*-> \(.*\)$'`
       if expr "$link" : '/.*' > /dev/null; then
          PRG="$link"
       else
          PRG=`dirname "$PRG"`"/$link"
       fi
    done

dir=$(dirname "$PRG")

python $dir/scripts/generate_indices.py -d $1 > $2_tmp_file
sort -k1 $2_tmp_file | gzip > $2_fast5_index.gz
rm $2_tmp_file