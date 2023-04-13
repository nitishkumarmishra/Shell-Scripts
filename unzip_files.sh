#!/bin/bash

for fq1 in ls *.txt.gz; do

        input=$(basename $fq1)
        output=$(basename $fq1 .gz)
        echo -e "$input"
        echo -e "$output"
        gunzip -c $input > $output
done

