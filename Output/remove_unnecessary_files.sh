#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
for i in {0..19..1}
do
    echo "Entering folder $i:"
    mv iteration_$i/u000101.vtu ./ && rm iteration_$i/* && mv u000101.vtu iteration_$i/ && ls iteration_$i    
done
