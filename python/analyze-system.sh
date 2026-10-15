#!/bin/bash

shopt -s nullglob

cd ../build

for dir in /scratch/user/cameron.parker/projects/runs/new/$1/*/*/;
    do echo "$dir";
    ./$1GammaAnalysis "$dir"
    ./$1GammaAnalysis "$dir" --abson

done

cd ../python