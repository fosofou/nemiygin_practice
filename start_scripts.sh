#!/bin/bash

output_file=results.txt
rm $output_file

g++ -fopenmp  main.cpp -o main.out

./main.out | tee -a $output_file

rm main.out
