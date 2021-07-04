# !/bin/bash

NFrom=880
NTo=20000

let "step=($NTo - $NFrom) / 10"
let "NTo=$NTo + $step"

gcc -fopenmp lab1_schedule_0.c -o lab3_schd0 -lm

for ((index=NFrom; index<NTo; ))
do

echo $index

./lab3_schd0 $index
cat result.txt >> result_s0.txt

let "index=index+(step)"

done