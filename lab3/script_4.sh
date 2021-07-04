# !/bin/bash

NFrom=880
NTo=20000

let "step=($NTo - $NFrom) / 10"
let "NTo=$NTo + $step"

gcc -fopenmp lab1_schedule_1.c -o lab3_schd1 -lm

for ((index=NFrom; index<NTo; ))
do

echo $index

./lab3_schd0 $index
cat result.txt >> result_s6.txt

let "index=index+(step)"

done