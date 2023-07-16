#! /usr/bin/bash
#path=../src/examples/
#db = primaryschool.csv highschool_2012.csv highschool_2011.csv hospital_ward.csv ht09.csv workplace_2013.csv
for args in $*
do
    timeout 7200s ./btwBenchmark -f ${args} -n 1 -b 1 -t 1
done

for args in $*
      do
         for j in 90 80 70 60 50 40 30 20 10
         do
             timeout 7200s ./btwBenchmark -f ${args} -n 1 -b 1 -p $j
         done
      done

