#! /usr/bin/bash
#path=../src/examples/
declare -a db=("shortest passive" "shortest active" "shortestrestless active" "shortestrestless passive" "shortestforemost passive")

for args in $*
do
    for i in "${db[@]}";
    do
        a=( $i );
        timeout 7200s ./btwBenchmark -f ${args} -n 1 -b 1 -t 1 -c "${a[0]}" -y "${a[1]}"
    done

    for args in $*
    do
        for j in 90 80 70 60 50 40 30 20 10
        do
            for i in "${db[@]}";
            do
                a=( $i );
                timeout 7200s ./btwBenchmark -f ${args} -n 1 -b 1 -p $j -c "${a[0]}" -y "${a[1]}"
            done

        done
    done
done


