#! /usr/bin/bash
#substitute 27 with the number of combinations that you want to run in parrallel
for((c=1;c<=27;c++))
do
    nohup root -l -b < ./starting_parameters/start_$c.cmd > &./out/log_$c.txt 2>&1 &
done