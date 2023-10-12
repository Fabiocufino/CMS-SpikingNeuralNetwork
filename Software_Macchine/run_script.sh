#! /usr/bin/bash
start_index=0
last_index=5
for((c=$start_index;c<=last_index;c++))
do
    nohup root -l -b < ./starting_parameters/start_$c.cmd >&./out/log_$c.txt 2>&1 &
done