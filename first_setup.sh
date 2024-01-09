#!/bin/bash

echo "Checking the directories structure"

if ! [[ -e "Code/MODE" ]];then
    mkdir Code/MODE
fi
if ! [[ -e "Code/MODE/SNNT" ]];then
    mkdir Code/MODE/SNNT
fi
if ! [[ -e "Code/MODE/CSV" ]];then
    mkdir Code/MODE/SNNT
fi
if ! [[ -e "Code/Data" ]];then
    mkdir Code/Data
fi
if ! [[ -e "Code/pdf" ]];then
    mkdir Code/pdf
fi
if ! [[ -e "Code/MODE/potentials.csv" ]];then
    touch "Code/MODE/potentials.csv"
fi

echo "Collecting the data files if missing... the process could take some minutes"

if ! [[ -e "Code/Data/100k_100br.root" ]];then
    wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1B6d3bxFBziWGIMmCp8uxxSF2As3k89m7' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1B6d3bxFBziWGIMmCp8uxxSF2As3k89m7" -O Code/Data/100k_100br.root && rm -rf /tmp/cookies.txt
fi

if ! [[ -e "Code/Data/ordered.root" ]];then
    wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=17YxLI8qBCfBDKPHVCa-fAKVGeLocOf09' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=17YxLI8qBCfBDKPHVCa-fAKVGeLocOf09" -O Code/Data/ordered.root && rm -rf /tmp/cookies.txt
fi

echo "Execution terminated"