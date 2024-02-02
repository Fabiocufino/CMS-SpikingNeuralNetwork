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
    wget -o - -O Code/Data/100k_100br.root "https://www.dropbox.com/scl/fi/dlauefynunvznx06kfol7/100k_100br.root?rlkey=qhmh2og38flpzo7spwmlsjm6w&dl=0"
fi

if ! [[ -e "Code/Data/ordered.root" ]];then
    wget -o - -O Code/Data/ordered.root "https://www.dropbox.com/scl/fi/2ipkrkxud5k9j7hglh64c/ordered.root?rlkey=8qrtubtqdedszypb6wd6o7d77&dl=0"
fi

echo "Execution terminated"