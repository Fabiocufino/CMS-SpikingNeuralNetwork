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
if ! [[ -e "Code/MODE/JSON" ]];then
    mkdir Code/MODE/JSON
fi
if ! [[ -e "Code/Data" ]];then
    mkdir Code/Data
fi

if ! [[ -e "Code/Data/delays.txt" ]];then
    touch Code/Data/delays.txt
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

if ! [[ -e "Code/Data/muons_100k_100br_new.root" ]];then
    wget -o - -O Code/Data/muons_100k_100br_new.root "https://www.dropbox.com/scl/fi/w7rgcgnv4pd9gusy70yt5/muons_100k_100br_new.root?rlkey=tmxrj8ckqjgcaidwdssph9w2e&st=zbch3kpb&dl=0"
fi

if ! [[ -e "Code/Data/muons_100k_100br.root" ]];then
    wget -o - -O Code/Data/muons_100k_100br.root "https://www.dropbox.com/scl/fi/i9o2525gz1t7buyxcl028/muons_100k_100br.root?rlkey=a8r54plfuy0r0gb8a0mc83vgc&dl=0"
fi

if ! [[ -e "Code/Data/01-02muons_100k_100br.root" ]];then
    wget -o - -O Code/Data/01-02muons_100k_100br.root "https://www.dropbox.com/scl/fi/atpq7ap1re04vvp38578r/01-02muons_100k_100br.root?rlkey=um54u3iz8kaeh77jpnjdbiab2&st=5zqb97m8&dl=0"
fi

if grep -q "SNN_PATH=" ~/.bashrc; then
    sed -i "s|export SNN_PATH=.*|export SNN_PATH=$(pwd)|g" ~/.bashrc
    source ~/.bashrc
    echo "SNN_PATH updated to current directory."
else
    echo "export SNN_PATH=$(pwd)" >> ~/.bashrc
    source ~/.bashrc
    echo "SNN_PATH is set to current directory."
fi


if ! [[ -e "Code/Data/muons_100k_100br_new.root" ]];then
    wget -o - -O Code/Data/muons_100k_100br_new.root "https://www.dropbox.com/scl/fi/w7rgcgnv4pd9gusy70yt5/muons_100k_100br_new.root?rlkey=tmxrj8ckqjgcaidwdssph9w2e&st=zbch3kpb&dl=0"
fi

echo "Execution terminated"