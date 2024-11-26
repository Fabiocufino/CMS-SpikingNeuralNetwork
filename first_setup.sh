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

if ! [[ -e "Code/Data/ordered.root" ]];then
    wget -o - -O Code/Data/ordered.root "https://www.dropbox.com/scl/fi/2ipkrkxud5k9j7hglh64c/ordered.root?rlkey=8qrtubtqdedszypb6wd6o7d77&dl=0"
fi

if ! [[ -e "Code/Data/muons_100k_100br_new.root" ]];then
    wget -o - -O Code/Data/muons_100k_100br_new.root "https://www.dropbox.com/scl/fi/w7rgcgnv4pd9gusy70yt5/muons_100k_100br_new.root?rlkey=tmxrj8ckqjgcaidwdssph9w2e&st=zbch3kpb&dl=0"
fi

if ! [[ -e "Code/Data/muons_amuons_100k_100br.root" ]];then
    wget -o - -O Code/Data/muons_amuons_100k_100br.root "https://www.dropbox.com/scl/fi/ag0ti9s9nl4l8lcwqufyn/muons_amuons_100k_100br.root?rlkey=gt6trjhc7ufviq1fmxldzbypl&st=40pz9xpa&dl=0"
fi

if ! [[ -e "Code/Data/muons_amuons_50k_200br.root" ]];then
    wget -o - -O Code/Data/muons_amuons_50k_200br.root "https://www.dropbox.com/scl/fi/64shia8q6k8wqosl66tez/muons_amuons_50k_200br.root?rlkey=l9h6q3aedib3hz5mvliwk4jnr&st=q7jqw5qa&dl=0"
fi

if ! [[ -e "Code/Data/muons_amuons_50k_300br.root" ]];then
    wget -o - -O Code/Data/muons_amuons_50k_300br.root "https://www.dropbox.com/scl/fi/prh95qicdutraa2pj79gu/muons_amuons_50k_300br.root?rlkey=62u2waok2o1yjmr4mxcmwdw6u&st=lkvim2o5&dl=0"
fi

if ! [[ -e "Code/Data/muons_amuons_50k_400br.root" ]];then
    wget -o - -O Code/Data/muons_amuons_50k_400br.root "https://www.dropbox.com/scl/fi/xnpf34o5ayw89c3f4xx9k/muons_amuons_50k_400br.root?rlkey=bni16mqwv55rv3d6rff7es0ej&st=mx8kgwf0&dl=0"
fi

if ! [[ -e "Code/Data/double_tracks_muons_amuons_50k_100br.root" ]];then
    wget -o - -O Code/Data/double_tracks_muons_amuons_50k_100br.root "https://www.dropbox.com/scl/fi/nv2db0d6ib1censt9xwkb/double_tracks_muons_amuons_50k_100br.root?rlkey=zlhanfdupsn1os6w6h3q4660c&st=zhtp1xfq&dl=0"
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


if ! [[ -e "Code/Data/muons_100k_100br_new.root" ]]; then
    wget -o - -O Code/Data/muons_100k_100br_new.root "https://www.dropbox.com/scl/fi/w7rgcgnv4pd9gusy70yt5/muons_100k_100br_new.root?rlkey=tmxrj8ckqjgcaidwdssph9w2e&st=zbch3kpb&dl=0"
fi

echo "Execution terminated"