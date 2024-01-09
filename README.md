# CMS-SpikingNeuralNetwork

Template commands used to compile:
Laptop:
c++ $(root-config --cflags) -I /home/ema/root/include -o SNNT13.out SNNT13.C ./Class/SNN.cc $(root-config --libs)

Home:
c++ $(root-config --cflags) -I /home/ema/anaconda3/envs/ROOT/bin -o SNNT13.out SNNT13.C ./Class/SNN.cc $(root-config --libs)

If you are using the program for the first time, please download these two files inside the Data folder:

wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1B6d3bxFBziWGIMmCp8uxxSF2As3k89m7' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1B6d3bxFBziWGIMmCp8uxxSF2As3k89m7" -O 100k_100br.root && rm -rf /tmp/cookies.txt

wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=17YxLI8qBCfBDKPHVCa-fAKVGeLocOf09' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=17YxLI8qBCfBDKPHVCa-fAKVGeLocOf09" -O ordered.root && rm -rf /tmp/cookies.txt