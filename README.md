# CMS-SpikingNeuralNetwork

Template commands used to compile:
Laptop:
`c++ $(root-config --cflags) -I /home/ema/root/include -I ./Class -o SNNT13.out SNNT13.C ./Class/SNN.cc $(root-config --libs)`
Mac-Os:
`c++ -g -I. `root-config --cflags` -o SNNT13.out SNNT13.C Class/SNN.cc `root-config --libs``

Home:
`g++ $(root-config --cflags) -I /home/ema/anaconda3/envs/ROOT/bin -I ./Class -o SNNT13.out SNNT13.C ./Class/SNN.cc $(root-config --libs)`

If you are using the program for the first time, please launch the first_setup.sh script.
On Linux systems you may have to use:
`chmod +x first_setup.sh`
`./first_setup.sh`

