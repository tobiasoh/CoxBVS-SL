#!/bin/bash
  # Open multiple screens, send a command to each screen
  #screen -dmS test1
  #screen -S test1 -p 0 -X stuff "Rscript Rtest.R 1\n"
  #screen -dmS test2
  #screen -S test2 -p 0 -X stuff "Rscript Rtest.R 2\n"
  #screen -dmS test3
  #screen -S test3 -p 0 -X stuff "Rscript Rtest.R 3\n"


#example run: chmod +x run_real_data.sh
#bash run_real_data.sh path/to/dataset path/to/graph name 1000 4 50
#run_real_data.sh ./RealData/march14clin/dataset_
#which means: 20 runs using true graph, 1000 MCMC iterations with thinning parameter 4.

if [ "$1" = "" ]; then
        echo "Missing command line argument. Needs type of graph (true, empty, noise, partial_uniform, partial_non-uniform) and number of MCMC iterations (default 30k)."
        exit -1
fi

MCMC=$4

if [ "$4" = "" ]; then
        MCMC=30000
fi

THIN=$5
if [ "$5" = "" ]; then
        THIN=1
fi

A=$7
if [ "$7" = "" ]; then
        A=3
fi

B=$8
if [ "$8" = "" ]; then
        B=0.01
fi


IFS="_"
# Read the input string into two variables
read -r first_part second_part <<< "$3"

# Access the second part
#echo "Second part: $second_part"

CPUNUM=$6
if [ "$6" = "" ]; then
        CPUNUM=50
fi

for i in {1..20} 
do 
NUMSCREEN=$(($CPUNUM+i))
screen -dmS "$second_part$i"
screen -S "$second_part$i" -p 0 -X stuff "nice -n 18 taskset -c $NUMSCREEN  Rscript /data/tobiasoh/run_real_data.R $1 $2 $3$i $MCMC server -t$THIN -b$B -a$A\n"
done

