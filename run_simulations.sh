#!/bin/bash
  # Open multiple screens, send a command to each screen
  #screen -dmS test1
  #screen -S test1 -p 0 -X stuff "Rscript Rtest.R 1\n"
  #screen -dmS test2
  #screen -S test2 -p 0 -X stuff "Rscript Rtest.R 2\n"
  #screen -dmS test3
  #screen -S test3 -p 0 -X stuff "Rscript Rtest.R 3\n"


#example run: chmod +x run_simulations.sh
#bash run_simulations.sh true 1000 4
#which means: 20 runs using true graph, 1000 MCMC iterations with thinning parameter 4.

if [ "$1" = "" ]; then
        echo "Missing command line argument. Needs type of graph (true, empty, noise, partial_uniform, partial_non-uniform) and number of MCMC iterations (default 30k)."
        exit -1
fi

MCMC=$2

if [ "$2" = "" ]; then
        MCMC=30000
fi

THIN=$3
if [ "$3" = "" ]; then
        THIN=1
fi



  for i in {1..20} 
  do 
  screen -dmS "$1$i"
  screen -S "$1$i" -p 0 -X stuff "nice -n 19 taskset -c 20-88 Rscript /data/tobiasoh/Run_Simulation_Study.R $1$i $MCMC server t$THIN\n"

  done
  
  #for i in {1..20}
  #do
  #echo "$1$i"
  #echo "R $1$i $MCMC"
  #done