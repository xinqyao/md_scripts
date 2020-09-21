#!/bin/bash

## First, you need to prepare a set of PDB files, in a folder named e.g. 'pdbs', as the input.
## These PDB files should be taken from your simulation trajectory. For example, using cpptraj you can extract 1% of the trajectory by skipping every 100 frames and save the frames as PDB.
##
## Then, you need a configuration file. In this folder, there is an example 'config.txt'. Check and make all necessary changes before proceeding.
## The variables that you most likely will change are:
##    last_frame
##    starting_point_atom
##
##
java -Xmx1200m -cp /software/caver_3.0/caver/lib -jar /software/caver_3.0/caver/caver.jar -conf config.txt -home /software/caver_3.0/caver -pdb pdbs -out results

## Results will be generated under the 'results' folder.
##
## You can also check out the 'plot_bottleneck_dynamics.r' for further analysis.

