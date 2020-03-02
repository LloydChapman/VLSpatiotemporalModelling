#!/bin/bash
cd /usr/local/MATLAB/R2018a/bin
./matlab -nosplash -nodisplay -nodesktop -r "RunMCMC(1); exit" &
./matlab -nosplash -nodisplay -nodesktop -r "RunMCMC(2); exit" &
./matlab -nosplash -nodisplay -nodesktop -r "RunMCMC(3); exit" &
./matlab -nosplash -nodisplay -nodesktop -r "RunMCMC(4); exit" &
./matlab -nosplash -nodisplay -nodesktop -r "RunMCMC(5); exit" &
./matlab -nosplash -nodisplay -nodesktop -r "RunMCMC(6); exit"
