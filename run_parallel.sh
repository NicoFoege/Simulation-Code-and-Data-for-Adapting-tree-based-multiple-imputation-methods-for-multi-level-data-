#!/bin/bash

# Starte 60 Screen-Sessions und setze in jeder einen Wert für `ii`
for ii in {1..60}; do
  screen -dmS "RunCode_$ii" bash -c "export OMP_NUM_THREADS=4; Rscript RunCode.R $ii"
done