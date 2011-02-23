#!/bin/bash

export FITDIR=$PWA_DATA_DIR/fits
export FITFILE=$FITDIR/fit.root
export PLOTFILE=${FITFILE/.root/.plots.root}

echo FITDIR=$FITDIR
echo FITFILE=$FITFILE
echo PLOTFILE=$PLOTFILE

cd $FITDIR
root -l -b "$ROOTPWA/generators/rootlogon.C" "$ROOTPWA/generators/plotGlobalWeightedEvts_3pin.C+(\"fit.plots.root\", \"weightedMC.globalplots.root\")"

exit

