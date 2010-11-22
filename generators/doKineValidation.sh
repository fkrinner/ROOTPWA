#!/bin/bash

export WORKDIR=/afs/e18/compass/analysis/sneubert/
export FITDIR=$WORKDIR/PWAFITS/LOWT/NEWACC/fit17GenSelected
#export FITDIR=$WORKDIR/PWAFITS/GENETICS/ltRUN21/gen39/set42
export FITFILE=$FITDIR/fit17.root
export PLOTFILE=${FITFILE/.root/.plots.root}
export BOOKY=${FITFILE/.root/.booky.pdf}
export DATADIR=$WORKDIR/5PiLTData3/

echo WORKDIR=$WORKDIR
echo FITDIR=$FITDIR
echo FITFILE=$FITFILE
echo PLOTFILE=$PLOTFILE
echo BOOKY=$BOOKY
echo DATADIR=$DATADIR


cd $DATADIR

rm $PLOTFILE;

for i in *; do
    echo "MassBin: $i";
    cd $i;
    # convert events to root tree if not already done
    test -s $i.root || cat $i.evt | evt2tree $i.root;
    # run evtweight on accepted events:
    cd ACCAMPS
    # cd PSPAMPS
    WEIGHTEDFILE=${FITFILE/.root/.kineval.$i.root}
    test -s $WEIGHTEDFILE || evtweight -e ../$i.acc.evt -o $WEIGHTEDFILE  -w $FITFILE -i accnorm.int -m $i  
    # produce nice plots
    root -b -q "$ROOTPWA/generators/doPlotWEvts.C(\"../$i.root\",\"$WEIGHTEDFILE\",\"$PLOTFILE\",\"$i\")"

    cd $DATADIR;
done;

echo "CREATING BOOKY $BOOKY ..."
# collect all plots into a booky:
cd $FITDIR
for i in *plots*.ps; do ps2pdf $i; rm -f $i; done;
pdftk *plots*.pdf cat output $BOOKY
for i in *plots*.pdf; do rm -f $i; done;


