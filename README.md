# dualtoy
dual calorimeter toy for HONR class




cmake -DGeant4_DIR=/cvmfs/geant4.cern.ch/geant4/10.5/x86_64-slc6-gcc63-opt/lib64/GEANT4-10.5.0

./dualtoy -c template.cfg -u Xm  (-m run.mac if batch mode) -o filename

make histograms with toyplot.cc