# dualtoy
dual calorimeter toy for HONR class


first clone with using
   git clone git@github.com:saraheno/dualtoy.git

cd into dualtoy and do
   sh ./g4env.sh

then compile it with
  cmake -DGeant4_DIR=/cvmfs/geant4.cern.ch/geant4/10.5/x86_64-slc6-gcc63-opt/lib64/GEANT4-10.5.0
  make

you can run it two ways.  to run interactively, do 
  ./dualtoy -c template.cfg -u Xm

while to run without the display do 
  ./dualtoy -c template.cfg -m run.mac -o filename
where filename is the name of the file you want GEANT to put the output into

Look inside of template.cfg.  You will see many options you can change about the detector geometry
You can also see other files like run*.mac.  these allow running different particles and energies
Compare them.

Once you have run a bunch of particles at different energies, you need to make histograms
from the geant output.  This is done with toyplot.cc.   go into root and type
   .x toyplot.cc
Do this once for each GEANT root file you created.  You should edit the two lines just after
"void toyplot() {" in this file to change the input and output file names.

Now, you want to take these histograms and use them to make resolution curve.
Use res.C to do this in root.


