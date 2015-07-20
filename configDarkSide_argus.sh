#!/bin/sh

if [ $DS_LNGS_CLUSTER ] ; then

 export G4DS=$PWD
 export DSDATA=$PWD/data/physics
 export PATH=$PATH:$PWD/tools
 export G4WORKDIR=$PWD/.g4ds

 unset G4VIS_USE_DAWN           
 unset G4VIS_USE_OPENGLQT
 unset G4VIS_USE_OPENGLX  						   
 unset G4VIS_USE_RAYTRACERX
 unset G4VIS_USE_VRML
 unset G4VIS_USE_OIX
 unset G4VIS_USE_OPENGLXM
 unset G4UI_USE_QT

 source /usr/local/share/geant4/geant4make/geant4make.sh

 export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$G4WORKDIR/tmp/Darwin-g++/g4ds

else

 unset ROOTSYS
 export G4DS=${PWD}
 export ROOTSYS=/Applications/Shared/root_5.34.00

 export DSDATA=$PWD/data/physics

 export PATH=${ROOTSYS}/bin:${PWD}/tools:${PATH}:/sw/lib
 export LD_LIBRARY_PATH=${ROOTSYS}/lib:${ROOTSYS}/lib/root:${LD_LIBRARY_PATH}
 export DYLD_LIBRARY_PATH=${ROOTSYS}/lib:${ROOTSYS}/lib/root:${DYLD_LIBRARY_PATH}

 source /Applications/Shared/geant4.9.5/share/Geant4-9.5.0/geant4make/geant4make.sh

 unset G4VIS_USE_DAWN           
 unset G4VIS_USE_OPENGLQT
 unset G4VIS_USE_OPENGLX  						   
 unset G4VIS_USE_RAYTRACERX
 
 unset G4VIS_USE_OIX
 unset G4VIS_USE_OPENGLXM
 unset G4UI_USE_QT
						
 export G4DATA=/Applications/Shared/geant4.9.5/share/Geant4-9.5.0/data/
 export G4LEDATA=$G4DATA/G4EMLOW6.23
 export G4LEVELGAMMADATA=$G4DATA/PhotonEvaporation2.2
 export G4RADIOACTIVEDATA=$G4DATA/RadioactiveDecay3.4
 export NeutronHPCrossSections=$G4DATA/G4NDL4.0
 export G4NEUTRONHPDATA=$G4DATA/G4NDL4.0
 export G4NEUTRONXSDATA=$G4DATA/G4NEUTRONXS1.1
 export G4REALSURFACEDATA=$G4DATA/RealSurface1.0
 export G4PIIDATA=$G4DATA/G4PII1.3
 export G4WORKDIR=$PWD/.g4ds

 #G4NEUTRONHPDATA=/work/GEANT4/geant4.9.6/share/Geant4-9.6.0/data/G4NDL4.2
 #G4LEDATA=/work/GEANT4/geant4.9.6/share/Geant4-9.6.0/data/G4EMLOW6.32
 #G4LEVELGAMMADATA=/work/GEANT4/geant4.9.6/share/Geant4-9.6.0/data/PhotonEvaporation2.3
 #G4RADIOACTIVEDATA=/work/GEANT4/geant4.9.6/share/Geant4-9.6.0/data/RadioactiveDecay3.6
 #G4NEUTRONXSDATA=/work/GEANT4/geant4.9.6/share/Geant4-9.6.0/data/G4NEUTRONXS1.2
 #G4PIIDATA=/work/GEANT4/geant4.9.6/share/Geant4-9.6.0/data/G4PII1.3
 #G4SAIDXSDATA=/work/GEANT4/geant4.9.6/share/Geant4-9.6.0/data/G4SAIDDATA1.1


 export DYLD_LIBRARY_PATH=${G4WORKDIR}/tmp/Darwin-g++/g4ds/:$XERCESCROOT/lib:${DYLD_LIBRARY_PATH}
 export LD_LIBRARY_PATH=${G4WORKDIR}/tmp/Darwin-g++/g4ds/:$XERCESCROOT/lib:${LD_LIBRARY_PATH}
fi
