# X-SCAPE-Photons-Analysis

## Follow directions on how to install X-SCAPE in singularity on an HPC grid, the one used in my analysis was Vanderbilt's ACCRE.
(Just in case, I will put down my steps here):
Operating system: LINUX
ssh into your HPC grid account
```
mkdir jetbox
cd jetbox
singularity build --jet.sif docker://jetscape/base:stable
singularity shell jet.sif
```

## X-SCAPE Installation
(in singularity shell)
```
mkdir myX-SCAPE
cd myX-SCAPE
git clone https://github.com/JETSCAPE/X-SCAPE.git
```

### External Packages
(cmake components will be in build directory after running scripts to install packages)
```
cd ../external_packages
./get_3dglauber.sh
./get_music.sh
./get_iSS.sh
./get_freestream-milne.sh
```
Download PYTHIA (in external_packages directory)
```
wget phttps://pythia.org/download/pythia83/pythia8309.tgz
tar zxf pythia8309.tgz
mv pythia8309.tgz pythia8309
cd pythia8309
./configure
make -j4 install
```
Download Eigen (in external_packages directory)
```
wget -O eigen.tar.bz2 https://gitlab.com/libeigen/eigen/-/archive/3.3.9/eigen-3.3.9.tar.bz2
tar xjf eigen.tar.bz2
mv eigen-3.3.9 eigen
mv eigen.tar.bz2 eigen
```

### Define Environment variables (in external_packages directory)
Create a new file 
```
emacs prepare.sh
```
Write this in the file:
```
export EIGEN_INSTALL_DIR=`readlink -f ./eigen`
export EIGEN3_ROOT=`readlink -f ./eigen`
export GSL=$(gsl-config --prefix)
export GSL_HOME=$(gsl-config --prefix)
export GSL_ROOT_DIR=$(gsl-config --prefix)
export XSCAPE_DIR=`readlink -f .`
export SMASH_DIR=/home/myX-SCAPE/X-SCAPE/external_packages/smash/smash_code
export PYTHIAINSTALLDIR=`readlink -f .`
export PYTHIA8DIR=${PYTHIAINSTALLDIR}/pythia8309
export PYTHIA8_ROOT_DIR=${PYTHIAINSTALLDIR}/pythia8309
export CC=gcc
export CXX=g++
export OpenMP_CXX=g++
```
Save and exit file, then source it
```
source prepare.sh
```
PYTHIA and Eigen are pre-requisites to SMASH, so now we can install SMASH
```
./get_smash.sh
```
And then source the environment file again
```
source prepare.sh
```
(The directory may need to be changed for the SMASH_DIR if it gives trouble in the next step)
### Compiling X-SCAPE using the external packages
```
cd ../build
cmake .. -DUSE_3DGlauber=ON -DUSE_MUSIC=ON -DUSE_FREESTREAM=ON -DUSE_ISS=ON -DUSE_SMASH=ON
```
Wait for cmake to compile, then make and install:
```
make -j4 install
```
### Running tests
MUSIC Test:
(in build directory)
```
mpirun -np 1 ./MUSICTest
```
Let it run, if there are root errors then try adding to environment file these lines:
```
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
```
There may also be an error that T_freeze = 0.3 is not physical, but I ignore it

SMASH Test:
```
./SMASHTest
```
Same here with the T_freeze error

### Optional Reader Test
```
./readerTest
```
(will save output under three files:
my_test.gv, my_test.gml, my_test.graphml

scp to local machine, then run command to convert to PDF
(from local machine terminal)
```
scp <username>@login.accre.vanderbilt.edu:/home/<username>/jetbox/myX-SCAPE/X-SCAPE/build/my_test.gv ~/
```
(Make sure you have the converter software "dot")
```
dot my_test.gv -Tdf -o outputPDF.pdf
```
When you open outputPDF.pdf, it should show parton splitting 

## Running JETSCAPE p+p results publication file
```
./runJetscape ../config/publications_config/arXiv_1910.05481/jetscape_user_PP_1910.05481.xml
```
Resulting files will be test_out.dat, test_out_final_state_hadrons.dat, and test_out_final_state_partons.dat


