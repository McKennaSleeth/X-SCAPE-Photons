# X-SCAPE-Photons-Analysis

## Follow directions on how to install X-SCAPE in singularity on an HPC grid, the one used in my analysis was Vanderbilt's ACCRE.
(Just in case, I will put down my steps here):
Operating system: LINUX
ssh into your HPC grid account
mkdir jetbox
cd jetbox
singularity build --jet.sif docker://jetscape/base:stable
singularity shell jet.sif

## X-SCAPE Installation
(in singularity shell)
mkdir myX-SCAPE
cd myX-SCAPE
git clone https://github.com/JETSCAPE/X-SCAPE.git

### External Packages
(cmake components will be in build directory after running scripts to install packages)
