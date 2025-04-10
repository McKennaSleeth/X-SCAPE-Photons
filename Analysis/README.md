# File Descriptions and Make Commands:

## Direct Photons:
The file under the name direct_photon_cross_sec.cpp calculates the photon cross section by taking final_state_hadrons.dat files from different pTHat runs of jetscape_pp_direct_photons.xml

## Fragmentation Photons:
The file under the name of frag_photon_cross_sec.cpp calculates the photon cross section by taking final_state_hadrons.dat files from different pTHat runs of jetscape_pp_direct_photons.xml, but with PromptPhoton:all = off

The file under the name of photonfrag_func_deltaR.cpp calculates the momentum fraction z and deltaR of photons found in jets by taking final_state_hadrons.dat file from one pTHat run of jetscape_pp_frag_photons.xml

## Charged and Neutral Pions:
The files under the name of pioncross_sec.cpp and pionpluscross_sec.cpp calculate the neutral pion and charged pion cross sections (respectively) by taking final_state_hadrons.dat files from different pTHat runs of jetscape_pp_pions.xml (not yet added to GitHub)

The file under the name of pionfrag_func_deltaR.cpp calculates the momentum fraction z and deltaR of pions found in jets by taking final_state_hadrons.dat file from one pTHat run of jetscape_pp_frag_pions.xml
