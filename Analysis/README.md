# File Descriptions and Make Commands:

## Direct Photons:
The file under the name direct_photon_cross_sec.cpp calculates the photon cross section by taking final_state_hadrons.dat files from different pTHat runs of jetscape_pp_direct_photons.xml

Make command: 
```
g++ -o direct_photon_cross_sec direct_photon_cross_sec.cpp -lgsl -lgslcblas -lm -std=c++17
```

## Fragmentation Photons:
The file under the name of frag_photon_cross_sec.cpp calculates the photon cross section by taking final_state_hadrons.dat files from different pTHat runs of jetscape_pp_direct_photons.xml, but with PromptPhoton:all = off

Make command:
```
g++ -o frag_photon_cross_sec frag_photon_cross_sec.cpp -lgsl -lgslcblas -lm -std=c++17 -I/usr/local/include/fastjet/include -Wl,-rpath,usr/local/lib -L/usr/local/lib -lfastjettools -lfastjet -lm
```

The file under the name of photonfrag_func_deltaR.cpp calculates the momentum fraction z and deltaR of photons found in jets by taking final_state_hadrons.dat file from one pTHat run of jetscape_pp_frag_photons.xml

Make command:
```
g++ -o photonfrag_func_deltaR photonfrag_func_deltaR.cpp -lgsl -lgslcblas -lm -std=c++17 -I/usr/local/include/fastjet/include -Wl,-rpath,usr/local/lib -L/usr/local/lib -lfastjettools -lfastjet -lm
```


## Charged and Neutral Pions:
The files under the name of pioncross_sec.cpp and pionpluscross_sec.cpp calculate the neutral pion and charged pion cross sections (respectively) by taking final_state_hadrons.dat files from different pTHat runs of jetscape_pp_pions.xml (not yet added to GitHub)

Make commands:
```
g++ -o pioncross_sec pioncross_sec.cpp -lgsl -lgslcblas -lm -std=c++17
```
```
g++ -o pionpluscross_sec pionpluscross_sec.cpp -lgsl -lgslcblas -lm -std=c++17
```



The file under the name of pionfrag_func_deltaR.cpp calculates the momentum fraction z and deltaR of pions found in jets by taking final_state_hadrons.dat file from one pTHat run of jetscape_pp_frag_pions.xml

Make command:
```
g++ -o pionfrag_func_deltaR pionfrag_func_deltaR.cpp -lgsl -lgslcblas -lm -std=c++17 -I/usr/local/include/fastjet/include -Wl,-rpath,usr/local/lib -L/usr/local/lib -lfastjettools -lfastjet -lm
```
