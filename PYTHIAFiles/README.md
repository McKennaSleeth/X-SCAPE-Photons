# File key:

## Isolated Photons:
main96.cc -- isolated photon cross sections (selection method = status codes to disable hadron decay and find specific initial state showers)

Make command:
```
g++ -o main96 main96.cc -lgsl -lgslcblas -lm -I../include -O2 -std=c+++11 -pedantic -W -Wall -Wshadow -fPIC -pthread -L../lib -WI,-rpath,../lib -lpythia8
```

## Fragmentation Photons:
main98.cc -- fragmentation photon cross section (selection method = status codes to disable hadron decay and find specific initial state showers)

Make command:
```
g++ -o main98 main98.cc -lgsl -lgslcblas -lm -I../include -O2 -std=c+++11 -pedantic -W -Wall -Wshadow -fPIC -pthread -L../lib -WI,-rpath,../lib -lpythia8 
```

main997.cc -- fragmentation photon cross section (selection method = photon found in jet with PromptPhoton:all = off and no hadron decay)

Make command:
```
g++ -o main997 main997.cc -lgsl -lgslcblas -lm -I../include -O2 -std=c+++11 -pedantic -W -Wall -Wshadow -fPIC -pthread -L../lib -WI,-rpath,../lib -lpythia8 -I/usr/local/include/fastjet/include -Wl,-rpath,usr/local/lib -L/usr/local/lib -lfastjettools -lfastjet -lm
```

## Direct Photons:
main99.cc -- direct photon cross section (selection method = PID and no hadron decay)

Make command:
```
g++ -o main99 main99.cc -lgsl -lgslcblas -lm -I../include -O2 -std=c+++11 -pedantic -W -Wall -Wshadow -fPIC -pthread -L../lib -WI,-rpath,../lib -lpythia8 
```

## Charged and Neutral Pions:
main100.cc -- pion cross section (selection method = 1. PID = |211| for pion plus and pion minus and 2. PID = 111 & no hadron decay for neutral pion)

Make command:
```
g++ -o main100 main100.cc -lgsl -lgslcblas -lm -I../include -O2 -std=c+++11 -pedantic -W -Wall -Wshadow -fPIC -pthread -L../lib -WI,-rpath,../lib -lpythia8 
```

