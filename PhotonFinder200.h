#ifndef PhotonF_h
#define PhotonF_h


// Stdlib header file for input and output.
#include <iostream>

// FastJet3 library.
#include "Pythia8/Pythia.h"
//#include "Pythia8Plugins/FastJet3.h"
//#include "/home/sleethmr/jetbox/myX-SCAPE/X-SCAPE/external_pa//ckages/fastjet-install/include/fastjet/PseudoJet.hh"
//#include "/home/sleethmr/jetbox/myX-SCAPE/X-SCAPE/external_packages/fastjet-install/include/fastjet/ClusterSequence.hh"

#include <vector>
#include <fstream>
#include <cstdio>
#include "iomanip"
#include "TMath.h"

// ROOT, for histogramming.
#include "TH1.h"
#include "TH2.h"

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraphErrors.h"

// ROOT, for saving file.
#include "TFile.h"

// Histogram coloring.
#include "TStyle.h"
#include "string"


using namespace std;

Pythia8::Pythia pythia ("",false);

int StartTime = time(NULL);                                               
int sqrtSNN = 200;
string Input;
string OutputDirName = "ResultsRoot";

int NpTHardBin = 0;
std::vector<int> pTHatMin;
std::vector<int> pTHatMax;

ostringstream Outfilename, Particlelistfile;
//Hardcrosssectionfile;
//char pTBinString[100];

//std::vector<TH1F*> HistFragPtHat;
//TH1F *HistFragTotal;

//std::vector<double> HardCrossSection;
//std::vector<double> HardCrossSectionError;

std::vector<int> nEvent;

int nEvent_k = 0;

void InputSettings(int argc, char* argv[]);
void GetInput();
void InitializeForThisPtHardBin(int k);
//void LoadHardCrossSectionForThisPtHardBin(int k);
void LoadParticlesFromList(int k);
//void PtHatSumHist();
void OutCalculationTime();

//---------------------------------------------------------------------------
int NpTHardBin200 = 23;
std::vector<int> pTHatMin200 = { 1, 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90};
std::vector<int> pTHatMax200 = { 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100};
//---------------------------------------------------------------------------
int NpTHardBin2760 = 54;
std::vector<int> pTHatMin2760 = { 1, 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900, 1000 };
std::vector<int> pTHatMax2760 = { 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900, 1000, 1380 };
//---------------------------------------------------------------------------
int NpTHardBin7000 = 69;
std::vector<int>pTHatMin7000 = { 1, 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2600, 2800, 3000 };
std::vector<int>pTHatMax7000 = { 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2600, 2800, 3000, 3500 };
//---------------------------------------------------------------------------


int NptBinATLAS2760 = 10;
std::vector<double>ptBinATLAS2760
= {1, 1.6, 2.5, 4.0, 6.3, 10, 16, 25, 40, 63, 100};

int NzBinATLAS2760 = 10;
std::vector<double>zBinATLAS2760
= {0.01, 0.016, 0.025, 0.04, 0.063, 0.1, 0.16, 0.25, 0.4, 0.63, 1};


#endif
