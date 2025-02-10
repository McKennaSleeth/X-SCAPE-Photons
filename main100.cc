// Copyright (C) 2022 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: analysis; root;

// This is a simple test program.
// It studies the charged multiplicity distribution at the LHC.
// Modified by Rene Brun, Axel Naumann and Bernhard Meirose
// to use ROOT for histogramming.

// Stdlib header file for input and output.
#include <iostream>
#include <vector>
#include <stdio.h>
#include <gsl/gsl_statistics.h>
#include <cmath>
// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"

using namespace Pythia8;

//int main(int argc, char* argv[]) {
int main (){

  // Create Pythia instance and set it up to generate hard QCD processes above pTHat = 20 GeV for pp collisions at 14 TeV.
  Pythia pythia;
  //Settings& settings = pythia.settings;
  const Info& info = pythia.info;
  const long int nEvent = 100000;
  const int pTHatbins[] = {5,15};
  const int npTHatbins = sizeof(pTHatbins)/sizeof(int);
  const double pTlist[] = {5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0};
  const int npT = sizeof(pTlist)/sizeof(double) - 1;
  std::vector<std::vector<double>> c(npT, std::vector<double>(nEvent, 0.0));

  //pythia.readString("PhaseSpace:pTHatMin = 0.");

  // Initialize for RHIC energies at 200 GeV, LHC energies at 2.76 TeV

  pythia.readString("Beams:eCM = 200.");
  pythia.readString("PromptPhoton:all = on");
  //pythia.readString("HadronLevel:Decay = off");
  pythia.readString("111:mayDecay = off");
  pythia.readString("HardQCD:all = on");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:Seed = 0");

  // Begin pTHatBin loop
  for (int ipTHat = 0; ipTHat < npTHatbins - 1; ipTHat++){
    const double pTHatMin = pTHatbins[ipTHat];
    const double pTHatMax = pTHatbins[ipTHat + 1];
    std::string min = "PhaseSpace:pTHatMin=";
    std::string max = "PhaseSpace:pTHatMax=";
    min += std::to_string(pTHatMin);
    max += std::to_string(pTHatMax);
    std::cout << min.c_str() << "\n";
    std::cout << max.c_str() << "\n";
    std::string pTHatMinString = "PhaseSpace:pTHatMin="+std::to_string(pTHatMin);
    std::string pTHatMaxString = "PhaseSpace:pTHatMax="+std::to_string(pTHatMax);
    pythia.readString(pTHatMinString.c_str());
    pythia.readString(pTHatMaxString.c_str());
    pythia.init();
    
    // Begin event loop. Generate event; skip if generation aborted.
    for (int iEvent = 0; iEvent < nEvent; iEvent++) {
      if (!pythia.next()) continue;

      // Begin event selection loop. Find number of all final charged particles.
      for (int i = 0; i < pythia.event.size(); ++i){

        //if (fabs(pythia.event[i].id())==211 && pythia.event[i].isFinal()){
	// Begin pion selection loop (neutral pions PID = 111, charged pions = |211|)
	if (pythia.event[i].id()==111){
	  // Begin eta selection loop (Data for pions has |eta| < 0.35)
	  if (fabs(pythia.event[i].eta())<0.35){

	    // Begin pTbin loop
	    for (int a = 0; a < npT; ++a){
	      
	      const double pTmin = pTlist[a];
	      const double pTmax = pTlist[a+1];
	      // Begin pT selection loop
	      if (pTmin < pythia.event[i].eT() && pythia.event[i].eT() < pTmax){
		c[a][iEvent] += 1*info.sigmaGen();
	      } // End of pion momentum selection loop
	      
	    } // End of pTbin for loop
	    
          } // End of eta selection loop
	  
        } // End of pion selection loop
	
      } // End event selection loop
      pythia.stat();
    } // End event number loop
    
  } // End pTHatbin loop
  
  
  
  double mean, stan_err_mean;
  // Begin second pTbin loop
  for (int a = 0; a < npT; ++a){
    mean = gsl_stats_mean(c[a].data(), 1, nEvent);
    stan_err_mean = (sqrt(gsl_stats_variance(c[a].data(), 1, nEvent))/sqrt(nEvent));
    const double pTmin = pTlist[a];
    const double pTmax = pTlist[a+1];
    double deltaE = pTmax - pTmin;
    double deltaeta = 0.7;
    double avgE = (pTmax + pTmin)/2;
    double cross_section = (1/(2*3.14*avgE)) *((mean)/(deltaE*deltaeta));
    
    double cross_sec_err = (cross_section/mean)*stan_err_mean;
    
    //Output
    
    cout << (pTmin+pTmax)/2 << " " << cross_section << " " << cross_sec_err << "\n";                                               
  } // End second pTbin loop

  // Done.
  return 0;
} // End main function
 
