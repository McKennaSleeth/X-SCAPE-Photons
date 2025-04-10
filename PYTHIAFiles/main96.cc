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

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"

using namespace Pythia8;

//int main(int argc, char* argv[]) {
int main (){
  
  // Create Pythia instance and set it up to generate hard QCD processes
  // above pTHat = 20 GeV for pp collisions at 14 TeV.
  Pythia pythia;
  //Settings& settings = pythia.settings;
  const Info& info = pythia.info;
  const long int nEvent = 1000000;
  //const int pTHatbins[3] = {5,15,25}; //2 pT bins
  const int pTHatbins[] = {3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20}; // pTHat bins
  const int npTHatbins = sizeof(pTHatbins)/sizeof(int);
  //cout << "npTHatbins = " << npTHatbins;
  //const int pTlist[8] = {2,5,7,10,12,15,17,20};
  const double pTlist[] = {3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10,12,14,16};
  const int npT = sizeof(pTlist)/sizeof(double);
  //cout << "npT = " << npT;
 
  //double c[npT][nEvent] = {0.0};
  std::vector<std::vector<double>> c(npT, std::vector<double>(nEvent, 0.0));
     
  // Initialize for RHIC energies at 200 GeV, LHC energies at 2.76 TeV
  pythia.readString("Beams:eCM = 200.");
  pythia.readString("PromptPhoton:all = on");
  pythia.readString("HardQCD:all = on");
  
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
      
      // Find number of all final charged particles.      
      for (int i = 0; i < pythia.event.size(); ++i){
	
	if (pythia.event[i].id()==22 && pythia.event[i].isFinal() ){
	  
	  if (fabs(pythia.event[i].eta())<0.35 && pythia.event[i].status()<90){ //pythia event status less than 90 to exclude decay photons

	    if(fabs(pythia.event[pythia.event[i].iTopCopy()].status())<40){ //pythia event status less than 40 to select for particles produced by initial state showers
		for (int a = 0; a < npT - 1; a++){
		  const double pTmin = pTlist[a];
		  const double pTmax = pTlist[a+1];
		  if (pTmin < pythia.event[i].eT() && pythia.event[i].eT() < pTmax){
		    
		    c[a][iEvent] += 1*info.sigmaGen();
		    
		  } // End of photon momentum selection loop
		  
		} // End of pTbin for loop
		
	    } // End of direct prompt photon selection loop
	    
	  } // End of eta selection and decay photon exclusion loop
	  
	} // End of photon selection loop
	
      } // End event selection loop
      
      //c.push_back(nPhoton);
      pythia.stat();
      
    } // End event number loop
    
  } // End pTHatbin loop
  
    
    //double* N = &c[0];
  double mean, stan_err_mean;
  ofstream RHIC_iso_4bins_hard_cross_sec;
  RHIC_iso_4bins_hard_cross_sec.open("RHIC_iso_4bins_hard_cross_sec.txt");
  //RHIC_iso_4bins_hard_cross_sec << "# Photon Cross Sections vs pT for RHIC energies (200 GeV) " << "\n";
  //RHIC_iso_4bins_hard_cross_sec << "# 100,000 events, pTHatbins = 5,15,25 GeV/c" << "\n";
  //RHIC_iso_4bins_hard_cross_sec << "# pTmin    pTmax    cs    sigmaGen     stan err mean      cs err" << "\n";
  
  for (int a = 0; a < npT - 1; a++){
    //mean = gsl_stats_mean(c[a], 1, nEvent);
    //stan_err_mean = (sqrt(gsl_stats_variance(c[a], 1, nEvent))/sqrt(nEvent));
    mean = gsl_stats_mean(c[a].data(), 1, nEvent);
    stan_err_mean = (sqrt(gsl_stats_variance(c[a].data(), 1, nEvent))/sqrt(nEvent));
    
    
    const double pTmin = pTlist[a];
    const double pTmax = pTlist[a+1];
    double deltaE = pTmax - pTmin;
    double deltaeta = 2;
    double avgE = (pTmax + pTmin)/2;
    double cross_section = (1/(2*3.14*avgE)) *((mean)/(deltaE*deltaeta)); //change this later
    double cross_sec_err = (cross_section/mean)*stan_err_mean;
    //RHIC_iso_4bins_hard_cross_sec <<  "pTmin=" << pTmin << " " << "pTmax=" << pTmax << "\n" << "cross_section =" << cross_section << " " << "sigmaGen=" << info.sigmaGen() << " " << "stan_err_mean=" << stan_err_mean << " " << "cross sec err=" << cross_sec_err << "\n" ;
    RHIC_iso_4bins_hard_cross_sec << (pTmin+pTmax)/2 << " " << cross_section << " " << cross_sec_err << "\n";
    
    //Output
    //cout << "pTmin=" << pTmin << " " << "pTmax=" << pTmax << "\n" << "cross_section =" << cross_section << " " << "sigmaGen=" << info.sigmaGen() << "\n";
    cout << (pTmin+pTmax)/2 << " " << cross_section << " " << cross_sec_err << "\n";
    //std::cout << "mean=" << mean << "\n";
    //std::cout << "stan_err_mean=" << stan_err_mean << "\n";
    //std::cout << "cross sec err=" << cross_sec_err << "\n";
  }
  
  RHIC_iso_4bins_hard_cross_sec.close();
  // Done.
  return 0;
  
}
  
