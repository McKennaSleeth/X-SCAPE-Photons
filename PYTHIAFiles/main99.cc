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
#include "TFile.h"
#include "TH2F.h"

using namespace Pythia8;

//int main(int argc, char* argv[]) {
int main (){
 
  // Create Pythia instance and set it up to generate hard QCD processes
  // above pTHat = 20 GeV for pp collisions at 14 TeV.
  Pythia pythia;
  //Settings& settings = pythia.settings;
  const Info& info = pythia.info;
  const long int nEvent = 100000;
  //const int pTHatbins[3] = {5,15,25};
  //const int pTHatbins[] = {5,20};
  const int pTHatbins[] = {3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
  //const int pTHatbins[] = {5,10,15,20,25};
  const int npTHatbins = sizeof(pTHatbins)/sizeof(int);
  //const int pTlist[8] = {2,5,7,10,12,15,17,20};
  const double pTlist[] = {3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 12.0, 14.0, 16.0};
  //const int pTlist[] = {5,10,15,20};
  const int npT = sizeof(pTlist)/sizeof(double) - 1;
  //double c[npT][nEvent] = {0};
  std::vector<std::vector<double>> c(npT, std::vector<double>(nEvent, 0.0));
  //double pionCountsPerBin[npT] = {0.0};

  TH2F* hPhotonStatus_vs_pT = new TH2F("hPhotonStatus_vs_pT", "Photon Status vs pT;Status;pT (GeV)", 100, 0, 100, 100, 0, 20);
  //pythia.readString("HardQCD:all = on");
  // pythia.readString("HardQCD:gg2gg=on"); // This interaction doesn't work
  //pythia.readString("HardQCD:gg2qqbar = on"); // This interaction works
  // pythia.readString("HardQCD:qg2qg = on"); // This interaction doesn't work
  // pythia.readString("HardQCD:qq2qq = on"); // This interaction doesn't work
  //pythia.readString("HardQCD:qqbar2gg = on"); // This interaction works
  //pythia.readString("HardQCD:qqbar2qqbarNew = on"); // This interaction works
  // pythia.readString("HardQCD:nQuarkNew = 3"); // This interaction works
  //pythia.readString("PhaseSpace:pTHatMin = 0.");
  
  // Initialize for RHIC energies at 200 GeV, LHC energies at 2.76 TeV
  
  pythia.readString("Beams:eCM = 200.");
  pythia.readString("PromptPhoton:all = on");
  pythia.readString("HadronLevel:Decay = off");
  /*
  pythia.readString("111:mayDecay = off");
  pythia.readString("221:mayDecay = off");
  pythia.readString("331:mayDecay = off");
  pythia.readString("223:mayDecay = off");
  pythia.readString("113:mayDecay = off");
  pythia.readString("213:mayDecay = off");
  //pythia.readString("-213:mayDecay = off");
  pythia.readString("333:mayDecay = off");
  pythia.readString("2114:mayDecay = off");
  pythia.readString("2214:mayDecay = off");
  pythia.readString("3212:mayDecay = off");
  pythia.readString("310:mayDecay = off");
  pythia.readString("130:mayDecay = off");
  pythia.readString("3122:mayDecay = off");
  */
  pythia.readString("HardQCD:all = on");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:Seed = 0");
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
    
 
    //std::vector<double> c;
    
  // Begin event loop. Generate event; skip if generation aborted.
    for (int iEvent = 0; iEvent < nEvent; iEvent++) {
      if (!pythia.next()) continue;
      
      // Find number of all final charged particles.      
      for (int i = 0; i < pythia.event.size(); ++i){
	
	//if (fabs(pythia.event[i].id())==211 && pythia.event[i].isFinal()){
	if (pythia.event[i].id()==22 && pythia.event[i].isFinal()){
	  if (fabs(pythia.event[i].eta())<0.35){
	    
	      for (int a = 0; a < npT; ++a){

		const double pTmin = pTlist[a];
		const double pTmax = pTlist[a+1];

		if (pTmin < pythia.event[i].eT() && pythia.event[i].eT() < pTmax){

		  //pionCountsPerBin[a] += 1*info.sigmaGen();
		  c[a][iEvent] += 1*info.sigmaGen();
		  hPhotonStatus_vs_pT->Fill(pythia.event[i].status(), pythia.event[i].eT());
		  //c[a][iEvent] += 1;
		} // End of pion momentum selection loop
	      
	      } // End of pTbin for loop

	  } // End of eta selection loop

	} // End of pion selection loop
	
      } // End event selection loop
      
      //c.push_back(nPhoton);
      pythia.stat();
      
    } // End event number loop

  } // End pTHatbin loop
  
  
  TFile* outFile = new TFile("PhotonStatus_vs_pT.root", "RECREATE");
  hPhotonStatus_vs_pT->Write();
  outFile->Close();
  
  double mean, stan_err_mean;
  ofstream DirectPhotons_pp_pTHat3_20_hadrondecayoff_cross_sec;
  DirectPhotons_pp_pTHat3_20_hadrondecayoff_cross_sec.open("DirectPhotons_pp_pTHat3_20_hadrondecayoff_cross_sec.txt");
  
  //PionPlus_pp_pTHat5_10__1Mil_final_counts << "# avgpT    cs    cs err" << "\n";
  //double cross_section[npT] = {0.0};
  //double cross_sec_err[npT] = {0.0};
  //double mean[npT] = {0.0};
  //double stan_err_mean[npT] = {0.0};
  for (int a = 0; a < npT; ++a){
    mean = gsl_stats_mean(c[a].data(), 1, nEvent);
    stan_err_mean = (sqrt(gsl_stats_variance(c[a].data(), 1, nEvent))/sqrt(nEvent));
    //mean[a] = pionCountsPerBin[a]/nEvent;
    //double sumSquaredDifferences[npT] = {0.0};
    //double difference[npT] = {0.0};
    //double variance[npT] = {0.0};
    //double standardDeviation[npT] = {0.0};
    //double pionCountInEvent[npT] = {0.0};
    for (int iEvent = 0; iEvent < nEvent; iEvent++){
      for (int i = 0; i < pythia.event.size(); ++i){
	
	//if (fabs(pythia.event[i].id())==211 && pythia.event[i].isFinal()){
	if (pythia.event[i].id()==22 && pythia.event[i].isFinal()){
	  if (fabs(pythia.event[i].eta())<0.35){
	    
	    //pionCountInEvent[a]++;
	  }
	}
      
	//difference[a] = pionCountInEvent[a] - mean[a];
	//sumSquaredDifferences[a] += difference[a]*difference[a];
      }
    }
    //variance[a] = sumSquaredDifferences[a] / nEvent;
    //standardDeviation[a] = std::sqrt(variance[a]);
    //stan_err_mean[a] = standardDeviation[a] / std::sqrt(nEvent);
    
    //int totalPhotons = 0;
    //for (int iEvent = 0; iEvent < nEvent; iEvent++){
    //totalPhotons += c[a][iEvent];
    //}
    // cout << "pT bin " << a << " photon count: " << totalPhotons << std::endl; 
         
    const double pTmin = pTlist[a];
    const double pTmax = pTlist[a+1];
    double deltaE = pTmax - pTmin;
    double deltaeta = 0.7;
    double avgE = (pTmax + pTmin)/2;
    double cross_section = (1/(2*3.14*avgE)) *((mean)/(deltaE*deltaeta)); 
    //cross_section[a] = mean[a] / (deltaE * deltaeta * 2 * M_PI * avgE);
    double cross_sec_err = (cross_section/mean)*stan_err_mean;
    //cross_sec_err[a] = (cross_section[a] / mean[a]) * (stan_err_mean[a]);
    
    //PionPlus_pp_pTHat5_10__1Mil_final_counts <<  "pTmin=" << pTmin << " " << "pTmax=" << pTmax << "\n" << "cross_section =" << cross_section << " " << "sigmaGen=" << info.sigmaGen() << " " << "stan_err_mean=" << stan_err_mean << " " << "cross sec err=" << cross_sec_err << "\n";
    DirectPhotons_pp_pTHat3_20_hadrondecayoff_cross_sec << (pTmin+pTmax)/2 << " " << cross_section << " " << cross_sec_err << "\n";
    
    //Output
    //cout << "pTmin= " << pTmin << " pTmax=" << pTmax << "\n" << "cross_section = " << cross_section << " sigmaGen= " << info.sigmaGen() << "\n";
    cout << (pTmin+pTmax)/2 << " " << cross_section << " " << cross_sec_err << "\n";
    
    //cout << (pTmin+pTmax)/2 << " " << cross_section[a] << " " << cross_sec_err[a] << "\n";
    
    
    //std::cout << "mean=" << mean << "\n";
    //std::cout << "stan_err_mean=" << stan_err_mean << "\n";
    //std::cout << "cross sec err=" << cross_sec_err << "\n";
    
  }
  
  DirectPhotons_pp_pTHat3_20_hadrondecayoff_cross_sec.close();

  delete hPhotonStatus_vs_pT;
  delete outFile;
  // Done.
  return 0;
}
