#include <vector>
#include <cmath>
#include <stdio.h>
#include <iostream>
#include <gsl/gsl_statistics.h>
#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

using namespace Pythia8;

double deltaR(const fastjet::PseudoJet& p1, const fastjet::PseudoJet& p2){
  double deta = p1.eta() - p2.eta();
  double dphi = std::abs(p1.phi() - p2.phi());
  if (dphi > M_PI) dphi = 2 * M_PI - dphi;
  return std::sqrt(deta * deta + dphi * dphi);
}

int main() {

  Pythia pythia;
  const Info& info = pythia.info;
  const long int nEvent = 1000000;
  const double pTlist[] = {3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 12.0, 14.0, 16.0};
  const int npT = sizeof(pTlist)/sizeof(double) - 1;
  
  const double pTHatbins[] = {3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
  const int npTHatbins = sizeof(pTHatbins)/sizeof(double) - 1;
  std::vector<std::vector<double>> c(npT, std::vector<double>(nEvent, 0.0));
  
  pythia.readString("Beams:eCM = 200.");
  //pythia.readString("PhaseSpace:pTHatMin = 3.");
  pythia.readString("PromptPhoton:all = off");
  pythia.readString("HardQCD:all = on");
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
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:Seed = 0");

  double R = 0.4;
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
  
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

    for (int iEvent = 0; iEvent < nEvent; iEvent++) {
      if (!pythia.next()) continue;
      
      std::vector<fastjet::PseudoJet> particles;
      std::vector<fastjet::PseudoJet> photons;
      
      for (int i = 0; i < pythia.event.size(); ++i) {
	if (pythia.event[i].isFinal()) {
	  fastjet::PseudoJet particle(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
	  particle.set_user_index(i);
	  particles.push_back(particle);
	  
	  if (pythia.event[i].id() == 22 && fabs(pythia.event[i].eta()) < 0.35) {
	    photons.push_back(particle);
	  }
	}
      }
      
      fastjet::ClusterSequence cs(particles, jet_def);
      std::vector<fastjet::PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
      
      for (const auto& jet : jets) {
	std::vector<fastjet::PseudoJet> jet_constituents = jet.constituents();
	if (jet_constituents.size() > 1) {
	  for (const auto& constituent : jet_constituents) {
	    if (pythia.event[constituent.user_index()].id() == 22) {
	      double photon_pt = constituent.pt();
	      double delta_R = deltaR(jet, constituent);
	      if (delta_R < R){
		for (int a = 0; a < npT; ++a) {
		  const double pTmin = pTlist[a];
		  const double pTmax = pTlist[a+1];
		  
		  if (pTmin < photon_pt && photon_pt < pTmax) {
		    c[a][iEvent] += 1 * info.sigmaGen();
		  }
		}
	      }
	    }
	  }
	}
      }
      pythia.stat();
    }
  }
  double mean, stan_err_mean;
  ofstream FragPhotons_pp_pTHat3_20_cross_sec;
  FragPhotons_pp_pTHat3_20_cross_sec.open("FragPhotons_pp_pTHat3_20_noprompt_deltaRcut_constituentcut_hadrondecayoff_cross_sec.txt");
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
    
    FragPhotons_pp_pTHat3_20_cross_sec << (pTmin+pTmax)/2 << " " << cross_section << " " << cross_sec_err << "\n";
    //Output
    
    cout << (pTmin+pTmax)/2 << " " << cross_section << " " << cross_sec_err << "\n";
    
  } // End second pTbin loop
  
  FragPhotons_pp_pTHat3_20_cross_sec.close();
  // Done.
  return 0;
} // End main function
