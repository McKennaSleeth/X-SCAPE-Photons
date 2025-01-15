#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdio>
#include <cmath>
#include <map>
#include <gsl/gsl_statistics.h>

struct Particle {
  int SN, PID, Status;
  double E, Px, Py, Pz, Eta, Phi;
};

double calculatePT(double px, double py){
  return std::sqrt(px*px+py*py);
}

std::map<int, std::vector<Particle>> LoadParticlesFromFile(const std::string& filename){
  double eventCount = 0.;
  std::map<int, std::vector<Particle>> eventParticles;
  
  std::ifstream ifs(filename);
  if (!ifs.is_open()) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return eventParticles;
  }
  std::string str;
  while (std::getline(ifs, str)) {
    if (str.find("Event") != std::string::npos) {
      sscanf(str.c_str(), "# Event %lf", &eventCount);
    
    } else {
      Particle p;
      sscanf(str.data(), "%d %d %d %lf %lf %lf %lf %lf %lf",
	     &p.SN, &p.PID, &p.Status, &p.E, &p.Px, &p.Py, &p.Pz, &p.Eta, &p.Phi);
      double pT = calculatePT(p.Px,p.Py);
      if(p.PID == 211 && pT > 4.5){
	eventParticles[eventCount].push_back(p);
      }
    }
  }
  return eventParticles;
}

int main() {
  std::string filename = "photon_pp_pTHat5_10_out_final_state_hadrons.dat";
  /*
  std::vector<std::string> filenames = {
					"photon_pp_pTHat5_10_out_final_state_hadrons.dat",
					"photon_pp_pTHat10_15_out_final_state_hadrons.dat",
					"photon_pp_pTHat15_20_out_final_state_hadrons.dat",
					"photon_pp_pTHat20_25_out_final_state_hadrons.dat"};
  */
  const double pTlist[] = {5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0};
  const int npT = sizeof(pTlist)/sizeof(double) - 1;
  //for (const auto& filename : filenames){
  std::map<int, std::vector<Particle>> eventSelectedParticles = LoadParticlesFromFile(filename);
  int totalPions = 0;
  std::vector<int> pionCountsPerBin(npT,0);
  for (const auto& event : eventSelectedParticles){
    //std::cout << "Event " << event.first << ": " << event.second.size() << " positive pions" << std::endl;
    totalPions += event.second.size();
    for (const Particle& pion : event.second){
      double pT = calculatePT(pion.Px, pion.Py);
      double absEta = std::fabs(pion.Eta);
      //std::cout << " Positive Pion pT: " << pT << std::endl;
      for (int i = 0; i < npT; ++i){
	if (pT >= pTlist[i] && pT <= pTlist[i+1] && absEta < 1.0){
	  pionCountsPerBin[i]++;
	}
      }
    }
  }
  std::cout << "Total number of Positive Pions in " << filename <<": " << totalPions << std::endl;
  for (int i = 0; i < npT; ++i){
    std::cout << (pTlist[i] + pTlist[i+1])/2 << " " << pionCountsPerBin[i] << std::endl;
  }
  //std::cout << "Total particles in " << filename << ": " << particles.size() << std::endl;
    /*
    for (int a = 0; a < npT; a++){
      double pTmin = pTlist[a];
      double pTmax = pTlist[a+1];
      int count = 0;
      //const long int nEvent = 100000;
      //int totalEvents = 0;
      //std::vector<std::vector<double>> c(npT, std::vector<double>(nEvent,0.0));
      //for (int iEvent = 0; iEvent < nEvent; iEvent++){	*/
    /*
    for (const auto& particle : particles){
	double pT = std::sqrt(particle.Px * particle.Px + particle.Py * particle.Py);
	//std::cout << "For event: " << bin << " with pTmin: " << pTmin << " and pTmax: " << pTmax << " Particle pT: " << pT << std::endl;
	double PID = particle.PID;
	//std::cout << "Particle PID: " << PID << std::endl;
	//if (pT >= pTmin && pT < pTmax){
	  if(PID == 211){
	    count++; //counting pions
	  }
	  //std::cout << "pT =" << pT << std::endl; 

	  //}
      }
    */
      //std::cout << "Number of particles with pT in (" << pTmin << ", " << pTmax << "): " << count << std::endl;
      
  //}
    //}	  
    
    
  //Photon_pp_pTHat0_25_tot_cross_sec.close();
  return 0;
}
