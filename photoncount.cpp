#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdio>
#include <cmath>
#include <gsl/gsl_statistics.h>

struct Particle {
  int SN, PID, Status;
  double E, Px, Py, Pz, Eta, Phi;
};
std::vector<Particle> LoadParticlesFromFile(const std::string& filename) {
  std::vector<Particle> particles;
  int eventCount = 0;
  
  //double sigmaGen = 0;
  std::ifstream ifs(filename);
  if (!ifs.is_open()) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return particles;
  }
  std::string str;
  while (std::getline(ifs, str)) {
    if (str.find("Event") != std::string::npos) {
      eventCount++;
    
    /*
    if (str.find("sigmaGen") != std::string::npos){
      double value;
      if (sscanf(str.c_str(), "# sigmaGen %lf", &value) == 1){
	sigmaGen = value;
	std::cout << "Extracted sigmaGen: " << sigmaGen << std::endl;
      }
    */
    } else {
      Particle p;
      sscanf(str.data(), "%d %d %d %lf %lf %lf %lf %lf %lf",
	     &p.SN, &p.PID, &p.Status, &p.E, &p.Px, &p.Py, &p.Pz, &p.Eta, &p.Phi);
      //if(p.PID == 22){
      particles.push_back(p);
      //}
    }
  }
  std::cout << eventCount << " events in " << filename << std::endl;
  return particles;
}
int main() {
  //std::string filename = "photon_pp_pTHat10_15_out_final_state_hadrons.dat";
  std::vector<std::string> filenames = {"photon_pp_pTHat0_5_out_final_state_hadrons.dat", "photon_pp_pTHat5_10_out_final_state_hadrons.dat", "photon_pp_pTHat10_15_out_final_state_hadrons.dat", "photon_pp_pTHat15_20_out_final_state_hadrons.dat", "photon_pp_pTHat20_25_out_final_state_hadrons.dat"};
  //std::vector<Particle> particles = LoadParticlesFromFile(filename);
  //std::cout << "Total particles: " << particles.size() << std::endl;
  //const double pTlist[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0};
  const double pTlist[] = {0.0,5.0,10.0,15.0,20.0};
  const int npT = sizeof(pTlist)/sizeof(double) - 1;
  //double mean, stan_err_mean;
  //std::ofstream Photon_pp_pTHat0_25_tot_cross_sec;
  //Photon_pp_pTHat0_25_tot_cross_sec.open("Photon_pp_pTHat0_25_tot_cross_sec.txt");
  for (const auto& filename : filenames){
    //double sigmaGen = 0.0;
    std::vector<Particle> particles = LoadParticlesFromFile(filename);
    std::cout << "Total photons in " << filename << ": " << particles.size() << std::endl;
    //std::cout << "Extracted sigmaGen from " << filename << ": " << sigmaGen << std::endl;
    for (int a = 0; a < npT; a++){
      double pTmin = pTlist[a];
      double pTmax = pTlist[a+1];
      int count = 0;
      //const long int nEvent = 100000;
      //int totalEvents = 0;
      //std::vector<std::vector<double>> c(npT, std::vector<double>(nEvent,0.0));
      //for (int iEvent = 0; iEvent < nEvent; iEvent++){	
      for (const auto& particle : particles){
	double pT = std::sqrt(particle.Px * particle.Px + particle.Py * particle.Py);
	//std::cout << "For event: " << bin << " with pTmin: " << pTmin << " and pTmax: " << pTmax << " Particle pT: " << pT << std::endl;
	double E = particle.E;
	//std::cout << "Particle E: " << E << std::endl;
	double PID = particle.PID;
	//std::cout << "Particle PID: " << PID << std::endl;
	//if (E >= pTmin && E < pTmax){
	if (pT >= pTmin && pT < pTmax){
	  if(PID == 22){
	    count++; //counting photons
	    //c[a][iEvent] += 1*sigmaGen;  //calculating cross section
	  }
	}
      }
	  //mean = gsl_stats_mean(c[a].data(), 1, nEvent);
	  /*
	    stan_err_mean = (std::sqrt(gsl_stats_variance(c[a].data(), 1, nEvent))/std::sqrt(nEvent));
	    double deltaE = pTmax - pTmin;
	    double deltaeta = 2;
	    double avgE = (pTmax+pTmin)/2;
	    double cross_section = (1/(2*3.14*avgE))*((mean)/(deltaE*deltaeta));
	    double cross_sec_err = (cross_section/mean)*stan_err_mean;
	    Photon_pp_pTHat0_25_tot_cross_sec << "pTmin " << pTmin << " pTmax " << pTmax << "\n" << "cross section " << cross_section << " cross sec err" << cross_sec_err << "\n"; 
	    */
      std::cout << "Number of particles with pT in (" << pTmin << ", " << pTmax << "): " << count << std::endl;
	    
	    //std::cout << "pTmin " << pTmin << " pTmax " << pTmax << "\n" << "cross section " << cross_section << " cross sec err " << cross_sec_err << "\n";
      
    }
  }	  
  
  
  //Photon_pp_pTHat0_25_tot_cross_sec.close();
  return 0;
}
