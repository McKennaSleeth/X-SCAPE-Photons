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

std::map<int, std::vector<Particle>> LoadParticlesFromFile(const std::string& filename, double& sigmaGen){
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
    
    }
    if (str.find("sigmaGen") != std::string::npos){
      double value;
      if (sscanf(str.c_str(), "# sigmaGen %lf", &value) == 1){
        sigmaGen = value;
        std::cout << "Extracted sigmaGen: " << sigmaGen << std::endl;
      }
    }
    else {
      Particle p;
      sscanf(str.data(), "%d %d %d %lf %lf %lf %lf %lf %lf",
	     &p.SN, &p.PID, &p.Status, &p.E, &p.Px, &p.Py, &p.Pz, &p.Eta, &p.Phi);
      double pT = calculatePT(p.Px,p.Py);
      if(std::fabs(p.PID) == 211){
	eventParticles[eventCount].push_back(p);
      }
    }
  }
  return eventParticles;
  
}

int main() {
  std::string filename = "photon_pp_pTHat5_25_out_colored_vir_fact_0.5_final_state_hadrons.dat";
  /*
  std::vector<std::string> filenames = {
					"photon_pp_pTHat5_10_out_colored_vir_fact_0.5_final_state_hadrons.dat",
					"photon_pp_pTHat10_15_out_colored_vir_fact_0.5_final_state_hadrons.dat",
					"photon_pp_pTHat15_20_out_colored_vir_fact_0.5_final_state_hadrons.dat",
					"photon_pp_pTHat20_25_out_colored_vir_fact_0.5_final_state_hadrons.dat"};
  */
  const double pTlist[] = {5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0};
  const int npT = sizeof(pTlist)/sizeof(double) - 1;
  //for (const auto& filename : filenames){
    double sigmaGen = 0.0;
    
    std::map<int, std::vector<Particle>> eventSelectedParticles = LoadParticlesFromFile(filename, sigmaGen);
    int totalPions = 0;
    double pionCountsPerBin[npT] = {0.0};
    double eventsize = 100000;
    //std::vector<std::vector<double>> c(npT, std::vector<double>(eventsize, 0.0));
    
    for (const auto& event : eventSelectedParticles){
      //std::cout << "Event " << event.first << ": " << event.second.size() << " positive pions" << std::endl;
      totalPions += event.second.size();
      int binCount = 0;
      for (const Particle& pion : event.second){
	double pT = calculatePT(pion.Px, pion.Py);
	double absEta = std::fabs(pion.Eta);
	//std::cout << " Positive Pion pT: " << pT << std::endl;
	for (int i = 0; i < npT; ++i){

	  if (pT >= pTlist[i] && pT <= pTlist[i+1] && absEta < 0.35){
	    //c[i][event.first] += 1*sigmaGen;
	    pionCountsPerBin[i] += 1*sigmaGen;
	  }
	}
	
      }
      
    }
    std::cout << "Total number of Positive Pions in " << filename <<": " << totalPions << std::endl;
    double cross_section[npT] = {0.0};
    double cross_sec_err[npT] = {0.0};
    double mean[npT] = {0.0};
    double stan_err_mean[npT] = {0.0};
    for (int i = 0; i < npT; ++i){
      //std::cout << (pTlist[i] + pTlist[i+1])/2 << " " << pionCountsPerBin[i] << std::endl;
      //double eventsize = 100000;

      mean[i] = pionCountsPerBin[i]/eventsize;

      double sumSquaredDifferences[npT] = {0.0};
      double difference[npT] = {0.0};
      double variance[npT] = {0.0};
      double standardDeviation[npT] = {0.0};
      for (const auto& event : eventSelectedParticles){
	double pionCountInEvent[npT] = {0.0};
	for (const Particle& pion : event.second){
	  double pT = calculatePT(pion.Px, pion.Py);
	  if (pT >= pTlist[i] && pT <= pTlist[i+1]){
	    pionCountInEvent[i]++;
      	  }
	}
      	difference[i] = pionCountInEvent[i] - mean[i];
	sumSquaredDifferences[i] += difference[i]*difference[i];
      }
      variance[i] = sumSquaredDifferences[i] / eventsize;
      standardDeviation[i] = std::sqrt(variance[i]);
      stan_err_mean[i] = standardDeviation[i] / std::sqrt(eventsize);
      //double mean = gsl_stats_mean(c[i].data(), 1, eventsize);
      //stan_err_mean[i] = gsl_stats_sd_m(pionCountsPerBin[i], 1, eventsize, mean[i]) / std::sqrt(eventsize);
      //double stan_err_mean = gsl_stats_sd_m(c[i].data(), 1, eventsize, mean) / std::sqrt(eventsize);
      double deltaE = pTlist[i+1] - pTlist[i];
      double deltaeta = 0.7;
      double avgE = (pTlist[i+1] + pTlist[i]) / 2;
      cross_section[i] = mean[i] / (deltaE * deltaeta * 2 * M_PI * avgE);
      //double cross_section = mean / (deltaE * deltaeta * 2 * M_PI * avgE);
      cross_sec_err[i] = (cross_section[i] / mean[i]) * (stan_err_mean[i]);
      //double cross_sec_err = (cross_section / mean) * (stan_err_mean);
      //pionplus_out << (pTlist[i]+pTlist[i+1])/2 << " " << cross_section << " " << cross_sec_err << "\n";
      //std::cout << (pTlist[i]+pTlist[i+1])/2 << " " << cross_section << " " << cross_sec_err << "\n";
      std::cout << (pTlist[i]+pTlist[i+1])/2 << " " << cross_section[i] << " " << cross_sec_err[i] << "\n";
    }
  //std::cout << "Total particles in " << filename << ": " << particles.size() << std::endl;
  //}	  
    
    
  //Photon_pp_pTHat0_25_tot_cross_sec.close();
  return 0;
}
