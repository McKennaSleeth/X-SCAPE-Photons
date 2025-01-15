#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <map>
#include <numeric>
#include <string>
#include <utility> // for std::pair
#include <gsl/gsl_statistics.h>

// Define the pion ID for filtering
const int PION_POSITIVE_ID = 211;

struct Particle {
  int N;
  int id;
  int status;
  double energy;
  double px;
  double py;
  double pz;
  double eta;
  double phi;
};

// Function to calculate pT
double calculatePT(double px, double py) {
  return std::sqrt(px * px + py * py);
}

// Function to read particles and sigmaGen from a file
bool readDataFromFile(const std::pair<int, int>& pTHatBin, std::vector<Particle>& particles, double& sigmaGen) {
  std::string filename = "photon_pp_pTHat" + std::to_string(pTHatBin.first) + "_" + std::to_string(pTHatBin.second) + "_out_final_state_hadrons.dat";

  //std::vector<Particle> particles;
  int eventCount = 0;
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Could not open file: " << filename << std::endl;
    return false;
  }
  std::string str;
  while (std::getline(file,str)){
    if (str.find("Event") != std::string::npos){
      eventCount++;
    }
    
    if (str.find("sigmaGen") != std::string::npos){
      double value;
      if (sscanf(str.c_str(), "# sigmaGen %lf", &value) == 1){
	sigmaGen = value;
	std::cout << "Extracted sigmaGen: " << sigmaGen << std::endl;
      }
    }
    else {
      Particle particle;
      sscanf(str.data(), "%d %d %d %lf %lf %lf %lf %lf %lf", &particle.N, &particle.id, &particle.status, &particle.energy, &particle.px, &particle.py, &particle.pz, &particle.eta, &particle.phi);
      particles.push_back(particle);
    }
  }
  
  //std::cout << eventCount << "events in " << file << std::endl;
  file.close();
  return true;
}

int main() {
    // Vector of pairs, each pair holds the min and max pTHat for each bin
    std::vector<std::pair<int, int>> pTHatBins = {
        {5, 10},
        {10, 15},
        {15, 20},
        {20, 25}
    };

    std::map<int, int> pTBinCounts;
    for (int i = 5; i <= 20; ++i){
      pTBinCounts[i];
    }
    // Open output file
    std::ofstream outputFile("pionplus_pp_pTHat5_25_cross_secv3.txt");
    if (!outputFile.is_open()) {
      std::cerr << "Could not open output file." << std::endl;
      return 1;
    }

    // Loop over all pTHat ranges
    //for (size_t i = 0; i < pTHatBins.size(); ++i) {
    for (const auto& pTHatBin : pTHatBins){
      std::vector<Particle> particles;
      double sigmaGen = 0.0;
      
      //const double pTlist[] = {5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0};
      //const int npT = sizeof(pTlist)/sizeof(double) - 1;

      // Read particles and sigmaGen from the current file
      if (!readDataFromFile(pTHatBin, particles, sigmaGen)) {
	continue; // Skip to the next file if reading fails
      }

      std::map<int, std::vector<double>> pTBinValues;
      
      // Process particles
      //std::vector<double> pTValues;
      //double pTSum = 0.0;
      //int nEvent = 0;
      //for (int a = 0; a < npT; a++){
      //double pTmin = pTlist[a];
      //double pTmax = pTlist[a+1];
      
      for (const Particle& p : particles) {
	// Calculate pT
	double pT = calculatePT(p.px, p.py);
	double absEta = std::fabs(p.eta);
	int pTBin = static_cast<int>(pT);
	if (pTBin >= 5 && pTBin <= 20 && absEta < 1){
	  //pTBinCounts[pTBin]++;
	  //pTValues.push_back(pT);
	  //pTSum += pT;
	  //nEvent++;
	  pTBinCounts[pTBin]++;
	  pTBinValues[pTBin].push_back(pT);
	}
      }
    
	
	    
    // Calculate invariant differential cross section
    // Placeholder for actual calculation, which should use sigmaGen and pTHat range
      for (int pTBin = 5; pTBin <= 20; ++pTBin){
	double pTSum = std::accumulate(pTBinValues[pTBin].begin(), pTBinValues[pTBin].end(), 0.0);
	int nEvent = pTBinValues[pTBin].size();
	double mean = (nEvent > 0) ? (gsl_stats_mean(pTBinValues[pTBin].data(), 1, nEvent)) : 0.0;
	double stan_err_mean = (nEvent > 1) ? (std::sqrt(gsl_stats_variance(pTBinValues[pTBin].data(), 1, nEvent)) / std::sqrt(nEvent)) : 0.0;
	double deltaE = 1.0;
	double deltaeta = 2.0;
	double avgE = pTBin + 0.5;
	double cross_section = sigmaGen * (1 / (2 * M_PI * avgE)) * ((mean) / (deltaE * deltaeta));
	double cross_sec_err = (cross_section / mean) * stan_err_mean;
	
	//double invariantCrossSection = (sigmaGen * pT); // Replace with actual formula
    
    // Write results to the output file
	//outputFile << "pTHatRange: (" << pTHatRanges[i].first << ", " << pTHatRanges[i].second << ") \n" << pT << cross_section <<  std::endl;
	if (nEvent > 0){
	  outputFile << pTBin + 0.5 << " " << cross_section << " " << cross_sec_err << "\n";
	  std::cout << pTBin + 0.5 << " " << cross_section << " " << cross_sec_err << "\n";
	}
      }
    }
    
    // Close output file
    outputFile.close();
    
    return 0;
}
