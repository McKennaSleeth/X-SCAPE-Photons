#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdio>
#include <cmath>
#include <utility>
#include <numeric>
#include <map>
#include <gsl/gsl_statistics.h>

struct Particle {
  int SN, PID, Status;
  double E, Px, Py, Pz, Eta, Phi;
};

double calculatePT(double px, double py) {
  return std::sqrt(px*px + py*py);
}

double calculateETA(double px, double py, double pz) {
  return std::atanh(pz / (std::sqrt(px*px + py*py + pz*pz)));
}

std::map<int, std::vector<Particle>> LoadParticlesFromFile(const std::string& filename, double& sigmaGen) {
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
    if (str.find("sigmaGen") != std::string::npos) {
      double value;
      if (sscanf(str.c_str(), "# sigmaGen %lf", &value) == 1) {
	sigmaGen = value;
	std::cout << "Extracted sigmaGen: " << sigmaGen << std::endl;
      }
    } else {
      Particle p;
      sscanf(str.data(), "%d %d %d %lf %lf %lf %lf", &p.SN, &p.PID, &p.Status, &p.E, &p.Px, &p.Py, &p.Pz);
      if (std::fabs(p.PID) == 211) {
	eventParticles[eventCount].push_back(p);
      }
    }
  }
  return eventParticles;
}

int main(int argc, char* argv[]) {
  
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <file1> <file2> ... <fileN>" << std::endl;
    return 1;
  }

  const double pTlist[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0};
  const int npT = sizeof(pTlist)/sizeof(double) - 1;

  std::vector<double> totalCounts(npT, 0);
  std::vector<double> totalCrossSections(npT, 0);
  std::vector<double> totalSquaredErrors(npT, 0);
  std::vector<double> totalWeights(npT, 0);
  double eventsize = 1000000;
  const double MIN_MEAN = 1e-10;
  
  //for (const auto& filename : filenames) {
  for (int i = 1; i < argc; ++i){
    double sigmaGen = 0.0;
    std::map<int, std::vector<Particle>> eventSelectedParticles = LoadParticlesFromFile(argv[i], sigmaGen);
    
    for (int j = 0; j < npT; ++j) {
      std::vector<double> pionCounts;
      pionCounts.reserve(eventsize);
      double sumCounts = 0;
      double sumSquaredCounts = 0;
      
      for (const auto& event : eventSelectedParticles) {
	double pionCount = 0;
	for (const Particle& pion : event.second) {
	  double pT = calculatePT(pion.Px, pion.Py);
	  double eta = calculateETA(pion.Px, pion.Py, pion.Pz);
	  if (pT >= pTlist[j] && pT < pTlist[j+1] && std::fabs(eta) < 0.35) {
	    pionCount += 1;
	  }
	}
	sumCounts += pionCount;
	sumSquaredCounts += pionCount * pionCount;
	//pionCounts.push_back(pionCount * sigmaGen);
      }
      
      double mean = sumCounts / eventsize;
      double variance = (sumSquaredCounts / eventsize) - (mean * mean);
      double standardDeviation = std::sqrt(variance);
      double standardError = standardDeviation / std::sqrt(eventsize) * sigmaGen;
      //double error = std::sqrt(variance / eventsize);
      
      //double mean = gsl_stats_mean(sumSquaredCounts.data(), 1, eventsize);
      //double standardDeviation = gsl_stats_sd(pionCounts.data(), 1, eventsize);
      
      
      
      double avgPT = (pTlist[j+1] + pTlist[j]) / 2;
      double deltaE = pTlist[j+1] - pTlist[j];
      double deltaEta = 0.7;
      
      //double crossSection = (mean*sigmaGen) / (deltaE * deltaEta * 2 * M_PI * avgPT);
      //double crossSectionError = (crossSection / mean) * standardError;
      
      
      //if (mean > 0) {
      //double weight = 1.0 / (crossSectionError * crossSectionError);
      double weight = 1.0 / (standardError * standardError);
      totalCounts[j] += mean * sigmaGen;
      //double weight = 1.0 / (crossSectionError * crossSectionError);
      totalSquaredErrors[j] += standardError * standardError;
      //totalCrossSections[i] += crossSection * weight;
      //totalWeights[i] += weight;
      
      
      
      //}
    }
  }
  
  for (int j = 0; j < npT; ++j) {
    //if (totalWeights[i] > 0) {
    double avgPT = (pTlist[j+1] + pTlist[j]) / 2;
    double deltaE = pTlist[j+1] - pTlist[j];
    double deltaEta = 0.7;
    
    //double finalCrossSection = totalCrossSections[i] / totalWeights[i];
    //double finalError = 1.0 / std::sqrt(totalWeights[i]);
    //std::cout << avgPT << " " << finalCrossSection << " " << finalError << "\n";
    
    //double combinedCount = totalCounts[i] / totalWeights[i];
    double combinedCount = totalCounts[j];
    //double combinedError = 1.0 / std::sqrt(totalWeights[i]);
    double combinedError = std::sqrt(totalSquaredErrors[j]);
    
    double crossSection = combinedCount / (deltaE * deltaEta * 2 * M_PI * avgPT);
    
    double crossSectionError = (crossSection / combinedCount) * combinedError;
    
    std::cout << avgPT << " " << crossSection << " " << crossSectionError << "\n";
    //std::cout << avgPT << " " << combinedCrossSection << " " << combinedError << "\n";
    //}
  }
  
  return 0;
}
