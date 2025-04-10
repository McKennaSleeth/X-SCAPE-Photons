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
            if (p.PID == 22) {
                eventParticles[eventCount].push_back(p);
            }
        }
    }
    return eventParticles;
}

int main() {
    std::vector<std::string> filenames = {
        "photon_pp_pTHat3_4_out_alldecayoff_1M_fixedPID_final_state_hadrons.dat",
        "photon_pp_pTHat4_5_out_alldecayoff_1M_fixedPID_final_state_hadrons.dat",
        "photon_pp_pTHat5_6_out_alldecayoff_1M_fixedPID_final_state_hadrons.dat",
	"photon_pp_pTHat6_7_out_alldecayoff_1M_fixedPID_final_state_hadrons.dat",
	"photon_pp_pTHat7_8_out_alldecayoff_1M_fixedPID_final_state_hadrons.dat",
	"photon_pp_pTHat8_9_out_alldecayoff_1M_fixedPID_final_state_hadrons.dat",
	"photon_pp_pTHat9_10_out_alldecayoff_1M_fixedPID_final_state_hadrons.dat",
	"photon_pp_pTHat10_11_out_alldecayoff_1M_fixedPID_final_state_hadrons.dat",
	"photon_pp_pTHat11_12_out_alldecayoff_1M_fixedPID_final_state_hadrons.dat",
	"photon_pp_pTHat12_13_out_alldecayoff_1M_fixedPID_final_state_hadrons.dat",
	"photon_pp_pTHat13_14_out_alldecayoff_1M_fixedPID_final_state_hadrons.dat",
	"photon_pp_pTHat14_15_out_alldecayoff_1M_fixedPID_final_state_hadrons.dat",
	"photon_pp_pTHat15_16_out_alldecayoff_1M_fixedPID_final_state_hadrons.dat"
    };

    const double pTlist[] = {3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 12.0, 14.0, 16.0};
    const int npT = sizeof(pTlist)/sizeof(double) - 1;

    std::vector<double> totalCounts(npT, 0);
    std::vector<double> totalCrossSections(npT, 0);
    std::vector<double> totalSquaredErrors(npT, 0);
    std::vector<double> totalWeights(npT, 0);
    double eventsize = 1000000;
    const double MIN_MEAN = 1e-10;

    for (const auto& filename : filenames) {
        double sigmaGen = 0.0;
        std::map<int, std::vector<Particle>> eventSelectedParticles = LoadParticlesFromFile(filename, sigmaGen);

	for (int i = 0; i < npT; ++i) {
	  std::vector<double> photonCounts;
	  photonCounts.reserve(eventsize);
	  double sumCounts = 0;
	  double sumSquaredCounts = 0;

            for (const auto& event : eventSelectedParticles) {
                double photonCount = 0;
                for (const Particle& photon : event.second) {
		  double pT = calculatePT(photon.Px, photon.Py);
		  double eta = calculateETA(photon.Px, photon.Py, photon.Pz);
		  if (pT >= pTlist[i] && pT < pTlist[i+1] && std::fabs(eta) < 0.35) {
		    photonCount += 1;
		  }
                }
                sumCounts += photonCount;
                sumSquaredCounts += photonCount * photonCount;
		//photonCounts.push_back(photonCount * sigmaGen);
            }
	    
            double mean = sumCounts / eventsize;
            double variance = (sumSquaredCounts / eventsize) - (mean * mean);
	    double standardDeviation = std::sqrt(variance);
	    double standardError = standardDeviation / std::sqrt(eventsize) * sigmaGen;
	    //double error = std::sqrt(variance / eventsize);

	    //double mean = gsl_stats_mean(sumSquaredCounts.data(), 1, eventsize);
	    //double standardDeviation = gsl_stats_sd(photonCounts.data(), 1, eventsize);
	    

	    
	    double avgPT = (pTlist[i+1] + pTlist[i]) / 2;
            double deltaE = pTlist[i+1] - pTlist[i];
            double deltaEta = 0.7;

	    //double crossSection = (mean*sigmaGen) / (deltaE * deltaEta * 2 * M_PI * avgPT);
	    //double crossSectionError = (crossSection / mean) * standardError;
	    
	    
            //if (mean > 0) {
	      //double weight = 1.0 / (crossSectionError * crossSectionError);
	      double weight = 1.0 / (standardError * standardError);
	      totalCounts[i] += mean * sigmaGen;
	      //double weight = 1.0 / (crossSectionError * crossSectionError);
	      totalSquaredErrors[i] += standardError * standardError;
	      //totalCrossSections[i] += crossSection * weight;
	      //totalWeights[i] += weight;
	    
	    
	    
	      //}
        }
    }

    for (int i = 0; i < npT; ++i) {
      //if (totalWeights[i] > 0) {
            double avgPT = (pTlist[i+1] + pTlist[i]) / 2;
            double deltaE = pTlist[i+1] - pTlist[i];
            double deltaEta = 0.7;

	    //double finalCrossSection = totalCrossSections[i] / totalWeights[i];
	    //double finalError = 1.0 / std::sqrt(totalWeights[i]);
	    //std::cout << avgPT << " " << finalCrossSection << " " << finalError << "\n";
	    
            //double combinedCount = totalCounts[i] / totalWeights[i];
	    double combinedCount = totalCounts[i];
            //double combinedError = 1.0 / std::sqrt(totalWeights[i]);
	    double combinedError = std::sqrt(totalSquaredErrors[i]);

            double crossSection = combinedCount / (deltaE * deltaEta * 2 * M_PI * avgPT);
	    
            double crossSectionError = (crossSection / combinedCount) * combinedError;

            std::cout << avgPT << " " << crossSection << " " << crossSectionError << "\n";
	    //std::cout << avgPT << " " << combinedCrossSection << " " << combinedError << "\n";
	    //}
    }

    return 0;
}
