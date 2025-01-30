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

//! adding fastjet headers
#include "fastjet/PseudoJet.hh"

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
            double eta = calculateETA(p.Px, p.Py, p.Pz);
            double pT = calculatePT(p.Px, p.Py);
            if (p.PID == 111) {
                eventParticles[eventCount].push_back(p);
            }
	}
    }
    return eventParticles;
}

int main() {

  //! adding a test pseudojet obj
  //fastjet::PseudoJet j1;


    std::string filename = "pion_pp_pTHat15_30_out_final_state_hadrons.dat";
    const double pTlist[] = {5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0};
        const int npT = sizeof(pTlist)/sizeof(double) - 1;

    double sigmaGen = 0.0;
    std::map<int, std::vector<Particle>> eventSelectedParticles = LoadParticlesFromFile(filename, sigmaGen);
    double eventsize = 100000;
    for (int i = 0; i < npT; ++i) {
      std::vector<double> pionCounts;
      pionCounts.reserve(eventsize);

        for (const auto& event : eventSelectedParticles) {
            double pionCount = 0;
            for (const Particle& pion : event.second) {
                double pT = calculatePT(pion.Px, pion.Py);
                double eta = calculateETA(pion.Px, pion.Py, pion.Pz);
                if (pT >= pTlist[i] && pT < pTlist[i+1] && std::fabs(eta) < 0.35) {
                    pionCount += 1* sigmaGen;
                }
            }
            pionCounts.push_back(pionCount);
        }

        double mean = gsl_stats_mean(pionCounts.data(), 1, eventsize);
        double standardDeviation = gsl_stats_sd(pionCounts.data(), 1, eventsize);
        double standardError = standardDeviation / std::sqrt(eventsize);
	
        double deltaE = pTlist[i+1] - pTlist[i];
        double deltaEta = 0.7;
        double avgPT = (pTlist[i+1] + pTlist[i]) / 2;

        double crossSection = (mean) / (deltaE * deltaEta * 2 * M_PI * avgPT );
        double crossSectionError = (crossSection / mean) * standardError;

        std::cout << avgPT << " " << crossSection << " " << crossSectionError << "\n";
    }

    return 0;
}
