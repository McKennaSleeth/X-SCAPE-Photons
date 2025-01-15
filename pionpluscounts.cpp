#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <string>
#include <utility>

// Define the pion ID for filtering
const int PION_POSITIVE_ID = 211;

struct Particle {
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

// Function to read particles and sigmaGen based on pTHat bin
void readDataFromBin(const std::pair<int, int>& pTHatBin, std::vector<Particle>& particles, double& sigmaGen) {
  int eventCount = 0;
  std::string filename = "photon_pp_pTHat" + std::to_string(pTHatBin.first) + "_" + std::to_string(pTHatBin.second) + "_out_final_state_hadrons.dat";
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Could not open file: " << filename << std::endl;
        return;
    }
    
    std::string str;
    while (std::getline(file,str)){
      /*
      if (str.find("Event") != std::string::npos){
	eventCount++;
      }
      */
      if (str.find("sigmaGen") != std::string::npos){
	double value;
	if (sscanf(str.c_str(), "# sigmaGen %lf", &value) == 1){
	  sigmaGen = value;
	  std::cout << "Extracted sigmaGen: " << sigmaGen << std::endl;
	}
      }
    }
      /*
      else {
	Particle p;
	sscanf(str.data(), "%d %d %d %lf %lf %lf %lf %lf %lf",
	       &p.SN, &p.PID, &p.Status, &p.E, &p.Px, &p.Py, &p.Pz, &p.Eta, &p.Phi);
	particles.push_back(p);
      }
    }
    std::cout << eventCount << " events in " << file << std::endl;
    */
    int particleCount = 0;
    int eventNumber, numParticles;
    while (file >> eventNumber >> numParticles) {
      eventCount++;
      for (int i = 0; i < numParticles; ++i) {
	Particle particle;
	file >> particle.id >> particle.status >> particle.energy
	     >> particle.px >> particle.py >> particle.pz
	     >> particle.eta >> particle.phi;
	particles.push_back(particle);
	particleCount++;
      }
    }
    
    // Read sigmaGen, which is expected to be at the end of the file
    //sigmaGen = 0.0;
    //if (file >> sigmaGen){
    //std::cout << "Extracted sigmaGen: " << sigmaGen << std::endl;
    //} else {
    //std::cerr << "Failed to read sigmaGen from file: " << filename << std::endl;
    //}

    file.close();
    //return particles;
}

int main() {
    // Vector of pairs, each pair holds the min and max pTHat for each bin
    std::vector<std::pair<int, int>> pTHatBins = {
        {5, 10},
        {10, 15},
        {15, 20},
        {20, 25}
    };

    // Map to hold the count of particles in each pT bin for each event
    std::map<int, std::map<int, int>> eventPTBinCounts;

    // Open output file
    std::ofstream outputFile("test_output.txt");
    if (!outputFile.is_open()) {
        std::cerr << "Could not open output file." << std::endl;
        return 1;
    }

    // Loop over all pTHat bins
    for (const auto& pTHatBin : pTHatBins) {
        std::vector<Particle> particles;
        double sigmaGen = 0.0;
	
        // Read particles and sigmaGen from the current pTHat bin
        readDataFromBin(pTHatBin, particles, sigmaGen);
        

		
        // Process particles
        //int eventNumber = -1;
	std::map<int, int> pTBinCounts;
	for (int i = 5; i <= 20; ++i){
	  pTBinCounts[i] = 0;
	}
        for (const Particle& p : particles) {
            // Calculate pT
            double pT = calculatePT(p.px, p.py);
	    double PID = p.id;
            // Find the corresponding pT bin (integer value)
            // and ensure it's within the range of 5 to 20 GeV
            int pTBin = static_cast<int>(pT);
            if (pTBin >= 5 && pTBin <= 20) {
	      // Increment the count for the pT bin for the current event
	      if (PID == 211){
		pTBinCounts[pTBin]++;
		//eventPTBinCounts[eventNumber][pTBin]++;
	      }
	    }   
	}

	for (const auto& binCount : pTBinCounts){
	  if (binCount.second > 0){
	    std::cout << pTHatBin.first << " " << pTHatBin.second << " " << binCount.first + 0.5 << " " << binCount.second << std::endl;
	  }
	}
    }
    // Close output file
    outputFile.close();
    
    return 0;
}
