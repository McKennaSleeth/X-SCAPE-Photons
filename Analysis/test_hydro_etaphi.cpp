#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <sstream>

const int NUM_MOMENTUM_FRACTION_BINS = 10; // Define number of bins for momentum fractions
const int NUM_DELTA_R_BINS = 10; // Define number of bins for delta R
const double MOMENTUM_FRACTION_BIN_WIDTH = 0.1; // Define width of momentum fraction bins
const double DELTA_R_BIN_WIDTH = 0.1; // Define width of delta R bins


struct JetInfo {
  int event_number;
  double eta;
  double phi;
  double pT;
};

struct Particle {
  int id;
  int pid;
  int status;
  double E, Px, Py, Pz;

  double pT() const {
    return std::sqrt(Px * Px + Py * Py);
  }

  double eta() const {
    double p = pT();
    return (p > 0) ? 0.5 * std::log((p + Pz) / (p - Pz)) : 0.0; // Simplified eta calculation
  }
  
  double phi() const {
    return std::atan2(Py, Px);
  }
};

struct MomentumFractionBin {
    double center;
    long long pion_count;
};

struct DeltaRBin {
    double center;
    long long pion_count;
};

std::vector<MomentumFractionBin> momentumFractionBins(NUM_MOMENTUM_FRACTION_BINS);
std::vector<DeltaRBin> deltaRBins(NUM_DELTA_R_BINS);

void initializeBins() {
    for (int i = 0; i < NUM_MOMENTUM_FRACTION_BINS; ++i) {
        momentumFractionBins[i].center = i * MOMENTUM_FRACTION_BIN_WIDTH + MOMENTUM_FRACTION_BIN_WIDTH / 2;
        momentumFractionBins[i].pion_count = 0;
    }
    for (int i = 0; i < NUM_DELTA_R_BINS; ++i) {
        deltaRBins[i].center = i * DELTA_R_BIN_WIDTH + DELTA_R_BIN_WIDTH / 2;
        deltaRBins[i].pion_count = 0;
    }
}

const double CONE_RADIUS = 0.4; // Define the radius of the cone for matching

double deltaR(double eta1, double phi1, double eta2, double phi2) {
    double deta = eta1 - eta2;
    double dphi = std::abs(phi1 - phi2);
    if (dphi > M_PI) dphi = 2 * M_PI - dphi;
    return std::sqrt(deta * deta + dphi * dphi);
}

bool isPion(const Particle& particle) {
  return std::abs(particle.pid) == 211; // Check if it is a pion (PID 211)
}

bool particleInCone(const Particle& particle, const JetInfo& jet) {
  return deltaR(particle.eta(), particle.phi(), jet.eta, jet.phi) < CONE_RADIUS;
}

void processFile(const std::string& jetFileName, const std::vector<std::string>& particleFiles) {
  initializeBins();
  // Read subtracted jet information
  std::ifstream jetFile(jetFileName);
  if (!jetFile.is_open()) {
    std::cerr << "Error opening jet file: " << jetFileName << std::endl;
    return;
  }
  
  std::vector<JetInfo> jets;
  std::string line;
  jets.reserve(54365);
  while (std::getline(jetFile, line)) {
    std::stringstream ss(line);
    JetInfo jet;
    if (ss >> jet.event_number >> jet.eta >> jet.phi >> jet.pT) {
      jets.push_back(jet);
    }
  }
  jetFile.close();
  
  // Process each particle file
  for (const auto& particleFileName : particleFiles) {
    std::ifstream particleFile(particleFileName);
    if (!particleFile.is_open()) {
      std::cerr << "Error opening particle file: " << particleFileName << std::endl;
      continue;
    }
    
    std::cout << "Processing particle file: " << particleFileName << std::endl;
    
    // Read particle information
    std::vector<Particle> particles;
    
    while (std::getline(particleFile, line)) {
      Particle particle;
      if (sscanf(line.c_str(), "%d %d %d %lf %lf %lf %lf", 
		 &particle.id, &particle.pid, &particle.status, 
		 &particle.E, &particle.Px, &particle.Py, &particle.Pz) == 7) {
	particles.push_back(particle);
      }
    }
    particleFile.close();
    
    // Find matches within cones
    for (const auto& jet : jets) {
      std::cout << "Event: " << jet.event_number << std::endl;
      for (const auto& particle : particles) {
	if (isPion(particle) && particleInCone(particle, jet)) {
	  double momentumFraction = particle.pT() / jet.pT; // Calculate momentum fraction
	  double dr = deltaR(particle.eta(), particle.phi(), jet.eta, jet.phi); // Calculate delta R
	  int momentumBinIndex = static_cast<int>(momentumFraction / MOMENTUM_FRACTION_BIN_WIDTH);
	  if (momentumBinIndex >= 0 && momentumBinIndex < NUM_MOMENTUM_FRACTION_BINS) {
	    momentumFractionBins[momentumBinIndex].pion_count++;
	  }
	  
	  // Bin delta R
	  int deltaRBinIndex = static_cast<int>(dr / DELTA_R_BIN_WIDTH);
	  if (deltaRBinIndex >= 0 && deltaRBinIndex < NUM_DELTA_R_BINS) {
	    deltaRBins[deltaRBinIndex].pion_count++;
	  }
	  //std::cout << "Event: " << jet.event_number << std::endl;
	  /*
	  std::cout << "Event " << jet.event_number 
		    << ": Found particle in cone with pT " << particle.pT()
		    << " (eta: " << particle.eta() 
		    << ", phi: " << particle.phi() << ")" 
		    << std::endl;
	  */
	}
      }
    }
  }
  std::cout << "Momentum Fraction Binning Results: " << std::endl;
  for (const auto& bin : momentumFractionBins) {
    std::cout << "Fraction Center: " << bin.center << ", Count: " << bin.pion_count << std::endl;
  }
  
  std::cout << "Delta R Binning Results: " << std::endl;
  for (const auto& bin : deltaRBins) {
    std::cout << "Delta R Center: " << bin.center << ", Count: " << bin.pion_count << std::endl;
  }
}

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <jet_file> <particle_file1> <particle_file2> ... <particle_fileN>" << std::endl;
    return 1;
  }
  
  std::string jetFileName = argv[1];
  std::vector<std::string> particleFiles;
  for (int i = 2; i < argc; ++i) {
    particleFiles.push_back(argv[i]);
  }
  
  processFile(jetFileName, particleFiles);
  return 0;
}
