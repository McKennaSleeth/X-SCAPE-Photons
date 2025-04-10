#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <random>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>

struct Particle {
    int id, pid, status;
    double E, Px, Py, Pz;
};

const double JET_RADIUS = 0.4;
const double MIN_JET_PT = 20.0;
const double MAX_JET_PT = 30.0;
const double ETA_CUT = 0.35;

const double pTlist[] = {3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 12.0, 14.0, 16.0};
const int npT = sizeof(pTlist)/sizeof(double) - 1;

double deltaR(const fastjet::PseudoJet& jet, const fastjet::PseudoJet& particle) {
  double deta = jet.eta() - particle.eta();
  double dphi = std::abs(jet.phi() - particle.phi());
  if (dphi > M_PI) dphi = 2 * M_PI - dphi;
  return std::sqrt(deta * deta + dphi * dphi);
}

double findSigmaGen(const std::string& filename) {
  std::ifstream file(filename, std::ios::ate);
  std::string line;
  double sigmaGen = 0.0;
  
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return sigmaGen;
  }
  
  long long pos = file.tellg();
  while (pos > 0) {
    pos = std::max(0LL, pos - 1000);
    file.seekg(pos);
    if (pos > 0) {
      std::getline(file, line); // Discard partial line
    }
    while (std::getline(file, line)) {
      if (line.find("sigmaGen") != std::string::npos) {
	sscanf(line.c_str(), "# sigmaGen %lf", &sigmaGen);
	file.close();
	return sigmaGen;
      }
    }
  }
  
  file.close();
  return sigmaGen;
}

std::pair<double, std::vector<double>> processFile(const std::string& filename, long long& total_jets, long long& total_photons) {
  std::ifstream file(filename);
  std::string line;
  std::vector<double> fileCounts(npT, 0);
  
  
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return {0.0, fileCounts};
  }
  
  long long processed_events = 0;
  while (std::getline(file, line) && processed_events < 1000000) {
    int event_number, n_hadrons;
    if (sscanf(line.c_str(), "#\tEvent\t%d\tweight\t%*f\tEPangle\t%*f\tN_hadrons\t%d", &event_number, &n_hadrons) == 2) {
      std::vector<Particle> particles;
      particles.reserve(n_hadrons);
      
      for (int i = 0; i < n_hadrons; ++i) {
	if (!std::getline(file, line)) {
	  std::cerr << "Unexpected end of file" << std::endl;
	  return {0.0, fileCounts};
	}
	Particle p;
	if (sscanf(line.c_str(), "%d %d %d %lf %lf %lf %lf", &p.id, &p.pid, &p.status, &p.E, &p.Px, &p.Py, &p.Pz) == 7) {
	  particles.push_back(p);
	}
      }
      
      std::vector<fastjet::PseudoJet> input_particles;
      for (const auto& p : particles) {
	fastjet::PseudoJet pj(p.Px, p.Py, p.Pz, p.E);
	pj.set_user_index(p.pid);
	input_particles.push_back(pj);
      }
      
      fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, JET_RADIUS);
      fastjet::ClusterSequence cs(input_particles, jet_def);
      std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
      
      for (const auto& jet : jets) {
	std::vector<fastjet::PseudoJet> constituents = jet.constituents();
	
	if (constituents.size() > 1){
	  total_jets++;
	  
	  for (const auto& constituent : constituents) {
	    if (constituent.user_index() == 22) { // PID 22 is for photons
	      double delta_r = deltaR(jet, constituent);
	      if (delta_r <= JET_RADIUS) {
		double pT = constituent.pt();
		
		double eta = constituent.eta();
		
		if (std::fabs(eta) < ETA_CUT) {
		  for (int i = 0; i < npT; ++i) {
		    if (pT >= pTlist[i] && pT < pTlist[i+1]) {
		      fileCounts[i]++;
		      total_photons++;
		      break;
		    }
		  }
		}
	      }  
	    }
	  }
	}
      }
      
      processed_events++;
    }
  }
  
  file.close();
  std::cout << "Processed " << processed_events << " events from " << filename << std::endl;
  
  double sigmaGen = findSigmaGen(filename);
  std::cout << "File: " << filename << " - Extracted sigmaGen: " << sigmaGen << std::endl;
  
  return {sigmaGen, fileCounts};
}

std::vector<double> calculateCrossSection(const std::vector<double>& counts, double sigmaGen, long long eventsize) {
  std::vector<double> crossSections(npT);
  double deltaEta = 0.7;
  
  for (int i = 0; i < npT; ++i) {
    double avgPT = (pTlist[i+1] + pTlist[i]) / 2;
    double deltaE = pTlist[i+1] - pTlist[i];
    
    double count = counts[i];
    double error = std::sqrt(count);
    
    crossSections[i] = (count * sigmaGen) / (deltaE * deltaEta * 2 * M_PI * avgPT * eventsize);
  }
  
  return crossSections;
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
					"photon_pp_pTHat14_15_out_alldecayoff_1M_fixedPID_final_state_hadrons.dat"
  };
  
  long long total_jets = 0;
  long long total_photons = 0;
  std::vector<double> totalCrossSections(npT, 0);
  std::vector<double> totalSquaredErrors(npT, 0);
  
  for (const auto& filename : filenames) {
    auto [sigmaGen, fileCounts] = processFile(filename, total_jets, total_photons);
    std::vector<double> fileCrossSections = calculateCrossSection(fileCounts, sigmaGen, 1000000);
    
    for (int i = 0; i < npT; ++i) {
      totalCrossSections[i] += fileCrossSections[i];
      totalSquaredErrors[i] += (fileCrossSections[i] * fileCrossSections[i]) / (fileCounts[i] > 0 ? fileCounts[i] : 1);
    }
  }
  
  for (int i = 0; i < npT; ++i) {
    double avgPT = (pTlist[i+1] + pTlist[i]) / 2;
    double crossSection = totalCrossSections[i];
    double crossSectionError = std::sqrt(totalSquaredErrors[i]);
    std::cout << avgPT << " " << crossSection << " " << crossSectionError << "\n";
  }
  
  std::cout << "Total jets: " << total_jets << std::endl;
  std::cout << "Total photons: " << total_photons << std::endl;
  
  return 0;
}
