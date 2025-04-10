#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>


struct Particle {
  int id;
  int pid;
  int status;
  double E, Px, Py, Pz;
};

const int NUM_BINS = 20;
const double BIN_WIDTH = 0.05;
const double DELTA_R_BIN_WIDTH = 0.02;
const int NUM_DELTA_R_BINS = 20;


const double MIN_JET_PT = 10.0;
const double MAX_JET_PT = 20.0;

const double JET_RADIUS = 0.4;

struct DNDZBin {
  double z_center;
  long long pion_count;
};

struct DNDDeltaRBin {
  double deltaR_center;
  long long pion_count;
};

double deltaR(const fastjet::PseudoJet& p1, const fastjet::PseudoJet& p2){
  double deta = p1.eta() - p2.eta();
  double dphi = std::abs(p1.phi() - p2.phi());
  if (dphi > M_PI) dphi = 2 * M_PI - dphi;
  return std::sqrt(deta * deta + dphi * dphi);
}

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
std::normal_distribution<double> jet_distribution(0.0, 0.15);
std::normal_distribution<double> pion_distribution(0.0, 0.01);

int main() {
    std::ifstream file("photon_pp_pTHat15_40_out_alldecayoff_noprompt_10M_final_state_hadrons.dat");
    std::string line;
    int event_number, n_hadrons;
    std::vector<Particle> particles;

    std::vector<DNDZBin> dndz_bins(NUM_BINS);
    std::vector<DNDDeltaRBin> dnd_deltaR_bins(NUM_DELTA_R_BINS);
    for (int i = 0; i < NUM_BINS; ++i){
      dndz_bins[i].z_center = (i + 0.5) * BIN_WIDTH;
      dndz_bins[i].pion_count = 0;
      
    }
    
    for (int i = 0; i < NUM_DELTA_R_BINS; ++i) {
      dnd_deltaR_bins[i].deltaR_center = (i + 0.5) * DELTA_R_BIN_WIDTH;
      dnd_deltaR_bins[i].pion_count = 0;
    }
    
    long long total_jets = 0;
    long long total_pions = 0;
    long long processed_events = 0;
    
    if (!file.is_open()) {
      std::cerr << "Error opening file" << std::endl;
      return 1;
    }
    
    while (std::getline(file, line) && processed_events < 1000000) {
      if (sscanf(line.c_str(), "#\tEvent\t%d\tweight\t%*f\tEPangle\t%*f\tN_hadrons\t%d", &event_number, &n_hadrons) == 2) {
	//std::cout << "Event " << event_number << ": " << n_hadrons << " hadrons" << std::endl;
        
	particles.clear();
	particles.reserve(n_hadrons);
	
	for (int i = 0; i < n_hadrons; ++i) {
	  if (!std::getline(file, line)) {
	    std::cerr << "Unexpected end of file" << std::endl;
	    return 1;
	  }
	  
	  Particle p;
	  if (sscanf(line.c_str(), "%d %d %d %lf %lf %lf %lf", &p.id, &p.pid, &p.status, &p.E, &p.Px, &p.Py, &p.Pz) == 7) {
	    particles.push_back(p);
	  } else {
	    std::cerr << "Error parsing particle data" << std::endl;
	    return 1;
	  }
	}
	
	// Here you can process the particles for each event
	std::vector<fastjet::PseudoJet> input_particles;
	for (const auto& p : particles) {
	  //std::cout << "Hadron ID: " << p.id << ", PID: " << p.pid << ", E: " << p.E << std::endl;
	  //input_particles.push_back(fastjet::PseudoJet(p.Px, p.Py, p.Pz, p.E));
	  fastjet::PseudoJet pj(p.Px, p.Py, p.Pz, p.E);
	  pj.set_user_index(p.pid);
	  input_particles.push_back(pj);
	}
	double R = 0.4;
	fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
	fastjet::ClusterSequence cs(input_particles, jet_def);
	std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
	
	//std::cout << "Jets found in event " << event_number << ":" << std::endl;
	for (const auto& jet : jets) {
	  //double jet_pt = jet.pt();
	  double original_jet_pt = jet.pt();
	  double smearing_jet_factor = 1.0 + jet_distribution(generator);
	  double smeared_jet_pt = original_jet_pt * smearing_jet_factor;
	  std::vector<fastjet::PseudoJet> constituents = jet.constituents();
	  
	  if (smeared_jet_pt >= MIN_JET_PT && smeared_jet_pt < MAX_JET_PT && constituents.size() > 1){
	    total_jets++;
	    
	    //for (size_t i = 0; i < jets.size(); i++) {
	    /*
	      std::cout << "Jet " << i << ": pt = " << jets[i].pt() 
	      << ", eta = " << jets[i].eta() 
	      << ", phi = " << jets[i].phi() << std::endl;
	    */
	    //int pion_count = 0;
	    //double pion_momentum_fraction = 0.0;
	    double total_jet_momentum = jet.modp();
	    //std::vector<double> pion_z_values;
	    for (const auto& constituent : constituents) {
	      if (constituent.user_index() == std::fabs(211)) {  // PID 22 is for pions
		double delta_r = deltaR(jet, constituent);
		if (delta_r <= JET_RADIUS){
		  //pion_count++;
		  //pion_momentum_fraction += constituent.modp() / total_jet_momentum;
		  double original_pion_pt = constituent.pt();
		  double smearing_pion_factor = 1.0 + pion_distribution(generator);
		  double smeared_pion_pt = original_pion_pt * smearing_pion_factor;
		  //double z = constituent.modp() / total_jet_momentum;
		  double z = smeared_pion_pt / smeared_jet_pt;
		  int z_bin = std::min(static_cast<int>(z / BIN_WIDTH), NUM_BINS - 1);
		  int delta_r_bin = std::min(static_cast<int>(delta_r / DELTA_R_BIN_WIDTH), NUM_DELTA_R_BINS - 1);
		  dndz_bins[z_bin].pion_count++;
		  dnd_deltaR_bins[delta_r_bin].pion_count++;
		  total_pions++;
		  //pion_z_values.push_back(z);
		  //std::cout << " Pion " << pion_count << ": z = " << z << std::endl;
		}
	      }
	    }
	  }
	  //std::cout << "  Pions in jet: " << pion_count << std::endl;
	  //std::cout << "Pion z: " << pion_momentum_fraction << std::endl;
	  // Calculate dN/dz (assuming the amount of bins is the amount of events)
	}
	processed_events++;
	if (processed_events % 10000 == 0){
	  std::cout << "Processed " << processed_events << " events" << std::endl;
	}
	
      }
    }
    
    
    file.close();
    
    std::cout << "Final dN/dz distribution: " << std::endl;
    for (const auto& bin : dndz_bins){
      double dndz = static_cast<double>(bin.pion_count) / (total_jets * BIN_WIDTH);
      double error = std::sqrt(bin.pion_count) / (total_jets * BIN_WIDTH);
      //std::cout << "z = " << bin.z_center << ": dN/dz = " << dndz << " +/- " << error << std::endl;
      std::cout << bin.z_center << " " << dndz << " " << error << std::endl;
    }
    
    std::cout << "Final dN/ddeltaR distribution: " << std::endl;
    for (const auto& bin : dnd_deltaR_bins){
      double dnd_deltaR = static_cast<double>(bin.pion_count) / (total_jets * DELTA_R_BIN_WIDTH);
      double error = std::sqrt(bin.pion_count) / (total_jets * DELTA_R_BIN_WIDTH);
      //std::cout << "z = " << bin.z_center << ": dN/dz = " << dndz << " +/- " << error << std::endl;
      std::cout << bin.deltaR_center << " " << dnd_deltaR << " " << error << std::endl;
    }
    
    return 0;
}
