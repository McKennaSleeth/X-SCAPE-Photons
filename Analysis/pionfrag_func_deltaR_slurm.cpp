#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <map>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/tools/BackgroundEstimatorBase.hh>
#include <fastjet/tools/JetMedianBackgroundEstimator.hh>
#include <fastjet/tools/GridMedianBackgroundEstimator.hh>
#include <fastjet/Selector.hh>
#include <fastjet/tools/Subtractor.hh>
//#include <fastjet/tools/BackgroundSubtractor.hh>
//#include <fastjet/tools/AreaWeightedBackgroundEstimator.hh>

struct Particle {
  int id;
  int pid;
  int status;
  double E, Px, Py, Pz;
};

const int NUM_BINS = 20;
const double BIN_WIDTH = 1.0 / NUM_BINS;
const int NUM_DELTA_R_BINS = 20;
const double DELTA_R_BIN_WIDTH = 0.4 / NUM_DELTA_R_BINS;

const double MIN_JET_PT = 20.0;
const double MAX_JET_PT = 30.0;

const int NUM_PT_BINS = 20;
const double PT_BIN_WIDTH = 100 / NUM_PT_BINS;

const double JET_RADIUS = 0.4;

const double JET_ABS_ETA_MAX = 0.7;
//const double JET_PT_MIN      = 10.0;
const double R_BKG           = 0.6;
const double ETA_MAX_GHOST   = 1.5;

//const double SOFT_PARTICLE_PT_THRESHOLD = 1.0;

long long total_jets = 0;
long long total_pions = 0;

struct DNDZBin {
  double z_center;
  long long pion_count;
};

struct DNDDeltaRBin {
  double deltaR_center;
  long long pion_count;
};

struct DNDPTBin {
  double pT_center;
  long long jet_count;
};

std::vector<DNDZBin> dndz_bins(NUM_BINS);
std::vector<DNDDeltaRBin> dnd_deltaR_bins(NUM_DELTA_R_BINS);
std::vector<DNDPTBin> dndpt_bins(NUM_PT_BINS);

double deltaR(const fastjet::PseudoJet& p1, const fastjet::PseudoJet& p2) {
  double deta = p1.eta() - p2.eta();
  double dphi = std::abs(p1.phi() - p2.phi());
  if (dphi > M_PI) dphi = 2 * M_PI - dphi;
  return std::sqrt(deta * deta + dphi * dphi);
}


unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
std::normal_distribution<double> jet_distribution(0.0, 0.15);
std::normal_distribution<double> pion_distribution(0.0, 0.10);

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <file1> <file2> ... <fileN>" << std::endl;
    return 1;
  }

  std::vector<std::string> filenames;
  for (int i = 1; i < argc; ++i) {
    filenames.push_back(argv[i]);
  }

  for (int i = 0; i < NUM_BINS; ++i) {
    dndz_bins[i].z_center = i * BIN_WIDTH + BIN_WIDTH / 2;
    dndz_bins[i].pion_count = 0;
  }
    
  for (int i = 0; i < NUM_DELTA_R_BINS; ++i) {
    dnd_deltaR_bins[i].deltaR_center = i * DELTA_R_BIN_WIDTH + DELTA_R_BIN_WIDTH / 2;
    dnd_deltaR_bins[i].pion_count = 0;
  }
    
  for (int i = 0; i < NUM_PT_BINS; ++i) {
    dndpt_bins[i].pT_center = i * PT_BIN_WIDTH + PT_BIN_WIDTH / 2;
    dndpt_bins[i].jet_count = 0;
  }
  
  long long processed_events = 0;
  std::map<double, long long> jet_pt_difference_counts;
  
  for (const auto& filename : filenames) {
    std::ifstream file(filename);
    if (!file.is_open()) {
      std::cerr << "Error opening file: " << filename << std::endl;
      continue;
    }

    std::cout << "Processing file: " << filename << std::endl;

    std::string line;
    int event_number, n_hadrons;
    std::vector<Particle> particles;
    
    while (std::getline(file, line)) {
      if (sscanf(line.c_str(), "#\tEvent\t%d\tweight\t%*f\tEPangle\t%*f\tN_hadrons\t%d", &event_number, &n_hadrons) == 2) {
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

	std::vector<fastjet::PseudoJet> input_particles;
	for (const auto& p : particles) {
	  fastjet::PseudoJet pj(p.Px, p.Py, p.Pz, p.E);
	  pj.set_user_index(p.pid);
	  input_particles.push_back(pj);

	}

	fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, JET_RADIUS);
	fastjet::GhostedAreaSpec area_spec(ETA_MAX_GHOST);
	fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, area_spec);
	fastjet::ClusterSequenceArea cs(input_particles, jet_def, area_def);

	//fastjet::Selector selector_jets = fastjet::SelectorAbsEtaMax(JET_ABS_ETA_MAX) && fastjet::SelectorPtMin(MIN_JET_PT);
			
	std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt((cs.inclusive_jets()));

	fastjet::JetDefinition jet_def_bkg(fastjet::kt_algorithm, R_BKG);
	fastjet::Selector selector_bkg = fastjet::SelectorAbsEtaMax(ETA_MAX_GHOST);
	fastjet::JetMedianBackgroundEstimator bge(selector_bkg, jet_def_bkg, area_def);
	bge.set_particles(input_particles);

	fastjet::Subtractor subtractor(&bge);
	//std::vector<fastjet::PseudoJet> subtracted_jets = subtractor(jets);
	
	for (const auto& jet : jets){
	  fastjet::PseudoJet subtracted_jet = subtractor(jet);
	  //if (subtracted_jet.pt() < MIN_JET_PT || std::fabs(subtracted_jet.eta()) > JET_ABS_ETA_MAX || subtracted_jet.E() < 0.0) continue;
	  std::vector<fastjet::PseudoJet> constituents = jet.constituents();
	  if (constituents.size() > 1){
	    double original_jet_pt = jet.pt();
	  
	  //std::cout << "Subtracted jet pT: " << subtracted_jet.pt() << " and Raw jet PT: " << jet.pt() << std::endl;
	  //std::cout << "Background density (rho): " << bge.rho(jet) << std::endl;
	  //std::cout << "Jet area: " << jet.area() << std::endl;

	    double smearing_jet_factor = 1.0 + jet_distribution(generator);
	    double smeared_jet_pt = original_jet_pt * smearing_jet_factor;
	    	  
	    
	    if (smeared_jet_pt >= MIN_JET_PT && smeared_jet_pt < MAX_JET_PT){
	      total_jets++;
	      int pt_bin = std::min(static_cast<int>(original_jet_pt / PT_BIN_WIDTH), NUM_PT_BINS - 1);
	      dndpt_bins[pt_bin].jet_count++;
	    
	      double total_jet_momentum = jet.modp();	    

	      double pt_difference = original_jet_pt - subtracted_jet.pt();
	      jet_pt_difference_counts[pt_difference]++;
	      
	      //int pt_bin = std::min(static_cast<int>(original_jet_pt / PT_BIN_WIDTH), NUM_PT_BINS - 1);
	      for (const auto& constituent : constituents) {
		if (std::fabs(constituent.user_index()) == 211) { // PID 22 is for photons, abs(211) is for charged pions
		  double delta_r = deltaR(jet, constituent);
		  if (delta_r <= JET_RADIUS) {
		    double original_pion_pt = constituent.pt();
		    double smearing_pion_factor = 1.0 + pion_distribution(generator);
		    double smeared_pion_pt = original_pion_pt * smearing_pion_factor;
		    double z = smeared_pion_pt / smeared_jet_pt;
		    int z_bin = std::min(static_cast<int>(z / BIN_WIDTH), NUM_BINS - 1);
		    int delta_r_bin = std::min(static_cast<int>(delta_r / DELTA_R_BIN_WIDTH), NUM_DELTA_R_BINS - 1);
		    dndz_bins[z_bin].pion_count++;
		    dnd_deltaR_bins[delta_r_bin].pion_count++;
		    total_pions++;
		  }
		}
	      }
	    }
	  }
	}
	processed_events++;
	if (processed_events % 1000000 == 0) {
	  std::cout << "Processed " << processed_events << " events in " << filename << std::endl;
	}
      }
    }
    
    file.close();
    std::cout << "Finished processing " << filename << std::endl;
  }
  
  std::cout << "Final dN/dz distribution: " << std::endl;
  for (const auto& bin : dndz_bins) {
    double dndz = static_cast<double>(bin.pion_count) / (total_jets * BIN_WIDTH);
    double error = std::sqrt(bin.pion_count) / (total_jets * BIN_WIDTH);
    std::cout << bin.z_center << " " << dndz << " " << error << std::endl;
  }
  
  std::cout << "Final dN/ddeltaR distribution: " << std::endl;
  for (const auto& bin : dnd_deltaR_bins) {
    double dnd_deltaR = static_cast<double>(bin.pion_count) / (total_jets * DELTA_R_BIN_WIDTH);
    double error = std::sqrt(bin.pion_count) / (total_jets * DELTA_R_BIN_WIDTH);
    std::cout << bin.deltaR_center << " " << dnd_deltaR << " " << error << std::endl;
  }
  
  std::cout << "Final dN/dpT distribution: " << std::endl;
  for (const auto& bin : dndpt_bins) {
    double dndpt = static_cast<double>(bin.jet_count) / (processed_events * PT_BIN_WIDTH);
    double error = std::sqrt(bin.jet_count) / (processed_events * PT_BIN_WIDTH);
    std::cout << bin.pT_center << " " << dndpt << " " << error << std::endl;
  }

  std::ofstream output_file("pion_AuAu_pTHat15_100_out_alldecayoff_prompt_hydro_jet_pt_differences.txt");
    if (!output_file.is_open()){
      std::cerr << "Error opening output file" << std::endl;
      return 1;
    }
    //output_file << "Jet pT differences (pT - subtracted pT) and counts: " << std::endl;
  for (const auto& entry : jet_pt_difference_counts){
    output_file << entry.first << " " << entry.second << std::endl;
  }
  output_file.close();
  
  std::cout << "Total processed events: " << processed_events << std::endl;
  std::cout << "Total jets: " << total_jets << std::endl;
  std::cout << "Total pions: " << total_pions << std::endl;
  
  
  return 0;
}
