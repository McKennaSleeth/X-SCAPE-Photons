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
#include <fastjet/tools/JetMedianBackgroundEstimator.hh>
#include <fastjet/tools/Subtractor.hh>
//#include <fastjet/tools/BackgroundSubtractor.hh>
#include <fastjet/AreaDefinition.hh>
//#include <fastjet/tools/AreaWeightedBackgroundEstimator.hh>
#include <fastjet/ClusterSequenceArea.hh>

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

long long total_jets = 0;
long long total_photons = 0;

struct DNDZBin {
  double z_center;
  long long photon_count;
};

struct DNDDeltaRBin {
  double deltaR_center;
  long long photon_count;
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
std::normal_distribution<double> photon_distribution(0.0, 0.10);

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
    dndz_bins[i].photon_count = 0;
  }
    
  for (int i = 0; i < NUM_DELTA_R_BINS; ++i) {
    dnd_deltaR_bins[i].deltaR_center = i * DELTA_R_BIN_WIDTH + DELTA_R_BIN_WIDTH / 2;
    dnd_deltaR_bins[i].photon_count = 0;
  }
    
  for (int i = 0; i < NUM_PT_BINS; ++i) {
    dndpt_bins[i].pT_center = i * PT_BIN_WIDTH + PT_BIN_WIDTH / 2;
    dndpt_bins[i].jet_count = 0;
  }
  
  long long processed_events = 0;

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

	double R = 0.4;

	//subtractor.set_background_estimator(background_estimator);
	//fastjet::AreaDefinition area_def(fastjet::active_area);
	fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
	//fastjet::AreaDefinition area_def(fastjet::active_area, fastjet::GhostedAreaSpec(R));
	//fastjet::ClusterSequenceArea cs(input_particles, jet_def, area_def);
	fastjet::ClusterSequence cs(input_particles, jet_def);
	std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());

	fastjet::JetMedianBackgroundEstimator bge(fastjet::SelectorIdentity());
	//background_estimator.set_particles(input_particles);

	fastjet::Subtractor subtractor(&bge);

	
	for (const auto& jet : jets){

	  PseudoJet subtracted_jet
	  //double jet_area = jet.area();

	  //double background_density = background_estimator.background(jet);
	  //double background_density = 0.0;
	  //for (const auto& particle : input_particles){
	  //  background_density += particle.pt();
	  //}
	  //background_density /= input_particles.size();
	  
	  //double background_energy = background_density * jet_area;
	  //fastjet::PseudoJet corrected_jet = jet;

	  //corrected_jet.reset_momentum(corrected_jet.px(), corrected_jet.py(), corrected_jet.pz(), corrected_jet.E() - background_energy);
	  //subtractor.subtract(jet);
	  //}
	
	//for (const auto& jet : jets) {
	  double corrected_jet_pt = corrected_jet.pt();
	  double smearing_jet_factor = 1.0 + jet_distribution(generator);
	  double smeared_jet_pt = corrected_jet_pt * smearing_jet_factor;

	  
	  std::vector<fastjet::PseudoJet> constituents = jet.constituents();

	  if (smeared_jet_pt >= MIN_JET_PT && smeared_jet_pt < MAX_JET_PT && constituents.size() > 1) {
	    total_jets++;
	    int pt_bin = std::min(static_cast<int>(corrected_jet_pt / PT_BIN_WIDTH), NUM_PT_BINS - 1);
            dndpt_bins[pt_bin].jet_count++;
	    double total_jet_momentum = corrected_jet.modp();

	    for (const auto& constituent : constituents) {
	      if (constituent.user_index() == 22) { // PID 22 is for photons
		double delta_r = deltaR(jet, constituent);
		if (delta_r <= JET_RADIUS) {
		  double original_photon_pt = constituent.pt();
		  double smearing_photon_factor = 1.0 + photon_distribution(generator);
		  double smeared_photon_pt = original_photon_pt * smearing_photon_factor;
		  double z = smeared_photon_pt / smeared_jet_pt;
		  int z_bin = std::min(static_cast<int>(z / BIN_WIDTH), NUM_BINS - 1);
		  int delta_r_bin = std::min(static_cast<int>(delta_r / DELTA_R_BIN_WIDTH), NUM_DELTA_R_BINS - 1);
		  dndz_bins[z_bin].photon_count++;
		  dnd_deltaR_bins[delta_r_bin].photon_count++;
		  total_photons++;
		}
	      }
	    }
	  }
	}

	processed_events++;
	if (processed_events % 10000000 == 0) {
	  std::cout << "Processed " << processed_events << " events in " << filename << std::endl;
	}
      }
    }

    file.close();
    std::cout << "Finished processing " << filename << std::endl;
  }

  std::cout << "Final dN/dz distribution: " << std::endl;
  for (const auto& bin : dndz_bins) {
    double dndz = static_cast<double>(bin.photon_count) / (total_jets * BIN_WIDTH);
    double error = std::sqrt(bin.photon_count) / (total_jets * BIN_WIDTH);
    std::cout << bin.z_center << " " << dndz << " " << error << std::endl;
  }

  std::cout << "Final dN/ddeltaR distribution: " << std::endl;
  for (const auto& bin : dnd_deltaR_bins) {
    double dnd_deltaR = static_cast<double>(bin.photon_count) / (total_jets * DELTA_R_BIN_WIDTH);
    double error = std::sqrt(bin.photon_count) / (total_jets * DELTA_R_BIN_WIDTH);
    std::cout << bin.deltaR_center << " " << dnd_deltaR << " " << error << std::endl;
  }

  std::cout << "Final dN/dpT distribution: " << std::endl;
    for (const auto& bin : dndpt_bins) {
      double dndpt = static_cast<double>(bin.jet_count) / (processed_events * PT_BIN_WIDTH);
      double error = std::sqrt(bin.jet_count) / (processed_events * PT_BIN_WIDTH);
      std::cout << bin.pT_center << " " << dndpt << " " << error << std::endl;
    }

  std::cout << "Total processed events: " << processed_events << std::endl;
  std::cout << "Total jets: " << total_jets << std::endl;
  std::cout << "Total photons: " << total_photons << std::endl;

  return 0;
}
