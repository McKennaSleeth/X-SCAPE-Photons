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

const double PT_DIFF_BIN_SIZE = 5.0;
const double RAW_PT_BIN_SIZE = 5.0;
const double SUB_PT_BIN_SIZE = 5.0;

const double PION_PT_BIN_SIZE = 0.5;
const double PARTICLE_PT_BIN_SIZE = 0.5;

const double RAW_CUT_PT_BIN_SIZE = 5.0;
const double SUB_CUT_PT_BIN_SIZE = 5.0;

const double RHO_BIN_SIZE = 1.0;

const double MIN_JET_PT = 20.0;
const double MAX_JET_PT = 30.0;

const int NUM_PT_BINS = 20;
const double PT_BIN_WIDTH = 100 / NUM_PT_BINS;

const double JET_RADIUS = 0.4;

const double JET_ABS_ETA_MAX = 1 - JET_RADIUS;
const double JET_PT_CUT      = 10.0;
const double R_BKG           = 0.4;
const double ETA_MAX_GHOST   = 1 - JET_RADIUS;

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

std::vector<fastjet::PseudoJet> filter_overlapping_jets(const std::vector<fastjet::PseudoJet>& jets, double min_delta_R) {
  std::vector<fastjet::PseudoJet> filtered_jets;

  for (size_t i = 0; i < jets.size(); ++i) {
    bool is_overlapping = false;
    for (size_t j = 0; j < filtered_jets.size(); ++j) {
      double delta_r = deltaR(jets[i], filtered_jets[j]);      
      // Check if the jets are overlapping
      if (delta_r < min_delta_R) {
	is_overlapping = true;
	break; // No need to check further
      }
    }
    if (!is_overlapping) {
      filtered_jets.push_back(jets[i]); // Keep this jet if not overlapping
    }
  }
  return filtered_jets;
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
  std::map<double, long long> raw_jet_pt_counts;
  std::map<double, long long> sub_jet_pt_counts;

  std::map<double, long long> raw_jet_cut_pt_counts;
  std::map<double, long long> sub_jet_cut_pt_counts;
  
  std::map<double, long long> pion_pt_counts;
  std::map<double, long long> particle_pt_counts;
  
  std::map<double, long long> background_density_counts;
  std::map<int, int> particles_in_event_and_jets;
  
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
	//int particle_count = 0;
	//int pion_count = 0;
	for (int i = 0; i < n_hadrons; ++i) {
	  if (!std::getline(file, line)) {
	    std::cerr << "Unexpected end of file" << std::endl;
	    return 1;
	  }
	  Particle p;
	  if (sscanf(line.c_str(), "%d %d %d %lf %lf %lf %lf", &p.id, &p.pid, &p.status, &p.E, &p.Px, &p.Py, &p.Pz) == 7) {
	    particles.push_back(p);
	    double particle_pt = std::sqrt(p.Px * p.Px + p.Py * p.Py); // Calculate pT
	    double binned_particle_pt = std::round(particle_pt / PARTICLE_PT_BIN_SIZE) * PARTICLE_PT_BIN_SIZE; // Bin the pT
	    particle_pt_counts[binned_particle_pt]++;
	    
	    if (std::fabs(p.pid) == 211){
	      double pion_pt = std::sqrt(p.Px * p.Px + p.Py * p.Py); // Calculate pT
	      double binned_pion_pt = std::round(pion_pt / PION_PT_BIN_SIZE) * PION_PT_BIN_SIZE; // Bin the pT
	      pion_pt_counts[binned_pion_pt]++;
	    }
	  } else {
	    std::cerr << "Error parsing particle data" << std::endl;
	    return 1;
	  }
	  //std::cout << "Total particles: " << particle_count << ", pions: " << pion_count << std::endl;
	  
	
	}
      
      int num_particles_in_event = particles.size();
      std::vector<fastjet::PseudoJet> input_particles;
      for (const auto& p : particles) {
	fastjet::PseudoJet pj(p.Px, p.Py, p.Pz, p.E);
	pj.set_user_index(p.pid);
	input_particles.push_back(pj);	
      }
	//std::cout << "Total fastjet particles: " << fastjet_particle_count << ", fastjet pions: " << fastjet_pion_count << std::endl;
	
	fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, JET_RADIUS);
	
	fastjet::GhostedAreaSpec area_spec(ETA_MAX_GHOST);
	fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, area_spec);
	
	fastjet::ClusterSequenceArea cs(input_particles, jet_def, area_def);
	
	fastjet::Selector sel_jets = fastjet::SelectorAbsEtaMax(JET_ABS_ETA_MAX);
	//fastjet::Selector sel_jets = fastjet::SelectorAbsEtaMax(JET_ABS_ETA_MAX) && fastjet::SelectorPtMin(JET_PT_CUT);
	std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(sel_jets(cs.inclusive_jets()));
	double min_delta_R = 0.4;
	std::vector<fastjet::PseudoJet> non_overlapping_jets = filter_overlapping_jets(jets, min_delta_R);			
	//std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt((cs.inclusive_jets()));

	fastjet::JetDefinition jet_def_bkg(fastjet::kt_algorithm, R_BKG);
	//fastjet::Selector selector_bkg = fastjet::SelectorAbsEtaMax(ETA_MAX_GHOST) && fastjet::SelectorPtMin(JET_PT_CUT);
	fastjet::Selector selector_bkg = fastjet::SelectorAbsEtaMax(ETA_MAX_GHOST);
	fastjet::JetMedianBackgroundEstimator bge(selector_bkg, jet_def_bkg, area_def);
	bge.set_particles(input_particles);
	
	double rho = bge.rho(); // Get the background density
        
	fastjet::Subtractor subtractor(&bge);
	subtractor.set_use_rho_m(true);
	//std::vector<fastjet::PseudoJet> subtracted_jets = subtractor(jets);
	
	int num_particles_in_jets = 0;
	for (const auto& jet: non_overlapping_jets){
	  num_particles_in_jets += jet.constituents().size();
	}
	particles_in_event_and_jets[num_particles_in_event]+= num_particles_in_jets;

	for (const auto& jet : non_overlapping_jets){
	  double jet_area = jet.area(); // Accumulate the area for each jet
	  double adjusted_rho = rho * jet_area;
	  double binned_rho = std::round((adjusted_rho) / RHO_BIN_SIZE) * RHO_BIN_SIZE;
	  background_density_counts[binned_rho]++; // Count occurrences of this rho
	  
	  fastjet::PseudoJet subtracted_jet = subtractor(jet);
	  
	  //if (subtracted_jet.pt() < MIN_JET_PT || std::fabs(subtracted_jet.eta()) > JET_ABS_ETA_MAX || subtracted_jet.E() < 0.0) continue;
	  std::vector<fastjet::PseudoJet> raw_constituents = jet.constituents();
	  std::vector<fastjet::PseudoJet> sub_constituents = subtracted_jet.constituents();
	  if (raw_constituents.size() > 1 && sub_constituents.size() > 1){

	    double original_jet_pt = jet.pt();
	    double subtracted_jet_pt = subtracted_jet.pt();
	    if (original_jet_pt > JET_PT_CUT && subtracted_jet_pt > JET_PT_CUT){
	      double raw_jet_pt = original_jet_pt;
	      double binned_raw_jet_pt = std::round(raw_jet_pt / RAW_PT_BIN_SIZE) * RAW_PT_BIN_SIZE;
	      raw_jet_pt_counts[binned_raw_jet_pt]++;
	      
	      double sub_jet_pt = subtracted_jet_pt;
	      double binned_sub_jet_pt = std::round(sub_jet_pt / SUB_PT_BIN_SIZE) * SUB_PT_BIN_SIZE;
	      sub_jet_pt_counts[binned_sub_jet_pt]++;
	      
	      double pt_difference = original_jet_pt - subtracted_jet.pt();
	      double binned_pt_difference = std::round(pt_difference / PT_DIFF_BIN_SIZE) * PT_DIFF_BIN_SIZE;
	      jet_pt_difference_counts[binned_pt_difference]++;

	      //std::cout << "Subtracted jet pT: " << subtracted_jet.pt() << " and Raw jet PT: " << jet.pt() << std::endl;
	      //std::cout << "Background density (rho): " << bge.rho(jet) << std::endl;
	      //std::cout << "Jet area: " << jet.area() << std::endl;
	      
	      //double smearing_jet_factor = 1.0 + jet_distribution(generator);
	      //double smeared_jet_pt = original_jet_pt * smearing_jet_factor;
	      
	      
	      if (subtracted_jet_pt >= MIN_JET_PT && subtracted_jet_pt < MAX_JET_PT && original_jet_pt >= MIN_JET_PT && original_jet_pt < MAX_JET_PT){
		total_jets++;
		int pt_bin = std::min(static_cast<int>(subtracted_jet_pt / PT_BIN_WIDTH), NUM_PT_BINS - 1);
		dndpt_bins[pt_bin].jet_count++;
		
		double total_jet_momentum = jet.modp();	    
		double total_subtracted_jet_momentum = subtracted_jet.modp();
		
		double raw_jet_cut_pt = original_jet_pt;
		double binned_raw_jet_cut_pt = std::round(raw_jet_cut_pt / RAW_CUT_PT_BIN_SIZE) * RAW_CUT_PT_BIN_SIZE;
		raw_jet_cut_pt_counts[binned_raw_jet_cut_pt]++;
		
		double sub_jet_cut_pt = subtracted_jet_pt;
		double binned_sub_jet_cut_pt = std::round(sub_jet_cut_pt / SUB_CUT_PT_BIN_SIZE) * SUB_CUT_PT_BIN_SIZE;
		sub_jet_cut_pt_counts[binned_sub_jet_cut_pt]++;
		
		
		//int pt_bin = std::min(static_cast<int>(original_jet_pt / PT_BIN_WIDTH), NUM_PT_BINS - 1);
		for (const auto& constituent : sub_constituents) {
		  if (std::fabs(constituent.user_index()) == 211) { // PID 22 is for photons, abs(211) is for charged pions
		    double delta_r = deltaR(jet, constituent);
		    if (delta_r <= JET_RADIUS) {
		      double original_pion_pt = constituent.pt();
		      //double smearing_pion_factor = 1.0 + pion_distribution(generator);
		      //double smeared_pion_pt = original_pion_pt * smearing_pion_factor;
		      double z = original_pion_pt / subtracted_jet_pt;
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
	}
	processed_events++;
	if (processed_events % 1000 == 0) {
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

  std::cout << "Total processed events: " << processed_events << std::endl;
  std::cout << "Total jets: " << total_jets << std::endl;
  std::cout << "Total pions: " << total_pions << std::endl;

  //std::ofstream particle_pt_output_file("particle_AuAu_pTHat15_100_out_alldecayoff_noprompt_hydro_iSS10_10k_pt_counts.txt");
  std::ofstream particle_pt_output_file("pion_AuAu_out_only_hydro_iSS10_10k_pt_counts.txt");
  if (!particle_pt_output_file.is_open()) {
    std::cerr << "Error opening particle pT counts output file" << std::endl;
    return 1;
  }

    //particle_pt_output_file << "Particle pT (GeV/c) and counts: " << std::endl;
  for (const auto& entry : particle_pt_counts) {
    particle_pt_output_file << entry.first << " " << entry.second << std::endl; // x = binned pT, y = count
  }

    // Close the particle pT output file
  particle_pt_output_file.close();
  
  //std::ofstream pion_pt_output_file("pion_AuAu_pTHat15_100_out_alldecayoff_noprompt_hydro_iSS10_10k_pt_counts.txt");
  std::ofstream pion_pt_output_file("pion_AuAu_out_only_hydro_iSS10_10k_pt_counts.txt");
  if (!pion_pt_output_file.is_open()) {
    std::cerr << "Error opening particle pT counts output file" << std::endl;
    return 1;
  }

    //particle_pt_output_file << "Particle pT (GeV/c) and counts: " << std::endl;
  for (const auto& entry : pion_pt_counts) {
    pion_pt_output_file << entry.first << " " << entry.second << std::endl; // x = binned pT, y = count
  }

    // Close the particle pT output file
  pion_pt_output_file.close();
  
  //std::ofstream jet_pt_diff_output_file("pion_AuAu_pTHat15_100_out_alldecayoff_noprompt_hydro_iSS10_10k_jet_pt_differences.txt");
  std::ofstream jet_pt_diff_output_file("pion_AuAu_out_only_hydro_iSS10_10k_jet_pt_differences.txt");
  if (!jet_pt_diff_output_file.is_open()){
      std::cerr << "Error opening output file" << std::endl;
      return 1;
    }
    //output_file << "Jet pT differences (pT - subtracted pT) and counts: " << std::endl;
  for (const auto& entry : jet_pt_difference_counts){
    jet_pt_diff_output_file << entry.first << " " << entry.second << std::endl;
  }
  jet_pt_diff_output_file.close();
  

  //std::ofstream raw_jet_pt_output_file("pion_AuAu_pTHat15_100_out_alldecayoff_noprompt_hydro_iSS10_10k_raw_jet_pt_counts.txt");
  std::ofstream raw_jet_pt_output_file("pion_AuAu_out_only_hydro_iSS10_10k_raw_jet_pt_counts.txt");
  if (!raw_jet_pt_output_file.is_open()){
      std::cerr << "Error opening output file" << std::endl;
      return 1;
    }
    //output_file << "Jet pT differences (pT - subtracted pT) and counts: " << std::endl;
  for (const auto& entry : raw_jet_pt_counts){
    raw_jet_pt_output_file << entry.first << " " << entry.second << std::endl;
  }
  raw_jet_pt_output_file.close();

  //std::ofstream raw_jet_cut_pt_output_file("pion_AuAu_pTHat15_100_out_alldecayoff_noprompt_hydro_iSS10_10k_raw_jet_cut_pt_counts.txt");
  std::ofstream raw_jet_cut_pt_output_file("pion_AuAu_out_only_hydro_iSS10_10k_raw_jet_cut_pt_counts.txt");
    if (!raw_jet_cut_pt_output_file.is_open()){
      std::cerr << "Error opening output file" << std::endl;
      return 1;
    }
    //output_file << "Jet pT differences (pT - subtracted pT) and counts: " << std::endl;
  for (const auto& entry : raw_jet_cut_pt_counts){
    raw_jet_cut_pt_output_file << entry.first << " " << entry.second << std::endl;
  }
  raw_jet_cut_pt_output_file.close();

  //std::ofstream sub_jet_pt_output_file("pion_AuAu_pTHat15_100_out_alldecayoff_noprompt_hydro_iSS10_10k_sub_jet_pt_counts.txt");
  std::ofstream sub_jet_pt_output_file("pion_AuAu_out_only_hydro_iSS10_10k_sub_jet_pt_counts.txt");
  if (!sub_jet_pt_output_file.is_open()){
      std::cerr << "Error opening output file" << std::endl;
      return 1;
    }
    //output_file << "Jet pT differences (pT - subtracted pT) and counts: " << std::endl;
  for (const auto& entry : sub_jet_pt_counts){
    sub_jet_pt_output_file << entry.first << " " << entry.second << std::endl;
  }
  sub_jet_pt_output_file.close();

  //std::ofstream sub_jet_cut_pt_output_file("pion_AuAu_pTHat15_100_out_alldecayoff_noprompt_hydro_iSS10_10k_sub_jet_cut_pt_counts.txt");
  std::ofstream sub_jet_cut_pt_output_file("pion_AuAu_out_only_hydro_iSS10_10k_sub_jet_cut_pt_counts.txt");
    if (!sub_jet_cut_pt_output_file.is_open()){
      std::cerr << "Error opening output file" << std::endl;
      return 1;
    }
    //output_file << "Jet pT differences (pT - subtracted pT) and counts: " << std::endl;
  for (const auto& entry : sub_jet_cut_pt_counts){
    sub_jet_cut_pt_output_file << entry.first << " " << entry.second << std::endl;
  }
  sub_jet_cut_pt_output_file.close();

  
  //std::ofstream background_output_file("pion_AuAu_pTHat15_100_out_alldecayoff_noprompt_hydro_iSS10_10k_background_density_area_counts.txt");
  std::ofstream background_output_file("pion_AuAu_out_only_hydro_iSS10_10k_background_density_area_counts.txt");
  if (!background_output_file.is_open()) {
    std::cerr << "Error opening background density output file" << std::endl;
    return 1;
  }

  //background_output_file << "Background density (rho) and counts: " << std::endl;
  for (const auto& entry : background_density_counts) {
    background_output_file << entry.first << " " << entry.second << std::endl; // x = rho, y = count
  }

  // Close the background density output file
  background_output_file.close();
  /*
  std::ofstream particles_output_file("pion_pp_pTHat15_100_out_alldecayoff_noprompt_particles_in_event_and_jets.txt");
  if (!particles_output_file.is_open()) {
    std::cerr << "Error opening particles output file" << std::endl;
    return 1;
  }

  //particles_output_file << "Number of particles in event, Number of particles in jets: " << std::endl;
  for (const auto& entry : particles_in_event_and_jets) {
    particles_output_file << entry.first << " " << entry.second << std::endl; // x = number of particles in event, y = number of particles in jets
  }

  // Close the particles output file
  particles_output_file.close();
  */
 
  
  
  return 0;
}
