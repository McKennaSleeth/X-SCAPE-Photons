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

const double PION_PT_BIN_SIZE = 1.0;
const double PARTICLE_PT_BIN_SIZE = 1.0;

const double RAW_CUT_PT_BIN_SIZE = 0.5;
const double SUB_CUT_PT_BIN_SIZE = 0.5;

const double RHO_BIN_SIZE = 1.0;

const double MIN_JET_PT = 20.0;
const double MAX_JET_PT = 30.0;

const int NUM_PT_BINS = 20;
const double PT_BIN_WIDTH = 100 / NUM_PT_BINS;

const double JET_RADIUS = 0.4;

const double JET_ABS_ETA_MAX = 1.1 - JET_RADIUS;
const double JET_PT_CUT      = 10.0;
const double R_BKG           = 0.4;
const double ETA_MAX_GHOST   = 1.1;

//const double CENTRALITY_0_1_THRESHOLD = 0.2;

long long total_raw_jets = 0;
long long total_subtracted_jets = 0;
long long total_raw_pions = 0;
long long total_subtracted_pions = 0;
long long total_pions = 0;
long long total_cones = 0;
long long total_cone_pions = 0;
int n_random_cones = 100;

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

std::vector<DNDZBin> raw_dndz_bins(NUM_BINS);
std::vector<DNDDeltaRBin> raw_dnd_deltaR_bins(NUM_DELTA_R_BINS);
std::vector<DNDPTBin> raw_dndpt_bins(NUM_PT_BINS);

std::vector<DNDZBin> sub_dndz_bins(NUM_BINS);
std::vector<DNDDeltaRBin> sub_dnd_deltaR_bins(NUM_DELTA_R_BINS);
std::vector<DNDPTBin> sub_dndpt_bins(NUM_PT_BINS);

std::vector<DNDZBin> cone_dndz_bins(NUM_BINS);
std::vector<DNDDeltaRBin> cone_dnd_deltaR_bins(NUM_DELTA_R_BINS);
std::vector<DNDPTBin> cone_dndpt_bins(NUM_PT_BINS);

double deltaR(const fastjet::PseudoJet& p1, const fastjet::PseudoJet& p2) {
  double deta = p1.eta() - p2.eta();
  double dphi = std::abs(p1.phi() - p2.phi());
  if (dphi > M_PI) dphi = 2 * M_PI - dphi;
  return std::sqrt(deta * deta + dphi * dphi);
}


//unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//std::default_random_engine generator(seed);
//std::normal_distribution<double> jet_distribution(0.0, 0.15);
//std::normal_distribution<double> pion_distribution(0.0, 0.10);

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
//std::uniform_real_distribution<double> pt_distribution(0.0, 100.0); // Random pT between 20 and 30 GeV
std::uniform_real_distribution<double> phi_distribution(-M_PI, M_PI);
std::uniform_real_distribution<double> eta_distribution(-1.1 + JET_RADIUS, 1.1 - JET_RADIUS);



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
    raw_dndz_bins[i].z_center = i * BIN_WIDTH + BIN_WIDTH / 2;
    raw_dndz_bins[i].pion_count = 0;
  }

  for (int i = 0; i < NUM_BINS; ++i) {
    sub_dndz_bins[i].z_center = i * BIN_WIDTH + BIN_WIDTH / 2;
    sub_dndz_bins[i].pion_count = 0;
  }

  for (int i = 0; i < NUM_BINS; ++i) {
    cone_dndz_bins[i].z_center = i * BIN_WIDTH + BIN_WIDTH / 2;
    cone_dndz_bins[i].pion_count = 0;
  }
    
  for (int i = 0; i < NUM_DELTA_R_BINS; ++i) {
    raw_dnd_deltaR_bins[i].deltaR_center = i * DELTA_R_BIN_WIDTH + DELTA_R_BIN_WIDTH / 2;
    raw_dnd_deltaR_bins[i].pion_count = 0;
  }

  for (int i = 0; i < NUM_DELTA_R_BINS; ++i) {
    sub_dnd_deltaR_bins[i].deltaR_center = i * DELTA_R_BIN_WIDTH + DELTA_R_BIN_WIDTH / 2;
    sub_dnd_deltaR_bins[i].pion_count = 0;
  }

  for (int i = 0; i < NUM_DELTA_R_BINS; ++i) {
    cone_dnd_deltaR_bins[i].deltaR_center = i * DELTA_R_BIN_WIDTH + DELTA_R_BIN_WIDTH / 2;
    cone_dnd_deltaR_bins[i].pion_count = 0;
  }
    
  for (int i = 0; i < NUM_PT_BINS; ++i) {
    raw_dndpt_bins[i].pT_center = i * PT_BIN_WIDTH + PT_BIN_WIDTH / 2;
    raw_dndpt_bins[i].jet_count = 0;
  }

  for (int i = 0; i < NUM_PT_BINS; ++i) {
    sub_dndpt_bins[i].pT_center = i * PT_BIN_WIDTH + PT_BIN_WIDTH / 2;
    sub_dndpt_bins[i].jet_count = 0;
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
  //std::map<double, long long> rhoA_vs_num_particles;
  std::map<double, long long> rho_vs_num_particles;
  //std::vector<std::tuple<double, int, long long>> rhoA_particle_event_counts; // (rho * A, number of particles, event index)
  std::vector<std::tuple<double, int, long long>> rho_particle_event_counts; // (rho, number of particles, event index)
  std::vector<int> particle_counts;

  std::ofstream subtracted_jet_etaphipt_output_file("pion_AuAu_pTHat15_100_hard_sector_hydro10_centr0_10_10k_events_subtracted_jet_trackpTcut_etaphipt_output.txt", std::ios::app);
  if (!subtracted_jet_etaphipt_output_file.is_open()) {
    std::cerr << "Error opening subtracted jet output file" << std::endl;
    return 1;
  }
  
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
    //std::vector<Particle> centr_0_1_particles;
    double EPangle;
    
    while (std::getline(file, line)) {
      if (sscanf(line.c_str(), "#\tEvent\t%d\tweight\t%*f\tEPangle\t%lf\tN_hadrons\t%d", &event_number, &EPangle, &n_hadrons) == 3) {
	//if (EPangle != 0){
	//  std::cout << "Event Number: " << event_number << ", EPangle: " << EPangle << ", Number of Hadrons: " << n_hadrons << std::endl;
	//}
	particles.clear();
	particles.reserve(n_hadrons);
        
	for (int i = 0; i < n_hadrons; ++i) {
	  if (!std::getline(file, line)) {
	    std::cerr << "Unexpected end of file" << std::endl;
	    return 1;
	  }
	  Particle p;
	  if (sscanf(line.c_str(), "%d %d %d %lf %lf %lf %lf", &p.id, &p.pid, &p.status, &p.E, &p.Px, &p.Py, &p.Pz) == 7) {
	    double particle_pt = std::sqrt(p.Px * p.Px + p.Py * p.Py); // Calculate pT
	    if (particle_pt > 0.5){
	      particles.push_back(p);
	    } //minimum track pT cut
	    double binned_particle_pt = std::round(particle_pt / PARTICLE_PT_BIN_SIZE) * PARTICLE_PT_BIN_SIZE; // Bin the pT
	    particle_pt_counts[binned_particle_pt]++;
	    
	    if (std::fabs(p.pid) == 211){
	      double pion_pt = std::sqrt(p.Px * p.Px + p.Py * p.Py); // Calculate pT
	      double binned_pion_pt = std::round(pion_pt / PION_PT_BIN_SIZE) * PION_PT_BIN_SIZE; // Bin the pT
	      total_pions++;
	      pion_pt_counts[binned_pion_pt]++;
	    }
	  } else { 
	    std::cerr << "Error parsing particle data" << std::endl;
	    return 1;
	    
	  } // scanning particle information and assigning variables condition
	} //hadron for loop
	 
	int num_particles_per_event = particles.size();
	
	double min_pt = std::numeric_limits<double>::max();
	double max_pt = std::numeric_limits<double>::lowest();
	  
	std::vector<fastjet::PseudoJet> input_particles;
	for (const auto& p : particles) {
	  fastjet::PseudoJet pj(p.Px, p.Py, p.Pz, p.E);
	  pj.set_user_index(p.pid);
	  input_particles.push_back(pj);
	  double pt = pj.pt();
	  if (pt < min_pt) min_pt = pt;
	  if (pt > max_pt) max_pt = pt;
	}
	  
	fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, JET_RADIUS);
	fastjet::GhostedAreaSpec area_spec(ETA_MAX_GHOST);
	fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, area_spec);
	//pick a random jet before clustering 
	fastjet::ClusterSequenceArea cs(input_particles, jet_def, area_def);
	
	//fastjet::Selector sel_jets = fastjet::SelectorAbsEtaMax(JET_ABS_ETA_MAX);
	fastjet::Selector sel_jets = fastjet::SelectorAbsEtaMax(JET_ABS_ETA_MAX) && fastjet::SelectorPtMin(JET_PT_CUT);
	std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(sel_jets(cs.inclusive_jets()));
	  
	//std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt((cs.inclusive_jets()));
	  
	fastjet::JetDefinition jet_def_bkg(fastjet::kt_algorithm, R_BKG);
	fastjet::Selector selector_bkg = fastjet::SelectorAbsEtaMax(ETA_MAX_GHOST) && fastjet::SelectorPtMin(JET_PT_CUT);
	fastjet::JetMedianBackgroundEstimator bge(selector_bkg, jet_def_bkg, area_def);
	bge.set_particles(input_particles);
	
	double rho = bge.rho(); // Get the background density
	
	
	fastjet::Subtractor subtractor(&bge);
	subtractor.set_use_rho_m(true);
	//std::vector<fastjet::PseudoJet> subtracted_jets = subtractor(jets);
	
	int num_particles_in_jets = 0;
	for (const auto& jet: jets){
	  num_particles_in_jets += jet.constituents().size();
	}
	particles_in_event_and_jets[num_particles_per_event]+= num_particles_in_jets;
	
	  
	//double total_rho = 0.0;
	//int jet_number = 0;
	//random cone loop before jet loop?
	for (const auto& jet : jets){
	  //double jet_area = jet.area(); // Accumulate the area for each jet
	  //double rhoA = rho * jet_area;
	  //total_rho += rho;
	  //jet_number++;
	  //double adjusted_rho = rho * jet_area;
	  double binned_rho = std::round((rho) / RHO_BIN_SIZE) * RHO_BIN_SIZE;
	  background_density_counts[binned_rho]++; // Count occurrences of this rho
	  //rhoA_vs_num_particles[binned_rho] += num_particles_in_event;
	  rho_particle_event_counts.emplace_back(rho, num_particles_per_event, processed_events);
	  //total_rhoA += binned_rho;
	  //jet_number;
	  //rhoA_vs_num_particles.emplace_back(binned_rho, num_particles_in_event);
	  fastjet::PseudoJet subtracted_jet = subtractor(jet);
	  
	  //if (subtracted_jet.pt() < MIN_JET_PT || std::fabs(subtracted_jet.eta()) > JET_ABS_ETA_MAX || subtracted_jet.E() < 0.0) continue;
	  std::vector<fastjet::PseudoJet> raw_constituents = jet.constituents();
	  std::vector<fastjet::PseudoJet> sub_constituents = subtracted_jet.constituents();
	  if (raw_constituents.size() > 1 && sub_constituents.size() > 1){
	    //total_jets++;
	    
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
	      
	      if (original_jet_pt >= MIN_JET_PT && original_jet_pt < MAX_JET_PT){
		total_raw_jets++;
		int raw_pt_bin = std::min(static_cast<int>(original_jet_pt / PT_BIN_WIDTH), NUM_PT_BINS - 1);
		//int pt_bin = std::min(static_cast<int>(original_jet_pt / PT_BIN_WIDTH), NUM_PT_BINS - 1);
		raw_dndpt_bins[raw_pt_bin].jet_count++;
		
		double total_jet_momentum = jet.modp();	    
		
		double raw_jet_cut_pt = original_jet_pt;
		double binned_raw_jet_cut_pt = std::round(raw_jet_cut_pt / RAW_CUT_PT_BIN_SIZE) * RAW_CUT_PT_BIN_SIZE;
		raw_jet_cut_pt_counts[binned_raw_jet_cut_pt]++;
		
		//int pt_bin = std::min(static_cast<int>(original_jet_pt / PT_BIN_WIDTH), NUM_PT_BINS - 1);
		for (const auto& constituent : raw_constituents) {
		  //for (const auto& constituent : raw_constituents) {
		  if (std::fabs(constituent.user_index()) == 211) { // PID 22 is for photons, abs(211) is for charged pions
		    double delta_r = deltaR(jet, constituent);
		    if (delta_r <= JET_RADIUS) {
		      double original_pion_pt = constituent.pt();
		      total_raw_pions++;
		      //double smearing_pion_factor = 1.0 + pion_distribution(generator);
		      //double smeared_pion_pt = original_pion_pt * smearing_pion_factor;
		      double z = original_pion_pt / original_jet_pt;
		      //double z = original_pion_pt / original_jet_pt;
		      int raw_z_bin = std::min(static_cast<int>(z / BIN_WIDTH), NUM_BINS - 1);
		      int raw_delta_r_bin = std::min(static_cast<int>(delta_r / DELTA_R_BIN_WIDTH), NUM_DELTA_R_BINS - 1);
		      raw_dndz_bins[raw_z_bin].pion_count++;
		      raw_dnd_deltaR_bins[raw_delta_r_bin].pion_count++;
		      //total_pions++;
		    } //deltaR condition
		  } //pion selection condition
		} //raw jet constituent loop
	      } //raw jet pT selection condition

	      //std::cout << "Subtracted jet pT: " << subtracted_jet.pt() << " and Raw jet PT: " << jet.pt() << std::endl;
	      //std::cout << "Background density (rho): " << bge.rho(jet) << std::endl;
	      //std::cout << "Jet area: " << jet.area() << std::endl;
	      
	      //double smearing_jet_factor = 1.0 + jet_distribution(generator);
	      //double smeared_jet_pt = original_jet_pt * smearing_jet_factor;
	      
	      
	      if (subtracted_jet_pt >= MIN_JET_PT && subtracted_jet_pt < MAX_JET_PT){
		//if (original_jet_pt >= MIN_JET_PT && original_jet_pt <=MAX_JET_PT){
		//std::cout << "Subtracted jet pT: " << subtracted_jet_pt << "Original jet pT: " << original_jet_pt << std::endl;
		double jet_eta = subtracted_jet.eta();
		double jet_phi = subtracted_jet.phi();
		double jet_pt = subtracted_jet.pt();
		double jet_event = processed_events;
		
		subtracted_jet_etaphipt_output_file << jet_event << " " << jet_eta << " " << jet_phi << " " << jet_pt << std::endl;
		
		total_subtracted_jets++;
		/*
		  double random_cone_eta = eta_distribution(generator);
		  double random_cone_phi = phi_distribution(generator);
		  std::uniform_real_distribution<double> pt_distribution(min_pt, max_pt);
		  double random_cone_pt = pt_distribution(generator);
		  
		  double E_cone = std::sqrt(random_cone_pt * random_cone_pt * std::cosh(random_cone_eta) * std::cosh(random_cone_eta));
		  double Px_cone = random_cone_pt * std::cos(random_cone_phi);
		  double Py_cone = random_cone_pt * std::sin(random_cone_phi);
		  double Pz_cone = random_cone_pt * std::sinh(random_cone_eta);
		  
		  fastjet::PseudoJet cone_center(Px_cone, Py_cone, Pz_cone, E_cone); // Create the cone center
		  
		  //std::cout << "Random cone pT (before selection): " << cone_center.pt() << std::endl;
		  total_cones++;
		  
		  //fastjet::Selector cone_selector = fastjet::SelectorCircle(JET_RADIUS);
		  
		  double cone_axis_x = std::cos(random_cone_phi) * std::cosh(random_cone_eta);
		  double cone_axis_y = std::sin(random_cone_phi) * std::cosh(random_cone_eta);
		  double cone_axis_z = std::sinh(random_cone_eta);
		  
		  // Normalize the cone axis
		  double cone_axis_length = std::sqrt(cone_axis_x * cone_axis_x + cone_axis_y * cone_axis_y + cone_axis_z * cone_axis_z);
		  cone_axis_x /= cone_axis_length;
		  cone_axis_y /= cone_axis_length;
		  cone_axis_z /= cone_axis_length;
		  
		  // Print the cone pT and axis
		  //std::cout << "Random Cone pT: " << cone_center.pt() << std::endl;
		  //std::cout << "Random Cone Axis: (" << cone_axis_x << ", " << cone_axis_y << ", " << cone_axis_z << ")" << std::endl;
		  
		  std::vector<fastjet::PseudoJet> cone_constituents;
		  for (const auto& particle : input_particles) {
		    // Check if the particle is within the cone
		    
		    if (deltaR(cone_center, particle) <= JET_RADIUS) {
		      //std::cout << "DeltaR of cone and particle: " << deltaR(cone_center, particle) << std::endl;
		      cone_constituents.push_back(particle);
		    }
		  }
		  
		  if (!cone_constituents.empty()){
		    for (const auto& constituent : cone_constituents) {
		      if (cone_constituents.size() > 1){
			double random_cone_pt = cone_center.pt();
			//if (random_cone_pt >= MIN_JET_PT && random_cone_pt < MAX_JET_PT){
			//std::cout << "Random Cone pT: " << random_cone_pt << std::endl; 
			if (std::fabs(constituent.user_index()) == 211) { // PID 22 is for photons, abs(211) is for charged pions
			  double delta_r_cone = deltaR(cone_center, constituent);
			  //std::cout << "DeltaR of cone and pion: " << deltaR(cone_center, constituent) << std::endl;
			  double cone_pion_pt = constituent.pt();
			  //std::cout << "Cone Pion pT: " << cone_pion_pt << std::endl;
			  double cone_z = cone_pion_pt / random_cone_pt;
			  //double z = original_pion_pt / original_jet_pt;
			  int cone_z_bin = std::min(static_cast<int>(cone_z / BIN_WIDTH), NUM_BINS - 1);
			  int cone_delta_r_bin = std::min(static_cast<int>(delta_r_cone / DELTA_R_BIN_WIDTH), NUM_DELTA_R_BINS - 1);
			  cone_dndz_bins[cone_z_bin].pion_count++;
			  total_cone_pions++;
			  cone_dnd_deltaR_bins[cone_delta_r_bin].pion_count++;
			} //pion selection condition
		      } //cone constituent condition 
		    } //cone constituent for loop
		  } //cone pT selection condition
	 
		  */
		int sub_pt_bin = std::min(static_cast<int>(subtracted_jet_pt / PT_BIN_WIDTH), NUM_PT_BINS - 1);
		//int pt_bin = std::min(static_cast<int>(original_jet_pt / PT_BIN_WIDTH), NUM_PT_BINS - 1);
		sub_dndpt_bins[sub_pt_bin].jet_count++;
		
		
		double total_subtracted_jet_momentum = subtracted_jet.modp();
		
		double sub_jet_cut_pt = subtracted_jet_pt;
		double binned_sub_jet_cut_pt = std::round(sub_jet_cut_pt / SUB_CUT_PT_BIN_SIZE) * SUB_CUT_PT_BIN_SIZE;
		sub_jet_cut_pt_counts[binned_sub_jet_cut_pt]++;
		
		
		
		//int pt_bin = std::min(static_cast<int>(original_jet_pt / PT_BIN_WIDTH), NUM_PT_BINS - 1);
		for (const auto& constituent : sub_constituents) {
		  //for (const auto& constituent : raw_constituents) {
		  if (std::fabs(constituent.user_index()) == 211) { // PID 22 is for photons, abs(211) is for charged pions
		    double delta_r = deltaR(subtracted_jet, constituent);
		    if (delta_r <= JET_RADIUS) {
		      double original_pion_pt = constituent.pt();
		      total_subtracted_pions++;
		      //double smearing_pion_factor = 1.0 + pion_distribution(generator);
		      //double smeared_pion_pt = original_pion_pt * smearing_pion_factor;
		      double z = original_pion_pt / subtracted_jet_pt;
		      //double z = original_pion_pt / original_jet_pt;
		      int sub_z_bin = std::min(static_cast<int>(z / BIN_WIDTH), NUM_BINS - 1);
		      int sub_delta_r_bin = std::min(static_cast<int>(delta_r / DELTA_R_BIN_WIDTH), NUM_DELTA_R_BINS - 1);
		      sub_dndz_bins[sub_z_bin].pion_count++;
		      sub_dnd_deltaR_bins[sub_delta_r_bin].pion_count++;
		      //total_pions++;
		    } //deltaR condition
		  } //pion selection condition
		} //subtracted jet constituent loop
	      } //subtracted jet pT selection condition
	    } //raw and subtracted jet pT cut condition
	  } //raw and subtracted jet constituent cut condition
	} //jet for loop

	  
	processed_events++;
	if (processed_events % 1000 == 0) {
	  std::cout << "Processed " << processed_events << " events in " << filename << std::endl;
	} //processed events condition
	
      } //event condition
    } //filename while loop
    
    file.close();
    std::cout << "Finished processing " << filename << std::endl;
  } //filename loop

  subtracted_jet_etaphipt_output_file.close();
	  
  std::cout << "Final Raw Jet dN/dz distribution: " << std::endl;
  for (const auto& bin : raw_dndz_bins) {
    double dndz = static_cast<double>(bin.pion_count) / (total_raw_jets * BIN_WIDTH);
    double error = std::sqrt(bin.pion_count) / (total_raw_jets * BIN_WIDTH);
    std::cout << bin.z_center << " " << dndz << " " << error << std::endl;
  }
  
  std::cout << "Final Raw dN/ddeltaR distribution: " << std::endl;
  for (const auto& bin : raw_dnd_deltaR_bins) {
    double dnd_deltaR = static_cast<double>(bin.pion_count) / (total_raw_jets * DELTA_R_BIN_WIDTH);
    double error = std::sqrt(bin.pion_count) / (total_raw_jets * DELTA_R_BIN_WIDTH);
    std::cout << bin.deltaR_center << " " << dnd_deltaR << " " << error << std::endl;
  }
  
  std::cout << "Final Raw dN/dpT distribution: " << std::endl;
  for (const auto& bin : raw_dndpt_bins) {
    double dndpt = static_cast<double>(bin.jet_count) / (processed_events * PT_BIN_WIDTH);
    double error = std::sqrt(bin.jet_count) / (processed_events * PT_BIN_WIDTH);
    std::cout << bin.pT_center << " " << dndpt << " " << error << std::endl;
  }

  std::cout << "Final Sub Jet dN/dz distribution: " << std::endl;
  for (const auto& bin : sub_dndz_bins) {
    double dndz = static_cast<double>(bin.pion_count) / (total_subtracted_jets * BIN_WIDTH);
    double error = std::sqrt(bin.pion_count) / (total_subtracted_jets * BIN_WIDTH);
    std::cout << bin.z_center << " " << dndz << " " << error << std::endl;
  }

  std::cout << "Final Sub dN/ddeltaR distribution: " << std::endl;
  for (const auto& bin : sub_dnd_deltaR_bins) {
    double dnd_deltaR = static_cast<double>(bin.pion_count) / (total_subtracted_jets * DELTA_R_BIN_WIDTH);
    double error = std::sqrt(bin.pion_count) / (total_subtracted_jets * DELTA_R_BIN_WIDTH);
    std::cout << bin.deltaR_center << " " << dnd_deltaR << " " << error << std::endl;
  }

  std::cout << "Final Sub dN/dpT distribution: " << std::endl;
  for (const auto& bin : sub_dndpt_bins) {
    double dndpt = static_cast<double>(bin.jet_count) / (processed_events * PT_BIN_WIDTH);
    double error = std::sqrt(bin.jet_count) / (processed_events * PT_BIN_WIDTH);
    std::cout << bin.pT_center << " " << dndpt << " " << error << std::endl;
  }
  /*
  std::cout << "Final Random Cone dN/dz distribution: " << std::endl;
  for (const auto& bin : cone_dndz_bins) {
    double cone_dndz = static_cast<double>(bin.pion_count) / (total_cones * BIN_WIDTH);
    double cone_error = std::sqrt(bin.pion_count) / (total_cones * BIN_WIDTH);
    std::cout << bin.z_center << " " << cone_dndz << " " << cone_error << std::endl;
  }

  std::cout << "Final Random Cone dN/ddeltaR distribution: " << std::endl;
  for (const auto& bin : cone_dnd_deltaR_bins) {
    double cone_dnd_deltaR = static_cast<double>(bin.pion_count) / (total_cones * DELTA_R_BIN_WIDTH);
    double cone_error = std::sqrt(bin.pion_count) / (total_cones * DELTA_R_BIN_WIDTH);
    std::cout << bin.deltaR_center << " " << cone_dnd_deltaR << " " << cone_error << std::endl;
  }
  */
  std::cout << "Total processed events: " << processed_events << std::endl;
  //std::cout << "Total cones: " << total_cones << std::endl;
  //std::cout << "Total cone pions: " << total_cone_pions << std::endl;
  std::cout << "Total raw jets: " << total_raw_jets << std::endl;
  std::cout << "Total subtracted jets: " << total_subtracted_jets << std::endl;
  std::cout << "Total pions: " << total_pions << std::endl;
  std::cout << "Total raw pions: " << total_raw_pions << std::endl;
  std::cout << "Total subtracted pions: " << total_subtracted_pions << std::endl;
  

  
  std::ofstream particle_pt_output_file("particle_JETSCAPE_AuAu_pTHat15_100_hard_sector_hydro10_centr0_10_10k_trackpTcut_pt_counts.txt");
  //std::ofstream particle_pt_output_file("particle_JETSCAPE_AuAu_out_only_hydro100_iSS1_1M_pt_counts.txt");
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
  
  
  std::ofstream pion_pt_output_file("pion_JETSCAPE_AuAu_pTHat15_100_hard_sector_hydro10_centr0_10_10k_trackpTcut_pt_counts.txt");
  //std::ofstream pion_pt_output_file("pion_JETSCAPE_AuAu_out_only_hydro100_iSS1_1M_pt_counts.txt");
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
  
  
  std::ofstream jet_pt_diff_output_file("pion_JETSCAPE_AuAu_pTHat15_100_hard_sector_hydro10_centr0_10_10k_trackpTcut_jet_pt_differences.txt");
  //std::ofstream jet_pt_diff_output_file("pion_JETSCAPE_AuAu_out_only_hydro100_iSS1_1M_jet_pt_differences.txt");
  if (!jet_pt_diff_output_file.is_open()){
      std::cerr << "Error opening output file" << std::endl;
      return 1;
    }
    //output_file << "Jet pT differences (pT - subtracted pT) and counts: " << std::endl;
  for (const auto& entry : jet_pt_difference_counts){
    jet_pt_diff_output_file << entry.first << " " << entry.second << std::endl;
  }
  jet_pt_diff_output_file.close();
  
  
  std::ofstream raw_jet_pt_output_file("pion_JETSCAPE_AuAu_pTHat15_100_hard_sector_hydro10_centr0_10_10k_trackpTcut_raw_jet_pt_counts.txt");
  //std::ofstream raw_jet_pt_output_file("pion_JETSCAPE_AuAu_out_only_hydro100_iSS1_1M_raw_jet_pt_counts.txt");
  if (!raw_jet_pt_output_file.is_open()){
      std::cerr << "Error opening output file" << std::endl;
      return 1;
    }
    //output_file << "Jet pT differences (pT - subtracted pT) and counts: " << std::endl;
  for (const auto& entry : raw_jet_pt_counts){
    raw_jet_pt_output_file << entry.first << " " << entry.second << std::endl;
  }
  raw_jet_pt_output_file.close();
  
  
  std::ofstream raw_jet_cut_pt_output_file("pion_JETSCAPE_AuAu_pTHat15_100_hard_sector_hydro10_centr0_10_10k_trackpTcut_raw_jet_cut_pt_counts.txt");
  //std::ofstream raw_jet_cut_pt_output_file("pion_JETSCAPE_AuAu_out_only_hydro100_iSS1_1M_raw_jet_cut_pt_counts.txt");
    if (!raw_jet_cut_pt_output_file.is_open()){
      std::cerr << "Error opening output file" << std::endl;
      return 1;
    }
    //output_file << "Jet pT differences (pT - subtracted pT) and counts: " << std::endl;
  for (const auto& entry : raw_jet_cut_pt_counts){
    raw_jet_cut_pt_output_file << entry.first << " " << entry.second << std::endl;
  }
  raw_jet_cut_pt_output_file.close();
  
  
  std::ofstream sub_jet_pt_output_file("pion_JETSCAPE_AuAu_pTHat15_100_hard_sector_hydro10_centr0_10_10k_trackpTcut_sub_jet_pt_counts.txt");
  //std::ofstream sub_jet_pt_output_file("pion_JETSCAPE_AuAu_out_only_hydro100_iSS1_1M_sub_jet_pt_counts.txt");
  if (!sub_jet_pt_output_file.is_open()){
      std::cerr << "Error opening output file" << std::endl;
      return 1;
    }
    //output_file << "Jet pT differences (pT - subtracted pT) and counts: " << std::endl;
  for (const auto& entry : sub_jet_pt_counts){
    sub_jet_pt_output_file << entry.first << " " << entry.second << std::endl;
  }
  sub_jet_pt_output_file.close();

  
  std::ofstream sub_jet_cut_pt_output_file("pion_JETSCAPE_AuAu_pTHat15_100_hard_sector_hydro100_centr0_10_trackpTcut_1k_sub_jet_cut_pt_counts.txt");
  //std::ofstream sub_jet_cut_pt_output_file("pion_JETSCAPE_AuAu_out_only_hydro100_iSS1_1M_sub_jet_cut_pt_counts.txt");
    if (!sub_jet_cut_pt_output_file.is_open()){
      std::cerr << "Error opening output file" << std::endl;
      return 1;
    }
    //output_file << "Jet pT differences (pT - subtracted pT) and counts: " << std::endl;
  for (const auto& entry : sub_jet_cut_pt_counts){
    sub_jet_cut_pt_output_file << entry.first << " " << entry.second << std::endl;
  }
  sub_jet_cut_pt_output_file.close();
  
  
  std::ofstream background_output_file("pion_JETSCAPE_AuAu_pTHat15_100_hard_sector_hydro10_centr0_10_10k_trackpTcut_background_density_area_counts.txt");
  //std::ofstream background_output_file("pion_JETSCAPE_AuAu_out_only_hydro100_iSS1_1M_background_density_counts.txt");
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
  
  
  std::ofstream rho_particle_output_file("particle_AuAu_pTHat15_100_hard_sector_hydro10_centr0_10_10k_trackpTcut_rho_particle_counts.txt");
  //std::ofstream rho_particle_output_file("pion_JETSCAPE_AuAu_out_only_hydro100_iSS1_1M_rho_particle_counts.txt");
  if (!rho_particle_output_file.is_open()) {
    std::cerr << "Error opening particle pT counts output file" << std::endl;
    return 1;
  }
  
  //for (const auto& entry : rhoA_vs_num_particles) {
  //rhoA_particle_output_file << entry.first << " " << entry.second << std::endl; // x = number of particles in event, y = number of particles in jets
  //}
  for (const auto& entry : rho_particle_event_counts) {
    rho_particle_output_file << std::get<0>(entry) << ", " << std::get<1>(entry) << ", " << std::get<2>(entry) << "\n"; // x = rho*A, y = number of particles, z = event index
  }

  // Close the particles output file
  rho_particle_output_file.close();
  
  
  std::ofstream particles_output_file("pion_pp_pTHat15_100_hard_sector_hydro10_centr0_10_10k_trackpTcut_particles_in_event_and_jets.txt");
  //std::ofstream particles_output_file("pion_JETSCAPE_AuAu_out_only_hydro100_iSS1_1M_particles_in_event_and_jets.txt");
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
  
  
  return 0;
}
