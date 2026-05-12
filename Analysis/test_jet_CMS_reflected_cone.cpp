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
#include <limits>
#include <sstream>
#include <complex>
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
const double BIN_CENTER_SHIFT = 0.25;
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
long long total_perp_cones = 0;
long long total_perp_cone_pions = 0;
long long total_reflected_cones = 0;
long long total_reflected_cone_pions = 0;

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

struct QVector {
  std::complex<double> Qn;
  std::complex<double> Q2n;
  double sum_w;
  double sum_w2;
  int M;
  
  QVector() : Qn(0, 0), Q2n(0, 0), sum_w(0), sum_w2(0), M(0) {}
  
  void addParticle(double phi, int n, double weight) {
    Qn += std::polar(weight, n * phi);
    Q2n += std::polar(weight, 2 * n * phi);
    sum_w += weight;
    sum_w2 += weight * weight;
    M += (weight > 0) ? 1 : 0;
  }
};

double calcEta(double Px, double Py, double Pz) {
  double pz = Pz;
  double p = std::sqrt(Px * Px + Py * Py + pz * pz);
  double cosTheta = pz / p;
  if (cosTheta > 1.) cosTheta = 1.;
  if (cosTheta < -1.) cosTheta = -1.;
  double theta = std::acos(cosTheta);
  return -std::log(std::tan(theta / 2));
}

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

  std::map<double, long long> raw_jet_pion_pt_counts;
  std::map<double, long long> sub_jet_pion_pt_counts;
  std::map<double, long long> perp_cone_pion_pt_counts;
  std::map<double, long long> perp_cone_pion_pt_weighted_counts;
  std::map<double, long long> reflected_cone_pion_pt_counts;
  std::map<double, long long> reflected_cone_pion_pt_weighted_counts;
  
  std::map<double, long long> pion_pt_counts;
  std::map<double, long long> particle_pt_counts;
  
  std::map<double, long long> background_density_counts;
  std::map<int, int> particles_in_event_and_jets;
  //std::map<double, long long> rhoA_vs_num_particles;
  std::map<double, long long> rho_vs_num_particles;
  //std::vector<std::tuple<double, int, long long>> rhoA_particle_event_counts; // (rho * A, number of particles, event index)
  std::vector<std::tuple<double, int, long long>> rho_particle_event_counts; // (rho, number of particles, event index)
  std::vector<int> particle_counts;

  //std::ofstream subtracted_jet_etaphipt_output_file("pion_AuAu_pTHat15_100_nohydro_10M_events_subtracted_jet_trackpTcut_etaphipt_output.txt", std::ios::app);
  //if (!subtracted_jet_etaphipt_output_file.is_open()) {
  //  std::cerr << "Error opening subtracted jet output file" << std::endl;
  //  return 1;
  //}
  
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

	// Calculate event Q-vectors with eta gap for v2 calculation from charged pions above pt cut
	const double pt_cut_v2 = 0.2;          // GeV/c
	const double eta_gap_v2 = 0.3;         // Gap around midrapidity
	
	// Build Q-vectors for two subevents separated by eta gap
	QVector Q2_A, Q2_B;
	
	for (const auto& p : particles) {
	  if (std::abs(p.pid) == 211) { // charged pions only for flow calc
	    double pt = std::sqrt(p.Px * p.Px + p.Py * p.Py);
	    if (pt < pt_cut_v2) continue;
	    double eta = calcEta(p.Px, p.Py, p.Pz);
	    double phi = std::atan2(p.Py, p.Px);
	    
	    if (eta < -eta_gap_v2 / 2.) {
	      Q2_A.addParticle(phi, 2, 1.0);
	    }
	    else if (eta > eta_gap_v2 / 2.) {
	      Q2_B.addParticle(phi, 2, 1.0);
	    }
	  }
	}
	
	// Calculate v2 from Q-vectors
	double v2 = 0.0;
	double Psi2 = 0.0;
	
	if (Q2_A.sum_w > 1 && Q2_B.sum_w > 1) {
	  double c2 = std::real(Q2_A.Qn * std::conj(Q2_B.Qn)) / (Q2_A.sum_w * Q2_B.sum_w);
	  if (c2 > 0) v2 = std::sqrt(c2);
	  
	  std::complex<double> Q_comb = Q2_A.Qn + Q2_B.Qn;
	  Psi2 = 0.5 * std::atan2(Q_comb.imag(), Q_comb.real());
	  if (Psi2 < 0) Psi2 += M_PI;  // Normalize Psi2 to (0, pi)
	}
	else {
	  // Handle low multiplicity events gracefully, e.g., keep v2=0 and Psi2=0
	  v2 = 0;
	  Psi2 = 0;
	}

	std::cout << "Event " << event_number << ": v2 = " << v2 << ", Psi2 = " << Psi2 << std::endl;
	
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
	
	double rho0 = bge.rho(); // Get the background density
	
	
	fastjet::Subtractor subtractor(&bge);
	subtractor.set_use_rho_m(true);
	//std::vector<fastjet::PseudoJet> subtracted_jets = subtractor(jets);
	
	//int num_particles_in_jets = 0;
	//for (const auto& jet: jets){
	// num_particles_in_jets += jet.constituents().size();
	//}
	//particles_in_event_and_jets[num_particles_per_event]+= num_particles_in_jets;
	std::vector<fastjet::PseudoJet> subtracted_jets = subtractor(jets);

	  
	//double total_rho = 0.0;
	//int jet_number = 0;
	//random cone loop before jet loop?
	for (const auto& jet : subtracted_jets){
	  //double jet_area = jet.area(); // Accumulate the area for each jet
	  //double rhoA = rho * jet_area;
	  //total_rho += rho;
	  //jet_number++;
	  //double adjusted_rho = rho * jet_area;
	  double binned_rho = std::round((rho0) / RHO_BIN_SIZE) * RHO_BIN_SIZE;
	  background_density_counts[binned_rho]++; // Count occurrences of this rho
	  //rhoA_vs_num_particles[binned_rho] += num_particles_in_event;
	  rho_particle_event_counts.emplace_back(rho0, num_particles_per_event, processed_events);
	  //total_rhoA += binned_rho;
	  //jet_number;
	  //rhoA_vs_num_particles.emplace_back(binned_rho, num_particles_in_event);
	  //fastjet::PseudoJet subtracted_jet = subtractor(jet);
	  
	  //if (subtracted_jet.pt() < MIN_JET_PT || std::fabs(subtracted_jet.eta()) > JET_ABS_ETA_MAX || subtracted_jet.E() < 0.0) continue;
	  //std::vector<fastjet::PseudoJet> raw_constituents = jet.constituents();
	  std::vector<fastjet::PseudoJet> sub_constituents = jet.constituents();
	  if (sub_constituents.size() > 1){
	    //total_jets++;
	    
	    
	    double subtracted_jet_pt = jet.pt();
	    if (subtracted_jet_pt > JET_PT_CUT){

	      
	      double sub_jet_pt = subtracted_jet_pt;
	      double binned_sub_jet_pt = std::round(sub_jet_pt / SUB_PT_BIN_SIZE) * SUB_PT_BIN_SIZE;
	      sub_jet_pt_counts[binned_sub_jet_pt]++;
	      
	      	      
	      if (subtracted_jet_pt >= MIN_JET_PT && subtracted_jet_pt < MAX_JET_PT){
		//if (original_jet_pt >= MIN_JET_PT && original_jet_pt <=MAX_JET_PT){
		//std::cout << "Subtracted jet pT: " << subtracted_jet_pt << "Original jet pT: " << original_jet_pt << std::endl;
		double jet_eta = jet.eta();
		double jet_phi = jet.phi();
		double jet_pt = jet.pt();
		double jet_event = processed_events;
		
		//subtracted_jet_etaphipt_output_file << jet_event << " " << jet_eta << " " << jet_phi << " " << jet_pt << std::endl;
		total_subtracted_jets++;

		double v2_jet = v2;
		double v2_cone = v2_jet;
		
		
		//Reflected Cones
		double reflected_cone_eta;
		double reflected_cone_phi = jet_phi;

		if (std::abs(jet_eta) > JET_RADIUS){
		  reflected_cone_eta = -jet_eta;
		} else if (jet_eta >= -JET_RADIUS && jet_eta < 0){
		  reflected_cone_eta = jet_eta + (2*JET_RADIUS);
		} else if (jet_eta >= 0 && jet_eta <= JET_RADIUS){
		  reflected_cone_eta = jet_eta - (2*JET_RADIUS);
		} else {
		  reflected_cone_eta = -jet_eta;
		}
		  
		  total_reflected_cones++;
		  // Build cone center pseudojet for searching cone particles
		  double E_reflected_cone = jet_pt * std::cosh(reflected_cone_eta);
		  double Px_reflected_cone = jet_pt * std::cos(reflected_cone_phi);
		  double Py_reflected_cone = jet_pt * std::sin(reflected_cone_phi);
		  double Pz_reflected_cone = jet_pt * std::sinh(reflected_cone_eta);
		  fastjet::PseudoJet reflected_cone_center(Px_reflected_cone, Py_reflected_cone, Pz_reflected_cone, E_reflected_cone);
		  
		  // Find particles inside cone radius
		  std::vector<fastjet::PseudoJet> reflected_cone_constituents;
		  for (const auto& p : input_particles) {
		    if (deltaR(reflected_cone_center, p) <= JET_RADIUS) {
		      reflected_cone_constituents.push_back(p);
		    }
		  }

		  if (reflected_cone_constituents.empty()) continue;
		  
		  for (const auto& c : reflected_cone_constituents) {
		    if (std::abs(c.user_index()) == 211) {  // charged pions only
		      const double reflected_cone_pion_pt = c.pt();
		      total_reflected_cone_pions++;
		      // CMS flow modulation correction:
		      // weight = 1 / [1 + 2 * v2_cone * v2_jet * cos(2(phi_cone - Psi2))]
		      double delta_phi = reflected_cone_phi - Psi2;
		      double denom = 1. + 2. * v2_cone * v2_jet * std::cos(2 * delta_phi);
		      if (std::abs(denom) < 1e-12) denom = 1e-12;
		      double flow_weight = 1. / denom;
		      
		      // Bin pion pt same way as for raw/subtracted jets
		      double shifted_reflected_cone_pion_pt = reflected_cone_pion_pt - BIN_CENTER_SHIFT;
		      int bin_index = static_cast<int>(std::floor(shifted_reflected_cone_pion_pt / PION_PT_BIN_SIZE + 0.5));
		      double binned_reflected_cone_pion_pt = BIN_CENTER_SHIFT + bin_index * PION_PT_BIN_SIZE;
		      
		      // Accumulate weighted counts rounded to nearest integer
		      reflected_cone_pion_pt_counts[binned_reflected_cone_pion_pt]++;
		      reflected_cone_pion_pt_weighted_counts[binned_reflected_cone_pion_pt] += flow_weight;
		    }
		  }
		  //}
		
		//Perpendicular Cones
		std::vector<double> perp_cone_phis = {jet_phi + M_PI / 2., jet_phi - M_PI / 2.};
		for (auto& ph : perp_cone_phis) {
		  while (ph > M_PI) ph -= 2 * M_PI;
		  while (ph <= -M_PI) ph += 2 * M_PI;
		}
		
		// Use v2_jets calculated event-wise, and assume cone v2 same as jet for modulation
		
		for (double perp_cone_phi : perp_cone_phis) {
		  total_perp_cones++;

		  // Build cone center pseudojet for searching cones particles
		  double E_perp_cone = jet_pt * std::cosh(jet_eta);
		  double Px_perp_cone = jet_pt * std::cos(perp_cone_phi);
		  double Py_perp_cone = jet_pt * std::sin(perp_cone_phi);
		  double Pz_perp_cone = jet_pt * std::sinh(jet_eta);
		  fastjet::PseudoJet perp_cone_center(Px_perp_cone, Py_perp_cone, Pz_perp_cone, E_perp_cone);
		  
		  // Find particles inside cone radius
		  std::vector<fastjet::PseudoJet> perp_cone_constituents;
		  for (const auto& p : input_particles) {
		    if (deltaR(perp_cone_center, p) <= JET_RADIUS) {
		      perp_cone_constituents.push_back(p);
		    }
		  }

		  if (perp_cone_constituents.empty()) continue;
		  
		  // Loop over cone pions: apply CMS flow correction and accumulate pT bins as you do for jets
		  for (const auto& c : perp_cone_constituents) {
		    if (std::abs(c.user_index()) == 211) {  // charged pions only
		      const double perp_cone_pion_pt = c.pt();
		      total_perp_cone_pions++;
		      // CMS flow modulation correction:
		      // weight = 1 / [1 + 2 * v2_cone * v2_jet * cos(2(phi_cone - Psi2))]
		      double delta_phi = perp_cone_phi - Psi2;
		      double denom = 1. + 2. * v2_cone * v2_jet * std::cos(2 * delta_phi);
		      if (std::abs(denom) < 1e-12) denom = 1e-12;
		      double flow_weight = 1. / denom;
		      
		      // Bin pion pt same way as for raw/subtracted jets
		      double shifted_perp_cone_pion_pt = perp_cone_pion_pt - BIN_CENTER_SHIFT;
		      int bin_index = static_cast<int>(std::floor(shifted_perp_cone_pion_pt / PION_PT_BIN_SIZE + 0.5));
		      double binned_perp_cone_pion_pt = BIN_CENTER_SHIFT + bin_index * PION_PT_BIN_SIZE;
		      
		      // Accumulate weighted counts rounded to nearest integer
		     perp_cone_pion_pt_counts[binned_perp_cone_pion_pt]++;
		     perp_cone_pion_pt_weighted_counts[binned_perp_cone_pion_pt] += flow_weight;
		    }
		  }
		}
		
		int sub_pt_bin = std::min(static_cast<int>(subtracted_jet_pt / PT_BIN_WIDTH), NUM_PT_BINS - 1);
		//int pt_bin = std::min(static_cast<int>(original_jet_pt / PT_BIN_WIDTH), NUM_PT_BINS - 1);
		sub_dndpt_bins[sub_pt_bin].jet_count++;
		
		
		double total_subtracted_jet_momentum = jet.modp();
		
		double sub_jet_cut_pt = subtracted_jet_pt;
		double binned_sub_jet_cut_pt = std::round(sub_jet_cut_pt / SUB_CUT_PT_BIN_SIZE) * SUB_CUT_PT_BIN_SIZE;
		sub_jet_cut_pt_counts[binned_sub_jet_cut_pt]++;
		
		
		
		//int pt_bin = std::min(static_cast<int>(original_jet_pt / PT_BIN_WIDTH), NUM_PT_BINS - 1);
		for (const auto& constituent : sub_constituents) {
		  //for (const auto& constituent : raw_constituents) {
		  if (std::fabs(constituent.user_index()) == 211) { // PID 22 is for photons, abs(211) is for charged pions
		    double delta_r = deltaR(jet, constituent);
		    if (delta_r <= JET_RADIUS) {
		      double original_pion_pt = constituent.pt();
		      total_subtracted_pions++;

		      double shifted_sub_pion_pt = original_pion_pt - BIN_CENTER_SHIFT;
		      int bin_index = static_cast<int>(std::floor(shifted_sub_pion_pt / PION_PT_BIN_SIZE + 0.5));
		      double binned_sub_pion_pt = BIN_CENTER_SHIFT + bin_index * PION_PT_BIN_SIZE;

		      //double binned_sub_pion_pt = std::round(original_pion_pt / PION_PT_BIN_SIZE) * PION_PT_BIN_SIZE;
		      sub_jet_pion_pt_counts[binned_sub_pion_pt]++;

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

  //subtracted_jet_etaphipt_output_file.close();
	  
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
  std::cout << "Total perpendicular cones: " << total_perp_cones << std::endl;
  std::cout << "Total reflected cones: " << total_reflected_cones << std::endl;
  std::cout << "Total pions: " << total_pions << std::endl;
  std::cout << "Total raw pions: " << total_raw_pions << std::endl;
  std::cout << "Total subtracted pions: " << total_subtracted_pions << std::endl;
  std::cout << "Total perp cone pions: " << total_perp_cone_pions << std::endl;
  std::cout << "Total reflected cone pions: " << total_reflected_cone_pions << std::endl;
  
  std::ofstream perp_cone_pion_pt_output_file("pion_JETSCAPE_AuAu_pTHat15_100_hard_soft_hydro10_SMASH_only_decay_100k_perpendicular_cone_pion_pt_counts.txt");
  if (!perp_cone_pion_pt_output_file.is_open()) {
    std::cerr << "Error opening output file for corrected perpendicular cone pions" << std::endl;
    return 1;
  }
  
  for (const auto& kv : perp_cone_pion_pt_counts) {
    perp_cone_pion_pt_output_file << kv.first << " " << kv.second << "\n";
  }
  perp_cone_pion_pt_output_file.close();

  std::ofstream perp_cone_pion_pt_weighted_output_file("pion_JETSCAPE_AuAu_pTHat15_100_hard_soft_hydro10_SMASH_only_decay_100k_corrected_perpendicular_cone_pion_pt_counts.txt");
  if (!perp_cone_pion_pt_weighted_output_file.is_open()) {
    std::cerr << "Error opening output file for corrected perpendicular cone pions" << std::endl;
    return 1;
  }
  
  for (const auto& kv : perp_cone_pion_pt_weighted_counts) {
    perp_cone_pion_pt_weighted_output_file << kv.first << " " << kv.second << "\n";
  }
  perp_cone_pion_pt_weighted_output_file.close();


    std::ofstream reflected_cone_pion_pt_output_file("pion_JETSCAPE_AuAu_pTHat15_100_hard_soft_hydro10_SMASH_only_decay_100k_Jussi_method_reflected_cone_pion_pt_counts.txt");
  if (!reflected_cone_pion_pt_output_file.is_open()) {
    std::cerr << "Error opening output file for corrected perpendicular cone pions" << std::endl;
    return 1;
  }
  
  for (const auto& kv : reflected_cone_pion_pt_counts) {
    reflected_cone_pion_pt_output_file << kv.first << " " << kv.second << "\n";
  }
  reflected_cone_pion_pt_output_file.close();

  std::ofstream reflected_cone_pion_pt_weighted_output_file("pion_JETSCAPE_AuAu_pTHat15_100_hard_soft_hydro10_SMASH_only_decay_100k_Jussi_method_corrected_reflected_cone_pion_pt_counts.txt");
  if (!reflected_cone_pion_pt_weighted_output_file.is_open()) {
    std::cerr << "Error opening output file for corrected perpendicular cone pions" << std::endl;
    return 1;
  }
  
  for (const auto& kv : reflected_cone_pion_pt_weighted_counts) {
    reflected_cone_pion_pt_weighted_output_file << kv.first << " " << kv.second << "\n";
  }
  reflected_cone_pion_pt_weighted_output_file.close();

  /*
  std::ofstream raw_jet_pion_pt_output_file("pion_JETSCAPE_AuAu_pTHat15_100_nohydro_10M_trackpTcut_raw_jet_pion_pt_counts.txt");
  if (!raw_jet_pion_pt_output_file.is_open()) {
    std::cerr << "Error opening raw jet pion pT output file" << std::endl;
    return 1;
  }
  for (const auto& entry : raw_jet_pion_pt_counts) {
    raw_jet_pion_pt_output_file << entry.first << " " << entry.second << std::endl;
  }
  raw_jet_pion_pt_output_file.close();
  
  std::ofstream sub_jet_pion_pt_output_file("pion_JETSCAPE_AuAu_pTHat15_100_nohydro_10M_trackpTcut_sub_jet_pion_pt_counts.txt");
  if (!sub_jet_pion_pt_output_file.is_open()) {
    std::cerr << "Error opening subtracted jet pion pT output file" << std::endl;
    return 1;
  }
  for (const auto& entry : sub_jet_pion_pt_counts) {
    sub_jet_pion_pt_output_file << entry.first << " " << entry.second << std::endl;
  }
  sub_jet_pion_pt_output_file.close();
  */
  
  return 0;
}
