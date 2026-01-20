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
#include <fastjet/Selector.hh>

struct Particle {
  int id;
  int pid;
  int status;
  double E, Px, Py, Pz;
};

const int NUM_BINS = 20;
const double BIN_WIDTH = 1.0 / NUM_BINS;

const double PION_PT_BIN_SIZE = 0.2;
const double PARTICLE_PT_BIN_SIZE = 0.2;

const double RHO_BIN_SIZE = 0.2;

const int NUM_PT_BINS = 20;
const double PT_BIN_WIDTH = 100 / NUM_PT_BINS;

const double ETA_MAX_GHOST = 1.1;

long long total_raw_jets = 0;
long long total_subtracted_jets = 0;
long long total_raw_pions = 0;
long long total_subtracted_pions = 0;
long long total_pions = 0;


double deltaR(const fastjet::PseudoJet& p1, const fastjet::PseudoJet& p2) {
  double deta = p1.eta() - p2.eta();
  double dphi = std::abs(p1.phi() - p2.phi());
  if (dphi > M_PI) dphi = 2 * M_PI - dphi;
  return std::sqrt(deta * deta + dphi * dphi);
}

double calculate_manual_background_density(const std::vector<fastjet::PseudoJet>& input_particles) {
  // Define acceptance region
  double eta_max = ETA_MAX_GHOST;
  
  // Accumulate transverse momentum in acceptance
  double total_pt = 0.0;
  double total_area = 2 * eta_max * (2 * M_PI);  // Full azimuthal coverage
    
  for (const auto& particle : input_particles) {
    // Check if particle is within acceptance
    if (std::abs(particle.eta()) < eta_max) {
      total_pt += particle.pt();
    }
  }
    
  return total_pt / total_area;
}


int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <file1> <file2> ... <fileN>" << std::endl;
    return 1;
  }

  std::vector<std::string> filenames;
  for (int i = 1; i < argc; ++i) {
    filenames.push_back(argv[i]);
  }

  long long processed_events = 0;
  std::map<double, long long> pion_pt_counts;
  std::map<double, double> pion_pt_sum;
  std::map<double, double> pion_pt_sum_sq;
  
  std::map<double, long long> particle_pt_counts;
  std::map<double, double> particle_pt_sum;
  std::map<double, double> particle_pt_sum_sq;
  
  std::map<double, long long> background_total_pt_counts;
  std::map<double, long long> background_density_counts;
  std::map<double, long long> rho_vs_num_particles;
  std::vector<std::tuple<double, int, long long>> rho_particle_event_counts; // (rho, number of particles, event index)
  std::vector<int> particle_counts;

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
	    double binned_particle_pt = std::round(particle_pt / PARTICLE_PT_BIN_SIZE) * PARTICLE_PT_BIN_SIZE; // Bin the pT
	    particle_pt_counts[binned_particle_pt]++;
	    particle_pt_sum[binned_particle_pt] += particle_pt;
	    particle_pt_sum_sq[binned_particle_pt] += particle_pt * particle_pt;

	    if (std::fabs(p.pid) == 211){
	      double pion_pt = std::sqrt(p.Px * p.Px + p.Py * p.Py); // Calculate pT
	      double binned_pion_pt = std::round(pion_pt / PION_PT_BIN_SIZE) * PION_PT_BIN_SIZE; // Bin the pT
	      total_pions++;
	      pion_pt_counts[binned_pion_pt]++;
	      pion_pt_sum[binned_pion_pt] += pion_pt;
	      pion_pt_sum_sq[binned_pion_pt] += pion_pt * pion_pt;

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

	double manual_rho = calculate_manual_background_density(input_particles);
        double binned_rho = std::round(manual_rho / RHO_BIN_SIZE) * RHO_BIN_SIZE;
        background_density_counts[binned_rho]++;
        rho_particle_event_counts.emplace_back(manual_rho, num_particles_per_event, processed_events);

	double background_total_pt = 0.0;
	for (const auto& particle : input_particles) {
	  if (std::abs(particle.eta()) < ETA_MAX_GHOST) {
	    background_total_pt += particle.pt();
	  }
	}
	// Bin the total background pT with the same bin size as rho
	double binned_total_pt = std::round(background_total_pt / PT_BIN_WIDTH) * PT_BIN_WIDTH;
	background_total_pt_counts[binned_total_pt]++;
	double total_area = 2 * ETA_MAX_GHOST * (2 * M_PI);
	std::cout << "Event " << processed_events 
		  << ": total_area = " << total_area
		  << ", manual_rho = " << manual_rho
		  << ", total_background_pt = " << background_total_pt << std::endl;
	
	processed_events++;
	if (processed_events % 1000 == 0) {
	  std::cout << "Processed " << processed_events << " events in " << filename << std::endl;
	} //processed events condition
	
      } //event condition
    } //filename while loop
    
    file.close();
    std::cout << "Finished processing " << filename << std::endl;
  } //filename loop
	  
  std::cout << "Total processed events: " << processed_events << std::endl;
  std::cout << "Total pions: " << total_pions << std::endl;

  
  std::ofstream particle_pt_output_file("particle_JETSCAPE_AuAu_soft_sector_only_hydro_centr30_40_10_events_pt_counts.txt");
  //std::ofstream particle_pt_output_file("particle_JETSCAPE_AuAu_out_only_hydro100_iSS1_1M_pt_counts.txt");
  if (!particle_pt_output_file.is_open()) {
    std::cerr << "Error opening particle pT counts output file" << std::endl;
    return 1;
  }

    //particle_pt_output_file << "Particle pT (GeV/c) and counts: " << std::endl;
  for (const auto& entry : particle_pt_counts) {
    double bin = entry.first;
    long long N = entry.second;
    double mean = particle_pt_sum[bin] / N;
    double mean_sq = particle_pt_sum_sq[bin] / N;
    double variance = mean_sq - (mean * mean);
    double sem = (N > 0) ? std::sqrt(variance) / std::sqrt(N) : 0.0;
    particle_pt_output_file << bin << " " << N << " " << sem << std::endl; // x = binned pT, y = count
  }

    // Close the particle pT output file
  particle_pt_output_file.close();
  
  
  std::ofstream pion_pt_output_file("pion_JETSCAPE_AuAu_soft_sector_only_hydro_centr30_40_10_events_pt_counts.txt");
  //std::ofstream pion_pt_output_file("pion_JETSCAPE_AuAu_out_only_hydro100_iSS1_1M_pt_counts.txt");
  if (!pion_pt_output_file.is_open()) {
    std::cerr << "Error opening particle pT counts output file" << std::endl;
    return 1;
  }

    //particle_pt_output_file << "Particle pT (GeV/c) and counts: " << std::endl;
  for (const auto& entry : pion_pt_counts) {
    double bin = entry.first;
    long long N = entry.second;
    double mean = pion_pt_sum[bin] / N;
    double mean_sq = pion_pt_sum_sq[bin] / N;
    double variance = mean_sq - (mean * mean);
    double sem = (N > 0) ? std::sqrt(variance) / std::sqrt(N) : 0.0;
    pion_pt_output_file << bin << " " << N << " " << sem << std::endl; // x = binned pT, y = count
  }

    // Close the particle pT output file
  pion_pt_output_file.close();
  
  
  std::ofstream background_output_file("pion_JETSCAPE_AuAu_soft_sector_only_hydro_centr30_40_10_events_background_density_area_counts.txt");
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

  std::ofstream background_total_pt_file("pion_JETSCAPE_AuAu_soft_sector_only_hydro_centr30_40_10_events_background_pt_counts.txt");
  if (!background_total_pt_file.is_open()) {
    std::cerr << "Error opening background total pT output file" << std::endl;
    return 1;
  }
  
  for (const auto& entry : background_total_pt_counts) {
    background_total_pt_file << entry.first << " " << entry.second << std::endl;
  }
  
  background_total_pt_file.close();

  
  std::ofstream rho_particle_output_file("particle_AuAu_soft_sector_only_hydro_centr30_40_10_events_rho_particle_counts.txt");
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
  
  
  return 0;
}
