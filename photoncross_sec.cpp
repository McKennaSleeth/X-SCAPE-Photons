#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdio>
#include <cmath>
#include <gsl/gsl_statistics.h>

struct Particle {
  int SN, PID, Status;
  double E, Px, Py, Pz, Eta, Phi;
};

std::vector<Particle> LoadParticlesFromFile(const std::string& filename, double& sigmaGen) {
  std::vector<Particle> particles;
  int eventCount = 0;

  std::ifstream ifs(filename);
  if (!ifs.is_open()) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return particles;
  }
  std::string str;
  while (std::getline(ifs, str)) {
    if (str.find("Event") != std::string::npos) {
      eventCount++;
    }
    if (str.find("sigmaGen") != std::string::npos) {
      double value;
      if (sscanf(str.c_str(), "# sigmaGen %lf", &value) == 1) {
        sigmaGen = value;
        std::cout << "Extracted sigmaGen: " << sigmaGen << std::endl;
      }
    } else {
      Particle p;
      sscanf(str.data(), "%d %d %d %lf %lf %lf %lf %lf %lf",
	     &p.SN, &p.PID, &p.Status, &p.E, &p.Px, &p.Py, &p.Pz, &p.Eta, &p.Phi);
      particles.push_back(p);
    }
    
  }
  std::cout << eventCount << " events in " << filename << std::endl;
  return particles;
}

int main() {
  std::string filename = "photon_pp_pTHat5_25_5Mil_out_final_state_hadrons.dat";
  /*
  std::vector<std::string> filenames = {
					"photon_pp_pTHat5_10_out_final_state_hadrons.dat",
					"photon_pp_pTHat10_15_out_final_state_hadrons.dat",
					"photon_pp_pTHat15_20_out_final_state_hadrons.dat",
					"photon_pp_pTHat20_25_out_final_state_hadrons.dat"
  };
  */
  //const double pTlist[] = {5.0, 10.0, 15.0, 20.0, 25.0};
  const double pTlist[] = {5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0};
  const int npT = sizeof(pTlist)/sizeof(double) - 1;
    
  //for (const auto& filename : filenames) {
  //std::string base_filename = filename.substr(0, filename.find("_out_final_state_hadrons.dat"));
  //std::string photon_outfile = base_filename + "combined_cross_sec.txt";
  std::string photon_outfile = "photon_pp_pTHat5_25_5Mil_cross_sec.txt";
    std::ofstream photon_out(photon_outfile);
    
    double sigmaGen = 0.0;
    std::vector<Particle> particles = LoadParticlesFromFile(filename, sigmaGen);
    std::cout << "Total particles in " << filename << ": " << particles.size() << std::endl;

    for (int a = 0; a < npT; a++) {
      double pTmin = pTlist[a];
      double pTmax = pTlist[a+1];
      std::vector<double> c(particles.size(), 0.0);
      for (const auto& particle : particles) {
	double pT = std::sqrt(particle.Px * particle.Px + particle.Py * particle.Py);
	double PID = particle.PID;
	double absEta = std::fabs(particle.Eta);
        if (PID == 22){
	  if (pT >= pTmin && pT < pTmax && absEta < 1.0) {
	    // Photon with PID 22 and within the pT range
	    c.push_back(sigmaGen);
	  } //end pT and eta selection loop
	} //end photon selection loop
      } //end particle loop
      
      double mean = gsl_stats_mean(c.data(), 1, c.size());
      double stan_err_mean = gsl_stats_sd_m(c.data(), 1, c.size(), mean) / std::sqrt(c.size());
      double deltaE = pTmax - pTmin;
      double deltaeta = 2;
      double avgE = (pTmax + pTmin) / 2;
      double cross_section = mean / (deltaE * deltaeta * 2 * M_PI * avgE);
      double cross_sec_err = (cross_section / mean) * stan_err_mean;
      photon_out << (pTmin+pTmax)/2 << " " << cross_section << " " << cross_sec_err << "\n";
      std::cout << (pTmin+pTmax)/2 << " " << cross_section << " " << cross_sec_err << "\n";
    } //end pT bin loop
    //} //end file loop
  photon_out.close();
  //}
  return 0;
} //end main function
