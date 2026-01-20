#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <complex>
#include <vector>
#include <cmath>
#include <numeric>

struct Particle {
  int id;
  int pid;
  int status;
  double E, Px, Py, Pz;
};

double calcPhi(double Px, double Py) {
  return std::atan2(Py, Px);
}

double calcPt(double Px, double Py) {
  return std::sqrt(Px*Px + Py*Py);
}

double calcEta(double Px, double Py, double Pz) {
  double pz = Pz;
  double p = std::sqrt(Px * Px + Py * Py + pz * pz);
  double cosTheta = pz / p;
  // Guard against floating point issues:
  if (cosTheta >  1.0) cosTheta =  1.0;
  if (cosTheta < -1.0) cosTheta = -1.0;
  double theta = std::acos(cosTheta);
  return -std::log(std::tan(theta / 2.0));
}

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


long long total_pions = 0;

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <file1> <file2> ... <fileN>" << std::endl;
    return 1;
  }

  constexpr double pt_cut = 0.2;       // GeV/c pT cut
  constexpr double eta_gap = 0.0;      // midrapidity gap width

  std::vector<double> c2_2_events;
  std::vector<double> c2_3_events;
    
  std::vector<std::string> filenames;
  for (int i = 1; i < argc; ++i) {
    filenames.push_back(argv[i]);
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
	    
	    if (std::fabs(p.pid) == 211){
	      total_pions++;
	    }
	  } else {
	    std::cerr << "Error parsing particle data" << std::endl;
	    return 1;
	  }
	}

	QVector Q2_A, Q2_B;
	QVector Q3_A, Q3_B;
	
	int pion_count = 0;
	double pion_pt_sum = 0.0;

	for(const auto& p : particles) {
	  if(std::abs(p.pid) == 211) {
	    double pt = calcPt(p.Px, p.Py);
	    if(pt < pt_cut) continue;
	    double eta = calcEta(p.Px, p.Py, p.Pz);
	    
	    if(eta < -eta_gap / 2.0) {
	      double phi = calcPhi(p.Px, p.Py);
	      Q2_A.addParticle(phi, 2, 1);
	      Q3_A.addParticle(phi, 3, 1);
	      pion_pt_sum += pt;
	      pion_count++;
	    }
	    else if(eta > eta_gap / 2.0) {
	      double phi = calcPhi(p.Px, p.Py);
	      Q2_B.addParticle(phi, 2, 1);
	      Q3_B.addParticle(phi, 3, 1);
	      pion_pt_sum += pt;
	      pion_count++;
	    }
	  }
	}
	
	double mean_pt = (pion_count > 0) ? (pion_pt_sum / pion_count) : 0.0;
			
	double c2_2 = 0.0;
	double c4_2 = 0.0;
	if(Q2_A.sum_w > 0 && Q2_B.sum_w > 0 && Q2_A.sum_w > 1 && Q2_B.sum_w >1) {
	  c2_2 = std::real(Q2_A.Qn * std::conj(Q2_B.Qn)) / (Q2_A.sum_w * Q2_B.sum_w);
	  
	  double numerator = (std::norm(Q2_A.Qn) - Q2_A.sum_w2) * (std::norm(Q2_B.Qn) - Q2_B.sum_w2) - (std::norm(Q2_A.Q2n) * std::norm(Q2_B.Q2n)) / (Q2_A.sum_w * Q2_B.sum_w);
	  
	  double denominator = Q2_A.sum_w * (Q2_A.sum_w - 1) * Q2_B.sum_w * (Q2_B.sum_w-1);
	  if(denominator > 0) {
	    c4_2 = numerator / denominator - 2.0 * c2_2 * c2_2;
	  }
	}
	
	// Weighted for v3
	double c2_3 = 0.0;
	double c4_3 = 0.0;
	if(Q3_A.sum_w > 0 && Q3_B.sum_w > 0 && Q3_A.sum_w > 1 && Q3_B.sum_w >1) {
	  c2_3 = std::real(Q3_A.Qn * std::conj(Q3_B.Qn)) / (Q3_A.sum_w * Q3_B.sum_w);
	  
	  double numerator = (std::norm(Q3_A.Qn) - Q3_A.sum_w2) * (std::norm(Q3_B.Qn) - Q3_B.sum_w2) - (std::norm(Q3_A.Q2n) * std::norm(Q3_B.Q2n)) / (Q3_A.sum_w * Q3_B.sum_w);

	  double denominator = Q3_A.sum_w * (Q3_A.sum_w - 1) * Q3_B.sum_w * (Q3_B.sum_w-1);
	  if(denominator > 0) {
	    c4_3 = numerator / denominator - 2.0 * c2_3 * c2_3;
	  }
	}
	
	// Store per-event c2 cumulants for mean and SEM calc
	if(c2_2 > 0) c2_2_events.push_back(c2_2);
	if(c2_3 > 0) c2_3_events.push_back(c2_3);
	
	// Calculate vn{2} and vn{4} per event for printing
	double v2_2 = (c2_2 > 0) ? std::sqrt(c2_2) : 0;
	double v2_4 = (c4_2 < 0) ? std::pow(-c4_2, 0.25) : 0;
	double v3_2 = (c2_3 > 0) ? std::sqrt(c2_3) : 0;
	double v3_4 = (c4_3 < 0) ? std::pow(-c4_3, 0.25) : 0;
	
	std::cout << "Event " << event_number
		  << " Pions: " << pion_count
		  << " Mean_pT: " << mean_pt
		  << " v2{2}: " << v2_2 << " v2{4}: " << v2_4
		  << " v3{2}: " << v3_2 << " v3{4}: " << v3_4
		  << std::endl;
	
	processed_events++;
	if (processed_events % 1000 == 0) {
	  std::cout << "Processed " << processed_events << " events in " << filename << std::endl;
	} //processed events condition
      } //event header condition
    } //filename while loop
    
    file.close();
    std::cout << "Finished processing " << filename << std::endl;
  } //filename loop

  auto computeMeanSEM = [](const std::vector<double>& vals, double& mean, double& sem) {
			  if(vals.empty()) { mean=0; sem=0; return;}
			  mean = std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
			  double sq_sum = 0.0;
			  for(auto v : vals) sq_sum += (v - mean)*(v - mean);
			  double variance = (vals.size() > 1) ? sq_sum / (vals.size() - 1) : 0;
			  sem = std::sqrt(variance/vals.size());
			};
  
  double mean_c2_2, sem_c2_2, mean_c2_3, sem_c2_3;
  computeMeanSEM(c2_2_events, mean_c2_2, sem_c2_2);
  computeMeanSEM(c2_3_events, mean_c2_3, sem_c2_3);
  
  double mean_v2_2 = (mean_c2_2 > 0) ? std::sqrt(mean_c2_2) : 0;
  double sem_v2_2 = (mean_c2_2 > 0) ? 0.5 * sem_c2_2 / std::sqrt(mean_c2_2) : 0;
  
  double mean_v3_2 = (mean_c2_3 > 0) ? std::sqrt(mean_c2_3) : 0;
  double sem_v3_2 = (mean_c2_3 > 0) ? 0.5 * sem_c2_3 / std::sqrt(mean_c2_3) : 0;
  
  std::cout << "\n=== Final averages with SEM ===" << std::endl;
  std::cout << "v2{2} = " << mean_v2_2 << " +/- " << sem_v2_2 << std::endl;
  std::cout << "v3{2} = " << mean_v3_2 << " +/- " << sem_v3_2 << std::endl;

  
  return 0;
  
}
