#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <random>
#include <sstream>

// Struct for particle info
struct Particle {
    int id;
    int pid;
    int status;
    double E, Px, Py, Pz;

    double pT() const { return std::sqrt(Px*Px + Py*Py); }

    double eta() const {
        double p = std::sqrt(Px*Px + Py*Py + Pz*Pz);
        if (p != Pz) {
            double arg = (p + Pz) / (p - Pz);
            if (arg <= 0) return 0.0;
            return 0.5 * std::log(arg);
        }
        return 0.0;
    }

    double phi() const { return std::atan2(Py, Px); }
};

// Struct for jets read from file
struct JetInfo {
    int event_number; // from jet file - not used for matching
    double eta;
    double phi;
    double pT;
};

// Returns deltaR distance between eta/phi coordinates
double deltaR(double eta1, double phi1, double eta2, double phi2) {
    double dEta = eta1 - eta2;
    double dPhi = std::fabs(phi1 - phi2);
    if(dPhi > M_PI) dPhi = 2*M_PI - dPhi;
    return std::sqrt(dEta*dEta + dPhi*dPhi);
}

bool isPion(const Particle& p) {
    return std::abs(p.pid) == 211;
}

// === Parameters for binning and cone ===
const double cone_radius = 0.4;
const double pion_pt_bin_width = 0.5; 
const int num_pion_pt_bins = 40;  // Cover 0 to 20 GeV pT bins

// Reads jets from a file: lines of "event_number eta phi pT"
bool readJets(const std::string& jet_filename, std::vector<JetInfo> &jets) {
    std::ifstream infile(jet_filename);
    if(!infile.is_open()) {
        std::cerr << "Failed to open jet file: " << jet_filename << std::endl;
        return false;
    }
    std::string line;
    while(std::getline(infile, line)) {
        if(line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        JetInfo jet;
        if(ss >> jet.event_number >> jet.eta >> jet.phi >> jet.pT) {
            jets.push_back(jet);
        }
    }
    return true;
}

// Reads hydro/particle event format files, returns vector of particle-vectors, one per event
bool readParticleEvents(const std::vector<std::string> &filenames, std::vector<std::vector<Particle>> &eventsParticles) {
    for(const auto &fname: filenames) {
        std::ifstream infile(fname);
        if(!infile.is_open()) {
            std::cerr << "Failed to open particle file " << fname << std::endl;
            continue;
        }

        std::string line;
        int event_number = -1, n_hadrons = 0;
        while(std::getline(infile, line)) {
            // Look for event header line e.g.
            // "#   Event 123 weight ... EPangle ... N_hadrons 20"
            if(sscanf(line.c_str(), "#\tEvent\t%d\tweight\t%*f\tEPangle\t%*f\tN_hadrons\t%d", &event_number, &n_hadrons) == 2) {
                std::vector<Particle> eventParticles;
                eventParticles.reserve(n_hadrons);

                for(int i=0; i<n_hadrons; ++i) {
                    if(!std::getline(infile, line)) {
                        std::cerr << "Unexpected EOF reading particles in event " << event_number << std::endl;
                        return false;
                    }
                    Particle p;
                    if(sscanf(line.c_str(), "%d %d %d %lf %lf %lf %lf", &p.id, &p.pid, &p.status, &p.E, &p.Px, &p.Py, &p.Pz) != 7) {
                        std::cerr << "Error reading particle line in event " << event_number << std::endl;
                        return false;
                    }
                    // Apply a minimum pT cut for tracks (adjust if needed)
                    if(p.pT() > 0.5) {
		      eventParticles.push_back(p);
		    }
                }
                eventsParticles.push_back(std::move(eventParticles));
            }
        }
    }
    return true;
}

int main(int argc, char* argv[]) {
    if(argc < 3) {
        std::cerr << "Usage: "<< argv[0] << " <jet_file> <particle_file1> [particle_file2 ...]" << std::endl;
        return 1;
    }

    // Read jets
    std::string jetFileName = argv[1];
    std::vector<JetInfo> jets;
    if(!readJets(jetFileName, jets)) {
        std::cerr << "Failed reading jet file." << std::endl;
        return 1;
    }
    if(jets.empty()) {
        std::cerr << "No jets found in file." << std::endl;
        return 1;
    }
    std::cout << "Read " << jets.size() << " jets." << std::endl;

    // Read particle events from hydro files
    std::vector<std::string> particleFiles;
    for(int i=2; i<argc; ++i) particleFiles.push_back(argv[i]);
    
    std::vector<std::vector<Particle>> eventsParticles;
    if(!readParticleEvents(particleFiles, eventsParticles)) {
        std::cerr << "Error reading particle event files." << std::endl;
        return 1;
    }
    if(eventsParticles.empty()) {
        std::cerr << "No background particle events loaded." << std::endl;
        return 1;
    }
    std::cout << "Read " << eventsParticles.size() << " particle events." << std::endl;

    // Set up pion pT bins accumulator
    std::vector<long long> pion_pt_counts(num_pion_pt_bins, 0);

    size_t total_events = eventsParticles.size();
    size_t bg_event_index = 0; // Cycle index over background events
    static int jet_counter = 0; // global or static to keep count

    // Loop over jets, assign background events cyclically (wrap-around)
    for(const auto &jet: jets) {
        const auto &bg_event_particles = eventsParticles[bg_event_index];
        bg_event_index = (bg_event_index + 1) % total_events;
	
	std::cout << "Jet " << jet_counter << " (eta=" << jet.eta << ", phi=" << jet.phi << ") assigned to background event " << bg_event_index << std::endl;
	jet_counter++;
        // Count pions in cone around jet eta, phi and bin by pT
        for(const auto &p: bg_event_particles) {
            if(!isPion(p)) continue;
            if(deltaR(p.eta(), p.phi(), jet.eta, jet.phi) <= cone_radius) {
                double pt = p.pT();
                int bin = static_cast<int>(pt / pion_pt_bin_width);
                if(bin >= 0 && bin < num_pion_pt_bins) {
                    pion_pt_counts[bin]++;
                }
            }
        }
    }

    // Write pion pT bin counts to output file
    std::ofstream outFile("pion_JETSCAPE_AuAu_pTHat15_100_hard_soft_hydro10_SMASH_only_decay_10k_trackpTcut_cone_pt_counts_v3.txt");
    if(!outFile.is_open()) {
        std::cerr << "Failed to open output file for writing." << std::endl;
        return 1;
    }
    outFile << "#pT_bin_center(GeV) pion_count\n";
    for(int i=0; i<num_pion_pt_bins; ++i) {
        double bin_center = (i + 0.5) * pion_pt_bin_width;
        outFile << bin_center << " " << pion_pt_counts[i] << "\n";
    }
    outFile.close();

    std::cout << "Pion counts per pT bin written to pion_counts_random_event_cones_pt.txt\n";

    return 0;
}
