#include "ShellModel.hh"

using namespace arma;

//Slater determinant vector contains sets of integers (enforces Pauli Principle)

int constructSlaterDeterminants(std::vector< std::vector<int> > slater_determinants, int num_particles,  std::vector<State> sp_states){
  //For right now I'm just manually making this
  std::vector<int> temp_sd; 
  if (num_particles % 2 != 0){
    std::cout << "ERROR: num_particles isn't even and can't be handled by this pairing model!" << std::endl;
    return -1;
  }
  for (unsigned int i = 1; i <= sp_states.size(); i++){
    for (unsigned int j = i+1; j <= sp_states.size(); j++){
      temp_sd.clear();
      temp_sd.push_back(i);
      temp_sd.push_back(j);
      slater_determinants.push_back(temp_sd);
    }
  }
  return slater_determinants.size();
}

int readSingleParticleStates(std::string const &file_name, std::vector<State> &single_particle_states){
  std::ifstream input_file;
  input_file.open(file_name.c_str());
  if (!input_file.is_open()){
    std::cout << "Failed to open matrix elements file" << std::endl;
    return 0;
  }
  int state_index =-1;
  int n = -1;
  int l = -1; 
  int j2 = -1;
  int mj2 = -1;
  int tz2 = -1;
  std::string line;
  while(std::getline(input_file,line)){
    sscanf(line.c_str(), "Orbit number: %d %d %d %d %d %d", &state_index,&n,&l,&j2,&mj2,&tz2); 
    State state(state_index,n,l,j2,mj2,tz2);
    state.Print();
    single_particle_states.push_back(state);
  }
  return single_particle_states.size();
}

int buildHamiltonian(arma::mat &hamiltonian, std::vector< std::vector<int> > slater_determinants, double interaction_strength, std::vector<State> const sp_states){
  int num_pairs = slater_determinants.at(0).size();
  int num_sp_states = sp_states.size();
  for (unsigned int alpha = 0; alpha < slater_determinants.size(); alpha++){
    std::vector<int> alpha_vec = slater_determinants.at(alpha);
    for (unsigned int beta = alpha; beta < slater_determinants.size(); beta++){
      std::vector<int> beta_vec = slater_determinants.at(beta);
      double h0 = 0; 
      double interaction = 0;

      if (beta == alpha){
        for (int i = 0; i < num_pairs; i++){
          h0 += beta_vec.at(i);
        }
        h0 -= num_pairs;
        h0 *= 2;
      }

      for (int p = 1; p <= num_sp_states; p++){
        for (int q = 1; q <= num_sp_states; q++){
          std::vector<int>::iterator q_pos = std::find(beta_vec.begin(), beta_vec.end(), q);
          if (q_pos == beta_vec.end()){
            continue;
          }

          if (p != q){
            beta_vec.erase(q_pos);
            if (std::find(beta_vec.begin(), beta_vec.end(), p) != beta_vec.end()){
              continue;
            }
            beta_vec.push_back(p); 
            if (beta_vec.size() != alpha_vec.size()){
              std::cout << "ERROR: Incompatible vector sizes! beta size ("<<beta_vec.size()<<") and alpha size("<<alpha_vec.size()<<")"<<std::endl;
              return -1;
            }
          }//p!=q

          //compare alpha & beta
          std::sort(alpha_vec.begin(), alpha_vec.end());
          std::sort(beta_vec.begin(), beta_vec.end());
          if (alpha_vec == beta_vec){
            interaction += interaction_strength; 
          }
        }//loop over q
      }//loop over p
      hamiltonian(alpha,beta) = hamiltonian(beta,alpha) = h0-interaction;
    }//loop over beta
  }//loop over alpha
  return 0;
}

int runShellModelCalculation(std::string single_particle_states_file_name, int num_particles, double interaction_strength){

  std::vector<State> single_particle_states;
  std::vector< std::vector<int> > slater_determinants;

  if (readSingleParticleStates(single_particle_states_file_name, single_particle_states) == 0){
    std::cout << "ERROR: Failed to read in any states! Exiting." << std::endl;
    return -1;
  };

  //currently doesn't do anything with the information given, just fills Slater 
  //Determinants with known values
  //Need way to know number of slater determinants; this means I think that a vector
  //of SlaterDeterminant objects is the best choice; will implement later
  int num_slater_determinants = constructSlaterDeterminants(slater_determinants, num_particles, single_particle_states);
  if (num_slater_determinants == 0){
    std::cout << "ERROR: No slater determinants found!" << std::endl;
    return -1;
  }
  mat hamiltonian = zeros(num_slater_determinants, num_slater_determinants);
  vec energies = zeros(num_slater_determinants,1);
  vec eigenvectors = zeros(num_slater_determinants,1);

  if (buildHamiltonian(hamiltonian, slater_determinants, interaction_strength, single_particle_states) != 0){
    std::cout << "ERROR: Failed to build Hamiltonian!" << std::endl;
    return -1;
  }

  eig_sym(energies, eigenvectors, hamiltonian);

  return 0; 
}


int main(int argc, char **argv){
  std::string usage("ShellModel [single particle states file name]  [number of particles] ");

  if (argc < 4){
    std::cout << usage << std::endl;
    return -1;
  }

//if (num_particles%2 != total_m%2){
//  std::cout << "Error! Incompatible number of particles ("<<num_particles<<") and total m value ("<<total_m<<")"<<std::endl;
//  return -1;
//}

  std::string sp_state_filename(argv[1]);
  int num_particles = atoi(argv[2]);
  double interaction_strength = 0.5;
  runShellModelCalculation(sp_state_filename, num_particles,interaction_strength);

}
