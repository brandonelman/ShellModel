#include "ShellModel.hh"

using namespace arma;

//Slater determinant vector contains sets of integers (enforces Pauli Principle)

int constructSlaterDeterminants(std::vector< std::vector<int> > &slater_determinants, int num_particles,  int num_sp_states){
  //For right now I'm just manually making this
  std::vector<int> temp_sd; 
  if (num_particles % 2 != 0){
    std::cout << "ERROR: num_particles isn't even and can't be handled by this pairing model!" << std::endl;
    return -1;
  }
  int num_pairs = num_particles/2;
  int num_pair_states = num_sp_states/2;
  if (num_pairs == 2){
    for (int i = 1; i < num_pair_states; i++){
      for (int j = i+1; j <= num_pair_states; j++){
        temp_sd.clear();
        temp_sd.push_back(i);
        temp_sd.push_back(j);
        slater_determinants.push_back(temp_sd);
      }
    }
  }

  else if (num_pairs == 3){
    for (int i = 1; i < num_pair_states; i++){
      for (int j = i+1; j < num_pair_states; j++){
        for (int k = j+1; k <=num_pair_states;k++){
          temp_sd.clear();
          temp_sd.push_back(i);
          temp_sd.push_back(j);
          temp_sd.push_back(k);
          slater_determinants.push_back(temp_sd);
        }
      }
    }
  }
  else if (num_pairs ==4){
    for (int i = 1; i < num_pair_states; i++){
      for (int j = i+1; j < num_pair_states; j++){
        for (int k = j+1; k <=num_pair_states;k++){
          for (int l = k+1; l <=num_pair_states;l++){
            temp_sd.clear();
            temp_sd.push_back(i);
            temp_sd.push_back(j);
            temp_sd.push_back(k);
            temp_sd.push_back(l);
            slater_determinants.push_back(temp_sd);
          }
        }
      }
    }
  }//num_Pairs == 4
  else{ 
    std::cout << "ERROR: Bad number of pairs: " << num_pairs << std::endl;
    return 0;
  }
  return slater_determinants.size();
}

int buildHamiltonian(arma::mat &hamiltonian, std::vector< std::vector<int> > slater_determinants, double interaction_strength, int num_sp_states){
  int num_pairs = slater_determinants.at(0).size();
  int num_pair_states = num_sp_states/2;
//std::cout << "Number of pairs: " << num_pairs << std::endl;
//std::cout << "Number of pair states: " << num_pair_states << std::endl;
//std::cout << "Slater determinant vector size: " << slater_determinants.size() <<  std::endl;
//std::cout << "Slater determinant size: " << slater_determinants.at(0).size() <<  std::endl;
  for (unsigned int alpha = 0; alpha < slater_determinants.size(); alpha++){
    std::vector<int> alpha_vec = slater_determinants.at(alpha);
    for (unsigned int beta = alpha; beta < slater_determinants.size(); beta++){
//    std::cout << "\n\nWorking on element ("<<alpha<<","<<beta<<")"<<std::endl;
      std::vector<int> beta_vec = slater_determinants.at(beta);
      double h0 = 0; 
      double interaction = 0;

//    std::cout << "\tComparing vector alpha with elements |"<<alpha_vec.at(0)<<" "<<alpha_vec.at(1)<<"> with vector beta"
//              << " with elements |"<<beta_vec.at(0)<<" "<<beta_vec.at(1)<<">"<<std::endl;
      
      if (beta == alpha){
        for (int i = 0; i < num_pairs; i++){
          h0 += beta_vec.at(i);
        }
        h0 -= num_pairs;
        h0 *= 2;
      }

      for (int p = 1; p <= num_pair_states; p++){
        for (int q = 1; q <= num_pair_states; q++){
          beta_vec = slater_determinants.at(beta);
//        std::cout << "\tIteration with p = " << p << " and q = " << q << std::endl;
          std::vector<int>::iterator q_pos = std::find(beta_vec.begin(), beta_vec.end(), q);
          std::vector<int>::iterator p_pos = std::find(beta_vec.begin(), beta_vec.end(), p);

          if (q_pos == beta_vec.end()){
            continue;
          }
          if (p != q && p_pos != beta_vec.end()){
            continue;
          }

          if (p != q){
            beta_vec.erase(q_pos);
            beta_vec.push_back(p); 
            if (beta_vec.size() != alpha_vec.size()){
              std::cout << "ERROR: Incompatible vector sizes! beta size ("<<beta_vec.size()<<") and alpha size("<<alpha_vec.size()<<")"<<std::endl;
              return -1;
            }
          }//p!=q

          //compare alpha & beta
          std::sort(alpha_vec.begin(), alpha_vec.end());
          std::sort(beta_vec.begin(), beta_vec.end());
//        std::cout << "\tVectors After Sorting" << std::endl;
//        std::cout << "\tComparing vector alpha with elements |"<<alpha_vec.at(0)<<" "<<alpha_vec.at(1)<<"> with vector beta"
//                  << " with elements |"<<beta_vec.at(0)<<" "<<beta_vec.at(1)<<">"<<std::endl;
          if (alpha_vec == beta_vec){
//          std::cout << "\tNon zero term!\n\n";
            interaction += interaction_strength; 
          }
        }//loop over q
      }//loop over p

      h0 = 0;
      hamiltonian(alpha,beta) = hamiltonian(beta,alpha) = h0-interaction;
    }//loop over beta
  }//loop over alpha
  return 0;
}

int runShellModelCalculation(int num_sp_states, int num_particles, double interaction_strength){

  std::vector< std::vector<int> > slater_determinants;

  std::cout << "Constructing slater determinants" << std::endl;
  int num_slater_determinants = constructSlaterDeterminants(slater_determinants, num_particles, num_sp_states);
  if (num_slater_determinants == 0){
    std::cout << "ERROR: No slater determinants found!" << std::endl;
    return -1;
  }
  std::cout << "Initializing Hamiltonian...." << std::endl;
  mat hamiltonian = zeros(num_slater_determinants, num_slater_determinants);
  vec energies = zeros(num_slater_determinants,1);
  mat eigenvectors = zeros(num_slater_determinants,num_slater_determinants);

  std::cout << "Building Hamiltonian...." << std::endl;
  if (buildHamiltonian(hamiltonian, slater_determinants, interaction_strength, num_sp_states) != 0){
    std::cout << "ERROR: Failed to build Hamiltonian!" << std::endl;
    return -1;
  }

  std::cout << "Diagonalizing..." << std::endl;
  eig_sym(energies, eigenvectors, hamiltonian);

  std::cout << "Saving files...." << std::endl;
  std::ofstream energy_file("energies.dat");
  energy_file << energies;
  energy_file.close();
  std::ofstream mat_file("out_mat_file.dat");
  mat_file << hamiltonian;
  mat_file.close();
  return 0; 
}


int main(int argc, char **argv){
  std::string usage("ShellModel [num single particle states]  [number of particles] ");

  if (argc < 3){
    std::cout << usage << std::endl;
    return -1;
  }

//if (num_particles%2 != total_m%2){
//  std::cout << "Error! Incompatible number of particles ("<<num_particles<<") and total m value ("<<total_m<<")"<<std::endl;
//  return -1;
//}

//std::string sp_state_filename(argv[1]);
  
  int num_sp_states = atoi(argv[1]);
  int num_particles = atoi(argv[2]);
  double interaction_strength = 0.5;
  runShellModelCalculation(num_sp_states, num_particles,interaction_strength);

}
