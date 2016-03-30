#ifndef __SHELL_MODEL_G__
#define __SHELL_MODEL_G__
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <armadillo>
#include "State.hh"


int readSingleParticleStates(std::string const &file_name, std::vector<State> &single_particle_states);
int runShellModelCalculation(std::string, int, int);
int constructSlaterDeterminants(std::vector<std::vector<int>> slater_determinants, int num_particles, std::vector<State> states);

//See notes on Full Configuration interaction, pg. 25) for meaning of variables
//p,q,r,s are single particle states; alpha is slate_determinant index
int buildHamiltonian(arma::mat &hamiltonian, std::vector<std::vector<int>> slater_determinants, double interaction_strength, std::vector<State> const sp_states);


#endif