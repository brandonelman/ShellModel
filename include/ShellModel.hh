#ifndef __SHELL_MODEL_G__
#define __SHELL_MODEL_G__
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <armadillo>
#include <fstream>
#include "State.hh"


int readSingleParticleStates(std::string const &file_name, int num_sp_states);
int runShellModelCalculation(int, int, int);
int constructSlaterDeterminants(std::vector< std::vector<int> > &slater_determinants, int num_particles, int num_sp_states);

//See notes on Full Configuration interaction, pg. 25) for meaning of variables
//p,q,r,s are single particle states; alpha is slate_determinant index
int buildHamiltonian(arma::mat &hamiltonian, std::vector< std::vector<int> > slater_determinants, double interaction_strength, int num_sp_states);


#endif
