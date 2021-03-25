#ifndef GEN_FUNCTIONS_H
#define GEN_FUNCTIONS_H

#include <cstdlib>
#include <ctime>
#include <random>

#include "LorentzVector.h"
#include "ThreeVector.h"


// generate uniform randorm value in [0, 1)
inline double doubleRand(){
  return double(rand()) / (double(RAND_MAX) + 1.0); 
}

// generate three vector with uniform distributions on hemisphere of radius p_mag
Hep3Vector generate_3vector(double p_mag, std::default_random_engine &generator);

// create four vector with t component en and
// with three vector with uniform distributions on hemisphere of radius p_mag
HepLorentzVector generate_lorentz(double p_mag, double en, 
				  std::default_random_engine &generator);

#endif
