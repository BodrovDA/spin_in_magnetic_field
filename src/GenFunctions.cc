#include "GenFunctions.h"

// generate three vector with uniform distributions on hemisphere of radius p_mag    
Hep3Vector generate_3vector(double p_mag, std::default_random_engine &generator){
  std::normal_distribution<double> distribution(0.0,1.0);
  double x = distribution(generator);
  double y = distribution(generator);
  double z = distribution(generator);
  double len = sqrt(x*x + y*y + z*z);
  if(len == 0){
    len == 1;
  }
  return Hep3Vector(x * p_mag / len, y * p_mag / len, z * p_mag / len);
}

// create four vector with t component en and                                       
// with three vector with uniform distributions on hemisphere of radius p_mag      
HepLorentzVector generate_lorentz(double p_mag, double en,
                                  std::default_random_engine &generator){
  std::normal_distribution<double> distribution(0.0,1.0);
  double x = distribution(generator);
  double y = distribution(generator);
  double z = distribution(generator);
  double len = sqrt(x*x + y*y + z*z);
  if(len == 0){
    len == 1;
  }
  return HepLorentzVector(x * p_mag / len, y * p_mag / len, z * p_mag / len, en);
}

