#ifndef EVENT_H
#define EVENT_H

#include <cstdlib>
#include <ctime>
#include <cmath>
#include <random>

#include "LorentzVector.h"
#include "constants.h"
#include "MichelFunctions.icc"
#include "GenFunctions.h"
#include "Particle.h"
#include "SpinParticleMagField.h"

class Event {
public:
  Event(HepLorentzVector beam, int charge, double time_limit, 
	double time_limit_lab, std::default_random_engine &generator);
  ~Event();
  
  void generate_time();
  void generate_tau_decay(void);
  void generate_tau();
  void generate_mu();
  HepLorentzVector generate_spin_mu(HepLorentzVector &p4_tau_mu);
  void generate_e();
  HepLorentzVector generate_spin_e(HepLorentzVector &p4_mu_e);
  double width_tau();
  double width_mu();
  
  int MagFieldRotation_mu(int flag);
  void generate_cascade(int flag = 0);

  //  void dumpToTree();

  // get methods
  
  inline double time_mu();
  inline double time_lab_mu();
  inline double rotation_angle();
  inline particle tau();
  inline particle mu();
  inline particle mu_out();
  inline particle e();

private:

  HepLorentzVector beam_;
  int charge_;
  double time_limit_;
  double time_limit_lab_;
  std::default_random_engine &generator_;

  double time_mu_;
  double time_lab_mu_;
  double rotation_angle_;
  particle tau_;
  particle mu_;
  particle mu_rotated_;
  particle e_;

};

inline double Event::time_mu() {
  return time_mu_;
}
inline double Event::time_lab_mu() {
  return time_lab_mu_;
}
inline double Event::rotation_angle() {
  return rotation_angle_;
}
inline particle Event::tau() {
  return tau_;
}
inline particle Event::mu() {
  return mu_;
}
inline particle Event::mu_out() {
  return mu_rotated_;
}
inline particle Event::e() {
  return e_;
}


#endif
