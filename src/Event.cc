#include "Event.h"

Event::Event(HepLorentzVector beam, int charge, std::default_random_engine &generator) 
  : beam_(beam), charge_(charge), generator_(generator) {


}

Event::~Event(){};

void Event::generate_time(double time_limit) {
  time_mu_ = 0;
  double func = 0;
  double rand = 0;
  do{
    time_mu_ = doubleRand() * time_limit;
    func = exp(-time_mu_ / decay_time);
    rand = doubleRand();
  } while (rand > func);
}

void Event::generate_tau_decay(void) {
  double width_tau_ = 0;
  double rand = 0;
  do{
    generate_tau();
    generate_mu();
    width_tau_ = width_tau();
    rand = doubleRand();
  } while (rand > width_tau_);
}

void Event::generate_tau() {
  double energy_cms_tau = beam_.mag() / 2;
  double p_cms_tau = sqrt(energy_cms_tau * energy_cms_tau - m_tau * m_tau);
  HepLorentzVector p4_cms_tau = generate_lorentz(p_cms_tau, energy_cms_tau, generator_);
  HepLorentzVector spin_tau = generate_lorentz(1, 0, generator_);
  tau_ = particle(p4_cms_tau, spin_tau);
}

void Event::generate_mu() {
  double energy_tau_mu = mu_energy_min + doubleRand() * (w_mu - mu_energy_min);
  double p_tau_mu = sqrt(energy_tau_mu * energy_tau_mu - m_mu * m_mu);
  HepLorentzVector p4_tau_mu = generate_lorentz(p_tau_mu, energy_tau_mu, generator_);
  HepLorentzVector spin_mu = generate_spin_mu(p4_tau_mu);
  mu_ = particle(p4_tau_mu, spin_mu);

}

HepLorentzVector Event::generate_spin_mu(HepLorentzVector &p4_tau_mu){
  Hep3Vector axis_z = p4_tau_mu.vect( ) * (1. / p4_tau_mu.vect().mag());
  Hep3Vector axis_y = p4_tau_mu.vect().cross(tau_.s().vect()) *
    (1. / p4_tau_mu.vect().cross(tau_.s().vect()).mag());
  Hep3Vector axis_x = axis_y.cross(axis_z);
  
  double cost_tau = cos(p4_tau_mu.vect().angle(tau_.s().vect()));
  double sint_tau = sqrt(1 - cost_tau * cost_tau);
  double x_mu = p4_tau_mu.e() / w_mu;
  double x0_mu = m_mu / w_mu;

  double spin_denominator = Fis(x_mu, x0_mu) + charge_ * Fas(x_mu, x0_mu) * cost_tau;
  double spin_mu_z = (charge_ * Fip(x_mu, x0_mu)  + Fap(x_mu, x0_mu) * cost_tau) / spin_denominator;
  double spin_mu_x = Ft1(x_mu, x0_mu) * sint_tau / spin_denominator;
  double spin_mu_y = Ft2(x_mu, x0_mu) * sint_tau / spin_denominator;

  Hep3Vector PL_mu_vect = spin_mu_x * axis_x + spin_mu_y * axis_y + spin_mu_z * axis_z;
  Hep3Vector r_tmp = PL_mu_vect;

  Hep3Vector spin_mu_vect_tmp(0,0,0);
  if(PL_mu_vect.mag() != 0) {

    while((r_tmp.angle(PL_mu_vect) < const_pi / 12) || (r_tmp.angle(PL_mu_vect) > 11 * const_pi / 12))
      r_tmp = generate_3vector(1, generator_);

    Hep3Vector PL_mu_perp = r_tmp - (r_tmp.dot(PL_mu_vect) / PL_mu_vect.mag2()) * PL_mu_vect;
    Hep3Vector spin_mu_perp = (sqrt(1 - PL_mu_vect.mag2()) / PL_mu_perp.mag()) * PL_mu_perp;

    spin_mu_vect_tmp = PL_mu_vect + spin_mu_perp;
  } else {
    spin_mu_vect_tmp = generate_3vector(1, generator_);
  }
  
  return HepLorentzVector(spin_mu_vect_tmp, 0);

}

void Event::generate_e() {

}

double Event::width_tau() {
  double cost_tau = cos(mu_.p().vect().angle(tau_.s().vect()));
  double x_mu = mu_.p().e() / w_mu;
  double x0_mu = m_mu / w_mu;
  double width_tau = coeff_tau * sqrt(x_mu * x_mu - x0_mu * x0_mu)
    * (Fis(x_mu, x0_mu) - Fas(x_mu, x0_mu) * cost_tau);
  return width_tau;
}

double Event::width_mu() {
  return 0;
}
