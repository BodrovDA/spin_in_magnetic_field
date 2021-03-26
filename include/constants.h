#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "ThreeVector.h"

// math constants
const double const_pi = 3.1415926;

// particles decay times
const double decay_time = 2.2E-6;

// particle masses
const double m_tau = 1.7769;
const double m_mu = 0.1057; 
const double m_e = 0.00051;

// kinematic parameters                                            
const double w_mu = (m_tau * m_tau + m_mu * m_mu) / 2 / m_tau;
const double w_e = (m_mu * m_mu + m_e * m_e) / 2 / m_mu;
const double mu_energy_min = m_mu;
const double e_energy_min = m_e;
const double coeff_tau =  w_mu*w_mu*w_mu*w_mu / (m_tau*m_tau*m_tau*m_tau) * 24 * 2;

// Michel parameters
// already measured
const double rho = 0.75;
const double delta = 0.75;
const double xi = 1.;
const double eta = 0;

// not measured yet
const double xi_p = 1;
const double xi_pp = xi * xi_p;   
const double eta_pp = 0;
const double alpha_p = 0;
const double beta_p = 0;
const double A = 16;

// old notation
const double rho_p = 3. / 4;
const double delta_p = 1. / 4;
const double eta_p = 0;
const double alpha = 0;
const double beta = 0;
const double a = 0;
const double b = 4;
const double c = 0;

// physics constants
const double magnb = 14E9 * m_e / m_mu;

// magnetic field
const Hep3Vector BfieldC(0., 0., 1.);
const Hep3Vector BfieldB(0., 0., 1.5);

// time step in particle propagation calculation
const double dt = 0.5 / 3. / 1.E10;


#endif 
