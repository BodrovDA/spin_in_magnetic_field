#include "SpinParticleMagField.h"


Hep3Vector spin_rotation(Hep3Vector &zeta, HepLorentzVector p, double time_dec){

  double time = 10000 * dt;
  // three vectors
  Hep3Vector p_vect = p.vect();
  Hep3Vector sperp1(zeta.x(), zeta.y(), 0);
  
  // Helix parameters
  double kappa = 1 / p_vect.perp();
  double tanlam = p_vect.z() * kappa;
  double alpha = 10000. / 3. / 10. / BfieldB.mag();
  double radius = alpha / kappa;

  // for tests
  double distance1 = 0;
               
  for(int i = 0; i < time_dec / dt; ++i){
    distance1 += p_vect.mag() / p.e() * 3E10 * dt;
    Hep3Vector dzeta = 2 * magnb * m_mu / p.e() * zeta.cross(BfieldB) * dt * 6.30253;
    zeta += dzeta;
    Hep3Vector dp = 2 * magnb * m_mu / p.e() * p_vect.cross(BfieldB) * dt * 6.30253;
    p_vect += dp;
       
  }

  // for tests
  Hep3Vector perp1(p.vect().x(), p.vect().y(), 0);
  Hep3Vector perp2(p_vect.x(), p_vect.y(), 0);
  Hep3Vector sperp2(zeta.x(), zeta.y(), 0);
  double distance2 = sqrt(radius * radius * (1 + tanlam * tanlam) * 
			  perp1.angle(perp2) * perp1.angle(perp2));
  return p_vect;
}
