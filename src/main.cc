//
// Created by Denis Bodrov on 2020-01-24.
//

// include standard c++ libraries
#include <string>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <random>

// include cern root libraries
#include <TH1F.h>
#include <TF1.h>
#include <TH2F.h>
#include <TApplication.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
//#include <TBranch.h>


// include CLHep
#include "Particle.h"
#include "LorentzVector.h"
#include "ThreeVector.h"

// include custom libraries
//#include "constants.h"
//#include "MichelFunctions.icc"
#include "GenFunctions.h"

//using namespace std;

const double const_pi = 3.1415926;
const double m_tau = 1.7769;
const double m_mu = 0.1057; //0.00051;//0.1057;
const double m_e = 0.00051;
const double w_mu = (m_tau * m_tau + m_mu * m_mu) / 2 / m_tau;
const double w_e = (m_mu * m_mu + m_e * m_e) / 2 / m_mu;
const double mu_energy_min = m_mu;
const double e_energy_min = m_e;

//const double bfieldc = 1.0 * 10E4;
//const double bfieldb = 1.5 * 10;
//const double balphac = 10000. / 3. / bfieldc;
//const double balphab = 10000. / 3. / bfieldb;
const Hep3Vector BfieldC(0., 0., 1.);
const Hep3Vector BfieldB(0., 0., 1.5);
const double dt = 0.5 / 3. / 1.E10;
const double magnb = 14E9 * m_e / m_mu;

//Measured Michel parameters
const double rho = 0.75;
const double delta = 0.75;
const double xi = 1.;
const double eta = 0;

// Not measured Michel parameters 
const double rho_p = 3. / 4;
const double delta_p = 1. / 4;
double xi_p = 1;
const double xi_pp = xi * xi_p;//1;
const double eta_p = 0;
const double eta_pp = 0;
const double alpha_p = 0;
const double alpha = 0;
const double beta = 0;
const double beta_p = 0;
const double a = 0;
const double b = 4;
const double c = 0;
const double A = 16;
const double decay_time = 2.2E-6;


inline double Fis(double x, double x0){
    return x*(1 - x) + 2./9. * rho * (4*x*x - 3*x - x0*x0) + eta*x0*(1 - x);
}
inline double Fas(double x, double x0){
  return (1. / 3.) * xi * sqrt(x*x - x0*x0) * ((1 - x) + 2. / 3. * delta * (4 * x - 4 + sqrt(1 - x0 * x0)));
}
inline double Fip(double x, double x0){
    return 1./54 * sqrt(x*x - x0*x0) * (9*xi_p*(-2*x + 2 + sqrt(1 - x0*x0)) + 4*xi*(delta - 3./4) 
					* (4*x - 4 + sqrt(1 - x0*x0)));
}
inline double Fap(double x, double x0){
      return 1./6 * (xi_pp*(2*x*x - x - x0*x0) + 4*(rho - 3./4)*(4*x*x - 3*x - x0*x0) + 
		     2*eta_pp*(1 - x)*x0);
}
inline double Ft1(double x, double x0){
  return 1./12 * (-2*(xi_pp + 12*(rho - 3./4))*(1 - x)*x0 - 3*eta*(x*x - x0*x0) + eta_pp*(-3*x*x + 4*x - x0*x0));
}
inline double Ft2(double x, double x0){
  return 1./3 * sqrt(x*x - x0*x0) * (3*(1 - x)*alpha_p/A + 2*sqrt(1 - x0*x0)*beta_p/A);
}

inline double gauss_width(double x, double width){
    return exp(- x * x / 2 / width / width);
}

Hep3Vector spin_rotation(Hep3Vector &zeta, HepLorentzVector p, double time_dec){
    Hep3Vector p_vect = p.vect();
    //std::cout << "yoshino " << cos(p.vect().angle(zeta));
    double kappa = 1 / p_vect.perp();
    double tanlam = p_vect.z() * kappa;
    double alpha = 10000. / 3. / 10. / BfieldB.mag();
    double radius = alpha / kappa;
    //std::cout << "homura, r: " << radius << " a: " << alpha << " k: " << kappa << " pz: " << p_vect.z() << " pt: " << p_vect.perp() << " p: " << p_vect.mag() << std::endl;
    double distance1 = 0;
    Hep3Vector sperp1(zeta.x(), zeta.y(), 0);
    double time = 10000 * dt;
   /* while(true){
        double tmp_time = decay_time * doubleRand();
        double rand = doubleRand();
        if(rand < exp(-tmp_time / decay_time * m_mu / p.e())){
            time = tmp_time;
            break;
        }
    }*/
    //time_dec = time;
    //int max_i = 0;
    for(int i = 0; i < time_dec / dt; ++i){
        distance1 += p_vect.mag() / p.e() * 3E10 * dt;
        Hep3Vector dzeta = 2 * magnb * m_mu / p.e() * zeta.cross(BfieldB) * dt * 6.30253;
        zeta += dzeta;
        Hep3Vector dp = 2 * magnb * m_mu / p.e() * p_vect.cross(BfieldB) * dt * 6.30253;
        p_vect += dp;
        //max_i = i;
        //std::cout << "ashigara " << dzeta.mag() / zeta.perp() << " " << dp.mag() / p_vect.perp() << " " << cos(p_vect.angle(zeta)) << std::endl;
    }
    Hep3Vector perp1(p.vect().x(), p.vect().y(), 0);
    Hep3Vector perp2(p_vect.x(), p_vect.y(), 0);
    Hep3Vector sperp2(zeta.x(), zeta.y(), 0);
    //    std::cout << "yugumo " << cos(perp1.angle(perp2)) << " " << cos(sperp1.angle(sperp2)) << std::endl;
    double distance2 = sqrt(radius * radius * (1 + tanlam * tanlam) * perp1.angle(perp2) * perp1.angle(perp2));
    //std::cout << "hitagi " << distance1 << " " << distance2 << " " << distance1 / distance2 << std::endl;
    //std::cout << "hanekawa " << max_i << std::endl;
    return p_vect;
}

int main(int argc, char **argv) {
  std::string file_name = "data/raw_belle_xip1_nB.root";
  TFile f(file_name.c_str(),"recreate");
    TTree t1("t1","tree with raw data from my MC");

    std::default_random_engine generator;
    /* initialize random seed: */
    srand(time(NULL));
    double cms_energy = 10.58; //3.56;//3.7;//10.58;
    double tau_energy = cms_energy / 2;
    double tau_p = sqrt(tau_energy * tau_energy - m_tau * m_tau);

    //prepare tree variables
    double bfield, time_mu, timel_mu;
    int ev_number = 0;
    //tau variables
    double p_x_tau, p_y_tau, p_z_tau, e_tau, pl_x_tau, pl_y_tau, pl_z_tau, el_tau;
    double s_x_tau, s_y_tau, s_z_tau, s_t_tau;
    //mu variables
    double p_x_mu, p_y_mu, p_z_mu, e_mu, pl_in_x_mu, pl_in_y_mu, pl_in_z_mu, el_in_mu;
    double pl_out_x_mu, pl_out_y_mu, pl_out_z_mu, el_out_mu;
    double s_in_x_mu, s_in_y_mu, s_in_z_mu, s_in_t_mu;
    double s_out_x_mu, s_out_y_mu, s_out_z_mu, s_out_t_mu;
    double st_in_x_mu, st_in_y_mu, st_in_z_mu, st_in_t_mu;
    //electron variables
    double p_x_e, p_y_e, p_z_e, e_e, pl_x_e, pl_y_e, pl_z_e, el_e;
    double s_x_e, s_y_e, s_z_e, s_t_e;
    double sm_x_e, sm_y_e, sm_z_e, sm_t_e;
    double cos_p_mu_s_tau, cos_p_e_s_mu, helicity_mu, helicity_e, cos_p_mu_s_mu, cos_p_e_p_mu, spin_mu_mag;
    double cos_s_mu_PL_mu, PL_mu_mag, test;

    //prepare branches
    t1.Branch("helicity_e",&helicity_e,"helicity_e/D");//
    t1.Branch("test",&test,"test/D");//
    t1.Branch("helicity_mu",&helicity_mu,"helicity_mu/D");//
    t1.Branch("spin_mu_mag",&spin_mu_mag,"spin_mu_mag/D");//
    t1.Branch("PL_mu_mag",&PL_mu_mag,"PL_mu_mag/D");//
    t1.Branch("cos_s_mu_PL_mu",&cos_s_mu_PL_mu,"cos_s_mu_PL_mu/D");//
    //    t1.Branch("spin_mu_mag",&spin_mu_mag,"spin_mu_mag/D");//
    t1.Branch("cms_energy",&cms_energy,"cms_energy/D");//
    t1.Branch("tau_energy",&tau_energy,"tau_energy/D");//
    t1.Branch("tau_p",&tau_p,"tau_p/D");//
    t1.Branch("bfield",&bfield,"bfield/D");
    t1.Branch("time_mu",&time_mu,"time_mu/D");//
    t1.Branch("timel_mu",&timel_mu,"timel_mu/D");//
    t1.Branch("ev_number",&ev_number,"ev_number/I");//

    t1.Branch("p_x_tau",&p_x_tau,"p_x_tau/D");//
    t1.Branch("p_y_tau",&p_y_tau,"p_y_tau/D");//
    t1.Branch("p_z_tau",&p_z_tau,"p_z_tau/D");//
    t1.Branch("e_tau",&e_tau,"e_tau/D");//
    t1.Branch("pl_x_tau",&pl_x_tau,"pl_x_tau/D");//
    t1.Branch("pl_y_tau",&pl_y_tau,"pl_y_tau/D");//
    t1.Branch("pl_z_tau",&pl_z_tau,"pl_z_tau/D");//
    t1.Branch("el_tau",&el_tau,"el_tau/D");//
    t1.Branch("s_x_tau",&s_x_tau,"s_x_tau/D");//
    t1.Branch("s_y_tau",&s_y_tau,"s_y_tau/D");//
    t1.Branch("s_z_tau",&s_z_tau,"s_z_tau/D");//
    t1.Branch("s_t_tau",&s_t_tau,"s_t_tau/D");//

    t1.Branch("p_x_mu",&p_x_mu,"p_x_mu/D");//
    t1.Branch("p_y_mu",&p_y_mu,"p_y_mu/D");//
    t1.Branch("p_z_mu",&p_z_mu,"p_z_mu/D");//
    t1.Branch("e_mu",&e_mu,"e_mu/D");//
    t1.Branch("cos_p_mu_s_tau",&cos_p_mu_s_tau,"cos_p_mu_s_tau/D");//
    t1.Branch("cos_p_mu_s_mu",&cos_p_mu_s_mu,"cos_p_mu_s_mu/D");//
    t1.Branch("cos_p_e_s_mu",&cos_p_e_s_mu,"cos_p_e_s_mu/D");//
    t1.Branch("cos_p_e_p_mu",&cos_p_e_p_mu,"cos_p_e_p_mu/D");//
    t1.Branch("pl_in_x_mu",&pl_in_x_mu,"pl_in_x_mu/D");//
    t1.Branch("pl_in_y_mu",&pl_in_y_mu,"pl_in_y_mu/D");//
    t1.Branch("pl_in_z_mu",&pl_in_z_mu,"pl_in_z_mu/D");//
    t1.Branch("el_in_mu",&el_in_mu,"el_in_mu/D");//
    t1.Branch("pl_out_x_mu",&pl_out_x_mu,"pl_out_x_mu/D");//
    t1.Branch("pl_out_y_mu",&pl_out_y_mu,"pl_out_y_mu/D");//
    t1.Branch("pl_out_z_mu",&pl_out_z_mu,"pl_out_z_mu/D");//
    t1.Branch("el_out_mu",&el_out_mu,"el_out_mu/D");//
    t1.Branch("s_in_x_mu",&s_in_x_mu,"s_in_x_mu/D");//
    t1.Branch("s_in_y_mu",&s_in_y_mu,"s_in_y_mu/D");//
    t1.Branch("s_in_z_mu",&s_in_z_mu,"s_in_z_mu/D");//
    t1.Branch("s_in_t_mu",&s_in_t_mu,"s_in_t_mu/D");//
    t1.Branch("s_out_x_mu",&s_out_x_mu,"s_out_x_mu/D");//
    t1.Branch("s_out_y_mu",&s_out_y_mu,"s_out_y_mu/D");//
    t1.Branch("s_out_z_mu",&s_out_z_mu,"s_out_z_mu/D");//
    t1.Branch("s_out_t_mu",&s_out_t_mu,"s_out_t_mu/D");//
    t1.Branch("st_in_x_mu",&st_in_x_mu,"st_in_x_mu/D");//
    t1.Branch("st_in_y_mu",&st_in_y_mu,"st_in_y_mu/D");//
    t1.Branch("st_in_z_mu",&st_in_z_mu,"st_in_z_mu/D");//
    t1.Branch("st_in_t_mu",&st_in_t_mu,"st_in_t_mu/D");//

    t1.Branch("p_x_e",&p_x_e,"p_x_e/D");//
    t1.Branch("p_y_e",&p_y_e,"p_y_e/D");//
    t1.Branch("p_z_e",&p_z_e,"p_z_e/D");//
    t1.Branch("e_e",&e_e,"e_e/D");//
    t1.Branch("pl_x_e",&pl_x_e,"pl_x_e/D");//
    t1.Branch("pl_y_e",&pl_y_e,"pl_y_e/D");//
    t1.Branch("pl_z_e",&pl_z_e,"pl_z_e/D");//
    t1.Branch("el_e",&el_e,"el_e/D");//
    t1.Branch("s_x_e",&s_x_e,"s_x_e/D");//
    t1.Branch("s_y_e",&s_y_e,"s_y_e/D");//
    t1.Branch("s_z_e",&s_z_e,"s_z_e/D");//
    t1.Branch("s_t_e",&s_t_e,"s_t_e/D");//
    t1.Branch("sm_x_e",&sm_x_e,"sm_x_e/D");//
    t1.Branch("sm_y_e",&sm_y_e,"sm_y_e/D");//
    t1.Branch("sm_z_e",&sm_z_e,"sm_z_e/D");//
    t1.Branch("sm_t_e",&sm_t_e,"sm_t_e/D");//


    //std::cout << Fis(1, 0) << " "  << Fas(1, 0) << std::endl; 
    //std::cout << Fis(0.118, 0.118) << " "  << Fas(0.118, 0.118) << std::endl; 
    //std::cout << Fis(1, 0.118) << " "  << Fas(1, 0.118) << std::endl; 

    const int charge = 1;
    int j = 0;
    int k = 0;

    double max_width = 0;
    double max_width_mu = 0;

    //histogram for time dependence
    const double time_limit = 1.E-6;
    const double time_limit_lab = 1.E-8;

    for (int i = 0; i < 1.E7; ++i) {
        double time = doubleRand() * time_limit;
        double func = exp(-time / decay_time);// / decay_time;
        double func_rand = doubleRand();// / decay_time;
        if(func_rand > func){
            continue;
        }
        time_mu = time;

        if (i % 100000 == 0) std::cout << "Event number " << i << std::endl;
        //generate tau lepton
        HepLorentzVector p_tau_cms = generate_lorentz(tau_p, tau_energy, generator);
        HepLorentzVector spin_tau = generate_lorentz(1, 0, generator);
        particle tau(p_tau_cms, spin_tau);
        
	//~~~~~~~~~~~~~~~~~
        p_x_tau = p_tau_cms.vect().x();
        p_y_tau = p_tau_cms.vect().y();
        p_z_tau = p_tau_cms.vect().z();
        e_tau = p_tau_cms.e();
        s_x_tau = spin_tau.vect().x();
        s_y_tau = spin_tau.vect().y();
        s_z_tau = spin_tau.vect().z();
        s_t_tau = spin_tau.t();
        //~~~~~~~~~~~~~~~~~

	double elec = 8.0, posi = 3.5;
	HepLorentzVector beam(elec*sin(0.022)/2., 0.,(elec*cos(0.022)-posi)/2., (elec+posi)/2.);
        //HepLorentzVector beam(0., 0., 0., (elec+posi)/2.);
        //HepLorentzVector beam(0., 0., 0., cms_energy);
        HepLorentzVector p_tau = p_tau_cms;
        p_tau.boost(beam.boostVector());
        //~~~~~~~~~~~~~~~~~
        pl_x_tau = p_tau.vect().x();
        pl_y_tau = p_tau.vect().y();
        pl_z_tau = p_tau.vect().z();
        el_tau = p_tau.e();
        //~~~~~~~~~~~~~~~~~

        //generate muon
        double mu_energy = mu_energy_min + doubleRand() * (w_mu - mu_energy_min);
        double mu_p = sqrt(mu_energy * mu_energy - m_mu * m_mu);
        HepLorentzVector p_mu_tau = generate_lorentz(mu_p, mu_energy, generator);
        //~~~~~~~~~~~~~~~~~
        p_x_mu = p_mu_tau.vect().x();
        p_y_mu = p_mu_tau.vect().y();
        p_z_mu = p_mu_tau.vect().z();
        e_mu = p_mu_tau.e();
        //s_x_mu = spin_mu.vect().x();
        //s_y_mu = spin_mu.vect().y();
        //s_z_mu = spin_mu.vect().z();
        //s_t_mu = spin_mu.t();
        //~~~~~~~~~~~~~~~~~

        Hep3Vector axis_z = p_mu_tau.vect( ) * (1. / p_mu_tau.vect().mag());
        Hep3Vector axis_y = p_mu_tau.vect().cross(spin_tau.vect()) * (1. / p_mu_tau.vect().cross(spin_tau.vect()).mag());
        Hep3Vector axis_x = axis_y.cross(axis_z);

        
        //double width_t = (q_t.mag2() * p_t.dot(k_t) + 2 * q_t.dot(p_t) * q_t.dot(k_t)) * p_mu_tau.vect().mag2()
        //        / p_mu_tau.e() / (m_tau * m_tau * m_tau * m_tau * m_tau);

        //calculate tau decay width
        double cost_tau = cos(p_mu_tau.vect().angle(spin_tau.vect()));
        double sint_tau = sqrt(1 - cost_tau * cost_tau);
        double x_mu = mu_energy / w_mu;
        double x0_mu = m_mu / w_mu;


	double spin_denominator = Fis(x_mu, x0_mu) - Fas(x_mu, x0_mu) * cost_tau;
	double spin_mu_z = (-Fip(x_mu, x0_mu)  + Fap(x_mu, x0_mu) * cost_tau) / spin_denominator;
	double spin_mu_x = Ft1(x_mu, x0_mu) * sint_tau / spin_denominator;
	double spin_mu_y = Ft2(x_mu, x0_mu) * sint_tau / spin_denominator;
	//	double spin_mu_mag = sqrt(spin_mu_x * spin_mu_x + spin_mu_y * spin_mu_y + spin_mu_z * spin_mu_z);

	Hep3Vector PL_mu_vect = spin_mu_x * axis_x + spin_mu_y * axis_y + spin_mu_z * axis_z;
	Hep3Vector r_tmp = PL_mu_vect;	
	
	Hep3Vector spin_mu_vect_tmp(0,0,0);

	if(PL_mu_vect.mag() != 0) {

	  while((r_tmp.angle(PL_mu_vect) < const_pi / 12) || (r_tmp.angle(PL_mu_vect) > 11 * const_pi / 12))
	    r_tmp = generate_3vector(1, generator);

	  Hep3Vector PL_mu_perp = r_tmp - (r_tmp.dot(PL_mu_vect) / PL_mu_vect.mag2()) * PL_mu_vect;
	  Hep3Vector spin_mu_perp = (sqrt(1 - PL_mu_vect.mag2()) / PL_mu_perp.mag()) * PL_mu_perp;
       
	  spin_mu_vect_tmp = PL_mu_vect + spin_mu_perp;
	} else {
	  spin_mu_vect_tmp = generate_3vector(1, generator);
	}
	spin_mu_mag = spin_mu_vect_tmp.mag();
	cos_s_mu_PL_mu = cos(spin_mu_vect_tmp.angle(PL_mu_vect));
	PL_mu_mag = PL_mu_vect.mag();

	//	test = cos(PL_mu_vect.angle(p_mu_tau.vect()));

	int spin_mu_projection = 1 - 2 * (rand() % 2);

	//	Hep3Vector spin_mu_vect_tmp = axis_x * spin_mu_projection * (spin_mu_x/spin_mu_mag) + 
	//  axis_y * spin_mu_projection * (spin_mu_y/spin_mu_mag) + axis_z * spin_mu_projection * (spin_mu_z/spin_mu_mag);

	/*	HepLorentzVector spin_mu(spin_mu_projection * p_mu_tau.x()/p_mu_tau.vect().mag(), 
				 spin_mu_projection * p_mu_tau.y()/p_mu_tau.vect().mag(), 
				 spin_mu_projection * p_mu_tau.z()/p_mu_tau.vect().mag(), 0);*/
	HepLorentzVector spin_mu(spin_mu_vect_tmp, 0);
	cos_p_mu_s_mu = cos(spin_mu.vect().angle(p_mu_tau));

	particle mu(p_mu_tau, spin_mu);
	//std::cout << spin_mu_projection << std::endl;
	helicity_mu = spin_mu_projection;

        HepLorentzVector spin_mu_tau(spin_mu.vect() + p_mu_tau.vect().dot(spin_mu.vect()) / m_mu /
        (m_mu + mu_energy) * p_mu_tau.vect(), p_mu_tau.vect().dot(spin_mu.vect()) / m_mu);
        double scal_mu1 = p_mu_tau.vect().dot(spin_mu.vect()) / p_mu_tau.vect().mag();
        double scal_mu2 = p_mu_tau.vect().dot(spin_mu_tau.vect()) 
	  / p_mu_tau.vect().mag() / spin_mu_tau.vect().mag();
        //~~~~~~~~~~~~~~~~~
        st_in_x_mu = spin_mu_tau.vect().x();
        st_in_y_mu = spin_mu_tau.vect().y();
        st_in_z_mu = spin_mu_tau.vect().z();
        st_in_t_mu = spin_mu_tau.t();
        //~~~~~~~~~~~~~~~~~


	
	double coeff_tau = w_mu*w_mu*w_mu*w_mu / (m_tau*m_tau*m_tau*m_tau) * 24 * 2;
        double width_tau = coeff_tau * sqrt(x_mu * x_mu - x0_mu * x0_mu) 
	  * (Fis(x_mu, x0_mu) - Fas(x_mu, x0_mu) * cost_tau);/* +
	     (-Fip(x_mu, x0_mu)  + Fap(x_mu, x0_mu) * cost_tau) * spin_mu.vect().dot(axis_z) 
	     + Ft1(x_mu, x0_mu) * sint_tau * spin_mu.vect().dot(axis_x)
	     + Ft2(x_mu, x0_mu) * sint_tau * spin_mu.vect().dot(axis_y));*/

        //generate decay spectrum
        if (width_tau > max_width) max_width = width_tau;
        double rand_t = doubleRand();
        if ((rand_t < width_tau)) {
	  cos_p_mu_s_tau = cost_tau;
	  
	  // **********
	  //muon parameters before rotation
            Hep3Vector spin_mu_vect = spin_mu.vect();
            HepLorentzVector p_mu_lab_in = p_mu_tau;
            p_mu_lab_in.boost(p_tau.boostVector());
            //~~~~~~~~~~~~~~~~~
            pl_in_x_mu = p_mu_lab_in.vect().x();
            pl_in_y_mu = p_mu_lab_in.vect().y();
            pl_in_z_mu = p_mu_lab_in.vect().z();
            el_in_mu = p_mu_lab_in.e();
            //~~~~~~~~~~~~~~~~~

            //muon rotation in magnetic field
            double time_dec = time * p_mu_lab_in.e() / p_mu_lab_in.mag() * 0;
            timel_mu = time_dec;
	    if(time_dec > time_limit_lab) continue;

            Hep3Vector p3_mu_lab_out = spin_rotation(spin_mu_vect, p_mu_lab_in, time_dec);
            HepLorentzVector p_mu_lab_out(p3_mu_lab_out, p_mu_lab_in.e());
	    test = cos(p3_mu_lab_out.angle(p_mu_lab_in));
	    //HepLorentzVector p_mu_lab_out = p_mu_lab_in;//(p3_mu_lab_out, p_mu_lab_in.e());
            //~~~~~~~~~~~~~~~~~
            pl_out_x_mu = p_mu_lab_out.vect().x();
            pl_out_y_mu = p_mu_lab_out.vect().y();
            pl_out_z_mu = p_mu_lab_out.vect().z();
            el_out_mu = p_mu_lab_out.e();
            //~~~~~~~~~~~~~~~~~

            //muon parameters after rotation
            HepLorentzVector p_mu_tau_out = p_mu_lab_out;
            p_mu_tau_out.boost(-p_tau.boostVector());
            // **********

            //generate electron
            double e_energy = e_energy_min + doubleRand() * (w_e - e_energy_min);
            double e_p = sqrt(e_energy * e_energy - m_e * m_e);
            HepLorentzVector p_e_mu = generate_lorentz(e_p, e_energy, generator);
	    int spin_e_projection = 1 - 2 * (rand() % 2);
            HepLorentzVector spin_e(spin_e_projection * p_e_mu.x()/p_e_mu.vect().mag(), 
				    spin_e_projection * p_e_mu.y()/p_e_mu.vect().mag(),
				    spin_e_projection * p_e_mu.z()/p_e_mu.vect().mag(), 0);
            particle e(p_e_mu, spin_e);
	    helicity_e = spin_e_projection;
            //~~~~~~~~~~~~~~~~~
            p_x_e = p_e_mu.vect().x();
            p_y_e = p_e_mu.vect().y();
            p_z_e = p_e_mu.vect().z();
            e_e = p_e_mu.e();
            s_x_e = spin_e.vect().x();
            s_y_e = spin_e.vect().y();
            s_z_e = spin_e.vect().z();
            s_t_e = spin_e.t();
            //~~~~~~~~~~~~~~~~~
            HepLorentzVector spin_mu_out(spin_mu_vect, spin_mu.t());
            //~~~~~~~~~~~~~~~~~
            s_out_x_mu = spin_mu_out.vect().x();
            s_out_y_mu = spin_mu_out.vect().y();
            s_out_z_mu = spin_mu_out.vect().z();
            s_out_t_mu = spin_mu_out.e();
            //~~~~~~~~~~~~~~~~~

            HepLorentzVector spin_e_mu(
                    spin_e.vect() + p_e_mu.vect().dot(spin_e.vect()) / m_e / (m_e + e_energy) * p_e_mu.vect(),
                    p_e_mu.vect().dot(spin_e.vect()) / m_e);
            //~~~~~~~~~~~~~~~~~
            sm_x_e = spin_e_mu.vect().x();
            sm_y_e = spin_e_mu.vect().y();
            sm_z_e = spin_e_mu.vect().z();
            sm_t_e = spin_e_mu.t();
            //~~~~~~~~~~~~~~~~~

            HepLorentzVector p_mu_rest(0, 0, 0, m_mu);
            HepLorentzVector mass_spin_mu(spin_mu_out.vect() * m_mu, spin_mu_out.t() * m_mu);
            HepLorentzVector p_m = p_mu_rest - mass_spin_mu;
            HepLorentzVector mass_spin_e_mu(spin_e_mu.vect() * m_e, spin_e_mu.t() * m_e);
            HepLorentzVector k_m = p_e_mu - mass_spin_e_mu;
            HepLorentzVector q_m = p_mu_rest - p_e_mu;

            //muon decay width calculation
            double cost_mu = cos(p_e_mu.vect().angle(spin_mu_out.vect()));
	    cos_p_e_s_mu = cost_mu;
            double scal_e1 = p_e_mu.vect().dot(spin_e.vect()) / p_e_mu.vect().mag();
            double scal_e2 = p_e_mu.vect().dot(spin_e_mu.vect()) / p_e_mu.vect().mag() / spin_e_mu.vect().mag();
            double x_e = e_energy / w_e;
            double x0_e = m_e / w_e;
            double width_mu = (q_m.mag2() * p_m.dot(k_m) + 2 * q_m.dot(p_m) * q_m.dot(k_m)) * p_e_mu.vect().mag2()
                    / p_e_mu.e() / (m_mu * m_mu * m_mu * m_mu * m_mu);

            //generate decay spectrum
            double rand_mu = doubleRand();
            if (width_mu > max_width_mu) max_width_mu = width_mu;
            if (rand_mu < width_mu) {

                ++j;
		cos_p_e_p_mu = cos(p_e_mu.vect().angle(p3_mu_lab_out));
                HepLorentzVector p_e_lab = p_e_mu;
                p_e_lab.boost(p_mu_lab_out.boostVector());
                //~~~~~~~~~~~~~~~~~
                pl_x_e = p_e_lab.vect().x();
                pl_y_e = p_e_lab.vect().y();
                pl_z_e = p_e_lab.vect().z();
                el_e = p_e_lab.e();
                //~~~~~~~~~~~~~~~~~
		++ev_number;
                if(ev_number % 10000 == 0)
                    std::cout << "ashigara " << ev_number << std::endl;
                t1.Fill();

	    }
        }
    }
    std::cout << max_width << " " << max_width_mu << std::endl;

    t1.Write();
    return 0;
}
