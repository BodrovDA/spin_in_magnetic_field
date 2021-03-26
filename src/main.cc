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
#include "constants.h"
//#include "MichelFunctions.icc"
#include "GenFunctions.h"
#include "SpinParticleMagField.h"
#include "Event.h"

//using namespace std;

//const double bfieldc = 1.0 * 10E4;
//const double bfieldb = 1.5 * 10;
//const double balphac = 10000. / 3. / bfieldc;
//const double balphab = 10000. / 3. / bfieldb;

inline double gauss_width(double x, double width){
    return exp(- x * x / 2 / width / width);
}


int main(int argc, char **argv) {
  std::string file_name = "data/raw_belle_xip1_nB.root";
  TFile f(file_name.c_str(),"recreate");
  TTree t1("t1","tree with raw data from my MC");

  double elec = 8.0, posi = 3.5;
  HepLorentzVector beam(elec*sin(0.022), 0.,(elec*cos(0.022)-posi), (elec+posi));

  
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



    const int charge = 1;
    int j = 0;
    int k = 0;

    double max_width = 0;
    double max_width_mu = 0;

    //histogram for time dependence
    const double time_limit = 1.E-7;
    const double time_limit_lab = 1.E-8;

    for (int i = 0; i < 1.E5; ++i) {
      Event ev(beam, -1, generator);
      ev.generate_time(time_limit);
      double time = ev.time_mu();
      time_mu = time;
      
      if (i % 100000 == 0) std::cout << "Event number " << i << std::endl;
      ev.generate_tau_decay();
      
        //generate tau lepton
      //	ev.generate_tau();
	particle tau = ev.tau();
	HepLorentzVector p_tau_cms = tau.p();
	HepLorentzVector spin_tau = tau.s();

        
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

        HepLorentzVector p_tau = p_tau_cms;
        p_tau.boost(beam.boostVector());
        //~~~~~~~~~~~~~~~~~
        pl_x_tau = p_tau.vect().x();
        pl_y_tau = p_tau.vect().y();
        pl_z_tau = p_tau.vect().z();
        el_tau = p_tau.e();
        //~~~~~~~~~~~~~~~~~

        //generate muon
	//ev.generate_mu();
        particle mu = ev.mu();
	HepLorentzVector p_mu_tau = mu.p();
	HepLorentzVector spin_mu = mu.s();
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


        //calculate tau decay width
        double cost_tau = cos(p_mu_tau.vect().angle(spin_tau.vect()));
        double sint_tau = sqrt(1 - cost_tau * cost_tau);
        double x_mu = p_mu_tau.e() / w_mu;
        double x0_mu = m_mu / w_mu;

	cos_p_mu_s_mu = cos(spin_mu.vect().angle(p_mu_tau.vect()));

        HepLorentzVector spin_mu_tau(spin_mu.vect() + p_mu_tau.vect().dot(spin_mu.vect()) / m_mu /
				     (m_mu + p_mu_tau.e()) * p_mu_tau.vect(), p_mu_tau.vect().dot(spin_mu.vect()) / m_mu);
        double scal_mu1 = p_mu_tau.vect().dot(spin_mu.vect()) / p_mu_tau.vect().mag();
        double scal_mu2 = p_mu_tau.vect().dot(spin_mu_tau.vect()) 
	  / p_mu_tau.vect().mag() / spin_mu_tau.vect().mag();
        //~~~~~~~~~~~~~~~~~
        st_in_x_mu = spin_mu_tau.vect().x();
        st_in_y_mu = spin_mu_tau.vect().y();
        st_in_z_mu = spin_mu_tau.vect().z();
        st_in_t_mu = spin_mu_tau.t();
        //~~~~~~~~~~~~~~~~~


	//	std::cout << Fas(x_mu, x0_mu) << std::endl;
	//double width_tau = ev.width_tau();
/* +
	     (-Fip(x_mu, x0_mu)  + Fap(x_mu, x0_mu) * cost_tau) * spin_mu.vect().dot(axis_z) 
	     + Ft1(x_mu, x0_mu) * sint_tau * spin_mu.vect().dot(axis_x)
	     + Ft2(x_mu, x0_mu) * sint_tau * spin_mu.vect().dot(axis_y));*/

        //generate decay spectrum
        //if (width_tau > max_width) max_width = width_tau;
        //double rand_t = doubleRand();
        //if ((rand_t < width_tau)) {
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

		//}
        }
    }
    std::cout << max_width << " " << max_width_mu << std::endl;

    t1.Write();
    return 0;
}
