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

// include CLHep
#include "Particle.h"
#include "LorentzVector.h"
#include "ThreeVector.h"

// include custom libraries
#include "constants.h"
#include "GenFunctions.h"
#include "SpinParticleMagField.h"
#include "Event.h"

double xi_p;
double xi_pp;

int main(int argc, char **argv) {
  
  if(argc != 4) return -1;
  
  xi_p = atof(argv[1]);
  xi_pp = xi * xi_p;

  int flag = atoi(argv[2]);
  int num_of_events = atoi(argv[3]);

  std::string file_name = "data/raw_belle_xip";
  if(flag) {
    std::cout << "Start program with xi_p = " << xi_p << ", and magnetic field on "  << std::endl;
    file_name = file_name + std::string(argv[1]) + "_B.root";
  } else {
    std::cout << "Start program with xi_p = " << xi_p << ", and magnetic field off " << std::endl;
    file_name = file_name + std::string(argv[1]) + "_nB.root";
  }
      
  // create file
  TFile f(file_name.c_str(),"recreate");
  TTree t1("t1","tree with raw data from my MC");

  // set beam parameters
  double elec = 8.0, posi = 3.5;
  HepLorentzVector beam(elec*sin(0.022), 0.,(elec*cos(0.022)-posi), (elec+posi));

  
  std::default_random_engine generator;
  srand(time(NULL));

  // prepare tree variables
  // event general variables
  double time_mu, timel_mu, rotation_angle;
  int ev_number = 0;
  // beam variables
  double beam4[4];
  // tau variables
  double pc_tau[4], pl_tau[4], s_tau[4];
  // mu variables
  double pt_mu[4], pc_mu[4], pl_mu[4], s_mu[4];
  double pt_mu_out[4], pc_mu_out[4], pl_mu_out[4], s_mu_out[4];
  // electron variables
  double pm_e[4], pt_e[4], pc_e[4], pl_e[4], s_e[4];

  // general variables branches
  t1.Branch("time_mu",&time_mu,"time_mu/D");//
  t1.Branch("timel_mu",&timel_mu,"timel_mu/D");//
  t1.Branch("rotation_angle",&rotation_angle,"rotation_angle/D");//
  t1.Branch("ev_number",&ev_number,"ev_number/I");//

  // beam branches
  t1.Branch("beam",&beam4,"beam4[4]/D");//      

  // tau branches
  t1.Branch("pc_tau",&pc_tau,"pc_tau[4]/D");//    
  t1.Branch("pl_tau",&pl_tau,"pl_tau[4]/D");//
  t1.Branch("s_tau",&s_tau,"s_tau[4]/D");//
  // mu branches
  t1.Branch("pt_mu",&pt_mu,"pt_mu[4]/D");//    
  t1.Branch("pc_mu",&pc_mu,"pc_mu[4]/D");//    
  t1.Branch("pl_mu",&pl_mu,"pl_mu[4]/D");//
  t1.Branch("s_mu",&s_mu,"s_mu[4]/D");//
  t1.Branch("pt_mu_out",&pt_mu_out,"pt_mu_out[4]/D");//    
  t1.Branch("pc_mu_out",&pc_mu_out,"pc_mu_out[4]/D");//    
  t1.Branch("pl_mu_out",&pl_mu_out,"pl_mu_out[4]/D");//
  t1.Branch("s_mu_out",&s_mu_out,"s_mu_out[4]/D");//
  // e branches
  t1.Branch("pm_e",&pm_e,"pm_e[4]/D");//    
  t1.Branch("pt_e",&pt_e,"pt_e[4]/D");//    
  t1.Branch("pc_e",&pc_e,"pc_e[4]/D");//    
  t1.Branch("pl_e",&pl_e,"pl_e[4]/D");//
  t1.Branch("s_e",&s_e,"s_e[4]/D");//
  
  const int charge = -1;  
  const double time_limit = 1.E-8;
  const double time_limit_lab = 1.E-8;
  
  for (int i = 0; i < num_of_events; ++i) {
    Event ev(beam, charge, time_limit, time_limit_lab, generator);
    ev.generate_cascade(flag);

    time_mu = ev.time_mu();
    timel_mu = ev.time_lab_mu();
    rotation_angle = ev.rotation_angle();

    //~~~~~~~~~~~~~~~~~
    beam4[0] = beam.x();
    beam4[1] = beam.y();
    beam4[2] = beam.z();
    beam4[3] = beam.e();
    //~~~~~~~~~~~~~~~~~

      
    particle tau = ev.tau();
    HepLorentzVector p_tau_cms = tau.p();
    HepLorentzVector spin_tau = tau.s();
    //~~~~~~~~~~~~~~~~~
    pc_tau[0] = p_tau_cms.x();
    pc_tau[1] = p_tau_cms.y();
    pc_tau[2] = p_tau_cms.z();
    pc_tau[3] = p_tau_cms.e();
    s_tau[0] = spin_tau.x();
    s_tau[1] = spin_tau.y();
    s_tau[2] = spin_tau.z();
    s_tau[3] = spin_tau.t();
    //~~~~~~~~~~~~~~~~~

    HepLorentzVector p_tau_lab = p_tau_cms;
    p_tau_lab.boost(beam.boostVector());
    //~~~~~~~~~~~~~~~~~ 
    pl_tau[0] = p_tau_lab.x();
    pl_tau[1] = p_tau_lab.y();
    pl_tau[2] = p_tau_lab.z();
    pl_tau[3] = p_tau_lab.e();
    //~~~~~~~~~~~~~~~~~

    particle mu = ev.mu();
    HepLorentzVector p_mu_tau = mu.p();
    HepLorentzVector spin_mu = mu.s();
    //~~~~~~~~~~~~~~~~~
    pt_mu[0] = p_mu_tau.x();
    pt_mu[1] = p_mu_tau.y();
    pt_mu[2] = p_mu_tau.z();
    pt_mu[3] = p_mu_tau.e();
    s_mu[0] = spin_mu.x();
    s_mu[1] = spin_mu.y();
    s_mu[2] = spin_mu.z();
    s_mu[3] = spin_mu.t();
    //~~~~~~~~~~~~~~~~~

    HepLorentzVector p_mu_cms = p_mu_tau;
    p_mu_cms.boost(p_tau_cms.boostVector());
    //~~~~~~~~~~~~~~~~~
    pc_mu[0] = p_mu_cms.x();
    pc_mu[1] = p_mu_cms.y();
    pc_mu[2] = p_mu_cms.z();
    pc_mu[3] = p_mu_cms.e();
    //~~~~~~~~~~~~~~~~~

    HepLorentzVector p_mu_lab = p_mu_cms;
    p_mu_lab.boost(beam.boostVector());
    //~~~~~~~~~~~~~~~~~
    pl_mu[0] = p_mu_lab.x();
    pl_mu[1] = p_mu_lab.y();
    pl_mu[2] = p_mu_lab.z();
    pl_mu[3] = p_mu_lab.e();
    //~~~~~~~~~~~~~~~~~

    particle mu_out = ev.mu_out();
    HepLorentzVector p_mu_lab_out = mu_out.p();
    HepLorentzVector spin_mu_out = mu_out.s();
    //~~~~~~~~~~~~~~~~~
    pl_mu_out[0] = p_mu_lab_out.x();
    pl_mu_out[1] = p_mu_lab_out.y();
    pl_mu_out[2] = p_mu_lab_out.z();
    pl_mu_out[3] = p_mu_lab_out.e();
    s_mu_out[0] = spin_mu_out.x();
    s_mu_out[1] = spin_mu_out.y();
    s_mu_out[2] = spin_mu_out.z();
    s_mu_out[3] = spin_mu_out.t();
    //~~~~~~~~~~~~~~~~~

    HepLorentzVector p_mu_cms_out = p_mu_lab_out;
    p_mu_cms_out.boost(-beam.boostVector());
    //~~~~~~~~~~~~~~~~~
    pc_mu_out[0] = p_mu_cms_out.x();
    pc_mu_out[1] = p_mu_cms_out.y();
    pc_mu_out[2] = p_mu_cms_out.z();
    pc_mu_out[3] = p_mu_cms_out.e();
    //~~~~~~~~~~~~~~~~~

    HepLorentzVector p_mu_tau_out = p_mu_cms_out;
    p_mu_tau_out.boost(-p_tau_cms.boostVector());
    //~~~~~~~~~~~~~~~~~
    pt_mu_out[0] = p_mu_tau_out.x();
    pt_mu_out[1] = p_mu_tau_out.y();
    pt_mu_out[2] = p_mu_tau_out.z();
    pt_mu_out[3] = p_mu_tau_out.e();
    //~~~~~~~~~~~~~~~~~

    particle e = ev.e();
    HepLorentzVector p_e_mu = e.p();
    HepLorentzVector spin_e = e.s();
    //~~~~~~~~~~~~~~~~~
    pm_e[0] = p_e_mu.x();
    pm_e[1] = p_e_mu.y();
    pm_e[2] = p_e_mu.z();
    pm_e[3] = p_e_mu.e();
    s_e[0] = spin_e.x();
    s_e[1] = spin_e.y();
    s_e[2] = spin_e.z();
    s_e[3] = spin_e.t();
    //~~~~~~~~~~~~~~~~~

    HepLorentzVector p_e_tau = p_e_mu;
    p_e_tau.boost(p_mu_tau.boostVector());
    //~~~~~~~~~~~~~~~~~
    pt_e[0] = p_e_tau.x();
    pt_e[1] = p_e_tau.y();
    pt_e[2] = p_e_tau.z();
    pt_e[3] = p_e_tau.e();
    //~~~~~~~~~~~~~~~~~

    HepLorentzVector p_e_cms = p_e_tau;
    p_e_cms.boost(p_tau_cms.boostVector());
    //~~~~~~~~~~~~~~~~~
    pc_e[0] = p_e_cms.x();
    pc_e[1] = p_e_cms.y();
    pc_e[2] = p_e_cms.z();
    pc_e[3] = p_e_cms.e();
    //~~~~~~~~~~~~~~~~~

    HepLorentzVector p_e_lab = p_e_cms;
    p_e_lab.boost(beam.boostVector());
    //~~~~~~~~~~~~~~~~~
    pl_e[0] = p_e_lab.x();
    pl_e[1] = p_e_lab.y();
    pl_e[2] = p_e_lab.z();
    pl_e[3] = p_e_lab.e();
    //~~~~~~~~~~~~~~~~~


    ++ev_number;
    if(ev_number % 10000 == 0)
      std::cout << "Event number: " << ev_number << std::endl;
    t1.Fill();
    
		
  }
  
  t1.Write();
  return 0;
}
