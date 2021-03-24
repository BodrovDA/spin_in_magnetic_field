//
// Created by Denis Bodrov on 2020-01-23.
//
#ifndef THEORY_PARTICLE_H
#define THEORY_PARTICLE_H
#include <string>
#include "LorentzVector.h"
#include "ThreeVector.h"

class particle {
public:
    particle(HepLorentzVector p, HepLorentzVector s) : p_(p), s_(s) {}
    ~particle(){}
    HepLorentzVector p(){
        return p_;
    }
    HepLorentzVector s(){
        return s_;
    }
    void setp(HepLorentzVector p){
        p_ = p;
    }
    void sets(HepLorentzVector s){
        s_ = s;
    }

private:
    HepLorentzVector p_;
    HepLorentzVector s_;
};

#endif //THEORY_PARTICLE_H
