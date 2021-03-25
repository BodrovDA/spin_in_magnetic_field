#ifndef SPIN_PARTICLE_MAG_FIELD_H
#define SPIN_PARTICLE_MAG_FIELD_H

#include "LorentzVector.h"
#include "ThreeVector.h"
#include "constants.h"

Hep3Vector spin_rotation(Hep3Vector &zeta, HepLorentzVector p, double time_dec);

#endif
