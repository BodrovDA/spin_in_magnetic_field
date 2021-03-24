// -*- C++ -*-
// $Id: LorentzVector.cc 10002 2007-02-26 06:56:17Z katayama $
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the implementation of the HepLorentzVector class.
//

#include <iostream>
#include "LorentzVector.h"
#include "Rotation.h"
#include "LorentzRotation.h"

HepLorentzVector &
HepLorentzVector::operator *= (const HepRotation & m) {
  pp *= m;
  return *this;
}

HepLorentzVector &
HepLorentzVector::transform(const HepRotation & m) {
  pp.transform(m);
  return *this;
}


void HepLorentzVector::boost(double bx, double by, double bz){
  double b2 = bx*bx + by*by + bz*bz;
  register double gamma = 1.0 / sqrt(1.0 - b2);
  register double bp = bx*x() + by*y() + bz*z();
  register double gamma2 = b2 > 0 ? (gamma - 1.0)/b2 : 0.0;

  setX(x() + gamma2*bp*bx + gamma*bx*t());
  setY(y() + gamma2*bp*by + gamma*by*t());
  setZ(z() + gamma2*bp*bz + gamma*bz*t());
  setT(gamma*(t() + bp));
}

HepLorentzVector &
HepLorentzVector::operator *= (const HepLorentzRotation & m) {
  return *this = m.vectorMultiplication(*this);
}

HepLorentzVector &
HepLorentzVector::transform(const HepLorentzRotation & m){
  return *this = m.vectorMultiplication(*this);
}

std::ostream & operator << (std::ostream & s, const HepLorentzVector & q)
{
  return s << "(" << q.x() << "," << q.y() << ","
           << q.z() << ";" << q.t() << ")";
}

