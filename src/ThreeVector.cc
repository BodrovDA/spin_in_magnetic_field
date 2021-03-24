// -*- C++ -*-
// $Id: ThreeVector.cc 10002 2007-02-26 06:56:17Z katayama $
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the implementation of the Hep3Vector class.
//


#include <iostream>

#include "ThreeVector.h"
#include "Rotation.h"




double Hep3Vector::operator () (int i) const {
  if (i == 0) {
    return x();
  }else if (i == 1) {
    return y();
  }else if (i == 2) {
    return z();
  }else{

    return 0.0;
  }
}

Hep3Vector & Hep3Vector::operator *= (const HepRotation & m){
  return *this = m * (*this);
}

Hep3Vector & Hep3Vector::transform(const HepRotation & m) {
  return *this = m * (*this);
}

void Hep3Vector::rotateX(double angle) {
  double s = sin(angle);
  double c = cos(angle);
  double yy = dy;
  dy = c*yy - s*dz;
  dz = s*yy + c*dz;
}

void Hep3Vector::rotateY(double angle) {
  double s = sin(angle);
  double c = cos(angle);
  double zz = dz;
  dz = c*zz - s*dx;
  dx = s*zz + c*dx;
}

void Hep3Vector::rotateZ(double angle) {
  double s = sin(angle);
  double c = cos(angle);
  double xx = dx;
  dx = c*xx - s*dy;
  dy = s*xx + c*dy;
}

void Hep3Vector::rotate(double angle, const Hep3Vector & axis){
  HepRotation trans;
  trans.rotate(angle, axis);
  operator*=(trans);
}

void Hep3Vector::rotateUz(Hep3Vector& NewUzVector){
  // NewUzVector must be normalized !

  double u1 = NewUzVector.x();
  double u2 = NewUzVector.y();
  double u3 = NewUzVector.z();
  double up = u1*u1 + u2*u2;

  if (up) {
      up = sqrt(up);
      double px = dx,  py = dy,  pz = dz;
      dx = (u1*u3*px - u2*py + u1*up*pz)/up;
      dy = (u2*u3*px + u1*py + u2*up*pz)/up;
      dz = (u3*u3*px -    px + u3*up*pz)/up;
    }
  else if (u3 < 0.) { dx = -dx; dz = -dz; }      // phi=0  teta=pi
  else {};
}

std::ostream & operator << (std::ostream & s, const Hep3Vector & q) {
  return s << "(" << q.x() << "," << q.y() << "," << q.z() << ")";
}

const Hep3Vector HepXHat(1.0, 0.0, 0.0);
const Hep3Vector HepYHat(0.0, 1.0, 0.0);
const Hep3Vector HepZHat(0.0, 0.0, 1.0);

const Hep3Vector xhat(1.0, 0.0, 0.0);
const Hep3Vector yhat(0.0, 1.0, 0.0);
const Hep3Vector zhat(0.0, 0.0, 1.0);