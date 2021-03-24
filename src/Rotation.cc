// -*- C++ -*-
// $Id: Rotation.cc 10002 2007-02-26 06:56:17Z katayama $
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the implementation of the HepRotation class.
//



#include <iostream>
#include "Rotation.h"

double HepRotation::operator () (int i, int j) const {
  if (i == 0) {
    if (j == 0) { return xx(); }
    if (j == 1) { return xy(); }
    if (j == 2) { return xz(); } 
  } else if (i == 1) {
    if (j == 0) { return yx(); }
    if (j == 1) { return yy(); }
    if (j == 2) { return yz(); } 
  } else if (i == 2) {
    if (j == 0) { return zx(); }
    if (j == 1) { return zy(); }
    if (j == 2) { return zz(); } 
  }
  return 0.0;
} 

HepRotation HepRotation::operator * (const HepRotation & b) const {
  return HepRotation(xx()*b.xx() + xy()*b.yx() + xz()*b.zx(),
		     xx()*b.xy() + xy()*b.yy() + xz()*b.zy(),
		     xx()*b.xz() + xy()*b.yz() + xz()*b.zz(),
		     yx()*b.xx() + yy()*b.yx() + yz()*b.zx(),
		     yx()*b.xy() + yy()*b.yy() + yz()*b.zy(),
		     yx()*b.xz() + yy()*b.yz() + yz()*b.zz(),
		     zx()*b.xx() + zy()*b.yx() + zz()*b.zx(),
		     zx()*b.xy() + zy()*b.yy() + zz()*b.zy(),
		     zx()*b.xz() + zy()*b.yz() + zz()*b.zz());
}

HepRotation & HepRotation::rotate(double psi, const Hep3Vector & axis) {
  HepRotation m;
  double phi = axis.phi();
  double theta = axis.theta();
  // First rotate so the axis coinsides with the z-axis
  m.rotateZ(-phi);
  m.rotateY(-theta);

  // Then rotate around the z-axis
  m.rotateZ(psi);

  // Finally rotate back to the original frame
  m.rotateY(theta);
  m.rotateZ(phi);
  transform(m);
  return *this;
}

HepRotation & HepRotation::rotateX(double psi) {
  double cp=cos(psi);
  double sp=sin(psi);
  HepRotation m(1.0, 0.0, 0.0, 0.0, cp, -sp, 0.0, sp, cp);
  transform(m);
  return *this;
}

HepRotation & HepRotation::rotateY(double theta){
  double st = sin(theta);
  double ct = cos(theta);
  HepRotation m(ct, 0.0, st, 0.0, 1.0, 0.0, -st, 0.0, ct);
  transform(m);
  return *this;
}

HepRotation & HepRotation::rotateZ(double phi) {
  double cp=cos(phi);
  double sp=sin(phi);
  HepRotation m(cp, -sp, 0.0, sp, cp, 0.0, 0.0, 0.0, 1.0);
  transform(m);
  return *this;
}

HepRotation & HepRotation::rotateAxes(const Hep3Vector &newX,
                                             const Hep3Vector &newY,
                                             const Hep3Vector &newZ) {
  double del = 0.001;
  Hep3Vector w = newX.cross(newY);

  if (abs(newZ.x()-w.x()) > del ||
      abs(newZ.y()-w.y()) > del ||
      abs(newZ.z()-w.z()) > del ||
      abs(newX.mag2()-1.) > del ||
      abs(newY.mag2()-1.) > del || 
      abs(newZ.mag2()-1.) > del ||
      abs(newX.dot(newY)) > del ||
      abs(newY.dot(newZ)) > del ||
      abs(newZ.dot(newX)) > del) {

    return *this;
  }else{
    return transform(HepRotation(newX.x(), newY.x(), newZ.x(),
                                 newX.y(), newY.y(), newZ.y(),
                                 newX.z(), newY.z(), newZ.z()));
  }
}

double HepRotation::phiX() const {
  return (yx() == 0.0 && xx() == 0.0) ? 0.0 : atan2(yx(),xx());
}

double HepRotation::phiY() const {
  return (yy() == 0.0 && xy() == 0.0) ? 0.0 : atan2(yy(),xy());
}

double HepRotation::phiZ() const {
  return (yz() == 0.0 && xz() == 0.0) ? 0.0 : atan2(yz(),xz());
}

double HepRotation::thetaX() const {
  return acos(zx());
}

double HepRotation::thetaY() const {
  return acos(zy());
}

double HepRotation::thetaZ() const {
  return acos(zz());
}

void HepRotation::getAngleAxis(double &angle, Hep3Vector &axis) const {
  double cosa  = 0.5*(xx()+yy()+zz()-1);
  double cosa1 = 1-cosa;
  if (cosa1 <= 0) {
    angle = 0;
    axis  = Hep3Vector(0,0,1);
  }else{
    double x=0, y=0, z=0;
    if (xx() > cosa) x = sqrt((xx()-cosa)/cosa1);
    if (yy() > cosa) y = sqrt((yy()-cosa)/cosa1);
    if (zz() > cosa) z = sqrt((zz()-cosa)/cosa1);
    if (zy() < yz())  x = -x;
    if (xz() < zx())  y = -y;
    if (yx() < xy())  z = -z;
    angle = acos(cosa);
    axis  = Hep3Vector(x,y,z);
  }
}

