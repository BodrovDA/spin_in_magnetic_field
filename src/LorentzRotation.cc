// -*- C++ -*-
// $Id: LorentzRotation.cc 10002 2007-02-26 06:56:17Z katayama $
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the implementation of the HepLorentzRotation class.

#include <iostream>
#include "LorentzRotation.h"

double HepLorentzRotation::operator () (int i, int j) const {
  if (i == 0) {
    if (j == 0) { return xx(); }
    if (j == 1) { return xy(); }
    if (j == 2) { return xz(); } 
    if (j == 3) { return xt(); } 
  } else if (i == 1) {
    if (j == 0) { return yx(); }
    if (j == 1) { return yy(); }
    if (j == 2) { return yz(); } 
    if (j == 3) { return yt(); } 
  } else if (i == 2) {
    if (j == 0) { return zx(); }
    if (j == 1) { return zy(); }
    if (j == 2) { return zz(); } 
    if (j == 3) { return zt(); } 
  } else if (i == 3) {
    if (j == 0) { return tx(); }
    if (j == 1) { return ty(); }
    if (j == 2) { return tz(); } 
    if (j == 3) { return tt(); } 
  } 

  return 0.0;
} 

void HepLorentzRotation::setBoost(double bx, double by, double bz) {
  double bp2 = bx*bx + by*by + bz*bz;
  double gamma = 1.0 / sqrt(1.0 - bp2);
  double bgamma = gamma * gamma / (1.0 + gamma);
  mxx = 1.0 + bgamma * bx * bx;
  myy = 1.0 + bgamma * by * by;
  mzz = 1.0 + bgamma * bz * bz;
  mxy = myx = bgamma * bx * by;
  mxz = mzx = bgamma * bx * bz;
  myz = mzy = bgamma * by * bz;
  mxt = mtx = gamma * bx;
  myt = mty = gamma * by;
  mzt = mtz = gamma * bz;
  mtt = gamma;
}

HepLorentzRotation
HepLorentzRotation::matrixMultiplication(const HepLorentzRotation & b) const {
  return HepLorentzRotation(
    xx()*b.xx() + xy()*b.yx() + xz()*b.zx() + xt()*b.tx(),
    xx()*b.xy() + xy()*b.yy() + xz()*b.zy() + xt()*b.ty(),
    xx()*b.xz() + xy()*b.yz() + xz()*b.zz() + xt()*b.tz(),
    xx()*b.xt() + xy()*b.yt() + xz()*b.zt() + xt()*b.tt(),

    yx()*b.xx() + yy()*b.yx() + yz()*b.zx() + yt()*b.tx(),
    yx()*b.xy() + yy()*b.yy() + yz()*b.zy() + yt()*b.ty(),
    yx()*b.xz() + yy()*b.yz() + yz()*b.zz() + yt()*b.tz(),
    yx()*b.xt() + yy()*b.yt() + yz()*b.zt() + yt()*b.tt(),

    zx()*b.xx() + zy()*b.yx() + zz()*b.zx() + zt()*b.tx(),
    zx()*b.xy() + zy()*b.yy() + zz()*b.zy() + zt()*b.ty(),
    zx()*b.xz() + zy()*b.yz() + zz()*b.zz() + zt()*b.tz(),
    zx()*b.xt() + zy()*b.yt() + zz()*b.zt() + zt()*b.tt(),

    tx()*b.xx() + ty()*b.yx() + tz()*b.zx() + tt()*b.tx(),
    tx()*b.xy() + ty()*b.yy() + tz()*b.zy() + tt()*b.ty(),
    tx()*b.xz() + ty()*b.yz() + tz()*b.zz() + tt()*b.tz(),
    tx()*b.xt() + ty()*b.yt() + tz()*b.zt() + tt()*b.tt());
}
