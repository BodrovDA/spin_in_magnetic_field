// -*- C++ -*-
// CLASSDOC OFF
// $Id: LorentzRotation.h 10030 2007-03-09 07:51:44Z katayama $
// ---------------------------------------------------------------------------
// CLASSDOC ON
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the definition of the HepLorentzRotation class for performing 
// Lorentz transformations (rotations and boosts) on objects of the
// HepLorentzVector class.
//
// .SS See Also
// ThreeVector.h, LorentzVector.h, Rotation.h
//
// .SS Author
// Leif Lonnblad. Modified by Evgueni Tcherniaev.

#ifndef BELLE_HEP_LORENTZROTATION_H
#define BELLE_HEP_LORENTZROTATION_H

#include "Rotation.h"
#include "LorentzVector.h"


class HepLorentzRotation {

public:

  inline HepLorentzRotation();
  // Default constructor. Gives a unit matrix.

  inline HepLorentzRotation(const HepRotation &);
  // Constructor for 3d rotations.

  inline HepLorentzRotation(const HepLorentzRotation &);
  // Copy constructor.

  inline HepLorentzRotation(double, double, double);
  inline HepLorentzRotation(const Hep3Vector &);
  // Constructors giving a Lorenz-boost.

  inline double xx() const;
  inline double xy() const;
  inline double xz() const;
  inline double xt() const;
  inline double yx() const;
  inline double yy() const;
  inline double yz() const;
  inline double yt() const;
  inline double zx() const;
  inline double zy() const;
  inline double zz() const;
  inline double zt() const;
  inline double tx() const;
  inline double ty() const;
  inline double tz() const;
  inline double tt() const;
  // Elements of the matrix.
  
  double operator () (int, int) const;
  // Returns (i,j) element of the matrix.

  inline HepLorentzRotation & operator = (const HepLorentzRotation &);
  inline HepLorentzRotation & operator = (const HepRotation &);
  // Assignment.

  inline bool operator == (const HepLorentzRotation &) const;
  inline bool operator != (const HepLorentzRotation &) const;
  // Comparisons.

  inline bool isIdentity() const;
  // Returns true if the Identity matrix.

  inline HepLorentzVector vectorMultiplication(const HepLorentzVector&) const;
  inline HepLorentzVector operator * (const HepLorentzVector &) const;
  // Multiplication with a Lorentz vector.

  HepLorentzRotation matrixMultiplication(const HepLorentzRotation &) const;
  inline HepLorentzRotation operator * (const HepLorentzRotation &) const;
  inline HepLorentzRotation & operator *= (const HepLorentzRotation &);
  inline HepLorentzRotation & transform(const HepLorentzRotation &);
  // Matrix multiplication.
  // Note: a *= b; <=> a = a * b; while a.transform(b); <=> a = b * a;

  inline HepLorentzRotation inverse() const;
  // Return the inverse.

  inline HepLorentzRotation & invert();
  // Inverts the LorentzRotation matrix.

  inline HepLorentzRotation & boost(double, double, double);
  inline HepLorentzRotation & boost(const Hep3Vector &);
  // Lorenz boost.

  inline HepLorentzRotation & rotateX(double);
  // Rotation around x-axis.

  inline HepLorentzRotation & rotateY(double);
  // Rotation around y-axis.

  inline HepLorentzRotation & rotateZ(double);
  // Rotation around z-axis.

  inline HepLorentzRotation & rotate(double, const Hep3Vector &);
  inline HepLorentzRotation & rotate(double, const Hep3Vector *);
  // Rotation around specified vector.

private:

  double mxx, mxy, mxz, mxt,
            myx, myy, myz, myt,
            mzx, mzy, mzz, mzt,
            mtx, mty, mtz, mtt;
  // The matrix elements.

  void setBoost(double, double, double);
  // Set elements according to a boost vector.

  inline HepLorentzRotation(double, double, double, double,
                            double, double, double, double,
                            double, double, double, double,
                            double, double, double, double);
  // Private constructor.
};


typedef HepLorentzRotation LRotation;

#include "LorentzRotation.icc"
#endif /* HEP_LORENTZROTATION_H */
