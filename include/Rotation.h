// -*- C++ -*-
// CLASSDOC OFF
// $Id: Rotation.h 10030 2007-03-09 07:51:44Z katayama $
// ---------------------------------------------------------------------------
// CLASSDOC ON
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the definition of the HepRotation class for performing rotations
// on objects of the Hep3Vector class.
//
// .SS See Also
// ThreeVector.h, LorentzVector.h, LorentzRotation.h
//
// .SS Author
// Leif Lonnblad

#ifndef BELLE_HEP_ROTATION_H
#define BELLE_HEP_ROTATION_H

#ifndef BELLE_CLHEP_THREEVECTOR_H
#include "ThreeVector.h"
#endif




class HepRotation {

public:

  inline HepRotation();
  // Default constructor. Gives a unit matrix.

  inline HepRotation(const HepRotation &);
  // Copy constructor.

  inline double xx() const;
  inline double xy() const;
  inline double xz() const;
  inline double yx() const;
  inline double yy() const;
  inline double yz() const;
  inline double zx() const;
  inline double zy() const;
  inline double zz() const;
  // Elements of the rotation matrix (Geant4).

  double operator () (int, int) const;
  // Returns (i,j) element of the rotation matrix.

  inline HepRotation & operator = (const HepRotation &);
  // Assignment.

  inline bool operator == (const HepRotation &) const;
  inline bool operator != (const HepRotation &) const;
  // Comparisons (Geant4).

  inline bool isIdentity() const;
  // Returns true if the identity matrix (Geant4).

  inline Hep3Vector operator * (const Hep3Vector &) const;
  // Multiplication with a Hep3Vector.

  HepRotation operator * (const HepRotation &) const;
  inline HepRotation & operator *= (const HepRotation &);
  inline HepRotation & transform(const HepRotation &);
  // Matrix multiplication.
  // Note a *= b; <=> a = a * b; while a.transform(b); <=> a = b * a;

  inline HepRotation inverse() const;
  // Returns the inverse.

  inline HepRotation & invert();
  // Inverts the Rotation matrix.

  HepRotation & rotateX(double);
  // Rotation around the x-axis.

  HepRotation & rotateY(double);
  // Rotation around the y-axis.

  HepRotation & rotateZ(double);
  // Rotation around the z-axis.

  HepRotation & rotate(double, const Hep3Vector &);
  inline HepRotation & rotate(double, const Hep3Vector *);
  // Rotation around a specified vector.

  HepRotation & rotateAxes(const Hep3Vector & newX,
                           const Hep3Vector & newY,
                           const Hep3Vector & newZ);
  // Rotation of local axes (Geant4).

  double phiX() const;
  double phiY() const;
  double phiZ() const;
  double thetaX() const;
  double thetaY() const;
  double thetaZ() const;
  // Return angles (RADS) made by rotated axes against original axes (Geant4).

  void getAngleAxis(double &, Hep3Vector &) const;
  // Returns the rotation angle and rotation axis (Geant4).

protected:

  inline HepRotation(double, double, double, double, double,
		     double, double, double, double);
  // Protected constructor.

  double rxx, rxy, rxz, ryx, ryy, ryz, rzx, rzy, rzz;
  // The matrix elements.
};


typedef HepRotation Rotation;


#include "Rotation.icc"



#endif /* HEP_ROTATION_H */

