// -*- C++ -*-
// CLASSDOC OFF
// $Id: LorentzVector.h 10030 2007-03-09 07:51:44Z katayama $
// ---------------------------------------------------------------------------
// CLASSDOC ON
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// HepLorentzVector is a Lorentz vector consisting of Hep3Vector and
// HepDouble components. Lorentz transformations (rotations and boosts)
// of these vectors are perfomed by multiplying with objects of
// the HepLorenzRotation class.
//
// .SS See Also
// ThreeVector.h, Rotation.h, LorentzRotation.h
//
// .SS Authors
// Leif Lonnblad and Anders Nilsson. Modified by Evgueni Tcherniaev.
//

#ifndef BELLE_HEP_LORENTZVECTOR_H
#define BELLE_HEP_LORENTZVECTOR_H

#ifndef BELLE_CLHEP_THREEVECTOR_H
#include "ThreeVector.h"
#endif

#include <iostream>

class HepLorentzRotation;
class HepRotation;


class HepLorentzVector {

public:

  inline HepLorentzVector(void);
  explicit inline HepLorentzVector(double x);
  inline HepLorentzVector(double x, double y,
			  double z = 0.0, double t = 0.0);
  // Constructor giving the components x, y, z, t.

  inline HepLorentzVector(const Hep3Vector &, double);
  // Constructor giving a 3-Vector and a time component.

  inline HepLorentzVector(const HepLorentzVector &);
  // Copy constructor.

  inline ~HepLorentzVector();
  // The destructor.

  inline operator Hep3Vector () const;
  inline operator Hep3Vector & ();
  // Conversion (cast) to Hep3Vector.

  inline double x() const;
  inline double y() const;
  inline double z() const;
  inline double t() const;
  // Get position and time.

  inline void setX(double);
  inline void setY(double);
  inline void setZ(double);
  inline void setT(double);
  // Set position and time.

  inline double px() const;
  inline double py() const;
  inline double pz() const;
  inline double e() const;
  // Get momentum and energy.

  inline void setPx(double);
  inline void setPy(double);
  inline void setPz(double);
  inline void setE(double);
  // Set momentum and energy.

  inline Hep3Vector vect() const;
  // Get spatial component. 

  inline void setVect(const Hep3Vector &);
  // Set spatial component. 

  inline double theta() const;
  inline double cosTheta() const;
  inline double phi() const;
  inline double rho() const;
  // Get spatial vector components in spherical coordinate system.

  inline void setTheta(double);
  inline void setPhi(double);
  inline void setRho(double);
  // Set spatial vector components in spherical coordinate system.

  inline double operator () (int) const;
  // Get components by index.

  inline HepLorentzVector & operator = (const HepLorentzVector &);
  // Assignment. 

  inline HepLorentzVector   operator +  (const HepLorentzVector &) const;
  inline HepLorentzVector & operator += (const HepLorentzVector &);
  // Additions.

  inline HepLorentzVector   operator -  (const HepLorentzVector &) const;
  inline HepLorentzVector & operator -= (const HepLorentzVector &);
  // Subtractions.

  inline HepLorentzVector operator - () const;
  // Unary minus.

  inline bool operator == (const HepLorentzVector &) const;
  inline bool operator != (const HepLorentzVector &) const;
  // Comparisons.

  inline double perp2() const;
  // Transverse component of the spatial vector squared.

  inline double perp() const;
  // Transverse component of the spatial vector (R in cylindrical system).

  inline double perp2(const Hep3Vector &) const;
  // Transverse component of the spatial vector w.r.t. given axis squared.

  inline double perp(const Hep3Vector &) const;
  // Transverse component of the spatial vector w.r.t. given axis.

  inline double angle(const Hep3Vector &) const;
  // Angle wrt. another vector.

  inline double mag2() const;
  inline double m2() const;
  // Invariant mass squared.

  inline double mag() const;
  inline double m() const;
  // Invariant mass. If mag2() is negative then -sqrt(-mag2()) is returned.

  inline double dot(const HepLorentzVector &) const;
  inline double operator * (const HepLorentzVector &) const;
  // Scalar product.

  inline double plus() const;
  inline double minus() const;
  // Returns the positive/negative light-cone component t +/- z.


  inline Hep3Vector boostVector() const ;
  // Returns the spatial components divided by the time component.

  void boost(double, double, double);
  inline void boost(const Hep3Vector &);
  // Lorentz boost.

  inline void rotateX(double);
  // Rotate the spatial component around the x-axis.

  inline void rotateY(double);
  // Rotate the spatial component around the y-axis.

  inline void rotateZ(double);
  // Rotate the spatial component around the z-axis.

  inline void rotateUz(Hep3Vector &);
  // Rotates the reference frame from Uz to newUz (unit vector).

  inline void rotate(double, const Hep3Vector &);
  // Rotate the spatial component around specified axis.

  HepLorentzVector & operator *= (const HepRotation &);
  HepLorentzVector & transform(const HepRotation &);
  // Transformation with HepRotation.

  HepLorentzVector & operator *= (const HepLorentzRotation &);
  HepLorentzVector & transform(const HepLorentzRotation &);
  // Transformation with HepLorenzRotation.

private:

  Hep3Vector pp;
  double  ee;

};

std::ostream & operator << (std::ostream &, const HepLorentzVector &);
// Output to a stream.

typedef HepLorentzVector VectorL;
typedef HepLorentzVector Vector4;
typedef HepLorentzVector DVectorL;
typedef HepLorentzVector DVector4;
typedef HepLorentzVector FVectorL;
typedef HepLorentzVector FVector4;

typedef HepLorentzVector HepLorentzVectorD;
typedef HepLorentzVector HepLorentzVectorF;

#include "LorentzVector.icc"



#endif /* HEP_LORENTZVECTOR_H */
