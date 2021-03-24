// -*- C++ -*-
// CLASSDOC OFF
// $Id: ThreeVector.h 10030 2007-03-09 07:51:44Z katayama $
// ---------------------------------------------------------------------------
// CLASSDOC ON
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// Hep3Vector is a general 3-vector class defining vectors in three
// dimension using double components. Rotations of these vectors are
// performed by multiplying with an object of the HepRotation class.
//
// .SS See Also
// LorentzVector.h, Rotation.h, LorentzRotation.h
//
// .SS Authors
// Leif Lonnblad and Anders Nilsson.
//

#ifndef BELLE_CLHEP_THREEVECTOR_H
#define BELLE_CLHEP_THREEVECTOR_H



#include <iosfwd>
#include <cmath>


class HepRotation;

class Hep3Vector {

public:

  inline Hep3Vector(void);
  inline explicit Hep3Vector(double x);
  inline Hep3Vector(double x, double y, double z = 0.0);
  // The constructor.

  inline Hep3Vector(const Hep3Vector &);
  // The copy constructor.

  inline ~Hep3Vector();
  // The destructor.

  inline double x() const;
  inline double y() const;
  inline double z() const;
  // The components in cartesian coordinate system.

  double operator () (int) const;
  // Get components by index (Geant4).

  inline void setX(double);
  inline void setY(double);
  inline void setZ(double);
  // Set the components in cartesian coordinate system.

  inline double phi() const;
  // The azimuth angle.

  inline double theta() const;
  // The polar angle.

  inline double cosTheta() const;
  // Cosine of the polar angle.

  inline double mag2() const;
  // The magnitude squared (rho^2 in spherical coordinate system).

  inline double mag() const;
  // The magnitude (rho in spherical coordinate system).

  inline void setPhi(double);
  // Set phi keeping mag and theta constant (BaBar).

  inline void setTheta(double);
  // Set theta keeping mag and phi constant (BaBar).

  inline void setMag(double);
  // Set magnitude keeping theta and phi constant (BaBar).

  inline double perp2() const;
  // The transverse component squared (R^2 in cylindrical coordinate system).

  inline double perp() const;
  // The transverse component (R in cylindrical coordinate system).

  inline double perp2(const Hep3Vector &) const;
  // The transverse component w.r.t. given axis squared.

  inline double perp(const Hep3Vector &) const;
  // The transverse component w.r.t. given axis.

  inline Hep3Vector & operator = (const Hep3Vector &);
  // Assignment.

  inline bool operator == (const Hep3Vector &) const;
  inline bool operator != (const Hep3Vector &) const;
  // Comparisons (Geant4). 

  inline Hep3Vector & operator += (const Hep3Vector &);
  // Addition.

  inline Hep3Vector & operator -= (const Hep3Vector &);
  // Subtraction.

  inline Hep3Vector operator - () const;
  // Unary minus.

  inline Hep3Vector & operator *= (double);
  // Scaling with real numbers.

  inline Hep3Vector unit() const;
  // Unit vector parallel to this.

  inline Hep3Vector orthogonal() const;
  // Vector orthogonal to this (Geant4).

  inline double dot(const Hep3Vector &) const;
  // Scalar product.

  inline Hep3Vector cross(const Hep3Vector &) const;
  // Cross product.

  inline double angle(const Hep3Vector &) const;
  // The angle w.r.t. another 3-vector.

  void rotateX(double);
  // Rotates the Hep3Vector around the x-axis.

  void rotateY(double);
  // Rotates the Hep3Vector around the y-axis.

  void rotateZ(double);
  // Rotates the Hep3Vector around the z-axis.

  void rotateUz(Hep3Vector&);
  // Rotates reference frame from Uz to newUz (unit vector) (Geant4).

  void rotate(double, const Hep3Vector &);
  // Rotates around the axis specified by another Hep3Vector.

  Hep3Vector & operator *= (const HepRotation &);
  Hep3Vector & transform(const HepRotation &);
  // Transformation with a Rotation matrix.

private:

  double dx, dy, dz;
  // The components.
};


std::ostream & operator << (std::ostream &, const Hep3Vector &);
// Output to a stream.


typedef Hep3Vector Vector3;
typedef Hep3Vector DVector3;
typedef Hep3Vector FVector3;

typedef Hep3Vector HepThreeVectorD;
typedef Hep3Vector HepThreeVectorF;

#include "ThreeVector.icc"




#endif /* CLHEP_THREEVECTOR_H */
