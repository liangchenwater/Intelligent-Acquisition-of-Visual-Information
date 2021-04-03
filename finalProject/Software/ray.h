#ifndef _RAY_H
#define _RAY_H

#include <iostream>
using namespace std;

#include "vectors.h"

// ====================================================================
// ====================================================================

// Ray class mostly copied from Peter Shirley and Keith Morley

class Ray {

public:

  // CONSTRUCTOR & DESTRUCTOR
  Ray () {}
  Ray (const Vec3f &orig, const Vec3f &dir) {
    origin = orig; 
    direction = dir; 
    }
  Ray (const Ray& r) {*this=r;}
  Ray (const Vec3f &orig, const Vec3f &dir,int tp) {
    origin = orig; 
    direction = dir; 
    type=tp;
 }

  // ACCESSORS
  const Vec3f& getOrigin() const { return origin; }
  const Vec3f& getDirection() const { return direction; }
  Vec3f pointAtParameter(float t) const {
    return origin+direction*t; }
  int getType() const{ return type;}

private:

  // REPRESENTATION
  Vec3f origin;
  Vec3f direction;
  int type=0;
};

inline ostream &operator<<(ostream &os, const Ray &r) {
  os << "Ray <o:" <<r.getOrigin()<<", d:"<<r.getDirection()<<">";
  return os;
}

// ====================================================================
// ====================================================================

#endif