#ifndef _MATERIAL_H_
#define _MATERIAL_H_

#include "vectors.h"
#include "ray.h"
#include "hit.h"
#include <cmath>
 #include <GL/gl.h>
 #include <GL/glu.h>
 #include <GL/glut.h>
class Material {

public:

  // CONSTRUCTORS & DESTRUCTOR
  Material(const Vec3f &d_color):diffuseColor(d_color){}
  virtual ~Material(){};
  virtual Vec3f Shade(const Ray &ray, const Hit &hit, const Vec3f &dirToLight, const Vec3f &lightColor) const=0;
  // ACCESSORS
  Vec3f getDiffuseColor() const { return diffuseColor; }
  virtual void glSetMaterial(void) const=0;
  virtual Vec3f getSpecularColor() const=0;
  virtual Vec3f getReflectiveColor() const=0;
  virtual Vec3f getTransparentColor() const=0;
  virtual float getExponent()=0;
  virtual float getIndexOfRefraction()=0;
protected:
  // REPRESENTATION
  Vec3f diffuseColor;
};

class PhongMaterial:public Material{
    public:
        PhongMaterial(const Vec3f & _diffuseColor, const Vec3f & _specularColor, float _exponent):Material(_diffuseColor),specularColor(_specularColor),exponent(_exponent){}
         PhongMaterial(const Vec3f &d_color):Material(d_color){}
         PhongMaterial(const Vec3f& _diffuseColor, const Vec3f& _specularColor, float _exponent, const Vec3f& _reflectiveColor, const Vec3f& _transparentColor,  float _indexOfRefraction):Material(_diffuseColor),specularColor(_specularColor),exponent(_exponent),reflectiveColor(_reflectiveColor),transparentColor(_transparentColor),indexOfRefraction(_indexOfRefraction){}
       virtual Vec3f Shade(const Ray &ray, const Hit &hit, const Vec3f &dirToLight, const Vec3f &lightColor) const;
       virtual void glSetMaterial(void) const;
        Vec3f getSpecularColor()const{return specularColor;}
        Vec3f getReflectiveColor() const{return reflectiveColor;}
         Vec3f getTransparentColor() const{return transparentColor;}
         float getExponent (){return exponent;}
        float getIndexOfRefraction(){return indexOfRefraction;}
    private:
        Vec3f specularColor;
        float exponent;
        Vec3f reflectiveColor;
        Vec3f transparentColor;
        float indexOfRefraction=1;
};
#endif