#include "material.h"
#include "glCanvas.h"  
#ifdef SPECULAR_FIX
extern int SPECULAR_FIX_WHICH_PASS;
#endif

Vec3f PhongMaterial::Shade(const Ray &ray, const Hit &hit, const Vec3f &dirToLight, const Vec3f &lightColor) const
{
    Vec3f dir_L=dirToLight,dir_R=-1.0f*ray.getDirection();
	 Vec3f normal=hit.getNormal();
	  dir_L.Normalize();
	 normal.Normalize();
   dir_R.Normalize();
	 Vec3f h=dir_R+dir_L;
	 h.Normalize();
	 if(normal.Dot3(dir_L)<0) return Vec3f(0,0,0);
	 Vec3f res=normal.Dot3(dir_L)*diffuseColor;
	#if SPECULAR_FIX
       res+=(float)pow(normal.Dot3(h),exponent)*specularColor*normal.Dot3(dir_L);    
    #else
        res+=(float)pow(normal.Dot3(h),exponent)*specularColor;
       #endif 

	 res.Set(res.x()*lightColor.x(),res.y()*lightColor.y(),res.z()*lightColor.z());
	return res;
}

void PhongMaterial::glSetMaterial(void) const {

  GLfloat one[4] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat zero[4] = { 0.0, 0.0, 0.0, 0.0 };
  GLfloat specular[4] = {
    getSpecularColor().r(),
    getSpecularColor().g(),
    getSpecularColor().b(),
    1.0};
  GLfloat diffuse[4] = {
    getDiffuseColor().r(),
    getDiffuseColor().g(),
    getDiffuseColor().b(),
    1.0};
  // NOTE: GL uses the Blinn Torrance version of Phong...      
  float glexponent = exponent;
  if (glexponent < 0) glexponent = 0;
  if (glexponent > 128) glexponent = 128;
  
#if !SPECULAR_FIX 
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, diffuse);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, &glexponent);
  
#else

  // OPTIONAL: 3 pass rendering to fix the specular highlight 
  // artifact for small specular exponents (wide specular lobe)

  if (SPECULAR_FIX_WHICH_PASS == 0) {
    // First pass, draw only the specular highlights
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, zero);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, zero);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, &glexponent);
    
  } else if (SPECULAR_FIX_WHICH_PASS == 1) {
    // Second pass, compute normal dot light 
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, one);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, zero);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, zero);
  } else {
    // Third pass, add ambient & diffuse terms
    assert (SPECULAR_FIX_WHICH_PASS == 2);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, zero);
  }

#endif
}
