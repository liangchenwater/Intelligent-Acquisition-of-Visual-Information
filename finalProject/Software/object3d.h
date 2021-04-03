#ifndef _OBJECT3D_H
#define _OBJECT3D_H
#include "ray.h"
#include "hit.h" 
#include "material.h"
 #include <GL/gl.h>
 #include <GL/glu.h>
 #include <GL/glut.h>
 #include<vector>
 using namespace std;
 
class Object3D{
	public:
		Object3D(){}
		Object3D(Material* _m):m(_m){}
 		virtual bool intersect (const Ray& r, Hit & h, float tmin)=0;
 		virtual bool intersectShadowRay(const Ray&, float)=0;
 		virtual void intersectTShadowRay(const Ray&, Vec3f&, vector<Hit>&, float)=0;
 		virtual void paint(void) = 0; 
 		virtual ~Object3D(){}
 	protected:
 		Material *m;	
 };
#endif
