#ifndef _GROUP_H
#define _GROUP_H
#include "object3d.h"
 #include <GL/gl.h>
 #include <GL/glu.h>
 #include <GL/glut.h>
 #include <vector>
 using namespace std;
 
class Group:public Object3D{
		public:
			Group(int _num);
			void addObject(int index, Object3D* obj);
			virtual bool intersect(const Ray&, Hit&,float);
			virtual bool intersectShadowRay(const Ray&, float);
			virtual void intersectTShadowRay(const Ray&, Vec3f&,vector<Hit>&,float);
			virtual void paint(void); 
			~Group();
		private:
			int num;
			Object3D** arr;
};

#endif