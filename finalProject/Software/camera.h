#ifndef _CAMERA_H
#define _CAMERA_H

#include "vectors.h"
#include "ray.h"
#include "matrix.h"
 #include <GL/gl.h>
 #include <GL/glu.h>
 #include <GL/glut.h>

class Camera{
	public:
		Camera(){}
		virtual ~Camera(){}
		virtual Ray generateRay(Vec2f point)=0;
		virtual float getTMin() const=0;
		virtual void glInit(int w, int h) = 0; 
		virtual void glPlaceCamera(void) = 0; 
		virtual void dollyCamera(float dist) = 0; 
		virtual void truckCamera(float dx, float dy) = 0; 
		virtual void rotateCamera(float rx, float ry) = 0; 
		virtual void moveAroundCamera(Vec3f,Vec3f,float)=0;
		virtual Vec3f getDir() const=0;
		virtual Vec3f getPos() const=0;
};

class OrthographicCamera:public Camera{
	public:
		OrthographicCamera(Vec3f& _center, Vec3f& _dir, Vec3f& _up,float _size);
		Ray generateRay(Vec2f point);
		virtual float getTMin() const;
		virtual void glInit(int w, int h) ; 
		virtual void glPlaceCamera(void) ; 
		virtual void dollyCamera(float dist) ; 
		virtual void truckCamera(float dx, float dy); 
		virtual void rotateCamera(float rx, float ry); 
		virtual void moveAroundCamera(Vec3f,Vec3f,float);
		Vec3f getDir()const{return direction;}
		Vec3f getPos() const{return center;}
	private:
		Vec3f center;
		Vec3f direction;
		Vec3f up;
		Vec3f screen_up;
		Vec3f horizontal;
		float size;
};

class PerspectiveCamera:public Camera{
	public:
		PerspectiveCamera(Vec3f&,Vec3f&,Vec3f&,float,float,float);
		Ray generateRay(Vec2f point);
		virtual float getTMin() const;
		virtual void glInit(int w, int h) ; 
		virtual void glPlaceCamera(void); 
		virtual void dollyCamera(float dist); 
		virtual void truckCamera(float dx, float dy); 
		virtual void rotateCamera(float rx, float ry) ; 
		virtual void moveAroundCamera(Vec3f,Vec3f,float);
		Vec3f getDir()const{return direction;}
		Vec3f getPos() const{return center;}
	private:
		Vec3f center;
		Vec3f direction;
		Vec3f up;
		Vec3f screen_up;
		Vec3f horizontal;
		float angle;
		float cx;
		float cy;
};
#endif