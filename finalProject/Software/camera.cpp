#include "camera.h"
#include<iostream>
using namespace std;
#define EPS 1e-6
#define ZERO(x) ((x)<=EPS&&(x)>=-EPS)
Ray OrthographicCamera::generateRay(Vec2f point)
{
	Vec3f orig=center+size*(point.x()-0.5)*horizontal+size*(point.y()-0.5)*screen_up;
	return Ray(orig,direction);
}

float OrthographicCamera::getTMin() const
{
	return EPS;
}

OrthographicCamera::OrthographicCamera(Vec3f& _center, Vec3f& _dir, Vec3f& _up,float _size):size(_size),center(_center),direction(_dir),up(_up){
		direction.Normalize();
	Vec3f::Cross3(horizontal,direction,up);
	horizontal.Normalize();
	Vec3f::Cross3(screen_up,horizontal,direction);
}

void OrthographicCamera::glInit(int w, int h) 
{ 
	glMatrixMode(GL_PROJECTION); 
	glLoadIdentity(); 
	if (w > h) glOrtho(-size/2.0, size/2.0, -size*(float)h/(float)w/2.0, size*(float)h/(float)w/2.0, 0.5, 40.0); 
	else glOrtho(-size*(float)w/(float)h/2.0, size*(float)w/(float)h/2.0, -size/2.0, size/2.0, 0.5, 40.0); 
} 

void OrthographicCamera::glPlaceCamera(void) 
{ 
	gluLookAt(center.x(), center.y(), center.z(), center.x()+direction.x(), center.y()+direction.y(), center.z()+direction.z(), up.x(), up.y(), up.z()); 
}

void OrthographicCamera::dollyCamera(float dist) { center += direction*dist; } 

void OrthographicCamera::truckCamera(float dx, float dy) 
{ 
	center += horizontal*dx + screen_up*dy;
} 

void OrthographicCamera::rotateCamera(float rx, float ry) 
{ 
	float tiltAngle = acos(up.Dot3(direction)); 
	if (tiltAngle-ry > 3.13)  ry = tiltAngle - 3.13; 
	else if (tiltAngle-ry < 0.01) ry = tiltAngle - 0.01; 
	Matrix rotMat = Matrix::MakeAxisRotation(up, rx); 
	rotMat *= Matrix::MakeAxisRotation(horizontal, ry); 
	rotMat.Transform(center); 
	rotMat.TransformDirection(direction); 
	rotMat.TransformDirection(up);
	direction.Normalize();
	Vec3f::Cross3(horizontal,direction,up);
	horizontal.Normalize();
	Vec3f::Cross3(screen_up,horizontal,direction);
}

void OrthographicCamera::moveAroundCamera(Vec3f axis,Vec3f r,float ro)
{
	Matrix rotMat= Matrix::MakeAxisRotation(axis,ro);
	Matrix tlMat=Matrix::MakeTranslation(-1.0f*r);
	rotMat=rotMat*tlMat;
	rotMat.Transform(center);
	rotMat.TransformDirection(direction);
	rotMat.TransformDirection(up);
	rotMat.TransformDirection(horizontal);
	rotMat.TransformDirection(screen_up);
	tlMat=-1.0f*tlMat;
	tlMat.Transform(center);
	tlMat.TransformDirection(direction);
	tlMat.TransformDirection(up);
	tlMat.TransformDirection(horizontal);
	tlMat.TransformDirection(screen_up);}


Ray PerspectiveCamera::generateRay(Vec2f point)
{
	Vec3f fake_center=center+direction*(1.0f/(2*tan(angle/2)));
	Vec3f fake_orig=fake_center+(point.x()-0.5-cx)*horizontal+(point.y()-0.5-cy)*screen_up;
	Vec3f dir=fake_orig-center;
	dir.Normalize();
	return Ray(center,dir);
}

float PerspectiveCamera::getTMin() const 
{
	return EPS;
}

PerspectiveCamera::PerspectiveCamera(Vec3f& _center, Vec3f& _dir, Vec3f& _up,float _anglex,float _cx,float _cy):angle(_anglex),center(_center),direction(_dir),up(_up),cx(_cx),cy(_cy)
{
	direction.Normalize();
	Vec3f::Cross3(horizontal,direction,up);
	horizontal.Normalize();
	Vec3f::Cross3(screen_up,horizontal,direction);
}

void PerspectiveCamera::glInit(int w, int h)
{ 
	glMatrixMode(GL_PROJECTION); 
	glLoadIdentity(); 
	gluPerspective(angle*180.0/3.14159, (float)w/float(h), 0.5, 40.0); 
} 

void PerspectiveCamera::glPlaceCamera(void)
 {
 	 gluLookAt(center.x(), center.y(), center.z(), center.x()+direction.x(), center.y()+direction.y(), center.z()+direction.z(), up.x(), up.y(), up.z()); 
} 

void PerspectiveCamera::dollyCamera(float dist) { center += direction*dist;}

void PerspectiveCamera::truckCamera(float dx, float dy) { center += horizontal*dx + screen_up*dy;}

void PerspectiveCamera::rotateCamera(float rx, float ry)
 { 
	float tiltAngle = acos(up.Dot3(direction)); 
	if (tiltAngle-ry > 3.13) ry = tiltAngle - 3.13; 
	else if (tiltAngle-ry < 0.01) ry = tiltAngle - 0.01; 
	Matrix rotMat = Matrix::MakeAxisRotation(up, rx); 
	rotMat *= Matrix::MakeAxisRotation(horizontal, ry);
	rotMat.Transform(center); 
	rotMat.TransformDirection(direction); 
	rotMat.TransformDirection(up);
	direction.Normalize();
	Vec3f::Cross3(horizontal,direction,up);
	horizontal.Normalize();
	Vec3f::Cross3(screen_up,horizontal,direction);
}

void PerspectiveCamera::moveAroundCamera(Vec3f axis,Vec3f r,float ro)
{
float tiltAngle = acos(up.Dot3(direction)); 
	if (tiltAngle-ro > 3.13) ro = tiltAngle - 3.13; 
	else if (tiltAngle-ro < 0.01) ro = tiltAngle - 0.01; 
	Matrix rotMat= Matrix::MakeAxisRotation(axis,ro);
	Matrix tlMat=Matrix::MakeTranslation(-1.0f*r);
	tlMat.Transform(center);
	rotMat.Transform(center);
	rotMat.TransformDirection(direction);
	rotMat.TransformDirection(up);
	rotMat.TransformDirection(horizontal);
	rotMat.TransformDirection(screen_up);
	tlMat=-1.0f*tlMat;
	tlMat.Transform(center);
}
