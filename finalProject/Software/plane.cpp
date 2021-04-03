#include "plane.h"
#define eps 1e-6
#define BIG 9*1e5

bool Plane::intersect(const Ray& r,Hit& h,float tmin)
{
		float judge=r.getDirection().Dot3(normal);
		if(judge<eps&&judge>-eps) return 0;
		float t=(d-r.getOrigin().Dot3(normal))/judge;
		if(t<tmin) return 0;
		if(t<h.getT())  h.set(t,m,normal,r);
		return 1;
}

void Plane::paint()
{
	Vec3f v(1,0,0),b1(0,0,0),b2(0,0,0),orig(d/normal.Length()/normal.Length()*normal);
	Vec3f::Cross3(b1,normal,v);
	if(b1.Length()<=eps&&b1.Length()>=-eps) b1.Set(0,1,0);
	Vec3f::Cross3(b2,normal,b1);
	b1.Normalize();
	b2.Normalize();
	Vec3f ver1(orig+b1*BIG),ver3(orig-b1*BIG),ver2(orig+b2*BIG),ver4(orig-b2*BIG);
	m->glSetMaterial();
	glBegin(GL_TRIANGLES);
   glNormal3f(normal.x(), normal.y(), normal.z());
   glVertex3f(ver1.x(), ver1.y(),ver1.z());
   glVertex3f(ver2.x(), ver2.y(),ver2.z());
   glVertex3f(ver3.x(), ver3.y(),ver3.z());
    glVertex3f(ver3.x(), ver3.y(),ver3.z());
   glVertex3f(ver4.x(), ver4.y(),ver4.z());
    glVertex3f(ver1.x(), ver1.y(),ver1.z());
   glEnd();
}

bool Plane::intersectShadowRay(const Ray& r, float tmin)
{
	float judge=r.getDirection().Dot3(normal);
	if(judge<=eps&&judge>=-eps) return 0;
	float t=(d-r.getOrigin().Dot3(normal))/judge;
	if(t<tmin) return 0;
	return 1;
}

void Plane::intersectTShadowRay(const Ray& r, Vec3f& term,vector<Hit>& v,float tmin)
{
	float judge=r.getDirection().Dot3(normal);
	if(judge<=eps&&judge>=-eps) return;
	float t=(d-r.getOrigin().Dot3(normal))/judge;
	if(t<tmin) return;
	Vec3f color=m->getTransparentColor();
	term.Set(term.x()*color.x(),term.y()*color.y(),term.z()*color.z());
	Hit h;
	h.set(t,m,normal,r);
	v.push_back(h);
}