#include "group.h"
#define eps 1e-6
Group::Group(int _num):num(_num){
	arr= new Object3D*[num];
}

Group::~Group(){
	delete[] arr;
	arr=NULL;
}

void Group::addObject(int index, Object3D* obj)
{
		arr[index]=obj;
}

bool Group::intersect(const Ray& r,Hit& h,float tmin)
{
	for(int i=0;i<num;i++) arr[i]->intersect(r,h,tmin);
	return true;
}

void Group::paint()
{
	for(int i=0;i<num;i++) arr[i]->paint();
}

bool Group::intersectShadowRay(const Ray& r, float tmin)
{
    for(int i=0;i<num;i++) 
        if(arr[i]->intersectShadowRay(r,tmin)) return true;
    return false;
}

void Group::intersectTShadowRay(const Ray& r, Vec3f& term,vector<Hit>& v ,float tmin)
{
    for(int i=0;i<num;i++) {
    	arr[i]->intersectTShadowRay(r,term,v,tmin);
    	if(term.Length()<=eps) break;
    }
}
