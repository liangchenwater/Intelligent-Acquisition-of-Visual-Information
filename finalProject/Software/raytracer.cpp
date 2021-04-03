#include "raytracer.h"
#define eps 1e-6
#include <vector>
using namespace std;
 
inline float gainTravelDist(const vector<Hit>& v,const Ray& r)
{
	if(v.size()<=1) return 1;
	Vec3f dir=r.getDirection();
	stack<float> temp;
	float d=v[v.size()-1].getT()-v[0].getT(),t_skip=0;
	for(int i=0;i<v.size();i++){
		if(temp.empty()&&i>0){ t_skip+=v[i].getT()-v[i-1].getT();  /*cout<<i<<endl;*/}
		Material* m=v[i].getMaterial();
		float idx=m->getIndexOfRefraction();
		if(!temp.empty()&&temp.top()==idx) temp.pop();
		else temp.push(idx);
	}
	d-=t_skip;
	if(d<=eps) return 1;
	return d;
}

Vec3f RayTracer::traceRay(Ray& ray, float tmin, int bounces, float weight,Hit& hit)
{	
 	bool shadeback=this->shadeback;
	Group* g=scene->getGroup();
	g->intersect(ray,hit,tmin);
	Vec3f nowPoint=hit.getIntersectionPoint(),dir_R=-1.0f*ray.getDirection(),n=hit.getNormal();
	dir_R.Normalize();n.Normalize();
	assert(!refractionStack.empty());
	float judge=n.Dot3(dir_R);
	int num_lights=scene->getNumLights();
	Material* bg=new PhongMaterial(scene->getBackgroundColor());
	
	Vec3f light_col,light_dir,res_col(0.0f,0.0f,0.0f);
	float dist_to_light=0;
	
	Vec3f reflectiveColor=hit.getMaterial()->getReflectiveColor();
	Vec3f transparentColor=hit.getMaterial()->getTransparentColor();
	Vec3f obj_col=hit.getMaterial()->getDiffuseColor();
	Vec3f ambient_col=scene->getAmbientLight();
	if(transparentColor.Length()>eps) shadeback=true;
		
	if(judge<-eps&&shadeback){
		 hit.set(hit.getT(),hit.getMaterial(),-1.0f*n,ray);
		 judge*=-1.0f;
		 n*=-1.0f;
	}
	
	
	if(judge<-eps)  return res_col;
	if(hit.getT()==inf) return obj_col;
	
	ambient_col.Set(obj_col.x()*ambient_col.x(),obj_col.y()*ambient_col.y(),obj_col.z()*ambient_col.z());
	res_col+=ambient_col;
	if(!shadows) for(int i=0;i<num_lights;i++){
		scene->getLight(i)->getIllumination(nowPoint,light_dir,light_col,dist_to_light);
		res_col+=hit.getMaterial()->Shade(ray,hit,light_dir,light_col);
		}
	else for(int i=0;i<num_lights;i++){
		scene->getLight(i)->getIllumination(nowPoint,light_dir,light_col,dist_to_light);
		light_dir.Normalize();
		Ray r(nowPoint,light_dir,1);
	   if(!transparent_shadow){
       	if(!g->intersectShadowRay(r,sqrt(epsilon))) res_col+=hit.getMaterial()->Shade(ray,hit,light_dir,light_col);
      	}
      else{
        Vec3f shadowTerm(1,1,1);
      	vector<Hit> v;
      	g->intersectTShadowRay(r,shadowTerm,v,sqrt(epsilon));
      	shadowTerm.Set((float)sqrt(shadowTerm.x()),(float)sqrt(shadowTerm.y()),(float)sqrt(shadowTerm.z()));
      	if(semi){
      	sort(v.begin(),v.end());
       	if(shadowTerm.Length()>eps) shadowTerm*=(semi_coe/(float)pow(gainTravelDist(v,r),0.15));
       	}
      	Vec3f render_col=hit.getMaterial()->Shade(ray,hit,light_dir,light_col);
      	render_col.Set(render_col.x()*shadowTerm.x(),render_col.y()*shadowTerm.y(),render_col.z()*shadowTerm.z());
      	res_col+=render_col;
      }
	}
	
	if(reflectiveColor.Length()>eps&&bounces>0&&weight>cutoff_weight){
		Vec3f out=2*n.Dot3(dir_R)*n-dir_R;
		out.Normalize();
		Ray r(nowPoint,out,2);
		Hit reflectHit(inf,bg,-1.0f*out);
		stack<float> tmp=refractionStack;
		Vec3f reflectiveTerm=traceRay(r,epsilon,bounces-1,weight*reflectiveColor.Length(),reflectHit);
		refractionStack=tmp;
		reflectiveTerm.Set(reflectiveTerm.x()*reflectiveColor.x(),reflectiveTerm.y()*reflectiveColor.y(),reflectiveTerm.z()*reflectiveColor.z());
		res_col+=reflectiveTerm;
	}
	
	if(transparentColor.Length()>eps&&bounces>0&&weight>cutoff_weight)
	{
		bool trans_flag=0;//in
		float index_obj=hit.getMaterial()->getIndexOfRefraction(),relativeIndex=refractionStack.top()/index_obj;
		if(index_obj==bg->getIndexOfRefraction()) index_obj*=1.0001;
		if(index_obj==refractionStack.top()) {
			trans_flag=1;//out
			refractionStack.pop();
			relativeIndex=index_obj/refractionStack.top();
			refractionStack.push(index_obj);
			}
		float delta=1-relativeIndex*relativeIndex*(1-judge*judge);
		if(delta>epsilon){
			Vec3f out=(relativeIndex*judge-(float)sqrt(delta))*n-relativeIndex*dir_R;
			out.Normalize();
			Ray r(nowPoint,out,3);
			Hit refractionHit(inf,bg,-1.0f*out);
			if(trans_flag==0) refractionStack.push(index_obj);
			else {
				assert(!refractionStack.empty());
				refractionStack.pop();
			}
			Vec3f refractiveTerm=traceRay(r,epsilon,bounces-1,weight*transparentColor.Length(),refractionHit);
			refractiveTerm.Set(refractiveTerm.x()*transparentColor.x(),refractiveTerm.y()*transparentColor.y(),refractiveTerm.z()*transparentColor.z());
			res_col+=refractiveTerm;
		}
	}
	if(bg) delete bg;
	bg=NULL;	
	return res_col;
}

void RayTracer::traceRay(Ray& ray,float tmin, int bounces,float weight,Hit& hit,int idx)
{
	bool shadeback=this->shadeback;
	
	Group* g=scene->getGroup();
	g->intersect(ray,hit,tmin);
	float tstop=hit.getT();
	Material* bg=new PhongMaterial(scene->getBackgroundColor());
	Vec3f reflectiveColor=hit.getMaterial()->getReflectiveColor();
	Vec3f transparentColor=hit.getMaterial()->getTransparentColor();

	if(hit.getT()==inf) return;
		
	Vec3f nowPoint=hit.getIntersectionPoint(),dir_R=-1.0f*ray.getDirection(),n=hit.getNormal();
	dir_R.Normalize();n.Normalize();
	float judge=n.Dot3(dir_R);
	
	if(transparentColor.Length()>eps) shadeback=true;
	if(judge<-eps&&shadeback){
		 hit.set(hit.getT(),hit.getMaterial(),-1.0f*n,ray);
		 judge*=-1.0f;
		 n*=-1.0f;
	}
	if(judge<-eps)  return;
	
	assert(idx>=0&&idx<=2);
	if(idx==0) RayTree::SetMainSegment(ray,0, tstop);
	else if(idx==1) RayTree::AddReflectedSegment(ray,0,tstop);
	else RayTree::AddTransmittedSegment(ray,0,tstop);

	int num_lights=scene->getNumLights();
	
	Vec3f light_col,light_dir;
	float dist_to_light=0;
	
	if(shadows)
	for(int i=0;i<num_lights;i++){
		scene->getLight(i)->getIllumination(nowPoint,light_dir,light_col,dist_to_light);
		light_dir.Normalize();
		Ray r(nowPoint,light_dir,1);
		Hit ShadowHit(inf,bg,light_dir);
		if(!transparent_shadow){
       	if(g->intersect(r,ShadowHit,sqrt(epsilon))) RayTree::AddShadowSegment(r,0,ShadowHit.getT()); 
       	}
       	else{
       		 Vec3f shadowTerm(1,1,1);
      		vector<Hit> v;
      		g->intersectTShadowRay(r,shadowTerm,v,sqrt(epsilon));
       		if(!v.empty()) RayTree::AddShadowSegment(r,0,v[v.size()-1].getT());
       	}
	}
	
	if(reflectiveColor.Length()>eps&&bounces>0&&weight>cutoff_weight){
		Vec3f out=2*n.Dot3(dir_R)*n-dir_R;
		Ray r(nowPoint,out,2);
		Hit reflectHit(inf,bg,-1.0f*out);
		stack<float> tmp=refractionStack;
		traceRay(r,epsilon,bounces-1,weight*reflectiveColor.Length(),reflectHit,1);
		refractionStack=tmp;
	}
	
	if(transparentColor.Length()>eps&&weight>cutoff_weight)
	{
		bool trans_flag=0;//in
		float index_obj=hit.getMaterial()->getIndexOfRefraction(),relativeIndex=refractionStack.top()/index_obj;
		if(index_obj==bg->getIndexOfRefraction()) index_obj*=1.0001;
		if(index_obj==refractionStack.top()) {
			trans_flag=1;//out
			refractionStack.pop();
			relativeIndex=index_obj/refractionStack.top();
			refractionStack.push(index_obj);
			}
		float delta=1-relativeIndex*relativeIndex*(1-judge*judge);
		if(delta>0){
			Vec3f out=(relativeIndex*judge-(float)sqrt(delta))*n-relativeIndex*dir_R;
			Ray r(nowPoint,out,3);
			Hit refractionHit(inf,bg,-1.0f*out);
			if(trans_flag==0) refractionStack.push(index_obj);
			else {
				assert(!refractionStack.empty());
				refractionStack.pop();
			}
			cout<<"Now the relative index is:"<<relativeIndex<<endl;
			traceRay(r,epsilon,bounces-1,weight*transparentColor.Length(),refractionHit,2);
			}
		
	}
	if(bg) delete bg;
	bg=NULL;	
}

