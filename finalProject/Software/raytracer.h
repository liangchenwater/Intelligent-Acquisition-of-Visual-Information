#include "image.h"
#include "scene_parser.h"
#include "hit.h"
#include "glCanvas.h"
#include "matrix.h"
#include "camera.h" 
#include "light.h"
#include "material.h"
#include "object3d.h"
#include "group.h" 
#include "sphere.h"
#include "plane.h"
#include "triangle.h"
#include "transform.h"
#include "rayTree.h"
#include<cmath>
#include<stack>
using namespace std;

class RayTracer{
	public:
		 RayTracer(SceneParser *s, int max_bounces, float cutoff_weight, bool shadows,bool _shadeback,bool transparent_shadow, bool _semi,float _eps,float semi_coe):scene(s),max_bounces(max_bounces),cutoff_weight(cutoff_weight),shadows(shadows),epsilon(_eps),shadeback(_shadeback),transparent_shadow(transparent_shadow),semi(_semi),semi_coe(semi_coe){ refractionStack.push(1);}
		 Vec3f traceRay(Ray& ray, float tmin, int bounces, float weight,Hit& hit);
		 void traceRay(Ray&,float,int,float,Hit&,int);
		 void clearStack(float now_index){while(!refractionStack.empty()) refractionStack.pop();refractionStack.push(now_index);}
		 ~RayTracer(){if(scene) delete scene; scene=NULL;}
	private:
		SceneParser *scene;
		int max_bounces;
		float cutoff_weight;
		bool shadows;
		float epsilon;
		bool shadeback;
		bool transparent_shadow;
		bool semi;
		float semi_coe;
		stack<float> refractionStack;
};