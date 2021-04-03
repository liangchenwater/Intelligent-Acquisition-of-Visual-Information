#include "vectors.h"
#include "matrix.h"
#include "image.h"
#include <map>
#include "hit.h"
#include <cmath>
#include <fstream>
#define inf  (0xff<<23)
#include "scene_parser.h"
#include "camera.h" 
#include "light.h"
#include "material.h"
#include "object3d.h"
#include "group.h" 
#include "sphere.h"
#include "plane.h"
#include "triangle.h"
#include "transform.h"
#include "raytracer.h"
#include <cstring>
#include <cstdlib>
#include<opencv2/opencv.hpp>

using namespace std;

//from (theta_light,phi_light,theta_camera,phi_camera) to (R,G,B)
map<Vec4f,Vec3f> BRDF;

/*input:
0. camera related: hx,hy
1. initial pic: extrinsic matrix, camera matrix,
2. another pics to fix the track: extrinsic matrix
3. light's position: (xL,yL) in two pics, extrinsic matrices of these two pics
4. object's position: (x1,y1),(x2,y2),(x3,y3) in three pics, extrinsic matrices of these two pics
*/
string input_file;
class InputGenerator{
	public:
		static InputGenerator* readFile()
		{
			ifstream in(input_file);
			int n,width,height;
			float exposure_time;
			float _focus_length,_cx,_cy,rx,ry,rz,tx,ty,tz,colx,coly,colz;
			in>>n>>width>>height;
			in>>exposure_time;
			in>>_focus_length>>_cx>>_cy;
			in>>rx>>ry>>rz>>tx>>ty>>tz;
			in>>colx>>coly>>colz;
			_cx/=width;
			_cy/=width;
			_cx-=0.5;
			_cy-=0.5;
			Vec3f _T(tx,ty,tz);
			Vec3f Rot(rx,ry,rz);
			Vec3f _light_col(colx,coly,colz);
			float alpha=Rot.Length();
			Rot.Normalize();
			Matrix _rotate=Matrix::MakeAxisRotation(Rot,alpha);
			InputGenerator* IG=new InputGenerator(_focus_length,_cx,_cy,_rotate,_T,_light_col,n,width,height,exposure_time);
			float tx1,ty1,tz1,tx2,ty2,tz2;
			in>>tx1>>ty1>>tz1>>tx2>>ty2>>tz2;
			Vec3f ex1(-tx1,-ty1,-tz1),ex2(-tx2,-ty2,-tz2);
			int dev=3;
			IG->cal_axisAndTrackC(ex1,ex2,dev);
			int light_pix1_x,light_pix1_y,light_pix2_x,light_pix2_y;
			float cx1,cy1,cx2,cy2;
			in>>light_pix1_x>>light_pix1_y>>light_pix2_x>>light_pix2_y;
			in>>cx1>>cy1>>cx2>>cy2;
			//in>>tx1>>ty1>>tz1>>tx2>>ty2>>tz2;
			int dev1,dev2;
			in>>dev1>>dev2;
			float alpha1=1.0f*dev1*exposure_time*IG->getV()/IG->getTrackR(),alpha2=1.0f*dev2*exposure_time*IG->getV()/IG->getTrackR();
			IG->cal_light_pos(light_pix1_x,light_pix1_y,light_pix2_x,light_pix2_y,cx1,cy1,cx2,cy2,alpha1,alpha2);
			int dot1_x,dot1_y,dot2_x,dot2_y,dot3_x,dot3_y,dot4_x,dot4_y,dot1_x1,dot1_y1,dot2_x1,dot2_y1,dot3_x1,dot3_y1,dot4_x1,dot4_y1;
			in>>dot1_x>>dot1_y>>dot2_x>>dot2_y>>dot3_x>>dot3_y>>dot4_x>>dot4_y>>dot1_x1>>dot1_y1>>dot2_x1>>dot2_y1>>dot3_x1>>dot3_y1>>dot4_x1>>dot4_y1;
			in>>cx1>>cy1>>cx2>>cy2;
			//in>>tx1>>ty1>>tz1>>tx2>>ty2>>tz2;
			in>>dev1>>dev2;
			alpha1=1.0f*dev1*exposure_time*IG->getV()/IG->getTrackR(),alpha2=1.0f*dev2*exposure_time*IG->getV()/IG->getTrackR();
			IG->cal_obj(dot1_x,dot1_y,dot2_x,dot2_y,dot3_x,dot3_y,dot4_x,dot4_y,dot1_x1,dot1_y1,dot2_x1,dot2_y1,dot3_x1,dot3_y1,dot4_x1,dot4_y1,cx1,cy1,cx2,cy2,alpha1,alpha2);
			cout<<dot1_x<<endl<<dot1_y<<endl;
			return IG;
		}
		float getTrackR(){
			return trackR;
		}
		float getV()
		{
			return velocity;
		}
		InputGenerator(float _focus_length, float _cx, float _cy, Matrix& _rotate,Vec3f& _T,Vec3f& _light_col,int _n,int _width, int _height,float _exposure_time):focus_length(_focus_length),cx(_cx),cy(_cy),rotate(_rotate),T(_T),light_col(_light_col),pic_num(_n),width(_width),height(_height),exposure_time(_exposure_time)
		{
				fovx=2*atan(width/2/focus_length);
				fovy=2*atan(width/2/focus_length);
				camera_dir.Set(-rotate.Get(2,0),-rotate.Get(2,1),-rotate.Get(2,2));
				camera_up.Set(-rotate.Get(1,0),-rotate.Get(1,1),-rotate.Get(1,2));
				camera_pos.Set(-T.x(),-T.y(),-T.z());
				camera_dir.Normalize();
				camera_up.Normalize();
				Vec3f::Cross3(horizontal,camera_dir,camera_up);			
		}
		Vec3f Rotate(Vec3f base,float r)
		{
			Matrix rotMat=Matrix::MakeAxisRotation(axis,r);
			Matrix tlMat=Matrix::MakeTranslation(-1.0f*trackC);
			rotMat=rotMat*tlMat;
			rotMat.Transform(base);
			tlMat=-1.0f*tlMat;
			tlMat.Transform(base);
			return base;
		}
		Vec3f RotateDir(Vec3f base,float r)
		{
			Matrix rotMat=Matrix::MakeAxisRotation(axis,r);
			rotMat.TransformDirection(base);
			return base;
		}
		void cal_axisAndTrackC(Vec3f ex1,Vec3f ex2,int dev1)
		{
			Vec3f pos1=ex1,pos2=ex2;
			Vec3f a=pos1-camera_pos,b=pos2-camera_pos;
			Vec3f::Cross3(axis,a,b);
			axis.Normalize();
			Vec3f d1,d2;
			Vec3f::Cross3(d1,a,axis);
			Vec3f::Cross3(d2,b,axis);
			d1.Normalize();
			d2.Normalize();
			Vec3f o1=(1.0f/2)*(pos1+camera_pos),o2=(1.0f/2)*(pos2+camera_pos);
			Vec3f o=o1-o2;
			float t=(o.x()*d2.y()-o.y()*d2.x())/(d2.x()*d1.y()-d1.x()*d2.y()),tu=(o.z()*d2.y()-o.y()*d2.z())/(d2.z()*d1.y()-d1.z()*d2.y());
			trackC=o1+t*d1;
			Vec3f trackRVec=trackC-camera_pos;
			trackR=trackRVec.Length();
			Vec3f dist=camera_pos-ex1;
			velocity=dist.Length()/exposure_time/dev1;
			cout<<trackR<<endl<<trackC;
		}
		Vec3f cal_dot_pos(int x1,int y1,int x2,int y2,float cx1,float cy1,float cx2,float cy2,float alpha1,float alpha2)
		{
			float pix_x1=-1.0f*(x1+width/2.0f)/width,pix_y1=1.0f*(y1+width/2.0f)/width,pix_x2=-1.0f*(x2+width/2.0f)/width,pix_y2=1.0f*(y2+width/2.0f)/width,ratio=tan(fovx/2)/tan(fovy/2);
			Vec3f orig1=Rotate(camera_pos,alpha1),orig3=Rotate(camera_pos,alpha1),orig2=Rotate(camera_pos,alpha2);
			Vec3f camera_dir1=RotateDir(camera_dir,alpha1),camera_dir3=RotateDir(camera_dir,alpha1),camera_dir2=RotateDir(camera_dir,alpha2);
			Vec3f camera_up1=RotateDir(camera_up,alpha1),camera_up3=RotateDir(camera_up,alpha1),camera_up2=RotateDir(camera_up,alpha2);
			Vec3f horizontal1=RotateDir(horizontal,alpha1),horizontal3=RotateDir(horizontal,alpha1),horizontal2=RotateDir(horizontal,alpha2);
			camera_dir1.Normalize();
			camera_dir2.Normalize();
			camera_up1.Normalize();
			horizontal1.Normalize();
			camera_up2.Normalize();
			horizontal2.Normalize();
			Vec3f fc1=orig1+camera_dir1*(1.0f/(2*tan(fovx/2))),fc2=orig2+camera_dir2*(1.0f/(2*tan(fovx/2)));
			Vec3f c1=fc1+(pix_x1-0.5-cx1)*horizontal1+(pix_y1-0.5*ratio-cy1*ratio)*camera_up1,c2=fc2+(pix_x2-0.5-cx2)*horizontal2+(pix_y2-0.5*ratio-cy2*ratio)*camera_up2;
			Vec3f rdir1=c1-orig1,rdir2=c2-orig2;
			rdir1.Normalize();
			rdir2.Normalize();
			Vec3f o=orig1-orig2;
			float t1=(o.x()*rdir2.y()-o.y()*rdir2.x())/(rdir2.x()*rdir1.y()-rdir1.x()*rdir2.y()),t2=(o.y()*rdir2.z()-o.z()*rdir2.y())/(rdir2.y()*rdir1.z()-rdir1.y()*rdir2.z()),t3=(o.x()*rdir2.z()-o.z()*rdir2.x())/(rdir2.x()*rdir1.z()-rdir1.x()*rdir2.z());
			float t=(t1+t2+t3)/3.0f;
			cout<<"t1:"<<t1<<endl;
			cout<<"t2:"<<t2<<endl;
			cout<<"t3:"<<t3<<endl;
			Vec3f dot=orig1+t*rdir1;
			return dot;
		}
		void cal_light_pos(int x1,int y1,int x2,int y2,float cx1,float cy1,float cx2,float cy2,float r1,float r2)
		{
			light_pos=cal_dot_pos(x1,y1,x2,y2,cx1,cy1,cx2,cy2,r1,r2);
		}
		
		/*int calculateCentreAndRadius(Vec3f q1,Vec3f q2, Vec3f q3, Vec3f q4, Vec3f &centre, float &radius)
		{
			vector<float> p1,p2,p3,p4;
			p1.push_back(q1.x());
			p1.push_back(q1.y());
			p1.push_back(q1.z());
			p2.push_back(q2.x());
			p2.push_back(q2.y());
			p2.push_back(q2.z());
			p3.push_back(q3.x());
			p3.push_back(q3.y());
			p3.push_back(q3.z());
			p4.push_back(q4.x());
			p4.push_back(q4.y());
			p4.push_back(q4.z());
    		float a = p1[0] - p2[0], b = p1[1] - p2[1], c = p1[2] - p2[2];
   			float a1 = p3[0] - p4[0], b1 = p3[1] - p4[1], c1 = p3[2] - p3[2];
    		float a2 = p2[0] - p3[0], b2 = p2[1] - p3[1], c2 = p2[2] - p3[2];
    		float A = p1[0] * p1[0] - p2[0] * p2[0];
    		float B = p1[1] * p1[1] - p2[1] * p2[1];
    		float C = p1[2] * p1[2] - p2[2] * p2[2];
    		float A1 = p3[0] * p3[0] - p4[0] * p4[0];
    		float B1 = p3[1] * p3[1] - p4[1] * p4[1];
    		float C1 = p3[2] * p3[2] - p4[2] * p4[2];
    		float A2 = p2[0] * p2[0] - p3[0] * p3[0];
    		float B2 = p2[1] * p2[1] - p3[1] * p3[1];
    		float C2 = p2[2] * p2[2] - p3[2] * p3[2];
    		float P = (A + B + C) /2;
    		float Q = (A1 + B1 + C1) / 2;
    		float R = (A2 + B2 + C2) / 2;

    		// D是系数行列式，利用克拉默法则
    		float D = a*b1*c2 + a2*b*c1 + c*a1*b2 - (a2*b1*c + a1*b*c2 + a*b2*c1);
    		float Dx = P*b1*c2 + b*c1*R + c*Q*b2 - (c*b1*R + P*c1*b2 + Q*b*c2);
    		float Dy = a*Q*c2 + P*c1*a2 + c*a1*R - (c*Q*a2 + a*c1*R + c2*P*a1);
    		float Dz = a*b1*R + b*Q*a2 + P*a1*b2 - (a2*b1*P + a*Q*b2 + R*b*a1);

    	if(D == 0){
        cerr << "四点共面" << endl;
        return -1;
   		 }else{
   		 	centre.Set(Dx/D,Dy/D,Dz/D);
        radius = sqrt((p1[0]-centre.x())*(p1[0]-centre.x()) +
                              (p1[1]-centre.y())*(p1[1]-centre.y()) +
                              (p1[2]-centre.z())*(p1[2]-centre.z()));
        return 0;
    }
}*/

void cal_sphere(Vec3f q1,Vec3f q2, Vec3f q3, Vec3f q4, Vec3f &centre, float &radius)
{
	Vec3f a1=q2-q1,b1=q3-q1,a2=q3-q2,b2=q4-q2;
	Vec3f n1,n2;
	Vec3f::Cross3(n1,a1,b1);
	Vec3f::Cross3(n2,a2,b2);
	n1.Normalize();
	n2.Normalize();
	Vec3f d1,e1,d2,e2;
	Vec3f::Cross3(d1,a1,n1);
	Vec3f::Cross3(e1,b1,n1);
	Vec3f::Cross3(d2,a2,n2);
	Vec3f::Cross3(e2,b2,n2);
	d1.Normalize();
	e1.Normalize();
	d2.Normalize();
	e2.Normalize();
	Vec3f o1=0.5f*(q2-q3),o2=0.5f*(q3-q4);//t1=(o.x()*rdir2.y()-o.y()*rdir2.x())/(rdir2.x()*rdir1.y()-rdir1.x()*rdir2.y())
	float t1=(o1.x()*e1.y()-o1.y()*e1.x())/(e1.x()*d1.y()-d1.x()*e1.y()),t2=(o2.x()*e2.y()-o2.y()*e2.x())/(e2.x()*d2.y()-d2.x()*e2.y());
	Vec3f c1=0.5f*(q1+q2)+t1*d1,c2=0.5f*(q2+q3)+t2*d2;
	Vec3f c=c1-c2;
	float t=(c.x()*n2.y()-c.y()*n2.x())/(n2.x()*n1.y()-n1.x()*n2.y());
	centre=c1+t*n1;
	Vec3f dist1=centre-q1,dist2=centre-q2,dist3=centre-q3,dist4=centre-q4;
	cout<<dist1.Length()<<endl;
	cout<<dist2.Length()<<endl;
	cout<<dist3.Length()<<endl;
	cout<<dist4.Length()<<endl;
	radius=dist1.Length();
}
		void cal_obj(int x11,int y11,int x21,int y21,int x31,int y31,int x41,int y41,int x12,int y12,int x22,int y22, int x32,int y32, int x42,int y42,float cx1,float cy1,float cx2,float cy2,float r1,float r2)
		{
			Vec3f dot1=cal_dot_pos(x11,y11,x12,y12,cx1,cy1,cx2,cy2,r1,r2),dot2=cal_dot_pos(x21,y21,x22,y22,cx1,cy1,cx2,cy2,r1,r2),dot3=cal_dot_pos(x31,y31,x32,y32,cx1,cy1,cx2,cy2,r1,r2),dot4=cal_dot_pos(x41,y41,x42,y42,cx1,cy1,cx2,cy2,r1,r2);
			/*calculateCentreAndRadius(dot1,dot2,dot3,dot4,objC,objR);
			cout<<"objC"<<objC<<"objR"<<objR<<endl;*/
			cal_sphere(dot1,dot2,dot3,dot4,objC,objR);
		}
		void writeFile()
		{
			ofstream out("scene.txt");
			out<<"PerspectiveCamera {"<<endl<<"	center "<<camera_pos<<"	direction "<<camera_dir<<"	up "<<camera_up<<"	angle "<<fovx*180/3.1415926<<endl<<"	cx "<<cx<<endl<<"	cy "<<cy<<endl<<"}"<<endl;
			out<<"Lights {"<<endl<<"	numLights 1"<<endl<<"	PointLight {"<<endl<<"		position "<<light_pos<<"		color "<<light_col<<"		attenuation 1 0 0"<<endl<<"	}"<<endl<<"}"<<endl;
			out<<"Background {"<<endl<<"	color 0 0 0"<<endl<<" ambientLight 0 0 0"<<endl<<"}"<<endl;
		out<<"Materials {"<<endl<<"	numMaterials 1"<<endl<<"	PhongMaterial {"<<endl<<"		diffuseColor 1 1 1"<<endl<<"		specularColor 0.6 0.6 0.6"<<endl<<"exponent 16"<<endl<<"	}"<<endl<<"}"<<endl;
			out<<"Group {"<<endl<<"	numObjects 1"<<endl<<"	MaterialIndex 0"<<endl<<"	Sphere {"<<endl<<"		center "<<objC<<"		radius "<<objR<<endl<<"	}"<<endl<<"}";
		}
	void getPicInfo(int& n,int& _width,int& _height)
	{
		n=pic_num;
		_width=width;
		_height=height;
	}
	void getTrackInfo(Vec3f& _trackC,Vec3f& _axis,float& _trackR,float& angle_velocity,float& _exposure_time)
	{
		_trackC=trackC;
		_axis=axis;
		_trackR=trackR;
		angle_velocity=velocity/trackR;
		_exposure_time=exposure_time;
	}
	private:
		int width,height;
		int pic_num;
		float exposure_time;
		float velocity;
		float fovx,fovy;
		Vec3f camera_dir;
		Vec3f camera_up;
		Vec3f camera_pos;
		Vec3f horizontal;
		Vec3f light_pos;
		Vec3f light_col;
		Vec3f objC;
		float objR;
		Vec3f axis,trackC;
		float trackR;
		Matrix rotate;
		Vec3f T;
		float focus_length;
		float cx,cy;
};

class InverseImage:public Image{
	public:
		InverseImage(int w,int h):Image(w,h)
		{
			isObject=new bool[w*h];
			BRDF_x=new Vec4f[w*h];
		}
		void SetBool(int x,int y,Hit& h,Ray& r,Light* l)
		{
			if(h.getT()!=inf){
				 isObject[y*width+x]=true;
				 Vec3f camera_dir=r.getDirection(),intersectionPoint=h.getIntersectionPoint();
				 Vec3f light_dir,col;
				 float dist_to_light;
				 l->getIllumination(intersectionPoint,light_dir,col,dist_to_light);
				 light_dir*=-1.0f; 
				 Vec3f xAxis(1,0,0),yAxis(0,1,0),zAxis(0,0,1);
				 float camera_phi=camera_dir.Dot3(zAxis)/camera_dir.Length(),light_phi=light_dir.Dot3(zAxis)/light_dir.Length();
				 camera_dir.Set(camera_dir.x(),camera_dir.y(),0);
				 light_dir.Set(light_dir.x(),light_dir.y(),0);
				 float camera_theta=camera_dir.Dot3(xAxis)/camera_dir.Length(),light_theta=light_dir.Dot3(xAxis)/light_dir.Length();
				 BRDF_x[y*width+x].Set(light_theta,light_phi,camera_theta,camera_phi);
			}
		}	
		bool IsObject(int x,int y){
			return isObject[y*width+x];
		}
		Vec4f getBRDFX(int x,int y){
			return BRDF_x[y*width+x];
		}
		~InverseImage(){
			if(isObject) delete[] isObject;
			if(BRDF_x) delete[] BRDF_x;
		}							
	private:
		bool* isObject;
		Vec4f* BRDF_x;
};


int main()
{
	cin>>input_file;
	InputGenerator* IG=InputGenerator::readFile();
	IG->writeFile();
	int width,height,pic_num;
	float angle_velocity,trackR,exposure_time;
	Vec3f trackC,axis;
	IG->getPicInfo(pic_num,width,height);
	IG->getTrackInfo(trackC,axis,trackR,angle_velocity,exposure_time);
	SceneParser* scene=new SceneParser("scene.txt");
	
	Camera *camera=scene->getCamera();
	
	Material* bg=new PhongMaterial(scene->getBackgroundColor());
	
	vector<InverseImage*> pics;
	
	for(int i=0;i<pic_num;i++){
		InverseImage *p=new InverseImage(width,width);
		 pics.push_back(p);
		}
	
	RayTracer* rt=new RayTracer(scene,0,0,0,1,0,0,1e-4,0);
	
	char** outputfile=new char*[pic_num];
	char** img_name=new char*[pic_num];
	for(int i=0;i<pic_num;i++){
		outputfile[i]=new char[10];
		img_name[i]=new char[10];
		sprintf(outputfile[i],"%d",i);
		sprintf(img_name[i],"%d",i);
		strcat(outputfile[i],".tga");
		strcat(img_name[i],".jpg");
	}
	cout<<"exp"<<exposure_time<<endl;
	cout<<"angle_velocity:"<<angle_velocity<<endl;
	//pre-rendering the images
	for(int i=0;i<pic_num;i++){
		//if(i%2==0){
		Vec3f R=trackC-camera->getPos();
		for(int x=0;x<width;x++)
		for(int y=0;y<width;y++){
			Vec2f pos(x*1.0/width,y*1.0/width);
			Ray r=camera->generateRay(pos);
			Hit h(inf,bg,-1.0*r.getDirection()); 
			rt->clearStack(1);
			Light* l=scene->getLight(0);
			pics[i]->SetPixel(x,y,rt->traceRay(r,camera->getTMin(),0,0,h));
			pics[i]->SetBool(x,y,h,r,l);
		}
		pics[i]->SaveTGA(outputfile[i]);
		//}
		camera->moveAroundCamera(axis,trackC,0.5*angle_velocity*exposure_time);
		camera->moveAroundCamera(axis,trackC,0.5*angle_velocity*exposure_time);
	}
	for(int i=0;i<pic_num;i++){
	cv::Mat img=cv::imread(img_name[i]);
	 for(int y=0;y<img.rows;y++)
	 for(int x=0;x<img.cols;x++){
	 if(pics[i]->IsObject(x,y+1920-1080)){
	 	for(int t=-300;t<300;t++)
	 	for(int k=-200;k<200;k++){
	 		x+=t;
	 		y+=k;
	 	float r=img.at<cv::Vec3b>(y,x)[2]/255.0,b=img.at<cv::Vec3b>(y,x)[1]/255.0,g=img.at<cv::Vec3b>(y,x)[0]/255.0;
	 	Vec3f BRDF_y(r,g,b);
	 	x-=t;
	 		y-=k;
	 	if(BRDF_y.Length()>0.15){
	 	Vec4f BRDF_x=pics[i]->getBRDFX(x,y+1920-1080);
	 	if(BRDF.find(BRDF_x)==BRDF.end()){
	 		BRDF.insert(pair<Vec4f,Vec3f>(BRDF_x,BRDF_y));
	 	}
	 	}
	 }
	 }
	}
	}
	if(BRDF.empty()) cout<<"e"<<endl;
	ofstream brdf("final_data.txt");
	for(auto it=BRDF.begin();it!=BRDF.end();it++){
		brdf<<"Theta1&2 and Phi1&2:"<<it->first.x()<<" "<<it->first.y()<<" "<<it->first.z()<<" "<<it->first.w()<<endl;
		brdf<<"RGB:"<<it->second.x()<<" "<<it->second.y()<<" "<<it->second.z()<<endl;
	}
}