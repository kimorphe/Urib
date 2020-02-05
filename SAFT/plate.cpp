#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"saft.h"


Plate::Plate(double z_bottom, double z_top){
	zb=z_bottom;
	zt=z_top;
	zs[0]=zb;
	zs[1]=zt;
	ht=zt-zb;
};
void Plate::set_src(double x, double y){
	xs[0]=x; xs[1]=y;
};
double Plate::mirror(
	double zcod, 
	int ud,	// -1:down, +1:up
	int nref // number of reflection
){
	if(nref==0) return(zcod);
	double zimg;
	if(ud==-1){	// downward mirror operationo
		zimg=2.*zb-zcod;
	}else if(ud==1){
		zimg=2.*zt-zcod;
	}
	nref--;
	if(nref >0) zimg=mirror(zimg,-ud,nref);
	return(zimg);
};

double Plate::path(double xs[2], double xr[2], int ud, int nref){

	double xsd[2],rx,ry;
		
	xsd[0]=xs[0];
	xsd[1]=mirror(xs[1],ud,nref);

	rx=xr[0]-xsd[0];
	ry=xr[1]-xsd[1];

	double xi,zi;
	int isgn,isrf;

	double eps=1.e-08;
	if(xsd[1]==zt) xsd[1]+=eps; 
	if(xsd[1]==zb) xsd[1]-=eps; 

	if(xsd[1]> zt){
		zi=zt;
		isgn=1;	
		isrf=1;
	}else{
		zi=zb;
		isgn=-1;
		isrf=0;
	};

	printf("%lf %lf\n",xr[0],xr[1]);
	double dz;
	for(int i=0;i<nref;i++){
		xi=(xsd[1]-zi)*xr[0]+(zi-xr[1])*xs[0];
		dz=xsd[1]-xr[1];
		if(abs(dz)<eps) dz=eps;
		xi/=dz;
		printf("%lf %lf\n",xi,zs[isrf%2]);
		zi+=(isgn*ht);
		isrf++;
	};
	printf("%lf %lf\n",xs[0],xs[1]);
	return(sqrt(rx*rx+ry*ry));
};

double Plate::TOF(double xs[2], double xr[2], int ud, int nref, double c){

	double xsd[2],rx,ry;
		
	xsd[0]=xs[0];
	xsd[1]=mirror(xs[1],ud,nref);

	rx=xr[0]-xsd[0];
	ry=xr[1]-xsd[1];

	return(sqrt(rx*rx+ry*ry)/c);
};
void Plate::TOFs(int ud, int nref, double c)
{
	double xsd[2],rx,ry;
	int i;

	xsd[0]=xs[0];
	xsd[1]=mirror(xs[1],ud,nref);
	for(i=0;i<rcv.nrec; i++){
		rx=xsd[0]-rcv.xr[i];
		ry=xsd[1]-rcv.yr[i];
		rcv.tf[i]+=sqrt(rx*rx+ry*ry)/c;
	}
};
void Plate::write_tx(char fname[128],char mode[2]){
	FILE *fp=fopen(fname,mode);
	for(int i=0;i<rcv.nrec; i++){
		fprintf(fp,"%lf %lf %lf\n",rcv.tf[i],rcv.xr[i],rcv.yr[i]);
	}
	//fprintf(fp,"\n");
	fclose(fp);
};

//---------------------------------------------------
void Recs::set_xr1(double x, double y){
	xr1[0]=x;
	xr1[1]=y;
};
void Recs::set_xr2(double x, double y){
	xr2[0]=x;
	xr2[1]=y;
};
void Recs::gen(int n){
	nrec=n;
	dxr[0]=0.0; dxr[1]=0.0;
	if(n>1){
		dxr[0]=(xr2[0]-xr1[0])/(nrec-1);
		dxr[1]=(xr2[1]-xr1[1])/(nrec-1);
	};

	xr=(double *)malloc(sizeof(double)*nrec);
	yr=(double *)malloc(sizeof(double)*nrec);
	tf=(double *)malloc(sizeof(double)*nrec);

	for(int i=0;i<nrec;i++){
		xr[i]=xr1[0]+dxr[0]*i;
		yr[i]=xr1[1]+dxr[1]*i;
		tf[i]=0.0;
	};
};
void Recs::print_cods(){
	puts("### Receiver Coordinates");
	for(int i=0;i<nrec;i++){
		printf("%lf %lf\n",xr[i],yr[i]);
	};
};
void Recs::set_tof(double val){
	for(int i=0;i<nrec;i++) tf[i]=val;
};
