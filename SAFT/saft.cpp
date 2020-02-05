#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"saft.h"


void Bscan::load(char fname[128]){
	FILE *fp=fopen(fname,"r");
	if(fp==NULL){
		printf("Cann't find file %s\n",fname);
		printf(" --> abort process.");
		exit(-1);
	}

	char cbff[128];

	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&Ny);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf, %lf, %d\n",&t0,&dt,&Nt);
	fgets(cbff,128,fp);
	printf("Ny,Nt=%d %d\n",Ny,Nt);

	Bscan::mem_alloc();

	int i,j;
	double tmp;
	for(i=0; i<Ny;i++){
	for(j=0; j<Nt;j++){
		//fscanf(fp,"%lf\n",amp[i]+j);
		fscanf(fp,"%lf\n",&tmp);
		amp[i][j]=tmp;
	}
	}

	fclose(fp);
};

void Bscan::mem_alloc(){
	double *pt=(double *)malloc(sizeof(double)*Nt*Ny);
	amp=(double **)malloc(sizeof(double *)*Ny);
	int i;
	for(i=0;i<Nt*Ny;i++) pt[i]=0.0;
	for(i=0;i<Ny;i++) amp[i]=pt+Nt*i;

	amp_sum=(double *)malloc(sizeof(double)*Nt);
	for(i=0;i<Nt;i++) amp_sum[i]=0.0;

};

void Bscan::stack_Ascans(double *tofs){
	int i,j,jd;
	double idly;
	tof=tofs[0];
	ntof=int((tof-t0)/dt);

	for(j=0;j<Nt;j++) amp_sum[j]=0.0;

	for(i=0;i<Ny;i++){
		idly=int((tofs[i]-tof)/dt);
		for(j=0;j<Nt;j++){
			jd=j+idly;
			if(jd<0) continue;
			if(jd>Nt-1) continue;
			amp_sum[j]+=amp[i][jd];
		}
//		amp_sum[i]/=Ny;
	}
};
double Bscan::get_Asum(int inc_ntof){
	int mtof=inc_ntof+ntof;
	if(mtof<0) return(0);
	if(mtof>Nt-1) return(0);
	return(amp_sum[mtof]);
};
void Bscan::fwrite_Ascan(char fname[128]){
	FILE *fp=fopen(fname,"w");
	int i;
	for(i=0;i<Nt;i++){
		fprintf(fp,"%lf %lf\n",t0+dt*i,amp_sum[i]);
	}
	fprintf(fp,"\n");
	fprintf(fp,"%lf %lf\n",tof,5.0);
	fprintf(fp,"%lf %lf\n",tof,-5.0);
	fclose(fp);
};
void Bscan::fwrite_Bscan(){
	FILE *fp=fopen("bscan.out","w");
	int i,j;

	fprintf(fp,"#\n");
	fprintf(fp,"%d\n",Ny);
	fprintf(fp,"#\n");
	fprintf(fp,"%lf, %lf, %d\n",t0,dt,Nt);
	fprintf(fp,"#\n");
	for(j=0;j<Ny;j++){
	for(i=0;i<Nt;i++){
		fprintf(fp,"%lf\n",amp[j][i]);
	}
	}
	fclose(fp);
};


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
//---------------------------------------------------
class IMG{
	public:
		double Xa[2],Xb[2],dx[2],Wd[2];
		int Ndiv[2],Ng[2]; 
		int npnt;
		double **V; 
		IMG();
		void set_Xa(double x, double y);
		void set_Wd(double W, double H);
		void set_Ng(int nx, int ny);
		double get_xcod(int l);
		double get_ycod(int l);
		void set_V(int l,double val);
		void fwrite_V(char fn[128]);
	private:
		void mem_alloc();
};
IMG::IMG(){
	Xa[0]=0.0; Xa[1]=0.0;
	Wd[0]=1.0; Wd[1]=1.0;
};
void IMG::set_Xa(double x, double y){
	Xa[0]=x; Xa[1]=y;
};
void IMG::set_Wd(double W, double H){
	Wd[0]=W; Wd[1]=H;
	Xb[0]=Xa[0]+W;
	Xb[1]=Xa[1]+H;
};
void IMG::set_Ng(int nx, int ny){
	Ng[0]=nx; Ng[1]=ny;
	for(int i=0;i<2;i++){
		Ndiv[i]=Ng[i]-1;
		dx[i]=0.0;
		if(Ndiv[i]>0) dx[i]=Wd[i]/Ndiv[i];
	}
	IMG::mem_alloc();
	npnt=Ng[0]*Ng[1];
};
void IMG::mem_alloc(){
	double *pt=(double *)malloc(sizeof(double)*Ng[0]*Ng[1]);
	V=(double **)malloc(sizeof(double *)*Ng[0]);

	int i;
	for(i=0;i<Ng[0]*Ng[1];i++) pt[i]=0.0;
	for(i=0;i<Ng[0];i++) V[i]=pt+Ng[1]*i;
};

double IMG::get_xcod(int l){
	int i=l/Ng[1];
	return(Xa[0]+i*dx[0]);
};
double IMG::get_ycod(int l){
	int j=l%Ng[1];
	return(Xa[1]+j*dx[1]);
};
void IMG::set_V(int l,double val){
	int i=l/Ng[1];
	int j=l%Ng[1];
	V[i][j]=val;
};
void IMG::fwrite_V(char fname[128]){
	FILE *fp=fopen(fname,"w");
	int i,j;
	fprintf(fp,"# Xa[0]  Xa[1]\n");
	fprintf(fp,"%lf %lf\n",Xa[0],Xa[1]);
	fprintf(fp,"# Xb[0]  Xb[1]\n");
	fprintf(fp,"%lf %lf\n",Xb[0],Xb[1]);
	fprintf(fp,"# Ng[0]  Ng[1]\n");
	fprintf(fp,"%d %d\n",Ng[0],Ng[1]);
	fprintf(fp,"# Values\n");
	for(i=0;i<Ng[0];i++){
	for(j=0;j<Ng[1];j++){
		fprintf(fp," %lf\n",V[i][j]);
	}
	}
	fclose(fp);
};

int main(){

//	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	char fdat_in[128]="wave_in4_2.dat";	// B-scan data (incident wave) 
	char fdat_sc[128]="wave_sc4_2.dat";	// B-scan data (scattered wave)
	char fnout[128]="saft4_2.out";		// Output file (SAFT image)
	char ftof_in[128]="tof_in.out",mode[2]="w"; // Travel time plot 
	char ftof_sc[128]="tof_sc.out"; // Travel time plot 
//	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	Bscan bwv_in;
	Bscan bwv_sc;
	bwv_in.load(fdat_in);
	bwv_sc.load(fdat_sc);

//		Plate Geometry 
	double ht=12.0;	// plate thickness
	double zb=-ht,zt=0.0;
	Plate PL(zb, zt);	// set plate surfaces

	double xin[2],xsc[2];
//	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//xin[0]=51.5; xin[1]= 0.0;	// source point
	xin[0]=52.5; xin[1]= 0.0;	// source point
	//xsc[0]=-7.0; xsc[1]=-4.0;	// scatterer location (used for TOF plotting)
	xsc[0]=-7.0; xsc[1]=-4.0;	// scatterer location (used for TOF plotting)
	//PL.rcv.set_xr1(20.0, 0.0); // receiver array (start)
	//PL.rcv.set_xr2( 0.0, 0.0); // receiver array (end)
	PL.rcv.set_xr1(22.5, 0.0); // receiver array (start)
	PL.rcv.set_xr2( 2.5, 0.0); // receiver array (end)
	PL.rcv.gen(81);	// generate reciver array (given number of receivers)
	double cT=3.036;	// wave velocity [km/s]
	double t_wedge=8.162;	// travel time in wedge [micro sec]
//	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	double tin;	// travel time from source to target
	int down=-1, up=1; // direction of wave upon excitation (up/down)
	int iref,jref;

//	----------- Draw Travel Time Curves -------------
	PL.set_src(xin[0],xin[1]);	// Incident point
	int nref=6;
	FILE *fp=fopen(ftof_in,"w");
	fprintf(fp,"# ncrv, nrec\n");
	fprintf(fp," %d, %d\n",nref/2+1,PL.rcv.nrec);	
	fclose(fp);

	sprintf(mode,"a");
	for(iref=0;iref<=nref;iref+=2){
		PL.rcv.set_tof(t_wedge);
		PL.TOFs(down,iref,cT);
		fp=fopen(ftof_in,"a");
		if(iref==0){
			fprintf(fp," %d\n",iref);	
		}else{
			fprintf(fp," %d\n",iref-1);	
		}
		fclose(fp);
		PL.write_tx(ftof_in,mode);
	}

	int ncrv;
	nref=3;
	ncrv=(nref+1)*(nref/2+1+(nref+1)/2);
	printf("ncrv=%d\n",ncrv);
	fp=fopen(ftof_sc,"w");
	fprintf(fp,"# ncrv, nrec\n");
	fprintf(fp," %d, %d\n",ncrv,PL.rcv.nrec);	
	fclose(fp);
	PL.set_src(xsc[0],xsc[1]);	// Crack tip 
	ncrv=0;
	for(iref=0;iref<=nref;iref++){
		tin=PL.TOF(xin,xsc,down,iref,cT);
		for(jref=0;jref<=nref;jref+=2){
			ncrv++;
			fp=fopen(ftof_sc,"a");
			fprintf(fp," %d, %d\n",iref,jref);	
			fclose(fp);
			PL.rcv.set_tof(tin+t_wedge);
			PL.TOFs(up,jref,cT);
			PL.write_tx(ftof_sc,mode);
		}
	}
	for(iref=0;iref<=nref;iref++){
		tin=PL.TOF(xin,xsc,down,iref,cT);
		for(jref=1;jref<=nref;jref+=2){
			ncrv++;
			fp=fopen(ftof_sc,"a");
			fprintf(fp," %d, %d\n",iref,jref);	
			fclose(fp);
			PL.rcv.set_tof(tin+t_wedge);
			PL.TOFs(down,jref,cT);
			PL.write_tx(ftof_sc,mode);
		}
	}
	printf("ncrv=%d\n",ncrv);

//	--------- Waveform Stacking ----------

	char fascn_in[128]="ascan_in.out";	// Reference Wave 
	PL.set_src(xin[0],xin[1]);	// Incident point
	PL.rcv.set_tof(t_wedge);
	PL.TOFs(down,1,cT);
	bwv_in.stack_Ascans(PL.rcv.tf);
	bwv_in.fwrite_Ascan(fascn_in);

	char fascn_sc[128]="ascan_sc.out";
	int dir_in=down; // wave direction upon generation
	int nref_in=2;	// number of reflection
	int dir_sc=up; // wave direction upon generation
	int nref_sc=0;	// number of reflection

	int k,npnt=11;
	double Wt=0.6;
	int Nwin=ceil(Wt/bwv_in.dt);
	double cor;

//	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	IMG IM;
	IM.set_Xa(-10.0-5.0,-12.0);	// Imaging area (lower left corner)
	IM.set_Wd(20.0,12.0);		// Width x Height [mm]
	IM.set_Ng(121,81);		// Number of pixels (width x height)
	
	//IM.set_Xa(xsc[0],xsc[1]);	// Imaging area (lower left corner)
	//IM.set_Wd(0.0,0.0);		// Width x Height [mm]
	//IM.set_Ng(1,1);		// Number of pixels (width x height)
//	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	double a_in,a_sc;
	for(k=0;k<IM.npnt;k++){
		xsc[0]=IM.get_xcod(k);
		xsc[1]=IM.get_ycod(k);

		PL.set_src(xsc[0],xsc[1]);	// Crack tip 
		tin=PL.TOF(xin,xsc,dir_in,nref_in,cT);
		PL.rcv.set_tof(tin+t_wedge);
		PL.TOFs(dir_sc,nref_sc,cT);
		bwv_sc.stack_Ascans(PL.rcv.tf);
	//	bwv_sc.fwrite_Ascan(fascn_sc);
		cor=0.0;
		for(int m=0;m<Nwin;m++){
			a_in=bwv_in.get_Asum(m);
			a_sc=bwv_sc.get_Asum(m);
			//cor+=(bwv_sc.get_Asum(m)*bwv_in.get_Asum(m));
			cor+=(a_in*a_sc);

		};
		printf("xs=(%lf, %lf), cor=%lf\n",xsc[0],xsc[1],cor);
		IM.set_V(k,cor);
	}
	IM.fwrite_V(fnout);
//	--------------------------------------

/* uncomment this section to draw ray paths
	int ud;
	double xs[2]={0.0,-5.0};
	double xr[2]={20.0,zt};
	nref=1; ud=1;
	PL.path(xs,xr,ud,nref);
	printf("\n");
	nref=3; ud=-1;
	PL.path(xs,xr,ud,nref);

	printf("\n");
	printf("%lf %lf\n",-5.0,PL.zt);
	printf("%lf %lf\n", 25.0,PL.zt);
	printf("\n");
	printf("%lf %lf\n",-5.0,PL.zb);
	printf("%lf %lf\n", 25.0,PL.zb);
*/

	return(0);
};
