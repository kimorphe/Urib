#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"saft.h"

int main(){

	char fgeom[128]="geom.inp";
	FILE *fp=fopen(fgeom,"r");

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
	double xin[2]={0.0,0.0};
	double xsc[2]={0.0,0.0};
	double xr1[2]={0.0,0.0};
	double xr2[2]={0.0,0.0};
	double cT=3.036;	// wave velocity [km/s]
	double t_wedge=8.162;	// travel time in wedge [micro sec]
	int nrec,nref=6;

	xsc[0]=-7.0; xsc[1]=-4.0;	// scatterer location (used for TOF plotting)

	char cbff[128];
	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",&ht);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",xin);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",xr1);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",xr2);
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&nrec);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",&cT);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",&t_wedge);
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&nref);
	fclose(fp);

	printf("t_wedge=%lf nref=%d\n",t_wedge,nref);

	double zb=-ht,zt=0.0;

	Plate PL(zb, zt);	// set plate surfaces

//	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//xin[0]=51.5; xin[1]= 0.0;	// source point
	//xin[0]=52.5; xin[1]= 0.0;	// source point

	//PL.rcv.set_xr1(20.0, 0.0); // receiver array (start)
	//PL.rcv.set_xr2( 0.0, 0.0); // receiver array (end)

	//PL.rcv.set_xr1(22.5, 0.0); // receiver array (start)
	//PL.rcv.set_xr2( 2.5, 0.0); // receiver array (end)
	//PL.rcv.gen(81);	// generate reciver array (given number of receivers)

	PL.rcv.set_xr1(xr1[0],xr1[1]); // receiver array (start)
	PL.rcv.set_xr2(xr2[0],xr2[1]); // receiver array (end)
	PL.rcv.gen(nrec);	// generate reciver array (given number of receivers)
//	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	double tin;	// travel time from source to target
	int down=-1, up=1; // direction of wave upon excitation (up/down)
	int iref,jref;

//	----------- Draw Travel Time Curves -------------
	PL.set_src(xin[0],xin[1]);	// Incident point
	fp=fopen(ftof_in,"w");
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
	fclose(fp);


//	--------- Waveform Stacking ----------

	char fascn_in[128]="ascan_in.out";	// Reference Wave 
	PL.set_src(xin[0],xin[1]);	// Incident point
	PL.rcv.set_tof(t_wedge);
	PL.TOFs(down,1,cT);
	double amax;
	amax=bwv_in.stack_Ascans(PL.rcv.tf);
	bwv_in.fwrite_Ascan(fascn_in);

	printf("amax=%lf\n",amax);



	char fascn_sc[128]="ascan_sc.out";
	int dir_in=down; // wave direction upon generation
	int nref_in=2;	// number of reflection
	int dir_sc=up; // wave direction upon generation
	int nref_sc=0;	// number of reflection

	int k,npnt=11;
	double Wt=0.6;
	int Nwin=ceil(Wt/bwv_in.dt);
	double cor;


	fp=fopen("saft.inp","r");
	fgets(cbff,128,fp);
	fgets(fnout,128,fp);
	puts(fnout);

	double Xa[2],Xb[2],Wd[2];
	int Ng[2];
	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf\n",Xa,Xa+1);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf\n",Xb,Xb+1);
	fgets(cbff,128,fp);
	fscanf(fp,"%d %d\n",Ng,Ng+1);
	Wd[0]=Xb[0]-Xa[0];
	Wd[1]=Xb[1]-Xa[1];
	printf("Ng=%d %d\n",Ng[0],Ng[1]);

	fclose(fp);
//	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	IMG IM;
	//IM.set_Xa(-10.0-5.0,-12.0);	// Imaging area (lower left corner)
	//IM.set_Wd(20.0,12.0);		// Width x Height [mm]
	//IM.set_Ng(121,81);		// Number of pixels (width x height)
	IM.set_Xa(Xa[0],Xa[1]);	// Imaging area (lower left corner)
	IM.set_Wd(Wd[0],Wd[1]);		// Width x Height [mm]
	IM.set_Ng(Ng[0],Ng[1]);		// Number of pixels (width x height)
	
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
