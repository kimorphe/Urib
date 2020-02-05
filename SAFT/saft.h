class Recs{
	public:
		int nrec; // number of receiver points
		double xr1[2]; // start
		double xr2[2]; // end --> equi-spaced nrec receivers will be generated
		void set_xr1(double x, double y);
		void set_xr2(double x, double y);
		void gen(int n);
		double *xr,*yr,*tf;
		double dxr[2];
		void print_cods();
		void set_tof(double val);
	private:
};
class Bscan{
	public:
		char fname[128];
		int Ny;
		double y0,dy;
		int Nt;
		double t0,dt;
		double **amp;
		double *amp_sum;
		void stack_Ascans(double* delay);
		void fwrite_Ascan(char fname[128]);
		void fwrite_Bscan();
		void load(char fn[128]);
		void mem_alloc();
		double tof;
		int ntof;
		double get_Asum(int inc_ntof);
	private:
};

class Plate{
	public:
		double zb;	//bottom surface
		double zt;	//top surface
		double ht;	//plate thickness
		double zs[2];
		Plate(double z_bottom, double z_top);
		double mirror(double zf, int ud, int nref);
		double path(double xs[2],double xr[2], int ud, int nref);
		double TOF(double xs[2],double xr[2], int ud, int nref,double c);
		void TOFs(int ud, int nref,double c);
		void write_tx(char fn[128],char mode[2]);
		Recs rcv;	// receiver points
		double xs[2];	// source point
		void set_src(double x,double y);
	private:
};

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
