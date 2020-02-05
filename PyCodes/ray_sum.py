#! /home/kazushi/anaconda3/bin/python

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import copy

import wv1d
import sample
import rays2

from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar
from mpl_toolkits.axes_grid1 import ImageGrid


class WV2D:
    def __init__(self,Nx,Ny,dx=1.0,dy=1.0,x0=0.0,y0=0):

        self.Nx=Nx;
        self.Ny=Ny;

        self.nfile=Nx*Ny;

        self.dx=dx;
        self.dy=dy;

        self.x0=x0;
        self.y0=y0;

        self.Wd=np.array([dx*(Nx-1),dy*(Ny-1)]);

        self.Xa=[x0,y0];
        self.Xb=self.Xa+self.Wd;

        self.xcod=self.Xa[0]+dx*np.arange(Nx);
        self.ycod=self.Xa[1]+dy*np.arange(Ny);

        self.AMAP=np.array([])


        self.tLim=[0,80];

    def load_1col(self,dirc="./",head="scope_",tail=".csv"):
        nums=np.arange(self.nfile);
        isum=0;
        for num in nums:
            fname=dirc+"/"+head+str(num)+tail
            print(fname)
            fp=open(fname,"r")
            dat=fp.readline()
            dat=dat.strip().split(",")
            t1=float(dat[0]);
            dt=float(dat[1]);
            Nt=int(fp.readline())

            amp=np.array(fp.readlines());
            amp=amp.astype(np.float32)
            fp.close()
            bias=np.mean(amp)
            if isum==0:
                Amp=np.zeros(len(nums)*Nt)

            Amp[isum*Nt:(isum+1)*Nt]=(amp-bias)
            isum+=1;

        self.Nt=Nt;
        time=np.array(range(Nt))*dt+t1;
        time=np.array(time)*1.e06
        self.t1=time[0];
        self.t2=time[-1];
        self.dt=dt*1.e06;
        
        Amp=Amp*2.0;    # vol --> m/s
        Amax=np.max(np.abs(Amp))
        Amp=np.reshape(Amp,[self.Ny,self.Nx,self.Nt])/Amax;

        self.time=time
        self.amp=Amp

    def load(self,dirc="./",head="scope_",tail=".csv"):
        nums=range(self.nfile);
        Amp=np.array([]);
        isum=0;
        for num in nums:
            fname=dirc+"/"+head+str(num)+tail
            print(fname)
            fp=open(fname,"r")

            
            Dat=fp.read()
            Dat=Dat.replace("\n",",")
            Dat=Dat.split(",")
            Dat=np.array(Dat[0:-1])
            Dat=Dat.astype(np.float32)

            Nt,=np.shape(Dat)
            Nt=int(Nt/2);
            indx=np.array(range(Nt))
            time=Dat[indx*2]
            amp=Dat[indx*2-1]
            if isum==0:
                Amp=np.zeros(Nt*len(nums))
            fp.close()
            dt=time[1]-time[0]
            bias=np.mean(amp)
            amp=amp-bias
            Amp[isum*Nt:(isum+1)*Nt]=np.array(amp)
            isum+=1;

        self.Nt=Nt;
        time=np.array(time)*1.e06
        self.t1=time[0];
        self.t2=time[-1];
        self.dt=dt*1.e06;
        
        Amp=Amp*2.0;    # vol. --> m/s
        Amax=np.max(np.abs(Amp))
        Amp=np.reshape(Amp,[self.Ny,self.Nx,self.Nt])/Amax;
        self.time=time
        self.amp=Amp

    def snapshot(self,ax,time):
        itm=np.argmin(np.abs(self.time-time))
        ext=[self.Xa[0],self.Xb[0],self.Xa[1],self.Xb[1]]
        Amax=np.max(np.abs(self.amp[:,:,:]))
        ax.imshow(self.amp[:,:,itm]/Amax,interpolation="bilinear",extent=ext,cmap="jet",vmin=-0.4,vmax=0.4,origin="lower");
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        txt="t="+str(time)+"[$\mu s$]"
        ax.set_title(txt)

    def multi_bscans(self,fig,x1,x2,nrows,ncols):
        grid=ImageGrid(fig, 111,
                nrows_ncols=(nrows,ncols),
                axes_pad=0.20,
                share_all=True,
                cbar_location="bottom", #right/left, top/bottom
                cbar_mode="single",
                cbar_size="8%",
                cbar_pad=0.25,
                )
        isum=0;
        xx=np.linspace(x1,x2,nrows*ncols)
        ext=[ self.t1,self.t2,self.Xa[1],self.Xb[1]];
        nfigs=nrows*ncols
        fsz=10
        for ax in grid:
            ix=np.argmin(np.abs(self.xcod-xx[isum]))
            #im=ax.imshow(self.amp[:,:,itm]/Amax,interpolation="bilinear",extent=ext,cmap="jet",vmin=-0.4,vmax=0.4,origin="lower");
            im=ax.imshow(self.amp[:,ix,:],cmap="jet",origin="lower",aspect="auto",interpolation="bilinear",extent=ext,vmin=-0.4,vmax=0.4)
            ax.set_xlim(self.tLim)
            txt='{:,.1f}'.format(self.xcod[ix])
            txt=" $x$="+txt+"mm"
            ax.text(self.tLim[0],self.Xb[1],txt,fontsize=fsz,horizontalalignment="left",verticalalignment="bottom")
            if isum >=nfigs-ncols:
                ax.set_xlabel("time [$\mu$s]",fontsize=fsz)
            if isum%ncols ==0:
                ax.set_ylabel("y [mm]",fontsize=fsz)
            isum+=1;
        ax.cax.colorbar(im)
        ax.cax.toggle_label(True)


    def multi_snapshots(self,fig,t1,t2,nrows,ncols):
        Amax=np.max(np.abs(self.amp[:,:,:]))
        ext=[self.Xa[0],self.Xb[0],self.Xa[1],self.Xb[1]]

        grid=ImageGrid(fig, 111,
                nrows_ncols=(nrows,ncols),
                axes_pad=0.20,
                share_all=True,
                cbar_location="bottom", #right/left, top/bottom
                cbar_mode="single",
                cbar_size="5%",
                cbar_pad=0.25,
                )
        isum=0;
        times=np.linspace(t1,t2,nrows*ncols)
        fsz=10; # fontsize
        nfig=nrows*ncols;
        for ax in grid:
            #im=ax.imshow(np.random.random((10,10)),vmin=0,vmax=1)
            itm=np.argmin(np.abs(self.time-times[isum]))
            im=ax.imshow(self.amp[:,:,itm]/Amax,interpolation="bilinear",extent=ext,cmap="jet",vmin=-0.4,vmax=0.4,origin="lower");
            txt='{:,.1f}'.format(self.time[itm])
            txt=" $t=$"+txt+"$\mu$s"
            ax.text(self.Xa[0],self.Xb[1],txt,fontsize=fsz,horizontalalignment="left",verticalalignment="bottom")
            if isum >= nfig-ncols:
                ax.set_xlabel("x[mm]",fontsize=fsz)
            if isum% ncols ==0:
                ax.set_ylabel("y[mm]",fontsize=fsz)
            isum+=1;
        ax.cax.colorbar(im)
        ax.cax.toggle_label(True)

    def mean_bscan(self,x1=0,x2=0):
        if x1*x1+x2*x2 == 0:
            i1=0;
            i2=self.Nx
        else:
            i1=np.argmin(abs(self.xcod-x1))
            i2=np.argmin(abs(self.xcod-x2))

        Bmp=np.mean(self.amp[:,i1:i2+1,:],1)

        #ext=[ self.t1,self.t2,self.Xa[1],self.Xb[1]];
        #im=ax.imshow(Bmp,cmap="jet",origin="lower",aspect="auto",interpolation="bilinear",extent=ext,vmin=-0.1,vmax=0.1)
        #ax.set_xlim(self.tLim)
        #ax.set_ylabel("y [mm]");
        #ax.set_xlabel("time [$\mu$s]");
        #ax.set_title("Travel-time plot of mean waveforms");

        #ax_div=make_axes_locatable(ax)
        #cax=ax_div.append_axes("right",size="5%",pad="2%")
        #cb=colorbar(im,cax=cax)

        bwv=wv1d.WV1D(self.Ny,dy=self.dy,y0=self.y0);
        bwv.amp=Bmp;
        bwv.time=self.time
        bwv.t1=bwv.time[0];
        bwv.t2=bwv.time[-1];
        bwv.dt=self.dt;
        bwv.Nt=self.Nt
        bwv.Ny=self.Ny

        return(bwv)

    def mean_bscan_FFT(self,ax,x1=0,x2=0):
        if x1*x1+x2*x2 == 0:
            i1=0;
            i2=self.Ny
        else:
            i1=np.argmin(abs(self.xcod-x1))
            i2=np.argmin(abs(self.xcod-x2))

        Bmp=np.mean(self.amp[:,i1:i2+1,:],1)
        BMP=np.fft.fft(Bmp,axis=1);
        fmax=1/self.dt
        ext=[0,fmax,self.Xa[1],self.Xb[1]]

        BMP=np.fft.ifft(BMP,axis=0);

        
        Bmax=np.max(np.abs(BMP))
        #im=ax.imshow(np.abs(BMP)/Bmax,cmap="jet",origin="lower",aspect="auto",interpolation="bilinear",extent=ext,vmin=0,vmax=0.5)

        BMP[0:int(self.Ny/2),0:int(self.Nt/2)]=0.0;
        BMP[int(self.Ny/2):-1,int(self.Nt/2):-1]=0.0;
        Bt=np.fft.fft(BMP,axis=0);
        Bt=np.fft.ifft(Bt,axis=1)
        Btmax=np.max(np.real(Bt))
        im=ax.imshow(np.real(Bt)/Btmax,cmap="jet",origin="lower",aspect="auto",interpolation="bilinear",extent=ext) #,vmin=-0.5,vmax=0.5)
        #ax.set_xlim([0,10])

        ax_div=make_axes_locatable(ax)
        cax=ax_div.append_axes("right",size="5%",pad="2%")
        cb=colorbar(im,cax=cax)
        ax.set_ylabel("y [mm]");
        ax.set_xlabel("frequency [MHz]");
        ax.set_title("Frequency spectrum")

    def bscan(self,ax,xcod):
        ix=int((xcod-self.x0)/self.dx);
        #Bmp=np.mean(self.amp[i1:i2+1,:,:],1)
        Bmp=np.mean(self.amp,1)
        ext=[ self.t1,self.t2,self.Xa[1],self.Xb[1]];
        ax.imshow(self.amp[:,ix,:]-Bmp,cmap="jet",origin="lower",aspect="auto",interpolation="bilinear",extent=ext,vmin=-0.7,vmax=0.5)
        ax.set_xlim(self.tLim);
        ax.set_xlabel("time [micro sec]");
        ax.set_ylabel("y [mm]");
        ax.text(100, self.Xb[1],"Difference waveform",verticalalignment="top",horizontalalignment="right");
        ax.text(100, self.Xa[1],"x="+str(self.xcod[ix])+"[mm]",verticalalignment="bottom",horizontalalignment="right");

    def ascans(self,ax):
        Bmp=np.sum(self.amp,1)
        Bmp/=np.max(np.abs(Bmp))
        ys=np.array([5,10,15,20,25,30])
        nums=ys/dy;
        A0=0.0;
        dA=1.0
        for k in nums.astype("int"):
            ax.plot(self.time,Bmp[k,:]+A0)
            A0+=dA
        ax.grid(True)
        t1=self.time[0];
        t2=self.time[-1];
        ax.set_xlim([t1,t2])

if __name__=="__main__":

    # Define sample geometry
    SMP=sample.SAMPLE(4)
    ntip=SMP.ntip
    nroot=SMP.nroot;
    xtip=np.array([SMP.PT[ntip,0], SMP.PT[ntip,1]])
    xroot=np.array([SMP.PT[nroot,0], SMP.PT[nroot,1]])
    h1=SMP.h1;  # top surface
    h2=SMP.h2;  # bottom surface
    ht=h1-h2;   # plate thickness [mm]
    cT=3.036    # Shear wave velocity [km/s]
    cR=2.910    # Rayleigh surface wave velocity [km/s]
    t0=8.162 # delay time

    xs=np.array([51.0,0.0])  # source point (incident point)
    SMP.set_src(xs[0],xs[1])

    recs=rays2.RECs() # use RECs class to handle mutiple receiver point
    Nrec=81;  # the number of receiver pointsｙ
    xr1=20.0;  # x-coordinate ( start )
    xr2=0.0;  # x-coordinate ( end )
    yr=h1;    # y-coordinate
    recs.setup(xr1,xr2,Nrec,h1) # setup reception point data
    SMP.set_rec(xr1,xr2,yr,Nrec)
    #------  Measurement Grid --------
    Nx=5; Ny=81;  
    dx=0.25; dy=-0.25;
    nfile=Nx*Ny;
    x0=-(Nx-1)*0.5*dx; y0=xr1;
    #dir_name="70_0_0211"; du=0.1
    #dir_name="70_2_0212"; du=0.1
    dir_name="70_4_0211"; du=0.1
    one_col=False
    #---------------------------------

    wv=WV2D(Nx=Nx,Ny=Ny,dx=dx,dy=dy,x0=x0,y0=y0)
    if one_col:
        wv.load_1col(dirc=dir_name)
    else:
        wv.load(dirc=dir_name)
    fig1=plt.figure(figsize=[8,4])
    ax1=fig1.add_subplot(111) 

    bwv=wv.mean_bscan()
    bwv.show_xt(ax1)

    wvp=copy.deepcopy(bwv)
    wvm=copy.deepcopy(bwv)
    (wvp.amp,wvm.amp)=bwv.split()

    #wvm.Lmax_normalize()
    #wvp.Lmax_normalize()
    #-----------------------------------
    fig1=plt.figure(figsize=[8,8])
    ax1=fig1.add_subplot(311) 
    ax2=fig1.add_subplot(312) 
    ax3=fig1.add_subplot(313) 
    fsz=12
    #cmap="bone"
    cmap="gray"
    bwv.show_xt(ax1,cmp=cmap,Vmin=-0.2,Vmax=0.2,fsz=fsz)
    wvp.show_xt(ax2,cmp=cmap,Vmin=-0.2,Vmax=0.2,fsz=fsz)
    wvm.show_xt(ax3,cmp=cmap,Vmin=-0.2,Vmax=0.2,fsz=fsz)
    ax3.set_xlabel("time [$\mu$s]",fontsize=fsz)
    ax1.set_title(dir_name,fontsize=fsz+2)
    #ax1.text(18,2.0,"total field",color="w",horizontalalignment="left",fontsize=fsz)
    #ax2.text(18,2.0,"incident field",color="w",horizontalalignment="left",fontsize=fsz)
    #ax3.text(18,2.0,"scattered field",color="w",horizontalalignment="left",fontsize=fsz)

    lwd=0.5
    tf_in=rays2.TOFs(xs,recs,h1,h2,0,"UP",cR)     # Incident Raylegh surface wave
    ax2.plot(tf_in+t0,recs.xrecs,"--w",linewidth=lwd)
    nrefs=[0,1,3]   # the number of reflections
    for nref in nrefs:
        tf_in=rays2.TOFs(xs,recs,h1,h2,nref,"UP",cT)  # Direct waves
        ax2.plot(tf_in+t0,recs.xrecs,"--r",linewidth=lwd)


    #-----------------------------------
    fig2=plt.figure(figsize=[6,8])
    bx1=fig2.add_subplot(111) 
    bx1.grid(True)

    lift=0
    for nsc in nrefs:
        tf_sc=rays2.TOFs(xtip,recs,h1,h2,nsc,"UP",cT) #  TOF of scattered waves
        for nin in nrefs:
            tf_in=rays2.TOF(xs,xtip,h1,h2,nin,"UP",cT) # TOF of incident waves
            tf=tf_in+tf_sc  # Total TOF
            ax3.plot(tf+t0,recs.xrecs,"w--",linewidth=lwd)
            lbl=str(nin)+"-"+str(nsc)
            ax3.text(tf[-1]+t0,recs.xrecs[-1],lbl,horizontalalignment="center",verticalalignment="bottom",color="k")

            tref=tf[0]+t0
            tf-=tf[0];
            ampm=wvm.delay_and_sum_ray(-tf,sgn=-1)
            bx1.plot(wvm.time,ampm+lift)
            bx1.text(wvm.time[0],lift,lbl)
            bx1.plot(tref,lift,"ok",markersize=3)
            lift+=du;
    """
    crs_slow=np.linspace(2.6,3.9,14)
    crs_fast=np.linspace(4.0,6.6,14)
    crss=[crs_slow,crs_fast]

    m=0;
    #du=0.1 # 70deg
    for crs in crss:
        lift=0
        ylabels=[]
        ax=axs[m]
        for cr in crs:
            #ampp=wvp.delay_and_sum_t(c=cr,sgn=1)
            #ampm=wvm.delay_and_sum_t(c=cr,sgn=-1)
            dly=-(wvm.ycod-wvm.y0)/cr
            ampm=wvm.delay_and_sum_ray(dly,sgn=-1)
            ax.plot(wvm.time,ampm+lift)
            lift+=du;
            ylabels.append("%2.1f"%cr)
        ax.set_xlim(wvm.time[0],wvm.time[-1])
        ax.set_title(dir_name,fontsize=fsz+2)
        ax.set_xlabel("time $t_0$ [$\mu$s]",fontsize=fsz)
        ax.set_yticks(np.arange(len(crs))*du)
        ax.set_yticklabels(ylabels)
        ax.set_ylabel("horizontal lwave speed $V$ [km/s]",fontsize=fsz)
        m+=1;
    """

    dir_out="./"
    fig1.savefig(dir_out+dir_name+"_bwv.png",bboxh_inches="tight")
    fig2.savefig(dir_out+dir_name+".png",bboxh_inches="tight")
    #fig3.savefig(dir_out+dir_name+"_fast.png",bboxh_inches="tight")
    """
    fig1.savefig(dir_out+dir_name+"_bwv.eps",bboxh_inches="tight")
    fig2.savefig(dir_out+dir_name+"_slow.eps",bboxh_inches="tight")
    fig3.savefig(dir_out+dir_name+"_fast.eps",bboxh_inches="tight")
    """
    plt.show()
