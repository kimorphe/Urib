#! /home/kazushi/anaconda3/bin/python

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import copy

from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar
from mpl_toolkits.axes_grid1 import ImageGrid 


class WV1D:
    def __init__(self,Ny,dy=1.0,y0=0):

        self.Ny=Ny;
        self.nfile=Ny;

        self.dy=dy;
        self.y0=y0;

        self.Wd=dy*(Ny-1);

        self.y1=y0+self.Wd;
        self.ycod=self.y0+dy*np.arange(Ny);
        self.tLim=[0,80];

        self.fftk_exist=False
        self.fftw_exist=False
        self.fftkw_exist=False

    def load_1col(self,dirc="./",head="scope_",tail=".csv",nums=[]):
        if len(nums)==0:
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

        time=np.array(range(Nt))*dt+t1;
        time=np.array(time)*1.e06
        
        Amp=Amp*2.0;    # vol --> m/s
        Amax=np.max(np.abs(Amp))
        Amp=np.reshape(Amp,[self.Ny,Nt])/Amax;

        self.Nt=Nt;
        self.t1=time[0];
        self.t2=time[-1];
        self.dt=dt*1.e06;
        self.time=time
        self.amp=Amp

    def load(self,dirc="./",head="scope_",tail=".csv",nums=[]):
        if len(nums)==0:
            nums=np.arange(self.nfile);
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

        time=np.array(time)*1.e06
        Amp=Amp*2.0;    # vol. --> m/s
        Amax=np.max(np.abs(Amp))
        Amp=np.reshape(Amp,[self.Ny,Nt])/Amax;

        self.t1=time[0];
        self.t2=time[-1];
        self.dt=dt*1.e06;
        self.Nt=Nt;
        self.time=time
        self.amp=Amp
    def fwrite(self,fname_out):
        fp=open(fname_out,"w")
        fp.write("# Ny\n")
        fp.write(str(self.Ny)+"\n")
        fp.write("# t1, dt, Nt\n")
        fp.write(str(self.t1)+", "+str(self.dt)+", "+str(self.Nt)+"\n")
        fp.write("# amp\n")
        for k in range(self.Ny):
            for l in range(self.Nt):
                fp.write(str(self.amp[k,l])+"\n")
        fp.close()

    def show_xt(self,ax,cmp="jet",Vmin=-0.2,Vmax=0.2,fsz=12):
        ext=[self.time[0],self.time[-1],self.y0,self.y1]
        im=ax.imshow(self.amp,vmin=Vmin,vmax=Vmax,extent=ext,origin="lower",interpolation="bilinear",cmap=cmp)
        #ax.set_xlabel("time [$\mu$s]",fontsize=fsz)
        ax.set_ylabel("y [mm]",fontsize=fsz)
        ax_div=make_axes_locatable(ax)
        cax=ax_div.append_axes("right",size="5%",pad="2%")
        cbar=colorbar(im,cax=cax)

    def show_xw(self,ax,cmp="jet",Vmin=0.0,Vmax=0.8):
        fsz=12
        if not self.fftw_exist:
            self.FFTw()
        ext=[self.freq[0],self.freq[-1],self.y0,self.y1]
        Amax=np.max(np.abs(self.AMPw))
        ax.imshow(np.abs(self.AMPw)/Amax,vmin=Vmin,vmax=Vmax,extent=ext,origin="lower",interpolation="bilinear",cmap=cmp,aspect="auto")
        ax.set_xlim([0,10])
        ax.set_xlabel("frequency [MHz]",fontsize=fsz)
        ax.set_ylabel("y [mm]",fontsize=fsz)

    def show_kt(self,ax,cmp="jet",Vmin=0.0,Vmax=0.8):
        fsz=12
        if not self.fftk_exist:
            self.FFTk()
        ext=[self.time[0],self.time[-1],self.kx[0],self.kx[-1]]
        Amax=np.max(np.abs(self.AMPk))
        ax.imshow(np.abs(self.AMPk)/Amax,vmin=Vmin,vmax=Vmax,extent=ext,origin="lower",interpolation="bilinear",cmap=cmp,aspect="auto")
        ax.set_ylim([self.kx[0],self.kx[-1]*0.5])
        ax.set_xlabel("time [$\mu$s]",fontsize=fsz)
        ax.set_ylabel("wave number [/mm]",fontsize=fsz)

    def show_kw(self,ax,cmp="jet",Vmin=0.0,Vmax=0.8):
        fsz=12
        if not self.fftkw_exist:
            self.FFTkw()
        ext=[self.freq[0],self.freq[-1],self.kx[0],self.kx[-1]]
        Amax=np.max(np.abs(self.AMPkw))
        ax.imshow(np.abs(self.AMPkw)/Amax,vmin=Vmin,vmax=Vmax,extent=ext,origin="lower",interpolation="bilinear",cmap=cmp,aspect="auto")
        ax.set_xlim([0,10])
        ax.set_ylim([self.kx[0],self.kx[-1]])
        ax.set_xlabel("frequency [MHz]",fontsize=fsz)
        ax.set_ylabel("wave number [/mm]",fontsize=fsz)


    def FFTk(self):
        self.AMPk=np.fft.ifft(self.amp,axis=0);
        kmax=1./self.dy
        dk=kmax/self.Ny
        kx=np.arange(self.Ny)*dk;
        self.dk=dk
        self.kx=kx

        self.fftk_exist=True

    def FFTw(self):
        self.AMPw=np.fft.fft(self.amp,axis=1);
        fmax=1./self.dt;
        df=fmax/self.Nt;
        freq=np.arange(self.Nt)*df;
        self.df=df
        self.freq=freq
        self.fftw_exist=True

    def FFTkw(self):
        AMP=np.fft.fft(self.amp,axis=1);
        self.AMPkw=np.fft.ifft(AMP,axis=0);

        kmax=1./self.dy
        dk=kmax/self.Ny
        kx=np.arange(self.Ny)*dk;
        fmax=1./self.dt;
        df=fmax/self.Nt;
        freq=np.arange(self.Nt)*df;

        self.df=df
        self.freq=freq
        self.dk=dk
        self.kx=kx
        self.fftkw_exist=True

    def split(self):
        Axt=self.mirror();
        Axw=np.fft.fft(Axt,axis=1)
        Akw=np.fft.ifft(Axw,axis=0)

        Akwp=np.array(np.zeros(np.shape(Akw),dtype=complex));
        Akwm=np.array(np.zeros(np.shape(Akw),dtype=complex));

        Nt2=int(self.Nt/2)
        Ny=self.Ny
        Akwp[0:Ny,0:Nt2]=Akw[0:Ny,0:Nt2]*2;
        Akwm[Ny:-1,0:Nt2]=Akw[Ny:-1,0:Nt2]*2;

        ampp=np.fft.fft(Akwp,axis=0)
        ampp=np.real(np.fft.ifft(ampp,axis=1))

        ampm=np.fft.fft(Akwm,axis=0)
        ampm=np.real(np.fft.ifft(ampm,axis=1))

        ampp=ampp[0:Ny,:];
        ampm=ampm[0:Ny,:];

        return((ampp,ampm))
        
    def mirror(self):

        tmp=np.zeros(np.shape(self.amp))
        indx=np.arange(self.Ny,0,-1)-1
        indx=indx.astype(np.int)
        tmp=self.amp[indx,:]
        A=np.vstack([self.amp,tmp])

        print(np.shape(self.amp))
        print(np.shape(A))

        return(A)

    def delay_and_sum_t(self,c=3.1,sgn=1):

        amp=np.zeros(self.Nt);
        for k in range(self.Ny):
            dly=(self.ycod[k]-self.y0)/c;
            ndly=int(dly/self.dt)
            if ndly>=self.Nt-1:
                break;
            if sgn==1:
                amp[0:-ndly-1]=amp[0:-ndly-1]+self.amp[k,ndly:-1]
            if sgn==-1:
                amp[ndly:-1]=amp[ndly:-1]+self.amp[k,0:-1-ndly]

        amp/=self.Ny

        return(amp)

    def delay_and_sum_ray(self,dly,sgn=1):
        amp=np.zeros(self.Nt);
        for k in range(self.Ny):
            #dly=(self.ycod[k]-self.y0)/c;
            ndly=int(dly[k]/self.dt)
            if ndly>=self.Nt-1:
                break;
            if sgn==1:
                amp[0:-ndly-1]=amp[0:-ndly-1]+self.amp[k,ndly:-1]
            if sgn==-1:
                amp[ndly:-1]=amp[ndly:-1]+self.amp[k,0:-1-ndly]

        amp/=self.Ny

        return(amp)

    def delay_and_sum(self,c=3.1,sgn=1):

        if not self.fftkw_exist:
            self.FFTkw()

        omg=2.*self.freq*np.pi
        indx=self.freq/(c*self.dk);
        indx=indx.astype(np.int)

        AMP=np.array(np.zeros(self.Nt),dtype=complex)

        if sgn==1: 
            for i in range(int(self.Nt/2)):
                k=indx[i];
                if k<0:
                    continue;
                #if k>=int(self.Ny/2):
                if k>=int(self.Ny-1):
                    continue;
                AMP[i]=np.conj(self.AMPkw[k,i])

        if sgn==-1: 
            for i in range(int(self.Nt/2)):
                k=self.Ny-1-indx[i]
                if k<int(self.Ny/2):
                #if k<0:
                    continue;
                if k>=(self.Ny-1):
                    continue;
                AMP[i]=np.conj(self.AMPkw[k,i])

        amp=np.fft.ifft(np.array(AMP,dtype=complex))

        return(np.real(amp))

    def Lmax_normalize(self):

        Lmax=np.max(np.abs(self.amp))
        self.amp/=Lmax
        return(Lmax)


if __name__=="__main__":

    #------  Measurement Grid --------
    Nx=5; Ny=81;  
    dx=0.25; dy=0.25;
    nfile=Nx*Ny;
    x0=-(Nx-1)*0.5*dx; y0=0.0;
    dir_name="90_2_0209"
    dir_name="90_0_0209"
    dir_name="90_4_0209"
    dir_name="70_0_0211"
    dir_name="70_4_0211"
    dir_name="70_2_0212"
    one_col=False
    #---------------------------------


    wv=WV1D(Ny=Ny,dy=dy,y0=y0)

    nums=np.arange(Ny)*5+2;
    nums=nums.astype(np.int)
    nums=nums.tolist()

    if one_col:
        wv.load_1col(dirc=dir_name,nums=nums)
    else:
        wv.load(dirc=dir_name,nums=nums)


    wvp=copy.deepcopy(wv)
    wvm=copy.deepcopy(wv)
    (wvp.amp,wvm.amp)=wv.split()

    wvm.Lmax_normalize()
    wvp.Lmax_normalize()
    #-----------------------------------
    fig1=plt.figure(figsize=[8,8])
    ax1=fig1.add_subplot(311) 
    ax2=fig1.add_subplot(312) 
    ax3=fig1.add_subplot(313) 
    fsz=12
    cmap="bone"
    wv.show_xt(ax1,cmp=cmap,Vmin=-0.3,Vmax=0.3,fsz=fsz)
    wvp.show_xt(ax2,cmp=cmap,Vmin=-0.3,Vmax=0.3,fsz=fsz)
    wvm.show_xt(ax3,cmp=cmap,Vmin=-0.3,Vmax=0.3,fsz=fsz)
    ax3.set_xlabel("time [$\mu$s]",fontsize=fsz)
    ax1.set_title(dir_name,fontsize=fsz+2)
    ax1.text(18,2.0,"total field",color="w",horizontalalignment="left",fontsize=fsz)
    ax2.text(18,2.0,"left-going wave",color="w",horizontalalignment="left",fontsize=fsz)
    ax3.text(18,2.0,"right-going  wave",color="w",horizontalalignment="left",fontsize=fsz)

    #-----------------------------------
    fig2=plt.figure(figsize=[6,8])
    bx1=fig2.add_subplot(111) 
    bx1.grid(True)

    fig3=plt.figure(figsize=[6,8])
    bx2=fig3.add_subplot(111) 
    bx2.grid(True)

    axs=[bx1,bx2]

    crs_slow=np.linspace(2.6,3.9,14)
    crs_fast=np.linspace(4.0,6.6,14)
    crss=[crs_slow,crs_fast]

    m=0;
    for crs in crss:
        lift=0
        du=0.2
        ylabels=[]
        ax=axs[m]
        for cr in crs:
            ampp=wvp.delay_and_sum_t(c=cr,sgn=1)
            ampm=wvm.delay_and_sum_t(c=cr,sgn=-1)
            ax.plot(wvm.time,ampm+lift)
            lift+=du;
            ylabels.append("%2.1f"%cr)
        ax.set_xlim(wvm.time[0],wvm.time[-1])
        ax.set_title(dir_name,fontsize=fsz+2)
        ax.set_xlabel("time [$\mu$s]",fontsize=fsz)
        ax.set_yticks(np.arange(len(crs))*du)
        ax.set_yticklabels(ylabels)
        ax.set_ylabel("wave speed [km/s]",fontsize=fsz)
        m+=1;

    fig1.savefig(dir_name+"_bwv.png",bboxh_inches="tight")
    fig2.savefig(dir_name+"_slow.png",bboxh_inches="tight")
    fig3.savefig(dir_name+"_fast.png",bboxh_inches="tight")
    plt.show()
