import numpy as np
import matplotlib.pyplot as plt


class TOFsc:
    def __init__(self,fname):
        fp=open(fname,"r")
        fp.readline()
        dat=fp.readline();
        dat=dat.strip().split(",")
        print(dat)
        ncrv=int(dat[0])
        Nx=int(dat[1])

        iin=[];
        isc=[];
        tof=[];
        xcod=[];
        for iref in range(ncrv):
            dat=fp.readline()
            dat=dat.strip().split(",")
            iin.append(int(dat[0]))
            isc.append(int(dat[1]))
            for i in range(Nx):
                dat=fp.readline().strip()
                dat=dat.split(" ")
                tof.append(float(dat[0]))
                xcod.append(float(dat[1]))
        fp.close()
        print(iin)
        print(isc)
        ncrv=len(iin);
        tof=np.array(tof);
        tof=np.reshape(tof,[ncrv,Nx])
        xcod=np.array(xcod);
        xcod=np.reshape(xcod,[ncrv,Nx])

        self.xcod=xcod;
        self.tof=tof;
        self.ncrv=ncrv;
        self.Nx=Nx;
        self.iin=iin;
        self.isc=isc;

    def show(self,ax):
        clrs=['b','g','r','c','m','y','k']
        nclrs=len(clrs)
        tb=["bottom","top"]
        x0=self.xcod[0,0]
        t0s=self.tof[:,0];
        indx=np.argsort(t0s);
        isum=0;
        for k in indx:
            lbl="("+str(self.iin[k])+","+str(self.isc[k])+")"
            t0=t0s[k]
            ax.plot(self.tof[k,:],self.xcod[k,:],color=clrs[isum%nclrs])
            #ax.text(t0,x0+(isum%2)*2.0,lbl,verticalalignment="bottom",horizontalalignment="center",color=clrs[isum%nclrs],rotation=70)
            ax.text(t0,x0,lbl,verticalalignment=tb[isum%2],horizontalalignment="center",color=clrs[isum%nclrs],rotation=75,fontsize=11)
            isum+=1;

class TOFin:
    def __init__(self,fname):
        fp=open(fname,"r")
        fp.readline()
        dat=fp.readline();
        dat=dat.strip().split(",")
        print(dat)
        ncrv=int(dat[0])
        Nx=int(dat[1])

        iin=[];
        tof=[];
        xcod=[];
        for iref in range(ncrv):
                dat=fp.readline()
                dat=dat.strip()
                iin.append(int(dat))
                for i in range(Nx):
                    dat=fp.readline().strip()
                    dat=dat.split(" ")
                    tof.append(float(dat[0]))
                    xcod.append(float(dat[1]))
        fp.close()
        print(iin)

        ncrv=len(iin);
        tof=np.array(tof);
        tof=np.reshape(tof,[ncrv,Nx])
        xcod=np.array(xcod);
        xcod=np.reshape(xcod,[ncrv,Nx])
        self.xcod=xcod;
        self.tof=tof;
        self.ncrv=ncrv;
        self.Nx=Nx;
        self.iin=iin;
    def show(self,ax):
        x0=self.xcod[0,0]
        t0s=self.tof[:,0];
        for k in range(self.ncrv):
            t0=t0s[k]
            ax.plot(self.tof[k,:],self.xcod[k,:],"k")
            ax.text(t0,x0,str(k),verticalalignment="bottom",horizontalalignment="center",color="k",fontsize=11)


if __name__=="__main__":
    fname1="tof_sc.out"
    fname2="tof_in.out"
    tf_sc=TOFsc(fname1)
    tf_in=TOFin(fname2)

    fig=plt.figure(figsize=(8,4))
    ax=fig.add_subplot(111)
    
    ax.grid(True)
    ax.set_xlim([15,50])
    ax.set_ylim((0,20))
    ax.set_xlabel("time [$\mu$sec]",fontsize=12)
    ax.set_ylabel("$x$ [mm]",fontsize=12)
    ax.tick_params(labelsize=12)

    tf_in.show(ax) 
    tf_sc.show(ax)
    fig.savefig("tof.png",bbox_inches="tight");
    plt.show()



