import numpy as np
import matplotlib.pyplot as plt

class SAMPLE:
    def __init__(self,num=0):
        if num==0: #0mm
            PT=[\
            [0,0],[0,12],[12,12],[12,13],[11,13],\
            [6.9,36],[14.9,36],[20,12],[180,12],[180,0]]
            nroot=2
            ntip=2
        elif num==2: #2mm
            PT=[[0,0],[0,12],[11.8,12],[12.3,10.1],\
            [12.5,10.1],[12,12],[12,13],[11,13],[6.9,36],\
            [14.9,36],[20,12],[180,12],[180,0]]
            nroot=5
            ntip=3
        elif num==4:
            PT=[[0,0],[0,12],[11.8,12],[12.8,8.1],[13.0,8.1],\
            [12,12],[12,13],[11,13],[6.9,36],[14.9,36],\
            [20,12],[180,12],[180,0]]
            nroot=5
            ntip=3
        else:
            print("No such sample of No.",num)
            exit()
        x0=20.0
        y0=12.0

        PT.append(PT[0])
        PT=np.array(PT)
        PT[:,0]-=x0
        PT[:,1]-=y0
        self.h1=PT[1,1]
        self.h2=PT[0,1]
        self.PT=PT
        self.nroot=nroot
        self.ntip=ntip

    def show(self,ax):
        ax.plot(self.PT[:,0],self.PT[:,1])
        try:
            ax.plot(self.xrs,self.yrs,".k")
        except:
            print("Reciver array has not been defined!");
        try:
            ax.plot(self.xs,self.ys,"bo")
        except:
            print("Source point has not been defined!");
    def show2(self,ax,clr="k"):
        ax.plot(self.PT[:,0],self.PT[:,1],color=clr,linewidth=1)

    def set_rec(self,xr1,xr2,yr,nrec):
        xrs=np.linspace(xr1,xr2,nrec);
        dxr=0.0
        if nrec >1:
            dxr=xrs[1]-xrs[0];
        yrs=np.ones(nrec)*yr;

        self.xrs=xrs;
        self.yrs=yrs;
        self.nrec=nrec;
        self.dxr=dxr;
    def set_src(self,xs,ys):
        self.xs=xs;
        self.ys=ys;
    
if __name__=="__main__":

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.grid(True)
    ax.set_aspect(1.0)
    ax.set_xlabel("x[mm]")
    ax.set_ylabel("y[mm]")

    SMP=SAMPLE(num=4)
    SMP.set_rec(5.,20.,0.0,41)
    SMP.set_src(56,0.0)

    SMP.show(ax)
    
    plt.show()


    def set_rec(self,xr1,xr2,yr,nrec):
        xrs=np.linspace(xr1,xr2,nrec);
        dxr=0.0
        if nrec >1:
            dxr=xrs[1]-xrs[0];
        yrs=np.ones(nrec)*yr;

        self.xrs=xrs;
        self.yrs=yrs;
        self.nrec=nrec;
        self.dxr=dxr;
    def set_src(self,xs,ys):
        self.xs=xs;
        self.ys=ys;
    
if __name__=="__main__":

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.grid(True)
    ax.set_aspect(1.0)
    ax.set_xlabel("x[mm]")
    ax.set_ylabel("y[mm]")

    SMP=SAMPLE(num=4)
    SMP.set_rec(5.,20.,0.0,41)
    SMP.set_src(56,0.0)

    SMP.show(ax)
    
    plt.show()

