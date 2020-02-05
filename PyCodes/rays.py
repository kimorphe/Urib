import numpy as np
import matplotlib.pyplot as plt
import sample 
import bscan


def mirror(ycod,y0):  # return reflection of ycod about y=y0
    return(2*y0-ycod)

# find the image point corresponding to the nref-th reflection waves
# (upward wave traveling toward the receiver)
def mirror_down(ycod,h_top,h_btm,nref):
    n=0
    yd=ycod;
    while n<np.abs(nref):
        yd=mirror(yd,h_btm)
        n+=1
        if n==np.abs(nref):
            break
        yd=mirror(yd,h_top)
        n+=1
    return(yd)

# find the image point corresponding to the nref-th reflection waves
# (donward wave traveling toward the receiver)
def mirror_up(ycod,h_top,h_btm,nref):
    n=0
    yd=ycod;
    while n<np.abs(nref):
        yd=mirror(yd,h_top)
        n+=1
        if n==np.abs(nref):
            break
        yd=mirror(yd,h_btm)
        n+=1
    return(yd)

# Get Time-of-Flight (TOF)
def TOF(xs,xr,h_top,h_btm,nref,UD,vel):

    X=xr[0]-xs[0];
    if UD=="UP":
        ycod=mirror_down(xr[1],h_top,h_btm,nref)
    if UD=="DOWN":
        ycod=mirror_up(xr[1],h_top,h_btm,nref)
    Y=ycod-xs[1];
    Lk=np.sqrt(X*X+Y*Y)

    return(Lk/vel)

# Receiver Class
class RECs:
    def setup(self,xr1,xr2,nr,yr):
        dxr=(xr2-xr1)/(nr-1)
        self.Nrec=nr;
        self.xrecs=np.arange(nr)*dxr+xr1
        self.yrecs=np.ones(nr)*yr;
        self.dxr=dxr

# TOFs gained for Receiver Class Data
def TOFs(xs,recs, h_top,h_btm,nref,UD,vel):
    xr=np.zeros(2)
    tf=[]
    for k in range(recs.Nrec):
        xr[0]=recs.xrecs[k]
        xr[1]=recs.yrecs[k]

        X=xr[0]-xs[0]
        if UD=="UP":
            ycod=mirror_down(xr[1],h_top,h_btm,nref)
        if UD=="DOWN":
            ycod=mirror_up(xr[1],h_top,h_btm,nref)

        Y=ycod-xs[1];
        Lk=np.sqrt(X*X+Y*Y)
        tf.append(Lk/vel)

    return(np.array(tf))


# Get travel path
def PATH(xs,xr,h_top,h_btm,nref,UD):

    X=xr[0]-xs[0];

    if UD=="UP":
        ycod=mirror_down(xr[1],h_top,h_btm,nref)
    if UD=="DOWN":
        ycod=mirror_up(xr[1],h_top,h_btm,nref)

    xk=[]
    xk.append(0)

    Y=ycod-xs[1];

    ht=h_top-h_btm
    hb=0.5*(h_top+h_btm)
    yref=np.arange(nref)
    yref=(yref%2)*2-1

    Xk=np.zeros(nref)
    if Y<0:
        tan_th=-X/Y;
        Xk=np.arange(nref)*tan_th*ht;
        Xk+=((xs[1]-h_btm)*tan_th);
        Yk=hb+yref*0.5*ht;
    elif Y>0:
        tan_th=X/Y;
        Xk=np.arange(nref)*tan_th*ht;
        Xk+=((h_top-xs[1])*tan_th);
        Yk=hb-yref*0.5*ht;
    else:
        tan_th=np.pi*0.5;
        Yk=hb-yref*0.5*ht
        #Xk=np.arange(nref)*tan_th*ht;
        #Xk+=((h_top-xs[1])*tan_th);
        Xk=np.arange(nref)*0.0
        Xk+=(h_top-xs[1])*0.0

    Xk+=xs[0];

    Xk=np.insert(Xk,0,xs[0])
    Yk=np.insert(Yk,0,xs[1])
    Xk=np.append(Xk,xr[0])
    Yk=np.append(Yk,xr[1])

    return([Xk,Yk])


#       Main Program
if __name__=="__main__":

    nums=range(3,404,5)
    Bscn=bscan.Bscan(nums)

    SMP=sample.SAMPLE(2)
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
    #--------- Travel Paths (example)------------
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.grid(True)
    ax.set_xlabel("y [mm]")
    ax.set_ylabel("x [mm]")
    Xmin=-10.; Xmax=70.
    ax.set_aspect("equal") # set aspect ration to unity

    xs=xtip
    xr=np.array([64.0,h1])      # source point

    nrefs=range(4);
    for nref in nrefs:  # the number of reflection varied
        [Xk,Yk]=PATH(xs,xr,h1,h2,nref,"UP") # get travel path data
        ax.plot(Xk,Yk)  # plot travel path

    SMP.show(ax)

    #---------- Travel Time Plot -----------------
    fig2=plt.figure()
    ax2=fig2.add_subplot(111)
    ax2.grid(True)
    ax2.set_xlabel("time [micro sec]")
    ax2.set_ylabel("receiver position x [mm]")

    xs=np.array([51.0,0.0])  # source point (incident point)
    SMP.set_src(xs[0],xs[1])

    recs=RECs() # use RECs class to handle mutiple receiver point
    Nrec=41;  # the number of receiver pointsｙ
    xr1=0.0;  # x-coordinate ( start )
    xr2=20.0;  # x-coordinate ( end )
    yr=h1;    # y-coordinate
    recs.setup(xr1,xr2,Nrec,h1) # setup reception point data
    SMP.set_rec(xr1,xr2,yr,Nrec)


    # Time-of-Flight (TOF) computation
    cax=Bscn.show(ax2,x1=xr2,dx=-0.25)

    tf_in=TOFs(xs,recs,h1,h2,0,"UP",cR)     # Incident Raylegh surface wave
    ax2.plot(tf_in+t0,recs.xrecs,"r")
    nrefs=[0,1,3,5]   # the number of reflections
    for nref in nrefs:
        tf_in=TOFs(xs,recs,h1,h2,nref,"UP",cT)  # Direct waves
        ax2.plot(tf_in+t0,recs.xrecs,"--b")

    for nsc in nrefs:
        tf_sc=TOFs(xtip,recs,h1,h2,nsc,"UP",cT) #  TOF of scattered waves
        for nin in nrefs:
            tf_in=TOF(xs,xtip,h1,h2,nin,"UP",cT) # TOF of incident waves
            tf=tf_in+tf_sc  # Total TOF
            ax2.plot(tf+t0,recs.xrecs,"g--")
            lbl=str(nin)+"-"+str(nsc)
            ax2.text(tf[-1]+t0,recs.xrecs[-1],lbl,horizontalalignment="center",verticalalignment="top",color="g")

    for nsc in nrefs:
        tf_sc=TOFs(xroot,recs,h1,h2,nsc,"UP",cT) #  TOF of scattered waves
        for nin in nrefs:
            tf_in=TOF(xs,xroot,h1,h2,nin,"UP",cT) # TOF of incident waves
            tf=tf_in+tf_sc  # Total TOF
            ax2.plot(tf+t0,recs.xrecs,"m--")
            lbl=str(nin)+"-"+str(nsc)
            ax2.text(tf[-1]+t0,recs.xrecs[-1],lbl,horizontalalignment="center",verticalalignment="bottom",color="m")

    plt.show()
