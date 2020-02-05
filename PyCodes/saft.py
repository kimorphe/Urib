#! /home/kazushi/anaconda3/bin/python
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
import sample

d=4 # crack length [mm]
fname="saft4_2.out";

SMP=sample.SAMPLE(2)

if len(sys.argv)>1:
    fname=sys.argv[1]

fp=open(fname,"r");

fp.readline(); # computational domain
Xa=list(map(float,fp.readline().lstrip().split(" ")));

fp.readline(); # computational domain
Xb=list(map(float,fp.readline().lstrip().split(" ")));

fp.readline(); # computational domain
Ndiv=list(map(int,fp.readline().lstrip().split(" ")));

fp.readline();	# Imaging area
K=list(map(float,fp.readlines()));	# Imaging area
K=np.transpose(np.reshape(K,Ndiv))


fig=plt.figure();
ax=fig.add_subplot(111)
divider=make_axes_locatable(ax)
cax=divider.append_axes("right",size="5%",pad=0.1)
#im=ax.imshow(K,origin="lower",extent=[Xa[0],Xb[0],Xa[1],Xb[1]],interpolation="none",vmin=0,vmax=2,cmap="gray");
#im=ax.imshow(-K,origin="lower",extent=[Xa[0],Xb[0],Xa[1],Xb[1]],interpolation="bilinear",cmap="jet",vmin=0,vmax=2000)
indx=np.unravel_index(np.argmax(-K),K.shape)
Xa=np.array(Xa)
Xb=np.array(Xb)
Ndiv=np.array(Ndiv)
dx=(Xb-Xa)/(Ndiv-1)
xmax=Xa[0]+indx[1]*dx[0];
ymax=Xa[1]+indx[0]*dx[1];
K/=np.max(-K)
im=ax.imshow(-K,origin="lower",extent=[Xa[0],Xb[0],Xa[1],Xb[1]],interpolation="bilinear",cmap="jet",vmin=0)
#ax.plot(xmax,ymax,"xw",markersize=10)
plt.colorbar(im,cax=cax)
fp.close();
fig.savefig("saft4.png",bbox_inches="tight");
SMP.show2(ax,clr="y")
ax.set_xlim([-20,10])
ax.set_ylim([-15,10])
ax.grid(True)
ax.tick_params(labelsize=12)
ax.set_xlabel("x[mm]",fontsize=12)
ax.set_ylabel("y[mm]",fontsize=12)

fnout=fname.replace(".out",".png")
print("-->"+fnout)
fig.savefig(fnout,bbox_inches="tight")
plt.show();


