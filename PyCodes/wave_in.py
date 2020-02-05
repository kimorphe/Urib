import numpy as np
import matplotlib.pyplot as plt

fname="wave_sc4.dat"
fname="wave_in4.dat"

fp=open(fname,"r")

fp.readline()
Ny=int(fp.readline())
fp.readline()
dat=fp.readline()
dat=dat.strip()
dat=dat.split(",");

t0=float(dat[0])
dt=float(dat[1])
Nt=int(dat[2])
print("dt=",dt)

fp.readline()
A=[];
for row in fp:
    A.append(float(row))    

A=np.array(A)
print(np.shape(A))
A=np.reshape(A,[Ny,Nt])
#A=np.transpose(A)

t2=Nt*dt+t0;
fig=plt.figure()
ax=fig.add_subplot(111)
ax.imshow(A,aspect="auto",origin="lower",cmap="jet",extent=[t0,t2,0,(Ny-1)*0.25],vmin=-0.2,vmax=0.2,interpolation="bilinear")
#ax.imshow(A,aspect="auto",origin="lower",cmap="jet",extent=[t0,t2,0,(Ny-1)*0.25],interpolation="bicubic")

fp.close()
plt.show()
