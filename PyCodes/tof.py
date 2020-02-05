import numpy as np
import matplotlib.pyplot as plt
print('input distance')
t=int(input())
x = np.arange(0, 100, 0.1)

ht=12.0     # plate thickness [mm]
t0=8.162    # delay time [micro sec]

cR=2.910    # surface wave velocity [km/s]
cT=3.036    # transverse wave velocity [km/s]

Ames=20.0   # measurement area [mm]


#1
y1=cR*(x-t0)
plt.plot(x,y1)

#2
nref=1   # the number of reflections 
m=nref+1
a=np.sqrt(((t/m)**2)+ht**2)
b=((m*a)/cT)+t0
c=np.sqrt((((t+Ames)/m)**2)+ht**2)
d=((m*c)/cT)+t0
plt.plot([b,d],[t,t+Ames])

#4
nref=3
m=nref+1
e=np.sqrt(((t/m)**2)+ht**2)
f=((m*e)/cT)+t0
g=np.sqrt((((t+Ames)/m)**2)+ht**2)
h=((m*g)/cT)+t0
plt.plot([f,h],[t,t+Ames])

#5
y2=-cR*x+((2*t)+89.75142)
plt.plot(x,y2)

#8
i=np.sqrt((((t+45)/3)**2)+ht**2)
j=np.sqrt((25**2)+ht**2)
k=(((3*i)+j)/cT)+t0
j2=np.sqrt((45**2)+ht**2)
k2=(((3*i)+j2)/cT)+t0
plt.plot([k2,k],[t,t+Ames])

#9
l=np.sqrt((((t+45)/2)**2)+ht**2)
m=np.sqrt(((25/2)**2)+ht**2)
n=(((2*l)+(2*m))/cT)+t0
m2=np.sqrt(((45/2)**2)+ht**2)
n2=(((2*l)+(2*m2))/cT)+t0
plt.plot([n2,n],[t,t+Ames])

#10
o=np.sqrt((((t+45)/4)**2)+ht**2)
p=np.sqrt(((25/2)**2)+ht**2)
q=(((4*o)+(2*p))/cT)+t0
p2=np.sqrt(((45/2)**2)+ht**2)
q2=(((4*o)+(2*p2))/cT)+t0
plt.plot([q2,q],[t,t+Ames])

#11
r=np.sqrt(((t+45)**2)+ht**2)
s=np.sqrt((25**2)+ht**2)
s2=np.sqrt((45**2)+ht**2)
u=((r+s)/cT)+t0
u2=((r+s2)/cT)+t0
plt.plot([u2,u],[t,t+Ames])

plt.xlabel("time(μs)")
plt.ｙlabel("y(mm)")
plt.xlim([15,65])
plt.ylim([t,t+Ames])
plt.show()
