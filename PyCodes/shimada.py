import numpy as np
import matplotlib.pyplot as plt
print('input distance')
t=int(input())
x = np.arange(0, 100, 0.1)

#1
y1=2.91*(x-8.162)
plt.plot(x,y1)

#2
a=np.sqrt(((t/2)**2)+12**2)
b=((2*a)/3.036)+8.162
c=np.sqrt((((t+20)/2)**2)+12**2)
d=((2*c)/3.036)+8.162
plt.plot([b,d],[t,t+20])

#4
e=np.sqrt(((t/4)**2)+12**2)
f=((4*e)/3.036)+8.162
g=np.sqrt((((t+20)/4)**2)+12**2)
h=((4*g)/3.036)+8.162
plt.plot([f,h],[t,t+20])

#5
y2=-2.91*x+((2*t)+89.75142)
plt.plot(x,y2)

#8
i=np.sqrt((((t+45)/3)**2)+12**2)
j=np.sqrt((25**2)+12**2)
k=(((3*i)+j)/3.036)+8.162
j2=np.sqrt((45**2)+12**2)
k2=(((3*i)+j2)/3.036)+8.162
plt.plot([k2,k],[t,t+20])

#9
l=np.sqrt((((t+45)/2)**2)+12**2)
m=np.sqrt(((25/2)**2)+12**2)
n=(((2*l)+(2*m))/3.036)+8.162
m2=np.sqrt(((45/2)**2)+12**2)
n2=(((2*l)+(2*m2))/3.036)+8.162
plt.plot([n2,n],[t,t+20])

#10
o=np.sqrt((((t+45)/4)**2)+12**2)
p=np.sqrt(((25/2)**2)+12**2)
q=(((4*o)+(2*p))/3.036)+8.162
p2=np.sqrt(((45/2)**2)+12**2)
q2=(((4*o)+(2*p2))/3.036)+8.162
plt.plot([q2,q],[t,t+20])

#11
r=np.sqrt(((t+45)**2)+12**2)
s=np.sqrt((25**2)+12**2)
s2=np.sqrt((45**2)+12**2)
u=((r+s)/3.036)+8.162
u2=((r+s2)/3.036)+8.162
plt.plot([u2,u],[t,t+20])

plt.xlabel("time(μs)")
plt.ｙlabel("y(mm)")
plt.xlim([15,65])
plt.ylim([t,t+20])
plt.show()
