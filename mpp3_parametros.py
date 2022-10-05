from cmath import pi
import matplotlib.pyplot as plt
import numpy as np
import math 


rho= 1.225 #kg/m3 density
eta= 1.81*10 ** -5 #kg/m*s dynamic viscosity
c= 343 #m/s
Z0=rho*c
f=np.arange(1, 10500, 100) #rango de frecuencias
w=2*pi*f
k=w/c

#panel 1
d1=0.15/1000 #diametro m
t1=0.74/1000 #espesor panel m
D1=1.75/100 #distancia cavidad m
p1=15/100 #ratio de perforacion %
r1=d1/2 #radio
e1=math.sqrt(p1) #epsilon
sigma1=32*eta/(p1*(d1 ** 2))
F1=(1-1.4092*e1+0.33818*e1 ** 3+0.06793*e1 ** 5-0.02287*e1 ** 6+0.03015*e1 ** 7-0.01641*e1 ** 8) ** -1
a1=1+0.85*d1/(t1*F1)
kd1=k*D1

#panel 2
d2=0.15/1000 
t2=1.5/1000 
D2=1.82/100 
p2=14.4/100 
r2=d2/2 
e2=math.sqrt(p2) 
sigma2=32*eta/(p2*d2 ** 2)
F2=(1-1.4092*e2+0.33818*e2 ** 3+0.06793*e2 ** 5-0.02287*e2 ** 6+0.03015*e2 ** 7-0.01641*e2 ** 8) ** -1
a2=1+0.85*d2/(t2*F2)
kd2=k*D2

# panel 3
d3=0.15/1000 
t3=1.5/1000 
D3=2.5/100 
p3=8/100 
r3=d3/2 
e3=math.sqrt(p3) 
sigma3=32*eta/(p3*d3 ** 2)
F3=(1-1.4092*e3+0.33818*e3 ** 3+0.06793*e3 ** 5-0.02287*e3 ** 6+0.03015*e3 ** 7-0.01641*e3 ** 8) ** -1
a3=1+0.85*d3/(t3*F3)
kd3=k*D3

Zm1=w*1j*rho*a1*t1/p1*(1+sigma1*p1/(1j*rho*w*a1)*(1+1j*4*w*rho*(a1 ** 2)*eta/((sigma1 ** 2)*(p1 ** 2)*(r1 ** 2))) ** 0.5)
Zm2=w*1j*rho*a2*t2/p2*(1+sigma2*p2/(1j*rho*w*a2)*(1+1j*4*w*rho*(a2 ** 2)*eta/((sigma2 ** 2)*(p2 ** 2)*(r2 ** 2))) ** 0.5)
Zm3=w*1j*rho*a3*t3/p3*(1+sigma3*p3/(1j*rho*w*a3)*(1+1j*4*w*rho*(a3 ** 2)*eta/((sigma3 ** 2)*(p3 ** 2)*(r3 ** 2))) ** 0.5)


Z3=Zm3-1j*Z0*np.cos(kd3)/np.sin(kd3)
Z2=Zm2+Z0*(Z3*np.cos(kd2)+1j*Z0*np.sin(kd2))/(Z0*np.cos(kd2)+1j*Z3*np.sin(kd2))
Z1=Zm1+Z0*(Z2*np.cos(kd1)+1j*Z0*np.sin(kd1))/(Z0*np.cos(kd1)+1j*Z2*np.sin(kd1))
R=(Z1-Z0)/(Z1+Z0)
alpha=1-(abs(R)) ** 2

#print('alpha',alpha)
plt.plot(f,alpha)
plt.grid()
plt.xlabel("frecuencia (Hz)")
plt.ylabel("\alpha_0")
plt.show()