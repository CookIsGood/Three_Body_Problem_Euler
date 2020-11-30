import numpy as np
import math as math
from scipy.integrate import odeint
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from sympy import *

f1=1
M0= 1.9885*10**8
M1= M0/332940
M2= M1*317.8
M=M0+M1+M2
a= Symbol('a', real=True)
eq1 = (M0+M1)*a**5+(3*M0+2*M1)*a**4+(3*M0+M1)*a**3-(M1+3*M2)*a**2-(2*M1+3*M2)*a-(M2+M1)
solve_a=solve(eq1, a)
A=solve_a[0]

nu1=f1*((M0+M1)-M2/(A)**2+M2/(A+1)**2)
c=200000
po=c**2/nu1
e=0.5
p=po/(1+e)




phi=0 #задано
vp=0 #скорость изменения стороны треугольника
F=0 #угол собственного вращения
params = [M0, M1, M2,f1,c]
t = np.linspace(0, 10, 1001)
def f(y, t, params): #формула 14.26
    vp,p,F= y
    M0,M1,M2, f1,c = params
    return [c**2/p**3-nu1/p**2,
            vp,
            c/(p**2)
            ]




y0 = [vp,p,F]
[vp,p,F] = odeint(f, y0, t, args=(params,), full_output=False).T
r1=p
r2=p*(A+1)
a11=np.cos(F)#направляющие косинусы из Маркеева с 50 51
a12=-np.sin(F)
a21=np.sin(F)
a22=np.cos(F)
#подвижная система координат с 741 14.13
x1=a11*r1 #vx1=da11*p+vp*a11
y1=a21*r1
x2=a11*r2*np.cos(phi)+a12*r2*np.sin(phi)
y2=a21*r2*np.cos(phi)+a22*r2*np.sin(phi)
w=c/p**2
da11=-np.sin(F)*w# #направляющие косинусы из Маркеева для скоростей(производные а)
da12=-np.cos(F)*w
da21=np.cos(F)*w
da22=-np.sin(F)*w
vx1=da11*r1+vp*a11 #производные координат
vy1=da21*r1+vp*a21
vx2=(da11*r2+vp*a11)*np.cos(phi)+np.sin(phi)*(da12*r2+vp*a12)
vy2=(da21*r2+vp*a21)*np.cos(phi)+np.sin(phi)*(da22*r2+vp*a22)
print('x1',x1)
print('x2',x2)
print('y1',y1)
print('y2',y2)
print('vx1',vx1)
print('vx2',vx2)
print('vy1',vy1)
print('vy2',vy2)

#главная система координат переход с 735 ниже 14.7'
xx0=-M1/M*x1-M2/M*x2
xx1=(M2+M0)/M*x1-M2/M*x2
xx2=-M1/M*x1+(M1+M0)/M*x2
yy0=-M1/M*y1-M2/M*y2
yy1=(M2+M0)/M*y1-M2/M*y2
yy2=-M1/M*y1+(M1+M0)/M*y2
fig = plt.figure(facecolor='white', figsize=(6, 6))
plt.plot(xx1, yy1, linewidth=1, color='red')
plt.plot(xx2, yy2, linewidth=1, color='blue')
plt.plot(xx0, yy0, linewidth=1, color='black')
plt.grid(True)
plt.axis('scaled')
plt.show()