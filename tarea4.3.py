from math import *

from scipy import optimize 
import pyfits #modulo para leer archivos fits
import matplotlib.pyplot as p #modulo para graficar
import numpy as n #este modulo es para trabajar con matrices como en matlab
import matplotlib as mp


G=6.67*1e-11
M=(1.0/6.67)*1e11
m=1.0

class Planeta(object):

    def __init__(self, condicion_inicial, alpha=0.0):

        self.y_actual = condicion_inicial
        self.alpha = alpha

    def ecuacion_de_movimiento(self, dx=0, dy=0, dvx=0, dvy=0):

        x, y, vx, vy = self.y_actual
        x += dx
        y += dy
        vx += dvx
        vy += dvy
        fx=G*M*x*((-n.sqrt(x**2+y**2)+2*self.alpha)/n.power(x**2+y**2,2))
        fy=G*M*y*((-n.sqrt(x**2+y**2)+2*self.alpha)/n.power(x**2+y**2,2))
        v=[vx,vy,fx,fy]


        return n.array(v)


    def avanza_verlet(self,dt):

        ux=self.ecuacion_de_movimiento()[0]
        uy=self.ecuacion_de_movimiento()[1]
        ax=self.ecuacion_de_movimiento()[2]
        ay=self.ecuacion_de_movimiento()[3]

        xyn= self.y_actual[0:2] + dt*(self.ecuacion_de_movimiento()[0:2])+ (self.ecuacion_de_movimiento()[2:4])*0.5*dt**2
        vn=self.y_actual[2:4]+0.5*dt*(self.ecuacion_de_movimiento()[2:4]+self.ecuacion_de_movimiento(xyn[0]-self.y_actual[0],xyn[1]-self.y_actual[1],0,0)[2:4])
        self.y_actual=[xyn[0], xyn[1], vn[0], vn[1]]
        pass

    def energia_total(self):

        x=self.y_actual[0]
        y=self.y_actual[1]
        vx=self.y_actual[2]
        vy=self.y_actual[3]

        E=(m/2)*(vx**2+vy**2) - G*M*m/n.sqrt(x**2+y**2) + self.alpha*G*M*m/(x**2+y**2)
        return E
    
    

condicion_inicial=[10.0,0.0,0.0,0.41]
ci=n.array(condicion_inicial)
Tierra=Planeta(ci,0.001054386896)  #incluye el valor de alfa

t_final=1110.0
N=25*10**5
dt=t_final/N


x=n.zeros(N+1)
y=n.zeros(N+1)
vx=n.zeros(N+1)
vy=n.zeros(N+1)
E=n.zeros(N)


x[0]=ci[0]
y[0]=ci[1]
vx[0]=ci[2]
vy[0]=ci[3]

for i in range(N):
    E[i]=Tierra.energia_total()
    Tierra.avanza_verlet(dt)
    x[i+1]=Tierra.y_actual[0]
    y[i+1]=Tierra.y_actual[1]
    vx[i+1]=Tierra.y_actual[2]
    vy[i+1]=Tierra.y_actual[3]

r2=x**2+y**2  #vectr del cuadrado de la distancia al foco


i_min=n.argmin(r2[2400000:])  #posicion del minimo (perihelio)

theta=atan(y[i_min+2400000]/x[i_min+2400000])  #angulo de precesion

print theta   #devuelve el valor del angulo de precesion 

