from math import *

from scipy import optimize 
import pyfits #modulo para leer archivos fits
import matplotlib.pyplot as p #modulo para graficar
import numpy as n #este modulo es para trabajar con matrices como en matlab
import matplotlib as mp


G=6.67*1e-11   #  Cte de gravitacion universal  
M=(1.0/6.67)*1e11   #Masa mayor
m=1.0                  #Masa menor

class Planeta(object):

    '''Planeta es una clase que crea un objeto (el cuerpo que orbita) y
        que puede almacenar informacion sobre su estado actual
        (posicion (x,y) y velocidad (vx,vy)). Se accede al estado mediante
        objeto.y_actual. Como input debe ingresarse la condicion inicial
        (que corresponde al primer estado del objeto) y el valor alfa

        ejemplo
        x=Planeta(n.array([1,2,3,4]), 0.01)
        estado_actual=x.y_actual   (=[1,2,3,4])
        '''

        

    def __init__(self, condicion_inicial, alpha=0.0):

        '''metodo especial que se usa para inicializar las instancias
            de una clase. Define los inputs: condicion inicial y alfa
            (0 por defecto)'''

        self.y_actual = condicion_inicial
        self.alpha = alpha

    def ecuacion_de_movimiento(self, dx=0, dy=0, dvx=0, dvy=0):

        '''Implementa la ecuacion de movimiento como sistema de ecuaciones
            de primer orden. Retorna un vector con la función que caracteriza
            al sistema de ecuaciones'''

        x, y, vx, vy = self.y_actual
        x += dx
        y += dy
        vx += dvx
        vy += dvy
        fx=G*M*x*((-n.sqrt(x**2+y**2)+2*self.alpha)/n.power(x**2+y**2,2))
        fy=G*M*y*((-n.sqrt(x**2+y**2)+2*self.alpha)/n.power(x**2+y**2,2))
        v=[vx,vy,fx,fy]


        return n.array(v)

    def avanza_euler(self,dt):

        '''Toma el estado actual del objeto y lo avanza al siguiente, en un
            intervalo de tiempo dt y utilizando el metodo de euler explicito.
            No retorna nada, solo actualiza el estado del objeto'''

        yn = self.y_actual + dt*(self.ecuacion_de_movimiento())
        self.y_actual=yn
        pass

    def avanza_rk4(self,dt):

        ''' Lo mismo que el anterior, pero utilizando el metodo de
            Runge-Kutta de orden 4'''

        k1=dt*(self.ecuacion_de_movimiento())
        k2=dt*(self.ecuacion_de_movimiento(k1[0]/2, k1[1]/2, k1[2]/2, k1[3]/2))
        k3=dt*(self.ecuacion_de_movimiento(k2[0]/2, k2[1]/2, k2[2]/2, k2[3]/2))
        k4=dt*(self.ecuacion_de_movimiento(k3[0], k3[1], k3[2], k3[3]))
        yn = self.y_actual + ((1.0/6)*(k1+2*k2+2*k3+k4))
        self.y_actual=yn
        pass

    def avanza_verlet(self,dt):

        ''' Lo mismo que el anterior pero con el metodo de verlet'''

        ux=self.ecuacion_de_movimiento()[0]
        uy=self.ecuacion_de_movimiento()[1]
        ax=self.ecuacion_de_movimiento()[2]
        ay=self.ecuacion_de_movimiento()[3]

        xyn= self.y_actual[0:2] + dt*(self.ecuacion_de_movimiento()[0:2])+ (self.ecuacion_de_movimiento()[2:4])*0.5*dt**2
        vn=self.y_actual[2:4]+0.5*dt*(self.ecuacion_de_movimiento()[2:4]+self.ecuacion_de_movimiento(xyn[0]-self.y_actual[0],xyn[1]-self.y_actual[1],0,0)[2:4])
        self.y_actual=[xyn[0], xyn[1], vn[0], vn[1]]
        pass

    def energia_total(self):

        '''Retorna la energía del sistema asociada al estado actual'''

        x=self.y_actual[0]
        y=self.y_actual[1]
        vx=self.y_actual[2]
        vy=self.y_actual[3]

        E=(m/2)*(vx**2+vy**2) - G*M*m/n.sqrt(x**2+y**2) + self.alpha*G*M*m/(x**2+y**2)
        return E
    
    

condicion_inicial=[10.0,0.0,0.0,0.41]  #condicion inicial
ci=n.array(condicion_inicial)
Tierra=Planeta(ci)     #creacion del objeto Tierra

t_final=5600.0     #tiempo final de integracion'''
N=10**6             #numero de divisiones del intervalo de int'''
dt=t_final/N          #paso de integracion'''


x=n.zeros(N+1)       #creacion de los vectores posicion, velocidad y energia'''
y=n.zeros(N+1)
vx=n.zeros(N+1)
vy=n.zeros(N+1)
E=n.zeros(N)


x[0]=ci[0]     #declarando condiciones iniciales'''
y[0]=ci[1]
vx[0]=ci[2]
vy[0]=ci[3]

for i in range(N):

    '''recurrencia que va llenando los vectores anteriores con los datos de
        cada estado actual con el metodo deseado. En este caso, verlet'''
    E[i]=Tierra.energia_total()
    Tierra.avanza_verlet(dt)
    x[i+1]=Tierra.y_actual[0]
    y[i+1]=Tierra.y_actual[1]
    vx[i+1]=Tierra.y_actual[2]
    vy[i+1]=Tierra.y_actual[3]


    
t=n.linspace(0, t_final -dt, N)  #vector tiempo para graficar la energia'''

p.plot(t,E)
p.xlabel('Tiempo [seg]')
p.ylabel('Energia [J]')
p.grid
p.title('Energia en funcion del tiempo')
p.show()    

p.plot(x, y)
p.xlabel('x [m]')
p.ylabel('y [m]')
p.title('Orbita planetaria')
p.grid
p.show()




