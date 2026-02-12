
import numpy as np
import scipy as sp
import matplotlib.pyplot as py

def f(x):
    return (np.sin(np.sqrt(100*x)))**2

def g(x):
    return (x**2)/np.sqrt(2-x)

def trapezoid(f,a,b,N):
    h = (b-a)/N #Interval size
    mysum = 0
    
    for i in range(1,N): #should go from 1 to N-1
        mysum = mysum + float(f(a+(i)*h)*h) 
    return mysum +(h/2)*(f(a) + f(b))

def TrapezoidTable(g,a,b):
    TrueVal = trapezoid(f,a,b,5000000) 

    print('\n',f"{'Intervals':<12} {'Approx Value':>16} {'Error' :>18}")
    for i in [2**n for n in range(14)]:
        val = trapezoid(g,a,b,i)
        print(f"{i:10d} {val:19.10f} {np.abs(TrueVal - val):18.10f}")
    print()
    return

#Legendre Polynomial plotting
def P(n,x):
    return sp.special.legendre(n)(x)

def PlotStuff():
    x = np.linspace(-1,1,100)
    fig, ax = py.subplots(4,4)

    for i in range(4):
        for j in range(4):
            ax[i,j].plot(x,P(i+1,x))
            ax[i,j].plot(x,P(j+1,x))
            ax[i,j].plot(x,P(i+1,x)*P(j+1,x))
    
    py.savefig('legendre')
    py.show()
    py.pause(5)
    return

#Gaussian Quadrature stuff

def F(g,a,b,u): #This is the transformed function that is to be integrated from -1 to 1
    return g(((b-a)/2)*u + (a+b)/2)*(b-a)/2

def GaussQuad(g,a,b,N):
    roots, weights = sp.special.roots_legendre(N)
    return np.sum([weights[i]*g(((b-a)/2)*roots[i] + (a+b)/2)*(b-a)/2 for i in range(N)])

def FindNumber(f,a,b,Alg,TrueVal,err):
    TrueVal = np.sqrt(8192)/15
    i = 2

    while np.abs(Alg(f,a,b,i) - TrueVal) > err:
        i *= 2

    upper = i
    lower = i/2
    guess = np.ceil((upper + lower)/2)

   #while np.abs(Alg(f,a,b,guess) - TrueVal - err) > 0.000000000001:
       #if 

    return
    
print(trapezoid(f,0,2,8192))
print(GaussQuad(f,0,2,16))
