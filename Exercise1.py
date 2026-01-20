import numpy as np

E = 1
m = 1 
def U(x):
        return 0

def VHO(x,omega):
        return (x**2)*m*(omega**2)/2

def EvolM(x):
	return np.array([[0,1],[-2*m*(U(x)-E),0]])

def der(x, y):
	return EvolM(x)@y

def RK4(x_i,x_f, y_i, h, der): #this f is the defivative function particularly
        Nits = int((x_f - x_i)/h)
        x_n = x_i
        y_n = np.array(y_i)
        for i in range(1, Nits + 1):
                k1 = der(x_n, y_n)
                k2 = der(x_n + h/2, y_n + h*k1/2)
                k3 = der(x_n + h/2, y_n + h*k2/2)
                k4 = der(x_n + h, y_n + h*k3)
                y_n = y_n + (h/6)*(k1 + 2*k2 + 2*k3 +k4)
                x_n += h
        return y_n

print("trial run:", RK4(0, 1, [0, 1], 0.01, der))

N = 1000 #number of steps
