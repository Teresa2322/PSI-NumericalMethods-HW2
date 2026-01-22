import numpy as np
import matplotlib.pyplot as plt

m = 1

def U(x):
        return 0 #(x**2)*m*(omega**2)/2 for Harmonic

def VHO(x,omega):
        return (x**2)*m*(omega**2)/2

def EvolM(x,E): #chose to not make it a function of U(x)
	return np.array([[0,1],[2*m*(U(x)-E),0]])

def der(x, y,E):
	return EvolM(x, E)@y

def RK4(x_i,x_f, y_i, h, der,E):
	Nits = int((x_f - x_i)/h)
	x_n = x_i
	y_n = np.array(y_i)
	x_array = []
	psi_array = []
	for i in range(1, Nits + 1):
		k1 = der(x_n, y_n, E)
		k2 = der(x_n + h/2, y_n + h*k1/2, E)	
		k3 = der(x_n + h/2, y_n + h*k2/2, E)
		k4 = der(x_n + h, y_n + h*k3, E)
		y_n = y_n + (h/6)*(k1 + 2*k2 + 2*k3 +k4)
		psi_array.append(y_n[0])
		x_n += h
		x_array.append(x_n)
	return x_array, psi_array

#E st psi(10) approx 0 
#I know there is some energy between 0.4 and 0.5 that satisfies this

tol = 10**(-5)
a = 0.4
b = 0.5

def RK4fE(E):
	return RK4(0, 10, [0, 1], 0.001, der, E)[1][-1]

while (b - a)/2 > tol:
	c = (a + b)/2
	if  RK4fE(c) == 0:
		break
	elif   RK4fE(a)* RK4fE(c) < 0:
		b = c
	elif   RK4fE(b)* RK4fE(c) < 0:
          	a = c
	print("Root found at:",c," with psi at that value being", RK4fE(c))

plt.plot(RK4(0, 10, [0, 1], 0.001, der, c)[0], RK4(0, 10, [0, 1], 0.001, der, c)[1])
plt.xlabel("x")
plt.ylabel("psi(x)")
plt.savefig('Exercise1Ca.png', dpi=300)
plt.close()

plt.show()

