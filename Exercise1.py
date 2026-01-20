import numpy as np
import matplotlib.pyplot as plt
E = 4
m = 1.79 
def U(x):
        return 0

def VHO(x,omega):
        return (x**2)*m*(omega**2)/2

def EvolM(x,E):
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

print("trial run:", RK4(0, 1, [0, 1], 0.01, der, E))

plt.plot(RK4(0, 10, [0, 1], 0.001, der, 1)[0], RK4(0, 10, [0, 1], 0.001, der, 1)[1])
plt.show()

y = 1
e_i = 1

while abs(y) > 10**(-4):
	y = RK4(0, 10, [0, 1], 0.001, der, e_i)[1][-1]
	e_i += 0.001
#this is a silly and inefficient approach of course
print("outputted y is", y, "and corresponding energy is", e_i)	
#Finding E st psi(10) approx 0 




