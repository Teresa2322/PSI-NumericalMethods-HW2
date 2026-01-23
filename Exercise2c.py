import numpy as np
import matplotlib.pyplot as plt

m = 1
omega = 1

def U(x):
	return (x**2)*m*(omega**2)/2
def V(x, a, b):
	if a <= x <= b:
		return U(x)
	else:
		return np.inf

def EvolM(x,E,a,b): #chose to not make it a function of U(x)
	return np.array([[0,1],[2*m*(V(x,a,b)-E),0]])

def der(x, y,E,a,b):
	return EvolM(x, E,a,b)@y

def RK4(x_i,x_f, y_i, h, der,E):
	Nits = int((x_f - x_i)/h)
	x_n = x_i
	y_n = np.array(y_i)
	x_array = []
	psi_array = []
	for i in range(1, Nits + 1):
		k1 = der(x_n, y_n, E, x_i,x_f)
		k2 = der(x_n + h/2, y_n + h*k1/2, E, x_i,x_f)	
		k3 = der(x_n + h/2, y_n + h*k2/2, E, x_i,x_f)
		k4 = der(x_n + h, y_n + h*k3, E, x_i,x_f)
		y_n = y_n + (h/6)*(k1 + 2*k2 + 2*k3 +k4)
		psi_array.append(y_n[0])
		x_n += h
		x_array.append(x_n)
	return x_array, psi_array

#E st psi(10) approx 0 
#I know there is some energy between 0.4 and 0.5 that satisfies this

L = 5.944244384765625

def RK4fE(E):
        return RK4(-L, L, [0, 1], 0.001, der, E)[1][-1]

tol = 10**(-10)

def root_finder(E_i, E_f):
	x0 = E_i
	roots_arr = []
	while x0 < E_f:
		x0 += 0.01
		if RK4fE(x0)*RK4fE(x0 + 0.01) < 0:
			a = x0 
			b = x0 + 0.01
			while abs(a-b) > tol:
				c = (a + b)/2
				if  RK4fE(c) == 0:
					break
				elif RK4fE(a)*RK4fE(c) < 0:
					b = c
				elif RK4fE(b)*RK4fE(c) < 0:
					a = c
			print("root", c)
			roots_arr.append(c)
	return roots_arr

#print("trial", root_finder(0,1))
result = root_finder(0,45)

def SHO_Elev(n):
	return omega*(n+(1/2))

n_arr = np.linspace(0,len(result)-1,len(result))#np.linspace(0,len(result), len(result))
print("narr is", n_arr)
plt.plot(n_arr, result,'o-', label = 'numerical')
plt.plot(n_arr, SHO_Elev(n_arr),'o-', label = 'analytical')
plt.legend()
#plt.plot(RK4(-L,L, [0, 1], 0.001, der, 24.762701482103346)[0], RK4(-L,L, [0, 1], 0.001, der, 24.762701482103346)[1])
plt.xlabel("n")
plt.ylabel("E_n")
plt.show()

