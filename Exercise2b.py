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
'''
E_arr = np.linspace(0,3*omega,100)
psi_barr = []
for e in E_arr:
	k = RK4fE(e)
	print("psi b:", k)
	psi_barr.append(k)


fig, (ax1,ax2) = plt.subplots(1,2,figsize = (10,4))

ax1.plot(E_arr , psi_barr, 'tab:blue')
ax1.set_xlim(0,3*omega)
ax1.set_title("Full Range")
ax1.set_xlabel("E")
ax1.set_ylabel("psi(b)")

ax2.plot(E_arr , psi_barr, 'tab:orange')
ax2.plot(E_arr, np.zeros(len(E_arr)), color = 'black')
ax2.set_title("Roots Region")
ax2.set_ylim(-1876987312.709363, 1876987312.709363)
ax2.set_xlim(0,3*omega)
ax2.set_xlabel("E")
ax2.set_ylabel("psi(b)")

fig.suptitle("psi(b) vs E [omega = 1]")
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()
'''
tol = 10**(-10)

def root_finder(E_i, E_f):
	x0 = E_i
	roots_arr = []
	while x0 < E_f:
		x0 += 0.01
		if RK4fE(x0)*RK4fE(x0 + 0.01) < 0:
			a = x0 
			b = x0 + 0.01
			c = 1
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

n_arr = [1,2,3,4,5,6]


plt.plot(n_arr,  root_finder(0, 3*omega), 'bo-')
plt.title("Infinite Square well Energy Levels")
plt.xlabel("nth Energy Level")
plt.ylabel("Energy")
plt.show()
