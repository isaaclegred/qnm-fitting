import numpy as np
import mylib
a = np.array([2.0,3.0,5.0,3.0])
c = np.array([5.0,2.0])
guess = np.array([2.2,2.8,4.9,2.8])
i = 100
j = 100
N = 100
MA = 4
mfit = MA
P = MA/2
pfit = P
SPREAD  = .01
x = np.linspace(0, 10,  100)
y = c[0]*np.exp(-((x - a[0]) / a[1])**2) + c[1]*np.exp(-((x - a[2]) / a[3])**2)
y *= (1 + SPREAD * np.random.randn())
marq  = mylib.get_marquardt_gauss(x, y,  np.ones((N)) * SPREAD, guess, 2)
print(marq.get_values())
print(marq.get_nl_params())
print(marq.get_lin_params())
#print(marq.get_times())
marq.fit()
