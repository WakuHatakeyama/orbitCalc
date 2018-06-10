from math import*
import matplotlib.pyplot as plt
import numpy as np  

if __name__ == '__main__':
    N = 15.5
    i = 97 * pi/180
    e = 0
    mu = 3.98603 * 10**14   # m^3/sec^2
    a_E = 6378160   # m
    J2 = 1.082645*10**(-3)
    omega_E = 7.292*10**(-5)    # rad/sec
    a = np.zeros((15,1))
    dT = np.zeros((15,1))

    # calculate initial condition
    T = 2*pi/N/omega_E
    a[0] = (mu*(T/(2*pi))**2)**(1/3)

    # print(T)
    # print(a)
    for i in range(len(a)-1):
        p = a[i]/a_E*(1-e**2)
        n0 = (mu/a[i]**3)**0.5
        print(n0)
        # print(p)

        # calculate n, Omega_dot, omega_dot, dT/da
        n = n0 * (1 + 3/2 * J2*(1-e**2)**0.5/p**2 * (1 - 3/2*sin(i)**2))
        Omega_dot = -3/2 * J2/p**2 * n * cos(i)
        omega_dot = 3/2 * J2/p**2 * n * (2-5/2*sin(i)**2)
        # print(n)
        # print(Omega_dot)
        # print(omega_dot)
        dT_da = -2*pi/(n+omega_dot)**2 \
                * (1-e**2)/a_E \
                * (-3*J2/p**3*(2-5/2*sin(i)**2)*n\
                +(1+3/2*J2/p**2*(2-5/2*sin(i)**2))
                *(-3*J2*(1-e**2)/p**3*(1-3/2*sin(i)**2)))
        print (1/dT_da)
        # print(a)
        

        T = 2*pi/(N*(omega_E-Omega_dot))
        Tn = 2*pi/(n+omega_dot)
        dT[i] = T-Tn
        print(T,Tn,a[i],dT[i])
        a[i+1] = a[i]+dT[i]/dT_da
    plt.plot((a-6380000)/1000)
    plt.show()
