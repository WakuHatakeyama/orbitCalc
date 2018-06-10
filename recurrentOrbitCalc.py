from math import*
import matplotlib.pyplot as plt
import numpy as np  

if __name__ == '__main__':
    N = 15
    i = 80 * pi/180
    e = 0
    mu = 3.98603 * 10**14   # m^3/sec^2
    a_E = 6378160   # m
    J2 = 1.082645*10**(-3)
    omega_E = 7.292*10**(-5)    # rad/sec
    iteration_num = 100000
    a = np.zeros((iteration_num,1))
    dT = np.zeros((iteration_num,1))

    # calculate initial condition
    T = 2*pi/N/omega_E
    a[0] = (mu*(T/(2*pi))**2)**(1/3)
    eps = 0.0001

    for j in range(len(a)-1):
        p = a[j]/a_E*(1-e**2)
        n0 = (mu/a[j]**3)**0.5

        # calculate n, Omega_dot, omega_dot, dT/da
        n = n0 * (1 + 3/2 * J2*(1-e**2)**0.5/p**2 * (1 - 3/2*sin(i)**2))
        Omega_dot = -3/2 * J2/p**2 * n * cos(i)
        omega_dot = 3/2 * J2/p**2 * n * (2-5/2*sin(i)**2)
        dT_da = -2*pi/(n+omega_dot)**2 \
                * (3/2*J2/p**2*n*(2-5/2*sin(i)**2)\
                -3/2*n/a[j]+n0*3/2*J2*(1-e**2)**0.5/p**2*(-2/a[j])*(1-3/2*sin(i)**2))
        

        T = 2*pi/(N*(omega_E-Omega_dot))
        Tn = 2*pi/(n+omega_dot)
        dT[j] = T-Tn
        # print(T,Tn)
        a[j+1] = a[j]+dT[j]/dT_da
        check = abs(dT[j]/T)+abs((a[j+1]-a[j])/a[j+1])
        if check < eps:
            print("height = "+ str((a[j]-a_E)/1000)+ " km")
            print("period = "+ str(T/60)+ " min")
            plt.plot(dT[0:j])
            break
    plt.show()