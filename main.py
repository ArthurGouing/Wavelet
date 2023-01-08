import os
import numpy as np
import matplotlib.pyplot as plt
import pywt
from math import *

# Init Constante
PI = 4*atan(1.)
J = 15             # Level of the Dyadic decomposition
N = 6             # Nb vanishing moment ( Db{N/2} )
M = 2**J          # Number of point on curves
K = N + 2**J      # Size of the matricies

alpha = 9.5*PI    # Equation variable  d²u + alpha*u = 0

def Write_line(T, i, j, Omega):
    if j<0 or j>=K:
        return
    if j==i:
        T[i][j] = alpha**2
    T[i][j] += Omega.get(j-i, 0)
    #print(f"update T[{i}][{j}] with Omega({j-i})")

def Mw(x, N):
    if x<0. or x>N-1:
        return 0.
    name = "db"+str(int(N/2))
    w = pywt.Wavelet(name)
    (phi, psi, Y) = w.wavefun(level=10)
    m = N
    for i, y in enumerate(Y):
        if abs(x-y)<m:
            m=abs(x-y)
            index=i
    #print(len(Y), Y)
    #plt.plot(Y, phi)
    #plt.show()
    return phi[index]

def main():

    # Initialize
    print(f"On cherche à résoudre l'EDO : d²u + alpha²*u = 0")
    print(f"On utilise 2^{J} soit {M} points")
    print(f"On aura {K} wavelet coefficients")
    X = np.linspace(0, 1, M)
    Sol_E = np.cos(alpha*X)-cos(alpha)/sin(alpha)*np.sin(alpha*X)
    #Sol_E = np.cos(alpha*X)

    # Wavelet Galerkine Mehtode
    # Compute the connective coefficient
    O_4=  8.777142857143141e1
    O_3=  1.872457142857136e3
    O_2= -1.435550476190484e4
    O_1=  5.554956190476204e4
    O_0= -8.630857142857112e4
    O1 =  5.554956190476204e4
    O2 = -1.435550476190484e4
    O3 =  1.872457142857136e3
    O4 =  8.777142857143141e1
    Omega_7 = {
        -4 :  8.777142857143141e1,
        -3 :  1.872457142857136e3,
        -2 : -1.435550476190484e4,
        -1 :  5.554956190476204e4,
        0  : -8.630857142857112e4,
        1  :  5.554956190476204e4,
        2  : -1.435550476190484e4,
        3  :  1.872457142857136e3,
        4  :  8.777142857143141e1 }

    # Build T
    I = O_0 + alpha
    T = np.zeros((K, K))
    for j in range(K): # bcp d'appel Mw inutile, il faut init que ~N-1 valeur
        k = j-N+1
        print(j, k)
        T[0][j]       = Mw(-k, N)
        T[K-1][j] = Mw(2**J-k, 6)
    for i in range(1, K-1): # Il faudra enlevé la 1er et dernière ligne qui coressondent au condition aux limites -> range(1, K-1)
        for j in range(i-4, i+5):
            Write_line(T, i, j, Omega_7)
			

    Ttest = np.array([[0,  Mw(4., 6), Mw(3., 6), Mw(2., 6), Mw(1., 6),  0.,  0.],
                 [O1,  I, O_1, O_2, O_3, O_4,  0.],
                 [O2, O1,   I, O_1, O_2, O_3, O_4],
                 [O3, O2,  O1,   I, O_1, O_2, O_3],
                 [O4, O3,  O2,  O1,   I, O_1, O_2],
                 [0., O4,  O3,  O2,  O1,   I, O_1],
                 [0., 0., Mw(4., 6), Mw(3., 6), Mw(2., 6), Mw(1., 6), 0.]])

    # Build B, le vecteur nul
    B = np.zeros(K)
    B[0]   = 1.
    B[K-1] = 0.

    #print("Tdiff : \n", abs(Ttest-T))
    print("T : \n", T)
    print("T1 : \n", T[0][:10])
    print("TK-1 : \n", T[K-1][K-10:])
    print("B : \n", B)
    
    
    C = np.linalg.solve(T, B)
    #C = [-0.9972, -0.8776, 0.1279, 1.0870, 0.2479, -0.5059]
    C_ex = pywt.wavedec(Sol_E, 'db3', level=0)
    Ca_ex = C_ex[0]
    print("Il y en a : ", len(C))
    print("On trouve les coefficients :", C)
    #print("On devrais avoir :", len(Ca_ex)," coefficients")
    print("On devrais avoir : ", len(Ca_ex), "coefficients")
    print("On devrais avoir :", Ca_ex)
    
    # Recreate the signal from the coefficients
    Sol_WGM = pywt.upcoef("a", C, 'db3', level=1)
    Sol_E_approx = pywt.upcoef("a", Ca_ex, 'db3', level=1)
    #print("On trouve la solution :", Sol_WGM)


    # Plot
    f = plt.figure()
    axe = f.add_axes([0.1, 0.1, 0.8, 0.8])
    plt.title("Cas test 1 : N="+str(N)+" et J="+str(J))
    print(f"taille : \n  Sol_E : {len(Sol_E)} \n Sol_E_approx : {len(Sol_E_approx)} \n Sol_WGM : {len(Sol_WGM)}")

    axe.plot(X, Sol_E, label="Solution exact")

    X_approx = np.linspace(0., 1., len(Sol_E_approx))
    axe.plot(X_approx, Sol_E_approx, label="Solution approx")

    X_WGM = np.linspace(0, 1, len(Sol_WGM))
    axe.plot(X_WGM, Sol_WGM, label="Solution WGM")
    axe.set_xlim([0,1])
    axe.set_ylim([-2, 2])
    axe.legend()
    plt.show()
    


def Om(d1, d2, l):
    return

if __name__ == '__main__':
    main()
