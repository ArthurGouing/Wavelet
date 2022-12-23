import os
import numpy as np
import matplotlib.pyplot as plt
import pywt
from math import *

PI = 4*atan(1.)

def main():
    print("début du code", PI)

    # Initialize
    j = 9
    N = 2**j
    print(f"On utilise {N} points")
    alpha = 9.5*PI
    X = np.linspace(0, 1, N)
    Sol_E = np.cos(alpha*X)-cos(alpha)/sin(alpha)*np.sin(alpha*X)
    Sol_E2 = np.cos(10*PI*X)

    # Compute with galerkine method
    # Build T
    """
    T = 

    # Build B, le vecteur null
    B = np.zeros(N)
    
    
    C = np.linalg.solve(T, B)
    # Recreate the signal from the coefficients
    Sol_WGM = pywt.wavrec(C, 'db6')
    """


    # Plot
    f = plt.figure()
    axe = f.add_axes([0.1, 0.1, 0.8, 0.8])
    plt.title("Cas test 1 : N=6 et j=9")
    axe.plot(X, Sol_E2, label="Solution exact")
    axe.set_xlim([0,1])
    axe.set_ylim([-2, 2])
    axe.legend()
    plt.show()
    

if __name__ == '__main__':
    main()
