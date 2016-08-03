import math as math

def calcparameters(EPS, LAM):
    G = 5  # Number of beads per monomer
    # Ree = (1-LAM)
    Ree = 2  # End-to-end distance of monomer
    L0 = Ree*EPS*((-0.5+0.5*math.exp(-2*EPS*G)+EPS*G)**(-0.5))  # Interbead distance

    print 'EPS = %.2f, LAM = %.2f, Ree = %.4f, L0 = %.4f' %(EPS, LAM, Ree, L0)
    
for EPS in [0.01,0.1,1.00]:
    for LAM in [-0.75,-0.50,-0.25,0.00,0.25]:
        calcparameters(EPS,LAM)
