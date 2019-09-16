import numpy as np
import sys
from numba import njit
import time


def contraction(p,q,r,s):
    status = False
    # heck if all particles are pairs
    if (q-p) == 1 and (s-r) == 1:
        # heck if pairs are legal pairs
        if p%2 == 0 and q%2 == 1 and r%2 == 0 and s%2 == 1:
            status = True
    return status

    
def initialize_all(p,h,delta,g):
    V_abcd = np.zeros((p, p, p, p))
    V_ijab = np.zeros((p, p, h, h))
    V_ijkl = np.zeros((h, h, h, h))
    f_ab =   np.zeros((p, p))
    f_ij =   np.zeros((h, h))
    t =      np.zeros((p, p, h, h))
    t_new =  np.zeros((p, p, h, h))
    H_ijab = np.zeros((p, p, h, h))
    D_ijab = np.zeros((p, p, h, h))

    v = -0.5*g # Interaction expectation value

    # Two-body (eq. 22)
    for a in range(p):
        for b in range(p):
            for c in range(p):
                for d in range(p):
                    if contraction(a,b,c,d) == True: 
                        V_abcd[a,b,c,d] = v
                        V_abcd[a,b,d,c] = -v
                        V_abcd[b,a,c,d] = -v
                        V_abcd[b,a,d,c] = v
                    if contraction(a,b,c,d) == True and c < h and d < h:
                        V_ijab[a,b,c,d] = v
                        V_ijab[a,b,d,c] = -v
                        V_ijab[b,a,c,d] = -v
                        V_ijab[b,a,d,c] = v
                    if contraction(a,b,c,d) == True and a < h and b < h and c < h and d < h: 
                        V_ijkl[a,b,c,d] = v
                        V_ijkl[a,b,d,c] = -v
                        V_ijkl[b,a,c,d] = -v
                        V_ijkl[b,a,d,c] = v

    # One-body (eq. 22)
    for a in range(p):
        for b in range(p):
            if a == b:
                if a <h and b < h:
                    f_ij[a,b] = (a - a%2)/2
                    for i in range(h):
                        f_ij[a,b] += V_ijkl[a,i,b,i]
                f_ab[a,b] = ((h+a) - (a+h)%2)/2

    # Cluster amplitudes (eq. 14)
    for i in range(h):
        for j in range(h):
            for a in range(p):
                for b in range(p):
                    D_ijab[a,b,i,j] = f_ij[i,i] + f_ij[j,j] - f_ab[a,a] - f_ab[b,b]
                    t[a,b,i,j] = V_ijab[a,b,i,j]/D_ijab[a,b,i,j]

    return V_abcd,V_ijkl,V_ijab,f_ab,f_ij,D_ijab,t,t_new,H_ijab

    
@njit
def update_H(p,h,
             V_abcd,V_ijkl,V_ijab,f_ab,f_ij,
             t,
             X_ijkl,X_il,X_ad,X_bcjk):
    """
    Eq. (35)
    """
    #H_ijab = V_ijab.copy()
    H_ijab = np.zeros((p,p,h,h))
    
    # Eq last
    for i in range(h):
        for j in range(h):
            for a in range(p):
                for b in range(p):
                    H_ijab[a,b,i,j] += V_ijab[a,b,i,j] # Term 1
                    for c in range(p):
                        H_ijab[a,b,i,j] += f_ab[b,c]*t[a,c,i,j] - f_ab[a,c]*t[b,c,i,j] # Term 2
                        H_ijab[a,b,i,j] += 0.5*(X_ad[a,c]*t[c,b,i,j] - X_ad[b,c]*t[c,a,i,j]) # Term 8 (c=d)
                        for d in range(p):
                            H_ijab[a,b,i,j] += 0.5*V_abcd[a,b,c,d]*t[c,d,i,j] # Term 4

                        #for k in range(h):
                            if d < h:
                                H_ijab[a,b,i,j] += 0.5*(t[a,c,i,d]*X_bcjk[b,c,j,d] # Term 6 (d=k)
                                                      - t[a,c,j,d]*X_bcjk[b,c,i,d]
                                                      - t[b,c,i,d]*X_bcjk[a,c,j,d]
                                                      + t[b,c,j,d]*X_bcjk[a,c,i,d])

                    for k in range(h):
                        H_ijab[a,b,i,j] -= (f_ij[k,j]*t[a,b,i,k] - f_ij[k,i]*t[a,b,j,k]) # Term 3
                        H_ijab[a,b,i,j] += 0.5*(X_il[i,k]*t[a,b,k,j] - X_il[j,k]*t[a,b,k,i]) # Term 7 k = l
                        for l in range(h):
                            H_ijab[a,b,i,j] += 0.5*V_ijkl[k,l,i,j]*t[a,b,k,l] # Term 5
                            H_ijab[a,b,i,j] += 0.25*t[a,b,k,l]*X_ijkl[k,l,i,j] # Term 9
                            
    return H_ijab

@njit
def update_X(p,h,V_ijab,t):
    X_ijkl = np.zeros((h,h,h,h))
    X_il = np.zeros((h,h))
    X_ad = np.zeros((p,p))
    X_bcjk = np.zeros((p,p,h,h))
    # hange a,b,i and j for n and m where it fits
    for c in range(p):
        for d in range(p):
            for k in range(h):
                for l in range(h):
                    for n in range(p):
                        for m in range(h):
                            if n < h:
                                X_ijkl[k,l,n,m] += V_ijab[k,l,c,d]*t[c,d,n,m] # n = i, m = j
                            X_bcjk[n,c,m,k] += V_ijab[k,l,c,d]*t[d,n,l,m] # n = b, m = j
                        if n < h:
                            X_il[n,l] += V_ijab[k,l,c,d]*t[c,d,n,k] # n = i
                        X_ad[n,d] += V_ijab[k,l,c,d]*t[n,c,k,l] # n = a
    return X_ijkl,X_il,X_ad,X_bcjk


@njit    
def update(p,h,
           V_abcd,V_ijkl,V_ijab,f_ab,f_ij,
           D_ijab,t,t_new,
           alpha):
    x1,x2,x3,x4 = update_X(p,h,V_ijab,t)
    H_ijab = update_H(p,h,V_abcd,V_ijkl,V_ijab,f_ab,f_ij,t,x1,x2,x3,x4)
    for i in range(h):
        for j in range(h):
            for a in range(p):
                for b in range(p):
                    t_new[a,b,i,j] = t[a,b,i,j] + H_ijab[a,b,i,j]/D_ijab[a,b,i,j] 

    t = alpha*t_new + (1-alpha)*t
    
    return t,t_new,H_ijab
    
@njit
def energy(p,h,V_ijab,t):
    """ 
    Eq. (34) 
    """
    E = 0
    for i in range(h):
        for j in range(h):
            for a in range(p):
                for b in range(p):
                    E += V_ijab[a,b,i,j]*t[a,b,i,j]
    return 0.25*E

def reference_energy(p,h,delta,g,V_ijkl):
    g = -0.5*g
    E_r = 0
    for i in range(h):
        E_r += delta*(i-(i%2))/2
        for j in range(h):
            E_r += 0.5*V_ijkl[i,j,i,j]
    return E_r
    
def run(p,h,delta,g,
        lim=1e-09, alpha=0.7,max_it=500):
    V_abcd,V_ijkl,V_ijab,f_ab,f_ij,D_ijab,t,t_new,H_ijab = initialize_all(p,h,delta,g)
    E_ref = reference_energy(p,h,delta,g,V_ijkl)
    energy1 = 0
    energy2 = 1000
    for i in range(max_it):
        t,t_new,H_ijab = update(p,h,V_abcd,V_ijkl,V_ijab,f_ab,f_ij,D_ijab,t,t_new,alpha)
        energy1 = energy(p,h,V_ijab,t)
        #print('It {}: {}, err={}'.format(i, energy2, np.abs(energy1-energy2)))
        #print(energy1,np.sum(np.abs(H_ijab)))
        if np.abs(energy1 - energy2) < lim:
            break
        energy2 = energy1
    
    return energy1,E_ref,i

if __name__ == '__main__':
    p_n = int(sys.argv[1])
    l_n = int(sys.argv[2])
    delta = float(sys.argv[3])
    g = float(sys.argv[4])
    alpha = float(sys.argv[5])

    h_n = l_n - p_n

    t1 = time.clock()
    E,E_ref,conv = run(p_n,h_n,delta,g)
    t2 = time.clock()
    print('-> Energy: ',E+E_ref)
    print('-> Time: ',t2-t1)
