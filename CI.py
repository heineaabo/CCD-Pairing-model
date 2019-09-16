from hamiltonian import *
import numpy as np
import pandas as pd
import time


#pairs = [1,2,4]
pairs = [2]
#gs = [-1,-0.5,0,0.5,1]
gs = [1]
delta = 1
E = []
basis = []
time_ = []
df = pd.DataFrame()

for n_pairs in pairs:
    df_dict = {}
    for g in gs:
        for j,n_basis in enumerate([i*n_pairs for i in range(2,30,2)]):
            print('{} / {}'.format(j+1,len(range(2,30,2))))
            #n_basis = int(8*n_pairs)
            t0 = time.process_time()
            H = hamiltonian(n_pairs,n_basis,delta,g)
            eigvals,eigvecs = np.linalg.eigh(H[0])
            t1 = time.process_time()
            time_.append(t1 - t0)

            eigval = eigvals[0]
            basis.append(n_basis)
            E.append(eigval)
            df_dict = {'E':eigval,'number_of_basis_states':n_basis,'Time':(t1 - t0)}
            df = df.append(df_dict,ignore_index=True)

df.to_csv('Data/FCITime_1.csv')



