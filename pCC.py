from w_inter import *
import numpy as np
import pandas as pd
import time


#pairs = [1,2,4]
pairs = [2]
#gs = [-1,-0.5,0,0.5,1]
gs = [-1,1]
delta = 1
E = []
basis = []
time_ = []
converge = []

alphas = [0.2,0.4,0.6,0.8]

for k,g in enumerate(gs):
    for a,alph in enumerate(alphas):
        df = pd.DataFrame()
        for n_pairs in pairs:
            p_n = 2*n_pairs
            df_dict = {}
            for j,n_basis in enumerate([i*n_pairs for i in range(2,50,2)]):
                print('{}/{} - {}/{} - {}/{}'.format(j+1,len(range(2,50,2)),a+1,len(alphas),k+1,len(gs)))
                t0 = time.process_time()
                Eig,E_ref,conv = run(n_basis,p_n,delta,g,alpha=alph)
                t1 = time.process_time()
                time_.append(t1 - t0)
                converge.append(conv)
                basis.append(n_basis)
                E.append(Eig+E_ref)
                df_dict = {'E':Eig+E_ref,'number_of_basis_states':n_basis,'Time':(t1 - t0),'conv':conv}
                df = df.append(df_dict,ignore_index=True)

        df.to_csv('Data/ny/CC_g={}_a={}.csv'.format(g,alph))
