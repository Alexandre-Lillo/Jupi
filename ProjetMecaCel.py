### FICHIER DE TRAVAIL 

#%% Constantes et Packages 

import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

m_jup = 1.8986e27 # masse de Jupiter en kg
m_sol = 1.9891e30 # masse du Soleil en kg
m_sat = 568.46e24 # masse de Saturne en kg
G = 1.4882142e-34 # cste gravitationnelle en AU3*kg-1*j-2
k = 30 # pas de temps en jours
N = 20 #int((365*5000 + 1250 + k)/k)


# Conditions initiales (01/01/2022) dans le plan de l'écliptique en AU, AU/j

q0_jup = [4.6580839119399, -1.7945574704081, -0.0967627432056]
q0_sol = [0, 0, 0]
q0_sat = [6.9600794880566, -7.0666123077577, -0.1541266614081]

p0_jup = [0.0026255239604, 0.0074041879836, -0.0000894925649]
p0_sol = [0, 0, 0]
p0_sat = [0.0036673089657, 0.0039101598352, -0.0002138877692]

#%% Dérivées pour 2 corps

q_jup = np.zeros((3, N))
q_jup[:,0] = q0_jup
p_jup = np.zeros((3, N))
p_jup[:,0] = p0_jup

q_sol = np.zeros((3, N))
p_sol = np.zeros((3, N))

dHdp_jup = p_jup/m_jup
dHdp_sol = p_sol/m_sol

dHdq_jup = np.zeros((3, N))
"""for i in range(N) :
    r = q_sol[:,i] - q_jup[:,i]
    dHdq_jup[:,i] = -G*m_jup*m_sol*r/(np.linalg.norm(r))**3"""

dHdq_sol = np.zeros((3, N))

#%% Heun à 2 corps



dHdqt_jup = np.zeros((3, N))
#for i in range(N) :
    
   # print(d)

dHdqt_sol = np.zeros((3, N))

for i in range(N-1) :
    r = q_sol[:,i] - q_jup[:,i]
    #print(r)
    d = np.sqrt(r[0]**2 + r[1]**2 + r[2]**2)
    dHdq_jup[:,i] = -G*m_jup*m_sol*r/(d)**3
    dHdq_sol[:,i] = G*m_jup*m_sol*r/(d)**3
    
    qt_jup = q_jup + k*dHdp_jup 
    pt_jup = p_jup - k*dHdq_jup
    
    qt_sol = q_sol + k*dHdp_sol 
    pt_sol = p_sol - k*dHdq_sol
    
    dHdpt_jup = pt_jup/m_jup
    dHdpt_sol = pt_sol/m_sol
    rt = qt_sol[:,i] - qt_jup[:,i]
    #print(r)
    dt = np.sqrt(rt[0]**2 + rt[1]**2 + rt[2]**2)
    dHdqt_jup[:,i] = -G*m_jup*m_sol*rt/(dt)**3
    dHdqt_sol[:,i] = G*m_jup*m_sol*rt/(dt)**3
    #print(d)
    q_jup[:,i+1] = q_jup[:,i] + (k/2)*(dHdp_jup[:,i] + dHdpt_jup[:,i+1])
    p_jup[:,i+1] = p_jup[:,i] - (k/2)*(dHdq_jup[:,i] + dHdqt_jup[:,i+1])
    q_sol[:,i+1] = q_sol[:,i] + (k/2)*(dHdp_sol[:,i] + dHdpt_sol[:,i+1])
    p_sol[:,i+1] = p_sol[:,i] - (k/2)*(dHdq_sol[:,i] + dHdqt_sol[:,i+1])
    
    
fig = plt.figure(figsize = (4, 4))
ax = fig.add_subplot(111, projection = '3d')

x_jup = q_jup[0,:]
y_jup = q_jup[1,:]
z_jup = q_jup[2,:]  
ax.plot(x_jup, y_jup, z_jup)

plt.show

print(q_jup)






























