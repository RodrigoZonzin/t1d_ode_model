import matplotlib.pyplot as plt
import numpy as np 
from scipy.integrate import solve_ivp
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as mcolors
from matplotlib import colormaps
import pandas as pd
import os, sys 

#nome arquivo de saida
#fout_name = sys.argv[1]

#valores iniciais
G_0 = 80
I_0 = 12
B_0 = 800 #980
T_0 = 1
y0 = [G_0, I_0, B_0, T_0]

t_range = (0, 500)
t_eval = np.linspace(*t_range, int(500/0.1))  #dt =0.5

def system(t, y, params): 
    #variaveis
    G = y[0]
    I = y[1]
    B = y[2]
    T = y[3]

    #parametros 
    RG, kG, muG, sI, muI, alphaB, muB, sigmaI, sigmaB, kB, alphaR, Tnaive, sE, muE = list(params.values())
     
    #equacoes 
    #dGdt = RG*( np.exp(- ((t - 50)/50)**2) + 1.0) - kG*I*G - muG*G
    dGdt = RG - kG*I*G - muG*G
    dIdt = (sI*B*G*G)/(1 + G*G) - muI*I
    dBdt = alphaB*G*B/(1 + sigmaB*B + sigmaI*I) - muB*B -(kB*B*T)/(1+alphaR*T)
    #print(t, sE, Tnaive, T, muE) if (t- 2.00 ) < 0.0000001 else 0+0
    dTdt = sE*(Tnaive - T) -muE*T

    return [dGdt, dIdt, dBdt, dTdt]

params = {
    #dGdt=RG-kG*I-muG*G
    'RG': 5.0,
    'kG': 0.005,
    'muG': 0.01125,

    #dIdt = (sI*B*G*G)/(1 + G*G) - muI*I
    'sI': 0.008,
    'muI': 0.6,
    
    #dBdt = alphaB*G*(1000-B) - muB*B -(kB*B*T)/(1+alphaR*T) 
    'alphaB': 0.25, #0.4, #0.4,
    'muB': 0.2, #0.3,
    'sigmaI': 0.2, 
    'sigmaB': 0.15, 
    'kB': 1,
    'alphaR': 0.001,

    #dTdt = sE*(Tnaive - T) -muE*T          #equação original dTdt = sE*(Tnaive - T) -muE*T*Treg, mas ainda nao temos TReg 
    'Tnaive': 370,
    'sE': 0.05, 
    'muE': 0.001
    
}

#parametros para texto
def params_to_str(p): 
    texto = ''
    for k, v in p.items(): 
        texto = texto+ f'{k}: {v}\n'
    return texto

menores_glicose = []
for kb in np.arange(0, 2.5, 0.1): 
    
    params['kB'] = kb

    #resolvendo o sistema

    sol = solve_ivp(
        fun=lambda t, y: system(t, y, params),
        t_span=t_range,
        y0=y0,
        t_eval=t_eval
    )


    results_glicose = sol.y[0]
    results_Tcells  = sol.y[3]

    menores_glicose.append(np.max(results_glicose))


plt.figure(figsize=(10, 8))
plt.scatter(y = menores_glicose, x = np.arange(0, 2.5, 0.1))
plt.savefig('results/testando_kbs.png', dpi = 400)
