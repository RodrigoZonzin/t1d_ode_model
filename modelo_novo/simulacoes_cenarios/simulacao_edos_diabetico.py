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
TReg_0 = 1
y0 = [G_0, I_0, B_0, T_0, TReg_0]

t_range = (0, 100)
t_eval = np.linspace(*t_range, int(100/0.1))  #dt =0.5

def system(t, y, params): 
    #variaveis
    G   = y[0]
    I   = y[1]
    B   = y[2]
    Te  = y[3]
    Treg = y[4]

    #parametros 
    RG, kG, muG, sI, muI, alphaB, muB, sigmaI, sigmaB, k_B, alpha_Te, k_Te, m_te, s_treg, sigma_treg, m_treg = list(params.values())
     
    #equacoes 
    #dGdt = RG*( np.exp(- ((t - 50)/50)**2) + 1) - kG*I*G - muG*G
    dGdt = RG -kG*I*G -muG*G
    dIdt = (sI*B*G*G)/(1 + G*G) - muI*I
    dBdt = alphaB*G*B/(1 + sigmaB*B + sigmaI*I) - muB*B - k_B*B*Te
    dTedt = alpha_Te*B/(1 + Treg) - k_Te*Te*Treg - m_te*Te
    dTregdt = (s_treg*Te)/(1 + sigma_treg*Treg) -m_treg*Treg

    return [dGdt, dIdt, dBdt, dTedt, dTregdt]

#parametros para texto
def params_to_str(p): 
    texto = ''
    for k, v in p.items(): 
        texto = texto+ f'{k}: {v}\n'
    return texto

params = {
    #dGdt=RG-kG*I-muG*G
    'RG': 5.0,
    'kG': 0.0049,
    'muG': 0.01125,

    #dIdt = (sI*B*G*G)/(1 + G*G) - muI*I
    'sI': 0.005,
    'muI': 0.34,
    
    #dBdt = alphaB*G*B/(1 + sigmaB*B + sigmaI*I) - muB*B - k_B*B*Te
    'alphaB': 0.39, #0.4,
    'sigmaB': 0.2,
    'sigmaI': 0.15, 
    'muB': 0.25, #0.3,
    'k_B': 0.15,

    #dTedt = alpha_Te*B/(1 + Treg) - k_Te*Te*Treg - m_te*Te
    'alpha_Te': 0.5,
    'k_Te': 0.2,
    'm_te': 0.01,

    #dTregdt = (s_treg*Te)/(1 + sigma_treg*Treg) -m_treg*Treg
    's_treg': 0.0, 
    'sigma_treg': 0.05,
    'm_treg': 0.1
}

plt.figure(figsize=(10, 6), dpi = 400)
linhas = ['-', ':', (0, (5, 2)), '-.', '--']
cores   = colormaps['viridis']
mcores = cores(np.linspace(0, 1, 10))

#for i, se in enumerate(np.arange(0.00001, 0.00005, 0.00001)):
#for i, kB in enumerate([0.01, 0.05, 0.15, 0.60, 0.98]):
for i, alpha in enumerate([0.01, 0.1, 0.25, 0.60, 0.98]):
    params['alpha_Te'] = alpha

    sol = solve_ivp(
        fun=lambda t, y: system(t, y, params),
        t_span=t_range,
        y0=y0,
        t_eval=t_eval
    )

    plt.plot(sol.t, sol.y[2], label = rf'$\alpha_{{TE}} = {alpha:.2f}$', color = mcores[i], ls = linhas[i])


plt.legend()
plt.tight_layout()
plt.savefig('impacto_alphaTe_em_B.png', dpi = 400)
print(np.min(sol.y[2]))

#plt.show()