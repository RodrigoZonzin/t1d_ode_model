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

t_range = (0, 800)
t_eval = np.linspace(*t_range, int(800/0.1))  #dt =0.5

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

#parametros para texto
def params_to_str(p): 
    texto = ''
    for k, v in p.items(): 
        texto = texto+ f'{k}: {v}\n'
    return texto

# Parâmetros fixos
params = {
    'RG': 5.0, 'kG': 0.005, 'muG': 0.01125,
    'sI': 0.008, 'muI': 0.6,
    'alphaB': 0.25, 'muB': 0.2,
    'sigmaI': 0.2, 'sigmaB': 0.15,
    'kB': 1, 'alphaR': 0.001,
    'Tnaive': 370, 'sE': 0.05, 'muE': 0.001
}

plt.figure(figsize=(10, 6), dpi = 400)
linhas = ['-', '--', (0, (5, 2)), '-.', ':']

#for i, se in enumerate(np.arange(0.00001, 0.00005, 0.00001)):
for i, se in enumerate([0.00001, 0.00005, 0.0001, 0.0005, 0.005]):
    params['sE'] = se

    sol = solve_ivp(
        fun=lambda t, y: system(t, y, params),
        t_span=t_range,
        y0=y0,
        t_eval=t_eval
    )

    plt.plot(sol.t, sol.y[3], label = rf'$s_E = {format(se, '.5f')}$', color = 'black', ls = linhas[i])

plt.legend()
plt.tight_layout()
plt.savefig('impacto_sE_em_T.png', dpi = 400)
print(np.min(sol.y[2]))

#plt.show()