import matplotlib.pyplot as plt
import numpy as np 
from scipy.integrate import solve_ivp
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as mcolors
from matplotlib import colormaps
import pandas as pd
import os, sys 

#nome arquivo de saida
fout_name = sys.argv[1]

#valores iniciais
G_0     = 80
I_0     = 12
B_0     = 800 #980
T_0     = 1
TReg_0  = 0
y0 = [G_0, I_0, B_0, T_0, TReg_0]

t_range = (0, 800)
t_eval = np.linspace(*t_range, int(800/0.1))  #dt =0.5

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
    'alpha_Te': 0.001,
    'k_Te': 0.2,
    'm_te': 0.01,

    #dTregdt = (s_treg*Te)/(1 + sigma_treg*Treg) -m_treg*Treg
    's_treg': 0.0, 
    'sigma_treg': 0.05,
    'm_treg': 0.1
}

#parametros para texto
def params_to_str(p): 
    texto = ''
    for k, v in p.items(): 
        texto = texto+ f'{k}: {v}\n'
    return texto

#resolvendo o sistema
sol = solve_ivp(
    fun=lambda t, y: system(t, y, params),
    t_span=t_range,
    y0=y0,
    t_eval=t_eval
)

#dataframe pra salvar os resultados
nomesVar = ['Glicose', 'Insulina', 'Beta', 'T', 'TReg']

df_results = pd.DataFrame(index = sol.t , columns=nomesVar)


for i, nome in enumerate(nomesVar):
    df_results[nome] = sol.y[i]

"""
    Aqui, apenas escolho algumas cores do viridis para manter a identidade visual no texto
"""
cmap_vir = colormaps['viridis']
mycolors = cmap_vir(np.linspace(0, 1, len(nomesVar)+2))


with PdfPages(f'results/{fout_name}.pdf') as pdf:
    #plotando cada curva numa pagina
    for i, nome in enumerate(nomesVar):
        plt.figure(figsize=(8,6)) 
        plt.plot(sol.t, sol.y[i], label=nome, color=mycolors[i], lw = 2)
        #plt.ylim(bottom=0)
        plt.xlabel('Tempo')
        plt.ylabel('Concentração')
        #plt.title(f'Evolução da {nome}')
        plt.tick_params(axis='both', which='major', labelsize=12)
        plt.legend()
        plt.tight_layout()

        pdf.savefig() 
        plt.savefig(f'results/{fout_name}_{nome}.png')
        plt.close()   

    #plotando todas as curvas na ultima pagina
    plt.figure(figsize=(15,5), dpi = 500)
    [plt.plot(sol.t, sol.y[j], label = nome, color = mycolors[j], lw = 2) for j, nome in enumerate(nomesVar)]
    plt.legend()
    plt.xlabel('Tempo')
    plt.ylabel('Concentração')
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.tight_layout()
    pdf.savefig()
    plt.savefig(f'results/{fout_name}.png', dpi = 500)
    plt.close
        
    plt.figure(figsize=(10,8), dpi = 500)
    plt.text(0.1, 0.9, params_to_str(params), fontsize=12, va='top', ha='left', wrap=True)
    plt.axis('off')
    pdf.savefig()
    plt.close()

#salvando os resultados 
df_results.to_csv(f'results/{fout_name}.csv', index=True) #se Index=True, salva o dt na primeira coluna do .csv