import matplotlib.pyplot as plt
import numpy as np 
from scipy.integrate import solve_ivp
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as mcolors
from matplotlib import colormaps
import pandas as pd
import os 


def system(t, y, params): 
    #variaveis
    G = y[0]

    #parametros 
    #RG, kG, muG, alphaI, muI, alphaB, muB = list(params.values())
    r, k = list(params.values())
     
    #equacoes 

    dGdt = r*G*(1-(G/k))

    return [dGdt]

#valores iniciais
G_0 = 2
y0 = [G_0]

t_range = (0, 250)
t_eval = np.linspace(*t_range, int(250/0.1))  #dt =0.5

params = {
    #dGdt=RG-kG*I-muG*G
    'r': 0.2,
    'k': 10
}

#resolvendo o sistema
sol = solve_ivp(
    fun=lambda t, y: system(t, y, params),
    t_span=t_range,
    y0=y0,
    t_eval=t_eval
)

#dataframe pra salvar os resultados
df_results = pd.DataFrame(index = sol.t , columns=['Y'])

nomesVar = ['Y']
for i, nome in enumerate(nomesVar):
    df_results[nome] = sol.y[i]

"""
    Aqui, apenas escolho algumas cores do viridis para manter a identidade visual no texto
"""
cmap_vir = colormaps['viridis']
mycolors = cmap_vir(np.linspace(0, 1, len(nomesVar)+2))


with PdfPages('results/resultados_EDO_teste.pdf') as pdf:
    #plotando cada curva numa pagina
    for i, nome in enumerate(nomesVar):
        plt.figure(figsize=(10,8)) 
        plt.plot(sol.t, sol.y[i], label=nome, color=mycolors[i])
        #plt.ylim(bottom=0)
        plt.xlabel('Tempo')
        plt.ylabel('Concentração')
        plt.title(f'Evolução da {nome}')
        plt.legend()
        plt.tight_layout()

        pdf.savefig() 
        plt.close()   

    #plotando todas as curvas na ultima pagina
    plt.figure(figsize=(10,3), dpi = 400)
    [plt.plot(sol.t, sol.y[j], label = nome, color = mycolors[j]) for j, nome in enumerate(nomesVar)]
    plt.legend()
    plt.tight_layout()
    pdf.savefig()
    plt.close
        
#salvando os resultados 
df_results.to_csv('results/resultados_EDO_teste.csv', index=True) #se Index=True, salva o dt na primeira coluna do .csv