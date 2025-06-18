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
    I = y[1]
    B = y[2]

    #parametros 
    RG, kG, muG, alphaI, muI, alphaB, muB = list(params.values())
     
    #equacoes 
    dGdt = RG-kG*I-muG*G
    dIdt = alphaI*B - muI*I
    dBdt = alphaB*G*(1000-B) - muB*B

    return [dGdt, dIdt, dBdt]

#valores iniciais
G_0 = 80
I_0 = 12
B_0 = 980
y0 = [G_0, I_0, B_0]

t_range = (0, 250)
t_eval = np.linspace(*t_range, int(250/0.1))  #dt =0.5

params = {
    #dGdt=RG-kG*I-muG*G
    'RG': 1,
    'kG': 0.008,
    'muG': 0.01125,

    #dIdt = alphaI*B - muI*I
    'alphaI': 0.01,
    'muI': 0.8,
    
    #dBdt = alphaB*G*(1000-B) - muB*B
    'alphaB': 0.4,
    'muB': 0.3
}

#resolvendo o sistema
sol = solve_ivp(
    fun=lambda t, y: system(t, y, params),
    t_span=t_range,
    y0=y0,
    t_eval=t_eval
)

#dataframe pra salvar os resultados
df_results = pd.DataFrame(index = sol.t , columns=['Glicose', 'Insulina', 'Beta'])

nomesVar = ['Glicose', 'Insulina', 'Beta']
for i, nome in enumerate(nomesVar):
    df_results[nome] = sol.y[i]

"""
    Aqui, apenas escolho algumas cores do viridis para manter a identidade visual no texto
"""
cmap_vir = colormaps['viridis']
mycolors = cmap_vir(np.linspace(0, 1, len(nomesVar)+2))


with PdfPages('results/resultados_EDO.pdf') as pdf:
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
df_results.to_csv('results/resultados_EDO.csv', index=True) #se Index=True, salva o dt na primeira coluna do .csv