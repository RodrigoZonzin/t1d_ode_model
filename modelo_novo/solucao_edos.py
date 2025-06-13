import matplotlib.pyplot as plt
import numpy as np 
from scipy.integrate import solve_ivp
from matplotlib.backends.backend_pdf import PdfPages
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
    dBdt = alphaB*G*B - muB*B

    return [dGdt, dIdt, dBdt]

#valores iniciais
G_0 = 80
I_0 = 20
B_0 = 1000
y0 = [G_0, I_0, B_0]

t_range = (0, 150)
t_eval = np.linspace(*t_range, 500)
params = {
    'RG': 8,
    'kG': 0.008,
    'muG': 0.001,
    'alphaI': 0.01,
    'muI': 0.8,
    'alphaB': 0.4,
    'muB': 0.3
}

# Resolver o sistema
sol = solve_ivp(
    fun=lambda t, y: system(t, y, params),
    t_span=t_range,
    y0=y0,
    t_eval=t_eval
)

nomesVar = ['Glicose', 'Insulina', 'Beta']


with PdfPages('resultados_EDO.pdf') as pdf:
    for i, nome in enumerate(nomesVar):
        plt.figure(figsize=(20,10))
        plt.plot(sol.t, sol.y[i], label=nome, color='tab:blue')
        plt.xlabel('Tempo')
        plt.ylabel('Concentração')
        plt.title(f'Evolução da {nome}')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        
        pdf.savefig() 
        plt.close()   
