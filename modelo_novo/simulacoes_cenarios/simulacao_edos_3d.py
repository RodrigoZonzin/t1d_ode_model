import matplotlib.pyplot as plt
import numpy as np 
from scipy.integrate import solve_ivp
from matplotlib import cm

# valores iniciais
G_0 = 80
I_0 = 12
B_0 = 800
T_0 = 1
TReg_0 = 1
y0 = [G_0, I_0, B_0, T_0, TReg_0]

t_range = (0, 800)
t_eval = np.linspace(*t_range, int(800/0.1))

# Sistema de EDOs
def system(t, y, params): 
    G, I, B, Te, Treg = y
    RG, kG, muG, sI, muI, alphaB, muB, sigmaI, sigmaB, k_B, alpha_Te, k_Te, m_te, s_treg, sigma_treg, m_treg = params.values()
    
    dGdt = RG -kG*I*G -muG*G
    dIdt = (sI*B*G*G)/(1 + G*G) - muI*I
    dBdt = alphaB*G*B/(1 + sigmaB*B + sigmaI*I) - muB*B - k_B*B*Te
    dTedt = alpha_Te*B/(1 + Treg) - k_Te*Te*Treg - m_te*Te
    dTregdt = (s_treg*Te)/(1 + sigma_treg*Treg) -m_treg*Treg

    return [dGdt, dIdt, dBdt, dTedt, dTregdt]

params = {
    'RG': 5.0, 
    'kG': 0.0049,
    'muG': 0.01125,
    'sI': 0.005, 
    'muI': 0.34,
    'alphaB': 0.39,
    'sigmaB': 0.2,
    'sigmaI': 0.15,
    'muB': 0.25,
    'k_B': 0.07,
    'alpha_Te': 0.5,
    'k_Te': 0.2,
    'm_te': 0.01,
    's_treg': 0.0,
    'sigma_treg': 0.05,
    'm_treg': 0.1
}

alpha_Te_values = np.linspace(0.1, 1, 20)
k_B_values = np.linspace(0.1, 1, 20)

B_final = np.zeros((len(alpha_Te_values), len(k_B_values)))

#fout = open('results3d.csv', 'w+')

for i, alpha in enumerate(alpha_Te_values):
    for j, kB in enumerate(k_B_values):
        params['alpha_Te'] = alpha
        params['k_B'] = kB

        sol = solve_ivp(lambda t, y: system(t, y, params), t_range, y0, t_eval=t_eval)
        #fout.write(f'{alpha}, {kB}, {np.min(sol.y[2])}\n')
        if np.min(sol.y[2]) < 0: 
            B_final[i, j] = 0
        else: 
            B_final[i, j] = np.max(sol.y[0]) 

#fout.close()
# Criar gráfico 3D
fig = plt.figure(figsize=(8, 6), dpi=300)
ax = fig.add_subplot(111, projection='3d')
A, K = np.meshgrid(k_B_values, alpha_Te_values)

surf = ax.plot_surface(A, K, B_final, cmap=cm.viridis)
ax.view_init(elev=30, azim=-45)
ax.set_xlabel(r'$s_I$')
ax.set_ylabel(r'$\alpha_{Te}$')
ax.set_zlabel(r'min I(t)')
#ax.set_title('Impacto de $\\alpha_{Te}$ e $k_{Te}$ em B')

fig.colorbar(surf, ax=ax, shrink=0.6, aspect=8)
plt.tight_layout()
plt.savefig('impacto_alphaTe_kB_em_G_3D.png', dpi=300)
#plt.show()
