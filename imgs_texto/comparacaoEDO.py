import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Parâmetros
k = 1.0
m = 2.0

# EDO: dy/dt = ky + m
def dydt(t, y):
    return k * y + m

# Solução analítica correta
def y_analitica(t):
    return (m / k) * (np.exp(k * t) - 1)

# Função fornecida pelo usuário (não é a solução da EDO original)
def y_usuario(t):
    return m * (1 - np.exp(k * t))

# Condição inicial
y0 = [0]
t_span = (0, 5)
t_eval = np.linspace(*t_span, 100)

# Solução numérica
sol = solve_ivp(dydt, t_span, y0, t_eval=t_eval)

# Plotando
plt.plot(sol.t, sol.y[0], label='Solução Numérica (solve_ivp)', linestyle='--')
plt.plot(t_eval, y_analitica(t_eval), label='Solução Analítica Correta', linewidth=2)
plt.plot(t_eval, y_usuario(t_eval), label='Função do Usuário (m(1 - e^{kt}))', linestyle=':')
plt.xlabel('t')
plt.ylabel('y(t)')
plt.legend()
plt.grid(True)
plt.title('Comparação entre soluções da EDO dy/dt = ky + m')
plt.savefig('result.png')
#plt.show()

