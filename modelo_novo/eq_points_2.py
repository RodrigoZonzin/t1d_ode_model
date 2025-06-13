import sympy as sp

# Definição das variáveis simbólicas
G, I, B = sp.symbols('G I B')
RG, kG, muG = sp.symbols('RG kG muG')
alphaI, muI = sp.symbols('alphaI muI')
alphaB, kB, muB = sp.symbols('alphaB kB muB')

# Equações diferenciais
def dGdt(G, I):
    return RG - kG * I * G - muG * G

def dIdt(B, I):
    return alphaI * B - muI * I

def dBdt(G, B):
    return alphaB * G * (1000 - B) - muB * B

# Substituições:
# 1. I = (alphaI / muI) * B
I_expr = (alphaI / muI) * B

# 2. B = (1000 * alphaB * G) / (alphaB * G + muB)
B_expr = (1000 * alphaB * G) / (alphaB * G + muB)

# 3. Substituir B em I
I_expr = I_expr.subs(B, B_expr)

# 4. Substituir B e I na equação dG/dt
G_eq = dGdt(G, I_expr)

# Resolver G
G_sol = sp.solve(sp.Eq(G_eq, 0), G, simplify=True)

# Para cada solução de G, encontrar os valores correspondentes de B e I
print("Pontos de equilíbrio encontrados:\n")
for g_sol in G_sol:
    g_sol_simpl = sp.simplify(g_sol)
    b_sol = sp.simplify(B_expr.subs(G, g_sol_simpl))
    i_sol = sp.simplify(I_expr.subs(G, g_sol_simpl))

    print(f"G = {g_sol_simpl}")
    print(f"B = {b_sol}")
    print(f"I = {i_sol}")
    print("-" * 40)
