params = {
    #dGdt=RG-kG*I-muG*G
    'RG': 5,
    'kG': 0.8,
    'muG': 0.01125,

    #dIdt = alphaI*B - muI*I
    'alphaI': 0.01,
    'sI': 0.02,
    'muI': 0.8,
    
    #dBdt = alphaB*G*(1000-B) - muB*B
    'alphaB': 0.4,
    'muB': 0.3
}

texto = 'p'
for k, v in params.items():
    texto = texto+ f'{k}: {v}\n'

print(texto)
