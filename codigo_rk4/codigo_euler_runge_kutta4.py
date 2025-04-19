# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from parameters_pops import *


def K1(G_param): 
    return G_param**2 / (G_param**2 + paper_parameters['Ghb']**2)

#trocar ** para multiplicacao 
def K2(E_param, R_param):
    return (paper_parameters['sE']*E_param)**2 / (1 +  (paper_parameters['sE']*E_param)**2 + (paper_parameters['sR']*R_param)**2)

#Apoptotic wave at 9 days
def W(B_param, t_param): 
    #print(B_param, t_param)
    return 0.1 * B_param * np.exp(-(((t_param - 9) / 9) ** 2))

def ode_system(t, u, constants):
    """
    system of first order differential equations
    t: discrete time step value
    y: state vector (vetor das variaveis/populacoes do modelo)
    """
    


    """
    parametros do modelo: 
    "J", "k", "b", "c", "fM", "e1", "e2", "alphaB", "sigmaB", "etha", "Ghb", "sE", "sR",
    "Bconv", "Qpanc", "d", "fMa", "ftD", "Dss", "fD", "R0", "G0", "SI", "sigmaI",
    "GI", "deltaI", "bDE", "muD", "bIR", "aE", "Tnaive", "Qspleen", "bp", "thetaD",
    "ram", "bE", "muE", "aR", "bR", "muR", "aEm"
    """
    J       = paper_parameters['J']
    k       = paper_parameters['k']
    b       = paper_parameters['b']
    c       = paper_parameters['c']
    fM      = paper_parameters['fM']
    e1      = paper_parameters['e1']
    e2      = paper_parameters['e2']
    alphaB  = paper_parameters['alphaB']
    sigmaB  = paper_parameters['sigmaB']
    etha    = paper_parameters['etha']
    Ghb     = paper_parameters['Ghb']
    sE      = paper_parameters['sE']
    sr      = paper_parameters['sR']
    Bconv   = paper_parameters['Bconv']
    Qpanc   = paper_parameters['Qpanc']
    d       = paper_parameters['d']
    fMa     = paper_parameters['fMa']
    ftD     = paper_parameters['ftD']
    Dss     = paper_parameters['Dss']
    fD      = paper_parameters['fD']
    R0      = paper_parameters['R0']
    G0      = paper_parameters['G0']
    SI      = paper_parameters['SI']
    sigmaI  = paper_parameters['sigmaI']
    GI      = paper_parameters['GI']
    deltaI  = paper_parameters['deltaI']
    bDE     = paper_parameters['bDE']
    muD     = paper_parameters['muD']
    bIR     = paper_parameters['bIR']
    aE      = paper_parameters['aE']
    Tnaive  = paper_parameters['Tnaive']
    Qspleen = paper_parameters['Qspleen']
    bp      = paper_parameters['bp']
    thetaD  = paper_parameters['thetaD']
    ram     = paper_parameters['ram']
    bE      = paper_parameters['bE']
    muE     = paper_parameters['muE']
    aR      = paper_parameters['aR']
    bR      = paper_parameters['bR']
    muR     = paper_parameters['muR']
    aEm     = paper_parameters['aEm']


    """Variaveis do modelo"""
    
    M   = u[0]
    Ma  = u[1]
    B   = u[2]
    Ba  = u[3]
    Bn  = u[4]
    G   = u[5]
    I   = u[6]
    D   = u[7]
    tD  = u[8]
    E   = u[9]
    R   = u[10]
    Em  = u[11]


    """Equacoes do modelo"""
    #Macrophage population (1)
    dMdt = J + (k+b)*Ma -c*M -fM*M*Ba -fM*M*Bn -e1*M*(M+Ma)

    #Activated Macrophage pop (2)
    dMadt = fM*M*Ba + fM*M*Bn -e2*Ma*(M+Ma)
    
    #Healthy beta cells (3)
    dBdt = alphaB*K1(G)*B -sigmaB*B -etha*K2(E,R)*B -W(B,t)

    #Apoptotic beta cells (4)
    conv_cte = Bconv/Qpanc
    dBadt = (sigmaB*conv_cte)*B + (etha*conv_cte)*K2(E,R)*B +(W(B,t)*conv_cte) -d*Ba -fM*M*Ba - fMa*Ma*Ba -ftD*(Dss - D)*Ba -fD*D*Ba

    #Necrotic beta cells (5)
    dBndt = d*Ba -fM*M*Bn -fMa*Ma*Bn -ftD*(Dss - D)*Bn -fD*D*Bn

    #Glucose (6)
    dGdt = R0 - (G0+SI*I)*G

    #Insulin (7)
    dIdt = deltaI * (G**2/(G**2 + G*I**2))*B -sigmaI*I

    #Immunogenic Dendritic Cells (8)
    dDdt = ftD*Bn*(Dss -D -tD) +ftD*Bn*t*D -bDE*E*D -muD*D

    #Tolerogenic Dendritic Cells (9)
    dtDdt = ftD*Ba*(Dss-D-t*D) -ftD*Bn*t*D -bIR*R*t*D -muD*t*D

    #Effector T-cells (10)
    dEdt = aE*(Tnaive/Qspleen - E) +bp*(D*E/(thetaD+D)) -ram*E +bE*D*Em -muE*E*R

    #Regulatory T-cells (11)
    dRdt = aR*(Tnaive/Qspleen - R) +bp*(t*D*R/(thetaD+t*D)) -ram*R +bR*t*D*Em - muR*E*R

    #Memory T-cells (12)
    dEmdt = ram*(E+R) -(aEm +bE*D +bR*t*D)*Em 
    

    return np.array([dMdt, dMadt, dBdt, dBadt, dBndt, dGdt, dIdt, dDdt, dtDdt, dEdt, dRdt, dEmdt])

def rk4(f, tk, _uk, _dt=0.01, **kwargs):
    """
    single-step fourth-order numerical integration (RK4) method
    f: system of first order ODEs
    tk: current time step
    _uk: current state vector [y1, y2, y3, ...]
    _dt: discrete time step size
    **kwargs: additional parameters for ODE system
    returns: y evaluated at time k+1
    """
    # evaluate derivative at several stages within time interval
    k1 = f(tk, _uk, **kwargs)
    k2 = f(tk + _dt / 2, _uk + (k1 * (_dt / 2)), **kwargs)
    k3 = f(tk + _dt / 2, _uk + (k2 * (_dt / 2)), **kwargs)
    k4 = f(tk + _dt, _uk + (k3 * _dt), **kwargs)

    # return an average of the derivative over tk, tk + dt
    return _uk + (_dt / 6) * (k1 + (2 * k2) + (2 * k3) + k4)

def euler(f, tk, _uk, _dt=0.001, **kwargs):
    """
    single-step explicit method
    func: system of first order ODEs
    tk: current time step
    _uk: current state vector [y1, y2, y3, ...]
    _dt: discrete time step size
    **kwargs: additional parameters for ODE system
    returns: y evaluated at time k+1
    """
    return _uk + _dt*f(tk, _uk, **kwargs)

def solveSystem(time, dt, y0, method):
    t = 0
    yk = y0
    state = []
    parameters = paper_parameters.values()
    #print(parameters)

    if method == "euler":
        for t in time:
            state.append(yk)
            yk = euler(ode_system, t, yk, dt, constants=parameters)

    elif method == "rk4":
        for t in time:
            state.append(yk)
            yk = rk4(ode_system, t, yk, dt, constants=parameters)

    return np.array(state)

def save(time, state, names):
    df = pd.DataFrame(state, columns = names)
    df.insert(0, 'time', time)
    df.to_csv('results.csv', float_format='%.5f', sep=',', index=False)

def plot(time, state, names):
    fig, ax = plt.subplots()
    fig.set_size_inches(8, 6)
    for i in range(len(names)):
        ax.plot(time, state[:, i], label=names[i], linewidth='2')

    #time2 = np.arange(0, tfinal + 0.001, 0.001)
    #sol = np.exp(r*time2)
    #ax.plot(time2, sol, label='solucao analitica', linewidth='2')
    ax.set(xlabel='dias', ylabel='populacao')
    plt.legend(loc='best')
    fig.savefig('ode_t1d.png', format='png')
    plt.show()

if __name__ == "__main__":
    names = ['M', 'Ma', 'B', 'Ba', 'Bn', 'G', 'I', 'D', 'tD', 'E', 'R', 'Em']
    dt = 0.00001
    tfinal = 10
    time = np.arange(0, tfinal + dt, dt)
    initial_condition = list(initial_values.values())
    
    result = solveSystem(time, dt, initial_condition, "rk4")
    save(time, result, names)
    plot(time, result, names)
