import argparse, contextlib, sys, os
import scipy
import numpy as np

from parameters_pops import *



def initial_values() -> np.ndarray:
    return list(paper_initial_values.values())

def constants() -> list:
    return paper_parameters


def K1(G_param):
    return G_param**2 / (G_param**2 + paper_parameters['Ghb']**2)

#trocar ** para multiplicacao 
def K2(E_param, R_param):
    return (paper_parameters['sE']*E_param)**2 / (1 +  (paper_parameters['sE']*E_param)**2 + (paper_parameters['sR']*R_param)**2)

#Apoptotic wave at 9 days
def W(B_param, t_param):
    #print(B_param, t_param)
    return 0.1 * B_param * np.exp(-(((t_param - 9) / 9) ** 2))


def constants_with_names() -> list:
    constants_list = [        
        ('J', 50000.0), 
        ('k', 0.4), 
        ('b', 0.09), 
        ('c', 0.1), 
        ('fM', 6.23e-07), 
        ('e1', 1e-08), 
        ('e2', 1e-08), 
        ('alphaB', 0.0334),
        ('sigmaB', 0.016666666666666666), 
        ('etha', 0.02), 
        ('Ghb', 90), 
        ('sE', 1), 
        ('sR', 36), 
        ('Bconv', 259000.0), 
        ('Qpanc', 0.194), 
        ('d', 0.5), 
        ('fMa', 6.23e-07), 
        ('ftD', 1.1899e-06), 
        ('Dss', 100000.0), 
        ('fD', 1.7101e-07), 
        ('R0', 864), 
        ('G0', 1.44), 
        ('SI', 0.72), 
        ('sigmaI', 43.2), 
        ('GI', 141.4214), 
        ('deltaI', 432), 
        ('bDE', 4.87e-06), 
        ('muD', 0.51), 
        ('bIR', 4.87e-06), 
        ('aE', 0.1199), 
        ('Tnaive', 370), 
        ('Qspleen', 0.1), 
        ('bp', 12), 
        ('thetaD', 212000.0), 
        ('ram', 0.01), 
        ('bE', 0.001), 
        ('muE', 2e-06), 
        ('aR', 0.1199), 
        ('bR', 0.001), 
        ('muR', 2e-06), 
        ('aEm', 0.01)]              
    return constants_list


def variable_names() -> list[str]:
    return [
        "M",
        "Ma",
        "B",
        "Ba",
        "Bn",
        "G",
        "I",
        "D", 
        "tD",
        "E",
        "R",
        "Em"
        ]

def system(t: np.float64, u: np.ndarray, *constants) -> np.ndarray:
    # populations
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

    #constants
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


    """Equacoes do modelo"""
    #Macrophage population (1)
    dMdt = J + (k+b)*Ma -c*M -fM*M*Ba -fM*M*Bn -e1*M*(M+Ma)

    #Activated Macrophage pop (2)
    dMadt = fM*M*Ba +fM*M*Bn -k*Ma -e2*Ma*(M+Ma)

    #Healthy beta cells (3)
    dBdt = alphaB*K1(G)*B -sigmaB*B -etha*K2(E,R)*B -W(B,t)

    #Apoptotic beta cells (4)
    conv_cte = Bconv/Qpanc
    dBadt = (sigmaB*conv_cte)*B + (etha*conv_cte)*K2(E,R)*B +(W(B,t)*conv_cte) -d*Ba -fM*M*Ba - fMa*Ma*Ba -ftD*(Dss - D -tD)*Ba -fD*D*Ba

    #Necrotic beta cells (5)
    dBndt = d*Ba -fM*M*Bn -fMa*Ma*Bn -ftD*(Dss - D -tD)*Bn -fD*D*Bn

     #Glucose (6)
    dGdt = R0 - (G0+SI*I)*G

    #Insulin (7)
    dIdt = deltaI * (G**2/(G**2 + G*I**2))*B -sigmaI*I

    #Immunogenic Dendritic Cells (8)
    dDdt = ftD*Bn*(Dss -D -tD) +ftD*Bn*t*D -bDE*E*D -muD*D

    #Tolerogenic Dendritic Cells (9)
    dtDdt = ftD*Ba*(Dss-D-tD) -ftD*Bn*tD -bIR*R*tD -muD*tD

    #Effector T-cells (10)
    dEdt = aE*(Tnaive/Qspleen - E) +bp*(D*E/(thetaD+D)) -ram*E +bE*D*Em -muE*E*R

    #Regulatory T-cells (11)
    dRdt = aR*(Tnaive/Qspleen - R) +bp*(t*D*R/(thetaD+t*D)) -ram*R +bR*tD*Em - muR*E*R

    #Memory T-cells (12)
    dEmdt = ram*(E+R) -(aEm +bE*D +bR*tD)*Em

    return np.array([dMdt, dMadt, dBdt, dBadt, dBndt, dGdt, dIdt, dDdt, dtDdt, dEdt, dRdt, dEmdt])


# includes! "ode-support.py"


def simulatiion_output_to_csv(sim_steps, simulation_output, write_to):
    if not simulation_output.success:
        print(simulation_output.message)
        return

    populatio_values_per_dt = simulation_output.y.T

    write_to.write(f"t,{','.join(variable_names())}\n")

    for dt, y in zip(sim_steps, populatio_values_per_dt):
        write_to.write(f"{dt},")
        write_to.write(",".join(f"{val:.4f}" for val in y))
        write_to.write("\n")


COLORS = [
    'tab:blue',
    'tab:orange',
    'tab:green',
    'tab:red',
    'tab:purple',
    'tab:brown',
    'tab:pink',
    'tab:gray',
    'tab:olive',
    'tab:cyan',
]

def plot_simulation(sim_steps, simulation_output, filename, x_label="time (days)", y_label="conc/ml"):
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    with PdfPages(filename) as pdf:
        # All
        all_fig, all_ax = plt.subplots()
        all_fig.set_size_inches(8, 6)
        all_ax.set(title="", xlabel=x_label, ylabel=y_label)

        # Individually
        for i, (variable_name, variable_line_data) in enumerate(zip(variable_names(), simulation_output.y)):
            fig, ax = plt.subplots()
            fig.set_size_inches(8, 6)
            ax.set(
                title=variable_name,
                xlabel=x_label, 
                ylabel=y_label, 
            )            
            ax.plot(simulation_output.t, variable_line_data, color=COLORS[i % len(COLORS)])
            all_ax.plot(simulation_output.t, variable_line_data)

            pdf.savefig(fig)
        all_ax.legend(variable_names(),loc="best")
        pdf.savefig(all_fig)


def file_or_stdout(filename: str | None):
    if filename:
        return open(filename, 'w')
    else:
        return sys.stdout


def update_constants_with_params(constants, params):
    updated_constants = constants.copy()

    constant_names = [constant[0] for constant in constants]

    for name, value in params.items():
        for idx, (const_name, const_value) in enumerate(updated_constants):
            if const_name == name:
                updated_constants[idx] = (const_name, value)

    return updated_constants


def simulate(filename, st=0, tf=100000, dt=0.1, plot=False, x_label="time (days)", y_label="populacao", params={}):
    sim_steps = np.arange(st, tf + dt, dt)

    constants_values = [value for _, value in update_constants_with_params(constants_with_names(), params)]

    simulation_output = scipy.integrate.solve_ivp(
        fun=system,
        t_span=(st, tf + dt * 2),
        y0=initial_values(),
        args=tuple(constants_values),
        t_eval=sim_steps,
    )

    if plot:
        plot_simulation(sim_steps, simulation_output, filename, x_label, y_label)
    else:
        with file_or_stdout(filename) as f:
            simulation_output_to_csv(sim_steps, simulation_output, f)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--st", type=float, default=0)
    parser.add_argument("--tf", type=float, default=50)
    parser.add_argument("--dt", type=float, default=0.01)
    parser.add_argument("-o", "--output", default=None)
    parser.add_argument("--csv", action=argparse.BooleanOptionalAction)
    parser.add_argument("--xlabel", type=str, default="time (days)")
    parser.add_argument("--ylabel", type=str, default="populacao")
    parser.add_argument("--params", type=str, default="")

    args = parser.parse_args()

    if args.params:
        params = {k: float(v) for k, v in (param.split('=') for param in args.params.split())}
    else:
        params = {}

    simulate(
        args.output,
        plot=not args.csv,
        st=args.st,
        tf=args.tf,
        dt=args.dt,
        x_label=args.xlabel,
        y_label=args.ylabel,
        params=params
    )
