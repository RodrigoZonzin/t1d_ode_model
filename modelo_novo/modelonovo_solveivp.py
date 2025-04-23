import argparse, contextlib, sys, os
import scipy
import numpy as np

from parameters_pops import *



def initial_values() -> np.ndarray:
    return list(paper_initial_values.values())

def constants() -> list:
    return paper_parameters



def constants_with_names() -> list:
    constants_list = []

    return constants_list


def variable_names() -> list[str]:
    return [
        "G",
        "I",
        "B", 
        "TE",
        "TReg",
        "APC",
        ]

def system(t: np.float64, u: np.ndarray, *constants) -> np.ndarray:
    # populations
    """Variaveis do modelo"""
    G   = u[0]
    I   = u[1]
    B   = u[2]
    TE  = u[3]
    TReg= u[4]
    APC = u[5]


    #constants
    k1          = paper_parameters['k1']
    alphaI      = paper_parameters['alphaI']
    muB         = paper_parameters['muB']
    alphaB      = paper_parameters['alphaB']
    muI         = paper_parameters['muI']
    alphaG      = paper_parameters['alphaG']
    alphaAPCTE  = paper_parameters['alphaAPCTE']
    muB         = paper_parameters['muB']
    alphaTReg   = paper_parameters['alphaTReg']
    muTReg      = paper_parameters['muTReg']
    k2          = paper_parameters['k2']
    k3          = paper_parameters['k3']


    """Equacoes do modelo"""
    #Macrophage population (1)
    dMdt = J + (k+b)*Ma -c*M -fM*M*Ba -fM*M*Bn -e1*M*(M+Ma)

    #Activated Macrophage pop (2)
    dMadt = fM*M*Ba +fM*M*Bn -k*Ma -e2*Ma*(M+Ma)

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
    dDdt = ftD*Bn*(Dss - D - tD) +ftD*Bn*t*D -bDE*E*D -muD*D

    #Tolerogenic Dendritic Cells (9)
    dtDdt = ftD*Ba*(Dss - D - tD) -ftD*Bn*tD -bIR*R*tD -muD*tD

    #Effector T-cells (10)
    dEdt = aE*(Tnaive/Qspleen - E) +bp*(D*E/(thetaD+D)) -ram*E +bE*D*Em -muE*E*R

    #Regulatory T-cells (11)
    dRdt = aR*(Tnaive/Qspleen - R) +bp*(t*D*R/(thetaD+t*D)) -ram*R +bR*tD*Em - muR*E*R

    #Memory T-cells (12)
    dEmdt = ram*(E+R) -(aEm +bE*D +bR*tD)*Em

    return np.array([dMdt, dMadt, dBdt, dBadt, dBndt, dGdt, dIdt, dDdt, dtDdt, dEdt, dRdt, dEmdt])


# includes! "ode-support.py"


def simulation_output_to_csv(sim_steps, simulation_output, write_to):
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


def simulate(filename, st=0, tf=100, dt=0.01, plot=False, x_label="time (days)", y_label="populacao", params={}):
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
    parser.add_argument("--dt", type=float, default=0.001)
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
