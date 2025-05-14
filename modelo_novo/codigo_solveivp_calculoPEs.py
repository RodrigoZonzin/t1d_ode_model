import argparse, contextlib, sys, os
import scipy
import numpy as np


def initial_values() -> np.ndarray:
    G_0   = (8/0.01)
    I_0   = 0    
    B_0   = 0 #300
    Te_0  = 0
    Treg_0= 0

    return np.array((G_0, I_0, B_0, Te_0, Treg_0))

def constants() -> list:
    RG      = 8. # 1.
    kG      = 0.008  #0.05 #0.72
    muG     = 0.001
    alphaI  = 0.01
    muI     = 0.8
    alphaB  = 0.5
    kB      = 0.4
    alphaE  = 0.01
    alpha1R = 0.01
    muB     = 0.3 #0.5
    sE      = 0.002*0
    Tnaive  = 1e4
    muE     = 0.02
    sR      = 0.2*0
    alpha2R = 0.1
    muR     = 0.001

    return (RG, kG, muG, alphaI, 
            muI, alphaB, kB, alphaE,
            alpha1R, muB, sE, Tnaive, 
            muE, sR, alpha2R, muR)


def constants_with_names() -> list:
    constants_list = [        
        ("RG", 0.003),
        ("kG", 0.1),
        ("muG", 0.003/100),
        ("alphaI", 0.1),
        ("muI", 0.001),
        ("alphaB", 0.01),
        ("kB", 0.04),
        ("alphaE", 0.01),
        ("alpha1R", 0.01),
        ("muB", 0.001),
        ("sE", 0.02),
        ("Tnaive", 1e3),
        ("muE", 0.02),
        ("sR", 0.2),
        ("alpha2R", 0.01),
        ("muR", 0.001),
        ]
    return constants_list


def variable_names() -> list[str]:
    return [
        "G",
        "I",
        "B",
        "Te",
        "Treg"
        ]

def system(t: np.float64, y: np.ndarray, *constants) -> np.ndarray:
    # populations
    G       = y[0]
    I       = y[1]
    B       = y[2]
    Te      = y[3]
    Treg    = y[4]
    
    # constants
    RG      = constants[0]
    kG      = constants[1]
    muG     = constants[2]
    alphaI  = constants[3]
    muI     = constants[4]
    alphaB  = constants[5]
    kB      = constants[6]
    alphaE  = constants[7]
    alpha1R = constants[8]
    muB     = constants[9]
    sE      = constants[10]
    Tnaive  = constants[11]
    muE     = constants[12]
    sR      = constants[13]
    alpha2R = constants[14]
    muR     = constants[15]


    #Sistema de EDOs
    #dGdt = RG -kG*I*G -muG*G      #atualiza o termo I*G
    #Gmax = 85
    #Bmax = 1000
    dGdt = RG -kG*I*G -muG*G
    dIdt = alphaI*B -muI*I
    dBdt = alphaB*G*(B) -((kB*B*Te)/((1+alphaE*Te) + (alpha1R*Treg))) -muB*B
    dTedt = sE*(Tnaive - Te) -muE*Te*Treg #retira o TE de  sE*Te
    dTregdt = ((sR*Te)/(1+alpha2R*Treg)) -muR*Treg

    return np.array([dGdt, dIdt, dBdt, dTedt, dTregdt])

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


def simulate(filename, st=0, tf=50, dt=0.1, plot=False, x_label="time (days)", y_label="conc/ml", params={}):
    sim_steps = np.arange(st, tf + dt, dt)

    #constants_values = [value for _, value in update_constants_with_params(constants_with_names(), params)]

    simulation_output = scipy.integrate.solve_ivp(
        fun=system,
        t_span=(st, tf + dt * 2),
        y0=initial_values(),
        args=constants(),
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
    parser.add_argument("--ylabel", type=str, default="conc/ml")
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