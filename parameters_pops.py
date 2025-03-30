

paper_parameters = {
    "J": 5e4,           # cells ml−1 d−1 - normal resting macrophage influx
    "k": 0.4,           # d−1 - Ma deactivation rate
    "b": 0.09,          # d−1 - macrophage recruitment rate by Ma
    "c": 0.1,           # d−1 - macrophage egress rate
    "fM": 0.0623e-5,    # ml cell−1 d−1 - basal phagocytosis rate per M
    "e1": 1e-8,         # cell−1 d−1 - anti-crowding term for macrophages
    "e2": 1e-8,         # cell−1 d−1 - anti-crowding term for macrophages
    "alphaB": 0.0334,   # d−1 - β-cell growth rate
    "deltaB": 1/60,     # d−1 - β-cell death rate
    "etha": 0.02,       # d−1 - rate at which T-cells eliminate β-cells
    "Ghb": 90,          # mg dl−1 - glucose level of half-max β-cell production
    "sE": 1,            # ml cells−1 - Relative impact of effector T-cells on β-cell death
    "sR": 36,           # ml cells−1 - Relative impact of regulatory T-cells on β-cell death
    "Bconv": 2.59e5,    # cell mg−1 - β-cells per milligram
    "Qpanc": 0.194,     # ml - volume of mouse pancreas
    "d": 0.50,          # d−1 - β-cell death rate
    "fMa": 0.0623e-5,   # ml cells−1 d−1 - activated phagocytosis rate per Ma
    "ftD": 1.1899e-6,   # ml cells−1 d−1 - rate naive DC engulf apoptotic β-cells
    "Dss": 1e5,         # cells ml−1 - steady-state DC population
    "fD": 1.7101e-7,    # ml cells−1 d−1 - rate naive DC engulf necrotic β-cells
    "R0": 864,          # mg dl−1 - basal rate glucose production
    "G0": 1.44,         # d−1 - rate of glucose decay
    "SI": 0.72,         # ml µU−1 d−1 - insulin rate of glucose elimination
    "sigmaI": 43.2,     # µU ml−1 d−1 mg−1 - maximum rate of insulin production by β-cells
    "GI": 141.4214,     # mg dl−1 - glucose level of half-max insulin production
    "deltaI": 432,      # d−1 - rate of insulin decay
    "bDE": 0.487e-5,    # ml cells−1 d−1 - DC elimination rate by effector T-cells
    "muD": 0.51,        # d−1 - DC death rate
    "bIR": 0.487e-5,    # ml cells−1 d−1 - DC elimination rate by regulatory T-cells
    "aE": 0.1199,       # d−1 - rate of initial expansion of naive T-cells to effector T-cells
    "Tnaive": 370,      # cells - number of naive T-cells contributing to initial production of effector and regulatory T-cells
    "Qspleen": 0.1,     # ml - volume of mouse spleen
    "bp": 12,           # d−1 - maximum expansion rate of effector T-cells due to DCs
    "thetaD": 2.12e5,   # d−1 - DC value for half-maximal effector T-cell expansion
    "ram": 0.01,        # d−1 - reversion rate of T-cells to memory T-cells
    "bE": 1e-3,         # ml d cells−1 - activation rate for effector T-cells from memory T-cells
    "muE": 2e-6,        # d−1 - effector T-cell removal rate due to competition
    "aR": 0.1199,       # d−1 - rate of initial expansion of naive T-cells to regulatory T-cells
    "bR": 1e-3,         # ml d cells−1 - activation rate for regulatory T-cells from Em T-cells
    "muR": 2e-6,        # d−1 - regulatory T-cell removal rate due to competition
    "aEm": 0.01         # d−1 - memory T-cell death rate
}

