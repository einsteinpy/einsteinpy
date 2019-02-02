from .RK4thOrder import RK4thOrder

try:
    from scipy.integrate import RK45 as RK45Scipy
except:
    pass