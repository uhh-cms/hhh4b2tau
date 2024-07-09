import numpy as np
import matplotlib.pyplot as plt

coupling_combinations = (
    # (c3, d4)
    (0, 0),
    (0, 99),
    (0, -1),
    (19, 19),
    (1, 0),
    (1, 2),
    (2, -1),
    (4, 9),
    (-1, 0),
    (-1, -1),
    (-1.5, -0.5),
)
coupling_combinations_xsec = np.array([
    # (c3, d4, xsecLO)
    [0, 0, 3.705e-05],
    [0, 99, 5.866e-03],
    [0, -1, 4.104e-05],
    [19, 19, 1.474e-01],
    [1, 0, 2.909e-05],
    [1, 2, 1.612e-05],
    [2, -1, 5.789e-05],
    [4, 9, 2.451e-04],
    [-1, 0, 1.130e-04],
    [-1, -1, 1.089e-04],
    [-1.5, -0.5, 1.938e-04],
])

def xscaler(c3, d4):
    return (
        1 - 0.79*c3 - 0.10*d4 + 0.81*c3**2 - 0.16*c3*d4 + 1.6e-2*d4**2
        - 0.23*c3**3 + 4.5e-2*c3**2*d4 + 3.5e-2*c3**4
    )
c3 = np.linspace(-15 ,15 , 61)
d4 = np.linspace(-99, 99, 397)
c3, d4 = np.meshgrid(c3, d4)
xsec = xscaler(c3, d4)
plt.pcolormesh(c3, d4, xsec)
plt.colorbar()
plt.grid()
plt.xlabel("c3")
plt.ylabel("d4")
plt.title("Cross-section scaling (theoretical values for 27 TeV)")
plt.show()
from IPython import embed; embed()
plt.pcolormesh(coupling_combinations_xsec[:,0], 
               coupling_combinations_xsec[:,1], 
               coupling_combinations_xsec[:,2]/3.705e-05)
plt.colorbar()
plt.grid()
plt.show()
