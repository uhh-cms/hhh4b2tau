from columnflow.util import maybe_import
import order as od
hist = maybe_import("hist")
sympy = maybe_import("sympy")
np = maybe_import("numpy")

morphing_coupling_combinations = (
    # (c3, d4) 
    # removed samples with high values of c3 
    # since only nine samples are neccessary
    (0, 0),
    (0, -1),
    (1, 0),
    (1, 2),
    (2, -1),
    (4, 9),
    (-1, 0),
    (-1, -1),
    (-1.5, -0.5),
)

new = (
    # couplings to be morphed into
    # (c3, d4) 
    (1, -1),
    (2, 0),
)

all_cc = morphing_coupling_combinations + new

def morph(kl,kc,s): # kl,kc couplings, s vector of cross-sections
        # define the matrix with nine scalings
        # box, two triangles, pentagon, interferences
        # kl=c3+1, kc=d4+1
        M = sympy.Matrix([
            [
                (c3+1)**4,
                (c3+1)**3,
                (c3+1)**2,
                (c3+1)**2 * (d4+1),
                (c3+1),
                (c3+1) * (d4+1),
                (d4+1)**2,
                (d4+1),
                1,
            ]
            for c3,d4 in morphing_coupling_combinations
        ])

        # the vector of couplings
        c = sympy.Matrix([
            [kl**4],
            [kl**3],
            [kl**2],  
            [kl**2 * kc],
            [kl],
            [kl*kc],
            [kc**2],
            [kc],
            [1],
        ])

        # actual computation,
        #  i.e., matrix inversion and multiplications with vectors
        M_inv = M.pinv()
        coeffs = c.transpose() * M_inv
        # equivalent to "sigma = coeffs * s" but for histograms
        sigma = sum([float(coeff_i) * s_i for coeff_i, s_i in zip(coeffs, s)])

        return sigma

def morphing_hook(
    task,
    hists: dict[od.Process, hist.Hist],
) -> dict[od.Process, hist.Hist]:
    # create hypothetical HHH_4b2tau sample by morphing given data 
    
    hhh_morph0 = od.Process(
          "hhh_c3_{c3}_d4_{d4}_4b2tau".format(
                      c3=str(new[0][0]).replace("-", "minus").replace(".", "p"),
                      d4=str(new[0][1]).replace("-", "minus").replace(".", "p"),
                      ), id="+",
            label="$(\kappa_3=2, \kappa_4=0)$", color1="#bd1f01")
    
    hhh_morph1 = od.Process(
          "hhh_c3_{c3}_d4_{d4}_4b2tau".format(
                      c3=str(new[1][0]).replace("-", "minus").replace(".", "p"),
                      d4=str(new[1][1]).replace("-", "minus").replace(".", "p"),
                      ), id="+",
            label="$(\kappa_3=3, \kappa_4=1)$", color1="#3f90da")
    
    data_histo = []
    for c3, d4 in morphing_coupling_combinations:
        name = "hhh_c3_{c3}_d4_{d4}_4b2tau".format(
                      c3=str(c3).replace("-", "minus").replace(".", "p"),
                      d4=str(d4).replace("-", "minus").replace(".", "p"),
                      )
        for proc, histo in hists.items():
            if proc.name == name:
                data_histo.append(histo)
                break
        else:
            raise Exception("Histograms are missing! Exactly nine are required for morphing.")

    # create the new HHH shape
    hists[hhh_morph0] = morph(kl=new[0][0]+1,kc=new[0][1]+1,s=data_histo)
    hists[hhh_morph1] = morph(kl=new[1][0]+1,kc=new[1][1]+1,s=data_histo)

    return hists