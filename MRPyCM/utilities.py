from vampyr import vampyr3d as vp
import matplotlib.pyplot as plt
import numpy as np
from qcelemental.physical_constants.context import PhysicalConstantsContext


def computeKappaOut(eps_out, Ionic_strength):
    """Computes the value of the Debye Huckel Screening parameter on
    the outside of the solute volume.

    Args:
        eps_out (float): Relative permittivity of the solvent in T=298.15K
        Ionic_strength (float): Ionic strength of the electrolyte in mol/L. This can also be considered the bulk concentration of the electrolyte.

    Returns:
        float: computed value of the Debye Huckel Screening parameter on the outside of the solute volume in atomic units.
    """
    
    context = PhysicalConstantsContext(context="CODATA2018")
    kb = context.kb
    e = context.get("elementary charge")
    e_0 = context.e0
    N_a = context.na
    m2au = 1/ context.bohr2m
    T = 298.15

    numerator = e_0*eps_out*kb*T
    denominator = 2.0*(e**2)*N_a*1000.0*Ionic_strength
    
    debye_length = ((numerator/denominator)**(1.0/2.0))*m2au
    kappa = 1.0/debye_length
    return kappa


def plotFunction(function, left=-1.0, right=1.0, npoints=100, title=""):
        x_plt = np.linspace(left, right, npoints)
        f_plt = [function([x, 0.0, 0.0]) for x in x_plt]
        plt.plot(x_plt, f_plt)
        plt.title(title)
        plt.show()


def constructChargeDensity(positions, charges, width_parameter):
    """Computes the charge density of a set of point charges using a Gaussian
    function. The Gaussian function is centered at the position of the charge
    and the width parameter is the standard deviation of the Gaussian.

    Args:
        positions (list of lists of floats): list of positions of the charges
        charges (list of floats): list of charges
        width_parameter (int, optional): Standard deviation of the Gaussian. Defaults to 1000.

    Returns:
        function: function that computes the charge density at a given position
    """
    charge_density = vp.GaussExp()
    for (pos, charge) in zip(positions, charges):
        beta = width_parameter
        alpha = (beta / np.pi)**(3.0/2.0)
        charge_density.append(vp.GaussFunc(beta=beta, alpha=alpha*charge, position=pos, poly_exponent=[0,0,0]))
    return charge_density