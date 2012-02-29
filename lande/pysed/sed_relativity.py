""" Functions to perform basic relativistic calculations.

    Author: Joshua Lande <joshualande@gmail.com>
"""
import numpy as np

def gamma_to_beta(gamma):
    r""" Note, 
            \gamma = {1-\beta^2)^{-1}, 
        so
            \beta = \sqrt{1-\gamma^{-2}}.
        """
    return np.sqrt(1-gamma**-2)
