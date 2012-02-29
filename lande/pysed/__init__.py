""" pysed

    An object oriented package for computing spectral energy
    distributions of astrophysical sources in python. 

    This code differs from other SED packages in valuing human readability
    and conceptual clarity over computational efficiency.

    Consequently, this code may not be fast enough for professional work
    and is designed primarily as a teaching aid.

    References:
    Notes: 
        * R&L is Rybicki and Lightman "Radiative Processes in Astrophysics

    Dependencies:
        * numpy, scipy, matplotlib, sympy

    General implementation notes:
        * Keep everything object oriented so that it is flexible.
        * As much as possible, vectorize the calculations using numpy and scipy
        * Try to perform most integrals using the simpson method so
          that the integrand evaluation can be vectorized.
        * When integrating over energy, sample points uniformly in
          log space (see sed_integrate.logsimps).

    This file is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


    TODO:
        * Try to make some of the spectrum calculators (like Bremsstrahlung)
          fully vectorized for increased speed. N.B., this probably
          requires some careful though about which axes logsimps integrates
          over.
        * General method for enforcing units of input to objects.
        * Improve documentation of methods + clarify units/formulas
        * write integrate function which applies to all spectrum objects!
        * Add unit testing for the integration methods


    Author: Joshua Lande <joshualande@gmail.com>
"""
