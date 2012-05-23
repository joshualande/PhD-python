import os
import sys

from skymaps import SkyDir

from uw.like.pointspec_helpers import PointSource
from uw.like.Models import PowerLaw

from . fit import fit_prefactor

def new_ps(roi, name, l, b, **tsmap_kwargs):
    """ A 'throw away' convenience function to add a new source to
        the ROI, localize it, and then print out a string which can be
        used to modify an ROI to add a new source. """

    skydir = SkyDir(l,b,SkyDir.GALACTIC)

    roi.print_summary(galactic=True)

    model=PowerLaw()

    ps = PointSource(
        name=name,
        skydir=skydir,
        model=model)

    roi.add_source(ps)

    roi.print_summary(galactic=True)

    fit_prefactor(roi, name, use_gradient=False)
    roi.fit(use_gradient=False)


    roi.print_summary(galactic=True)

    roi.localize(which=name, update=True)
    roi.fit(use_gradient=False)

    roi.print_summary(galactic=True)

    print roi

    path=os.path.abspath(sys.argv[0])
    print """
Code to recreate point source:

    # Parameters came from %s
    model=PowerLaw(norm=%g, index=%g, e0=%g)
    skydir=SkyDir(%s,%s,SkyDir.GALACTIC)
    ps=PointSource(name="%s", model=model, skydir=skydir)
    roi.add_source(ps)
""" % (path,model['norm'],model['index'],model.e0, ps.skydir.l(),ps.skydir.b(),name)

    roi.save('roi.dat')

    roi.plot_tsmap(filename='residual_tsmap.pdf', **tsmap_kwargs)
