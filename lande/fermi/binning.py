from GtApp import GtApp

def gtbin_from_binfile(evfile, outfile, scfile, binfile):
    """ Run gtbin on ft1 file evfile to create outfile outfile.
        but make the outfile have a binning consistent with the reference
        binfile. """

    GtApp(evfile=ft1, 
                outfile=binfile,
                scfile=ft2,
                algorithm='CCUBE',
                nxpix=
                nypix=
                binsz=
                xref=
                yref=
                axisrot=
                proj=

