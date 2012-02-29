""" Old SED code.

    Author: Joshua Lande <joshualande@gmail.com>
"""
def plot_lat(file,ax,
             plot_kwargs = dict(color='black', marker='o', 
                                linestyle='none', capsize=0, markersize=4, 
                                label='LAT'),
             ul_kwargs = dict(color='black', marker='o', 
                              linestyle='none', markersize=0),
             xlabel='Energy (MeV)',
             ylabel=r'$\mathrm{E}^2$ dN/dE $(\mathrm{MeV}\ \mathrm{cm}^{-2}\ \mathrm{s}^{-1})$'
            ):
    """ Plot onto a matplotlib axes ax a file created by the SED object. 
        fmt = format of points

        X axes is assumed to be 
    """
    lat_lines=open(file).readlines()

    lat_lines=[line.strip() for line in lat_lines]
    lat_lines=[line for line in lat_lines if line != '' and line[0] != '#' ]

    # First, make sure header is good
    if lat_lines[0].split() != ['Lower_Energy', 'Upper_Energy', 'Energy','dN/dE','dN/dE_Err']:
        raise Exception("Unable to parse LAT header. Bad names: %s" % lat_lines[0])
    if lat_lines[1].split() != ['[MeV]', '[MeV]', '[MeV]','[ph/cm^2/s/MeV]', '[ph/cm^2/s/MeV]']:
        raise Exception("Unable to parse LAT header. Bad units: %s" % lat_lines[1])
    
    # load in data
    lat_lines = lat_lines[2:]
    lower_energy, upper_energy, energy, flux, flux_err = zip(*[line.split() if len(line.split())==5 else line.split()+[''] for line in lat_lines])
    energy = np.asarray(map(float,energy))

    ul = np.asarray([True if '<' in i else False for i in flux])
    flux = np.asarray([float(i.replace('<','')) for i in flux])

    flux_err = np.asarray([float(i) if i is not '' else 0 for i in flux_err])

    # plot data points which are not upper limits
    ax.errorbar(energy[~ul],(energy**2*flux)[~ul],
                yerr=(energy**2*flux_err)[~ul],
                **plot_kwargs
                )

    if sum(ul)>0:
        # Plot upper limits
        ax.errorbar(energy[ul],(energy**2*flux)[ul],
                    yerr=[ (0.4*energy**2*flux)[ul], np.zeros_like(energy[ul]) ],
                    lolims=True,
                    **ul_kwargs
                   )

    if xlabel is not None: ax.set_xlabel(xlabel)
    if ylabel is not None: ax.set_ylabel(ylabel)

    ax.set_xscale('log')
    ax.set_yscale('log')

def plot_hess(file,ax, 
              xlabel='Energy (MeV)',
              ylabel=r'$\mathrm{E}^2$ dN/dE $(\mathrm{MeV}\ \mathrm{cm}^{-2}\ \mathrm{s}^{-1})$',
              **kwargs
             ):
    """ Plot HESS data points taken
        from HESS Auxilary website pages. 
    """

    plot_kwargs=dict(capsize=0, marker='+', color='black', markersize=4, linestyle='none', label='H.E.S.S')
    plot_kwargs.update(kwargs)


    hess_lines=open(file).readlines()

    hess_lines=[line.strip() for line in hess_lines]
    hess_lines=[line for line in hess_lines if line != '']

    if hess_lines[0].split() != ['Energy','Flux','Flux','Error_low','Flux','Error_high']:
        raise Exception("Unable to parse hess header")
    if hess_lines[1].split() != ['[TeV]','[/TeV','cm^2','s]']:
        raise Exception("Unable to parse hess header")

    hess_lines = hess_lines[2:]
    energy, flux, flux_low, flux_high = map(np.asarray,zip(*[map(float,line.split()) for line in hess_lines]))

    # convert energy to MeV
    energy *= 1e6
    # convert TeV poitns to MeV cm^-2 s^-1
    # TeV^-1 cm^-2 s^-1 * (1TeV/1e6 MeV) * MeV**2 = MeV cm^-1 s^-1
    flux *= 1e-6
    flux_low *= 1e-6
    flux_high *= 1e-6

    # not sure if there are

    ax.errorbar(energy,energy**2*flux,
                yerr=[energy**2*flux_low,energy**2*flux_high],
                **plot_kwargs)

    if xlabel is not None: ax.set_xlabel(xlabel)
    if ylabel is not None: ax.set_ylabel(ylabel)

    ax.set_xscale('log')
    ax.set_yscale('log')

