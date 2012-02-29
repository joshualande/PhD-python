

def plot_gtlike_cutoff_test(cutoff_results, sed_results, filename=None, title=None, 
                            model_0_kwargs=dict(color='red', zorder=0),
                            model_1_kwargs=dict(color='blue', zorder=0),
                            sed_kwargs=dict(),
                            plot_kwargs=dict(),
                           ):
    """ Plots the cutoff test performed of a spectrum using the function
        gtlike_test_cutoff.

        Input:
            cutoff_dict: created by gtlike_test_cutoff
            sed_dict: created by LandeSED.todict(). Can also be a yaml
              file created by LandeSED.save().

            model_0_kwargs: kwargs for model_0's plot 
            model_1_kwargs: kwargs for model_0's plot 
            sed_kwargs: kwargs to pass into LandeSED
              E.G. flux_units, flux_units, figsize, ...
    """
    sed=LandeSED(sed_results,**sed_kwargs)

    axes=sed.plot(plot_spectral_fit=False, **plot_kwargs)
    axes.autoscale(enable=False)

    model_0 = LandeSED.dict_to_spectrum(cutoff_results['model_0'])
    model_1 = LandeSED.dict_to_spectrum(cutoff_results['model_1'])
    sed.plot_spectrum(model_0, **model_0_kwargs)
    sed.plot_spectrum(model_1, **model_1_kwargs)

    if title is None:
        axes.set_title('Gtlike Cutoff test for %s' % sed.name)
    else:
        axes.set_title(title)


    if filename is not None: P.savefig(filename)
    return axes
