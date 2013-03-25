from pandas import DataFrame

def HDUtoPandas(hdu, key=None):
    """ Convert a pyfits hdu to a pandas data structure where the indices in the DataFrame
        are one of the columns in the fits file. """
    names = list(hdu.data.names)

    if key is not None: assert key in names
    
    # stip out vector columns (for now)
    names = [str(name) for name in names if len(hdu.data[name].shape)==1]
    df=DataFrame({n:hdu.data[n] for n in names}, columns=names)
        
    if key is not None: df=df.set_index([key])
    return df
