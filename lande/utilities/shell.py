

def format_command(*args, **kwargs):
    r""" Create a string suitable for running a shell program command where 
        *args are the positional arguments for the command and
        **kwargs are the keyword arguments for the script

        For example:

            >>> print format_command('ls','-al', '--author')
            ls \
                -al \
                --author
            >>> print format_command('gtlike', evfile='ft1.fits')
            gtlike \
                evfile=ft1.fits

        If you need parameters with dashes, you can pass in a dictionary:

            >>> print format_command('du', '-h', {'--max-depth':3})
            du \
                -h \
                --max-depth=3

        This function is not (yet) very robust, but does what I need.
    """
    line_break = ' \\' # slash
    tab='    '
    sep = '\n'.join([line_break,tab])

    args=list(args)
    for i,v in enumerate(args):
        if isinstance(v,dict):
            kwargs.update(args.pop(i))

    if args < 1: raise Exception("Command name must be passed into script")

    return sep.join(map(str,args) + ['%s=%s' % (a,b) for a,b in kwargs.items()])
    

if __name__ == "__main__":
    import doctest
    doctest.testmod()
