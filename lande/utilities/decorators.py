import types
import inspect

def modify_defaults(**kwargs):
    """ Function decorator which replaces the default values
        in the object whith the defaults specified by kwargs.
        Useful if you disagree about what a meaningful
        default is. """
    def decorator(func):
        spec=inspect.getargspec(func)
        names,defaults=list(spec.args),list(spec.defaults)
        # Defaults is generally shorter because not all parameters have a default.
        offset=len(defaults)-len(names)

        for k,v in kwargs.items():
            index=names.index(k)
            defaults[offset+index]=v

        if isinstance(func, types.MethodType):
            func.im_func.func_defaults=tuple(defaults)
        else:
            func.func_defaults=tuple(defaults)

        return func
    return decorator

def select_quiet(func):
    def new(self,*args,**kwargs):
        old_quiet=self.quiet
        if kwargs.has_key('quiet'):
            self.quiet=kwargs.pop('quiet')
        if kwargs.has_key('verbose'):
            self.quiet=not kwargs.pop('verbose')

        ret=func(self,*args,**kwargs)
        self.quiet=old_quiet
        return ret
    return new
    return decorator
