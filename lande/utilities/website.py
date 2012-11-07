""" Code to help making websites. """
import os


def t2t(lines,filename,flags='--style color.css --css-sugar'):
    """ create the HTML for a given t2t file. """

    # b
    if isinstance(lines,list) and isinstance(filename,str):
        website='\n'.join(lines)
    elif isinstance(lines,str) and isinstance(filename,list):
        lines,filename=filename,lines
        website='\n'.join(lines)
    elif isinstance(lines,str) and isinstance(filename,str):
        website=lines
    else:
        raise Exception("...")


    file=open(filename,'w')
    file.write(website)
    file.close()

    os.system('txt2tags --target html %s %s' % (flags,filename))

