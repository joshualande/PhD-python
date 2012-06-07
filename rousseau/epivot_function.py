import math
import sys

def epivot(K,g,CKg,Cgg,E1,E2):
    #if len(sys.argv) != 7:
    #printUsage %s K g CKg Cgg E1 E2

    """Calculate the pivot energy from the parameters of a PowerLaw2 and
    elements of the error matrix. See Jean Ballet's memo:

    http://tinyurl.com/epivot for further details

    for futher details. This program implements Eqns 7 and 7bis.

    Required parameters:

    K   - Integral flux from PL2
    g   - Spectral index (defined in the sense of the ST, i.e. less than zero)
    CKg - Covarience between K and g
    Cgg - Varience of g
    E1  - Low energy bound of PL2 energy range
    E2  - High energy bound of PL2 energy range"""
    

    if g == -1:
        loge = (math.log(E1)+math.log(E2))/2-CKg/K/Cgg;
    else:
        epsilon=(E2/E1)**(1+g);
        loge=(math.log(E1)-epsilon*math.log(E2))/(1-epsilon)-1/(g+1)-CKg/K/Cgg;
        
    e=math.exp(loge);
    print e
    return e
