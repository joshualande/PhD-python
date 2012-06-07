import pyfits, pywcsgrid2
import matplotlib.pyplot as plt
import numpy as np
import pylab
from skymaps import SkyDir
from uw.utilities import region_writer
import pyregion
import pywcs


IMAGE_SIZE = 30

def angle_equatorial_from_galactic(center, angle_galactic):
    """ For a given position in the sky,
        and a given angle defined as a rotation east of galactic
        north, comute the angle defined as a rotation east of
        celestial north. """
    l,b=center.l(),center.b()
    ra,dec=center.ra(),center.dec()

    # define a galacitc projection coordiante system centered on our position in the sky
    wcs = pywcs.WCS(naxis=2)

    # n.b. (l,b) has an image coordiante of (0,0), which is convenient
    wcs.wcs.crpix = [0, 0]
    wcs.wcs.cdelt = np.array([-0.1, 0.1])
    wcs.wcs.crval = [l, b]
    wcs.wcs.ctype = ["GLON-ZEA", "GLAT-ZEA"]

    # compute the axes of increasing 'dec' and 'b'

    gal_axis = wcs.wcs_sky2pix(np.array([[l,b + 0.01]], float),1)[0]
    # find the image vector of increasing galactic latitude
    gal_axis /= np.sqrt(gal_axis.dot(gal_axis))
    
    cel_dir = SkyDir(ra,dec+0.01)
    # find the image vector of increasing equatorial latitude
    cel_axis = wcs.wcs_sky2pix(np.array([[cel_dir.l(),cel_dir.b()]], float),1)[0]
    cel_axis /= np.sqrt(cel_axis.dot(cel_axis))

    # compute the angle between the two axes

    # N.B. use compute angle of each axes using atan2 to get correct coordinate.
    # The difference of these two angles should be the absolute angle
    # with the correct minus sign.
    rotation = np.degrees(
        np.arctan2(gal_axis[1],gal_axis[0]) - np.arctan2(cel_axis[1],cel_axis[0])
    )

    # the new angle is just the old angle + the difference between axes

    # The + sign was determined emperically
    angle_equatorial = angle_galactic +  rotation
    return angle_equatorial 


def new_angle(center, angle):
    """ here, angel is measured with respect to galactic origin. """
    l,b = center.l(), center.b()
    ra,dec = center.ra(), center.dec()

    celn= center(ra, dec+0.1)
    celn_l = celn.l()
    celn_b = celn.b()

    

def gal_header(center, major_axes, galactic=True):

    values = [
            ["NAXIS",  2,          ],

            ["NAXIS1", IMAGE_SIZE,       ],
            ["NAXIS2", IMAGE_SIZE,       ],

            ["CTYPE1", 'GLON-ZEA' if galactic else 'RA---ZEA' ],
            ["CTYPE2", 'GLAT-ZEA' if galactic else 'DEC--ZEA' ],

            ["CRPIX1", IMAGE_SIZE/2. + 0.5,       ],
            ["CRPIX2", IMAGE_SIZE/2. + 0.5,       ],

            ["CRVAL1", center.l() if galactic else center.ra(),        ],
            ["CRVAL2", center.b() if galactic else center.dec(),       ],

            ["CDELT1", -3.*major_axes/IMAGE_SIZE,       ],
            ["CDELT2", 3.*major_axes/IMAGE_SIZE,        ],
    ]

    if galactic is False:
        values += [
            ['RADECSYS','FK5'],
            ['EQUINOX',2000],
        ]


    cards = [pyfits.Card(*i) for i in values]

    header=pyfits.Header(cards=cards)

    return header

equ_header = lambda *args: gal_header(*args, galactic=False)

def draw_gal(center, major_axes):
    ra,dec=center.ra(),center.dec()

    header=gal_header(center, major_axes)

    global gal_ax

    gal_ax=pywcsgrid2.subplot(121, header=header)
    gal_ax.set_ticklabel_type("absdeg")
    gal_ax.imshow(np.zeros((IMAGE_SIZE,IMAGE_SIZE)), origin="lower", cmap=plt.cm.gray_r)
    gal_ax.grid()

    # add floating equatorial axes on top of galactic axes
    gal_ax.axis["ra=%.2f" % ra] = gal_ax["fk5"].new_floating_axis(0, ra)
    gal_ax.axis["dec=%.2f" % dec] = gal_ax["fk5"].new_floating_axis(1, dec)
    gal_ax.add_compass(loc=1)

    return gal_ax


def draw_cel(center, major_axes):

    header = equ_header(center, major_axes)

    l,b=center.l(),center.b()
    cel_ax=pywcsgrid2.subplot(122, header=header)
    cel_ax.set_ticklabel_type("absdeg")
    cel_ax.imshow(np.zeros((IMAGE_SIZE,IMAGE_SIZE)), origin="lower", cmap=plt.cm.gray_r)
    cel_ax.grid()


    # add floating galacitc axes on top of equatorial axes
    cel_ax.axis["l=%.2f" % l] = cel_ax["gal"].new_floating_axis(0, l)
    cel_ax.axis["b=%.2f" % b] = cel_ax["gal"].new_floating_axis(1, b)
    cel_ax.add_compass(loc=1)

    return cel_ax

def draw_ellipse(ax, center, region_string):

    header = ax.projection._pywcs.to_header()

    cel_ax['fk5'].plot([center.ra()],[center.dec()],'r+', zorder=10, markersize=20)

    reg = pyregion.parse(region_string).as_imagecoord(header)
    patch_list, artist_list = reg.get_mpl_patches_texts()
    for p in patch_list:
        p.set_zorder(10)
        ax.add_patch(p)

def plot(center, major_axes, region_string):
    global gal_ax, cel_ax


    fig = pylab.figure(None,(15,10))
    fig.subplots_adjust(wspace=0.4)


    gal_ax=draw_gal(center, major_axes)
    cel_ax=draw_cel(center, major_axes)

    draw_ellipse(cel_ax, center, region_string)
    draw_ellipse(gal_ax, center, region_string)

    pylab.suptitle('(l,b)=(%.2f,%.2f)\n(ra,dec)=(%.2f,%.2f)' % \
             (center.l(),center.b(),center.ra(),center.dec()))

    pylab.savefig('ellipse.png')


def plot_galactic(center, major_axes, minor_axes, angle_galactic):
    angle_equatorial = angle_equatorial_from_galactic(center, angle_galactic)
    print 'equatorial angle',angle_equatorial 

    region_string="\n".join(
        ["fk5; ellipse(%s, %s, %s, %s, %s) # color=red" % \
         (center.ra(),center.dec(),minor_axes,major_axes,angle_equatorial),
         "galactic; ellipse(%s, %s, %s, %s, %s) # color=blue dash=1" % \
         (center.l(),center.b(),minor_axes,major_axes,angle_galactic)
        ])
    print region_string

    plot(center, major_axes, region_string)

def test_around_poles():

    for center in [
        # points around equatorial north pole

        SkyDir(160,-10,SkyDir.GALACTIC),
        SkyDir(80,-10,SkyDir.GALACTIC),
        SkyDir(80,80,SkyDir.GALACTIC),
        SkyDir(160,50,SkyDir.GALACTIC),

        # points around equatorial south pole

        SkyDir(280,-10,SkyDir.GALACTIC),
        SkyDir(320,-10,SkyDir.GALACTIC),
        SkyDir(320,-40,SkyDir.GALACTIC),
        SkyDir(280,-40,SkyDir.GALACTIC),
    ]:

        major_axes= .2
        minor_axes= 0.005
        angle_galactic = 10

        plot_galactic(center, major_axes, minor_axes, angle_galactic)
        pylab.show()


if __name__ == '__main__':

    test_around_poles()
