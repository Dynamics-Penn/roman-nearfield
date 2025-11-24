#adapted by Robyn Sanderson from code written by Saurabh Jha and Javier Sanchez
import numpy as np
import healpy as hp
import astropy.coordinates
from astropy import units as u

def change_coord(m, coord):   
    """ Change coordinates of a HEALPIX mapc

    Parameters
    ----------
    m : map or array of maps
      map(s) to be rotated
    coord : sequence of two character
      First character is the coordinate system of m, second character
      is the coordinate system of the output map. As in HEALPIX, allowed
      coordinate systems are 'G' (galactic), 'E' (ecliptic) or 'C' (equatorial)

    Example
    -------
    The following rotate m from galactic to equatorial coordinates.
    Notice that m can contain both temperature and polarization.
    >>>> change_coord(m, ['G', 'C'])
    """
    # Basic HEALPix parameters
    npix = m.shape[-1]
    nside = hp.npix2nside(npix)
    ang = hp.pix2ang(nside, np.arange(npix))

    # Select the coordinate transformation
    rot = hp.Rotator(coord=reversed(coord))

    # Convert the coordinates
    new_ang = rot(*ang)
    new_pix = hp.ang2pix(nside, *new_ang)

    return m[..., new_pix]

 	# from https://stackoverflow.com/questions/44443498/how-to-convert-and-save-healpy-map-to-different-coordinate-system

def read_galaxies():
    """ reads list of all galaxies within 10 Mpc 
    (requires nbg.cat--- file with table of galaxies from Karachentsev et al)
      
    Returns
    -------
    nbg : pandas dataframe with table in it
    nbgs : list of skycoords for each galaxy 
    
    Example
    -------
    Plot all the galaxies on a predefined skyproj
    >>>> plot_galaxies(sp)
    """

    import pandas as pd
    #read and process the table
    nbg=pd.read_csv('nbg.cat',sep="|",usecols=range(1,13))
    ras=[]
    decs=[]
    hms_str = ['h','m','s']
    dms_str = ['d','m','s']
    for i in range(len(nbg)):
        ra_str=nbg['RA J2000  '][i].split(' ')
        dec_str = nbg['DEC J2000'][i].split(' ')
        ras.append("".join(x+y for x,y in zip(ra_str,hms_str)))
        decs.append("".join(x+y for x,y in zip(dec_str,dms_str)))
    nbgs= astropy.coordinates.SkyCoord(ras,decs)

    return nbg, nbgs

# Read maps with the "tiles" covered -- these are built from APT's outputs
HLWAS_wide = hp.read_map('HLWAS_WIDE_0328_APT994_1024.fits.gz')
HLWAS_medium = hp.read_map('HLWAS_MEDIUM_0328_APT994_1024.fits.gz')
maps = [HLWAS_wide, HLWAS_medium]

#combine maps
combined_map = np.zeros_like(HLWAS_wide)
for i, m in enumerate(maps):
    m[m==hp.UNSEEN] = 0.
    combined_map[m!=0] = (i+1)*1.0

#find galaxies in HLWAS fields
nbgs, nbg_coords = read_galaxies() #replace this with your own read function
nbgs_in_hlwas_north=[]
nbgs_in_hlwas_south=[]
for i,sat in enumerate(nbg_coords):
    ipix = hp.ang2pix(hp.get_nside(combined_map),sat.galactic.l.value, sat.galactic.b.value,lonlat=True)
    if not combined_map[ipix]==hp.UNSEEN:
        if sat.galactic.b.value > 0:
            nbgs_in_hlwas_north.append(i)
        else:
            nbgs_in_hlwas_south.append(i)

print(nbgs_in_hlwas_north)
print(nbgs_in_hlwas_south)
