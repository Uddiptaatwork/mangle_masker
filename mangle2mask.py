#!/usr/bin/python3
#
# DESCRIPTION:
#  Converts polygon mangle masks (ply) into a binary mask in fits format.
#
# USAGE:
#  ra1       : R.A. center of the mask [deg]
#  dec1      : Dec. center of the mask [deg]
#  fov_x     : field of view [deg]
#  fov_y     : field of view [deg]
#  pixSize   : resolution of the pixels of the output mask [deg]
#  mangle_fn : name of the input mangle file
#  output_fn : name of the output mask fits file
#
#  python mangle2mask (ra1, dec1, fov_x, fov_y, pixSize, AMICO_mask_file)
#
# AUTHORS:
#  Uddipta Bhardwaj (ITA, Heidelberg 2020)
#  Matteo Maturi (ITA, Heideberg 2020)

# Import modulus
import sys
import math
import numpy as np
import pymangle
from astropy.io import fits
from astropy import wcs
    
def mask_creator (ra_ctr, dec_ctr, fov_x, fov_y, pixSize, mangle_fn):  

    # Define the number of pixels of the output mask (AMICO style)
    nx = int(math.floor(fov_x/pixSize + 1))
    ny = int(math.floor(fov_y/pixSize + 1))

    # Set WCS system
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = [int(nx/2.0),int(ny/2.0)]
    w.wcs.crval = [ra_ctr,dec_ctr]
    w.wcs.cdelt = np.array([pixSize, pixSize])
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    
    # Read the mangle polygons describing regions to be masked using pymangle
    mask = pymangle.Mangle(mangle_fn)

    # Loop over the mask image
    m = np.arange(nx*ny).reshape(ny,nx)
    for i in range(0, nx):
        for j in range(0, ny):
            # Transform to celestial coordinates
            this_ra, this_dec = w.wcs_pix2world(i, j, 0)
            m[j,i] = int(1 - mask.contains(this_ra,this_dec) + 0.5)

    return m, w.to_header()
    
#
# RUN INSTRUCTIONS
#

# Read parameters and sanity check
total = len(sys.argv)
cmdargs = str(sys.argv)
if total != 7+1:
    print ("Error: wrong number of arguments")
    print ('  usage: mangle2mask <ra> <dec> <fov_x> <fov_y> <pixSize> <mangle_fn> <output_fn>')
    sys.exit(1)
    
ra = float(sys.argv[1])
dec = float(sys.argv[2]) 
fov_x = float(sys.argv[3])
fov_y = float(sys.argv[4])
pixSize = float(sys.argv[5])
mangle_fn = str(sys.argv[6])
output_fn = str(sys.argv[7])
print ("  executing: %s " % cmdargs)

# Execute
m, h = mask_creator (ra, dec, fov_x, fov_y, pixSize, mangle_fn)

# Write into a binary fits file            
hdu = fits.PrimaryHDU (m, header=h)
hdu.writeto (output_fn)


print ('Done!')        
sys.exit(0)

