import os
import time
from fnmatch import fnmatch

import xarray as xr


def compress(infile,outfile):
    """Compress netcdf file and delete original

    Parameters
    ----------
    infile : str
        Path to netcdf file
    outfile : str
        Path to output netcdf file

    Returns
    -------
    None
    """    
    # Open dataset using xarray
    ds = xr.open_dataset(infile,mask_and_scale=False,decode_times=False)

    # Compress
    comp = dict(zlib=True, complevel=1)
    encoding = {var: comp for var in ds.data_vars}
    ds.to_netcdf(outfile, encoding=encoding)

    # Delete original netcdf file
    os.remove(infile)
    
    # Rename compressed file to original name
    os.rename(outfile,infile)
    return None
    
root = '/Volumes/MATHAM_2/mid_lat'
pattern = "*_MOV.nc"

paths = []
# Get paths of all netcdf files in root directory
for path, subdirs, files in os.walk(root):
    for name in files:
        if fnmatch(name, pattern):
            paths.append(os.path.join(path,name))

# Loop through all netcdf files and compress
for idx, path in enumerate(paths):
    print("Progress: " + str(idx) +"/" + str(len(paths)))
    compress(path,path+'compress.nc')

