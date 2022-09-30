import xarray as xr
import os
import time
from fnmatch import fnmatch 

def compress(infile,outfile):
    ds = xr.open_dataset(infile,mask_and_scale=False,decode_times=False)

    comp = dict(zlib=True, complevel=1)
    encoding = {var: comp for var in ds.data_vars}
    ds.to_netcdf(outfile, encoding=encoding)

    os.remove(infile)
    
    os.rename(outfile,infile)
    
root = '/Volumes/MATHAM_2/mid_lat'
pattern = "*_MOV.nc"

paths = []
for path, subdirs, files in os.walk(root):
    for name in files:
        if fnmatch(name, pattern):
            paths.append(os.path.join(path,name))
# tic = time.perf_counter()
# fn = '/Users/tylerpaladino/Downloads/mid-lat_75m_100ms_20ms/atham_netCDF_MOV.nc'

for idx, path in enumerate(paths):
    print("Progress: " + str(idx) +"/" + str(len(paths)))
    compress(path,path+'compress.nc')

