import os
from pathlib import Path
from tqdm import tqdm
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
    
scratch_path = Path('/scratch/palatyle/')
pattern = "*_MOV.nc"


nc_files = list(scratch_path.glob('tropical_step*/*_MOV.nc'))

for file in tqdm(nc_files):
    size = file.stat().st_size/1e9 # get size and convert from b -> Gb
    if size < 16.0:
        compress(file,file.parent/'atham_netCDF_MOV_compress.nc')
    else:   
        print(file)

