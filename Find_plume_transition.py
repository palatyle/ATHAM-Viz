from intersect import intersection
import xarray as xr
from pathlib import Path
import argparse
import sys
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
from matplotlib import font_manager


plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams.update({'font.sans-serif':'Myriad Pro'})
xticks_font = font_manager.FontProperties(stretch='condensed', size=14)

def parse():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--atmosphere_dir",
                        default=None,
                        help="Full path of directory containing ATHAM output files for a particular atmosphere",
                        required=True)
    
    inputs = parser.parse_args()
    return inputs

def main():
    inputs = parse()
    
    # Read in atmosphere directory path
    atmosphere_dir = Path(inputs.atmosphere_dir)


    trop_u_wind = pd.read_csv('/Users/tylerpaladino/Documents/ISU/Thesis/ATHAM_wind/ATHAM-Viz/u_profile_norm_trop.csv', 
                            header=None,
                            names=['u_wind', 'altitude'])
    mid_lat_u_wind = pd.read_csv('/Users/tylerpaladino/Documents/ISU/Thesis/ATHAM_wind/ATHAM-Viz/u_profile_norm_mid_lat.csv', 
                            header=None,
                            names=['u_wind', 'altitude'])
    polar_u_wind = pd.read_csv('/Users/tylerpaladino/Documents/ISU/Thesis/ATHAM_wind/ATHAM-Viz/u_profile_norm_polar.csv', 
                            header=None,
                            names=['u_wind', 'altitude'])

    # Get all 0 ms wind speed netcdf MOV files
    netcdf_files = atmosphere_dir.glob('*/*/*_0ms/*_MOV.nc')
    
    
    transition_height = []
    atmos = []
    vent_rad = []
    vent_vel = []
    
    fig,ax = plt.subplots(1,2,figsize=(5.98,5.48))
    for file in netcdf_files:
        # Open netcdf file
        ds = xr.open_dataset(file,mask_and_scale=True,decode_times=False)
        
        # Grab dimensions of x/y
        x_size = ds.sizes['x']
        y_size = ds.sizes['y']
        
        # Grab z vector and divide by 1000 to convert to km
        z_full = ds.z/1000
        
        # Read in vertical profile of bulk density at center of simulation
        # 0th index gives ambient atmospheric density while last index gives 
        # plume density  
        atmos_den_full = ds.density[0,:,int(y_size/2),int(x_size/2)].to_numpy()
        plume_den_full = ds.density[-1,:,int(y_size/2),int(x_size/2)].to_numpy()
        

            
        # # Shrink z vector and density vectors based on nans (volcanic edifice)
        # z = z_full[~np.isnan(plume_den_full)]
        # atmos_den = atmos_den_full[~np.isnan(atmos_den_full)]
        # plume_den = plume_den_full[~np.isnan(plume_den_full)]
                
        if all(np.isnan(plume_den_full)) == False:    
            # Find all intersection points between plume and atmosphere profiles
            xi, yi = intersection(plume_den_full,z_full,atmos_den_full,z_full)
            ax[0].plot(plume_den_full,z_full,'-',label='Plume Density')
            ax[0].plot(atmos_den_full,z_full,'-',label='Atmosphere Density')
            ax[0].plot(xi[0],yi[0],'r*',markersize=8,label='Transition Point')
            ax[0].set_xlabel('Bulk Density (kg/m^3)', stretch = 'condensed', fontsize=16)
            ax[0].set_ylabel('Altitude (km)', stretch = 'condensed', fontsize=16)
            ax[0].set_xlim(0,1.8)
            ax[0].set_ylim(0,12)
            ax[0].legend(loc='upper right', prop={'stretch':'condensed','size':14})
            ax[0].grid()

            
            if atmosphere_dir.parts[3] == 'tropical':
                ax[1].plot(abs(trop_u_wind['u_wind']),trop_u_wind['altitude'],label='Tropical')
                ax[1].plot(abs(mid_lat_u_wind['u_wind']),mid_lat_u_wind['altitude'],'--',label='Mid-Lat', alpha=0.4)
                ax[1].plot(abs(polar_u_wind['u_wind']),polar_u_wind['altitude'],'--',label='Polar', alpha=0.4)
            elif atmosphere_dir.parts[3] == 'mid_lat':
                ax[1].plot(abs(trop_u_wind['u_wind']),trop_u_wind['altitude'],'--',label='Tropical', alpha=0.4)
                ax[1].plot(abs(mid_lat_u_wind['u_wind']),mid_lat_u_wind['altitude'],label='Mid-Lat')
                ax[1].plot(abs(polar_u_wind['u_wind']),polar_u_wind['altitude'],'--',label='Polar', alpha=0.4)
            elif atmosphere_dir.parts[3] == 'polar':
                ax[1].plot(abs(trop_u_wind['u_wind']),trop_u_wind['altitude'],'--',label='Tropical', alpha=0.4)
                ax[1].plot(abs(mid_lat_u_wind['u_wind']),mid_lat_u_wind['altitude'],'--',label='Mid-Lat', alpha=0.4)
                ax[1].plot(abs(polar_u_wind['u_wind']),polar_u_wind['altitude'],label='Polar')
            
            ax[1].grid()
            ax[1].set_xlabel('Normalized Wind Speed (m/s)', stretch = 'condensed', fontsize=16)
            ax[1].set_ylabel('Altitude (km)', stretch = 'condensed', fontsize=16)
            ax[1].yaxis.set_label_position("right")
            ax[1].yaxis.tick_right()
            ax[1].set_ylim(0,12)
            ax[1].legend(loc='upper right', prop={'stretch':'condensed','size':14})
            
            for tick in ax[0].get_xticklabels():
                tick.set_fontproperties(xticks_font)
            for tick in ax[1].get_xticklabels():
                tick.set_fontproperties(xticks_font)
            for tick in ax[0].get_yticklabels():
                tick.set_fontproperties(xticks_font)
            for tick in ax[1].get_yticklabels():
                tick.set_fontproperties(xticks_font)
            fig.suptitle(f'Transition Point & Wind\n Vent Radius = {file.parts[4][:-1]} m | Vent Velocity = {file.parts[5][:-2]} m/s', fontsize = 18, weight=700, stretch = 'condensed')

            fig.savefig((atmosphere_dir / 'density_plots' / file.parts[6]).with_suffix('.pdf'), format='pdf',bbox_inches=None, dpi=300)
            ax[0].clear()
            ax[1].clear()
        else:
            yi = [np.nan]
        
        # First index of yi is the very first cross over, i.e. transition between 
        # jet-thrust and convective rise. 
        transition_height.append(yi[0])
        
        # Atmosphere type is 3rd index
        atmos.append(file.parts[3])
        # Grab vent radius from file object. 4th index is vent rad directory name. 
        if file.parts[4] == '22_5m':
            vent_rad.append(22.5)
        else:
            # regex line is to remove non numeric characters
            vent_rad.append(int(re.sub("[^0-9]","",file.parts[4])))
            
        # Grab vent velocity from file object. 5th index is vent velocity directory name. 
        vent_vel.append(int(re.sub("[^0-9]","",file.parts[5])))

        # Append vent + atmos + transition height information to dataframe
        # df = df.append({'Vent Radius (m)': vent_rad ,
        #                 'Vent Velocity (m/s)': vent_vel,
        #                 'Atmosphere': atmos,
        #                 'Transition Height (km)': transition_height})

    df = pd.DataFrame({'Vent Radius (m)': vent_rad,
                        'Vent Velocity (m/s)': vent_vel,
                        'Atmosphere': atmos,
                        'Transition Height (km)': transition_height})

    df.to_csv(atmosphere_dir / 'out.csv')
    sys.exit(0)
    
if __name__ == '__main__':
    main()