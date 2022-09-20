#%%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


# plt.style.use('seaborn-talk')
# plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams.update({'font.sans-serif':'Myriad Pro'})
plt.rcParams.update({'xtick.labelsize': 16, 
                     'ytick.labelsize': 16,
                     'axes.titlesize': 22,
                     'figure.titlesize': 16,
                     'axes.labelsize': 16,
                     'axes.labelsize': 16,
                     'legend.fontsize': 14,
                     'legend.title_fontsize': 14,
                     'figure.facecolor':(240/255,240/255,240/255),
                     'savefig.facecolor':(240/255,240/255,240/255)})

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
# Define vent geometry and environmental parameters (for use later in for loops)
vent_vel = [50, 70, 100, 150, 300] # m/s
vent_rad = [20, 30, 75, 127, 303] # m
wind_speed = range(0,55,5) # m/s
wind_speed_str = [format(x, '02d') for x in wind_speed]
lats = ['tropical', 'polar', 'mid_lat']
lats_output = ['Tropical', 'Polar', 'Mid-Lat']
bulk_density = 4.37 #kg/m^3

input_dir = '/Users/tylerpaladino/Documents/ISU/Thesis/ATHAM_wind/ATHAM_output/v7_stability_calc/'

# Read in all data as multiple dataframes in a python dictionary
dataframes_dict = {}
for lat in lats:
    for vent_idx,vent in enumerate(vent_rad):
        if vent == 127:
            vent = '127_5'
        # Create filename
        fn = input_dir + lat + '_' + str(vent) + 'm.txt'

        # Read in data directly into dict
        dataframes_dict[lat + '_' + str(vent)] = pd.read_csv(fn) 

        # Replace 'ms' string with nothing in each dataframe, then convert resulting strings to integers
        dataframes_dict[lat + '_' + str(vent)] = dataframes_dict[lat + '_' + str(vent)].replace('ms','',regex=True)
        dataframes_dict[lat + '_' + str(vent)]['Wind Speed (m/s)'] = dataframes_dict[lat + '_' + str(vent)]['Wind Speed (m/s)'].astype(int)
        dataframes_dict[lat + '_' + str(vent)]['Vent speed (m/s)'] = dataframes_dict[lat + '_' + str(vent)]['Vent speed (m/s)'].astype(int)
        dataframes_dict[lat + '_' + str(vent)]['stability mean'] = dataframes_dict[lat + '_' + str(vent)]['stability mean'] * 100

        dataframes_dict[lat + '_' + str(vent)]['Atmosphere'] = [lat] * len(dataframes_dict[lat + '_' + str(vent)])
        
        dataframes_dict[lat + '_' + str(vent)]['MER'] = dataframes_dict[lat + '_' + str(vent)]['Vent speed (m/s)'] * bulk_density * np.pi * vent_rad[vent_idx]**2
#%% Stability vs. Wind Speed
for lat_idx, lat in enumerate(lats):
# Create Figure and axes. Create 4x1 subplot grid
    fig1, ax1 = plt.subplots(len(vent_rad),1,figsize=(10,8))
    # Turn autolayout off since we're doing some sick, custom shit
    fig1.autolayout = False
    for idx, vent in enumerate(vent_rad):
        if vent == 127:
            vent = '127_5'
        vent_speed_group = dataframes_dict[lat + '_' + str(vent)].groupby('Vent speed (m/s)') # group by vent speed

        colors_vent_speed = plt.cm.viridis(np.linspace(0,1,num=vent_speed_group.ngroups)) # get colors for each group

        for count, value in enumerate(vent_vel):
            vent_speed_group.get_group(value).sort_values(by=['Wind Speed (m/s)']).plot(
            kind='line',
            style = '-o', 
            x='Wind Speed (m/s)',
            y='stability mean', 
            legend = False,  
            grid = True, 
            color=colors_vent_speed[count],
            ax=ax1[idx])

            ax1[idx].errorbar(vent_speed_group.get_group(value)['Wind Speed (m/s)'],
                vent_speed_group.get_group(value)['stability mean'],
                yerr=vent_speed_group.get_group(value)['stability SD']*100,
                fmt = 'none',
                ecolor = 'black',
                barsabove=True)
            ax1[idx].set_ylim(10,100)
    
        # Add ticks and labels based on wind_speed vector
        ax1[idx].set_xticks(wind_speed)
        ax1[idx].minorticks_on()
        # Add grid minor gird lines with a dashed, semi transparent appearance. 
        ax1[idx].grid(which='minor', linestyle='--', linewidth='0.25', color='grey', alpha=0.5)

        if vent == '127_5':
            text_str = '127.5' + ' m Vent Radius'
        else:
            text_str = str(vent) + ' m Vent Radius'

        ax1[idx].text(0.00, 1.04, text_str, transform=ax1[idx].transAxes, fontsize=16, horizontalalignment='left', weight=700, stretch = 'condensed')

        # Turn off xlabels and tick labels for all but bottom subplot
        if vent != vent_rad[-1]:
            ax1[idx].set(xlabel=None)
            ax1[idx].set_xticklabels([])

    # Add legend to bottom of figure
    fig1.legend(vent_vel, ncol = len(vent_vel), borderpad = 0.3, frameon=True, fancybox=True, title='Vent Speed (m/s)', loc=8)
    fig1.tight_layout()
    fig1.subplots_adjust(top=0.95,bottom = 0.18)   
    fig1.suptitle(lats_output[lat_idx] +' ' + 'Stability vs. Wind Speed',y=.99, fontsize = 22, weight=700, stretch = 'condensed')
    fig1.supylabel('Stability (%)', x=0)


    fig1.savefig(input_dir+lat+"_"+ "StabilityVsWind.pdf", format = 'pdf',bbox_inches=None, dpi=300)



#%% Stability vs. Vent Speed

for lat_idx, lat in enumerate(lats):
# Create Figure and axes. Create 4x1 subplot grid
    fig2, ax2 = plt.subplots(len(vent_rad),1,figsize=(10,8))
    # Turn autolayout off since we're doing some custom shit
    fig2.autolayout = False
    for idx,vent in enumerate(vent_rad):
        if vent == 127:
            vent = '127_5'
        wind_speed_group = dataframes_dict[lat + '_' + str(vent)].groupby('Wind Speed (m/s)') # group by wind speed

        colors_wind_speed = plt.cm.viridis(np.linspace(0,1,num=wind_speed_group.ngroups)) # get colors for each group

        for count, value in enumerate(wind_speed):
            wind_speed_group.get_group(value).sort_values(by=['Vent speed (m/s)']).plot(
                kind='line',
                style='-o',
                x='Vent speed (m/s)',
                y='stability mean', 
                grid = True,
                legend = False, 
                color=colors_wind_speed[count],
                ax=ax2[idx])
            
        ax2[idx].set_ylim(0,100)

        # Add ticks and labels based on vent velocities
        ax2[idx].set_xticks(vent_vel)
        ax2[idx].minorticks_on()
        # Add grid minor gird lines with a dashed, semi transparent appearance. 
        ax2[idx].grid(which='minor', linestyle='--', linewidth='0.25', color='grey', alpha=0.5)

        if vent == '127_5':
            text_str = '127.5' + ' m Vent Radius'
        else:
            text_str = str(vent) + ' m Vent Radius'

        ax2[idx].text(0.00, 1.04, text_str, transform=ax2[idx].transAxes, fontsize=16, horizontalalignment='left', weight=700, stretch = 'condensed')


        # Turn off xlabels and tick labels for all but bottom subplot
        if vent != vent_rad[-1]:
            ax2[idx].set(xlabel=None)
            ax2[idx].set_xticklabels([])

        # Add legend to bottom of figure
        fig2.legend(wind_speed_str, ncol = len(wind_speed), borderpad = 0.3, handletextpad = 0.3, columnspacing = .25, frameon=True, fancybox=True, title='Wind Speed (m/s)', loc=8)
        fig2.tight_layout()
        fig2.subplots_adjust(top=0.95,bottom = 0.18)   
        fig2.suptitle(lats_output[lat_idx] +' ' + 'Stability vs. Vent Velocity',y=.99, fontsize = 22, weight=700, stretch = 'condensed')
        fig2.supylabel('Stability (%)')
        
        fig2.savefig(input_dir+lat+"_"+ "StabilityVsVentVelocity.pdf", format = 'pdf',bbox_inches=None, dpi=300)

        # fig2.legend(wind_speed_str, ncols = len(wind_speed), center=True, loc = 'b', frame=True, fancybox=True, title='Wind Speed (m/s)')

print('Done')


# %%
trop_max_plume_heights = []
mid_lat_max_plume_heights = []
polar_max_plume_heights = []


flattened_df = pd.concat(dataframes_dict, ignore_index=True, axis = 0)


# for lat in lats:
#     for vent in vent_rad:
#         if vent == 127:
#             vent = '127_5'
#         if lat == 'tropical':
            
#             trop_max_plume_heights.extend(dataframes_dict[lat + '_' + str(vent)]['Max plume height (km)'])
#         elif lat == 'mid_lat':
#             mid_lat_max_plume_heights.extend(dataframes_dict[lat + '_' + str(vent)]['Max plume height (km)'])
#         elif lat == 'polar':
#             polar_max_plume_heights.extend(dataframes_dict[lat + '_' + str(vent)]['Max plume height (km)'])


# flattened_df = np.vstack((trop_max_plume_heights, mid_lat_max_plume_heights, polar_max_plume_heights))
import seaborn as sns

fig3, ax3 = plt.subplots(figsize=(10,8))

sns.violinplot(x = 'Atmosphere', y = 'Max plume height (km)', data = flattened_df, order = ["tropical", "mid_lat", "polar"], bw=0.25 , ax = ax3)
ax3.set_axisbelow(True)
ax3.set_title('Max Plume Height vs. Atmosphere', fontsize = 22, weight = 700, stretch = 'condensed')
fig3.savefig(input_dir+"plume_heightvsAtmos.pdf", format = 'pdf',bbox_inches=None, dpi=300)


fig4, ax4 = plt.subplots(figsize=(10,8))

sns.violinplot(x = 'Atmosphere', y = 'stability mean', data = flattened_df, order = ["tropical", "mid_lat", "polar"], bw=0.25 , ax = ax4)
ax4.set_axisbelow(True)
ax4.set_title('Average Stability vs. Atmosphere', fontsize = 22, weight = 700, stretch = 'condensed')
fig4.savefig(input_dir+"stabilityvsAtmos.pdf", format = 'pdf',bbox_inches=None, dpi=300)

# fig5, ax5 = plt.subplots(figsize=(10,8))

# sns.violinplot(x = 'Wind Speed (m/s)', y = 'stability mean', data = flattened_df, bw=0.25 , ax = ax5)
# ax5.set_axisbelow(True)
# ax5.set_title('Average Stability vs. Atmosphere', fontsize = 22, weight = 700, stretch = 'condensed')


print("hello")


# %%
for lat in lats:
    # Create Figure and axes. Create 4x1 subplot grid
    fig5, ax5 = plt.subplots(len(vent_rad),1,figsize=(10,8))
    # Turn autolayout off since we're doing some custom shit
    fig5.autolayout = False
    for idx, vent in enumerate(vent_rad):
        if vent == 127:
            vent = '127_5'

        
        vent_df = dataframes_dict[lat + '_' + str(vent)]
        
        # vent_df.plot(kind='scatter',x='MER',y='Wind Speed (m/s)', c='stability mean', s=100, logx=True, colorbar = True, cmap='viridis', ax=ax5[idx])
        sc = ax5[idx].scatter(vent_df['MER'], vent_df['Wind Speed (m/s)'], c=vent_df['stability mean'], s=100, cmap='viridis')
        ax5[idx].set_xscale('log')
        ax5[idx].set(ylabel=None)
        ax5[idx].set_yticks(np.arange(min(wind_speed),max(wind_speed)+5,20))
        ax5[idx].set_xlim([10**5,10**9])
        ax5[idx].grid(which='minor', linestyle='--', linewidth='0.25', color='grey', alpha=0.5)
        ax5[idx].set_axisbelow(True)

        plt.colorbar(sc, label = 'Stability (%)',location='right', aspect = 10, ax=ax5[idx])

        # ax5[idx].colorbar(label='Stability (%)',location='right')


        if vent != vent_rad[-1]:
            ax5[idx].set(xlabel=None)
            ax5[idx].set_xticklabels([])
    fig5.tight_layout()
    fig5.supylabel('Wind Speed (m/s)',x=0)
    fig5.subplots_adjust(top=0.95,bottom = 0.18)   
    fig5.suptitle(lats_output[lat_idx] +' ' + 'Wind Speed vs. Mass Eruption Rate',y=.99, fontsize = 22, weight=700, stretch = 'condensed')

# %%
from matplotlib.colors import Normalize
import matplotlib.cm as cm
cmap = cm.get_cmap('viridis')
normalizer = Normalize(0,100)
im = cm.ScalarMappable(norm=normalizer, cmap=cmap)
for lat in lats:
    
    # Create Figure and axes. Create 4x1 subplot grid
    fig6, ax6 = plt.subplots(len(vent_rad),1,figsize=(10,8))
    # Turn autolayout off since we're doing some custom shit
    fig6.autolayout = False
    for idx, vent in enumerate(vent_rad):
        if vent == 127:
            vent = '127_5'

        
        vent_df = dataframes_dict[lat + '_' + str(vent)]
        
        # vent_df.plot(kind='scatter',x='MER',y='Wind Speed (m/s)', c='stability mean', s=100, logx=True, colorbar = True, cmap='viridis', ax=ax5[idx])
        ax6[idx].scatter(vent_df['MER'], vent_df['Wind Speed (m/s)'], c=vent_df['stability mean'], s=100, cmap=cmap, norm=normalizer)
        ax6[idx].set_xscale('log')
        ax6[idx].set(ylabel=None)
        ax6[idx].set_yticks(np.arange(min(wind_speed),max(wind_speed)+5,20))
        ax6[idx].set_xlim([10**5,10**9])
        ax6[idx].grid(which='minor', linestyle='--', linewidth='0.25', color='grey', alpha=0.5)
        ax6[idx].set_axisbelow(True)

        # plt.colorbar(label = 'Stability (%)',location='right', ax=ax5[idx])

        # ax5[idx].colorbar(label='Stability (%)',location='right')


        if vent != vent_rad[-1]:
            ax5[idx].set(xlabel=None)
            ax5[idx].set_xticklabels([])
    fig6.tight_layout()
    fig6.supylabel('Wind Speed (m/s)',x=0)
    fig6.subplots_adjust(top=0.95,bottom = 0.18)   
    fig6.suptitle(lats_output[lat_idx] +' ' + 'Wind Speed vs. Mass Eruption Rate',y=.99, fontsize = 22, weight=700, stretch = 'condensed')
    fig6.colorbar(im,ax=ax6.ravel().tolist())

print('hi')
# %%

for lat_idx, lat in enumerate(lats):
# Create Figure and axes. Create 4x1 subplot grid
    fig7, ax7 = plt.subplots(len(vent_rad),1,figsize=(10,8))
    # Turn autolayout off since we're doing some sick, custom shit
    fig7.autolayout = False
    for idx, vent in enumerate(vent_rad):
        if vent == 127:
            vent = '127_5'
        vent_speed_group = dataframes_dict[lat + '_' + str(vent)].groupby('Vent speed (m/s)') # group by vent speed

        colors_vent_speed = plt.cm.viridis(np.linspace(0,1,num=vent_speed_group.ngroups)) # get colors for each group

        for count, value in enumerate(vent_vel):
            vent_speed_group.get_group(value).sort_values(by=['Wind Speed (m/s)']).plot(
            kind='line',
            style = '-o', 
            x='Wind Speed (m/s)',
            y='Max plume height (km)', 
            legend = False,  
            grid = True, 
            color=colors_vent_speed[count],
            ax=ax7[idx])
            
            # ax7[idx].set_ylim(0,50)
    
        # Add ticks and labels based on wind_speed vector
        ax7[idx].set_xticks(wind_speed)
        ax7[idx].minorticks_on()
        # Add grid minor gird lines with a dashed, semi transparent appearance. 
        ax7[idx].grid(which='minor', linestyle='--', linewidth='0.25', color='grey', alpha=0.5)

        if vent == '127_5':
            text_str = '127.5' + ' m Vent Radius'
        else:
            text_str = str(vent) + ' m Vent Radius'

        ax7[idx].text(0.00, 1.04, text_str, transform=ax7[idx].transAxes, fontsize=16, horizontalalignment='left', weight=700, stretch = 'condensed')

        # Turn off xlabels and tick labels for all but bottom subplot
        if vent != vent_rad[-1]:
            ax7[idx].set(xlabel=None)
            ax7[idx].set_xticklabels([])

    # Add legend to bottom of figure
    fig7.legend(vent_vel, ncol = len(vent_vel), borderpad = 0.3, frameon=True, fancybox=True, title='Vent Speed (m/s)', loc=8)
    fig7.tight_layout()
    fig7.subplots_adjust(top=0.95,bottom = 0.18)   
    fig7.suptitle(lats_output[lat_idx] +' ' + 'Max Plume Height vs. Wind Speed',y=.99, fontsize = 22, weight=700, stretch = 'condensed')
    fig7.supylabel('Height (km)', x=0)


    fig7.savefig(input_dir+lat+"_"+ "MaxPlumeHeightVsWind.pdf", format = 'pdf',bbox_inches=None, dpi=300)



# %% Flat runs
fn = '/Users/tylerpaladino/Documents/ISU/Thesis/ATHAM_wind/ATHAM_output/v7_stability_calc/tropical_flat_75m.txt'
df = pd.read_csv(fn)
df=df.replace('ms','',regex=True)
df['Wind Speed (m/s)'] = df['Wind Speed (m/s)'].astype(int)
df['Vent speed (m/s)'] = df['Vent speed (m/s)'].astype(int)
df['stability mean'] = df['stability mean'] * 100

vent_speed_group = df.groupby('Vent speed (m/s)') # group by vent speed
fig8, ax8 = plt.subplots(figsize=(10,8))
colors_vent_speed = plt.cm.viridis(np.linspace(0,1,num=vent_speed_group.ngroups)) # get colors for each group

for count,value in enumerate([50,70,100]):
    vent_speed_group.get_group(value).sort_values(by=['Wind Speed (m/s)']).plot(
    kind='line',
    style = '-o', 
    x='Wind Speed (m/s)',
    y='stability mean', 
    legend = False,  
    grid = True, 
    color=colors_vent_speed[count],
    ax=ax8)
    

    ax8.errorbar(vent_speed_group.get_group(value)['Wind Speed (m/s)'],
        vent_speed_group.get_group(value)['stability mean'],
        yerr=vent_speed_group.get_group(value)['stability SD']*100,
        fmt = 'none',
        ecolor = colors_vent_speed[count],
        capsize=5,
        barsabove=True)

ax8.set_ylim(10,100)
# Add ticks and labels based on wind_speed vector
ax8.set_xticks(wind_speed)
ax8.minorticks_on()
# Add grid minor gird lines with a dashed, semi transparent appearance. 
ax8.grid(which='minor', linestyle='--', linewidth='0.25', color='grey', alpha=0.5)
fig8.legend([50,70,100], ncol = len([50,70,100]), borderpad = 0.3, frameon=True, fancybox=True, title='Vent Speed (m/s)', loc=8)

fig8.savefig(input_dir+ "75m_tropical_flat_StabilityVsWind.pdf", format = 'pdf',bbox_inches=None, dpi=300)