import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.cm import ScalarMappable

# Plotting parameters
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams.update({'font.sans-serif':'Myriad Pro'})
xticks_font = font_manager.FontProperties(stretch='condensed', size=14)
colors_vent_speed = ([0/255, 119/255, 187/255],[51/255, 187/255, 238/255],[0/255, 153/255, 136/255],[238/255, 119/255, 51/255],[204/255, 51/255, 17/255])
shapes_vent_speed = ['-o','-o','-s','-s','-^']
face_colors = [colors_vent_speed[0],'white',colors_vent_speed[2],'white',colors_vent_speed[4]]

# Transition height data directory
data_dir = Path('/Users/tylerpaladino/Documents/ISU/Thesis/ATHAM_wind/ATHAM-Viz/transition_height_out')

# Stability output directory 
stab_out = '/Users/tylerpaladino/Documents/ISU/Thesis/ATHAM_wind/ATHAM-Viz/v8_stability_calc/'
# stab_files = stab_out.glob('*m.txt')





trop_u_wind = pd.read_csv('/Users/tylerpaladino/Documents/ISU/Thesis/ATHAM_wind/ATHAM-Viz/u_profile_norm_trop.csv', 
                        header=None,
                        names=['u_wind', 'altitude'])
mid_lat_u_wind = pd.read_csv('/Users/tylerpaladino/Documents/ISU/Thesis/ATHAM_wind/ATHAM-Viz/u_profile_norm_mid_lat.csv', 
                        header=None,
                        names=['u_wind', 'altitude'])
polar_u_wind = pd.read_csv('/Users/tylerpaladino/Documents/ISU/Thesis/ATHAM_wind/ATHAM-Viz/u_profile_norm_polar.csv', 
                        header=None,
                        names=['u_wind', 'altitude'])

means_dict = {}
std_devs_dict = {}

vent_rad = [15, 22.5, 75, 135, 315] # m
vent_vel = [50, 70, 100, 150, 300] # m/s
lats = ['tropical', 'mid_lat', 'polar']
bulk_density = 4.37 #kg/m^3

# Read in transition height data for each atmosphere
tropical_df = pd.read_csv(data_dir / "tropical.csv")
mid_lat_df = pd.read_csv(data_dir / "mid_lat.csv")
polar_df = pd.read_csv(data_dir / "polar.csv")

total_df = pd.concat([tropical_df,mid_lat_df,polar_df])
total_df['MER (kg/s)'] = total_df['Vent speed (m/s)'] * bulk_density * np.pi * (total_df['Vent Radius (m)']**2)


fig1,ax1 = plt.subplots(1,len(vent_rad)+1, figsize=(10,8),sharey=True)

df_list = []
for lat in lats:
    for vent in vent_rad:
        # Convert the .5 case to string with underscore for reading file
        if vent == 22.5:
            vent = '22_5'
        fn = stab_out + lat + '_' + str(vent) + 'm.txt'
        df_temp = pd.read_csv(fn)
        # Replace 'ms' string with nothing in each dataframe, then convert resulting strings to integers
        df_temp = df_temp.replace('ms','',regex=True)
        df_temp['Wind Speed (m/s)'] = df_temp['Wind Speed (m/s)'].astype(int)
        df_temp['Vent speed (m/s)'] = df_temp['Vent speed (m/s)'].astype(int)
        df_temp['stability mean'] = df_temp['stability mean'] * 100
        
        # Add atmosphere column
        df_temp['Atmosphere'] = [lat] * len(df_temp)
        # Convert back to float for comparisons 
        if vent == '22_5':
            vent = 22.5
        # Add vent radius column
        df_temp['Vent Radius (m)'] = [float(vent)] * len(df_temp)
        
        # Add in Transition height column. Fill with nans for now
        df_temp['Transition Height (km)'] = np.nan
        
        # Loop through vent velocities
        for vspeed in vent_vel:
            try:
                # Get the transition height at specified vent radius, atmosphere, and vent speed. 
                trans_height = float(total_df[(total_df['Vent Radius (m)'] == vent) 
                                              & (total_df['Atmosphere'] == lat) 
                                              & (total_df['Vent speed (m/s)'] == vspeed)]['Transition Height (km)'])
            except:
                trans_height = np.nan
            # Set transition height of all values in stability dataframe to value found above at specific vent radius, atmosphere, and vent speed. 
            df_temp.loc[(df_temp['Vent Radius (m)'] == vent) 
                        & (df_temp['Atmosphere'] == lat) 
                        & (df_temp['Vent speed (m/s)'] == vspeed),
                        'Transition Height (km)'] = trans_height
        # append dataframe to list
        df_list.append(df_temp)
        
stab_df = pd.concat(df_list)

trans_height_max = np.max(total_df['Transition Height (km)'])
trans_height_min = np.min(total_df['Transition Height (km)'])
for rad_count, rad_value in enumerate(vent_rad):
    means = []
    std_devs = []
    vent_group = total_df[total_df['Vent Radius (m)'] == rad_value].groupby('Vent speed (m/s)')
    for vel_count, vel_value in enumerate(vent_vel):
        try:
            vent_group.get_group(vel_value).plot(kind='line',
                                            style=shapes_vent_speed[vel_count],
                                            x='Atmosphere',
                                            y='Transition Height (km)',
                                            legend=False,
                                            color = colors_vent_speed[vel_count],
                                            fillstyle='full',
                                            markerfacecolor=face_colors[vel_count],
                                            ax = ax1[rad_count])
            means.append(np.mean(vent_group.get_group(vel_value)['Transition Height (km)']))
            std_devs.append(np.std(vent_group.get_group(vel_value)['Transition Height (km)']))
        except:
            ax1[rad_count].plot(['tropical',
                                 'mid-lat',
                                 'polar'],
                                [0,0,0],
                                shapes_vent_speed[vel_count],
                                fillstyle='full',
                                markerfacecolor=face_colors[vel_count],
                                color=colors_vent_speed[vel_count])
            means.append(0)
            std_devs.append(0)

            print("Missing data for vent speed: " + str(vel_value) + " and vent radius: " + str(rad_value))
            continue
    ax1[rad_count].grid()
    ax1[rad_count].set_title('R = '+str(rad_value)+' m')
    ax1[rad_count].minorticks_on()
    ax1[rad_count].set_ylim(0,8)
    ax1[rad_count].tick_params(axis='both',which='both',direction='in',top=True, right=True)
    ax1[rad_count].tick_params(axis='x', rotation=45)
    ax1[rad_count].set(xlabel=None)


    means_dict[rad_value] = means
    std_devs_dict[rad_value] = std_devs
ax1[5].plot(abs(trop_u_wind['u_wind']),trop_u_wind['altitude'],label='Tropical')
ax1[5].plot(abs(mid_lat_u_wind['u_wind']),mid_lat_u_wind['altitude'],label='Mid-Lat')
ax1[5].plot(abs(polar_u_wind['u_wind']),polar_u_wind['altitude'],label='Polar')
ax1[5].grid()
ax1[5].legend(loc='upper right', prop={'stretch':'condensed'})
ax1[5].set_xlabel('Normalized Wind Speed (m/s)', stretch = 'condensed')
ax1[5].tick_params(axis='both',which='both',direction='in',top=True, right=True,labelright=True)


fig1.legend(vent_vel, ncol = len(vent_vel), frameon=True, fancybox=True, title='Vent Speed (m/s)', loc=8,prop={'stretch':'condensed'})
# fig1.tight_layout()
# fig1.subplots_adjust(top=0.95,bottom = 0.18) 
fig1.suptitle('Transition Height vs. Atmosphere Type', fontsize = 22, weight=700, stretch = 'condensed')
ax1[0].set_ylabel('Altitude (km)')
fig1.savefig(data_dir / 'Transition_height_summary.pdf', format='pdf',bbox_inches=None, dpi=300)

# fig1.savefig(input_dir+ "Atmos_comparison.pdf", format = 'pdf',bbox_inches=None, dpi=300)

# fig2,ax2 = plt.subplots()
# for vel_count, vel_value in enumerate(vent_vel):
#     vent_group = total_df.groupby('Vent speed (m/s)')
#     try:
#         vent_group.get_group(vel_value).plot(kind='scatter',
#                                     style='-o',
#                                     x='MER (kg/s)',
#                                     y='Transition Height (km)',
#                                     legend=False,
#                                     color = colors_vent_speed[vel_count],
#                                     ax = ax2)
#     except:
#         continue

# # ax2.plot(total_df['MER (kg/s)'],total_df['Transition Height (km)'],'-o')
# fig2.legend(vent_vel, ncol = len(vent_vel), borderpad = 0.3, frameon=True, fancybox=True, title='Vent Speed (m/s)', loc=8)
# # fig1.tight_layout()
# # fig1.subplots_adjust(top=0.95,bottom = 0.18) 
# fig2.suptitle('Transition Height vs. Mass Eruption Rate',y=.99, fontsize = 22, weight=700, stretch = 'condensed')
# fig2.supylabel('Transition Height (km)', x=0)

# fig3,ax3 = plt.subplots()
# legend_list = []
# for rad_count, rad_value in enumerate(vent_rad):
#     for vel_count, vel_value in enumerate(vent_vel):
#         plot, = ax3.plot(vent_vel,means_dict[rad_value],'-o',color=colors_vent_speed[rad_count],label=rad_value)
#         ax3.errorbar(vent_vel,means_dict[rad_value],yerr=std_devs_dict[rad_value],fmt='none')
#     legend_list.append(plot)
# # fig3.legend(vent_rad, ncol = len(vent_rad), frameon=True, fancybox=True, title='Vent Radius (m)', loc=8)
# fig3.legend(handles=legend_list,ncol = len(vent_vel), frameon=True, fancybox=True, title='Vent Radius (m)', loc=8)
# fig3.suptitle('Transition Height vs. Vent Radius',y=.99, fontsize = 22, weight=700, stretch = 'condensed')
# fig3.supylabel('Transition Height (km)', x=0)
# fig3.supxlabel('Vent Radius (m)', x=0)


atmos_group = stab_df[stab_df['Vent Radius (m)'] <= 75 ].groupby('Atmosphere')

fig2,ax2 = plt.subplots(1,3,sharey=True)
# cmap = plt.get_cmap('cividis')
cmap = plt.get_cmap('viridis')

norm = plt.Normalize(vent_vel[0], vent_vel[-1])

# atmos_group.get_group('tropical').plot.scatter(x='Transition Height (km)',
#                                        y='stability mean',
#                                        s='Vent Radius (m)',
#                                        c='Vent speed (m/s)',
#                                        colormap=cmap,
#                                        colorbar=False,
#                                        grid=True,
#                                        alpha=0.5,
#                                        ax=ax2[0])

# atmos_group.get_group('tropical').plot.scatter(x='Transition Height (km)',
#                                        y='stability mean',
#                                        s='Vent Radius (m)',
#                                        c='none',
#                                        colormap=cmap,
#                                        colorbar=False,
#                                        grid=True,
#                                        edgecolors='black',
#                                        ax=ax2[0])

stab_df_tropical = atmos_group.get_group('tropical')
order = np.argsort(-stab_df_tropical['Vent Radius (m)'])

ax2[0].scatter(x=np.take(stab_df_tropical['Transition Height (km)'],order),
               y=np.take(stab_df_tropical['stability mean'],order),
               s=np.take(stab_df_tropical['Vent Radius (m)'],order)*4,
               c=np.take(stab_df_tropical['Vent speed (m/s)'],order),
               cmap=cmap,
               edgecolors='black',
               alpha=0.8)

ax2[0].grid()
ax2[0].axis('tight')
ax2[0].set_xlabel('Transition Height (km)')
ax2[0].set_ylabel('Stability (%)')
ax2[0].set_title('Tropical')
ax2[0].set_xlim(1,5)
ax2[0].set_ylim(0,100)
ax2[0].tick_params(axis='both',which='both',direction='in',top=True, right=True)

# atmos_group.get_group('mid_lat').plot.scatter(x='Transition Height (km)',
#                                        y='stability mean',
#                                        s='Vent Radius (m)',
#                                        c='Vent speed (m/s)',
#                                        colormap=cmap,
#                                        colorbar=False,
#                                        grid=True,
#                                        alpha=0.5,
#                                        ax=ax2[1])

# atmos_group.get_group('mid_lat').plot.scatter(x='Transition Height (km)',
#                                        y='stability mean',
#                                        s='Vent Radius (m)',
#                                        c='none',
#                                        colormap=cmap,
#                                        colorbar=False,
#                                        grid=True,
#                                        edgecolors='black',
#                                        ax=ax2[1])

stab_df_mid_lat = atmos_group.get_group('mid_lat')
order = np.argsort(-stab_df_mid_lat['Vent Radius (m)'])

ax2[1].scatter(x=np.take(stab_df_mid_lat['Transition Height (km)'],order),
               y=np.take(stab_df_mid_lat['stability mean'],order),
               s=np.take(stab_df_mid_lat['Vent Radius (m)'],order)*4,
               c=np.take(stab_df_mid_lat['Vent speed (m/s)'],order),
               cmap=cmap,
               edgecolors='black',
               alpha=0.8)

# ax2[1].scatter(x=stab_df_polar['Transition Height (km)'],
#                y=stab_df_polar['stability mean'],
#                s=stab_df_polar['Vent Radius (m)']*4,
#                c=stab_df_polar['Vent speed (m/s)'],
#                cmap=cmap,
#                alpha=0.5)

# ax2[1].scatter(x=stab_df_polar['Transition Height (km)'],
#                y=stab_df_polar['stability mean'],
#                s=stab_df_polar['Vent Radius (m)']*4,
#                c='none',
#                cmap=cmap,
#                edgecolors='black'
#                )

ax2[1].grid()
ax2[1].axis('tight')

ax2[1].set_xlabel('Transition Height (km)')
ax2[1].set_title('Mid-lat')
ax2[1].set_xlim(1,5)
ax2[1].set_ylim(0,100)
ax2[1].tick_params(axis='both',which='both',direction='in',top=True, right=True)

# atmos_group.get_group('polar').plot.scatter(x='Transition Height (km)',
#                                        y='stability mean',
#                                        s='Vent Radius (m)',
#                                        c='Vent speed (m/s)',
#                                        colormap=cmap,
#                                        colorbar=False,
#                                        grid=True,
#                                        alpha=0.5,
#                                        ax=ax2[2])

# atmos_group.get_group('polar').plot.scatter(x='Transition Height (km)',
#                                        y='stability mean',
#                                        s='Vent Radius (m)',
#                                        c='none',
#                                        colormap=cmap,
#                                        colorbar=False,
#                                        grid=True,
#                                        edgecolors='black',
#                                        ax=ax2[2])
stab_df_polar = atmos_group.get_group('polar')
order = np.argsort(-stab_df_polar['Vent Radius (m)'])

ax2[2].scatter(x=np.take(stab_df_polar['Transition Height (km)'],order),
               y=np.take(stab_df_polar['stability mean'],order),
               s=np.take(stab_df_polar['Vent Radius (m)'],order)*4,
               c=np.take(stab_df_polar['Vent speed (m/s)'],order),
               cmap=cmap,
               edgecolors='black',
               alpha=0.8)


# ax2[2].scatter(x=stab_df_polar['Transition Height (km)'],
#                y=stab_df_polar['stability mean'],
#                s=stab_df_polar['Vent Radius (m)']*4,
#                c=stab_df_polar['Vent speed (m/s)'],
#                cmap=cmap,
#                alpha=0.5)

# ax2[2].scatter(x=stab_df_polar['Transition Height (km)'],
#                y=stab_df_polar['stability mean'],
#                s=stab_df_polar['Vent Radius (m)']*4,
#                c='none',
#                cmap=cmap,
#                edgecolors='black'
#                )

ax2[2].grid()
ax2[2].axis('tight')
ax2[2].set_xlabel('Transition Height (km)')
ax2[2].set_title('Polar')
ax2[2].set_xlim(1,5)
ax2[2].set_ylim(0,100)
ax2[2].tick_params(axis='both',which='both',direction='in',top=True, right=True)

sm = ScalarMappable(norm=norm, cmap=cmap)
fig2.colorbar(sm, alpha=0.8, label= 'Vent Velocity (m/s)' ,ax=ax2)
print("pause")