# ATHAM-Viz

## Scripts used in Paladino et al. (2024)
`ATHAM_viz_multi.m` -- Script to visualize batches of ATHAM simulations using `ATHAM_viz_ts.m`

`ATHAM_viz_ts.m` -- Script to visualize individual ATHAM simulation.

`Stability_stats.py` -- Python script to visualize stats from ATHAM viz scripts. Creats many of the figures in Paladino et al. (2024)

`netcdf_compression.py` -- Python script to compress multiple ATHAM netCDF files at once.

`plane_dir_viz.m` -- Script to visualize the directionality output from `ATHAM_viz_ts.m`

## Geometry files used in `plane_dir_viz.m`
`coarse_geo.mat` -- .mat file containing geometry of the ATHAM grid for the coarse grid used for larger vent radii.

`fine_geo.mat` -- .mat file containing geometry of the ATHAM grid for the fine grid used for smaller vent radii.

## Inputs
`IO_ref/` -- Input files for ATHAM incluidng atmospheric profiles and volcano inputs. 

`PBS_files/` -- Directory containing all PBS submission files used in simulations. These submission files reference inputs in `IO_ref/`

## Visualization output
`direction_out/` -- Directory containing the diretionality .mat files for each simulation

`v8_stability_calc/` -- Directory containing the stability output from `ATHAM_viz_multi.m`

## Helper shell scripts
`batch_pbs_scripts/batch_pbs_create.sh` -- Main shell script used to create ATHAM simulations for vent sizes of 75, 135, and 315 meters. Also submits to PBS queue.

`batch_pbs_scripts/batch_pbs_create_flat.sh` -- Shell script used to create ATHAM simulations for the flat(uniform) atmospheric profile. Also submits to PBS queue.

`batch_pbs_scripts/batch_pbs_create_step.sh` -- Shell script used to create ATHAM simulations for the step function atmospheric profile. Also submits to PBS queue.

`batch_pbs_scripts/batch_small_vent.sh` -- Shell script used to create ATHAM simulations for vent sizes of 15 and 22.5 meters. Also submits to PBS queue.

`batch_pbs_scripts/batch_small_vent_144.sh` -- Shell script used to create ATHAM simulations for vent sizes of 15 and 22.5 meters running on 144 cores. Also submits to PBS queue.

