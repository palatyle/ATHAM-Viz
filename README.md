# ATHAM-Viz

## Scripts used in Paladino et al. (2023)
`ATHAM_viz_multi.m` -- Script to visualize batches of ATHAM simulations using `ATHAM_viz_ts.m`

`ATHAM_viz_ts.m` -- Script to visualize individual ATHAM simulation.

`Stability_stats.py` -- Python script to visualize stats from ATHAM viz scripts. 

`netcdf_compression.py` -- Python script to compress multiple ATHAM netCDF files at once.

`plane_dir_viz.m` -- Script to visualize the directionality output from `ATHAM_viz_ts.m`

## Geometry files
`coarse_geo.mat` -- .mat file containing geometry of the ATHAM grid for the coarse grid used for larger vent radii.

`fine_geo.mat` -- .mat file containing geometry of the ATHAM grid for the fine grid used for larger vent radii.

## Visualization output
`direction_out/` -- Directory containing the diretionality .mat files for each simulation

`v8_stability_calc` -- Directory containing the stability output from `ATHAM_viz_multi.m`