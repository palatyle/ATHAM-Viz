#!/bin/bash
# Sets up and submits a batch of ATHAM runs to the PBS queue for the uniform flat wind profile

# Declare parameters to loop through
declare -a lats=("tropical")
declare -a ventsizes=("30")
declare -a velocities=("50" "70" "100" "150" "300")
declare -a windspeeds=("0" "5" "10" "15" "20" "25" "30" "35" "40" "45" "50")

# Define constant flags
i_flag="-i IO_ref "
f_flag=" -f atham "
a_flag="-a INPUT_atham_setup_small_grid "
d_flag=" -d INPUT_dynamic_setup"

for lat in ${lats[@]}; do
    for ventsize in ${ventsizes[@]}; do
        for velocity in ${velocities[@]}; do
            for windspeed in ${windspeeds[@]}; do
                # Set name of pbs file based off current parameters
                pbs_name=${lat}_flat${ventsize}m_${velocity}ms_${windspeed}ms.pbs

                # Copy master pbs file
                cp pbsfile.pbs $pbs_name

                # Edit -N flag in new pbs file
                N_var='8s/.*/''#PBS -N '${lat}_flat${ventsize}m_${velocity}ms_${windspeed}ms'/'
                sed -i "$N_var" $pbs_name

                # Create string to edit run line

                # set output directory and create it 
                out_dir='/scratch/palatyle/'${lat}_flat${ventsize}m_${velocity}ms_${windspeed}ms
                mkdir $out_dir
                o_flag='-o '$out_dir
                
                # Set atmospheric profile flag based off which latitude loop we're in
                if [ "${lat}" = "tropical" ]; then
                    p_flag='-p INPUT_profile_'${windspeed}mps'_trop_flat '
                elif [ "${lat}" = "mid-lat" ]; then
                    p_flag='-p INPUT_profile_'${windspeed}mps'_mid_lat '
                elif [ "${lat}" = "polar" ]; then
                    p_flag='-p INPUT_profile_'${windspeed}mps'_polar '
                fi

                # Set volcanic input flag based off vent size and velocity
                v_flag='-v INPUT_volcano_'${ventsize}_${velocity}

                # Set run line 
                run_line='34s|.*|mpirun ./exec/atham_st '${i_flag}${o_flag}${f_flag}${a_flag}${p_flag}${v_flag}${d_flag}'|'
                sed -i "$run_line" $pbs_name


                # Put into pbs queue 
                qsub $pbs_name

                echo $pbs_name
            done
        done
    done
done
