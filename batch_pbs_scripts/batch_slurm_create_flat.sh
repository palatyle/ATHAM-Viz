#!/bin/bash
# Sets up and submits a batch of ATHAM runs to the PBS queue for the step function wind profile


# Declare parameters to loop through
declare -a lats=("tropical")
declare -a ventsizes=("22_5")
declare -a velocities=("50" "70" "100")
declare -a windspeeds=("0" "5" "10" "15" "20" "25" "30" "35" "40" "45" "50")

# Define constant flags
i_flag="-i /lfs/palatyle.isu/ATHAM-Viz/IO_ref "
f_flag=" -f atham "
a_flag="-a INPUT_atham_setup "
d_flag=" -d INPUT_dynamic_setup"

for lat in ${lats[@]}; do
    for ventsize in ${ventsizes[@]}; do
        for velocity in ${velocities[@]}; do
            for windspeed in ${windspeeds[@]}; do
                # Set name of pbs file based off current parameters
                pbs_name=${lat}_step${ventsize}m_${velocity}ms_${windspeed}ms.batch

                # Copy master pbs file
                cp ../PBS_files/slurmfile.batch $pbs_name

                # Edit -N flag in new pbs file
                N_var='5s/.*/''#SBATCH --job-name='${lat}_step${ventsize}m_${velocity}ms_${windspeed}ms'/'
                sed -i '' "$N_var" $pbs_name

                # Create string to edit run line

                # set output directory and create it 
                out_dir='/lfs/palatyle.isu/output/'${lat}_step${ventsize}m_${velocity}ms_${windspeed}ms
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
                run_line='34s|.*|mpirun -n 64 /lfs/palatyle.isu/ATHAM_IO/exec/atham '${i_flag}${o_flag}${f_flag}${a_flag}${p_flag}${v_flag}${d_flag}'|'
                sed -i '' "$run_line" $pbs_name


                # Put into pbs queue 
                sbatch $pbs_name

                echo $pbs_name
            done
        done
    done
done
