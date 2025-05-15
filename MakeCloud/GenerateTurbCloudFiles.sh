#!/bin/bash

# This script runs the MakeCloud program to generate the parameter files for Gizmo
# It takes in the following arguments:
# 1. the file with the relavents masses and virial alpha values
# 2. Where to save the parameter files

# Define the radius of the cloud
radius=1.5 # pc
mass=600 # solar masses

file=$1
savePath=$2

# define the path for glass and turbulence files
glass_path=/work2/08770/kp32595/stampede3/StarFormProj/MakeCloud/glass_orig.npy
turb_path=/work2/08770/kp32595/stampede3/StarFormProj/MakeCloud/turbulence # folder where to store the turbulence files
# Define the number of particles
N=3.3e4
>Gizmo_jobs_file
module rm python2/2.7.15
source activate Cloud
# NewTurbSeed=$((42 + $RANDOM % 30))
# run the MakeCloud program with awk
awk -F" " '{print "python MakeCloud.py --R='$radius' --M='$mass' --N='$N' --turb_sol=0.5 --bturb=0.01 --boxsize=15 --alpha_turb="$1" --turb_seed="$2" --makebox  --glass_path='$glass_path' --turb_path='$turb_path'}' $file | bash
# python MakeCloud.py --R=1.5 --M=600 --N=3.3e4 --turb_sol=0.5 --bturb=0.01 --boxsize=10 --alpha_turb=2 --glass_path=/work2/08770/kp32595/stampede2/StarFormProj/MakeCloud/glass_orig.npy --turb_path=/work2/08770/kp32595/stampede2/StarFormProj/MakeCloud/turbulence
# find all the files that were just created

files=$(find . -name "*R${radius}*hdf5*")
# move the files to the save path
for file in $files
do
    # check if the file is a text file
    if [[ $file == *.txt ]]
    then
        # grep the name of the output directory
        dir=$(grep 'OutputDir' $file | awk -F" " '{print $2}')
        # check if the directory exists
        mkdir -p $dir
        # grep the mass from the file name
        mass=$(echo $file | awk -F"_" '{print $2}')
        alpha=$(echo $file | awk -F"_" '{print $6}')
        # grep the turb seed from the file name
        turbSeed=$(echo $file | awk -F"_" '{print $12}')
        # remove ./ from the file name
        tmp=$(echo $file | awk -F"./" '{print $2}')
        # write the GIZMO command to the job file
        echo "ibrun -n 2 ./GIZMO ./VP012_ParamFiles/$tmp 2 1>GizmoLogs/RunLogs/gizmoC_Seed$turbSeed.out 2>GizmoLogs/RunLogs/gizmoC_Seed$turbSeed.err & wait" >> Gizmo_jobs_file
        # grep IC file
        IC=$(grep 'InitCondFile' $file | awk -F" " '{print $2}')
        # add savepath to InitCondFile
        sed -i "s|$IC|./VP012_ParamFiles/$IC|g" $file

        # # replace the random number
        sed -i "s|TurbDrive_RandomNumberSeed          42|TurbDrive_RandomNumberSeed          $NewTurbSeed|g" $file
        # # find the time between snapshots
        # timeBetween=$(grep 'TimeBetSnapshot' $file | awk -F" " '{print $2}')
        # exp=$(($(echo $timeBetween | cut -d 'e' -f2)*1))
        # # increase by a factor of 10
        # exp=$(($exp+1))
        # exp=$(printf e%03d $exp)
        # # recombine the number
        # newTimeBet=$(echo $timeBetween | cut -d 'e' -f1)$exp
        # # replace the time between snapshots
        # sed -i "s|$timeBetween|$newTimeBet|g" $file

        # find the crossing time
        tCross=$(grep 'TurbDrive_CoherenceTime' $file | awk -F" " '{print $2}')
        # Correct by a factor of 2
        tCross=$(bc -l <<< "$tCross * 2")
        # define the new end time
        newEndTime=$(bc -l <<< "$tCross * 10")
        # find max time
        timeMax=$(grep 'TimeMax' $file | awk -F" " '{print $2}')
        # replace the end time
        sed -i "s|$timeMax|0$newEndTime|g" $file
    fi
    mv -vn $file $savePath
done
