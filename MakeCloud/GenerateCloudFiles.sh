#!/bin/bash

# This script runs the MakeCloud program to generate the parameter files for Gizmo
# It takes in the following arguments:
# 1. the file with the relavents masses and virial alpha values
# 2. Where to save the parameter files

# Define the radius of the cloud
radius=0.025392 # pc
file=$1
savePath=$2

# define the path for glass and turbulence files
glass_path=/path/to/MakeCloud/glass_orig.npy
turb_path=/path/to/MakeCloud/turbulence # folder where to store the turbulence files

# Define the number of particles
N=10000
>Gizmo_jobs_file
module rm python2/2.7.15
source activate Cloud
# run the MakeCloud program with awk
awk -F" " '{print "python MakeCloud.py --M="$2" --R='$radius' --N='$N' --alpha_turb=0 --bturb=0 --warmgas --bfixed=0 --spin=0 --nsnap=300 --tmax=1 --glass_path='$glass_path' --turb_path='$turb_path'"}' $file | bash

# find all the files that were just created
files=$(find . -name "*R${radius}*")
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
        # remove ./ from the file name
        tmp=$(echo $file | awk -F"./" '{print $2}')
        # write the GIZMO command to the job file
        echo "./GIZMO ./ParamFiles/$tmp 2 1>GizmoLogs/RunLogs/gizmoC_$mass.out 2>GizmoLogs/RunLogs/gizmoC_$mass.err & wait" >> Gizmo_jobs_file
        # grep IC file
        IC=$(grep 'InitCondFile' $file | awk -F" " '{print $2}')
        # add savepath to InitCondFile
        sed -i "s|$IC|./ParamFiles/$IC|g" $file
    fi
    mv -vn $file $savePath
done
