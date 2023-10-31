# Hybrid-Eddy-detection
This is the hyrbrid eddy detection project for "A Hybrid 3D Eddy Detection Technique Based on Sea Surface Height and Velocity Field" in EnvirVis 2023 by Vizlab at Rutgers University.

![plot](./method_comparison.jpg)

![plot](./more_dataset.jpg)


See our previous [feature tracking](https://github.com/VizlabRutgers/Feature_Tracking) project for a more general usage.

Our paper is accepted by EnvirVis 2023. If you use our code, please cite our paper:  

>Hua, Weiping, et al. "A Hybrid 3D Eddy Detection Technique Based on Sea Surface Height and Velocity Field." arXiv preprint arXiv:2305.08229 (2023).

We've used it in SciViz Contest 2020 with Red Sea data https://kaust-vislab.github.io/SciVis2020/.

Please contact Rutgers Vizlab or the author if you have any problem.
https://vizlab.rutgers.edu/ or huaweiping0@gmail.com

# Eddy Visualization

The visualization program of the eddy is here: https://github.com/VizlabRutgers/Hybrid-Eddy-Visualization-and-Comparison

# Installation
## With Singularity

We recommend you to use the Singularity to avoid installation and environments configuration (especially in HPC, which might not support OpenGL/VTK). If so, you could get the image by the command below and jump to the next section.
```
singularity pull --arch amd64 library://huaweiping/hybrid_eddy_env/v1.6:latest
```
## Without Singularity
This program is a C++ program developed in linux system. We strongly recommend you to compile the code in linux system. You may need to change the CMakelist file if you're using other platforms.

Validated environment:  
Ubuntu-18.04  
CMake-3.14  
gcc-7.5.0  
Hdf5-1.14.1  
NetCDF-4.3.1  
VTK-8.2.0  
OpenCV-3.4.13  

IDE (Optional):  
QT-5.14.2 

Installation Steps:
1. Install Cmake/CMake-gui.(not use the latest one as qt may not support that currently. Version 3.14 verified).
2. (optional)Install qt5.
3. Install Hdf5 and Netcdf first.
4. Install VTK library.  
    a)You may meet some problem because of the incompatible version of cmake or VTK library.  
    b)You need to change the library path in CMakeList file.  
    c)You may need to use apt-file command to find necessary file (for QT).  
    d)If you are working in QT, you may need to compile VTK with QT.  
6. Compile the source code by Cmake(or Do it in QT).
7. Make file.
8. You could see the excutable file "FT" in your folder. 

# Use our software 
## With Singularity
If you are using Singularity, you need to modify the path in the configuration file first and then you could run our program in one command:
```
singularity exec (--bind <The directory of your data and output path>) <The path to your Singularity image> /eddy_hybrid_build/FT <The path to your FreatureTrack.Conf file>
```

## Without Singularity
If you are not using singularity, please follow the instruction below:
1. Modify the configuration file (FeatureTrack.conf) to match the path of the dataset
2. Copy the FeatureTrack.Conf to folder where you compile the code to.
3. Run the executable file with configuration file. E.g.
```
weiping@Precision: .\FT .\FeatureTrack.Conf
```



# Configuration File
The program support 2D and 3D dataset with variables including temperature and salinity.
Pleases modify at least the first three parameters in the Configuration file to a proper path for the dataset and generated result. 
Here's the detail of the configuration file (FeatureTrack.conf) for our hybrid detection approach.

```
# This is the data file path for splitted nc file (each nc file includes only 1 frame). Ignored when STACKED_NC_DATA_PATH exist.
DATA_FILES_PATH:  /home/weiping/data/ft_changes/sixty_frame/source/

# This is the path for generated file after detection.    
GENERATED_FILES_PATH:  /home/weiping/data/ft_changes/FT_result/10/

# This is the data file path for merged nc file (only one nc file for whole dataset).
STACKED_NC_DATA_PATH:  /home/weiping/data/SciViz/SciVisContest2020/ensembles/0010/COMBINED_2011013100.nc

# This is the base name for the input and output data.
FILE_BASE_NAME: red_sea_

# This is the file extension for the input data. Support .nc/.nc4/.vtk/.hd5 file. Only verified by .nc file.
FILE_EXTENSION: .nc

# This is initial time frame.
INITIAL_TIME_STEP: 1

# This is final time frame.
FINAL_TIME_STEP: 60

# This is time frame increment.
TIME_STEP_INCREMENT: 1

# This is time frame precision.
TIME_STEP_PRECISION: 1

# Used for value-based detection approach (such as OW). Ignored in hybrid detection.
VARIABLE_NAMES: omega
THRESHOLD1: -1.0
STARTRADIUS: 3
DELTA_X_THRESHOLD: 0.01
DELTA_Y_THRESHOLD: 0.01
DELTA_Z_THRESHOLD: 0.01

# This is the volume filter threshold. Will filter out small objects below this volumn.
SMALLEST_OBJECT_VOLUME_TO_TRACK: 20

# Dimension information for the dataset.  
X_Dim: 500
Y_Dim: 500
Z_Dim: 50

# Start/end voxel coordinate.
X1_Dim: 499
Y1_Dim: 499
Z1_Dim: 49
X0_Dim: 0
Y0_Dim: 0
Z0_Dim: 0

# The name of variables in the dataset (Netcdf dataset in this example)
xCoord_Name: XC
yCoord_Name: YC
zCoord_Name: Z_MIT40
Temp_Name: TEMP
Salinity_Name: SALT
Velocity_U_Name: U
Velocity_V_Name: V
Welocity_W_Name: V
SSH_Name: ETA
```
Except the configuration file. You may also need to change the variable name in the code depending on different dataset. The names of the value are located in Line 2575-2588 in mainFeatureTracki.cpp (if you're using 3D volumn data in stacked data file). This problem will be replaced by variable name in configuration file.


# Output
The program will automatically generate the ./Seperated Structure, ./Seperated Structure/clockwise and ./Seperated Structure/counterclockwise folder under the root path of the assigned GENERATED_FILES_PATH.

## root folder
.uocd file and BASE_FILE_NAME_clean.uocd file under ./(GENERATED_FILES_PATH) directory will include the information of all voxels detected for each frame.
The prior one will seperate by objects. The latter only includes voxel data without any other information, which benefits the file read in visualization program.
Both file output same information as follow:

`Voxel_ID | Voxel_xCoord | Voxel_yCoord | Voxel_zCoord | Voxel_Data1(OW value by default) | Voxel_Data1(u_velocity by default) | Voxel_Data1(v_velocity by default) | Voxel_Data1(w_velocity by default) | Voxel_Data1(salinity by default) | Voxel_Data1(temperature by default)`

## ./Seperated Structures folder
The information for each seperated structures will be generated here with one structure per file.
The BASE_FILE_NAME_center.uocd file output the center coordinates from sea surface to bottom for each structure as follow:

`center_xCoord | center_yCoord | center_zCoord | radius_at_this_depth | velocity_minimum_xCoord | velocity_minimum_yCoord`

The ./Seperated Structure/clockwise and ./Seperated Structure/counterclockwise folder includes the information of every voxels for each structures as follow:
`eddy_center_on_surface_xCoord | eddy_center_on_surface_yCoord | voxel_xCoord voxel_yCoord | voxel_zCoord | Voxel_Data1(OW value by default) | Voxel_Data1(u_velocity by default) | Voxel_Data1(v_velocity by default) | Voxel_Data1(w_velocity by default) | Voxel_Data1(salinity by default) | Voxel_Data1(temperature by default) | 0 | radius_of_this_voxel | clockwiseFlag(clockwise=1 counterclockwise=0) | timeFrame | Object_ID`


