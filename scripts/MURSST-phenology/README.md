# Python Module for Phenology

This repository helps organize and maintain the records of the Phenology Project in PO.DAAC and is intended to be used by the phenology community.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

```
* [python 3.0 and above]
* [numpy]
* [netCDF4]
```

### File List

* subset_dataset.py - L3/L4 Dataset Subsetting Module
* data_info.py - PODAAC dataset infomation module
* phenologyalg.py - Phenology main algorithm module
* phenologyplt.py - Phenology main plotting module
* phenology.cfg - Phenology configuration file
* examples.py - example file

## Config file description

[DEFAULT]
* Dataset_Name = CMC
* Root_Name = CMC
* Parameter_Name = analysed_sst
* Parameter_Error_Name = analysis_error
* Start_Year    = 1993
* Start_Month    = 1
* Start_Day     = 1
* End_Year      = 2018
* End_Month    = 1
* End_Day       = 1
* Lat_Min       =  25.00
* Lat_Max       =  65.00
* Lon_Min       = -83.00
* Lon_Max       =   0.00
* Grid_Skip     = 1
* Climatology_File = sst_cmc_climatology.nc
* Climatology_File_Box = sst_cmc_climatology_box.nc
* Lat_Boxsize = 5.0
* Lon_Boxsize = 5.0

[OUTPUT]
* Output_Dir = Output 

## Running the code

```
python examples.py -c phenology.cfg
```

## Authors

* **Yibo JIang, Ed Armstrong** - *Initial work*

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Phenology project support
* Inspiration
* etc

