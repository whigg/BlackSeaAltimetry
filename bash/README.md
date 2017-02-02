Bash script to perform various operation on the data files containing altimetry measurements.
* *blacksea_netcdf_ascii*: transform the netcdf files obtained from AVISO FTP into ascii files directly readable by DIVA.
* *blacksea_concat_days*: concatenate the daily files (obtained with * *) according to the time of measurement.
* *blacksea_concat_missions*: concatenate the files corresponding to common dates but with different missions.
* *blacksea_apply_weight_time*: compute the weights according to the time difference, using a gaussian distribution.
* *blacksea_altimetry_count_data*: count the number of data for each daily file
* *blacksea_analysis_V1*: loop over all the files corresponding to a given year and perform analysis.
* *blacksea_seminorm_V1*: same of before but perform semi-normed analysis.

As the file format (netCDF) might have changed since 2013, the first script *blacksea_netcdf_ascii* shall be adapted in order to provide the same ascii files.

