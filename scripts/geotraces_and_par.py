from netCDF4 import *
import numpy as np
import numpy.ma as ma
import seawater as sw
import matplotlib.pyplot as plt
import datetime as dt
import os
import pandas as pd
from scipy import interpolate

ctd_data = Dataset('../data/oceanographic_data/GEOTRACES_IDP2017_v2/ctd_sensor_data/netcdf/GEOTRACES_IDP2017_v1_CTD_Sensor_Data.nc')

discrete_data = Dataset('../data/oceanographic_data/GEOTRACES_IDP2017_v2/discrete_sample_data/netcdf/GEOTRACES_IDP2017_v2_Discrete_Sample_Data.nc')

# getting variables from the discrete data
pressure_disc = discrete_data.variables['var1'][:]
temp_disc = discrete_data.variables['var7'][:]
salinity_disc = discrete_data.variables['var8'][:]
depth_disc = discrete_data.variables['var2'][:]

lon_geo_disc = discrete_data.variables['longitude']
lat_geo_disc = discrete_data.variables['latitude']

cruise_name_disc = discrete_data.variables['metavar1'][:]
station_name_disc = discrete_data.variables['metavar2'][:]

# calculate seawater density
density_disc = sw.pden(s=salinity_disc, t=temp_disc, p = pressure_disc)

## getting variables from the ctd data
pressure = ctd_data.variables['var1'][:]
temp = ctd_data.variables['var2'][:]
salinity = ctd_data.variables['var3'][:]
lon_geo = ctd_data.variables['longitude']
lat_geo = ctd_data.variables['latitude']

cruise_name = ctd_data.variables['metavar1'][:]
station_name = ctd_data.variables['metavar2'][:]

# calculate seawater density
density = sw.pden(s=salinity, t=temp, p = pressure)


def get_geo_names(geo_names_masked_array):
    
    ### extract names from the cruise name variable
    ### this function gets the string and decodes it and then appends to a list
    
    list_of_names = []
    for string_element in range(len(geo_names_masked_array)):
        blah = geo_names_masked_array[string_element]
        blank_str = ''
        for str_i in range(len(blah)):
            blank_str = blank_str + blah[str_i].decode()
        list_of_names.append(blank_str)
    
    return(list_of_names)

def get_mixed_layer_depth(density_nc):
    
    # gets the mixed layer depth from the ctd data
    
    mld_density_list = []
    # loop through each station
    for station in range(len(density_nc)):
        # get the density at 10 meters (the 9th index) and add 0.03
        mld_density_raw = density_nc[station][9]
        
        if not isinstance(mld_density_raw, np.float32):
            mld_density_raw = density_nc[station][10] + 0.03
        
        if not isinstance(mld_density_raw, np.float32):
            mld_density_raw = density_nc[station][8] + 0.03
            
        else:
            print('argh')
            
        mld_density = mld_density_raw + 0.03

        # add this density to the mld_density list
        mld_density_list.append(mld_density)
    
    # go through each station and through the depth profile
    depth_val = []
    depth_val_na = []
    for station_i in range(len(mld_density_list)):
        # the target density layer you are looking for the corresponding depth
        target_density = mld_density_list[station_i]
        # density profile at the given station
        den_0 = density[station_i]
        # go through the depth profile
        # for each density value, is this value less than the mld density?
        for density_value_at_depth in range(9, len(den_0)):
        #     if yes, continue to next metre density
            if den_0[density_value_at_depth] < target_density:
                continue

        #     if not stop and append the depth value
            elif den_0[density_value_at_depth] > target_density:
                depth_val.append(density_value_at_depth + 1)
                break
            
            else:
                depth_val.append('Na')
                depth_val_na.append(station_i)
                break
    return([depth_val,depth_val_na])

def get_month_from_geotraces(geotraces_ctd_data, ctd_or_disc = 'ctd'):
    
    ## gets the month associated with a date and time from the ctd data or from the discrete data
    
    # getting date time
    ctd_dates = geotraces_ctd_data.variables['date_time'][:].tolist()
    # the date_time var is Decimal Gregorian Days of the station	days since 2007-01-01 00:00:00 UTC
    if ctd_or_disc == 'ctd':
        first_sample_ref_date = dt.date(2007, 1, 1)
    if ctd_or_disc == 'disc':
        first_sample_ref_date = dt.date(6, 1, 1)
    
    geotraces_months = []
    for date_i in range(len(ctd_dates)):
        target_date = first_sample_ref_date + dt.timedelta(days=ctd_dates[date_i])
        geotraces_months.append(target_date.month)
    
    return(geotraces_months)

def discrete_data_mld(density_profile, depth_profile):
    
    ## gets mixed layer depth from the discrete data
    
    ## linear interpolation of density from density profile
    ## quality control that two data points must be above 25 meters
    ## from linear interpolation calculate the density at 10m
    ## figure out the depth at density 10m + 0.03 rho
    ## return MLD per station

    mld_depths = []
    
    for profile_i in range(len(density_profile)):

        density_profile_i = density_profile[profile_i]
        depth_profile_i = depth_profile[profile_i]

        ## check that there are at least two density measurements above 25m
        # print(depth_profile_i  < 25)
#         print(profile_i)
        if (depth_profile_i[depth_profile_i < 25].count() < 3) or (density_profile_i.count() < 3):
            mld_depths.append('NA')
            continue

        # making the surface density equal to the least dense water
        append_depth = ma.append(depth_profile_i, ma.masked_values([0.01], 0))
        append_density = ma.append(density_profile_i, ma.masked_values([density_profile_i.min()], 0))
        

        f = interpolate.interp1d(append_depth, append_density)
        rho_at_ten = f(10) + 0.03
        
        if rho_at_ten > density_profile_i.max():
            mld_depths.append('NA')
            continue
        
        f2 = interpolate.interp1d(append_density, append_depth)
        mld_output = f2(rho_at_ten)
        mld_depths.append(mld_output.item(0))
    
    return(mld_depths)

geo_names_formatted_disc = get_geo_names(ma.getdata(cruise_name_disc))
geo_station_names_formatted_disc = get_geo_names(ma.getdata(station_name_disc))

geo_names_formatted = get_geo_names(ma.getdata(cruise_name))
geo_station_names_formatted = get_geo_names(ma.getdata(station_name))

mld_disc = discrete_data_mld(density_profile=density_disc, depth_profile=depth_disc)
geotraces_months_disc = get_month_from_geotraces(geotraces_ctd_data=discrete_data, ctd_or_disc='disc')

mld = get_mixed_layer_depth(density_nc=density)
geotraces_months = get_month_from_geotraces(geotraces_ctd_data=ctd_data, ctd_or_disc ='ctd')


base_dir = '../data/oceanographic_data/ocean_colour_data/'
file_list = os.listdir(base_dir)
par_file_list = ["PAR" in i for i in file_list]
k_file_list = ["490" in i for i in file_list]

month_order_1 = ['2','7','11','5','10','4','12','9','6','3','1','8']
month_order_3 = ['12','1','4','7','5','11','3','8','6','10','9','2']

def get_nc_list(file_list):
    
    list_of_par_nc = []
    list_of_k_nc = []
    k_file_list_names = []
    par_file_list_names = []

    for un_file in range(len(file_list)):
        print(file_list[un_file])
    
        if par_file_list[un_file]:
            nc_temp = Dataset(base_dir + file_list[un_file])
            list_of_par_nc.append(nc_temp)
            par_file_list_names.append(file_list[un_file])
        
        if k_file_list[un_file]:
            nc_temp = Dataset(base_dir + file_list[un_file])
            list_of_k_nc.append(nc_temp)
            k_file_list_names.append(file_list[un_file])
            
    return([list_of_par_nc, par_file_list_names, list_of_k_nc, k_file_list_names])

# longitude for the par data is from -180 to 180, but longitudes for geotraces data is from 0-360
# must convert the par data longitude by adding *-1 and adding 180 to all negative numbers
def convert_long(long_180_m180):
    for long_i in range(len(long_180_m180)):
        if long_180_m180[long_i] < 0:
            long_180_m180[long_i] = long_180_m180[long_i]*-1 + 180
    return(long_180_m180)

def get_closest_par(get_nc_list_output, month_order_1, month_order_2,
                    long_stations, lat_stations, geotrace_month_per_st):
    
    # this function takes the above list of nc files for par and kd
    # along with the np array of long and lats from the geotraces data
    # and the month sampled for each station
    # the output is a 2 element list, each with a list of par and Kd_490 values for each station
    # in the GEOTRACES dataset
    
    station_kd_vars = []
    station_par_vars = []
    
    # go through each station
    for station_i in range(len(long_stations)):
        
        target_lon = long_stations[station_i]
        target_lat = lat_stations[station_i]
        
        ## something here that figures out the month from geotraces
        temp_month = str(geotrace_month_per_st[station_i])
        
        par_file_index = month_order_1.index(temp_month)
        k_file_index = month_order_2.index(temp_month)
        
        print(par_file_index)
        print(k_file_index)
        
        
        # get the longitude and latitude from those files
        lat_var = get_nc_list_output[0][par_file_index].variables['lat'][:]
        long_var = get_nc_list_output[0][par_file_index].variables['lon'][:]
        # convert the long to 0 -- 360
        lon_converted = convert_long(long_var)
        
        # get the par values from the file
        par = get_nc_list_output[0][par_file_index].variables['par'][:]
        k_val = get_nc_list_output[2][k_file_index].variables['Kd_490'][:]
        
        # find which long and lat are closest to the stations long and lat
        closest_lon = min(lon_converted, key=lambda x:abs(x-target_lon))
        closest_lat = min(lat_var, key=lambda x:abs(x-target_lat))
        
        # find the index of those long and lat
        closest_lat_index = np.where(lat_var == closest_lat)
        closest_lon_index = np.where(long_var == closest_lon)
        
        # get the corresponding par and k vals from that index
        target_kd_val = k_val[closest_lat_index, closest_lon_index]
        target_par_val = par[closest_lat_index, closest_lon_index]

        station_kd_vars.append(target_kd_val)
        station_par_vars.append(target_par_val)
    
    return([station_kd_vars, station_par_vars])


master_par_k = get_nc_list(file_list = file_list)

station_vals_disc = get_closest_par(geotrace_month_per_st=geotraces_months_disc, 
                               long_stations = lon_geo_disc, lat_stations = lat_geo_disc,
                              get_nc_list_output = master_par_k,
                               month_order_1 = month_order_1,
                              month_order_2 = month_order_3)

def get_median_light(get_closest_par_output_test, 
                     station_names_test, 
                     cruise_names_test, 
                     lats_test, 
                     longs_test, 
                     mld_function_test):
    
    par_list = []
    k_list = []
    
    for un_val in range(len(get_closest_par_output_test[0])):
        k_list.append(get_closest_par_output_test[0][un_val].tolist()[0][0])
        par_list.append(get_closest_par_output_test[1][un_val].tolist()[0][0])
    
    # convert par units from einstein m^-2 day^-1 into microE m^-2 s^-1
    station_surface_par_converted = par_list#*(1e6)/86400
    
    d = {'cruise':cruise_names_test,
         'station':station_names_test,
        'kd_490':k_list,
        'par':par_list,
        'lat':lats_test,
        'long':longs_test,
        'mld':mld_function_test}
#     print(d)
    df = pd.DataFrame(d)
    
    return(df)

df = get_median_light(get_closest_par_output_test = station_vals_disc,
                     station_names_test = geo_station_names_formatted_disc, 
                     cruise_names_test = geo_names_formatted_disc, 
                     lats_test = lat_geo_disc, 
                     longs_test = lon_geo_disc, 
                     mld_function_test = mld_disc)
    
df.to_csv('../data/oceanographic_data/geotraces_mld_light_by_station_discrete_data.csv')
