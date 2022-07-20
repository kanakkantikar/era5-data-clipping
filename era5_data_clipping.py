"""
Topic: Watershed based ERA5 Data Clipping
Date of last editing: 07/20/2022
"""

import geopandas
import xarray
from osgeo import gdal, ogr
import netCDF4
import os

shape_file = 'D:/LOUP shp/Loup_basin.shp'                           # Put the shape file location
nc_file_path = 'F:/ERA5_Data/Temperature/2m temperature/'           # Put nc directory

nc_file_list = os.listdir(nc_file_path)
sf = geopandas.read_file(shape_file)


def makeMask(lon,lat,res):
    source_ds = ogr.Open(shape_file)
    source_layer = source_ds.GetLayer()
 
    # Create high res raster in memory
    mem_ds = gdal.GetDriverByName('MEM').Create('', lon.size, lat.size, gdal.GDT_Byte)
    mem_ds.SetGeoTransform((lon.min(), res, 0, lat.max(), 0, -res))
    band = mem_ds.GetRasterBand(1)
 
    # Rasterize shapefile to grid
    gdal.RasterizeLayer(mem_ds, [1], source_layer, burn_values=[1])
 
    # Get rasterized shapefile as numpy array
    array = band.ReadAsArray()
 
    # Flush memory file
    mem_ds = None
    band = None
    return array

full_processed_nc_list = []
for file in nc_file_list:
    nc_file = nc_file_path + file
    print(nc_file)
    nc = netCDF4.Dataset(nc_file, 'r')
    sf = geopandas.read_file(shape_file)
   
    lons = nc.variables['longitude'][:]
    # get the latitude information
    lats = nc.variables['latitude'][:]
    # calculate the cellsize
    cellsize = lons[:][1] - lons[:][0]
    mask = makeMask(lons, lats, cellsize)

    ## cropping mask to smaller area
    mask = mask[15:25, 15:35]   ## [latitude, longitude]
    mask = xarray.DataArray(mask, dims=['latitude','longitude'])

    ## cropping whole data to exact mask area
    nc = xarray.open_dataset(nc_file)
    test = nc.isel(latitude=slice(15, 25), longitude=slice(15, 35))

    test.where(mask == 1, drop=True)
    full_processed_nc_list.append(test)

full_masked_nc = xarray.concat(full_processed_nc_list, dim='time')

mean_nc = full_masked_nc.mean(dim=['latitude', 'longitude'])    # For "mean" here (it can be "sum" for other variables)
mean_nc_c = mean_nc - 273.0                                     # Fahrenheit to Celsius conversion

mean_nc_c.to_netcdf('t2m_avg.nc')

print(mean_nc_c)
