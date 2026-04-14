#!/usr/bin/python
# -*- coding: utf-8 -*-

import os, sys
import pygplates
import numpy as np

from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
from scipy import ndimage

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pandas as pd

import cartopy.crs as ccrs

from shapely.geometry import Polygon

ROTATION_DIR = os.path.dirname(os.path.realpath(__file__)) + '/data'
DATA_DIR = './output_ATOM_Precipitation'
BATHYMETRY_SUFFIX = 'Ma_smooth.xyz'

script_dir = os.path.dirname(os.path.realpath(__file__)) 
ATOM_HOME = script_dir + '/../'
TOPO_DIR = ATOM_HOME + '/data/topo_grids/'                              #  path of the grid files
                                                                        #assign plate id to the points_on_land features
                                                                        #we only reconstruct the points on continents

def convert_atom_to_gmt(data):
    # 1. Reshape the flattened data to its original 2D structure
    # Based on your indexing [j1*181+i], the original shape is (361, 181)
    reshaped = data.flatten().reshape(361, 181)
    
    # 2. Transpose to get (181, 361) to match your new_data dimensions
    transposed = reshaped.T  # Now shape is (181, 361)
    
    # 3. Perform the longitude shift (equivalent to your j > 179 logic)
    # np.roll shifts the elements. Shifting by 180 indices swaps the hemispheres.
    new_data = np.roll(transposed, shift=180, axis=1)
    
    return new_data



def convert_gmt_to_atom(data):
    """
    Shifts longitude hemispheres and transposes GMT (lat, lon) data 
    to ATOM (lon, lat) format.
    """
    # 1. Swap hemispheres: lon < 180 moves to the end, lon >= 180 moves to the start.
    # This is equivalent to shifting the array by 180 positions along the longitude axis.
    shifted_data = np.roll(data, shift=-180, axis=1)
    
    # 2. Transpose from (181, 361) to (361, 181)
    # This matches your new_data[lon][90-lat] indexing.
    return shifted_data.T



def add_lon_lat_to_gmt_data(data):
    # 1. Define coordinates matching the -Rd -I1. grid (181 lats, 361 lons)
    lat = np.linspace(90, -90, 181)
    lon = np.linspace(-180, 180, 361)
    
    # 2. Create 2D grids of coordinates
    lon_2d, lat_2d = np.meshgrid(lon, lat)
    
    # 3. Ensure data is the correct 2D shape
    z_2d = data.reshape((181, 361))
    
    # 4. Flatten and stack into (N, 3) columns: [Lon, Lat, Z]
    # This format is required for writing to .xyz files for GMT
    return np.column_stack((lon_2d.ravel(), lat_2d.ravel(), z_2d.ravel()))



def add_lon_lat_to_atom_data(data):
    # 1. Define coordinate vectors
    lat = np.linspace(90, -90, 181)
    lon = np.linspace(0, 360, 361)
    
    # 2. Create 2D grids (Note: ATOM format uses lon as the first dimension)
    # Using indexing='ij' ensures Y (lat) aligns with the second dimension of 'data'
    lon_2d, lat_2d = np.meshgrid(lon, lat, indexing='ij')
    
    # 3. Flatten and stack into (N, 3) columns: [Lon, Lat, Z]
    # This is the standard format for GMT .xyz files
    return np.column_stack((lon_2d.ravel(), lat_2d.ravel(), data.ravel()))



def interp_grid(a):                                                     # obviously not applied
    d = np.array(a)
    x, y, z = d[:, 0], d[:, 1], d[:, 2]
    
    # Define the 1-degree global grid
    lat = np.linspace(90, -90, 181)
    lon = np.linspace(-180, 180, 361)
    X, Y = np.meshgrid(lon, lat)
    
    # Perform cubic interpolation
    # Note: fill_value=0 might create artificial "cliffs" at data edges
    grid_data = griddata((x, y), z, (X, Y), method='cubic', fill_value=0)

    # Efficiently clip values to original min/max bounds
    return np.clip(grid_data, z.min(), z.max())



def interp_grid_gmt(a):

    data = np.array(a)

    zmax = data[:,2].max()
    zmin = data[:,2].min()

    # Create the directory if it doesn't exist
    if not os.path.isdir("/tmp/atom/"):
        os.makedirs('/tmp/atom/', exist_ok=True)

    # Define paths and grid parameters
    input_xyz =   '/tmp/atom/no_land.xyz'
    cleaned_xyz = '/tmp/atom/cleaned_data.xyz'
    output_nc =   '/tmp/atom/fill_land_gap.nc'
    output_xyz =  '/tmp/atom/fill_land_gap.xyz'

    # Write data using f-strings for formatting
    with open(input_xyz, 'w') as of:
        for line in data:
            # line is expected to be [x, y, z]
            of.write(f"{line[0]} {line[1]} {line[2]}\n")

    # 1. Preprocess: Combine data points within the same 1-degree block
    # -Rd (Global), -I1. (1-degree resolution)
    os.system(f'gmt blockmean {input_xyz} -Rd -I1. > {cleaned_xyz}')

    # 2. Interpolate: Create the grid using the cleaned data
    # -Q0.8 applies tension to the spherical surface
    os.system(f'gmt sphinterpolate {cleaned_xyz} -G{output_nc} -Rd -I1. -Q0.8')

    # 3. Clip values: Anything below zmin or above zmax becomes 'NaN' (no-data)
    os.system(f'gmt grdclip {output_nc} -G{output_nc} -Sb{zmin}/NaN -Sa{zmax}/NaN')

    # 4. Convert grid to XYZ: Outputs x, y, z columns for every grid node
    os.system(f'gmt grd2xyz {output_nc} > {output_xyz}')

    # 5. Load the XYZ data into a NumPy file. Floats are still the standard for coordinates.
    data = np.genfromtxt(output_xyz, dtype=float)   
    
    return data[:,2] 



def gmt_filter(a, zmin=None, zmax=None):                                # obviously not used for reconstruction
    data_array = np.array(a)
    # Ensure directory exists
    os.makedirs('/tmp/atom/', exist_ok=True)
    
    # 1. Write input data to XYZ
    input_xyz = '/tmp/atom/result.xyz'
    with open(input_xyz, 'w') as of:
        for line in data_array:
            of.write(f"{line[0]} {line[1]} {line[2]}\n")

    # 2. Grid the data (MANDATORY: grdfilter requires a grid, not an XYZ file)
    # Using -Rd for global and -I1. for 1-degree resolution
    input_nc = '/tmp/atom/input_to_filter.nc'
    os.system(f'gmt sphinterpolate {input_xyz} -G{input_nc} -Rd -I1. -Q0.8')

    # 3. Apply the filter
    # -Fm7: Median filter with 7km (or degrees, depending on -D) width
    filtered_nc = '/tmp/atom/result.nc'
    os.system(f'gmt grdfilter {input_nc} -G{filtered_nc} -Fm7 -Dp -Vl')

    # 4. Optional: Clip if bounds were provided
    if zmin is not None and zmax is not None:
        os.system(f'gmt grdclip {filtered_nc} -G{filtered_nc} -Sb{zmin}/{zmin} -Sa{zmax}/{zmax}')

    # 5. Convert back to XYZ and load
    output_xyz = '/tmp/atom/result_filtered.xyz'
    os.system(f'gmt grd2xyz {filtered_nc} > {output_xyz}')
    
    # Using genfromtxt to handle any NaNs produced
    result_data = np.genfromtxt(output_xyz)
    
    # Return the Z-column (3rd column)
    return result_data[:, 2]



def get_coastline_polygons_from_topography(filename):
    # Load and prepare data (standard 181x361 grid)
    data = np.genfromtxt(filename)   
    lons = data[:, 0].reshape((181, 361))
    lats = data[:, 1].reshape((181, 361))
    topo = data[:, 2].reshape((181, 361))
    
    # Binary land/sea mask (Land=1, Sea=0)
    topo_mask = np.where(topo >= 0, 1, 0)

    # Use a hidden Matplotlib plot to extract contour paths
    # PlateCarree is the direct replacement for Basemap's 'cyl' projection
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    cs = ax.contourf(lons, lats, topo_mask, levels=[0.5, 1.5])
    plt.close(fig) # We only need the paths, not the plot

    polygons = []
    # cs.collections[0] contains the paths for the land areas (topo=1)
    for collection in cs.collections:
        for path in collection.get_paths():
            # Convert paths to Shapely polygons (handles holes automatically)
            for poly_coords in path.to_polygons():
                if len(poly_coords) < 3:
                    continue
                # Create Shapely Polygon
                poly = Polygon(poly_coords)
                if poly.is_valid:
                    polygons.append(poly)

    # Convert Shapely Polygons to pyGPlates Feature objects
    polygon_features = []
    for p in polygons:
        # Extract exterior coordinates (y=lat, x=lon)
        # pyGPlates expects (lat, lon) for its sphere geometries
        lons_ext, lats_ext = p.exterior.xy
        
        f = pygplates.Feature()
        # Create PolygonOnSphere from (latitude, longitude) tuples
        f.set_geometry(pygplates.PolygonOnSphere(zip(lats_ext, lons_ext)))
        polygon_features.append(f)
        
    return polygon_features



def remove_data_nearby_coastline_fast(data, topo, len_to_remove=5):
    # 1. Reshape
    data_2d = data.reshape((181, 361))
    topo_2d = topo.reshape((181, 361))
    
    # 2. Create a land mask (1 for land, 0 for ocean)
    land_mask = (topo_2d != 0)
    
    # 3. Dilate the land mask to create the "coastal buffer"
    # This expands land into the ocean by 'len_to_remove' pixels
    struct = ndimage.generate_binary_structure(2, 2) # Includes diagonals
    dilated_land = ndimage.binary_dilation(land_mask, structure=struct, iterations=len_to_remove)
    
    # 4. Identify pixels that are in the buffer but NOT original land
    ocean_buffer = dilated_land & ~land_mask
    
    # 5. Apply NaN to the buffer
    data_2d[ocean_buffer] = np.nan

    return data_2d.flatten()



def reconstruct_grid(from_time, input_grid, to_time, output_grid, reconstruction_dir=script_dir):
    
    data = np.genfromtxt(input_grid)
    data = convert_atom_to_gmt(data[:,2])
    
    if True:
        topo = np.genfromtxt(TOPO_DIR + '/{}Ma_smooth.xyz'.format(from_time))
        topo = topo[:,2]
        topo[topo>0] = True
        topo[topo<=0] = False
        remove_data_nearby_coastline_fast(data, topo, 5)
        
    data = add_lon_lat_to_gmt_data(data)

    static_polygons = pygplates.FeatureCollection(
        reconstruction_dir + '/data/Muller_etal_AREPS_2016_StaticPolygons.gpmlz' )
        
    rotation_files = [reconstruction_dir + '/data/Rotations/Global_EarthByte_230-0Ma_GK07_AREPS.rot']
    rotation_model = pygplates.RotationModel(rotation_files)

    #use matplotlib contour function to extract polygons from topography data
    coastline_polygons_to_time = get_coastline_polygons_from_topography(TOPO_DIR+'/{}Ma_smooth.xyz'.format(to_time))
    coastline_polygons_from_time = get_coastline_polygons_from_topography(TOPO_DIR+'/{}Ma_smooth.xyz'.format(from_time))
   
    #turn grid data into point feature
    #use subductionZoneDepth and subductionZoneSystemOrder properties to keep grid data and point index
    #the reason of using these two properties is I don't know how to store data in a feature in other way
    #if you know better way to store data in a feature, you may change the code below
    points_in_ocean = []
    points_on_land = []
    data = data.reshape((181*361,3))
    for idx, point in enumerate(data):        
        on_land = False
        for p in coastline_polygons_from_time: 
            if p.get_geometry().is_point_in_polygon((float(point[1]),float(point[0]))):
                f = pygplates.Feature()
                f.set_geometry(pygplates.PointOnSphere(float(point[1]),float(point[0])))
                f.set_double(
                    pygplates.PropertyName.create_gpml('subductionZoneDepth'),
                    point[2])
                f.set_integer(
                    pygplates.PropertyName.create_gpml('subductionZoneSystemOrder'),
                    idx)
                points_on_land.append(f)
                on_land=True
                break
        if not on_land and (not np.isnan(point[2])):
            points_in_ocean.append(point)


    #assign plate id to the points_on_land features
    #we only reconstruct the points on continents

    partitioned_features, unpartitioned_features = pygplates.partition_into_plates(
            static_polygons,
            rotation_files,
            points_on_land,
            properties_to_copy = [
                pygplates.PartitionProperty.reconstruction_plate_id],
            reconstruction_time = from_time,
            partition_return = pygplates.PartitionReturn.separate_partitioned_and_unpartitioned
            )

    #use equivalent stage rotation to reconstruct the points and save the data 
    new_points_on_land=[]

    fs = sorted(partitioned_features, key=lambda x: x.get_reconstruction_plate_id())


    from itertools import groupby


    for key, group in groupby(fs, lambda x: x.get_reconstruction_plate_id()):

        fr = rotation_model.get_rotation(to_time, key)

        for f in group:
            ll = (fr * f.get_geometry()).to_lat_lon()
            v = f.get_double(pygplates.PropertyName.create_gpml('subductionZoneDepth'))
            new_points_on_land.append([ll[1], ll[0], v])


    grid_data = interp_grid_gmt(points_in_ocean)                        #fill the gaps in the grid

    grid_data = grid_data.reshape((181,361))

    points_in_ocean_to_time = []
    data = add_lon_lat_to_gmt_data(grid_data)
    data = data.reshape((181*361,3))


    for idx, point in enumerate(data):        
        on_land = False
        for p in coastline_polygons_to_time: 
            if p.get_geometry().is_point_in_polygon((float(point[1]),float(point[0]))):
                on_land=True
                break
        if not on_land:
            points_in_ocean_to_time.append(point)
    
    new_grid_data = points_in_ocean_to_time
    #only keep points which are inside the to_time polygons. 
    #prevent info on continent leaking into oceans
    for row in new_points_on_land: 
        for f in coastline_polygons_to_time:
            if f.get_geometry().is_point_in_polygon((row[1], row[0])):
                new_grid_data.append(row)
                continue
    
    grid_data = interp_grid_gmt(new_grid_data)                          #fill the gaps in the grid
        
    grid_data = grid_data.reshape((181,361))
        
    output_data = convert_gmt_to_atom(grid_data)
    output_data = add_lon_lat_to_atom_data(output_data)
    output_data = output_data.reshape((361*181,3))

    np.savetxt(output_grid,output_data,fmt='%1.2f') 



def reconstruct_velocity_grid(from_time, input_grid, to_time, output_grid, reconstruction_dir=script_dir):
    
    data = np.genfromtxt(input_grid)
    data = convert_atom_to_gmt(data[:,2])
    
    if True:
        topo = np.genfromtxt(TOPO_DIR+'/{}Ma_smooth.xyz'.format(from_time))   
        topo = topo[:,2]
        topo[topo>0] = True
        topo[topo<=0] = False
        remove_data_nearby_coastline_fast(data, topo, 5)
        
    data = add_lon_lat_to_gmt_data(data)

    #use matplotlib contour function to extract polygons from topography data
    coastline_polygons_to_time = get_coastline_polygons_from_topography(TOPO_DIR+'/{}Ma_smooth.xyz'.format(to_time))
    coastline_polygons_from_time = get_coastline_polygons_from_topography(TOPO_DIR+'/{}Ma_smooth.xyz'.format(from_time))
   
    #get points in ocean
    points_in_ocean = []
    data = data.reshape((181*361,3))
    for idx, point in enumerate(data):        
        on_land = False
        for p in coastline_polygons_from_time: 
            if p.get_geometry().is_point_in_polygon((float(point[1]),float(point[0]))):
                on_land=True
                break
        if not on_land and (not np.isnan(point[2])):
            points_in_ocean.append(point)

    grid_data = interp_grid_gmt(points_in_ocean)                        #fill the gaps in the grid
    grid_data = grid_data.reshape((181,361))

    new_grid_data = []
    data = add_lon_lat_to_gmt_data(grid_data)
    data = data.reshape((181*361,3))
    for idx, point in enumerate(data):        
        for p in coastline_polygons_to_time: 
            if p.get_geometry().is_point_in_polygon((float(point[1]),float(point[0]))):
                point[2]=0
                break
        new_grid_data.append(point[2])
        
    grid_data = np.array(new_grid_data).reshape((181,361))
        
    output_data = convert_gmt_to_atom(grid_data)
    output_data = add_lon_lat_to_atom_data(output_data)
    output_data = output_data.reshape((361*181,3))

    np.savetxt(output_grid,output_data,fmt='%1.2f')     



def reconstruct_temperature_fast(time_0, time_1, suffix='Ma_smooth.xyz'):
    # Load data: Column 0=Lon, 1=Lat, 6=Temp
    input_path = f"{DATA_DIR}/{time_0}{suffix.removesuffix('.xyz')}_PlotData_Atm.xyz"
    df = pd.read_csv(input_path, sep='\s+', skiprows=1, usecols=[0, 1, 6], header=None)
    
    # Sort: Latitude descending (-data[:,1]), Longitude ascending
    df = df.sort_values(by=[0, 1], ascending=[True, False])
    
    # Save formatted XYZ
    temp_file = f"{DATA_DIR}/{time_0}Ma_Atm_Temperature.xyz"
    df.to_csv(temp_file, sep=' ', header=False, index=False, float_format='%1.2f')

    # Execute the reconstruction
    reconstruct_grid(
        time_0,
        temp_file,
        time_1,
        f"{DATA_DIR}/{time_1}Ma_Reconstructed_Temperature.xyz"
    )



def reconstruct_precipitation(time_0, time_1, suffix='Ma_smooth.xyz'):
    input_path = f"{DATA_DIR}/{time_0}{suffix.removesuffix('.xyz')}_PlotData_Atm.xyz"
    
    # Load Longitude (0), Latitude (1), and Precipitation (8)
    df = pd.read_csv(input_path, sep='\s+', skiprows=1, usecols=[0, 1, 8], header=None)
    
    # Sort: Lon ascending, Lat descending
    df = df.sort_values(by=[0, 1], ascending=[True, False])
    
    # Save intermediate file
    temp_file = f"{DATA_DIR}/{time_0}Ma_Atm_Precipitation.xyz"
    df.to_csv(temp_file, sep=' ', header=False, index=False, float_format='%1.2f')

    # Run tectonic reconstruction
    reconstruct_grid(
        time_0,
        temp_file,
        time_1,
        f"{DATA_DIR}/{time_1}Ma_Reconstructed_Precipitation.xyz"
    )



def reconstruct_salinity(time_0, time_1, suffix='Ma_smooth.xyz'):
    # Load Hydrological data (Salinity is typically column 7)
    input_path = f"{DATA_DIR}/{time_0}{suffix.removesuffix('.xyz')}_PlotData_Hyd.xyz"
    df = pd.read_csv(input_path, sep='\s+', skiprows=1, usecols=[0, 1, 7], header=None)
    
    # Sort: Longitude ascending, Latitude descending
    df = df.sort_values(by=[0, 1], ascending=[True, False])
    
    # Save intermediate file for the reconstruction function
    temp_file = f"{DATA_DIR}/{time_0}Ma_Hyd_Salinity.xyz"
    df.to_csv(temp_file, sep=' ', header=False, index=False, float_format='%1.2f')

    # Execute tectonic reconstruction
    reconstruct_grid(
        time_0,
        temp_file,
        time_1,
        f"{DATA_DIR}/{time_1}Ma_Reconstructed_Salinity.xyz"
    )



def reconstruct_wind_v(time_0, time_1, suffix='Ma_smooth.xyz'):
    input_path = f"{DATA_DIR}/{time_0}{suffix.removesuffix('.xyz')}_PlotData_Atm.xyz"
    
    # Load Lon(0), Lat(1), and V-Wind(3)
    df = pd.read_csv(input_path, sep='\s+', skiprows=1, usecols=[0, 1, 3], header=None)
    
    # Sort: Longitude ascending, Latitude descending
    df = df.sort_values(by=[0, 1], ascending=[True, False])
    
    # Save intermediate file
    temp_file = f"{DATA_DIR}/{time_0}Ma_Atm_v.xyz"
    df.to_csv(temp_file, sep=' ', header=False, index=False, float_format='%1.2f')

    # Use the velocity-specific reconstruction function
    reconstruct_velocity_grid(
        time_0,
        temp_file,
        time_1,
        f"{DATA_DIR}/{time_1}Ma_Reconstructed_wind_v.xyz"
    )



def reconstruct_wind_w(time_0, time_1, suffix='Ma_smooth.xyz'):
    input_path = f"{DATA_DIR}/{time_0}{suffix.removesuffix('.xyz')}_PlotData_Atm.xyz"
    
    # Load Lon(0), Lat(1), and W-Wind(4)
    df = pd.read_csv(input_path, sep='\s+', skiprows=1, usecols=[0, 1, 4], header=None)
    
    # Sort: Longitude ascending, Latitude descending
    df = df.sort_values(by=[0, 1], ascending=[True, False])
    
    # Save intermediate file
    temp_file = f"{DATA_DIR}/{time_0}Ma_Atm_w.xyz"
    df.to_csv(temp_file, sep=' ', header=False, index=False, float_format='%1.2f')

    # Execute velocity-aware reconstruction
    reconstruct_velocity_grid(
        time_0,
        temp_file,
        time_1,
        f"{DATA_DIR}/{time_1}Ma_Reconstructed_wind_w.xyz"
    )



def test(filename):                                                     #  test file to proove available data
    from_time = 0
    from_file = filename

    for t in range(time_start,time_end,time_step):

        to_file = f"{t}Ma.xyz"

        reconstruct_grid(from_time, from_file, t, to_file)

        # Update for the next iteration
        from_time = t
        from_file = to_file



def main():                                                             #  this python program is called in FileIO_Atm.cpp and executed by int ret = system(cmd_str.c_str());
    try:
#        if len(sys.argv) == 2:                                          #  if two arguments exist (program name and first user argument
        if len(sys.argv) == 3:                                          #  if two arguments exist (program name and first user argument
            test(sys.argv[1])                                           #  preceding timestep, user provided argument in function def test(filename)
            return

        global DATA_DIR                                                 #  DATA_DIR = './output_ATOM_Precipitation'
        global BATHYMETRY_SUFFIX                                        #  BATHYMETRY_SUFFIX = 'Ma_smooth.xyz'

        time_0 = int(sys.argv[1])                                       #  preceding timestep
        time_1 = int(sys.argv[2])                                       #  actual timestep

        DATA_DIR = sys.argv[3]                                          #  DATA_DIR = './output_ATOM_Precipitation'
        BATHYMETRY_SUFFIX = sys.argv[4]                                 #  BATHYMETRY_SUFFIX = 'Ma_smooth.xyz'

        atm_or_hyd = sys.argv[5]                                        #  atmosphere or ocean code running

        print()
        print('main in reconstruct_atom_data.py, preceding timestep = time_0 compares to  Ma = ', time_0)
        print()
        print('main in reconstruct_atom_data.p, actual timestep = time_1 compares to  Ma = ', time_1)
        print()


        if atm_or_hyd == 'atm':
            reconstruct_temperature_fast(time_0, time_1, BATHYMETRY_SUFFIX)
#            reconstruct_precipitation(time_0, time_1, BATHYMETRY_SUFFIX)
#            reconstruct_wind_v(time_0, time_1, BATHYMETRY_SUFFIX)
#            reconstruct_wind_w(time_0, time_1, BATHYMETRY_SUFFIX)
        else:
            reconstruct_salinity(time_0, time_1, BATHYMETRY_SUFFIX)

    except:
        print("Error message: Usage: python reconstruct_atom_data.py 0 10 ./output Ma_smooth.xyz atm/hyd") 
        import traceback
        traceback.print_exc()



if __name__ == "__main__":
    main()


