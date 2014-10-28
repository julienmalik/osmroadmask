#!/usr/bin/env python

#  Copyright 2014 Julien Malik
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

import sys
import argparse
import os

# Check if GDAL is installed
try:
  from osgeo import ogr, osr, gdal
except:
  sys.exit('ERROR: cannot import GDAL/OGR python modules. Did you install python-gdal ?')

# Enable exceptions instead of error codes
gdal.UseExceptions()

# Check if requests is installed
try:
  import requests
except:
  sys.exit('ERROR: cannot import requests python modules. Did you install python-requests ?')

# Now provided directly from source dir 
#requests__version_info__ = tuple([ int(num) for num in requests.__version__.split('.')])
#if requests__version_info__ < (2,4,3):
#except:
#  sys.exit('ERROR: use >= 2.4.3 for requests')

# Check if gdal_translate is installed
try:
  from distutils import spawn
  gdalbin_path = spawn.find_executable("gdal_translate")
  if not gdalbin_path :
    raise ValueError()
except:
  sys.exit('ERROR: cannot find gdal_translate. Did you install gdal-bin ?')

def apply_geotransform(pixel, geotransform):
    x, y = (pixel[0], pixel[1])
    X = geotransform[0] + x*geotransform[1] + y*geotransform[2]
    Y = geotransform[3] + x*geotransform[4] + y*geotransform[5]
    return (X, Y)

def get_lat_lon_bbox(image):
    print "get_lat_lon_bbox"
    src_ds = gdal.Open( image )
    
    crs, geotransform = (src_ds.GetProjectionRef(), src_ds.GetGeoTransform())
    if not crs:
      raise ValueError("Input image is not georeferenced")

    south,west,north,east=(0,0,0,0)
    # http://stackoverflow.com/questions/13439357/extract-point-from-raster-in-gdal
    
    # get the WGS84 spatial reference
    wgs84_sr = osr.SpatialReference()
    wgs84_sr.ImportFromEPSG(4326)
    
    # get the spatial reference of input raster
    raster_sr = osr.SpatialReference(crs);
    
    # instantiate transform
    ct = osr.CoordinateTransformation(raster_sr, wgs84_sr)
    
    
    # get lat/lon of 4 corners
    w, h = (src_ds.RasterXSize, src_ds.RasterYSize)
    ul = apply_geotransform( [0,0], geotransform )
    ur = apply_geotransform( [0,w], geotransform )
    ll = apply_geotransform( [h,0], geotransform )
    lr = apply_geotransform( [h,w], geotransform )
    
    ul_lonlat = ct.TransformPoint(ul[0], ul[1], 0)
    ur_lonlat = ct.TransformPoint(ur[0], ur[1], 0)
    ll_lonlat = ct.TransformPoint(ll[0], ll[1], 0)
    lr_lonlat = ct.TransformPoint(lr[0], lr[1], 0)
    
    # compute bbox
    south = min( [ul_lonlat[1], ur_lonlat[1], ll_lonlat[1], lr_lonlat[1]] )
    north = max( [ul_lonlat[1], ur_lonlat[1], ll_lonlat[1], lr_lonlat[1]] )
    west  = min( [ul_lonlat[0], ur_lonlat[0], ll_lonlat[0], lr_lonlat[0]] )
    east  = max( [ul_lonlat[0], ur_lonlat[0], ll_lonlat[0], lr_lonlat[0]] )
    
    bbox = (south,west,north,east)
    return bbox

def get_highway_types():
  # http://wiki.openstreetmap.org/wiki/FR:Map_Features#Route_.28highway.29
  
  return [
    "motorway",
    "trunk",
    "primary",
    "secondary",
    "tertiary",
    "unclassified",
    "residential",
    "service",
    
    "motorway_link",
    "trunk_link",
    "primary_link",
    "secondary_link",
    "tertiary_link",
    
    "living_street",
    "pedestrian",
    "track",
    "bus_guideway",
    "raceway",
    "road",
    
    "footway",
    "cycleway",
    "bridleway",
    "steps",
    "path",
    
    "cycleway",
    ]
    
def make_overpass_query(bbox, filename):
    print "make_overpass_query"
    # http://wiki.openstreetmap.org/wiki/Overpass_API/Overpass_QL
  
    template = """
    [timeout:25][out:xml][bbox:{south},{west},{north},{east}];
    (
      way["highway"~"{highway_types}"];
      >;
    );
    out body;
    """
    
    query = template.format(south=bbox[0], west=bbox[1], north=bbox[2], east=bbox[3],
                              highway_types="|".join(get_highway_types()))
    
    from contextlib import closing
    with closing(requests.get("http://overpass-api.de/api/interpreter",
                              params={"data":query}, stream=True)) as r:

      if r.status_code != 200:
        if r.status_code == 400:
          message = 'Query syntax error'
        else:
          message = 'Error from Overpass API'
        raise RuntimeError(r.status_code, message)

      # success : dump output to filename
      print "Retrieving roads data..."
      chunk_size = 10240 # bytes
      with open(filename, 'wb') as fd:
        for chunk in r.iter_content(chunk_size):
          fd.write(chunk)

def reproject_vector(args):
    print "reproject_vector"
  
    src_ds = gdal.Open(args.image)
    target_sr = osr.SpatialReference( src_ds.GetProjectionRef() )
    proj4=target_sr.ExportToProj4()
    
    input_layername="lines"
    output_layername="osmroads"
    
    #command_template = "ogr2ogr -t_srs '{proj4}' -f 'SQLite' -dsco SPATIALITE=YES -gt 65536 -nln {output_layername} {reproj_vector} {input_osm_vector} {input_layername}"
    #command_template = "ogr2ogr -t_srs '{proj4}' -f 'SQLite' -dsco SPATIALITE=YES -gt 65536 -lco SPATIAL_INDEX=NO -nln {output_layername} {reproj_vector} {input_osm_vector} {input_layername}"
    command_template = "ogr2ogr -overwrite -t_srs '{proj4}' -f 'SQLite' -gt 65536 -nln {output_layername} {reproj_vector} {input_osm_vector} {input_layername}"
    #command_template = "ogr2ogr -t_srs '{proj4}' -f 'ESRI Shapefile' -nln  {output_layername} {reproj_vector} {input_osm_vector} {input_layername}"
    
    command = command_template.format(proj4=proj4, reproj_vector=args.reproj, input_osm_vector=args.osm, output_layername=output_layername, input_layername=input_layername)
    os.system(command)
    
    #print "creating index"
    #command_template = 'ogrinfo {sqlite} -sql "SELECT CreateSpatialIndex(\'{layername}\',\'GEOMETRY\')"'
    #command = command_template.format(sqlite=args.reproj, layername=output_layername)
    #os.system(command)
    

def create_output_raster(args):
    print "create_output_raster"
    
    command_template = "GDAL_PAM_ENABLED=NO gdal_translate {input} {output_mask} -of GTiff -b 1 -ot Byte -scale 0 65525 0 0"
    #command_template = "GDAL_PAM_ENABLED=NO gdal_translate {input} {output_mask} -of ENVI -b 1 -ot Byte -scale 0 65525 0 0"
    command = command_template.format(input=args.image, output_mask=args.mask)
    os.system(command)
    
def rasterize(args):
    print "rasterizing"
    command_template = "gdal_rasterize -b 1 -burn 255 -l {output_layername} {reproj_vector} {output_mask}"
    command = command_template.format(output_layername="osmroads", reproj_vector=args.reproj, output_mask=args.mask)
    os.system(command)
  
def validate_args(args):
    if os.path.splitext(args.osm)[1] != ".osm":
      raise ValueError("use .osm for osm vector data filename extension")

def main(args):
    validate_args(args)
    
    bbox = get_lat_lon_bbox(args.image)
    make_overpass_query(bbox, args.osm)
    reproject_vector(args)
    create_output_raster(args)
    rasterize(args)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
      description="""
                  Given an input georeferenced image, 
                  computes a superimposed mask image 
                  containing road network extracted from OSM
                  """,
       epilog="""
       Example : osm_road.py input.tif rawosmdata.osm vectortorasterize.sqlite roadmask.tif
       """
       )
    parser.add_argument('image', metavar='image',
                    help='the input georeferenced image')
    parser.add_argument('osm', metavar='output_osm',
                    help='the roads vector data in osm XML format')
    parser.add_argument('reproj', metavar='output_proj',
                    help='the roads vector data reprojected in the input raster coordinate system')
    parser.add_argument('mask', metavar='mask',
                    help='the output road mask')
    args = parser.parse_args()

    main(args)