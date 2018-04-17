
from sentinel1_pre import *
from sentinel2_pre import *

if __name__ == '__main__':
    
	data_dir = '/home/azalazar/data/S1_Ibague/'
	out_dir = '/home/azalazar/data/pre/'
	area_of_int = '/home/azalazar/data/spatial_ref/ibagueextent.geojson'
	ref_raster = '/home/azalazar/data/spatial_ref/ibague.dim'
	unzipfiles(data_dir)
	pre_process_s1(data_dir, out_dir, area_of_int, ref_raster, polarizations=['VV','VH'])
	
	