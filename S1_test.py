
from sentinel1_pre import *

if __name__ == '__main__':
    
	data_dir = '/home/azalazar/data/'
	out_dir = '/home/azalazar/data/pre/'
	area_of_int = '/home/azalazar/data/spatial_ref/ibagueextent.geojson'
	ref_raster = '/home/azalazar/data/spatial_ref/ibagueraster_10m_S2_B3.img'
	
	pre_process_s1(data_dir, out_dir, area_of_int, ref_raster, polarizations=['VV','VH'])
	
	