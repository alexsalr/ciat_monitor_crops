import math
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import datetime as dt

@xr.register_dataarray_accessor('eotemp')
class EOTempArray(object):
    
    def __init__(self, xarray_obj):
        self._obj = xarray_obj
        
    def plot_RGB(self, bands=['red', 'green', 'blue'], factor=0.0001):
        
        crs = ccrs.UTM('18N')
        av_dates = self._obj.coords['time'].data.tolist()
        
        plt.figure(figsize=(15,3*math.ceil(len(av_dates)/5.0)))
        try:
            rgb = self._obj.sel(band=['red', 'green', 'blue']).sortby('time')
        except:
            print('The EOTempArray does not contain rgb bands')
            return None
        for ix, date in enumerate(av_dates):
            ap = rgb.isel(time=ix)*factor
            ax = plt.subplot(math.ceil(len(av_dates)/5.0),5,1+ix, projection=crs)
            ap.plot.imshow(rgb='band', transform=crs)
        
        #plt.savefig("/home/azalazar/data/pre/stacks/rgb.pdf", dpi=300)
        plt.show()
        
#xr.register_accessor('spec', EOTempArray)