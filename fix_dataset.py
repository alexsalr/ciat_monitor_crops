
import os
import xarray as xr

if __name__ == '__main__':
    
location = os.environ['WIN_SVR_DATA']+'Saldana/'

files_names = ['S2_201602.nc', 'S2_201601.nc', 'S2_201512.nc']

files = list(map(lambda x: xr.open_dataset(location+x), files_names))

print(files)

aligned = xr.align(files[0], files[1], files[2], join='left', exclude={'time'})

print(aligned)

for idx, dataset in enumerate(aligned):
    print('saving {}', str(idx))
    dataset.to_netcdf(location+'stack/'+files_names[idx])