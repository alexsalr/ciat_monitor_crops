
mkdir tmp
cd tmp

## Install Anaconda2
# Download anaconda (for Python 2.7.14) ## version is hard-coded. TODO check for further updates
curl -O https://repo.continuum.io/archive/Anaconda2-5.1.0-Linux-x86_64.sh
# Install anaconda, create environment snappy, activate environment
bash Anaconda2-5.1.0-Linux-x86_64.sh

## Create pre-processing conda environment
# prepend the Anaconda2 install location to PATH? yes
# install Ms Vcode? no
# activate anaconda
source ~/.bashrc
# Create environment for preprocessing
conda create --name pre-process python=2
# activate the environment for futher installation
source activate pre-process

## Download SNAP
curl -O http://step.esa.int/downloads/6.0/installers/esa-snap_sentinel_unix_6_0.sh
# Install SNAP
bash esa-snap_sentinel_unix_6_0.sh
## Configure snappy to pre-process conda environment
cd ~/.snap/snap-python/snappy/
python setup.py install
## Update modules in SNAP
bash snap --nosplash --nogui --modules --update-all

## Install sen2cor
cd ~/
mkdir sen2cor
cd sen2cor
curl -O http://step.esa.int/thirdparties/sen2cor/2.5.5/Sen2Cor-02.05.05-Linux64.run
bash Sen2Cor-02.05.05-Linux64.run

## Install sentinelsat and parmap
pip install sentinelsat
conda install -c marufr parmap

## Deactivate conda environment
source deactivate

## You can call L2A processor with '~/sen2cor/Sen2Cor-02.05.05-Linux64/bin/L2A_Process'
## Default configuration file is '~/sen2cor/2.5/cfg/L2A_GIPP.xml'

## To download eo products:
## Install sentinelsat
## pip install sentinelsat
## espa-bulk-downloader: to download Landsat orders
## pip install git+https://github.com/USGS-EROS/espa-bulk-downloader.git
## download_espa_order.py -h
## use example:
## python ./download_espa_order.py -e alex.salazarr@gmail.com -d /mnt/workspace_cluster_6/TRANSVERSAL_PROJECTS/MADR/COMPONENTE_2/Imagenes_Satelitales/Temp/asalazar_tests/ -u asalazarr



## Create a second environment for load and analysis
conda create --name read-explore python=2

source activate read-explore
conda install -c anaconda rasterio xarray netCDF4 cartopy geopandas dask

#conda install -c anaconda scipy pyproj  shapely joblib
#pip install salem

## A third environment for ml

conda install -c anaconda h2o h2o-py

### For Landsat files preprocessing
#cd ~/tmp
#wget https://www.ars.usda.gov/ARSUserFiles/80420510/software/364/public_AROP_v2.2.9.tar.gz ~/tmp

## Configure location for data download and processing
echo "# Added by ciat_monitor_crops" >> ~/.bashrc 
echo "export WIN_SVR_DATA='/mnt/workspace_cluster_6/TRANSVERSAL_PROJECTS/MADR/COMPONENTE_2/Imagenes_Satelitales/Temp/asalazar_tests/data/'" >> ~/.bashrc 
echo "export LOCAL_DATA='/home/azalazar/data/'" >> ~/.bashrc 

source ~/.bashrc
