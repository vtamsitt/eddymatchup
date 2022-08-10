# eddymatchup

This package colocates any oceanographic position data (Latitude, Longitude, and Time) to eddy trajectories identified in eddy databases available from AVISO. Data url: https://www.aviso.altimetry.fr/en/data/products/value-added-products/global-mesoscale-eddy-trajectory-product.html

Future versions of this package will add additional publicly available eddy tracking databases to the eddy matchups, including near-real time products.

We thank Jennifer Bonin and Don Chambers at the University of South Florida for providing their original fortran eddy match up code on which this package is based.

Typical usage looks like:
	
	import xarray as xr
	import eddymatchup as ed
        
 	ds = xr.open_dataset('your_favourite_oceanographic_data.nc')
       
	matches = ed.match(ds.lon,ds.lat,ds.datetime,database='META3.2_allsat',latmin=-90,latmax=-35,hourrange=24,radiusrange=1)


To Do: 
* Add Environment file
