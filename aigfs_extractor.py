import os
import requests
from datetime import datetime
import xarray as xr
import pymap3d as pm
from math import *


def get_latest_model_run():
    utc_hour = datetime.utcnow().hour
    print(utc_hour)
    # more because model takes time to generate
    if utc_hour >= 22:
        return ("18", utc_hour-18)
    elif utc_hour >= 16:
        return ("12", utc_hour-12)
    elif utc_hour >= 10:
        return ("06", utc_hour-6)
    else:
        return ("00", utc_hour-0)

def get_prediction_grib2_link():
    today = datetime.utcnow().strftime("%Y%m%d")
    # prediction model is run 4 times daily, 00, 06, 12, 18
    model_run, hours_after = get_latest_model_run()
    # model predicts future weather in intervals of 6 hours, see how many units have passed since model run
    nearest_prediction = str(int(round(hours_after/6) * 6)).zfill(3)
    #print(nearest_prediction)
    link = f"https://nomads.ncep.noaa.gov/pub/data/nccf/com/aigfs/prod/aigfs.{today}/{model_run}/model/atmos/grib2/aigfs.t{model_run}z.pres.f{nearest_prediction}.grib2"
    return link

def download_grib2_to_data():
    os.makedirs("data", exist_ok=True)

    url = get_prediction_grib2_link()
    filename = os.path.join("data", url.split("/")[-1])

    if os.path.exists(filename):
        print(f"Already downloaded dataset {filename}")
        return filename
    print(f"Downloading {url} ...")
    response = requests.get(url, stream=True)
    response.raise_for_status()

    with open(filename, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)

    print(f"Downloaded to {filename}")
    return filename





class WindPredictor():
    def __init__(self):
        self.local_file = download_grib2_to_data()
        self.ds = xr.open_dataset(self.local_file, engine='cfgrib')
        self.ds = self.ds.sortby("latitude")
        print(self.ds)
    def getWindLatLong(self, pressure, lat, lon):
        if lon < 0:
            lon = 360 + lon
        pressure = max(min(pressure, float(self.ds.isobaricInhPa.max())),
                   float(self.ds.isobaricInhPa.min()))
        """
        u_val = self.ds['u'].sel(isobaricInhPa=pressure, latitude=lat, longitude=lon, method="nearest").values
        v_val = self.ds['v'].sel(isobaricInhPa=pressure, latitude=lat, longitude=lon, method="nearest").values"""
        u_val = self.ds['u'].interp(
            isobaricInhPa=pressure,
            latitude=lat,
            longitude=lon
        ).values

        v_val = self.ds['v'].interp(
            isobaricInhPa=pressure,
            latitude=lat,
            longitude=lon
        ).values
        
        #print(self.ds['u'].sel(isobaricInhPa=pressure, latitude=lat, longitude=lon, method="nearest"))
        return u_val, v_val
    def getWindEC(self, x, y, z, altitude, pressure, debug=True):
        lat, lon, alt = pm.ecef2geodetic(x, y, z)
        if lon<0:
            lon = 360+lon
        # pressure from pascals to millibars
        u_val, v_val = self.getWindLatLong(pressure/100, lat, lon)
        if debug:
            print(f"Calc wind at {lat, lon, altitude, pressure/100} is {u_val, v_val}")
        return pm.enu2ecefv(u_val, v_val, 0, lat, lon, altitude)
    
