


import os
import requests
from datetime import datetime, timedelta, timezone
import xarray as xr
import pymap3d as pm
from math import *
from scipy.interpolate import RegularGridInterpolator
import numpy as np
import pandas

def get_model_run(target_time):
    """
    Gets the model run most appropriate for a target time
    """
    # this is the time when the model made its predictions
    model_run_time = target_time
    if model_run_time > datetime.now(tz=datetime.timezone.utc):
        model_run_time = datetime.now(tz=datetime.timezone.utc)
    model_run = 0
    # if target time in last four hours we need prev model run since model takes 4 hours to generate
    if target_time > (datetime.now(tz=timezone.utc) - timedelta(hours=4)):
        model_run_time-=timedelta(hours=4)
        model_run = floor(model_run_time.hour/6)*6
    else:
        # if target time more than 4 hours in past
        model_run = floor(model_run_time.hour/6)*6
    model_run_time = model_run_time.replace(hour=model_run, minute=0, second=0, microsecond=0)

    model_day_str = model_run_time.strftime("%Y%m%d")
    return model_run_time

def get_datasets(target_time):
    model_initial_time = get_model_run(target_time=target_time)
    model_inlet_time = target_time - model_initial_time
    print(model_inlet_time.total_seconds())
    # since model only does predictions every 6 hours, we need two datasets to interpolate between them
    prediction_hour_lower = floor(model_inlet_time.total_seconds()/(3600*6)) * 6
    prediction_hour_upper = ceil(model_inlet_time.total_seconds()/(3600*6)) * 6
    model_initial_time_str = model_initial_time.strftime("%Y%m%d")
    model_run_str = str(model_initial_time.hour).zfill(2)
    prediction_hour_lower_str = str(prediction_hour_lower).zfill(3)
    prediction_hour_upper_str = str(prediction_hour_upper).zfill(3)

    # the links for the two data sets on either side of our current target time
    ds1_link = f"https://nomads.ncep.noaa.gov/pub/data/nccf/com/aigfs/prod/aigfs.{model_initial_time_str}/{model_run_str}/model/atmos/grib2/aigfs.t{model_run_str}z.pres.f{prediction_hour_lower_str}.grib2"
    ds2_link = f"https://nomads.ncep.noaa.gov/pub/data/nccf/com/aigfs/prod/aigfs.{model_initial_time_str}/{model_run_str}/model/atmos/grib2/aigfs.t{model_run_str}z.pres.f{prediction_hour_upper_str}.grib2"

    os.makedirs("data", exist_ok=True)
    ds_dir = os.path.join("data", model_initial_time_str)
    os.makedirs(ds_dir, exist_ok=True)
    ds1_filename = os.path.join(ds_dir, ds1_link.split("/")[-1])
    ds2_filename = os.path.join(ds_dir, ds2_link.split("/")[-1])

    # if the dataset isn't already downloaded, download it
    if not os.path.exists(ds1_filename):
        print(f"Downloading {ds1_filename} ...")
        response = requests.get(ds1_link, stream=True)
        response.raise_for_status()

        with open(ds1_filename, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
    if not os.path.exists(ds2_filename):
        print(f"Downloading {ds2_filename} ...")
        response = requests.get(ds2_link, stream=True)
        response.raise_for_status()

        with open(ds2_filename, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
    return ds1_filename, ds2_filename


class WindPredictor():
    def update_datasets(self, target_time):
        # the two wind datasets, six hours apart, we'll interpolate between
        self.ds1_filename, self.ds2_filename = get_datasets(target_time)
        self.ds1 = xr.open_dataset(self.ds1_filename, engine='cfgrib')
        self.ds1 = self.ds1.sortby("latitude")
        self.ds2 = xr.open_dataset(self.ds2_filename, engine='cfgrib')
        self.ds2 = self.ds2.sortby("latitude")

        self.d1_pressures = self.ds1["isobaricInhPa"].values
        self.d1_lats = self.ds1["latitude"].values
        self.d1_lons = self.ds1["longitude"].values
        self.d1_u_data = self.ds1["u"].values
        self.d1_v_data = self.ds1["v"].values

        self.d2_pressures = self.ds2["isobaricInhPa"].values
        self.d2_lats = self.ds2["latitude"].values
        self.d2_lons = self.ds2["longitude"].values
        self.d2_u_data = self.ds2["u"].values
        self.d2_v_data = self.ds2["v"].values

        self.ds1_valid_time = pandas.to_datetime(self.ds1["valid_time"].values).tz_localize("UTC")
        self.ds2_valid_time = pandas.to_datetime(self.ds2["valid_time"].values).tz_localize("UTC")
        print(self.ds1_valid_time, self.ds1_valid_time)

        self.d1_u_interp = RegularGridInterpolator(
            (self.d1_pressures, self.d1_lats, self.d1_lons), 
            self.d1_u_data, 
            bounds_error=False, fill_value=None
        )
        self.d1_v_interp = RegularGridInterpolator(
            (self.d1_pressures, self.d1_lats, self.d1_lons), 
            self.d1_v_data, 
            bounds_error=False, fill_value=None
        )
        self.d2_u_interp = RegularGridInterpolator(
            (self.d2_pressures, self.d2_lats, self.d2_lons), 
            self.d2_u_data, 
            bounds_error=False, fill_value=None
        )
        self.d2_v_interp = RegularGridInterpolator(
            (self.d2_pressures, self.d2_lats, self.d2_lons), 
            self.d2_v_data, 
            bounds_error=False, fill_value=None
        )
    def getWindLatLong(self, pressure, lat, lon, target_time):
        if lon < 0:
            lon = 360 + lon
        # makes sure pressure is in range
        pressure = max(min(pressure, self.d1_pressures.max()), self.d1_pressures.min())
        print(pressure)
        point = np.array([[pressure, lat, lon]])

        # two wind vecs for dataset before us in time and foward in time,
        u_val_ds1 = self.d1_u_interp(point)[0]
        v_val_ds1 = self.d1_v_interp(point)[0]

        u_val_ds2 = self.d2_u_interp(point)[0]
        v_val_ds2 = self.d2_v_interp(point)[0]

        # ratio through the 6 hour timeframe between the two datasets
        print(type(self.ds1_valid_time), self.ds1_valid_time)
        r_through = ((target_time - self.ds1_valid_time).total_seconds())/ (6 * 3600)
        print(u_val_ds1, v_val_ds1, u_val_ds2, v_val_ds2, r_through)
        #print(f"WIND1: {u_val_ds1}, {v_val_ds1} WIND2: {u_val_ds2}, {v_val_ds2}")

        interpolated_u = ((1-r_through) * u_val_ds1) + (r_through * u_val_ds2)
        interpolated_v = ((1-r_through) * v_val_ds1) + (r_through * v_val_ds2)


        return interpolated_u, interpolated_v
    def getWindEC(self, x, y, z, altitude, pressure, target_time, debug=False):
        if target_time > self.ds2_valid_time:
            self.update_datasets(target_time=target_time)
        lat, lon, alt = pm.ecef2geodetic(x, y, z)
        if lon<0:
            lon = 360+lon
        # pressure from pascals to millibars
        u_val, v_val = self.getWindLatLong(pressure/100, lat, lon, target_time)
        if debug:
            print(f"Calc wind at {lat, lon, altitude, pressure/100} is {u_val, v_val}")
        return pm.enu2ecefv(u_val, v_val, 0, lat, lon, altitude)
    

    def __init__(self, sim_start=datetime.now(tz=timezone.utc)):
        self.update_datasets(sim_start)

#W = WindPredictor()
#W.getWindLatLong(100, 50, 1, datetime.now(tz=timezone.utc))