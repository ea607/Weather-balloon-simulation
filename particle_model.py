import numpy as np
import math
# atmosphere model from https://www.pdas.com/atmos.html
from atmos_model import SimpleAtmosphere, TZERO, PZERO, RHOZERO
from aigfs_extracter2 import WindPredictor
import pymap3d as pm
from datetime import datetime, timezone, timedelta

EARTH_MASS = 5.9722 * 10**24
GRAVITATIONAL_CONST = 6.6743 * 10**-11

EARTH_CENT_POS = np.array([0, 0, 0], dtype=float)
EARTH_RADIUS = 6371e3
SEA_LEVEL_PRESSURE = 101325

windPredictor = WindPredictor()

class Particle:
    running_id = 0
    def __init__(self):
        self.position = np.array([0.0, 6371 * 1000, 0.0])
        self.velocity = np.array([0.0, 0.0, 0.0])
        self.acceleration = np.array([0.0, 0.0, 0.0])
        self.mass = 1
        self.forces = {}
        self.id = self.running_id
        Particle.running_id+=1
    def setPosition(self, position):
        self.position = np.array(position)
    def setVelocity(self, velocity):
        self.velocity = np.array(velocity)
    def setAcceleration(self, acceleration):
        self.acceleration = np.array(acceleration)
    
    # force can be 3d array or a function returning 3d array
    def setForce(self, label, force=[0.0, 0.0, 0.0]):
        if callable(force):
            self.forces[label] = force
        else:
            self.forces[label] = np.array(force, dtype=float)
    def removeForce(self, label):
        self.forces.pop(label, None)

    # abstract method for subclasses
    def onUpdate(self, dt):
        pass

    def update(self, dt=1):
        self.onUpdate(dt)
        resultant = np.array([0.0, 0.0, 0.0])
        for label, force in self.forces.items():
            if callable(force):
                force = force()
            resultant += force
        self.acceleration = resultant / self.mass
        self.velocity += (self.acceleration * dt)
        self.position += (self.velocity * dt)
    def getStatus(self):
        return f"P#{self.id}, POS={self.position.round(3).tolist()}, VEL={self.velocity.round(3).tolist()}, ACC={self.acceleration.round(3).tolist()}"
    
class WeatherBalloon(Particle):
    def __init__(self, lat=90, lon=0, mass=0.4, radius=0.5, burst_radius=1.5, sim_time=datetime.now(tz=timezone.utc)):
        super().__init__()

        self.mass = mass
        self.radius = radius
        # radius at bursting point
        self.burst_radius = burst_radius
        self.sim_time = sim_time

        self.volume = (4/3) * math.pi * (self.radius **3)
        self.initial_radius = self.radius
        self.initial_volume = self.volume
        self.burst_volume = (4/3) * math.pi * (self.burst_radius **3)
        self.drag_coeffient = 0.35
        self.area = math.pi* (self.radius ** 2)

        self.hasBurst = False
        
        # wind data comes from large database so dont recalculate if we havent moved signifcantly
        self._wind_cache = np.zeros(3)
        self._last_wind_pos = np.array(self.position)
        self._wind_threshold = 500
        self._alt = 0
        self._lat = lat
        self._lon = lon
        # location is north pole
        x, y, z = pm.geodetic2ecef(lat, lon, 0)
        print([x,y,z])
        self.setPosition([x,y,z])
        self.setForce("weight", self.calculateWeight)
        self.setForce("lift", self.calculateLift)
        self.setForce("drag", self.calculateDrag)
    def calculateG(self):
        
        dist_from_centre = np.linalg.norm(EARTH_CENT_POS - self.position)
        # this is magnitude of force, we need vector acting towards earth centre
        force_mag = GRAVITATIONAL_CONST * ((EARTH_MASS) / dist_from_centre**2)
        # unit vector direction
        force_unit_vec = (EARTH_CENT_POS - self.position) / np.linalg.norm(EARTH_CENT_POS - self.position)
        return force_mag * force_unit_vec
    def calculateWeight(self):
        return self.calculateG() * self.mass

    def calculateStatsAtAlt(self, altitude):
        density_rel, pressure_rel, temp_rel = SimpleAtmosphere(altitude/1000)
        return (density_rel * RHOZERO, pressure_rel*PZERO, temp_rel* TZERO)
    def calculatePressureAtAlt():
        pass
    def calculateLift(self):
        # calculate g at current altitude
        g = self.calculateG()
        # calculate density of air at our current altitude
        air_density, pressure, temp = self.calculateStatsAtAlt(self._alt)
        # calculate weight of displaced air
        # to get lift direction is opposite to g
        lift = (air_density * self.volume) * g * -1
        return lift
    def calculateDrag(self):
        air_density, pressure, temp = self.calculateStatsAtAlt(self._alt)
        # https://www.grc.nasa.gov/www/k-12/VirtualAero/BottleRocket/airplane/drageq.html

        wind_vec = self.calculateWindVec()
        relative_velocity = self.velocity - wind_vec
        speed = np.linalg.norm(relative_velocity)
        if speed == 0:
            return np.zeros(3)
        opposite_vel_unit = (relative_velocity / speed) * -1
        drag = self.drag_coeffient * air_density * ((speed**2)/2) * self.area * opposite_vel_unit
        return drag
    def calculateWindVec(self):
        dist_moved = np.linalg.norm(self.position - self._last_wind_pos)
        if dist_moved > self._wind_threshold:
            air_density, pressure, temp = self.calculateStatsAtAlt(self._alt)
            x, y, z = self.position
            self._wind_cache = windPredictor.getWindEC(x, y, z, self._alt, pressure=pressure, target_time=self.sim_time)
            self._last_wind_pos = np.array(self.position)
        return self._wind_cache
    def updateVolume(self):
        #print("ALT!!!!", self.calculateAltitude())
        air_density, pressure, temp = self.calculateStatsAtAlt(self._alt)
        #print("PRESSURE!!!!", pressure)
        self.volume = (self.initial_volume) * (SEA_LEVEL_PRESSURE / pressure) * (temp/TZERO)
        
        if self.volume >= self.burst_volume:
            self.burst()
    def burst(self):
        print(f"BALLOON BURST at ALT={self._alt}")
        self.hasBurst = True
    def onUpdate(self, dt):
        self.sim_time += timedelta(seconds=dt)
        self._lat, self._lon, self._alt = pm.ecef2geodetic(*self.position)
        self.updateVolume()
    def getStatus(self):
        output = "WEATHER BALLOON\n"
        output += f"P#{self.id}, POS={self.position.round(3).tolist()}, VEL={self.velocity.round(3).tolist()}, ACC={self.acceleration.round(3).tolist()}"
        output += f"\nALT={self._alt}, WIND={self.calculateWindVec()}\n"
        output += f"VOL={self.volume}"
        output += f"\nLAT, LON={self._lat}, {self._lon}"
        return output


class superPressureWeatherBalloon(WeatherBalloon):
    # super pressure balloon has constant volume
    def updateVolume(self):
        pass