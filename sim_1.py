from particle_model import WeatherBalloon
import datetime
import numpy as np
import pymap3d as pm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

lat = float(input("Enter latitude "))
lon = float(input("Enter longitude "))
p = WeatherBalloon(lat=lat, lon=lon)

dt = 0.1
positions = []
start_time = datetime.datetime.now()
simulation_time_passed = 0
for i in range(50000):
    if i % 500 == 0:
        print(i)
        print(np.linalg.norm(p.acceleration))
        
        print(p.getStatus())

    # try changing time step for efficieny
    if np.linalg.norm(p.acceleration)<0.01:
        dt = 2
    else:
        dt = 0.1
    p.update(dt)
    simulation_time_passed += dt
    if p.hasBurst:
        break
    positions.append(p.position.copy())
print(f"SIMULATION TIME PASSED: {simulation_time_passed}")
print(f"EXECUTION TIME: {datetime.datetime.now() - start_time}")

lats, lons, alts = [], [], []
for pos in positions:
    lat, lon, alt = pm.ecef2geodetic(*pos)
    lats.append(lat)
    lons.append(lon)
    alts.append(alt)

lats = np.array(lats)
lons = np.array(lons)
alts = np.array(alts)

fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

ax.plot(lons, lats, alts, color='blue', label='Balloon path')
ax.scatter(lons[0], lats[0], alts[0], color='green', label='Start')
ax.scatter(lons[-1], lats[-1], alts[-1], color='red', label='End')

ax.set_xlabel('Longitude (deg)')
ax.set_ylabel('Latitude (deg)')
ax.set_zlabel('Altitude (m)')
ax.legend()
plt.show()
