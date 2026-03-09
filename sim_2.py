from particle_model import superPressureWeatherBalloon
import datetime
import numpy as np
import pymap3d as pm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

lat = float(input("Enter latitude "))
lon = float(input("Enter longitude "))
p = superPressureWeatherBalloon(lat=lat, lon=lon)

dt = 0.1
positions = []
start_time = datetime.datetime.now()
simulation_time_passed = 0
simulation_time = 3600 * 24 * 10
prev_status_report = 0
while simulation_time_passed < simulation_time:
    if simulation_time_passed - prev_status_report > 2000:
        print(simulation_time_passed)
        print(np.linalg.norm(p.acceleration))
        
        print(p.getStatus())
        prev_status_report = simulation_time_passed
        positions.append(p.position.copy())

    # try changing time step for efficieny
    mag_acc = np.linalg.norm(p.acceleration)
    
    
    p.update(dt)
    simulation_time_passed += dt
    if p.hasBurst:
        break
    if mag_acc<0.01:
        dt = 2
    elif mag_acc <0.005:
        dt=20
    elif mag_acc<0.0005:
        dt=200
    else:
        dt = 0.1
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
