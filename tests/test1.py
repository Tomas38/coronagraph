import numpy as np
from coronagraph import Coronagraph


theta_sun = np.radians(32 / 60) / 2  # Angular size of the Sun in radians
theta_v0 = 1.2 * theta_sun  # Maximal angle of 100 % vignetting
theta_v1 = 15 * theta_sun  # Minimal angle of 0 % vignetting
theta_m = 10 * theta_sun  # Maximal angle of the field of view

cor = Coronagraph(Ra=13, theta_v0=theta_v0, theta_v1=theta_v1, theta_m=theta_m, la=5, ld=5, lL=5, f1_=150, f2_=100, f3_=96.62)
