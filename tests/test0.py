import matplotlib.pyplot as plt
import numpy as np

from coronagraph.coronagraph_paraxial import Coronagraph

# The Sun's angular size varies from 31' 27.7" to 32' 31.9"
# http://sun.stanford.edu/~sasha/PHYS780/SOLAR_PHYSICS/L2/Lecture_02_PHYS780.pdf
theta_sun = np.radians(32 / 60) / 2  # Angular size of the Sun in radians
theta_v0 = 1.2 * theta_sun  # Maximal angle of 100 % vignetting
theta_v1 = 15 * theta_sun  # Minimal angle of 0 % vignetting
theta_m = 10 * theta_sun  # Maximal angle of the field of view

cor = Coronagraph(Ra=13, theta_v0=theta_v0, theta_v1=theta_v1, theta_m=theta_m,
                    la=5, ld=5, lL=5, f1_=150, f2_=100, f3_=96.62)

print(f"Ra': {cor.Ra_:.4f} mm")
print(f"R0: {cor.R0:.2f} mm")
print(f"d0: {cor.d0:.4f} mm")
print(f"Rd: {cor.Rd:.4f} mm")
print(f"Rf: {cor.Rf:.4f} mm")
print(f"RL: {cor.RL:.4f} mm")
print(f"Rf': {cor.Rf_:.4f} mm")
print(f"Ri: {cor.Ri:.4f} mm")

print(f"F1: {cor.f1_:.2f} mm")
print(f"F2: {cor.f2_:.2f} mm")
print(f"F3: {cor.f3_:.2f} mm")

print(f"beta2: {cor.beta2:.4f}")
print(f"beta3: {cor.beta3:.4f}")

print(f"dd: {cor.dd:.4f} mm")
print(f"L1: {cor.L1:.4f} mm")
print(f"a2: {cor.a2:.4f} mm")
print(f"a2': {cor.a2_:.4f} mm")
print(f"la': {cor.la_:.4f} mm")
print(f"dL: {cor.dL:.4f} mm")
print(f"L2: {cor.L2:.4f} mm")
print(f"a3: {cor.a3:.4f} mm")
print(f"a3': {cor.a3_:.4f} mm")
print(f"L3: {cor.L3:.4f} mm")

print()
print("Thicknesses:")
print(f"t1: {cor.t1:.4f} mm")
print(f"t2: {cor.t2:.4f} mm")
print(f"t3: {cor.t3:.4f} mm")
print(f"t4: {cor.t4:.4f} mm")
print(f"t5: {cor.t5:.4f} mm")
print(f"t6: {cor.t6:.4f} mm")
print(f"t7: {cor.t7:.4f} mm")
print(f"t8: {cor.t8:.4f} mm")

print(f"EFFL of the system: {cor.f_c:.4f} mm")
print(f"f23_: {cor.f23_:.4f} mm")
print(f"p1F: {cor.p1F:.4f} mm")
print(f"p2_F_: {cor.p2_F_:.4f} mm")
print(f"p1P: {cor.p1P:.4f} mm")
print(f"p2_P_: {cor.p2_P_:.4f} mm")
print(f"p1P = p1F - f23: {cor.p1F - cor.f23:.4f} mm")
print(f"p2_P_ = p2_F_ - f23_: {cor.p2_F_ - cor.f23_:.4f} mm")
print(f"EFFL using the combination of lens 1 and lens 2+3: {cor.f_c2:.4f} mm")
print(f"Magnification of lens 2 + lens 3: {cor.beta23:.4f}")


fig, ax = plt.subplots(figsize=(12, 6))

one = np.array([1, 1])
sym = np.array([1, -1])
ax.plot(one * 0, sym * cor.R0, '-', label='External Occulter')
ax.plot(one * cor.t1, sym * cor.Ra, '-', label='Entrance Aperture')
ax.plot(one * (cor.t1 + cor.t2), sym * cor.Ra, '-', label='Objective Lens')
ax.plot(one * (cor.t1 + cor.t2 + cor.t3), sym * cor.Rf, '-', label='Field Stop')
ax.plot(one * (cor.t1 + cor.t2 + cor.t3 + cor.t4), sym * cor.Rd, 
        '-', label='Internal Occulting Disc')
ax.plot(one * (cor.t1 + cor.t2 + cor.t3 + cor.t4 + cor.t5), sym * cor.Ra, 
        '-', label='Field Lens')
ax.plot(one * (cor.t1 + cor.t2 + cor.t3 + cor.t4 + cor.t5 + cor.t6), sym * cor.RL, 
        '-', label='Lyot Stop')
ax.plot(one * (cor.t1 + cor.t2 + cor.t3 + cor.t4 + cor.t5 + cor.t6 + cor.t7), sym * cor.Ra, 
        '-', label='Relay Lens')
ax.plot(one * (cor.t1 + cor.t2 + cor.t3 + cor.t4 + cor.t5 + cor.t6 + cor.t7 + cor.t8), 
        sym * cor.Ri, 
        '-', label='Image Plane')
ax.set_xlabel('Distance (mm)')
ax.set_ylabel('Radius (mm)')
ax.set_title('Paraxial Coronagraph Layout')
ax.set_ylim(-100, 100)
ax.set_aspect('equal', adjustable='box')
fig.legend(ncol=4)
ax.grid()
#fig.tight_layout()
plt.show()
