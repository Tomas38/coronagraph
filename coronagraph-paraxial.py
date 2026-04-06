import numpy as np
import matplotlib.pyplot as plt


# Lens 1 = Objective lens
# Lens 2 = Field lens (not exactly field lens to be correct)
# Lens 3 = Relay lens

Ra = 13  # Entrance pupil/aperture radius

# The Sun's angular size varies from 31' 27.7" to 32' 31.9"
# http://sun.stanford.edu/~sasha/PHYS780/SOLAR_PHYSICS/L2/Lecture_02_PHYS780.pdf
theta_sun = np.radians(32 / 60) / 2  # Angular size of the Sun in radians
theta_v0 = 1.2 * theta_sun  # Maximal angle of 100 % vignetting
theta_v1 = 15 * theta_sun  # Minimal angle of 0 % vignetting
theta_m = 10 * theta_sun  # Maximal angle of the field of view

t = np.tan(theta_v0) / np.tan(theta_v1)
R0 = Ra * (1 + t) / (1 - t)  # External occulter radius
d0 = (R0 + Ra) / np.tan(theta_v1) # Distance between the external occulter and the entrance pupil/aperture

la = 5  # Distance between the entrance (a)perture/pupil and the objective lens
ld = 5  # Distance between the internal occulting (d)isc and the field lens
lL = 5  # Distance between the (L)yot stop and the relay lens
f1_ = 150  # Focal length of the objective lens
f2_ = 100  # Focal length of the field lens
f3_ = 96.62  # Focal length of the relay lens
Rf = f1_ * np.tan(theta_m)  # Field stop radius

a1 = -np.inf  # Object distance for lens 1
a1_ = f1_  # Image distance for lens 1

dd = 1 / (1 / a1_ + 1 / (-(d0 + la)))  # Distance from lens 1 to the occulting disc
L1 = dd + ld  # Distance between lens 1 and lens 2
a2 = -(L1 - a1_)  # Object distance for lens 2
a2_ = 1 / (1 / f2_ + 1 / a2)  # Image distance for lens 2

la_ = 1 / (1 / f1_ + 1 / (-la))  # Position of the virtual image of the entrance aperture through lens 1
dL = 1 / (1 / f2_ + 1 / (-(L1 - la_)))  # Distance from Lyot stop to lens 3
L2 = dL + lL  # Distance between lens 2 and lens 3
a3 = -(L2 - a2_)  # Object distance for lens 3

a3_ = 1 / (1 / f3_ + 1 / a3)  # Image distance for lens 3
L3 = a3_  # Distance between lens 3 and the image plane

beta2 = a2_ / a2  # Lateral magnification for lens 2
beta3 = a3_ / a3  # Lateral magnification for lens 3
beta23 = beta2 * beta3  # Total lateral magnification for the combination of lens 2 and lens 3

Rd = R0 * dd / (d0 + la)  # Radius of the internal occulting disc
Ra_ = Ra * (-la_) / (la)  # Radius of the (virtual) image of the entrance aperture through lens 1
RL = Ra_ * dL / (L1 - la_)  # Radius of the Lyot stop

Rf_ = Rf * (-beta2)  # Radius of the field stop
Ri = Rf_ * (-beta3)  # Radius of the field on the image plane

t1 = d0  # Distance from the external occulter to the entrance aperture
t2 = la  # Distance from the entrance aperture to lens 1
t3 = f1_  # Distance from lens 1 to the field stop
t4 = dd - f1_  # Distance from the field stop to the internal occulting disc
t5 = ld  # Distance from the occulting disc to lens 2
t6 = dL  # Distance from lens 2 to the Lyot stop
t7 = lL  # Distance from the Lyot stop to lens 3
t8 = L3  # Distance from lens 3 to the image plane

f_c = (a1_ * a2_ * a3_) / (a2 * a3)  # Total focal length

f23_ = (f2_ * f3_) / (f2_ + f3_ - L2)  # Focal length of the combination of lens 2 and lens 3
f23 = - f23_
p1F = - (f3_ - L2) * f2_ / (f2_ + f3_ - L2)  # Position of the object focal plane of the combination of lens 2 and lens 3 with respect to lens 2
p2_F_ = (f2_ - L2) * f3_ / (f2_ + f3_ - L2)  # Position of the image focal plane of the combination of lens 2 and lens 3 with respect to lens 3
p1P = + (f2_ * L2) / (f2_ + f3_ - L2)  # Position of the object principal plane of the combination of lens 2 and lens 3 with respect to lens 2
p2_P_ = - (f3_ * L2) / (f2_ + f3_ - L2)  # Position of the image principal plane of the combination of lens 2 and lens 3 with respect to lens 3

f_c2 = (f1_ * f23_) / (f1_ + f23_ - (L1 + p1P))  # Total focal length using the combination of lens 1 and lens 2

print(f"Ra': {Ra_:.4f} mm")
print(f"R0: {R0:.2f} mm")
print(f"d0: {d0:.4f} mm")
print(f"Rd: {Rd:.4f} mm")
print(f"Rf: {Rf:.4f} mm")
print(f"RL: {RL:.4f} mm")
print(f"Rf': {Rf_:.4f} mm")
print(f"Ri: {Ri:.4f} mm")

print(f"F1: {f1_:.2f} mm")
print(f"F2: {f2_:.2f} mm")
print(f"F3: {f3_:.2f} mm")

print(f"beta2: {beta2:.4f}")
print(f"beta3: {beta3:.4f}")

print(f"dd: {dd:.4f} mm")
print(f"L1: {L1:.4f} mm")
print(f"a2: {a2:.4f} mm")
print(f"a2': {a2_:.4f} mm")
print(f"la': {la_:.4f} mm")
print(f"dL: {dL:.4f} mm")
print(f"L2: {L2:.4f} mm")
print(f"a3: {a3:.4f} mm")
print(f"a3': {a3_:.4f} mm")
print(f"L3: {L3:.4f} mm")

print()
print("Thicknesses:")
print(f"t1: {t1:.4f} mm")
print(f"t2: {t2:.4f} mm")
print(f"t3: {t3:.4f} mm")
print(f"t4: {t4:.4f} mm")
print(f"t5: {t5:.4f} mm")
print(f"t6: {t6:.4f} mm")
print(f"t7: {t7:.4f} mm")
print(f"t8: {t8:.4f} mm")

print(f"EFFL of the system: {f_c:.4f} mm")
print(f"f23_: {f23_:.4f} mm")
print(f"p1F: {p1F:.4f} mm")
print(f"p2_F_: {p2_F_:.4f} mm")
print(f"p1P: {p1P:.4f} mm")
print(f"p2_P_: {p2_P_:.4f} mm")
print(f"p1P = p1F - f23: {p1F - f23:.4f} mm")
print(f"p2_P_ = p2_F_ - f23_: {p2_F_ - f23_:.4f} mm")
print(f"EFFL using the combination of lens 1 and lens 2+3: {f_c2:.4f} mm")
print(f"Magnification of lens 2 + lens 3: {beta23:.4f}")


fig, ax = plt.subplots(figsize=(12, 6))

one = np.array([1, 1])
sym = np.array([1, -1])
ax.plot(one * 0, sym * R0, '-', label='External Occulter')
ax.plot(one * t1, sym * Ra, '-', label='Entrance Aperture')
ax.plot(one * (t1 + t2), sym * Ra, '-', label='Objective Lens')
ax.plot(one * (t1 + t2 + t3), sym * Rf, '-', label='Field Stop')
ax.plot(one * (t1 + t2 + t3 + t4), sym * Rd, '-', label='Internal Occulting Disc')
ax.plot(one * (t1 + t2 + t3 + t4 + t5), sym * Ra, '-', label='Field Lens')
ax.plot(one * (t1 + t2 + t3 + t4 + t5 + t6), sym * RL, '-', label='Lyot Stop')
ax.plot(one * (t1 + t2 + t3 + t4 + t5 + t6 + t7), sym * Ra, '-', label='Relay Lens')
ax.plot(one * (t1 + t2 + t3 + t4 + t5 + t6 + t7 + t8), sym * Ri, '-', label='Image Plane')
ax.set_xlabel('Distance (mm)')
ax.set_ylabel('Radius (mm)')
ax.set_title('Paraxial Coronagraph Layout')
ax.set_ylim(-100, 100)
ax.set_aspect('equal', adjustable='box')
fig.legend(ncol=4)
ax.grid()
#fig.tight_layout()
plt.show()
