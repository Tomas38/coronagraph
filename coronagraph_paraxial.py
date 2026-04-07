import numpy as np
import matplotlib.pyplot as plt


class Coronagraph:
    # This class calculates axial distances of optical components in Lyot configuration with external occulter.
    def __init__(
            self,
            Ra,  # Entrance pupil/aperture radius
            theta_v0,  # Maximal angle of 100 % vignetting
            theta_v1,  # Minimal angle of 0 % vignetting
            theta_m,  # Maximal angle of the field of view
            la=0.0,  # Distance between the entrance (a)perture/pupil and the objective lens
            ld=0.0,  # Distance between the internal occulting (d)isc and the field lens
            lL=0.0,  # Distance between the (L)yot stop and the relay lens
            f1_=150.0,  # Focal length of the objective lens
            f2_=100.0,  # Focal length of the field lens
            f3_=96.62  # Focal length of the relay lens
            ):
        # Lens 1 = Objective lens
        # Lens 2 = Field lens (not exactly field lens to be correct)
        # Lens 3 = Relay lens

        # EO = External Occulter
        # AP = Entrance aperture/pupil
        # L1 = Lens 1
        # FP = Focal Plane (where the field stop is located)
        # IO = Internal Occulting disc
        # L2 = Lens 2
        # LS = Lyot Stop
        # L3 = Lens 3
        # IM = Image plane

        # (EO) <----d0----> (AP) <-la-> (L1) <-------------L1-------------> (L2) <--------L2------> (L3) <-L3-> (IM)
        # (EO) <----d0----> (AP) <-la-> (L1) <-------dd-------> (IO) <-ld-> (L2) <-dL-> (LS) <-lL-> (L3) <-L3-> (IM)
        # (EO) <----t1----> (AP) <-t2-> (L1) <-t3-> [FP] <-t4-> (IO) <-t5-> (L2) <-t6-> (LS) <-t7-> (L3) <-t8-> (IM)

        self.Ra = Ra
        self.theta_v0 = theta_v0
        self.theta_v1 = theta_v1
        self.theta_m = theta_m
        self.la = la
        self.ld = ld
        self.lL = lL
        self.f1_ = f1_
        self.f2_ = f2_
        self.f3_ = f3_

        t = np.tan(self.theta_v0) / np.tan(self.theta_v1)
        self.R0 = Ra * (1 + t) / (1 - t)  # External occulter radius
        self.d0 = (self.R0 + Ra) / np.tan(self.theta_v1) # Distance between the external occulter and the entrance pupil/aperture

        self.f1_max = self.d0 + self.la  # Maximal focal length of the objective lens to ensure that the image of the external occulter is formed behind the objective lens
        if self.f1_max <= self.f1_:
            raise ValueError("The focal length of the objective lens is larger than the distance between the external occulter and the objective lens.")
        
        self.Rf = self.f1_ * np.tan(self.theta_m)  # Field stop radius

        self.a1 = -np.inf  # Object distance for lens 1
        self.a1_ = self.f1_  # Image distance for lens 1

        self.dd = 1 / (1 / self.a1_ + 1 / (-(self.d0 + self.la)))  # Distance from lens 1 to the occulting disc
        self.L1 = self.dd + self.ld  # Distance between lens 1 and lens 2

        self.a2 = -(self.L1 - self.a1_)  # Object distance for lens 2
        self.a2_ = 1 / (1 / self.f2_ + 1 / self.a2)  # Image distance for lens 2

        self.la_ = 1 / (1 / self.f1_ + 1 / (-self.la))  # Position of the virtual image of the entrance aperture through lens 1

        self.f2_max = self.L1 - self.la_  # Maximal focal length of the field lens to ensure that the image of the entrance aperture through lens 2 is formed behind the lens 2
        if self.f2_max <= self.f2_:
            # equivalent to self.dL < 0.0 condition
            raise ValueError("The focal length of the field lens is larger than the distance between the virtual image of the entrance aperture and lens 2.")

        self.dL = 1 / (1 / self.f2_ + 1 / (-(self.L1 - self.la_)))  # Distance from Lyot stop to lens 3
        self.L2 = self.dL + self.lL  # Distance between lens 2 and lens 3
        self.a3 = -(self.L2 - self.a2_)  # Object distance for lens 3

        self.a3_ = 1 / (1 / self.f3_ + 1 / self.a3)  # Image distance for lens 3
        self.f3_max = self.L2 - self.a2_  # Maximal focal length of the relay lens to ensure that the image of the entrance aperture through lens 3 is formed behind the lens 3
        if self.f3_max <= self.f3_:
            raise ValueError("The focal length of the relay lens is larger than the distance of the virtual object to be displayed and the relay lens.")
        self.L3 = self.a3_  # Distance between lens 3 and the image plane

        self.beta2 = self.a2_ / self.a2  # Lateral magnification for lens 2
        self.beta3 = self.a3_ / self.a3  # Lateral magnification for lens 3
        self.beta23 = self.beta2 * self.beta3  # Total lateral magnification for the combination of lens 2 and lens 3

        self.Rd = self.R0 * self.dd / (self.d0 + self.la)  # Radius of the internal occulting disc
        self.Ra_ = self.Ra * (-self.la_) / (self.la)  # Radius of the (virtual) image of the entrance aperture through lens 1
        self.RL = self.Ra_ * self.dL / (self.L1 - self.la_)  # Radius of the Lyot stop

        self.Rf_ = self.Rf * (-self.beta2)  # Radius of the field stop
        self.Ri = self.Rf_ * (-self.beta3)  # Radius of the field on the image plane

        self.t1 = self.d0  # Distance from the external occulter to the entrance aperture
        self.t2 = self.la  # Distance from the entrance aperture to lens 1
        self.t3 = self.f1_  # Distance from lens 1 to the field stop
        self.t4 = self.dd - self.f1_  # Distance from the field stop to the internal occulting disc
        self.t5 = self.ld  # Distance from the occulting disc to lens 2
        self.t6 = self.dL  # Distance from lens 2 to the Lyot stop
        self.t7 = self.lL  # Distance from the Lyot stop to lens 3
        self.t8 = self.L3  # Distance from lens 3 to the image plane

        self.f_c = (self.a1_ * self.a2_ * self.a3_) / (self.a2 * self.a3)  # Total focal length

        self.f23_ = (self.f2_ * self.f3_) / (self.f2_ + self.f3_ - self.L2)  # Focal length of the combination of lens 2 and lens 3
        self.f23 = - self.f23_
        self.p1F = - (self.f3_ - self.L2) * self.f2_ / (self.f2_ + self.f3_ - self.L2)  # Position of the object focal plane of the combination of lens 2 and lens 3 with respect to lens 2
        self.p2_F_ = (self.f2_ - self.L2) * self.f3_ / (self.f2_ + self.f3_ - self.L2)  # Position of the image focal plane of the combination of lens 2 and lens 3 with respect to lens 3
        self.p1P = + (self.f2_ * self.L2) / (self.f2_ + self.f3_ - self.L2)  # Position of the object principal plane of the combination of lens 2 and lens 3 with respect to lens 2
        self.p2_P_ = - (self.f3_ * self.L2) / (self.f2_ + self.f3_ - self.L2)  # Position of the image principal plane of the combination of lens 2 and lens 3 with respect to lens 3

        self.f_c2 = (self.f1_ * self.f23_) / (self.f1_ + self.f23_ - (self.L1 + self.p1P))  # Total focal length using the combination of lens 1 and lens 2


if __name__ == "__main__":
    # The Sun's angular size varies from 31' 27.7" to 32' 31.9"
    # http://sun.stanford.edu/~sasha/PHYS780/SOLAR_PHYSICS/L2/Lecture_02_PHYS780.pdf
    theta_sun = np.radians(32 / 60) / 2  # Angular size of the Sun in radians
    theta_v0 = 1.2 * theta_sun  # Maximal angle of 100 % vignetting
    theta_v1 = 15 * theta_sun  # Minimal angle of 0 % vignetting
    theta_m = 10 * theta_sun  # Maximal angle of the field of view

    cor = Coronagraph(Ra=13, theta_v0=theta_v0, theta_v1=theta_v1, theta_m=theta_m, la=5, ld=5, lL=5, f1_=150, f2_=100, f3_=96.62)

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
    ax.plot(one * (cor.t1 + cor.t2 + cor.t3 + cor.t4), sym * cor.Rd, '-', label='Internal Occulting Disc')
    ax.plot(one * (cor.t1 + cor.t2 + cor.t3 + cor.t4 + cor.t5), sym * cor.Ra, '-', label='Field Lens')
    ax.plot(one * (cor.t1 + cor.t2 + cor.t3 + cor.t4 + cor.t5 + cor.t6), sym * cor.RL, '-', label='Lyot Stop')
    ax.plot(one * (cor.t1 + cor.t2 + cor.t3 + cor.t4 + cor.t5 + cor.t6 + cor.t7), sym * cor.Ra, '-', label='Relay Lens')
    ax.plot(one * (cor.t1 + cor.t2 + cor.t3 + cor.t4 + cor.t5 + cor.t6 + cor.t7 + cor.t8), sym * cor.Ri, '-', label='Image Plane')
    ax.set_xlabel('Distance (mm)')
    ax.set_ylabel('Radius (mm)')
    ax.set_title('Paraxial Coronagraph Layout')
    ax.set_ylim(-100, 100)
    ax.set_aspect('equal', adjustable='box')
    fig.legend(ncol=4)
    ax.grid()
    #fig.tight_layout()
    plt.show()
