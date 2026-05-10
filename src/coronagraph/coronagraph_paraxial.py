from typing import Any, overload

import numpy as np
import numpy.typing as npt


class Coronagraph:
    def __init__(self, Ra,
                 theta_v0, theta_v1,
                 theta_m, la=0.0,
                 ld=0.0, lL=0.0,
                 f1_=150.0, f2_=100.0, f3_=96.62):
        """This class calculates axial distances of optical components
        in Lyot configuration with external occulter.

        Parameters
        ----------
        Ra : _type_
            Entrance pupil/aperture radius
        theta_v0 : _type_
            Maximal angle of 100 % vignetting
        theta_v1 : _type_
            Minimal angle of 0 % vignetting
        theta_m : _type_
            Maximal angle of the field of view
        la : float, optional
            Distance between the entrance (a)perture/pupil and the objective lens, by default 0.0
        ld : float, optional
            Distance between the internal occulting (d)isc and the field lens, by default 0.0
        lL : float, optional
            Distance between the (L)yot stop and the relay lens, by default 0.0
        f1_ : float, optional
            Focal length of the objective lens, by default 150.0
        f2_ : float, optional
            Focal length of the field lens, by default 100.0
        f3_ : float, optional
            Focal length of the relay lens, by default 96.62

        Raises
        ------
        ValueError
            If the focal length of the objective lens is larger than the distance
            between the external occulter and the objective lens.
        ValueError
            If the focal length of the field lens is larger than the distance
            between the virtual image of the entrance aperture and lens 2.
        ValueError
            If the focal length of the relay lens is larger than the distance
            of the virtual object to be displayed and the relay lens.
        """
        
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

        # (EO) <--d0--> (AP) <la> (L1) <-----L1-----> (L2) <-----------L2---------> (L3) <L3> (IM)
        # (EO) <--d0--> (AP) <la> (L1) <-----dd-----> (IO) <ld> (L2) <dL> (LS) <lL> (L3) <L3> (IM)
        # (EO) <--t1--> (AP) <t2> (L1) <t3> [FP] <t4> (IO) <t5> (L2) <t6> (LS) <t7> (L3) <t8> (IM)

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
        # Distance between the external occulter and the entrance pupil/aperture
        self.d0 = (self.R0 + Ra) / np.tan(self.theta_v1)

        # Maximal focal length of the objective lens to ensure that the image
        # of the external occulter is formed behind the objective lens
        self.f1_max = self.d0 + self.la
        if self.f1_max <= self.f1_:
            raise ValueError("The focal length of the objective lens is larger " \
            "than the distance between the external occulter and the objective lens.")
        
        self.Rf = self.f1_ * np.tan(self.theta_m)  # Field stop radius

        self.a1 = -np.inf  # Object distance for lens 1
        self.a1_ = self.f1_  # Image distance for lens 1

        # Distance from lens 1 to the occulting disc
        self.dd = 1 / (1 / self.a1_ + 1 / (-(self.d0 + self.la)))
        self.L1 = self.dd + self.ld  # Distance between lens 1 and lens 2

        self.a2 = -(self.L1 - self.a1_)  # Object distance for lens 2
        self.a2_ = 1 / (1 / self.f2_ + 1 / self.a2)  # Image distance for lens 2

        if self.la == 0.0:
            self.la_ = 0.0
        else:
            # Position of the virtual image of the entrance aperture through lens 1
            self.la_ = 1 / (1 / self.f1_ + 1 / (-self.la))

        # Maximal focal length of the field lens to ensure that the image
        # of the entrance aperture through lens 2 is formed behind the lens 2
        self.f2_max = self.L1 - self.la_
        if self.f2_max <= self.f2_:
            # equivalent to self.dL < 0.0 condition
            raise ValueError("The focal length of the field lens is larger " \
            "than the distance between the virtual image of the entrance aperture and lens 2.")
        
        # Distance from Lyot stop to lens 3
        self.dL = 1 / (1 / self.f2_ + 1 / (-(self.L1 - self.la_)))
        self.L2 = self.dL + self.lL  # Distance between lens 2 and lens 3
        self.a3 = -(self.L2 - self.a2_)  # Object distance for lens 3

        self.a3_ = 1 / (1 / self.f3_ + 1 / self.a3)  # Image distance for lens 3
        # Maximal focal length of the relay lens to ensure that the image of
        # the entrance aperture through lens 3 is formed behind the lens 3
        self.f3_max = self.L2 - self.a2_
        if self.f3_max <= self.f3_:
            raise ValueError("The focal length of the relay lens is larger than " \
            "the distance of the virtual object to be displayed and the relay lens.")
        self.L3 = self.a3_  # Distance between lens 3 and the image plane

        self.beta2 = self.a2_ / self.a2  # Lateral magnification for lens 2
        self.beta3 = self.a3_ / self.a3  # Lateral magnification for lens 3
        # Total lateral magnification for the combination of lens 2 and lens 3
        self.beta23 = self.beta2 * self.beta3

        self.Rd = self.R0 * self.dd / (self.d0 + self.la)  # Radius of the internal occulting disc
        if self.la == 0.0:
            self.Ra_ = self.Ra
        else:
            # Radius of the (virtual) image of the entrance aperture through lens 1
            self.Ra_ = self.Ra * (-self.la_) / (self.la)
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

        self.R1 = self.Ra + self.la * np.tan(self.theta_m)  # Radius (minimal) of the Lens 1
        r2 = (((self.t4 + self.t5) / (self.t3)) * 
              (self.Rf + self.Ra - self.la * np.tan(self.theta_m)))
        self.R2 = self.Rf + r2 # Radius (minimal) of the Lens 2

        self.f_c = (self.a1_ * self.a2_ * self.a3_) / (self.a2 * self.a3)  # Total focal length

        # Focal length of the combination of lens 2 and lens 3
        self.f23_ = (self.f2_ * self.f3_) / (self.f2_ + self.f3_ - self.L2)
        self.f23 = - self.f23_
        # Position of the object focal plane of the combination of 
        # lens 2 and lens 3 with respect to lens 2
        self.p1F = - (self.f3_ - self.L2) * self.f2_ / (self.f2_ + self.f3_ - self.L2)
        # Position of the image focal plane of the combination of 
        # lens 2 and lens 3 with respect to lens 3
        self.p2_F_ = (self.f2_ - self.L2) * self.f3_ / (self.f2_ + self.f3_ - self.L2)
        # Position of the object principal plane of the combination 
        # of lens 2 and lens 3 with respect to lens 2
        self.p1P = + (self.f2_ * self.L2) / (self.f2_ + self.f3_ - self.L2)
        # Position of the image principal plane of the combination 
        # of lens 2 and lens 3 with respect to lens 3
        self.p2_P_ = - (self.f3_ * self.L2) / (self.f2_ + self.f3_ - self.L2)

        # Total focal length using the combination of lens 1 and lens 2
        self.f_c2 = (self.f1_ * self.f23_) / (self.f1_ + self.f23_ - (self.L1 + self.p1P))

    @overload
    def vignetting(self, theta: float) -> float:
        ...

    @overload
    def vignetting(self, theta: npt.NDArray[np.floating[Any]]) -> npt.NDArray[np.floating[Any]]:
        ...

    def vignetting(
        self,
        theta: float | npt.NDArray[np.floating[Any]],
    ) -> float | npt.NDArray[np.floating[Any]]:
        """Calculate the vignetting for a given angle. Returns 1 for angles larger than theta_m.

        Parameters
        ----------
        theta : float or numpy.ndarray
            Angle(s) in radians

        Returns
        -------
        float or numpy.ndarray
            Vignetting factor(s) between 0 and 1
        """
        theta_arr = np.asarray(theta)
        is_scalar_input = np.ndim(theta_arr) == 0
        theta_arr = np.atleast_1d(theta_arr).astype(float)

        vign = np.empty_like(theta_arr, dtype=float)
        mask_full = theta_arr < self.theta_v0
        mask_none = theta_arr > self.theta_v1
        mask_large = theta_arr > self.theta_m
        mask_partial = ~(mask_full | mask_none)

        vign[mask_full] = 1.0
        vign[mask_none] = 0.0
        vign[mask_large] = 1.0

        if np.any(mask_partial):
            x_theta = self.d0 * np.tan(theta_arr[mask_partial])

            alpha_arg = (self.R0**2 + x_theta**2 - self.Ra**2) / (2 * self.R0 * x_theta)
            beta_arg = (x_theta**2 + self.Ra**2 - self.R0**2) / (2 * x_theta * self.Ra)
            alpha = np.arccos(np.clip(alpha_arg, -1.0, 1.0))
            beta = np.arccos(np.clip(beta_arg, -1.0, 1.0))

            area_aper_vignetted = ((np.pi - beta) * self.Ra**2
                                   + x_theta * self.R0 * np.sin(alpha)
                                   - alpha * self.R0**2)
            area_aper_unvignetted = np.pi * self.Ra**2
            vign[mask_partial] = 1.0 - area_aper_vignetted / area_aper_unvignetted

        if is_scalar_input:
            return float(vign[0])
        return vign
