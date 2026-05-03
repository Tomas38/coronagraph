import numpy as np

from coronagraph import Coronagraph

theta_sun = np.radians(32 / 60) / 2  # Angular size of the Sun in radians
theta_v0 = 1.2 * theta_sun  # Maximal angle of 100 % vignetting
theta_v1 = 10 * theta_sun  # Minimal angle of 0 % vignetting
theta_m = 15 * theta_sun  # Maximal angle of the field of view

cor = Coronagraph(Ra=13,
                  theta_v0=theta_v0, theta_v1=theta_v1, theta_m=theta_m,
                  la=5, ld=5, lL=5,
                  f1_=150, f2_=100, f3_=96.62)


def test_vignetting_scalar_boundaries():
    assert cor.vignetting(theta_v0 * 0.5) == 1.0
    assert cor.vignetting(theta_v1 * 1.5) == 0.0


def test_vignetting_exact_boundaries():
    assert np.isclose(cor.vignetting(theta_v0), 1.0, rtol=0.0, atol=1e-6)
    assert np.isclose(cor.vignetting(theta_v1), 0.0, rtol=0.0, atol=1e-6)


def test_vignetting_vectorized_matches_scalar():
    theta = np.array([
        theta_v0 * 0.5,
        0.5 * (theta_v0 + theta_v1),
        theta_v1 * 1.5,
    ])
    vector_result = cor.vignetting(theta)
    scalar_result = np.array([cor.vignetting(float(t)) for t in theta])

    assert isinstance(vector_result, np.ndarray)
    np.testing.assert_allclose(vector_result, scalar_result, rtol=1e-12, atol=1e-12)


def test_vignetting_vectorized_shape_and_range():
    theta = np.linspace(0, theta_m, 1000)
    vignetting = cor.vignetting(theta)

    assert vignetting.shape == theta.shape
    assert np.all(vignetting >= 0.0)
    assert np.all(vignetting <= 1.0)

# To run the tests, use the following command in the terminal:
# python -m pytest tests/test_vignetting.py -q
