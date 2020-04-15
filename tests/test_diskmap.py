import os

import pytest
import numpy as np

from astropy.io import fits

import diskmap


class TestDiskmap:

    def setup_class(self) -> None:

        np.random.seed(1)
        image = np.random.normal(loc=0, scale=1., size=(51, 51))
        fits.writeto('image.fits', image)

        self.limit = 1e-10

    def teardown_class(self) -> None:

        os.remove('image.fits')
        os.remove('diskmap_deprojected.fits')
        os.remove('diskmap_phase_function.dat')
        os.remove('diskmap_radius.fits')
        os.remove('diskmap_r2_scaled.fits')
        os.remove('diskmap_scat_angle.fits')
        os.remove('diskmap_total_intensity.fits')

    def test_diskmap(self) -> None:

        mapping = diskmap.DiskMap(fitsfile='image.fits',
                                  pixscale=1e-2,
                                  inclination=30.,
                                  pos_angle=70.,
                                  distance=100.)

        mapping.map_disk(power_law=(0., 0.01, 1.),
                         radius=(1., 200., 100),
                         surface='power-law',
                         filename=None)

        mapping.deproject_disk()

        mapping.r2_scaling(r_max=30.)

        mapping.total_intensity(pol_max=1.)

        mapping.phase_function(radius=(10., 20.),
                               n_phase=30)

        mapping.write_output(filename='diskmap')

        data = fits.getdata('diskmap_deprojected.fits')
        assert np.nansum(data) == pytest.approx(62.24650962576587, rel=self.limit, abs=0.)
        assert data.shape == (51, 51)

        data = fits.getdata('diskmap_radius.fits')
        assert np.sum(data) == pytest.approx(54749.27298678206, rel=self.limit, abs=0.)

        data = fits.getdata('diskmap_r2_scaled.fits')
        assert np.sum(data) == pytest.approx(27930.225422640284, rel=self.limit, abs=0.)

        data = fits.getdata('diskmap_scat_angle.fits')
        assert np.sum(data) == pytest.approx(232491.5391211969, rel=self.limit, abs=0.)

        data = fits.getdata('diskmap_total_intensity.fits')
        assert np.sum(data) == pytest.approx(38264.57782663277, rel=self.limit, abs=0.)

        data = np.loadtxt('diskmap_phase_function.dat')
        assert np.sum(data) == pytest.approx(998.046913483814, rel=self.limit, abs=0.)
        assert data.shape == (11, 5)
