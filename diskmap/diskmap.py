"""
Module with mapping functionalities for protoplanetary disks.
"""

import math

from typing import Optional, Tuple

import numpy as np

from astropy.io import fits
from scipy.interpolate import griddata, interp1d
from typeguard import typechecked


class DiskMap:
    """
    Class for mapping a surface layer of a protoplanetary disk.
    """

    @typechecked
    def __init__(self,
                 fitsfile: str,
                 pixscale: float,
                 inclination: float,
                 pos_angle: float,
                 distance: float) -> None:
        """
        Parameters
        ----------
        fitsfile : str
            FITS file with the scattered light imaged.
        pixscale : float
            Pixel scale (arcsec).
        inclination : float
            Inclination of the disk (deg). Include a minus sign to exchange the near and far side
            in the mapping of the disk.
        pos_angle : float
            Position angle of the disk (deg). Defined in counterclockwise direction with respect
            to the vertical axis (i.e. east of north).
        distance : float
            Distance (pc).

        Returns
        -------
        NoneType
            None
        """

        self.image = fits.getdata(fitsfile)
        self.image = np.nan_to_num(self.image)

        if self.image.ndim != 2:
            raise ValueError('DiskMap requires a 2D image.')

        if self.image.shape[0] != self.image.shape[1]:
            raise ValueError('The dimensions of the image should have the same size.')

        self.pixscale = pixscale  # (arcsec)
        self.incl = math.radians(inclination)  # (rad)
        self.pos_ang = math.radians(pos_angle)  # (rad)
        self.distance = distance  # (pc)

        self.grid = 501  # should be odd

        self.radius = None
        self.azimuth = None
        self.opening = None
        self.scatter = None
        self.im_deproj = None
        self.im_scaled = None
        self.stokes_i = None
        self.phase = None

        # sum_before = np.sum(self.image)

        # im_scale = rescale(image=np.asarray(self.image, dtype=np.float64),
        #                    scale=(10., 10.),
        #                    order=5,
        #                    mode='reflect',
        #                    anti_aliasing=True,
        #                    multichannel=False)

        # sum_after = np.sum(im_scale)

        # self.image = im_scale * (sum_before / sum_after)
        # self.pixscale /= 10.

        self.npix = self.image.shape[0]

    @typechecked
    def map_disk(self,
                 power_law: Tuple[float, float, float],
                 radius: Tuple[float, float, int] = (1., 500., 100),
                 surface: str = 'power-law',
                 filename: Optional[str] = None) -> None:
        """
        Function for mapping a scattered light image to a power-law disk surface.

        Parameters
        ----------
        power_law : tuple(float, float, float)
            The argument for the power-law function, provided as (a, b, c) with
            f(x) = a + b*x^c, with ``a`` and ``b`` in au. Set all values to zero for the mapping
            and deprojection of a geometrically flat disk, in which case only the inclination is
            used for the deprojection.
        radius : tuple(float, float, int)
            Radius points that are sampled, provided as (r_in, r_out, n_r), with ``r_in`` and
            ``r_out`` in au. The outer radius should be set large enough such that a radius is
            sampled for each pixel in the field of view. To check if any NaNs are present, have
            a look at the `_radius.fits` output.
        surface : str
            Parameterization type for the disk surface ('power-law' or 'file').
        filename : star, None
            Filename which contains the radius in au (first column) and the height of the disk
            surface in au (second column).

        Returns
        -------
        NoneType
            None
        """

        # Create geometric disk model

        if surface == 'power-law':

            # Power-law disk height

            @typechecked
            def power_law_height(x_power: np.ndarray,
                                 a_power: float,
                                 b_power: float,
                                 c_power: float) -> np.ndarray:

                return a_power + b_power*x_power**c_power

            # midplane radius (au)
            disk_radius = np.linspace(radius[0], radius[1], radius[2])

            # disk height (au)
            disk_height = power_law_height(disk_radius, power_law[0], power_law[1], power_law[2])

            # opening angle (rad)
            disk_opening = np.arctan2(disk_height, disk_radius)

        elif surface == 'file':

            # Read disk height from ASCII file

            data = np.loadtxt(filename)

            # midplane radius (au)
            disk_radius = np.linspace(radius[0], radius[1], radius[2])

            # disk height (au)
            height_interp = interp1d(data[:, 0], data[:, 1])
            disk_height = height_interp(disk_radius)

            # opening angle (rad)
            disk_opening = np.arctan2(disk_height, disk_radius)  # (au)

        # Project disk height to image plane

        disk_phi = np.linspace(0., 359., 360)  # (deg)
        disk_phi = np.radians(disk_phi)  # (rad)

        x_im = []
        y_im = []
        r_im = []
        o_im = []
        s_im = []
        p_im = []

        for i, r_item in enumerate(disk_radius):
            for j, p_item in enumerate(disk_phi):

                x_tmp = r_item * np.sin(p_item)

                y_tmp = disk_height[i]*math.sin(self.incl) - \
                    r_item*np.cos(p_item)*math.cos(self.incl)

                x_rot = x_tmp*math.cos(math.pi-self.pos_ang) - \
                    y_tmp*math.sin(math.pi-self.pos_ang)

                y_rot = x_tmp*math.sin(math.pi-self.pos_ang) + \
                    y_tmp*math.cos(math.pi-self.pos_ang)

                x_im.append(x_rot)
                y_im.append(y_rot)

                r_im.append(math.sqrt(r_item**2+disk_height[i]**2))
                p_im.append(p_item)
                o_im.append(disk_opening[i])

                ang_tmp = math.pi/2.+disk_opening[i]

                par1 = math.sin(ang_tmp)*math.cos(math.pi+p_item)*math.sin(self.incl)
                par2 = math.cos(ang_tmp)*math.cos(self.incl)

                s_im.append(math.pi - math.acos(par1+par2))

        # Sort image plane points along x-axis

        x_index = np.argsort(x_im)

        y_sort = np.zeros(len(y_im))
        r_sort = np.zeros(len(y_im))
        p_sort = np.zeros(len(y_im))
        o_sort = np.zeros(len(y_im))
        s_sort = np.zeros(len(y_im))

        for i in range(len(y_im)):
            y_sort[i] = y_im[x_index[i]]
            r_sort[i] = r_im[x_index[i]]
            p_sort[i] = p_im[x_index[i]]
            o_sort[i] = o_im[x_index[i]]
            s_sort[i] = s_im[x_index[i]]

        grid_xy = np.zeros((self.grid**2, 2))

        count = 0

        for i in range(-(self.grid-1)//2, (self.grid-1)//2+1):
            for j in range(-(self.grid-1)//2, (self.grid-1)//2+1):
                grid_xy[count, 0] = float(i)
                grid_xy[count, 1] = float(j)

                count += 1

        image_xy = np.zeros((len(x_im), 2))

        for i, _ in enumerate(x_im):
            image_xy[i, 0] = x_im[i]
            image_xy[i, 1] = y_im[i]

        # Interpolate image plane

        if self.npix % 2 == 0:
            x_grid = np.linspace(-self.npix/2+0.5, self.npix/2-0.5, self.npix)
            y_grid = np.linspace(-self.npix/2+0.5, self.npix/2-0.5, self.npix)

        elif self.npix % 2 == 1:
            x_grid = np.linspace(-(self.npix-1)/2, (self.npix-1)/2, self.npix)
            y_grid = np.linspace(-(self.npix-1)/2, (self.npix-1)/2, self.npix)

        x_grid *= self.pixscale*self.distance  # (au)
        y_grid *= self.pixscale*self.distance  # (au)

        grid = np.zeros((self.npix**2, 2))

        count = 0

        for i in range(self.npix):
            for j in range(self.npix):
                grid[count, 0] = x_grid[i]
                grid[count, 1] = x_grid[j]

                count += 1

        fit_radius = griddata(image_xy, r_im, grid, method='linear')
        fit_azimuth = griddata(image_xy, p_im, grid, method='linear')
        fit_opening = griddata(image_xy, o_im, grid, method='linear')
        fit_scatter = griddata(image_xy, s_im, grid, method='linear')

        self.radius = np.zeros((self.npix, self.npix))
        self.azimuth = np.zeros((self.npix, self.npix))
        self.opening = np.zeros((self.npix, self.npix))
        self.scatter = np.zeros((self.npix, self.npix))

        count = 0

        for i in range(self.npix):
            for j in range(self.npix):
                self.radius[i, j] = fit_radius[count]
                self.azimuth[i, j] = fit_azimuth[count]
                self.opening[i, j] = fit_opening[count]
                self.scatter[i, j] = fit_scatter[count]

                count += 1

    @typechecked
    def deproject_disk(self) -> None:
        """
        Function for deprojecting a disk surface based on the mapping of ``map_disk``.

        Returns
        -------
        NoneType
            None
        """

        # Deproject disk surface

        if self.radius is None or self.azimuth is None:
            raise ValueError('Please run \'map_disk\' and before using \'deproject_disk\'.')

        x_disk = []
        y_disk = []
        im_disk = []

        for i in range(self.npix):
            for j in range(self.npix):
                x_tmp = self.radius[i, j]*np.cos(self.azimuth[i, j])
                y_tmp = self.radius[i, j]*np.sin(self.azimuth[i, j])

                x_disk.append(x_tmp*math.cos(math.pi/2.-self.pos_ang) -
                              y_tmp*math.sin(math.pi/2.-self.pos_ang))

                y_disk.append(x_tmp*math.sin(math.pi/2.-self.pos_ang) +
                              y_tmp*math.cos(math.pi/2.-self.pos_ang))

                im_disk.append(self.image[i, j])

        # Sort disk plane points along x-axis

        x_index = np.argsort(x_disk)

        y_sort = np.zeros(len(y_disk))
        im_sort = np.zeros(len(y_disk))

        for i in range(len(y_disk)):
            y_sort[i] = y_disk[x_index[i]]
            im_sort[i] = im_disk[x_index[i]]

        grid_xy = np.zeros((self.grid**2, 2))

        count = 0

        for i in range(-(self.grid-1)//2, (self.grid-1)//2+1):
            for j in range(-(self.grid-1)//2, (self.grid-1)//2+1):
                grid_xy[count, 0] = float(i)
                grid_xy[count, 1] = float(j)

                count += 1

        image_xy = np.zeros((len(y_disk), 2))

        for i, _ in enumerate(y_disk):
            image_xy[i, 0] = x_disk[i]
            image_xy[i, 1] = y_disk[i]

        # Interpolate disk plane

        if self.npix % 2 == 0:
            x_grid = np.linspace(-self.npix/2+0.5, self.npix/2-0.5, self.npix)
            y_grid = np.linspace(-self.npix/2+0.5, self.npix/2-0.5, self.npix)

        elif self.npix % 2 == 1:
            x_grid = np.linspace(-(self.npix-1)/2, (self.npix-1)/2, self.npix)
            y_grid = np.linspace(-(self.npix-1)/2, (self.npix-1)/2, self.npix)

        x_grid *= self.pixscale*self.distance  # (au)
        y_grid *= self.pixscale*self.distance  # (au)

        grid = np.zeros((self.npix**2, 2))

        count = 0

        for i in range(self.npix):
            for j in range(self.npix):
                grid[count, 0] = x_grid[i]
                grid[count, 1] = x_grid[j]

                count += 1

        try:
            fit_im = griddata(image_xy, im_disk, grid, method='linear')

        except ValueError:
            raise ValueError('The radius sampling should cover the complete field of view of the '
                             'image. Try increasing the outer \'radius\' value in \'map_disk\' '
                             'and have a look at the \'_radius.fits\' output to check for NaNs.')

        self.im_deproj = np.zeros((self.npix, self.npix))

        count = 0

        for i in range(self.npix):
            for j in range(self.npix):
                self.im_deproj[i, j] = fit_im[count]

                count += 1

    @typechecked
    def r2_scaling(self,
                   r_max: float) -> None:
        """
        Function for correcting a scattered light image for the r^2 decrease of the stellar
        irradiation of the disk surface.

        Parameters
        ----------
        r_max : float
            Maximum disk radius (au) for the r^2-scaling. Beyond this distance, a constant
            r^2-scaling is applied of value ``r_max``.

        Returns
        -------
        NoneType
            None
        """

        # r^2 scaling

        if self.radius is None:
            raise ValueError('Please run \'map_disk\' before using \'r2_scaling\'.')

        self.im_scaled = np.zeros((self.npix, self.npix))

        for i in range(self.npix):
            for j in range(self.npix):
                if self.radius[i, j] < r_max:
                    self.im_scaled[i, j] = self.radius[i, j]**2 * self.image[i, j]

                else:
                    self.im_scaled[i, j] = r_max**2 * self.image[i, j]

    @typechecked
    def total_intensity(self,
                        pol_max: 1.) -> None:
        """
        Function for estimating the total intensity image when ``fitsfile`` contains a polarized
        light image. A Rayleigh phase function is assumed and effects of multiple are ignored.

        Parameters
        ----------
        pol_max : float
            The peak of the Rayleigh phase function, which normalizes the estimated total
            intensity image.

        Returns
        -------
        NoneType
            None
        """

        # Total intensity

        if self.scatter is None or self.im_scaled is None:
            raise ValueError('Please run \'map_disk\' before using \'total_intensity\'.')

        alpha = np.cos(self.scatter)
        deg_pol = -pol_max*(alpha*alpha-1.)/(alpha*alpha+1.)

        self.stokes_i = self.im_scaled / deg_pol

    @typechecked
    def phase_function(self,
                       radius: Tuple[float, float],
                       n_phase: int):
        """
        Function for estimating the polarized and total intensity phase function when ``fitsfile``
        contains a polarized light image. A Rayleigh phase function is assumed and effects of
        multiple are ignored. Pixel values are scaled by r^2 such that irradiation effects do not
        affect the result, which is extracted in arbitrary units.

        Parameters
        ----------
        radius : tuple(float, float)
            Inner and outer radius (au) between which pixels are selected for estimating the phase
            function.
        n_phase : int
            Number of sampling points for the phase function between 0 and 180 deg.

        Returns
        -------
        NoneType
            None
        """

        # Phase function

        # Phase function is normalizedf so pol_max has not effect
        pol_max = 1.

        if self.radius is None or self.scatter is None or self.im_scaled is None:
            raise ValueError('Please run \'map_disk\' and \'r2_scaling\' before using '
                             '\'phase_function\'.')

        scat_select = []
        im_select = []

        for i in range(self.npix):
            for j in range(self.npix):
                if self.radius[i, j] > radius[0] and self.radius[i, j] < radius[1]:
                    scat_select.append(math.degrees(self.scatter[i, j]))

                    # use im_scaled to correct for differences across the selected radius
                    im_select.append(self.im_scaled[i, j])

        phase_bins = np.linspace(0., 180., num=n_phase)
        bin_index = np.digitize(scat_select, phase_bins)

        im_bins = []
        scat_bins = []

        for i in range(n_phase):
            im_bins.append([])
            scat_bins.append([])

        for i, _ in enumerate(im_select):
            im_bins[bin_index[i]-1].append(im_select[i])
            scat_bins[bin_index[i]-1].append(scat_select[i])

        angle = []
        pol_flux = []
        pol_error = []

        for i in range(n_phase):
            if len(im_bins[i]) > 0:
                angle.append(np.nanmean(scat_bins[i]))
                pol_flux.append(np.nanmean(im_bins[i]))
                pol_error.append(np.nanstd(im_bins[i])/math.sqrt(len(im_bins[i])))

        # Degree of polarization

        alpha = np.cos(np.array(angle)*np.pi/180.)
        deg_pol = -pol_max*(alpha*alpha-1.)/(alpha*alpha+1.)

        tot_flux = pol_flux/deg_pol
        tot_error = pol_error/deg_pol

        # Normalization

        flux_norm = max(pol_flux)
        pol_flux = np.array(pol_flux)/flux_norm
        pol_error = np.array(pol_error)/flux_norm

        flux_norm = max(tot_flux)
        tot_flux = np.array(tot_flux)/flux_norm
        tot_error = np.array(tot_error)/flux_norm

        self.phase = np.column_stack([angle, pol_flux, pol_error, tot_flux, tot_error])

    @typechecked
    def write_output(self,
                     filename: str) -> None:
        """
        Function for writing the available results to FITS files.

        Parameters
        ----------
        filename : str
            Filename start that is used for all the output file.

        Returns
        -------
        NoneType
            None
        """

        # Write FITS output

        if self.radius is not None:
            fits.writeto(f'{filename}_radius.fits', self.radius, overwrite=True)

        if self.scatter is not None:
            fits.writeto(f'{filename}_scat_angle.fits', np.degrees(self.scatter), overwrite=True)

        if self.im_deproj is not None:
            fits.writeto(f'{filename}_deprojected.fits', self.im_deproj, overwrite=True)

        if self.im_scaled is not None:
            fits.writeto(f'{filename}_r2_scaled.fits', self.im_scaled, overwrite=True)

        if self.stokes_i is not None:
            fits.writeto(f'{filename}_total_intensity.fits', self.stokes_i, overwrite=True)

        if self.phase is not None:
            header = 'Scattering angle (deg) - Polarized intensity (a.u.) - Uncertainty (a.u.) ' \
                     '- Total intensity (a.u.) - Uncertainty (a.u.)'

            np.savetxt(f'{filename}_phase_function.dat', self.phase, header=header)
