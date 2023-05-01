"""
Module with mapping functionalities for protoplanetary disks.
"""

import math
import warnings

from typing import List, Optional, Tuple, Union, Callable

import numpy as np
import numpy.ma as ma

from astropy.io import fits
from astropy.stats import sigma_clip
from scipy.interpolate import griddata, interp1d
from scipy.ndimage import gaussian_filter, median_filter
from typeguard import typechecked


class DiskMap:
    """
    Class for mapping a surface layer of a protoplanetary disk.
    """

    @typechecked
    def __init__(self,
                 fitsfile: Union[str, np.ndarray],
                 pixscale: float,
                 inclination: float,
                 pos_angle: float,
                 distance: float,
                 image_type: str = 'polarized') -> None:
        """
        Parameters
        ----------
        fitsfile : str, np.ndarray
            Name of the FITS file with the scattered light image.
            Alternatively, a 2D ``numpy`` array with the image can
            be directly provided.
        pixscale : float
            Pixel scale of the image (arcsec per pixel).
        inclination : float
            Inclination of the disk (deg). The convention is such that
            the near side of the disk is on the right side of the image
            when using an inclination between 0 and 90 deg and using a
            `pos_angle` of 0 deg. The near and far side of the disk
            mapping can be exchanged by using a minus sign for the
            inclination. To be certain about using the correct near and
            far side, it is best to check the `_radius.fits` file with
            a somewhat large inclination. The near side will show more
            strongly compressed radii compared to the far side of the
            disk.
        pos_angle : float
            Position angle of the disk (deg). Defined in
            counterclockwise direction with respect to the vertical
            axis (i.e. east of north).
        distance : float
            Distance between observer and star (pc).
        image_type : str
            Image type ('polarized' or 'total'). This parameter affects
            the output that will be stored. For example, the conversion
            from polarized to total intensity phase function is only
            done with `image_type='polarized'`.

        Returns
        -------
        NoneType
            None
        """

        if isinstance(fitsfile, str):
            self.image = fits.getdata(fitsfile)
        else:
            self.image = fitsfile

        self.image = np.nan_to_num(self.image)

        if self.image.ndim == 3:
            warnings.warn(
                "The FITS file contains a 3D data cube so using the first image."
            )

            self.image = self.image[
                0,
            ]

        elif self.image.ndim != 2:
            raise ValueError("DiskMap requires a 2D image.")

        if self.image.shape[0] != self.image.shape[1]:
            raise ValueError("The dimensions of the image should have the same size.")
        
        if self.image.dtype == np.float32 or self.image.dtype == np.dtype('>f4'):
            warnings.warn(
                "The FITS file data is of type float32, this will be converted to float64"
            )
            self.image = self.image.astype(np.float64)
            
        if self.image.dtype != np.float64 and self.image.dtype != np.dtype('>f8'):
            raise ValueError(
                f"The FITS file data should be either of type float32 or float64"
            )

        if image_type not in ["polarized", "total"]:
            raise ValueError(
                "The argument of 'image_type' should be set "
                "to 'polarized' or 'total'."
            )

        self.pixscale = pixscale  # (arcsec)
        self.incl = math.radians(inclination)  # (rad)
        self.pos_ang = math.radians(pos_angle)  # (rad)
        self.distance = distance  # (pc)
        self.image_type = image_type

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
    def map_disk(
        self,
        power_law: Tuple[float, float, float],
        radius: Tuple[float, float, int] = (1.0, 500.0, 100),
        surface: str = "power-law",
        height_func: Optional[Callable[[np.ndarray],np.ndarray]] = None,
        filename: Optional[str] = None,
    ) -> None:
        """
        Function for mapping a scattered light image to a height
        profile (i.e. height of the scattering surface as function
        of radius in the disk midplane) that is either parameterized
        with a power-law function or read from an input file (for
        example returned by a radiative transfer code).

        Parameters
        ----------
        power_law : tuple(float, float, float)
            The height of the scattering surface is a power-law
            function, :math:`h(r) = a + b*(r/1\\,\\mathrm{au})^c`,
            with :math:`a`, :math:`b`, :math:`r`, and :math:`h(r)` in
            au. The argument of ``power_law`` should be provided as
            ``(a, b, c)``. Set all values to zero for the mapping
            a geometrically flat disk, in which case only the
            inclination is used for the deprojection.
        radius : tuple(float, float, int)
            Radius points that are sampled, provided as (r_in, r_out,
            n_r), with ``r_in`` and ``r_out`` in au. The outer radius
            should be set large enough such that a radius is sampled for
            each pixel in the field of view. To check if any NaNs are
            present, have a look at the `_radius.fits` output.
        surface : str
            Parameterization type for the disk surface ('power-law' or
            'function' or 'file').
        height_func : callable, None
            Function that returns the height of the scattering surface
            as a function of radius. The radii and returned height must
            be in au. Only used if surface='function'.
        filename : str, None
            Filename which contains the radius in au (first column) and
            the height of the disk surface in au (second column).

        Returns
        -------
        NoneType
            None
        """

        if surface == "power-law":

            # Power-law disk height

            @typechecked
            def power_law_height(
                x_power: np.ndarray, a_power: float, b_power: float, c_power: float
            ) -> np.ndarray:

                return a_power + b_power * x_power ** c_power

            # midplane radius (au)
            disk_radius = np.linspace(radius[0], radius[1], radius[2])

            # disk height (au)
            disk_height = power_law_height(
                disk_radius, power_law[0], power_law[1], power_law[2]
            )

            # opening angle (rad)
            disk_opening = np.arctan2(disk_height, disk_radius)
            
        elif surface == "function":
            
            if height_func is None:
                raise ValueError("If using surface=='function', you must specify height_func")
            
            # midplane radius (au)
            disk_radius = np.linspace(radius[0], radius[1], radius[2])

            # disk height (au)
            disk_height = height_func(disk_radius)

            # opening angle (rad)
            disk_opening = np.arctan2(disk_height, disk_radius)

        elif surface == "file":

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

        disk_phi = np.linspace(0.0, 359.0, 360)  # (deg)
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

                y_tmp = disk_height[i] * math.sin(self.incl) - r_item * np.cos(
                    p_item
                ) * math.cos(self.incl)

                x_rot = x_tmp * math.cos(math.pi - self.pos_ang) - y_tmp * math.sin(
                    math.pi - self.pos_ang
                )

                y_rot = x_tmp * math.sin(math.pi - self.pos_ang) + y_tmp * math.cos(
                    math.pi - self.pos_ang
                )

                x_im.append(x_rot)
                y_im.append(y_rot)

                r_im.append(math.sqrt(r_item ** 2 + disk_height[i] ** 2))
                p_im.append(p_item)
                o_im.append(disk_opening[i])

                ang_tmp = math.pi / 2.0 + disk_opening[i]

                par1 = (
                    math.sin(ang_tmp) * math.cos(math.pi + p_item) * math.sin(self.incl)
                )

                par2 = math.cos(ang_tmp) * math.cos(self.incl)

                s_im.append(math.pi - math.acos(par1 + par2))

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

        grid_xy = np.zeros((self.grid ** 2, 2))

        count = 0

        for i in range(-(self.grid - 1) // 2, (self.grid - 1) // 2 + 1):
            for j in range(-(self.grid - 1) // 2, (self.grid - 1) // 2 + 1):
                grid_xy[count, 0] = float(i)
                grid_xy[count, 1] = float(j)

                count += 1

        image_xy = np.zeros((len(x_im), 2))

        for i, _ in enumerate(x_im):
            image_xy[i, 0] = x_im[i]
            image_xy[i, 1] = y_im[i]

        # Interpolate image plane

        if self.npix % 2 == 0:
            x_grid = np.linspace(-self.npix / 2 + 0.5, self.npix / 2 - 0.5, self.npix)
            y_grid = np.linspace(-self.npix / 2 + 0.5, self.npix / 2 - 0.5, self.npix)

        elif self.npix % 2 == 1:
            x_grid = np.linspace(-(self.npix - 1) / 2, (self.npix - 1) / 2, self.npix)
            y_grid = np.linspace(-(self.npix - 1) / 2, (self.npix - 1) / 2, self.npix)

        x_grid *= self.pixscale * self.distance  # (au)
        y_grid *= self.pixscale * self.distance  # (au)

        grid = np.zeros((self.npix ** 2, 2))

        count = 0

        for i in range(self.npix):
            for j in range(self.npix):
                grid[count, 0] = x_grid[i]
                grid[count, 1] = x_grid[j]

                count += 1

        fit_radius = griddata(image_xy, r_im, grid, method="linear")
        fit_azimuth = griddata(image_xy, p_im, grid, method="linear")
        fit_opening = griddata(image_xy, o_im, grid, method="linear")
        fit_scatter = griddata(image_xy, s_im, grid, method="linear")

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
    def sort_disk(self) -> Tuple[List[np.float64], np.ndarray]:
        """
        Function for creating a list with pixel values and creating a
        2D array with the x and y pixel coordinates.

        Returns
        -------
        NoneType
            None
        """

        # Create lists with x and y coordinates and pixel values

        x_disk = []
        y_disk = []
        im_disk = []

        for i in range(self.npix):
            for j in range(self.npix):
                x_tmp = self.radius[i, j] * np.cos(self.azimuth[i, j])
                y_tmp = self.radius[i, j] * np.sin(self.azimuth[i, j])

                x_disk.append(
                    x_tmp * math.cos(math.pi / 2.0 - self.pos_ang)
                    - y_tmp * math.sin(math.pi / 2.0 - self.pos_ang)
                )

                y_disk.append(
                    x_tmp * math.sin(math.pi / 2.0 - self.pos_ang)
                    + y_tmp * math.cos(math.pi / 2.0 - self.pos_ang)
                )

                im_disk.append(self.image[i, j])

        # Sort disk plane points along x-axis

        x_index = np.argsort(x_disk)

        y_sort = np.zeros(len(y_disk))
        im_sort = np.zeros(len(y_disk))

        for i in range(len(y_disk)):
            y_sort[i] = y_disk[x_index[i]]
            im_sort[i] = im_disk[x_index[i]]

        # count = 0
        #
        # grid_xy = np.zeros((self.grid**2, 2))
        #
        # for i in range(-(self.grid-1)//2, (self.grid-1)//2+1):
        #     for j in range(-(self.grid-1)//2, (self.grid-1)//2+1):
        #         grid_xy[count, 0] = float(i)
        #         grid_xy[count, 1] = float(j)
        #
        #         count += 1

        image_xy = np.zeros((len(y_disk), 2))

        for i, _ in enumerate(y_disk):
            image_xy[i, 0] = x_disk[i]
            image_xy[i, 1] = y_disk[i]

        return im_disk, image_xy

    @typechecked
    def deproject_disk(self) -> None:
        """
        Function for deprojecting a disk surface based on the mapping
        of :meth:`~diskmap.diskmap.DiskMap.map_disk`.

        Returns
        -------
        NoneType
            None
        """

        if self.radius is None or self.azimuth is None:
            raise ValueError(
                "The disk has not been mapped yet with the 'map_disk' method."
            )

        im_disk, image_xy = self.sort_disk()

        # Interpolate disk plane

        if self.npix % 2 == 0:
            x_grid = np.linspace(-self.npix / 2 + 0.5, self.npix / 2 - 0.5, self.npix)
            y_grid = np.linspace(-self.npix / 2 + 0.5, self.npix / 2 - 0.5, self.npix)

        elif self.npix % 2 == 1:
            x_grid = np.linspace(-(self.npix - 1) / 2, (self.npix - 1) / 2, self.npix)
            y_grid = np.linspace(-(self.npix - 1) / 2, (self.npix - 1) / 2, self.npix)

        x_grid *= self.pixscale * self.distance  # (au)
        y_grid *= self.pixscale * self.distance  # (au)

        grid = np.zeros((self.npix ** 2, 2))

        count = 0

        for i in range(self.npix):
            for j in range(self.npix):
                grid[count, 0] = x_grid[i]
                grid[count, 1] = x_grid[j]

                count += 1

        try:
            fit_im = griddata(image_xy, im_disk, grid, method="linear")

        except ValueError:
            raise RuntimeError(
                "The radius sampling should cover the "
                "complete field of view of the image. Try "
                "increasing the outer 'radius' value in "
                "'map_disk' and have a look at the "
                "'_radius.fits' output to check for NaNs."
            )

        self.im_deproj = np.zeros((self.npix, self.npix))

        count = 0

        for i in range(self.npix):
            for j in range(self.npix):
                self.im_deproj[i, j] = fit_im[count]

                count += 1

    @typechecked
    def r2_scaling(self,
                   r_max: float,
                   mask_planet: Optional[Tuple[int, int, float, float]] = None) -> None:
        """
        Function for correcting a scattered light image for the r^2
        decrease of the stellar irradiation of the disk surface.

        Parameters
        ----------
        r_max : float
            Maximum disk radius (au) for the r^2-scaling. Beyond this
            distance, a constant r^2-scaling is applied of value
            ``r_max``.
        mask_planet : tuple(int, int, float, float), None
            Mask for a planet such that it will not be scaled. The
            tuple should have the following format:
            ``(x_planet, y_planet, r_mask, scaling)``. Here,
            ``x_planet`` and ``y_planet`` are the central pixel of
            the planet position, ``r_mask`` is the size (in pixels)
            of the planet signal, and ``scaling`` is the scaling
            factor of the planet flux. This parameter is a bit
            experimental and typically not used by setting the
            argument to ``None``.

        Returns
        -------
        NoneType
            None
        """

        if self.radius is None:
            raise ValueError("Please run 'map_disk' before using 'r2_scaling'.")

        if mask_planet is not None:
            x_grid = np.linspace(-mask_planet[0], self.npix-mask_planet[0], self.npix)
            y_grid = np.linspace(-mask_planet[1], self.npix-mask_planet[1], self.npix)

            xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)
            rr_grid = np.sqrt(xx_grid**2 + yy_grid**2)

        # self.low_snr = np.zeros((self.npix, self.npix))
        #
        # for i in range(3, self.npix-3):
        #     for j in range(3, self.npix-3):
        #         aperture = np.nansum(self.image[i-2:i+3, j-2:j+3])
        #
        #         if aperture > 90. and j < 157:
        #             self.low_snr[i, j] = np.nan
        #         else:
        #             self.low_snr[i, j] = self.image[i, j].copy()
        #
        # std_low = np.nanstd(self.low_snr)
        #
        # for i in range(3, self.npix-3):
        #     for j in range(3, self.npix-3):
        #         if not np.isnan(self.low_snr[i, j]) and np.abs(self.image[i, j]) > 3.*std_low:
        #             select = self.low_snr[i-3:i+4, j-3:j+4].copy()
        #             select[3, 3] = np.nan
        #
        #             self.image[i, j] = np.nanmedian(select)

        # self.low_snr = sigma_clip(self.low_snr, 5., maxiters=10)

        # for k in range(2):
        #     for i in range(3, self.npix-3):
        #         for j in range(3, self.npix-3):
        #             if not np.isnan(self.low_snr[i, j]):
        #                 select = self.low_snr[i-3:i+4, j-3:j+4].copy()
        #                 select[3, 3] = np.nan
        #
        #                 if self.image[i, j] > 0.5*np.nanstd(select):
        #                     self.image[i, j] = np.nanmedian(select)

        # for i in range(3, self.npix-3):
        #     for j in range(3, self.npix-3):
        #         if not np.isnan(self.high_snr[i, j]):
        #             self.high_snr[i, j] = self.low_snr[i, j]
        #         # elif type(self.low_snr[i, j]) is ma.core.MaskedConstant:
        #         #     self.high_snr[i, j] = np.nansum(self.low_snr[i-2:i+3, j-2:j+3])
        #         # else:
        #         #     self.high_snr[i, j] = self.image[i, j]
        #         else:
        #             self.high_snr[i, j] = self.low_snr[i, j].copy()

        # for k in range(2):
        # for i in range(3, self.npix-3):
        #     for j in range(3, self.npix-3):
        #         if np.isnan(self.filt[i, j]):
        #             select = self.image[i-2:i+3, j-2:j+3].copy()
        #             select[2, 2] = np.nan
        #
        #             if np.abs(self.image[i, j]) > 1.5*np.nanstd(select):
        #                 self.image[i, j] = np.nanmedian(select)

        # self.filt[i, j] = median_filter(self.image, 2)
        # self.image = self.high_snr.copy()

        self.im_scaled = np.zeros((self.npix, self.npix))

        for i in range(self.npix):
            for j in range(self.npix):
                if self.radius[i, j] < r_max:
                    if mask_planet is None or rr_grid[i, j] > mask_planet[2]:
                        self.im_scaled[i, j] = self.radius[i, j]**2 * self.image[i, j]

                    else:
                        self.im_scaled[i, j] = mask_planet[3] * self.image[i, j]

                else:
                    if mask_planet is None or rr_grid[i, j] > mask_planet[2]:
                        self.im_scaled[i, j] = r_max**2 * self.image[i, j]

                    else:
                        self.im_scaled[i, j] = mask_planet[3] * self.image[i, j]

    @typechecked
    def total_intensity(self, pol_max: float = 1.0) -> None:
        """
        Function for estimating the (stellar irradiation corrected)
        total intensity image when ``fitsfile`` contains a polarized
        light image and ``image_type='polarized'``. A bell-shaped
        degree of polarized is assumed and effects of multiple
        scattering are ignored.

        Parameters
        ----------
        pol_max : float
            The peak of the bell-shaped degree of polarization, which
            effectively normalizes the estimated total intensity image.

        Returns
        -------
        NoneType
            None
        """

        if self.image_type != "polarized":
            raise ValueError(
                "The 'total_intensity' method should only be "
                "used if the input image is a polarized light "
                "image (i.e. image_type='polarized')."
            )

        if self.scatter is None or self.im_scaled is None:
            raise ValueError("Please run 'map_disk' before using 'total_intensity'.")

        alpha = np.cos(self.scatter)
        deg_pol = -pol_max * (alpha ** 2 - 1.0) / (alpha ** 2 + 1.0)

        self.stokes_i = self.im_scaled / deg_pol

    @typechecked
    def phase_function(self, radius: Tuple[float, float], n_phase: int):
        """
        Function for extracting the phase function. If
        ``image_type='polarized'``, the polarized phase function is
        extracted and the total intensity phase function is estimated
        by assuming a bell-shaped degree of polarization. If
        ``image_type='polarized'``, the total intensity phase function
        is extracted. The extracting is done on the r$^2$-scaled pixel
        values such that the phase function is not biased by
        irradiation effects. The phase functions are have been
        normalized by their maximum value.

        Parameters
        ----------
        radius : tuple(float, float)
            Inner and outer radius (au) between which pixels are
            selected for estimating the phase function.
        n_phase : int
            Number of sampling points for the phase function between
            0 and 180 deg.

        Returns
        -------
        NoneType
            None
        """

        # Phase function is normalizedf so pol_max has not effect
        pol_max = 1.0

        if self.radius is None or self.scatter is None or self.im_scaled is None:
            raise ValueError(
                "Please run 'map_disk' and 'r2_scaling' "
                "before using 'phase_function'."
            )

        scat_select = []
        im_select = []

        for i in range(self.npix):
            for j in range(self.npix):
                if self.radius[i, j] > radius[0] and self.radius[i, j] < radius[1]:
                    scat_select.append(math.degrees(self.scatter[i, j]))

                    # use im_scaled to correct for differences across the
                    # selected radius
                    im_select.append(self.im_scaled[i, j])

        phase_bins = np.linspace(0.0, 180.0, num=n_phase)
        bin_index = np.digitize(scat_select, phase_bins)

        im_bins = []
        scat_bins = []

        for i in range(n_phase):
            im_bins.append([])
            scat_bins.append([])

        for i, _ in enumerate(im_select):
            im_bins[bin_index[i] - 1].append(im_select[i])
            scat_bins[bin_index[i] - 1].append(scat_select[i])

        angle = []
        pol_flux = []
        pol_error = []

        for i in range(n_phase):
            if len(im_bins[i]) > 0:
                angle.append(np.nanmean(scat_bins[i]))
                pol_flux.append(np.nanmean(im_bins[i]))
                pol_error.append(np.nanstd(im_bins[i]) / math.sqrt(len(im_bins[i])))

        if self.image_type == "polarized":
            # Degree of polarization

            alpha = np.cos(np.array(angle) * np.pi / 180.0)
            deg_pol = -pol_max * (alpha * alpha - 1.0) / (alpha * alpha + 1.0)

            tot_flux = pol_flux / deg_pol
            tot_error = pol_error / deg_pol

            # Normalization

            flux_norm = max(pol_flux)
            pol_flux = np.array(pol_flux) / flux_norm
            pol_error = np.array(pol_error) / flux_norm

            flux_norm = max(tot_flux)
            tot_flux = np.array(tot_flux) / flux_norm
            tot_error = np.array(tot_error) / flux_norm

            self.phase = np.column_stack(
                [angle, pol_flux, pol_error, tot_flux, tot_error]
            )

        else:
            self.phase = np.column_stack([angle, pol_flux, pol_error])

    @typechecked
    def write_output(self, filename: str) -> None:
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

        if self.radius is not None:
            fits.writeto(f"{filename}_radius.fits", self.radius, overwrite=True)

        if self.scatter is not None:
            fits.writeto(
                f"{filename}_scat_angle.fits", np.degrees(self.scatter), overwrite=True
            )

        if self.im_deproj is not None:
            fits.writeto(f"{filename}_deprojected.fits", self.im_deproj, overwrite=True)

        if self.im_scaled is not None:
            fits.writeto(f"{filename}_r2_scaled.fits", self.im_scaled, overwrite=True)

        if self.stokes_i is not None:
            fits.writeto(
                f"{filename}_total_intensity.fits", self.stokes_i, overwrite=True
            )

        if self.phase is not None:
            if self.image_type == "polarized":
                header = (
                    "Scattering angle (deg) - Normalized polarized "
                    + "flux - Error - Normalized total flux - Error"
                )

            else:
                header = "Scattering angle (deg) - Normalized total flux - Error"

            np.savetxt(f"{filename}_phase_function.dat", self.phase, header=header)
