import math
import numpy as np

from astropy.io import fits
from scipy.interpolate import griddata, interp1d


class DiskMap:
    def __init__(self,
                 fitsfile,
                 pixscale,
                 inclination,
                 pos_angle,
                 distance):

        self.image = fits.getdata(fitsfile)
        self.image = np.nan_to_num(self.image)

        if self.image.ndim != 2:
            raise ValueError('DiskMap requires a 2D image.')

        if self.image.shape[0] != self.image.shape[1]:
            raise ValueError('The dimensions of the image should have the same size.')

        self.pixscale = pixscale  # [arcsec]
        self.incl = inclination  # [deg]
        self.pos_ang = pos_angle  # [deg]
        self.distance = distance  # [pc]

        self.npix = self.image.shape[0]

        self.incl = math.radians(self.incl)  # [rad]
        self.pos_ang = math.radians(self.pos_ang-90.)  # [rad]

        self.grid = 500

        self.radius = None
        self.opening = None
        self.scatter = None
        self.im_scaled = None
        self.stokes_i = None
        self.phase = None

    def map_disk(self,
                 surface,
                 power_law,
                 radius,
                 filename=None):

        # Create geometric disk model

        if surface == 'power-law':

            # Power-law disk height

            def power_law_height(x_power, a_power, b_power, c_power):
                return a_power + b_power*x_power**c_power

            radius = np.linspace(radius[0], radius[1], radius[2])  # [au]
            height = power_law_height(radius, power_law[0], power_law[1], power_law[2])  # [au]
            opening = np.arctan2(height, radius)

        elif surface == 'file':

            # Read disk height from ASCII file

            data = np.loadtxt(filename)
            radius = data[:, 0]  # [au]
            height = data[:, 1]  # [au]

            radius = np.linspace(radius[0], radius[-1], 100)  # [au]

            height_interp = interp1d(radius, height)
            height = height_interp(radius)  # [au]

            opening = np.arctan2(height, radius)  # [au]

        # Project disk height to image plane

        phi = np.linspace(0., 359., 360)  # [deg]
        phi = np.radians(phi)  # [rad]

        x_im = []
        y_im = []
        r_im = []
        o_im = []
        s_im = []
        p_im = []

        for i, r_item in enumerate(radius):
            for j, p_item in enumerate(phi):

                x_temp = r_item * np.sin(p_item)
                y_temp = height[i]*math.sin(self.incl) - r_item*np.cos(p_item)*math.cos(self.incl)

                x_rot = x_temp*math.cos(self.pos_ang) - y_temp*math.sin(self.pos_ang)
                y_rot = x_temp*math.sin(self.pos_ang) + y_temp*math.cos(self.pos_ang)

                x_im.append(x_rot)
                y_im.append(y_rot)

                r_im.append(math.sqrt(r_item**2+height[i]**2))
                p_im.append(p_item)
                o_im.append(opening[i])

                par1 = math.sin(math.pi/2.+opening[i])*math.cos(math.pi+p_item)*math.sin(self.incl)
                par2 = math.cos(math.pi/2.+opening[i])*math.cos(self.incl)

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
        for i in range(-self.grid // 2, self.grid // 2):
            for j in range(-self.grid // 2, self.grid // 2):
                grid_xy[count, 0] = float(i)
                grid_xy[count, 1] = float(j)
                count += 1

        image_xy = np.zeros((len(x_im), 2))
        for i, _ in enumerate(x_im):
            image_xy[i, 0] = x_im[i]
            image_xy[i, 1] = y_im[i]

        # Interpolate data

        fit_radius = griddata(image_xy, r_im, grid_xy)
        fit_phi = griddata(image_xy, p_im, grid_xy)
        fit_opening = griddata(image_xy, o_im, grid_xy)
        fit_scatter = griddata(image_xy, s_im, grid_xy)

        grid_size = int(math.sqrt(np.size(fit_radius)))
        radius = np.zeros((grid_size, grid_size))
        phi = np.zeros((grid_size, grid_size))
        opening = np.zeros((grid_size, grid_size))
        scatter = np.zeros((grid_size, grid_size))

        count = 0
        for i in range(grid_size):
            for j in range(grid_size):
                radius[j, i] = fit_radius[count]
                phi[j, i] = fit_phi[count]
                opening[j, i] = fit_opening[count]
                scatter[j, i] = fit_scatter[count]
                count += 1

        # Map image plane

        if self.npix % 2 == 0:
            x_grid = y_grid = np.linspace(-self.npix/2+0.5, self.npix/2-0.5, self.npix)
        elif self.npix % 2 == 1:
            x_grid = y_grid = np.linspace(-(self.npix-1)/2, (self.npix-1)/2, self.npix)

        x_grid *= self.pixscale*self.distance  # [au]
        y_grid *= self.pixscale*self.distance  # [au]

        grid = np.zeros((self.npix**2, 2))

        count = 0

        for i in range(self.npix):
            for j in range(self.npix):
                grid[count, 0] = x_grid[i]
                grid[count, 1] = x_grid[j]

                count += 1

        fit_radius = griddata(image_xy, r_im, grid)
        fit_opening = griddata(image_xy, o_im, grid)
        fit_scatter = griddata(image_xy, s_im, grid)

        self.radius = np.zeros((self.npix, self.npix))
        self.opening = np.zeros((self.npix, self.npix))
        self.scatter = np.zeros((self.npix, self.npix))

        count = 0

        for i in range(self.npix):
            for j in range(self.npix):
                self.radius[i, j] = fit_radius[count]
                self.opening[i, j] = fit_opening[count]
                self.scatter[i, j] = fit_scatter[count]

                count += 1

    def r2_scaling(self,
                   r_max):

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

    def total_intensity(self,
                        r_max,
                        pol_max=0.5):

        # Total intensity

        if self.scatter is None or self.im_scaled is None:
            raise ValueError('Please run \'map_disk\' before using \'total_intensity\'.')

        alpha = np.cos(self.scatter)
        deg_pol = -pol_max*(alpha*alpha-1.)/(alpha*alpha+1.)
        self.stokes_i = self.im_scaled / deg_pol

    def phase_function(self,
                       radius,
                       n_phase,
                       pol_max=0.5):

        # Phase function

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

    def write_output(self,
                     filename):

        # Write FITS output

        if self.radius is not None:
            fits.writeto(f'{filename}_radius.fits', self.radius, overwrite=True)

        if self.scatter is not None:
            fits.writeto(f'{filename}_scatter.fits', np.degrees(self.scatter), overwrite=True)

        if self.im_scaled is not None:
            fits.writeto(f'{filename}_scaled.fits', self.im_scaled, overwrite=True)

        if self.stokes_i is not None:
            fits.writeto(f'{filename}_stokes_i.fits', self.stokes_i, overwrite=True)

        if self.phase is not None:
            header = 'Scattering angle [deg] - Polarized intensity [a.u.] - Uncertainty [a.u.] ' \
                     '- Total intensity [a.u.] - Uncertainty [a.u.]'

            np.savetxt(f'{filename}_phase_function.dat', self.phase, header=header)
