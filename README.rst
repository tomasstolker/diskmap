*diskmap*
=========

.. image:: https://badge.fury.io/py/diskmap.svg
    :target: https://badge.fury.io/py/diskmap

.. image:: https://img.shields.io/badge/Python-3.6%2C%203.7-yellow.svg?style=flat
    :target: https://pypi.python.org/pypi/diskmap

.. image:: https://img.shields.io/badge/MIT-blue.svg
    :target: https://github.com/tomasstolker/diskmap/blob/master/LICENSE

.. image:: http://img.shields.io/badge/arXiv-1609.09505-orange.svg?style=flat
    :target: https://arxiv.org/abs/1609.09505

Python tool for scattered light mapping of protoplanetary disks. The disk surface is parameterized with a power law profile or read from an input file. The projected radius and scattering angle is then calculated at each pixel. From this, a deprojected image, a stellar irradiation corrected image, and a scattering phase function can be extracted.

Installation
------------

Installation from PyPI:

.. code-block:: console

    $ pip install diskmap

Installation from Github:

.. code-block:: console

    $ git clone git@github.com:tomasstolker/diskmap.git
    $ python setup.py install

Example
-------

.. code-block:: python

   >>> import diskmap

   >>> mapping = diskmap.DiskMap(fitsfile='image.fits',
                                 pixscale=1e-2,
                                 inclination=40.,
                                 pos_angle=70.,
                                 distance=100.)

   >>> mapping.map_disk(power_law=(0., 0.1, 1.15),
                        radius=(1., 500., 100))

   >>> mapping.deproject_disk()

   >>> mapping.r2_scaling(r_max=200.)

   >>> mapping.total_intensity(pol_max=1.)

   >>> mapping.phase_function(radius=(50., 70.),
                              n_phase=30)

   >>> mapping.write_output(filename='diskmap')


Attribution
-----------

Please cite `Stolker et al. (2016) <https://ui.adsabs.harvard.edu/abs/2016A%26A...596A..70S>`_ whenever results from *diskmap* are used in a publication.

Contributing
------------

Contributions are welcome, please consider forking the repository and creating a pull request. Bug reports can be provided by creating an `issue <https://github.com/tomasstolker/diskmap/issues>`_ on the Github page.

License
-------

Copyright 2020 Tomas Stolker

*diskmap* is distributed under the MIT License. See the LICENSE file for the terms and conditions.
