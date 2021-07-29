.. _installation:

Installation
============

*diskmap* is compatible with Python 3.6/3.7/3.8/3.9 and available through `PyPI <https://pypi.org/project/diskmap/>`_ and `Github <https://github.com/tomasstolker/diskmap>`_.

Installation from PyPI
----------------------

*diskmap* can be installed with the `pip package manager <https://packaging.python.org/tutorials/installing-packages/>`_:

.. code-block:: console

    $ pip install diskmap

Or to update the package to the most recent version:

.. code-block:: console

   $ pip install --upgrade diskmap

Installation from Github
------------------------

Installation from Github is done by cloning the repository:

.. code-block:: console

    $ git clone git@github.com:tomasstolker/diskmap.git

And running the setup script to install *diskmap* and its dependencies:

.. code-block:: console

    $ python setup.py install

Once a local copy of the repository exists, new commits can be pulled from Github with:

.. code-block:: console

    $ git pull origin master

Do you want to makes changes to the code? Then please fork the `diskmap` repository on the Github page and clone your own fork instead of the main repository. Contributions and pull requests are welcome (see :ref:`about` section).

Testing `diskmap`
-----------------

The installation can be tested by starting Python in interactive mode, importing *diskmap*, and printing the installed version number:

.. code-block:: python

    >>> import diskmap
    >>> diskmap.__version__
