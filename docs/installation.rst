.. _installation:

Installation
============

``diskmap`` is available on `PyPI <https://pypi.org/project/diskmap/>`_ and `Github <https://github.com/tomasstolker/diskmap>`_.

Installation from PyPI
----------------------

``diskmap`` can be installed from `PyPI <https://pypi.org/project/diskmap/>`_  with the `pip package manager <https://packaging.python.org/tutorials/installing-packages/>`_:

.. code-block:: console

    $ pip install diskmap

Or to update the package to the most recent version:

.. code-block:: console

   $ pip install --upgrade diskmap

Installation from Github
------------------------

Using pip
^^^^^^^^^

Installation from Github is also possible with ``pip``:

.. code-block:: console

   $ pip install git+https://github.com/tomasstolker/diskmap.git

Cloning the repository
^^^^^^^^^^^^^^^^^^^^^^

Alternatively, the Github repository can be cloned, which is in particular useful if you want to look into and/or make changes to the code

.. code-block:: console

    $ git clone https://github.com/tomasstolker/diskmap.git

The package is installed by running ``pip`` in the local cloned repository:

.. code-block:: console

    $ pip install -e .

Once a local copy of the repository exists, new commits can be pulled from Github with:

.. code-block:: console

    $ git pull origin main

Do you want to makes changes to the code? Then please fork the `diskmap` repository on the Github page and clone your own fork instead of the main repository. Contributions and pull requests are welcome (see :ref:`about` section).

Testing `diskmap`
-----------------

The installation can be tested by starting Python in interactive mode, importing ``diskmap``, and printing the installed version number:

.. code-block:: python

    >>> import diskmap
    >>> diskmap.__version__
