.. _installation:

Installation
============

*diskmap* is compatible with Python 3.6/3.7/3.8 and is available in the |pypi| and on |github|.

Installation from PyPI
----------------------

*diskmap* can be installed with the |pip|:

.. code-block:: console

    $ pip install diskmap

And to update to the most recent version:

.. code-block:: console

   $ pip install --upgrade diskmap

Installation from Github
------------------------

Installation from Github is done by cloning the repository:

.. code-block:: console

    $ git clone git@github.com:tomasstolker/diskmap.git

And running the setup script to install the package and its dependencies:

.. code-block:: console

    $ python setup.py install

Once a local copy of the repository exists, new commits can be pulled from Github with:

.. code-block:: console

    $ git pull origin master

Do you want to makes changes to the code? Then please fork the `diskmap` repository on the Github page and clone your own fork instead of the main repository. Contributions and pull requests are very welcome (see :ref:`about` section).

Testing `diskmap`
-----------------

The installation can be tested by starting Python in interactive mode and printing the `diskmap` version:

.. code-block:: python

    >>> import diskmap
    >>> diskmap.__version__

.. |pypi| raw:: html

   <a href="https://pypi.org/project/diskmap/" target="_blank">PyPI repository</a>

.. |github| raw:: html

   <a href="https://github.com/tomasstolker/diskmap" target="_blank">Github</a>

.. |pip| raw:: html

   <a href="https://packaging.python.org/tutorials/installing-packages/" target="_blank">pip package manager</a>
