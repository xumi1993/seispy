Installation
=================

Dependencies
---------------

The current version has been integration tested using **Python 3**. 

Following Python modules are required for the Seispy, **but do not have to be installed manually.** Once the Seispy has been installed, dependencies are automatically installed to your environment.

- `Numpy <https://numpy.org/>`_
- `Scipy <https://www.scipy.org/scipylib/index.html>`_
- `Obspy <https://docs.obspy.org/>`_
- `Matplotlib <https://matplotlib.org/>`_
- `PySide6 <https://doc.qt.io/qtforpython/index.html>`_


Install and update via `Anaconda <https://www.anaconda.com/>`_ 
-------------------------------------------------------------------

We recommend to install the `Anaconda <https://www.anaconda.com/>`_ as the Python environment. 

Installing ``seispy`` from the ``conda-forge`` channel can be achieved by adding ``conda-forge`` to your channels with:

.. code-block:: shell

    conda config --add channels conda-forge
    conda config --set channel_priority strict

.. note::

    For Chinese users, we recommend to change the source of the ``conda-forge`` to `清华大学 Anaconda 镜像 <https://mirrors.tuna.tsinghua.edu.cn/help/anaconda/>`_

Once the ``conda-forge`` channel has been enabled, ``seispy`` can be installed with:

.. code-block:: shell

    conda create -n seispy matplotlib-base obspy pandas pyside6 scikits-bootstrap pyproj numba
    conda activate seispy
    pip install python-seispy

``Seispy`` can be updated with:

.. code-block:: shell

    pip install python-seispy -U


Install and update from source codes
--------------------------------------

Download source codes
^^^^^^^^^^^^^^^^^^^^^^^^

Clone the source code from `Github <https://github.com/xumi1993/seispy.git>`_ to any directory.

.. code-block:: shell

    git clone https://github.com/xumi1993/seispy.git

To access developing version, users can clone the source codes with 

.. code-block:: shell

    git clone --branch=dev https://github.com/xumi1993/seispy.git

For **Chinese users**, try to clone the source code from `Gitlab repository <https://gitlab.com/xumi1993/seispy.git>`_

.. code-block:: shell

    git clone https://gitlab.com/xumi1993/seispy.git

or

.. code-block:: shell

    git clone --branch=dev https://gitlab.com/xumi1993/seispy.git

Install Seispy to the Python environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Change path to where the source code was cloned into, and install the module via `Python pip <https://pip.pypa.io/>`_ command

.. code-block:: shell

    cd path/to/seispy
    pip install .

Update Seispy
^^^^^^^^^^^^^^^^

To update the Seispy, please change to directory of the source code, and execute following commands.

.. code-block:: shell

    cd path/to/seispy
    git pull
    pip install .


FAQ
-------
- When old users update seispy to v1.3.0, An error probably raised because of incompatible Qt library
    .. code-block:: shell

        qt.core.plugin.loader: In /Users/zhangxiaoqing/miniconda3/plugins/platforms/libqwebgl.dylib:
        Plugin uses incompatible Qt library (5.15.0) [release]
        qt.core.plugin.loader: In /Users/zhangxiaoqing/miniconda3/plugins/platforms/libqoffscreen.dylib:
        Plugin uses incompatible Qt library (5.15.0) [release]
        qt.core.plugin.loader: In /Users/zhangxiaoqing/miniconda3/plugins/platforms/libqminimal.dylib:
        Plugin uses incompatible Qt library (5.15.0) [release]
        qt.core.plugin.loader: In /Users/zhangxiaoqing/miniconda3/plugins/platforms/libqcocoa.dylib:
        Plugin uses incompatible Qt library (5.15.0) [release]
        qt.qpa.plugin: Could not find the Qt platform plugin "cocoa" in ""
        This application failed to start because no Qt platform plugin could be initialized. Reinstalling the application may fix this problem.

    Please uninstall PyQt5 with ``pip`` or ``conda``

    .. code-block:: shell

        pip uninstall pyqt5

    or 

    .. code-block:: shell

        conda remove pyqt