Batch calculation and CCP stacking
=======================================

Calculate PRFs of a seismic array in batch
-------------------------------------------

This example provide a procedure for batch calculating PRFs of NCISP-III array (`Zheng et al., 2007 <https://doi.org/10.1016/j.pepi.2007.01.004>`_)

Download this example
^^^^^^^^^^^^^^^^^^^^^^^^^

Download link: `ex-batrf.tar.gz <https://osf.io/xghrk/download>`_

Download and unzip the file to any directory.

.. code-block:: shell

    wget https://osf.io/xghrk/download -O ex-batrf.tar.gz
    tar -xzf ex-batrf.tar.gz

The package include ``Data.ZX`` with 47 sub folders named as ``network.station``, which include ``sac`` files cut out from 300 s after origin time to 1000 s after origin time seismic events. The ``SACHeader` of station information (``netwk``, ``stnm``, ``stla`` and ``stlo``) had been written into ``sac`` files.

Write a script for batch calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Seispy has provided a command as ``setpar``, which will modify values in the ``.cfg`` file. Combining with a list of station information with 3 columns (``staname``, ``stla``, ``stlo``), We can use the ``setpar`` and ``prf`` command in shell script (``run.sh``) for RF calculation in batch. In this example, please execute ``run.sh`` in ``ex-batRF`` directory for the batch calculation.

.. code-block:: bash

    while read staname stla stlo
    do
        setpar rf.cfg path datapath ./Data.ZX/$staname
        setpar rf.cfg path rfpath ./RFresult/$staname
        prf rf.cfg
    done < sta.lst


- The ``sta.lst`` includes station information.
- The ``rf.cfg`` is the configure file. In addition to the ``datapath`` and ``rfpath``, other parameters can be set in advance.
- If set up a local catalog, the calculation time will be saved.

CCP stacking along a linear seismic array
------------------------------------------

Download this example
^^^^^^^^^^^^^^^^^^^^^^^

Download link: `ex-ccp.tar.gz <https://osf.io/hzq2x/download>`_

Download and unzip the file to any directory.

.. code-block:: shell

    wget https://osf.io/hzq2x/download -O ex-ccp.tar.gz
    tar -xzf ex-ccp.tar.gz

This package include:

- A Folder with PRFs (``RFresult``)
- A configure file with ``cfg`` format (``ccp.cfg``).

Time-to-depth conversion
^^^^^^^^^^^^^^^^^^^^^^^^^

specify following parameters in the ``ccp.cfg`` for time-to-depth conversion

- ``rfpath``: Path to PRFs
- ``depthdat``: Output data structure after time-to-depth conversion
- ``stalist``: The list of station information with 3 columns (``staname``, ``stla``, ``stlo``)

Then execute follow command to do the conversion

.. code-block:: shell

    rf2depth ccp.cfg

After the time-to-depth conversion, the ``RFdepth.mat`` will be generated following the path set as ``depthdat``.

CCP stacking
^^^^^^^^^^^^^^^^

specify following parameters in the ``ccp.cfg`` for CCP stacking. See :doc:`../usage/ccp` in detail.

- ``stackfile``: Output file for the CCP stacking
- ``stack_sta_list``: Stations used for the CCP stacking. In this example, we use the same stations as the ``stalist``.
- Section ``[bin]`` for parameters of stacking bins
- Section ``[line]`` for location of the profile.
- Section ``[stack]`` for Depth samples.

Then execute following command for CCP stacking

.. code-block:: shell

    ccp_profile ccp.cfg -t


An text file set by ``stackfile`` will be generated including 6 columns: ``bin_lat``, ``bin_lon``, ``bin_distance``, ``depth``, ``amplitude`` and ``stack_num``.

Plot the stacking image
^^^^^^^^^^^^^^^^^^^^^^^^^^

Seispy dose not provide any functions or script for plotting stacking image. Here we provide a script with GMT6 format for plotting the image. Just run in ``ex-ccp`` directory.

.. code-block:: shell

    sh ps_profile.sh

.. figure:: ../_static/profile_ZX.png
    :alt: profile
    :align: center

    CCP stacking along the seismic array

Reference
^^^^^^^^^^^

Zheng T, Chen L, Zhao L, et al. Crustal structure across the Yanshan belt at the northern margin of the North China Craton[J]. Physics of the Earth and Planetary Interiors, 2007, 161(1-2): 36-49.