Calculate S wave Receiver Functions (SRF) for single station
---------------------------------------------------------------

Version after 1.2.13 has involved functions for SRF calculation. The workflow is similar to PRF calculation. 

Materials
==========

Users should prepare following materials.

- SAC files of teleseismic events including S wave. Refer to :doc:`../notes/rf_cfg`.
- A configure file with the same content with that of PRF. Refer to :doc:`./rf-sta`.

Run in command line
===================

Similar to the workflow for PRFs, we have provided a command ``srf`` for quickly calculating SRFs.

::

    usage: srf [-h] [-l] [-r N|E|NE] [-s] [-b [BAZ]] [-w] [-p] [-i] cfg_file

    Calculating RFs for single station

    positional arguments:
    cfg_file    Path to RF configure file

    optional arguments:
    -h, --help  show this help message and exit
    -l          use local catalog, defaults to false
    -r N|E|NE   Reverse components: N, E or NE
    -s          Switch the East and North components
    -b [BAZ]    Correct back-azimuth. If "baz" is specified, the corr_baz = raw_baz + baz. If there is no argument, the back-azimuth will be corrected with
                minimal energy of T component. The searching range is raw_baz +/- 90
    -w          Write project to localfile
    -p          Wether or not manually pick arrival time and waveforms arround S phase with a GUI.
    -i          Wether grid search incidence angle

.. note::

    Both ``TRZ`` and ``LQT`` are available for SRF calculation. the Z or L component of RF is deconvolved by R or Q component.

GUI for visual checking S wave
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Compared with ``prf``, ``srf`` has a new option of ``-p`` for visual checking seismic waveforms and picking S phase arrival with a GUI:

.. figure:: ../_static/srfui.png
    :align: center

    GUI for checking SRFs. Black lines denote manually picked arrival time. Red lines denote theoretical arrival time.

- This GUI is open after arrival time calculation and rotation. The seismic waveforms in 3 components are cut around theoretical S wave arrival time (the red line). A pre-picking is preformed with a STA/LTA method in T component (the black line).

- If the current event has a low quality, users can mark the event to "poor" via button ``mark to poor`` or hotkey ``d``.

- Directly click is used to manually pick S arrival time. The pre-picking line will move to where you click on.

Save checked events to a local file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Specified ``-w`` allow users to save events, parameters and picked S arrival time into a local ``pkl`` file. For re-calculation please load the file in Python script with:

.. code:: Python

    from seispy.rf import RF

    rf = RF.loadpjt('path/to/xxx.pkl')
    rf.detrend()
    rf.filter()
    rf.rotate()
    rf.trim()
    rf.deconv()
    rf.saverf()

Visual checking SRFs
========================

After calculation, the same step of visual checking as ``prf`` can be preformed via ``pickrf`` command. See :doc:`../usage/pickrf` in detail.


.. note::

    In this scheme, files of ``*_L.sac`` in ``rfpath`` will be searched in prior, when ``comp = LQT`` is set in configure file. If ``comp = TRZ``, files of ``*_Z.sac`` will be searched.
