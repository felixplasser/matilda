Vibronic Coupling
_________________

Tools for analysing vibronic coupling effects.

Huang-Rhys factors
~~~~~~~~~~~~~~~~~~

The ``huang_rhys.py`` script is used for computing Huang-Rhys (HR) factors, which quantify the electron-phonon coupling between 
the ground and excited state geometries of a molecule along vibrational normal modes.
It also calculates the reorganization energy, representing the nuclear reconfiguration cost during electronic transitions.

.. code:: text

    huang_rhys.py <gs_struc> <es_struc> <vib_mld>

* ``gs_struc``: Ground state geometry file in .xyz format
* ``es_struc``: Excited-state geometry file in .xyz format
* ``vib_mld``: Vibrational normal modes and frequencies in molden format.
  These should normally be computed at the ground state geometry.

Example 

.. code:: text

    huang_rhys.py GS.xyz ES.xyz freq.molden

Normal mode displacement
~~~~~~~~~~~~~~~~~~~~~~~~

Displace a molecular structure along a normal mode.

.. code:: text

    add_normal_mode.py <in_struc> <out_struc> <vib_file> <nm_ind> <disp> [<type>]
