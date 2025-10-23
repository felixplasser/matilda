Vibronic Coupling
_________________

Tools for analysing vibronic coupling effects.

Huang-Rhys factors
~~~~~~~~~~~~~~~~~~

The ``huang_rhys.py`` script is used for computing Huang-Rhys (HR) factors, which quantify the electron-phonon coupling between 
the ground and excited state geometries of a molecule along vibrational normal modes.
It also calculates the reorganization energy, representing the nuclear reconfiguration cost during electronic transitions.
The computed data are saved to an output file, which can be subsequently used as input for the ``fc_mini.py`` script.

To get detailed information about how to use this script, one can run the help command in the terminal:

.. code:: text

   huang_rhys.py -h

This command will display a comprehensive usage guide, including:

A brief description of the script’s purpose (computing Huang-Rhys parameters and reorganization energies).

Required Arguments: Ground state,  Excited-state geometry file, Selection of vibrational molden files for ground or excited state (-g or -e).

Other options: File format type (-t), Display options for top vibrational modes by Huang-Rhys factors or reorganization energies, plotting and saving output options.

.. code:: text

   usage: huang_rhys.py [-h] (-g GS_VIB | -e ES_VIB) [-t FILETYPE] [-s N] [-l N] [-p PNGFILE] gs_file es_file

   Compute Huang-Rhys parameters and reorganization energies.

   positional arguments:
    gs_file               Ground state geometry file (e.g., gs.xyz)
    es_file               Excited-state geometry file (e.g., es.xyz)


   options:
     -h, --help            show this help message and exit
     -g GS_VIB, --vib_gs GS_VIB
                        Ground-state vibrational Molden file (default)
     -e ES_VIB, --vib_es ES_VIB
                        Excited-state vibrational Molden file
     -t FILETYPE, --filetype FILETYPE
                        Filetype of structure files (default: xyz)
     -s N, --topS N        Show top N modes with highest Huang-Rhys factors
     -l N, --topLambda N   Show top N modes with highest reorganization energies
     -p PNGFILE, --plot PNGFILE
                        Plot Huang-Rhys spectrum and save to file

#NOTE: You are using excited-state vibrational modes (-e). If results seem unreasonable or unphysical, consider using ground-state modes (-g) instead, 
especially if excited-state frequencies are less reliable in your calculations.

Example 

.. code:: text

    huang_rhys.py GS.xyz ES.xyz -e freq.molden -s 10 -l 10 -p huang_rhys.png 

Normal mode displacement
~~~~~~~~~~~~~~~~~~~~~~~~

Displace a molecular structure along a normal mode.

.. code:: text

    add_normal_mode.py <in_struc> <out_struc> <vib_file> <nm_ind> <disp> [<type>]
