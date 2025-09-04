Internal Coordinates
____________________

MATILDA provides a variety of tools for analysing internal coordinates.

Basic analysis
~~~~~~~~~~~~~~

Analyse internal coordinates of a molecule.
This also works for a file containing multiple geometries of one molecule (e.g. coming from dynamics).

::

    int_coor_multi.py

    usage: int_coor_multi.py [-h] [-f FILENAME [FILENAME ...]] [--filetype FILETYPE] [-d DIST DIST]
                         [-b BEND BEND BEND] [-t TORS TORS TORS TORS] [-dg DIGITS]

    Outputs the value of the speciefied internal coordinate

    options:
    -h, --help            show this help message and exit
    -f FILENAME [FILENAME ...], --filename FILENAME [FILENAME ...]
    --filetype FILETYPE
    -d DIST DIST, --dist DIST DIST
                        The indices of atoms for distance
    -b BEND BEND BEND, --bend BEND BEND BEND
                        The indices of atoms for bend angle
    -t TORS TORS TORS TORS, --tors TORS TORS TORS TORS
                        The indices of atoms for dihedral angle
    -dg DIGITS, --digits DIGITS
                        Number of decimal points

	
Example

::

    int_coor_multi.py -d 1 2 -b 1 2 4 -t 5 1 2 6 -f dyn.xyz

Lanthanide complexes
~~~~~~~~~~~~~~~~~~~~

The ``lanth_geo.py`` script provides functionality specifically for lanthanide complexes.
Aside from bond lengths it computes the SAP/tSAP twisting angle for the main eight coordination sites.

To use this functionality, one needs to specify four `BOTTOM_INDS` defining the bottom square and
four `TOP_INDS` defining the top square.
Alternatively, one can specify 3 indices each to analyse a trigonal (anti)prism geometry.

.. code:: text

    lanth_geo.py

    usage: lanth_geo.py [-h] [-f FILENAME] [--filetype FILETYPE] [-l LANTH_IND] [-b BOTTOM_INDS] [-t TOP_INDS]

    Determine SAP or tSAP geomtry for lanthanide complex

    options:
    -h, --help            show this help message and exit
    -f FILENAME, --filename FILENAME
    --filetype FILETYPE
    -l LANTH_IND, --lanth_ind LANTH_IND
                        index defining the lanthanide
    -b BOTTOM_INDS, --bottom_inds BOTTOM_INDS
                        indices defining the bottom square
    -t TOP_INDS, --top_inds TOP_INDS
                        indices defining the top square

Example

.. code:: text

    lanth_geo.py -f nucfid.xyz -b 2 23 16 9 -t 33 83 77 71
