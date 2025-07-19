Internal Coordinates
____________________

MATILDA provides a variety of tools for analysing internal coordinates.

Basic analysis
~~~~~~~~~~~~~~

Analyse internal coordinates of a molecule.
This also works for a file containing multiple geometries of one molecule (e.g. coming from dynamics).

::

    int_coor_multi.py

Example

::

    int_coor_multi.py dist 1 2 bend 1 2 4 tors 5 1 2 6 dyn.xyz

Lanthanide complexes
~~~~~~~~~~~~~~~~~~~~

The ``lanth_geo.py`` script provides functionality specifically for lanthanide complexes.
Aside from bond lengths it computes the SAP/tSAP twisting angle for the main eight coordination sites.

To use this functionality, one needs to specify four `BOTTOM_INDS` defining the bottom square and
four `TOP_INDS` defining the top square.
Alternatively, one can specify 3 indices each to analyse a trigonal (anti)prism geometry.

.. code:: text

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

    lanth_geo.py -f struc.xyz -l 1 -b '2 3 4 5' -t '6 7 8 9'
