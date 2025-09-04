Utilities
_________

Format conversion
~~~~~~~~~~~~~~~~~

File format conversion utility, based on and similar to the babel program but with additional support for Columbus/Newton-X ("col") and Tinker ("txyz2") formats.

::

    babel.py <in_format> <in_file> <out_format> <out_file>

Example (to convert from Turbomole format to Columbus)

.. code:: text

    babel.py tmol coord col geom

Renumber file
~~~~~~~~~~~~~

Renumbering of an xyz file.
Renumber the xyz-file `<in_file>`.
The file `<renumber_file>` contains the old numbers of the atoms in the new order
(if optionally some of the atoms are left out, they are added at the end).
The result is written into `<out_file>`.

::

    renumber_xyz.py <in_file> <renumber_file> <out_file>

Example

.. code:: text

    renumber_xyz.py struc.xyz renumber.txt struc_renumber.xyz
