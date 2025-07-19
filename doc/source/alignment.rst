Structure alignment
___________________

MATILDA provides a variety of tools for structure alignment and computing
root-mean-square deviations (RMSD).

Superimpose
~~~~~~~~~~~

Superposition (alignment) of different structures of one molecule using a quaternion fit.

.. code:: text

    superimpose.py

    superimpose.py <ref_struc> <struc1> <struc2> ... [-t<type=xyz>] [-mwp<mass_weight_power=1>]

Example

.. code:: text

    superimpose.py ref.xyz struc1.xyz struc2.xyz -txyz

Distance/RMSD
~~~~~~~~~~~~~

Computation of the distances (i.e. norm of the difference vector, RMSD) between a set of structures of one molecule.

.. code:: text

    distances.py

    distances.py <struc1> <struc2> ... [-t<type>] [-mwp<mass_weight_power>] [-fit<fitting>] [-dig<digits>]

Example

.. code:: text

    distances.py struc1.xyz struc2.xyz -txyz

Insert molecule
~~~~~~~~~~~~~~~

Insert molecules into a larger template file and optionally combine the molecular vibrations.

.. code:: text

    combine_vib.py

*Requirements*: Structure files for the template and the molecules to insert.
Normal modes (if required) are given as molden format files.

*Execution*: Follow the explanations on the screen.
The non-trivial task is to identify the corresponding atoms between the structures.
For this it is suggested to open both molecules at the same time in a molecular structure viewing program to compare the indices.
Note that the program asks for the indices of the atoms in the (larger) template file corresponding to indices in the molecule inserted (putting the indices in the opposite order will give a wrong result).

*Output*:
* `combined_struc.xyz`: Structure file with inserted structures.
* `combined_vib.mld`: Molden file with combined normal modes.
