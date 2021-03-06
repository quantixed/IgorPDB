# IgorPDB

PDB file load and display in IGOR Pro

Read a PDB file into Igor for display in gizmo. The reader is built on code written by Jan Ilavsky, John Weeks and tcaus.

- Keep alpha carbon backbone
- Chains as separate objects for easy manipulation
- Biological assembly using symmetry operators (BIOMT only)
- Downsample structures for simplified view

Limited testing using 3iyv.pdb on IGOR Pro 9 beta.

## Example

Using `IgorPDB()` a PDB file can be loaded (alpha carbon only) and the biological assmbly generated and displayed in a gizmo

![img](img/p_pdb3iyv0_basic.png?raw=true "image")

The key feature of this code is that a simplified view of a biological assembly can be generated by "down-sampling" the backbone. An example is shown here using `DownsampleStructure(5,25)` as shown in the ipf file `Figure3IYV.ipf` which generates some images I was interested in from the PDB file of clathrin assembled into a cage.

![img](img/p_pdb3iyv0.png?raw=true "image")
![img](img/p_pdb3iyv1.png?raw=true "image")

Some other views are available in the `img` directory.