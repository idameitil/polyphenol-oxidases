Go to 

# Periodic boundary condition (getting out of Pacman world), center
printf "Protein\nSystem\n" | gmx trjconv -s md.tpr -f md.xtc -o md_noPBC.xtc -pbc mol -center

# Orient the protein the same way in every frame (rotation)
printf "Protein\nSystem\n" | gmx trjconv -s md.tpr -f md_noPBC.xtc -o md_fit.xtc -fit progressive

# Download md_fit_AoCO4_wt_CuPerGua.xtc (trajectory)

# Open in pymol
pymol md_AoCO4_wt_CuPerGua.gro
# file -> open -> md_fit_AoCO4_wt_CuPerGua.xtc
sele resi 102+110+119+284+288+312+384+385

# Start tmux
tmux a -t gmxmd

# Detach from tmux: ctrl b + d

From carp-rtp, copy S-lignin block into lignin.rtp

From top_lignin.top, copy the molecule figure into lignin.rtp

Modify the molecule into syringol (replace alcohol part with a hydrogen, H1)

Remove every atom line with atoms that were removed and add the new hydrogen. Remove bonds that were removed and add the new hydrogen bond.
Set charge of the new H to 0.115 and the charge of the C to -0.115.

Change name to syringol and syro.

.rtp is Gromacs format.

lignin.hdb, hydrogens to be added. Gromacs removes protons. We have to add them again.

1st column, how many protons. 2nd column is H type. 3rd column is which one it is connected to. 4th and 5th columns are which ones they oriented by.

# Docking
conda activate vina

Docking output pymol, choose the lowest energy state, export molecule

# MD
module swap gromacs/2021.2-gpu

Change atom names of syringol pdb (mol), C1, C2 etc. Change UNK to Syro

Run ./prepare_md.sh