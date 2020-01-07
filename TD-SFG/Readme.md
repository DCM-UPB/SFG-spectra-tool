# simulate-SFG
This program allows one to compute the time-dependent vibrational SFG spectrum for a given planar interface system. Using this approach, one can estimate the transistion time from one type of stretching mode to another. Theory of this method can be found in the following references [1,2,3]. 

## How this program works?
**System requirements - C++11, Open MPI libraries**


Copy your xyz and vibration.dat files to the Trajectory folder. 

```
cp -r your_files Trajectory/ 
```
**Important note: The xyz file should contain only water molecules. The water molecules should be sorted in an order, where the oxygen atom coordinate followed by its two neighbour hydrogen atoms. For reference
please check this sample file [planar.xyz](https://github.com/DCM-UPB/SFG-spectrum/blob/master/SFG-planar-interface/Trajectory/planar.xyz).**


Define the trajectory length, number of atoms, trajectory file name and location,
length of the simulation box, time step, time correlation function length and the interface regions in the SRC/input file.
```
vim SRC/input
```
Compile and Execute the program.  
```
chmod +x simulate-SFG.sh
```
```
./simulate-SFG.sh
```
Plot the SFG spectrum (.SFG), density profile (.DP) simulated in the Trajectory folder.

## References

1. [T. Ohto, K. Usui, T. Hasegawa, M Bonn, and Y Nagata
(2015)
Toward ab initio molecular dynamics modeling for sum-frequency generation spectra; an efficient algorithm based on surface-specific velocity-velocity correlation function,
Chemical Physics,
DOI: https://doi.org/10.1063/1.4931106](https://aip.scitation.org/doi/citedby/10.1063/1.4931106)


2. [N. K. Kaliannan, A. H. Aristizabal, H. Wiebeler, F. Zysk, T. Ohto, Y. Nagata and T. D. Kühne 
(2019)
Impact of intermolecular vibrational coupling effects on the sum-frequency generation spectra of the water/air interface,
Molecular Physics,
DOI: 10.1080/00268976.2019.1620358](https://www.tandfonline.com/doi/full/10.1080/00268976.2019.1620358)

3. [D. Ojha, N. K. Kaliannan and T. D. Kühne 
(2019)
Time-dependent vibrational sum-frequency generation spectroscopy of the air-water interface,
Communications Chemistry,
DOI: 10.1038/s42004-019-0220-6](https://doi.org/10.1038/s42004-019-0220-6)


