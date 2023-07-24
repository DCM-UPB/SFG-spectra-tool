rm bin/*
cd SRC
module load mpi
make
cd ../bin 
./SFG_Mu_Alpha.exe 
mpirun -np 1 ./SFG_Spectra.exe 




