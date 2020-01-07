/**
 * This program computes the sum frequency generation (SFG)
 * spectrum (both imaginary and real terms) for the given 
 * interfacial water molecules. 
 * 
 *
 * @Author  Naveen Kumar Kaliannan
 * Reach me via naveenkumar5892@gmail.com
 * 
 */



#include<mpi.h>
#include<iostream>
#include<fstream>
#include<string>
#include<cmath>
#include<typeinfo>
#include<cstdlib>
#include<vector>

using namespace std;


// Constants
const float PI      = 3.14592; // PI constant value
const int d         = 3;       // Dimension
const int Nfreq     = 5000;
const double light_vel = 299792458E-13;


double min_distance(double r, double l){  return r - l * round(r/l);}
double norm(double x,double y,double z){  return sqrt(x*x+y*y+z*z);}
double dp(double x,double y,double z,double xx,double yy,double zz){  return x*xx+y*yy+z*zz;}
double mindis(double dx,double dy,double dz,double a,double b){  return norm(dx - a * round(dx/a),dy - b * round(dy/b), dz);}
void FFT(vector<double> &vvacf,vector<double> & vvacf_r_fft,vector<double> & vvacf_i_fft, float tcfl, float dt)
{
  for(unsigned int freq = 0; freq < Nfreq;++freq)
    {
      for(unsigned int t_ = 0;t_ < tcfl;++t_)
        {
          vvacf_r_fft[freq] = vvacf_r_fft[freq] + vvacf[t_] * cos(2.0*PI*light_vel*t_*dt*freq) * dt;
          vvacf_i_fft[freq] = vvacf_i_fft[freq] + vvacf[t_] * sin(2.0*PI*light_vel*t_*dt*freq) * dt;
        }
    }
  
}


// Main implementation
int main(int argc, char** argv)
{

  unsigned int  Traj_len = 0, N = 0, tcfl = 0,  N_strings;                   // Traj_len = trajectory len, N - number of atoms in the given system, tcfl - time correlation function length
  float dt = 0, T1 = 0, T2 = 0, B1 = 0, B2 = 0, T_w = 0;       // dt - time step in femtosecond, rcut - cross-correlation cutoff T, B - surface region
  vector<double>  L(d, 0.0);                                     // r - atomic position, L - length of the simulation box, v - atomic velocity
  string dummy, filename1, filename2, filename3, filename4;
  char  dummy_;


  ifstream input("../SRC/input");
  input >> dummy >> Traj_len  ; 
  input >> dummy >> N         ; 
  input >> dummy >> filename1 ; 
  input >> dummy >> N_strings ; 
  input >> dummy >> filename2 ; 
  input >> dummy >> L[0] >> L[1] >> L[2] ; 
  input >> dummy >> dt        ; 
  input >> dummy >> filename4 ; 
  input >> dummy >> tcfl    ;   
  input >> dummy >> filename3 ; 
  input >> dummy >> T1 >> dummy >> T2 ;
  input >> dummy >> B1 >> dummy >> B2 ; 
  input >> dummy >> T_w ; 
  input.close();
  input.clear();


  int size, rank, root = 0; 
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  Traj_len = Traj_len - 1 ;
  Traj_len = size * floor(Traj_len / size) ;                       // Trajectory length is slightly adjusted based on the number of processors
  N  = ( N * 2 ) / 3;                                              // Number of OH group is determined

  vector<double> vvacf(tcfl,0.0);
  vector<double> vvacf_r_fft(Nfreq,0);
  vector<double> vvacf_i_fft(Nfreq,0);
  vector<double> vz(N * Traj_len, 0);
  vector<double> v_proj(N * Traj_len, 0); 
  vector<double> H_coord_x(N * Traj_len, 0); 
  vector<double> H_coord_y(N * Traj_len, 0);   
  vector<double> H_coord_z(N * Traj_len, 0);     
  vector<double> O_coord_x(N * Traj_len, 0); 
  vector<double> O_coord_y(N * Traj_len, 0);   
  vector<double> O_coord_z(N * Traj_len, 0);  
  vector<int> Layer(N * Traj_len, 0);   
  vector<double> vibration(N * Traj_len, 0);   


  unsigned int from, to, divide; 
  divide =  tcfl + T_w + 1000 ;
  from = (rank * (Traj_len - divide ) / size) ;
  to = ((rank + 1) * (Traj_len - divide) / size);
  if(rank == root)
    {
      ifstream infile(filename2);
      for(unsigned int t = 0;t < Traj_len;++t)
        {
          infile >> N ;
          for(unsigned int i = 0;i < N;++i)
            {
              infile >> vz[N*t+i] >> v_proj[N*t+i] >> H_coord_x[N*t+i] >> H_coord_y[N*t+i] >> H_coord_z[N*t+i] >> O_coord_x[N*t+i] >> O_coord_y[N*t+i] >> O_coord_z[N*t+i] >> Layer[N*t+i] >> vibration[N*t+i];
            }           
        }
      infile.close();
      infile.clear();   
    }


  MPI_Barrier(MPI_COMM_WORLD);  
  MPI_Bcast(&vvacf[0], vvacf.size(), MPI_DOUBLE, root, MPI_COMM_WORLD); 
  MPI_Bcast(&vz[0], vz.size(), MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&v_proj[0], v_proj.size(), MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&H_coord_x[0], H_coord_x.size(), MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&H_coord_y[0], H_coord_y.size(), MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&H_coord_z[0], H_coord_z.size(), MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&O_coord_x[0], O_coord_x.size(), MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&O_coord_y[0], O_coord_y.size(), MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&O_coord_z[0], O_coord_z.size(), MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&Layer[0], Layer.size(), MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(&vibration[0], vibration.size(), MPI_DOUBLE, root, MPI_COMM_WORLD);  


  double mean_ = 0;
  for(unsigned int t = from;t < to;++t)
    { 
      mean_ = 0;
      for(unsigned int i = 0;i < N;++i)
        { 
          double vibration_cut = 3744.5; // This value changes depending on the water model used.
          if( O_coord_z[N*t+i] >  T1  && O_coord_z[N*t+i] < T2 && vibration[N*t+i] < vibration_cut ) //- for xyz file format
          ///if( Layer[N*t+i] == 1 && O_coord_z[N*t+i] > 23 && vibration[N*t+i] > vibration_cut ) // - for pdb file format - top side water molecules
            { 
              mean_ += 1; 
              for(unsigned int t_ = 0;t_ < tcfl;++t_)
                {
                  vvacf[t_] =  vvacf[t_] + vz[N*(t+T_w)+i] *  v_proj[N*(t+t_+T_w)+i];          //---> Auto correlation
                }
            }  
          else if( O_coord_z[N*t+i] >  B1  && O_coord_z[N*t+i] < B2 && vibration[N*t+i] < vibration_cut ) // - for xyz file format
           //else if( Layer[N*t+i] == 1 && O_coord_z[N*t+i] < 23 && vibration[N*t+i] > vibration_cut ) //- for pdb file format - bottom side water molecules
            {
              mean_ += 1;
              for(unsigned int t_ = 0;t_ < tcfl;++t_)
                {
                  vvacf[t_] =  vvacf[t_] - vz[N*(t+T_w)+i] *  v_proj[N*(t+t_+T_w)+i];          //---> Auto correlation
                }
            }  
         } 
    }


  MPI_Barrier(MPI_COMM_WORLD);  
  double mean_global = 0;
  MPI_Reduce(&mean_, &mean_global, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
  mean_global = mean_global / size;
  for(unsigned int t_ = 0;t_ < tcfl; ++t_)
    {
      double local = 0, global = 0;
      local = vvacf[t_]; 
      MPI_Reduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
      if(rank == root)
        {
          vvacf[t_] = global;
        }
    }


  if(rank == root)
    {
      filename2.append(".vvcf");
      ofstream outfile_1(filename2);
      for(unsigned int t_ = 0;t_ < tcfl;++t_)
        {
          outfile_1 << t_ * dt << "  " << vvacf[t_] << "  " << mean_global << "  ";
          vvacf[t_] = vvacf[t_] * pow(cos((PI) * t_ * dt / (tcfl * dt * 2.0)),2) / mean_global; 
          outfile_1 << vvacf[t_] << "\n";
        }
      outfile_1.close();
      outfile_1.clear();

      // the imaginary and real terms are exchanged due to the imaginary unit i infront of the resonant nonlinear susceptibility
      FFT(vvacf,vvacf_i_fft,vvacf_r_fft, tcfl, dt);

      ofstream outfile_2(filename3);
      outfile_2 << " #Freq   imaginary-term   Real-term \n";
      for(unsigned int freq = 0;freq < Nfreq;++freq)
        {  //         frequency         imaginary term                real term (negative sign is due to the imaginary unit i infront of the susceptibility)
          outfile_2 << freq << "   " << vvacf_i_fft[freq] << "   " << (-1) * vvacf_r_fft[freq] << "\n";
          // Quantum correction factor and non-condon approximation are not included here. 
        }
      outfile_2.close();
      outfile_2.clear();
    }


  vz.clear();vz.shrink_to_fit(); 
  v_proj.clear();v_proj.shrink_to_fit();  
  H_coord_x.clear();H_coord_x.shrink_to_fit();  
  O_coord_x.clear();O_coord_x.shrink_to_fit();
  H_coord_y.clear();H_coord_y.shrink_to_fit();  
  O_coord_y.clear();O_coord_y.shrink_to_fit();
  H_coord_z.clear();H_coord_z.shrink_to_fit();  
  O_coord_z.clear();O_coord_z.shrink_to_fit();
  Layer.clear();Layer.shrink_to_fit();
  vibration.clear();vibration.shrink_to_fit();
  vvacf.clear();vvacf.shrink_to_fit();
  vvacf_r_fft.clear();vvacf_r_fft.shrink_to_fit();
  vvacf_i_fft.clear();vvacf_i_fft.shrink_to_fit();


  MPI_Finalize();  
  return 0;
}
