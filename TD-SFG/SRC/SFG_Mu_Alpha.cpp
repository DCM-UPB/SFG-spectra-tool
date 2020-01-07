/**
 * This program computes the dipole moment and
 * bond polarisability for given water molecules 
 * using the OH stretching velocities
 *
 * @Author  Naveen Kumar Kaliannan
 * Reach me via naveenkumar5892@gmail.com
 * 
 */

#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<cmath>
#include<typeinfo>
#include<cstdlib>


const float PI      = 3.14592; // PI constant value
const int d         = 3;       // Dimension


using namespace std;



double dp(double x,double y,double z,double xx,double yy,double zz){  return x*xx+y*yy+z*zz;}
double min_distance(double r, double l){  return r - l * round(r/l);}
double norm(double x,double y,double z){  return sqrt(x*x+y*y+z*z);}


//----------------------------- Main implementation -----------------------------//
int main(int argc, char** argv)
{


  unsigned int  Traj_len = 0, N = 0, N_strings;                               // N - number of atoms in the given system
  float dt = 0;                                                    // dt - time step in femtosecond
  vector<double>  L(d, 0.0);                                       // r - atomic position, L - length of the simulation box, v - atomic velocity
  string dummy, filename1, filename2, filename3;
  char  dummy_;


  ifstream input("../SRC/input");
  input >> dummy >> Traj_len  ; 
  input >> dummy >> N         ; 
  input >> dummy >> filename1 ; 
  input >> dummy >> N_strings ; 
  input >> dummy >> filename2 ; 
  input >> dummy >> L[0] >> L[1] >> L[2] ; 
  input >> dummy >> dt        ; 
  input >> dummy >> filename3        ;   
  input.close();
  input.clear();


  vector<double> r (N*d*Traj_len, 0.0), v(N*d*Traj_len, 0.0);       // r - atomic position, L - length of the simulation box, v - atomic velocity
  vector<double> comx (Traj_len, 0.0), comy (Traj_len, 0.0), comz (Traj_len, 0.0);       // com - center of mass
  vector<char> atom(N*Traj_len);                                    // atom type - atom name O or H
  vector<int> Layer(N*Traj_len, 0.0);                                  // Layer - location of the molecule from the instantaneous water surface


//-----pdb file format-----///
/**
  ifstream infile(filename1);
  for(unsigned int words = 0; words < 10; ++words)
    {
      infile >> dummy; 
    }
  for(unsigned int t = 0;t < Traj_len;++t)
    {

      infile >> dummy >> dummy >> dummy >> dummy;
      for(unsigned int words = 0; words < 8; ++words)
        {
          infile >> dummy;
        }
     for(unsigned int i = 0; i <  N; ++i)
       {
         infile >>  dummy; 
         infile >>  dummy; 
         infile >> atom[N*t+i] ;  
         if(atom[N*t+i] == 'O')
           {
            infile >> dummy >> dummy >> dummy ;
           }
         else if (atom[N*t+i] == 'H')
           {
             infile >> dummy >> dummy >> dummy ;
           } 
         infile >>r[N*t*d+i*d+0] >> r[N*t*d+i*d+1] >> r[N*t*d+i*d+2] >> dummy; //cout << r[N*t*d+i*d+0]  << " " <<  r[N*t*d+i*d+1] << " " <<  r[N*t*d+i*d+2] << endl;

      infile >> Layer[N*t+i] >> dummy >> dummy >> dummy; 
    } 
  infile >> dummy ; 
    }
  infile.close();
  infile.clear();
**/

//-----xyz file format-----///

  ifstream infile(filename1);
  for(unsigned int t = 0;t < Traj_len;++t)
    {
      for (unsigned int words = 0 ; words < N_strings + 1 ; words++)
        {
      	  infile  >> dummy; 
        }
      for(unsigned int i = 0; i <  N; ++i)
        {
          infile >>  atom[N*t+i] >> r[N*t*d+i*d+0] >> r[N*t*d+i*d+1] >> r[N*t*d+i*d+2];
 
        }  
    }
  infile.close();
  infile.clear();



  // center of mass 
  for(unsigned int t = 0;t < Traj_len;++t)
    {
      for(unsigned int i = 0;i < N;++i)
        {
          comx[t] = comx[t] + r[N*t*d+i*d+0];
          comy[t] = comy[t] + r[N*t*d+i*d+1]; 
          comz[t] = comz[t] + r[N*t*d+i*d+2];
        }
      comx[t] = comx[t] / N ;
      comy[t] = comy[t] / N ;
      comz[t] = comz[t] / N ; 
      // bringing the system to the origin in the surface direction, i.e., the z-direction
      for(unsigned int i = 0;i < N;++i)
        {
          r[N*t*d+i*d+2] = r[N*t*d+i*d+2] - comz[t] ;
        }
    }

  // radial density profie calculation
  unsigned int RDP_size = 1000, tframe = 30;
  double RDP_h = 0.3; //(h-stepsize) 
  vector<double> z (RDP_size, 0.0), Density_total(RDP_size, 0.0);
  for(unsigned int i = 0;i < RDP_size;++i)
    {
      z[i]          =  i * 0.1 - 60.0;
    }  
  for(unsigned int t = 0;t < Traj_len;++t)
    {
      if(t%tframe == 1){
      for(unsigned int i = 1;i < RDP_size;++i)
        {
          unsigned int N_H = 0,N_O = 0;
          double rij = 0;
          for(unsigned int j = 0;j < N;++j)
            {
              rij = r[2+j*d];
              if(rij <= z[i] + (RDP_h/2) && rij > z[i] - (RDP_h/2) && atom[j] == 'H') // counting H
                {
                  N_H = N_H + 1;
                }
              else if(rij <= z[i] + (RDP_h/2) && rij > z[i] - (RDP_h/2) && atom[j] == 'O')// counting O
                {
                  N_O = N_O + 1;
                }
            }
          Density_total[i] += (N_H * 1.00  + N_O * 15.99) * 1.66  / (L[0] * L[1] * RDP_h);
        }}
    }

  // determining the velocity 
  for(unsigned int t = 0;t < (Traj_len - 1);++t)
    {
      for(unsigned int i = 0;i < N;++i)
        {
          //Minimum image convention is applied only in x- and y-directions. 
          //Free boundary condition is applied in the surface direction, i.e., the z-direction
          v[N*t*d+i*d+0] = min_distance(r[N*(t+1)*d+i*d+0] - r[N*t*d+i*d+0], L[0]) / dt;
          v[N*t*d+i*d+1] = min_distance(r[N*(t+1)*d+i*d+1] - r[N*t*d+i*d+1], L[1]) / dt;
          v[N*t*d+i*d+2] = ( r[N*(t+1)*d+i*d+2] - r[N*t*d+i*d+2] ) / dt;
        }           
    }

  // determining the dipole moment and polarisability term
  ofstream outfile(filename2);
  infile.open(filename3);
  double vibration = 0;
  for(unsigned int t = 0;t < (Traj_len - 1);++t)
    {
      outfile << (2*N)/3 << endl;
      for(unsigned int i = 0;i < N;i = i + 3)
        {
          double dist[4], velo[3];
          dist[0] =  r[N*t*d+(i+1)*d+0] - r[N*t*d+i*d+0];dist[1] = r[N*t*d+(i+1)*d+1] - r[N*t*d+i*d+1];dist[2] = r[N*t*d+(i+1)*d+2] - r[N*t*d+i*d+2];
          velo[0] =  v[N*t*d+(i+1)*d+0] - v[N*t*d+i*d+0];velo[1] = v[N*t*d+(i+1)*d+1] - v[N*t*d+i*d+1];velo[2] = v[N*t*d+(i+1)*d+2] - v[N*t*d+i*d+2];  
          //Minimum image convention is applied only in x- and y- directions
          dist[0] =  min_distance(dist[0], L[0]);dist[1] = min_distance(dist[1], L[1]);dist[2] = dist[2];
          dist[3] =  norm(dist[0],dist[1],dist[2]);
          dist[0] = dist[0] / dist[3]; dist[1] = dist[1] / dist[3]; dist[2] = dist[2] / dist[3];
          infile >> vibration;           //reading the vibration of the first OH group in the water molecules
          outfile << velo[2] << "  " << dp(dist[0],dist[1],dist[2],velo[0],velo[1],velo[2]) << "  " << r[N*t*d+(i+1)*d+0] << "  " << r[N*t*d+(i+1)*d+1] << "  " <<  r[N*t*d+(i+1)*d+2] << "  " <<r[N*t*d+i*d+0] << "  " << r[N*t*d+i*d+1] << "  " << r[N*t*d+i*d+2] << "  " << Layer[N*t+i] << "  " << vibration <<  endl;


          dist[0] =  r[N*t*d+(i+2)*d+0] - r[N*t*d+i*d+0];dist[1] = r[N*t*d+(i+2)*d+1] - r[N*t*d+i*d+1];dist[2] = r[N*t*d+(i+2)*d+2] - r[N*t*d+i*d+2];
          velo[0] =  v[N*t*d+(i+2)*d+0] - v[N*t*d+i*d+0];velo[1] = v[N*t*d+(i+2)*d+1] - v[N*t*d+i*d+1];velo[2] = v[N*t*d+(i+2)*d+2] - v[N*t*d+i*d+2];  
          //Minimum image convention is applied only in x- and y- directions
          dist[0] =  min_distance(dist[0], L[0]);dist[1] = min_distance(dist[1], L[1]);dist[2] = dist[2];
          dist[3] =  norm(dist[0],dist[1],dist[2]);
          dist[0] = dist[0] / dist[3]; dist[1] = dist[1] / dist[3]; dist[2] = dist[2] / dist[3]; 
          infile >> vibration;           //reading the vibration of the second OH group in the water molecules
          outfile << velo[2] << "  " << dp(dist[0],dist[1],dist[2],velo[0],velo[1],velo[2]) << "  " << r[N*t*d+(i+2)*d+0] << "  " << r[N*t*d+(i+2)*d+1] << "  " <<  r[N*t*d+(i+2)*d+2] << "  " <<r[N*t*d+i*d+0] << "  " << r[N*t*d+i*d+1] << "  " << r[N*t*d+i*d+2] << "  " << Layer[N*t+i] << "  " << vibration <<  endl;
        }          
    } 
  infile.close();
  infile.clear();
  outfile.close();
  outfile.clear();

  filename2.append(".DP");
  ofstream outfile_(filename2);
  for(unsigned int i = 0; i < RDP_size; ++i)
    { 
      outfile_ << z[i] <<" \t "<< Density_total[i] /(Traj_len/tframe) << "\n";
    }
  outfile_.close();
  outfile_.clear(); 


  v.clear();v.shrink_to_fit();  
  r.clear();r.shrink_to_fit();  
  atom.clear();atom.shrink_to_fit();  
  L.clear();L.shrink_to_fit();  
  Layer.clear();Layer.shrink_to_fit(); 
 
 
  return 0;
}



