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
double angle_btwn_3points(vector<double> const &r,int i,int t, double a,double b,double c, unsigned int N)
{
  double x1,x2;
  double y1,y2;
  double z1,z2;

  x1 = r[N*t*d+(i+0)*d +0] - 0;
  x2 = 0 - 0;
  y1 = 0;
  y2 = 0 - 0;
  z1 = r[N*t*d+(i+0)*d +2] - 0;
  z2 = 0 - 1; 

  x1 = x1 - a * round(x1/a);
  x2 = x2 - a * round(x2/a);
  y1 = y1 - b * round(y1/b);
  y2 = y2 - b * round(y2/b);
  z1 = z1 - c * round(z1/c);
  z2 = z2 - c * round(z2/c);

  double top =  (x1*x2+y1*y2+z1*z2) - 0.000001 ;
  double bot =  norm(x1,y1,z1) * norm(x2,y2,z2);

  if( top == bot )
    {
      return 180; 
    }
  else
    {
      return  180 - (acos(top/bot) * 57.296);
    }
}


//----------------------------- Main implementation -----------------------------//
int main(int argc, char** argv)
{


  unsigned int  Traj_len = 0, N = 0, N_strings;                               // N - number of atoms in the given system
  float dt = 0;                                                    // dt - time step in femtosecond
  vector<double>  L(d, 0.0);                                       // r - atomic position, L - length of the simulation box, v - atomic velocity
  string dummy, filename1, filename2;
  char  dummy_;


  ifstream input("../SRC/input");
  input >> dummy >> Traj_len  ; 
  input >> dummy >> N         ; 
  input >> dummy >> filename1 ; 
  input >> dummy >> N_strings ; 
  input >> dummy >> filename2 ; 
  input >> dummy >> L[0] >> L[1] >> L[2] ; 
  input >> dummy >> dt        ;  
  input.close();
  input.clear();


  vector<double> r (N*d*Traj_len, 0.0), v(N*d*Traj_len, 0.0);       // r - atomic position, v - atomic velocity
  vector<double> r_ (N*d*Traj_len, 0.0), v_(N*d*Traj_len, 0.0);     // r_, v_ transfored atomic coordinates and velocities
  vector<double> comx (Traj_len, 0.0), comy (Traj_len, 0.0), comz (Traj_len, 0.0);       // com - center of mass
  vector<char> atom(N*Traj_len);                                    // atom type - atom name O or H
  vector<int> Layer(N*Traj_len, 0.0);                               // Layer - location of the molecule from the instantaneous water surface

  ifstream infile(filename1);
  for(unsigned int t = 0;t < Traj_len;++t)
    {
      for (unsigned int words = 0 ; words <  N_strings + 1  ; words++)
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
      // Normalized to the origin
      for(unsigned int i = 0;i < N;++i)
        {
          r[N*t*d+i*d+0] = r[N*t*d+i*d+0] - comx[t] ;
          r[N*t*d+i*d+1] = r[N*t*d+i*d+1] - comy[t] ; 
          r[N*t*d+i*d+2] = r[N*t*d+i*d+2] - comz[t] ;
        }
    }

  // radial density profie calculation
  unsigned int RDP_size = 1000, tframe = 40;
  double RDP_h = 0.3; //(h-stepsize) 
  vector<double> z (RDP_size, 0.0), Density_total(RDP_size, 0.0);       // r - atomic position, v - atomic velocity
  for(unsigned int i = 1;i < RDP_size;++i)
    {
      z[i]          = i * 0.5;
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
              rij = norm(r[N*t*d+j*d+0],r[N*t*d+j*d+1]-r[N*t*d+j*d+1],r[N*t*d+j*d+2]);;
              if(rij <= z[i] + (RDP_h/2) && rij > z[i] - (RDP_h/2) && atom[j] == 'H') //Hydrogen
                {
                  N_H = N_H + 1;
                }
              else if(rij <= z[i] + (RDP_h/2) && rij > z[i] - (RDP_h/2) && atom[j] == 'O')//Oxygen
                {
                  N_O = N_O + 1;
                }
            }
          Density_total[i] += (N_H * 1.00  + N_O * 15.99) * 1.66  / (PI * L[1] * (pow(z[i]+RDP_h,2) - pow(z[i],2)) );
        }}
    }


  // velocity is determined using minimum image convention
  for(unsigned int t = 0;t < (Traj_len - 1);++t)
    {
      for(unsigned int i = 0;i < N;++i)
        {
          v[N*t*d+i*d+0] = min_distance(r[N*(t+1)*d+i*d+0] - r[N*t*d+i*d+0], L[0]) / dt;
          v[N*t*d+i*d+1] = min_distance(r[N*(t+1)*d+i*d+1] - r[N*t*d+i*d+1], L[1]) / dt;
          v[N*t*d+i*d+2] = min_distance(r[N*(t+1)*d+i*d+2] - r[N*t*d+i*d+2], L[2]) / dt;
        }           
    }



  // transformed coordinates and velocity
  double angle = 0;
  for(unsigned int t = 0;t < (Traj_len - 1);++t)
    {
      for(unsigned int i = 0;i < N;i = i + 3)
        {
          double x1,x2;
          double y1,y2;
          double z1,z2;

          x1 = r[N*t*d+i*d+0] - 0;
          x2 = 0 - 0.;
          z1 = r[N*t*d+i*d+2] - 0;
          z2 = 0 - 1;  

          x1 = x1 - L[0] * round(x1/L[0]);
          x2 = x2 - L[0] * round(x2/L[0]);
          z1 = z1 - L[2] * round(z1/L[2]);
          z2 = z2 - L[2] * round(z2/L[2]);

          double top =  (x1*x2+z1*z2) - 0.0001 ;
          double bot =  norm(x1,0,z1) * norm(x2,0,z2)  ;            
          if(r[N*t*d+i*d+0] >=  0){ angle = 180 - acos(top/bot) * 57.296 ; }
          else if (r[N*t*d+i*d+0] < 0){angle = acos(top/bot) * 57.296 - 180; }
          if(isnan(angle)) {angle = 180;}

          for(unsigned int i_ = 0; i_ < 3 ; i_ = i_ + 1)
            {
              r_[N*t*d+(i+i_)*d+2] = r[N*t*d+(i+i_)*d+2] * cos ( angle * PI / 180.0 )       + r[N*t*d+(i+i_)*d+0] * sin ( angle * PI / 180.0 );
              r_[N*t*d+(i+i_)*d+0] = r[N*t*d+(i+i_)*d+2] * -1 * sin ( angle * PI / 180.0 )  + r[N*t*d+(i+i_)*d+0] * cos ( angle * PI / 180.0 );
              r_[N*t*d+(i+i_)*d+1] = r[N*t*d+(i+i_)*d+1];

              v_[N*t*d+(i+i_)*d+2] = v[N*t*d+(i+i_)*d+2] * cos ( angle * PI / 180.0 )       + v[N*t*d+(i+i_)*d+0] * sin ( angle * PI / 180.0 ); 
              v_[N*t*d+(i+i_)*d+0] = v[N*t*d+(i+i_)*d+2] * -1 * sin ( angle * PI / 180.0 )  + v[N*t*d+(i+i_)*d+0] * cos ( angle * PI / 180.0 );
              v_[N*t*d+(i+i_)*d+1] = v[N*t*d+(i+i_)*d+1];
            }
        }           
    }


  ofstream outfile(filename2);
  for(unsigned int t = 0;t < (Traj_len - 1);++t)
    {
      outfile << (2*N)/3 << endl;
      for(unsigned int i = 0;i < N;i = i + 3)
        {
          double dist[4], velo[3];
          dist[0] =  r_[N*t*d+(i+1)*d+0] - r_[N*t*d+i*d+0];dist[1] = r_[N*t*d+(i+1)*d+1] - r_[N*t*d+i*d+1];dist[2] = r_[N*t*d+(i+1)*d+2] - r_[N*t*d+i*d+2];
          velo[0] =  v_[N*t*d+(i+1)*d+0] - v_[N*t*d+i*d+0];velo[1] = v_[N*t*d+(i+1)*d+1] - v_[N*t*d+i*d+1];velo[2] = v_[N*t*d+(i+1)*d+2] - v_[N*t*d+i*d+2];  
          dist[0] =  dist[0];dist[1] = dist[1];dist[2] = dist[2];
          dist[3] =  norm(dist[0],dist[1],dist[2]); 
          dist[0] = dist[0] / dist[3]; dist[1] = dist[1] / dist[3]; dist[2] = dist[2] / dist[3];
          outfile << velo[2] << "  " << dp(dist[0],dist[1],dist[2],velo[0],velo[1],velo[2]) << "  " << r[N*t*d+(i+1)*d+0] << "  " << r[N*t*d+(i+1)*d+1] << "  " <<  r[N*t*d+(i+1)*d+2] << "  " <<r[N*t*d+i*d+0] << "  " << r[N*t*d+i*d+1] << "  " << r[N*t*d+i*d+2] << "  " << Layer[N*t+i] <<  endl;


          dist[0] =  r_[N*t*d+(i+2)*d+0] - r_[N*t*d+i*d+0];dist[1] = r_[N*t*d+(i+2)*d+1] - r_[N*t*d+i*d+1];dist[2] = r_[N*t*d+(i+2)*d+2] - r_[N*t*d+i*d+2];
          velo[0] =  v_[N*t*d+(i+2)*d+0] - v_[N*t*d+i*d+0];velo[1] = v_[N*t*d+(i+2)*d+1] - v_[N*t*d+i*d+1];velo[2] = v_[N*t*d+(i+2)*d+2] - v_[N*t*d+i*d+2];  
          dist[0] =  dist[0];dist[1] = dist[1];dist[2] = dist[2];
          dist[3] =  norm(dist[0],dist[1],dist[2]);
          dist[0] = dist[0] / dist[3]; dist[1] = dist[1] / dist[3]; dist[2] = dist[2] / dist[3]; 
          outfile << velo[2] << "  " << dp(dist[0],dist[1],dist[2],velo[0],velo[1],velo[2]) << "  " << r[N*t*d+(i+2)*d+0] << "  " << r[N*t*d+(i+2)*d+1] << "  " <<  r[N*t*d+(i+2)*d+2] << "  " <<r[N*t*d+i*d+0] << "  " << r[N*t*d+i*d+1] << "  " << r[N*t*d+i*d+2] << "  " << Layer[N*t+i] <<  endl;
        }          
    } 
  outfile.close();
  outfile.clear();

  filename2.append(".DP");
  ofstream outfile_(filename2);
  for(unsigned int i = 1; i < RDP_size; ++i)
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




