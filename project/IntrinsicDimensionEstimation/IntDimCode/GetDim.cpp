#include <stdlib.h>
#include <math.h>
#include <iostream.h>
#include <string.h>

long num;
long dim;
long dim_right;
float CORR_DIM;
float DIM;
float TAK_DIM;
long** sp_index;
long* nr_elements;

long ctr=4;
long ctr2=0;

/* uniform [0,1] random number generator
   developed by Pierre Lecuyer based on a clever
   and tested combination of two linear congruential
   sequences 
*/

/*
s1 and s2 are the seeds (nonnegative integers)
*/

double uni()
{
static long s1 = 55555;
static long s2 = 99999;
static double factor = 1.0/2147483563.0;
register long k,z;
	k= s1 /53668;
	s1 =40014*(s1%53668)-k*12211;
	if (s1 < 0) s1 += 2147483563;
	k=s2/52774;
	s2=40692*(s2%52774)-k*3791;
	if (s2 < 0) s2 += 2147483399;

	/*
	z = abs(s1 ^ s2);
	*/
	z= (s1 - 2147483563) + s2;
	if (z < 1) z += 2147483562;

	return(((double)(z))*factor);
}


float** generate_data(long data_set)
{
  //this part generates toy training data
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // choice of dataset
  // 0: a circle in 2d
  // 1: a sphere in 3d
  // 2: a 3d-affine space in 5d
  // 3: a strange 4-dim figure in 5d (but very concentrated so it is ess 3d)
  // 4: a 4-dim manifold in 8 dimensions
  // 5: a 2d helix in 3d
  // 6: a 6-dim manifold in 36 dimensions
  // 7: "swiss roll"
  // 8: 12-dim manifold in 72 dimensions
  // 9: A 20-dim affine space in 100-dimensions
  //10: k-hypercube, uniformly sampled
  //11: Moebius band
  //12: Multivariate Gaussian
  //13: one-dimensional curve in dim Dimensions


  // the data vector
  float** data;

  long i,j,k;
  long num_read;
  FILE* trainfile;
  long file_length;
  unsigned char* image;
  float* image_fl;

  float* radii;
  float* phi;
  float* theta;

  float norm;

  long* labeled;

  long nr_par;
  float** vectors;
  float** para;
  float* para1; float* para2; float* para3;
  double** para4;

  switch(data_set)
  {
    case 0:  //generates a (slightly distorted) sinusoid over the circle in R^3
      dim=3;
      dim_right=1;

      data = new float*[num];
      for(i=0;i<num;i++)
       data[i]=new float[dim];

      // generate random angles and radii for the two moons
      radii=new float[num];
      phi  =new float[num];

      for(i=0;i<num;i++)
      {
        radii[i]=1;//+0.05*uni();
        phi[i]  =uni()*2*M_PI;  // random numbers in [0,2pi]
      }

      // the circle
      for(i=0;i<num;i++)
      {
        data[i][0]=radii[i]*cos(phi[i]);
        data[i][1]=radii[i]*sin(phi[i]);
        data[i][2]=0.1*sin(150*phi[i]);
        //data[i][2]=0.2*uni()-0.1;
      }
      delete radii;
      delete phi;
      break;

   case 1:  // generates a k-sphere  (uniformly sampled)
      dim_right=dim-1;

      nr_par=dim;
      if(nr_par %2 ==1) nr_par+=1;

      data = new float*[num];
      for(i=0;i<num;i++)
       data[i]=new float[dim];

      para=new float*[nr_par];
      for(i=0;i<nr_par;i++)
       para[i]=new float[num];

      for(j=0;j<nr_par;j++)
       for(i=0;i<num;i++)
       {
         para[j][i]=uni();
         while(para[j][i]==0)
          para[j][i]=uni();
       }


      for(i=0;i<num;i++)
       { for(j=0;j<dim;j++) data[i][j]=0; }
      k=dim;
      for(i=0;i<num;i++)
      {
        for(j=0;j<dim-2;j+=2)
        {
          data[i][j  ]=sqrt(-2*log(para[j][i]))*cos(2*M_PI*para[j+1][i]);
          data[i][j+1]=sqrt(-2*log(para[j][i]))*sin(2*M_PI*para[j+1][i]);
        }
        if(dim % 2==0)
        {
          data[i][j  ]=sqrt(-2*log(para[j][i]))*cos(2*M_PI*para[j+1][i]);
          data[i][j+1]=sqrt(-2*log(para[j][i]))*sin(2*M_PI*para[j+1][i]);
        }
        else
          data[i][j  ]=sqrt(-2*log(para[j][i]))*cos(2*M_PI*para[j+1][i]);

        // now normalize these Gaussian variables
        norm=0;
        for(j=0;j<dim;j++)
         norm+=data[i][j]*data[i][j];
        norm=sqrt(norm);
        for(j=0;j<dim;j++)
         data[i][j]/=norm;
      }

      for(i=0;i<dim-1;i++)
       delete para[i];
      delete para;
      break;

    case 2: // generatesa 3d affine space in 5 dimensions
      dim=5;
      dim_right=3;

      data = new float*[num];
      for(i=0;i<num;i++)
       data[i]=new float[dim];

      para1=new float[num];
      para2=new float[num];
      para3=new float[num];

      for(i=0;i<num;i++)
      {
        para1[i]=4*uni();
        para2[i]=4*uni();
        para3[i]=4*uni();
      }

      // the affine space
      for(i=0;i<num;i++)
      {
        data[i][0]= 1.2*para1[i]- 0.5*para2[i]+3;
        data[i][1]= 0.5*para1[i]+ 0.9*para3[i]-1;
        data[i][2]=-0.5*para1[i]- 0.2*para2[i] +   para3[i];
        data[i][3]= 0.4*para1[i]- 0.9*para2[i] -  0.1*para3[i];
        data[i][4]= 1.1*para1[i]- 0.3*para3[i]+8;
      }

      delete para1,para2,para3;
      break;

    case 3: // generates a strange 4-dim figure in 6 dimensions
      dim=6;
      dim_right=4;

      data = new float*[num];
      for(i=0;i<num;i++)
       data[i]=new float[dim];

      para=new float*[4];
      for(i=0;i<4;i++)
        para[i]=new float[num];

      for(i=0;i<4;i++)
       for(j=0;j<num;j++)
        para[i][j]=uni();

      //generate the figure
      for(i=0;i<num;i++)
      {
        data[i][0]= para[1][i]*para[1][i]*cos(2*M_PI*para[0][i]);
        data[i][1]= para[2][i]*para[2][i]*sin(2*M_PI*para[0][i]);
        data[i][2]= para[1][i]+para[2][i]+pow(para[1][i]-para[3][i],2);
        data[i][3]= para[1][i]-2*para[2][i]+pow(para[0][i]-para[3][i],2);
        data[i][4]=-para[1][i]-2*para[2][i]+pow(para[2][i]-para[3][i],2);
        data[i][5]= para[0][i]*para[0][i]-para[1][i]*para[1][i]
                   +para[2][i]*para[2][i]-para[3][i]*para[3][i];
      }

      for(i=0;i<4;i++)
       delete para[i];
      delete para;

      break;

     case 4: // generates a 4-dim manifold in 8 dimensions
       dim=8;
       dim_right=4;

       data = new float*[num];
       for(i=0;i<num;i++)
         data[i]=new float[dim];

       para=new float*[4];
       for(i=0;i<4;i++)
        para[i]=new float[num];

       for(i=0;i<4;i++)
        for(j=0;j<num;j++)
         para[i][j]=uni();

       //generate the figure
       for(i=0;i<num;i++)
       {
         data[i][0]= para[1][i]*cos(2*M_PI*para[0][i]);
         data[i][1]= para[1][i]*sin(2*M_PI*para[0][i]);
         data[i][2]= para[2][i]*cos(2*M_PI*para[1][i]);
         data[i][3]= para[2][i]*sin(2*M_PI*para[1][i]);
         data[i][4]= para[3][i]*cos(2*M_PI*para[2][i]);
         data[i][5]= para[3][i]*sin(2*M_PI*para[2][i]);
         data[i][6]= para[0][i]*cos(2*M_PI*para[3][i]);
         data[i][7]= para[0][i]*sin(2*M_PI*para[3][i]);
       }

       for(i=0;i<4;i++)
        delete para[i];
       delete para;
       break;

     case 5: // generates a 2d helix in 3d
      dim=3;
      dim_right=2;

      data = new float*[num];
      for(i=0;i<num;i++)
       data[i]=new float[dim];

      radii=new float[num];
      phi  =new float[num];

      for(i=0;i<num;i++)
      {
        radii[i]=2*M_PI*uni()-M_PI;
        phi[i]  =(uni())*2*M_PI;
      }

      // the helix
      for(i=0;i<num;i++)
      {
        data[i][0]=radii[i]*sin(phi[i]);
        data[i][1]=radii[i]*cos(phi[i]);
        data[i][2]=phi[i];
      }
      delete radii;
      delete phi;
      break;

    case 6: // generates a 6-dim manifold in 36- dimensions
       dim=36;
       dim_right=6;

       data = new float*[num];
       for(i=0;i<num;i++)
         data[i]=new float[dim];

       para=new float*[6];
       for(i=0;i<6;i++)
        para[i]=new float[num];

       for(i=0;i<6;i++)
        for(j=0;j<num;j++)
         para[i][j]=uni();

       //generate the figure
       for(i=0;i<num;i++)
       {
         for(j=0;j<10;j+=2)
         {
           data[i][j  ]= para[j/2+1][i]*cos(2*M_PI*para[j/2][i]);
           data[i][j+1]= para[j/2+1][i]*sin(2*M_PI*para[j/2][i]);
         }
         data[i][10]= para[0][i]*cos(2*M_PI*para[5][i]);
         data[i][11]= para[0][i]*sin(2*M_PI*para[5][i]);
       }
       for(i=0;i<num;i++)
       {
         for(j=0;j<12;j++)
         {
           data[i][j+12]=data[i][j];
           data[i][j+24]=data[i][j];
         }
       }

       for(i=0;i<6;i++)
        delete para[i];
       delete para;
       break;

    case 7: //swiss roll
      dim=3;
      dim_right=2;

      data = new float*[num];
      for(i=0;i<num;i++)
       data[i]=new float[dim];

      radii=new float[num];
      phi  =new float[num];

      for(i=0;i<num;i++)
      {
        radii[i]=2*M_PI*uni()-M_PI;
        phi[i]  =(uni())*2*M_PI;
      }

      // the swiss roll
      for(i=0;i<num;i++)
      {
        data[i][0]=(phi[i])*sin(2.5*phi[i]);
        data[i][1]=radii[i];
        data[i][2]=(phi[i])*cos(2.5*phi[i]);

      }
      delete radii;
      delete phi;
      break;


    case 8: // generates a 12-dim  in 72- dimensions
       dim=72;
       dim_right=12;

       data = new float*[num];
       for(i=0;i<num;i++)
         data[i]=new float[dim];

       para4=new double*[12];
       for(i=0;i<12;i++)
        para4[i]=new double[num];

       for(i=0;i<12;i++)
        for(j=0;j<num;j++)
         para4[i][j]=uni();

       //generate the figure
       for(i=0;i<num;i++)
       {
         for(j=0;j<22;j+=2)
         {
           data[i][j  ]= para4[j/2+1][i]*cos(2*M_PI*para4[j/2][i]);
           data[i][j+1]= para4[j/2+1][i]*sin(2*M_PI*para4[j/2][i]);
         }
         data[i][22]= para4[0][i]*cos(2*M_PI*para4[11][i]);
         data[i][23]= para4[0][i]*sin(2*M_PI*para4[11][i]);
       }
       for(i=0;i<num;i++)
       {
         for(j=0;j<24;j++)
         {
           data[i][j+24]=data[i][j];
           data[i][j+48]=data[i][j];
         }
       }

       for(i=0;i<12;i++)
        delete para4[i];
       delete para4;
       break;

     case 9: // generates a 20-dim affine space  in 20- dimensions
       dim=20;
       dim_right=20;

       data = new float*[num];
       for(i=0;i<num;i++)
         data[i]=new float[dim];

       para=new float*[20];
       for(i=0;i<20;i++)
        para[i]=new float[num];

       for(i=0;i<20;i++)
        for(j=0;j<num;j++)
         para[i][j]=5.0*uni()-2.5;

       vectors=new float*[20];
       for(i=0;i<20;i++)
        vectors[i]=new float[dim];

       for(i=0;i<20;i++)
        for(j=0;j<dim;j++)
        {
          if(i==j)
           vectors[i][j]=1;
          else
           vectors[i][j]=0;//0.25*(uni()-0.125);
        }
         //vectors[i][j]=10*(uni()-0.5);

       //generate the figure
       for(i=0;i<num;i++)
       {
         for(j=0;j<dim;j++)
         {
           data[i][j]=0;
           for(k=0;k<20;k++)
            data[i][j]+=para[k][i]*vectors[k][j];
         }
       }

       for(i=0;i<20;i++)
        delete para[i];
       delete para;
       break;



    case 10: // generates a k-hypercube
      nr_par=dim-1;
      dim_right=dim-1;

      data = new float*[num];
      for(i=0;i<num;i++)
       data[i]=new float[dim];

      para=new float*[nr_par];
      for(i=0;i<nr_par;i++)
       para[i]=new float[num];

      for(j=0;j<nr_par;j++)
       for(i=0;i<num;i++)
       {
         para[j][i]=uni();
         while(para[j][i]==0)
          para[j][i]=uni();
       }

      for(i=0;i<num;i++)
      {
        for(j=0;j<dim-1;j++)
        {
          data[i][j]=para[j][i];
        }
        data[i][dim-1]=0;
      }

      for(j=0;j<dim-1;j++)
       for(i=0;i<num;i++)
        para[j][i]  =(uni());


      for(i=0;i<dim-1;i++)
       delete para[i];
      delete para;
      break;

    case 11: //generates a moebius-band in 3d
      dim=3;
      dim_right=2;

      data = new float*[num];
      for(i=0;i<num;i++)
       data[i]=new float[dim];

      para = new float*[2];
      for(i=0;i<2;i++)
       para[i]=new float[num];

      for(i=0;i<num;i++)
      {
        para[0][i]=(uni())*2*M_PI;
        para[1][i]=2.0*(uni())-1.0;
      }

      // the moebius band
      for(i=0;i<num;i++)
      {
        data[i][0]=(1+0.5*para[1][i]*cos(0.5*10.0*para[0][i]))*cos(para[0][i]);
        data[i][1]=(1+0.5*para[1][i]*cos(0.5*10.0*para[0][i]))*sin(para[0][i]);
        data[i][2]=0.5*para[1][i]*sin(0.5*10.0*para[0][i]);
      }
      for(i=0;i<2;i++)
       delete para[i];
      delete para;
      break;


    case 12: // generate a multivariate Gaussian
      data = new float*[num];
      for(i=0;i<num;i++)
       data[i]=new float[dim];

      para2=new float[dim];
      para2[0]=1.0;
      for(i=1;i<dim;i++)
       para2[i]=1.0;

      nr_par=dim;
      if(nr_par %2 ==1) nr_par+=1;

      para=new float*[nr_par];
      for(i=0;i<nr_par;i++)
       para[i]=new float[num];

      for(j=0;j<nr_par;j++)
       for(i=0;i<num;i++)
       {
         para[j][i]=uni();
         while(para[j][i]==0)
          para[j][i]=uni();
       }

      for(i=0;i<num;i++)
      {
        for(j=0;j<dim-2;j+=2)
        {
          data[i][j  ]=para2[j]  *sqrt(-2*log(para[j][i]))*cos(2*M_PI*para[j+1][i]);
          data[i][j+1]=para2[j+1]*sqrt(-2*log(para[j][i]))*sin(2*M_PI*para[j+1][i]);
        }
        if(dim % 2==0)
        {
          data[i][j  ]=para2[j]  *sqrt(-2*log(para[j][i]))*cos(2*M_PI*para[j+1][i]);
          data[i][j+1]=para2[j+1]*sqrt(-2*log(para[j][i]))*sin(2*M_PI*para[j+1][i]);
        }
        else
          data[i][j  ]=para2[j]*sqrt(-2*log(para[j][i]))*cos(2*M_PI*para[j+1][i]);

      }
      for(i=0;i<nr_par;i++)
       delete para[i];
      delete para;
      delete para2;
      break;

    case 13:  //generates a one-dimensional curve in dim dimensions
      data = new float*[num];
      for(i=0;i<num;i++)
       data[i]=new float[dim];

      // parameter of the curve
      phi  =new float[num];
      for(i=0;i<num;i++)
       phi[i]  =uni()*2.0*M_PI;  // random numbers in [0,2pi]

      // the curve
      for(i=0;i<num;i++)
      {
        for(j=0;j<dim;j++)
        {
          data[i][j]=phi[i]/(2.0*M_PI);
          for(k=0;k<j;k++)
           data[i][j]+=sin((k+1)*phi[i]);
          data[i][j]/=(1.0*(j+1));
        }
      }
      delete phi;
      break;

  }
  return data;
}



float comp_slope(float* x,float* y,long num)
{
  double avg_x=0, avg_y=0, slope=0,enu=0,deno=0, w_tot=0;
  long i;
  double* w=new double[num];

  for(i=0;i<num;i++)
   { w[i]=1.0/(1.0*(num-i)); w_tot+=w[i]; }
  for(i=0;i<num;i++)
   { avg_x+=w[i]*x[i]; avg_y+=w[i]*y[i];}

  avg_x/=w_tot; avg_y/=w_tot;

  for(i=0;i<num;i++)
  {
    enu  +=w[i]*x[i]*(y[i]-avg_y);
    deno +=w[i]*x[i]*(x[i]-avg_x);
  }
  slope=enu/deno;
  delete w;
  return slope;
}



// estimates the intrinsic dimension of the dim-dimensional data from num samples
// max_dim: is the maximal possible dimension to be estimated
// if TAKENS==1, CORR==1 then also the Takens estimator and the correlation dimension estimator are calculated
// SPARSE=1 if SPARSE data is used

void Calc_Dimension(float** data, long dim, long max_dim, long num, char TAKENS,char CORR,char SPARSE)
{
  long i,j,k,l,m;

  float dist;
  float* min_dist=new float[num];
  for(i=0;i<num;i++)
   min_dist[i]=1E10;

  double* norm2;
  if(SPARSE)
  {
    norm2 = new double[num];
    for(i=0;i<num;i++)
    {
      norm2[i]=0;
      for(k=0;k<nr_elements[i];k++) norm2[i]+=data[i][sp_index[i][k]]*data[i][sp_index[i][k]];
    }
  }

  /*time_t anf; time_t end;
  anf=time(NULL);
  struct time t1,t2;
  gettime(&t1);  */

  // compute the one-nearest neighbor distances of each point
  for(i=0;i<num;i++)
  {
    for(j=i+1;j<num; j++)
    {
      // compute the squared distance to this point
      dist=0;
      if(!SPARSE)
      {
        for(k=0;k<dim;k++)
        {
          dist += (data[i][k]-data[j][k])*(data[i][k]-data[j][k]);
          if(dist>min_dist[i] && dist>min_dist[j])
           { k=dim; }
        }
      }
      else
      {
        dist = norm2[i]+norm2[j];
        if(nr_elements[i]<nr_elements[j])
        {
          for(k=0;k<nr_elements[i];k++)
           if(data[j][sp_index[i][k]]!=0) dist-=2.0*data[i][sp_index[i][k]]*data[j][sp_index[i][k]];
        }
        else
        {
          for(k=0;k<nr_elements[j];k++)
           if(data[i][sp_index[j][k]]!=0) dist-=2.0*data[i][sp_index[j][k]]*data[j][sp_index[j][k]];
        }
        if(dist<0) dist=0;
      }
      if(dist < min_dist[i])  // check if it is smaller than the maximal minimal distance
         min_dist[i]=dist;
      if(dist < min_dist[j])
         min_dist[j]=dist;
    }
  }

  // compute mean minimal distance
  float avg_min=0;
  float min_min=1E+10;
  float max_min=0;
  for(i=0;i<num;i++)
  {
    avg_min+=sqrt(min_dist[i]);
    if(min_dist[i]<min_min)
     min_min=min_dist[i];
    if(min_dist[i]>max_min)
     max_min=min_dist[i];
  }
  max_min=sqrt(max_min);
  min_min=sqrt(min_min);
  avg_min/=num;

  // compute variance of the minimal distance
  float std_min=0;
  for(i=0;i<num;i++)
   std_min+=pow(sqrt(min_dist[i])-avg_min,2);
  std_min/=num;
  std_min=sqrt(std_min);

  long NR_EST=5;

  float TAKENS_SCALE=0;                 // the scale (upper bound on the distances) used in the Takens estimator
  float* CORR_SCALE=new float[NR_EST];

  TAKENS_SCALE=max_min*max_min;
  for(k=0;k<5;k++)
  {
    CORR_SCALE[k] =avg_min + k*0.2*std_min;
    CORR_SCALE[k]*=CORR_SCALE[k];
  }
  TAKENS_SCALE=avg_min+std_min;
  TAKENS_SCALE*=TAKENS_SCALE;

  // cut_dist is defined as the mean minimal distance
  float cut_dist=avg_min;  //avg_min+0*std_min;
  float cut_dist_sq=pow(cut_dist,2);


  // now the final algorithm starts
  // we assume a form h ~ (log n/n)^1/d, this guarantees nh^d -> infinity and h -> 0

  float scale=0;
  double** h=new double*[max_dim];
  for(i=0;i<max_dim;i++)
   h[i]=new double[NR_EST];

  double** h_squared=new double*[max_dim];
  for(i=0;i<max_dim;i++)
   h_squared[i]=new double[NR_EST];

  double*** est=new double**[max_dim];
  for(i=0;i<max_dim;i++)
  {
    est[i]=new double*[NR_EST];
    for(j=0;j<NR_EST;j++)
     est[i][j]=new double[NR_EST*NR_EST];
  }

  float TAKENS_EST=0;
  long TAKENS_NR=0;
  float* CORR_EST=new float[NR_EST];
  for(i=0;i<NR_EST;i++) CORR_EST[i]=0;

  long* divisions=new long[NR_EST];
  for(i=0;i<NR_EST;i++)
   divisions[i]=NR_EST-i;

  float* size=new float[NR_EST];
  for(i=0;i<NR_EST;i++)
   size[i]=floor(num/(1.0*divisions[i]));


  // initialize the scale of each dimension and the corresponding h
  for(m=1;m<=max_dim;m++)
  {
    for(k=0;k<NR_EST;k++)                      // we have NR_EST partitions of the data
    {
      h[m-1][k]=cut_dist*powl((num*log(size[k]))/(size[k]*log(num)),1.0/m);
      h_squared[m-1][k]=h[m-1][k]*h[m-1][k];
      for(j=0;j<NR_EST*NR_EST;j++)
       est[m-1][k][j]=0;
    }
  }

  float min_scale=h[0][0]*h[0][0];
  long index_k_x; long index_k_y;

  for(i=0;i<num;i++)
  {
    for(j=i+1;j<num;j++)
    {
      dist=0;
      if(!SPARSE)
      {
        for(l=0;l<dim;l++)
        {
          dist += (data[i][l]-data[j][l])*(data[i][l]-data[j][l]);
          if(dist>min_scale && (!TAKENS || dist > TAKENS_SCALE) && (!CORR || dist > CORR_SCALE[4]))
           l=dim;
        }
      }
      else
      {
        dist = norm2[i]+norm2[j];
        if(nr_elements[i]<nr_elements[j])
        {
          for(k=0;k<nr_elements[i];k++)
            if(data[j][sp_index[i][k]]!=0) dist-=2.0*data[i][sp_index[i][k]]*data[j][sp_index[i][k]];
        }
        else
        {
          for(k=0;k<nr_elements[j];k++)
            if(data[i][sp_index[j][k]]!=0) dist-=2.0*data[i][sp_index[j][k]]*data[j][sp_index[j][k]];
        }
        if(dist<0) dist=0;
      }
      if(TAKENS && dist<TAKENS_SCALE && dist>0 )
        { TAKENS_EST+=0.5*log(dist/TAKENS_SCALE);  TAKENS_NR++; }
      if(dist<min_scale)
      {
        if(CORR)
        {
          for(k=0;k<5;k++)
          {
            if(dist<CORR_SCALE[k]) CORR_EST[k]++;
          }
        }
        for(k=0;k<NR_EST;k++)
        {
          index_k_x=i % divisions[k];
          index_k_y=j % divisions[k];
          for(m=1;m<=max_dim;m++)  // try all dimensions from 1 to max_dim
          {
            if(dist < h_squared[m-1][k])
            {
              if(index_k_x<index_k_y)
                est[m-1][k][index_k_x*divisions[k]+index_k_y]+=1.0-dist/h_squared[m-1][k];
              else
                est[m-1][k][index_k_y*divisions[k]+index_k_x]+=1.0-dist/h_squared[m-1][k];
            }
            else m=max_dim;  // note that h_squared is monotonically decreasing with m
                             // which means that if the condition dist < h_squared[m-1][k]  does not
                             // hold for m it will also not hold for any larger m
          }
        }
      }
    }
  }
  long double total=0;
  for(m=1;m<=max_dim;m++)
  {
    for(k=0;k<NR_EST;k++)
    {
      total=0;
      for(i=0;i<divisions[k];i++)
      {
        for(j=0;j<divisions[k];j++)
        {
          if(i==j)
          {
            if(i < num % divisions[k])
             total+=est[m-1][k][i*divisions[k]+j]/(0.5*(size[k]+1)*size[k]);//*pow(h[m-1][k],m));
            else
             total+=est[m-1][k][i*divisions[k]+j]/(0.5*size[k]*(size[k]-1));//*pow(h[m-1][k],m));
          }
          else
          {
            if(i < num % divisions[k] && j < num % divisions[k])
             total+=est[m-1][k][i*divisions[k]+j]/2.0/(0.5*(size[k]+1)*(size[k]+1));//*pow(h[m-1][k],m));
            else
            {
               if(i>=num % divisions[k] && j >= num % divisions[k])
                total+=est[m-1][k][i*divisions[k]+j]/2.0/(0.5*size[k]*size[k]);//*pow(h[m-1][k],m));
               else
                total+=est[m-1][k][i*divisions[k]+j]/2.0/(0.5*(size[k]+1)*size[k]);//*pow(h[m-1][k],m));
            }
          }
          /*if(i==j)
           total+=est[m-1][k][i*divisions[k]+j]/(0.5*size[k]*(size[k]-1)*pow(h[m-1][k],m));
          else
           total+=est[m-1][k][i*divisions[k]+j]/2.0/(0.5*size[k]*size[k]*pow(h[m-1][k],m)); */
        }
      }
      est[m-1][k][0]=2.0*total/(divisions[k]*(divisions[k]+1));
    }
  }

  if(TAKENS)
  {
    TAKENS_EST/=TAKENS_NR;//size[4]*(size[4]-1);
    TAKENS_EST=-1/TAKENS_EST;
    TAK_DIM=TAKENS_EST;
    //printf("Takens estimate: %f\n",TAKENS_EST);
  }
  if(CORR)
  {
    for(k=0;k<5;k++)
    {
      CORR_EST[k]/=0.5*size[NR_EST-1]*(size[NR_EST-1]-1);
      CORR_EST[k]=log(CORR_EST[k]);
      CORR_SCALE[k]=log(sqrt(CORR_SCALE[k]));
    }
    CORR_DIM=comp_slope(CORR_SCALE,CORR_EST,5);
    //printf("Correlation dimension: %f\n",CORR_DIM);
  }
  float minimum=100000;
  float* slope_est=new float[max_dim];
  for(m=1;m<=max_dim;m++)
  {
    for(k=0;k<NR_EST;k++)
    {
      CORR_EST[k]=log(est[m-1][k][0]) - m*log(h[m-1][k]);
      CORR_SCALE[k]=log(h[m-1][k]);
    }
    slope_est[m-1]=comp_slope(CORR_SCALE,CORR_EST,NR_EST);
    if(fabs(slope_est[m-1])<minimum) { minimum=slope_est[m-1]; DIM=m; }
    //printf("Slope of Dimension: %d, %f, %f\n",m,slope_est[m-1],slope_est[m-1]+m);
  }
  //printf("Intrinsic Dim.  estimate: %d\n",min_est);

  /*end=time(NULL);
  gettime(&t2);
  float est_time;

  if(t2.ti_hund - t1.ti_hund < 0) { t1.ti_sec=t1.ti_sec+1; est_time=(100.0+(t2.ti_hund - t1.ti_hund))/100.0;}
  else {est_time=(t2.ti_hund - t1.ti_hund)/100.0;}

  if(t2.ti_sec - t1.ti_sec < 0) { t1.ti_min=t1.ti_min+1; est_time+=60.0+(t2.ti_sec - t1.ti_sec);}
  else {est_time+=t2.ti_sec - t1.ti_sec;}

  if(t2.ti_min - t1.ti_min < 0) { est_time+=60.0*(60.0-(t2.ti_min-t1.ti_min)); }
  else {est_time +=60.0*(t2.ti_min-t1.ti_min);}
  printf("Minimal Scale: %f,   Time: %1.2f seconds \n",h[0][4],est_time); */
  printf("Minimal Scale: %f \n",h[0][4]);


  // free memory
  delete min_dist;
}

int main(int argc, char * argv[])
{
  // arguments for the intrinsic dimensionality estimator
  // Two options:

  // First option - Data is loaded from a file
  // argument 1: name of the file with the data
  //             where the file should have the following format:
  //             Binary, float, and data points are written sequentially
  //             that is p[0][0],...,p[0][dim-1],p[1][0],....
  // argument 2: datatype  - c, uc, l , ul, f, d
  //                       - char, unsigned char, long, unsigned long, float, double
  // argument 3: dimension of the data
  // argument 4: upper bound on the dimension to be estimated (<=dim)
  //             can usually limited to something below ~100 (default)

  // Second option - Data is generated by the program
  // argument 1: number of the dataset
  // argument 2: dimension of the dataset
  // argument 3: number of points
  // argument 4: number of runs

  FILE* trainfile;
  float mean_dim;
  char TAKENS=1; // 1: calculate also the Takens estimate
  char CORR  =1; // 1: calculate the correlation dimension
  long file_length,i,j;

  long max_dim, data_set, runs, data_type;

  // the data vector
  float** data;

  long num_read;

  float* dim_est, *dim_est_corr, *dim_est_tak;

  char SPARSE=0, OPTION;
  char* image_c; unsigned char* image_uc; long* image_l; unsigned long* image_ul; float* image_f; double* image_d;

  if(argc!=5)
  {
     char answer;
     printf("Wrong number of Arguments - Usage either :\n 1) DIM_EST_PUB Filename Dimension_of_the_Data Data_Type Maximal_Estimated_Dimension\n 2) DIM_EST_PUB Number_of_Dataset Dimension_of_Dataset Number_of_Dataset Number_of_Runs\n Press a Key");
     scanf("%s\n",&answer);
     exit(EXIT_FAILURE);
  }

  if(argc==5 && strlen(argv[1])>3) // First option - data is loaded from a file
  {
    OPTION=0; runs=1;
    // first get the filesize
    trainfile=fopen(argv[1],"rb");
    if(trainfile==NULL)
    {
      char answer;
      printf("File could not be opened - usage: DIM_EST_PUB Filename Dimension_of_the_Data Data_Type Maximal_Estimated_Dimension\n Press a Key");
      scanf("%s",&answer);
      exit(EXIT_FAILURE);
    }

    // third determine dimension of the data
    data_type=-1;
    if(strcmp(argv[2],"c")==0) data_type=0;
    if(strcmp(argv[2],"uc")==0) data_type=1;
    if(strcmp(argv[2],"l")==0) data_type=2;
    if(strcmp(argv[2],"ul")==0) data_type=3;
    if(strcmp(argv[2],"f")==0) data_type=4;
    if(strcmp(argv[2],"d")==0) data_type=5;

    if(data_type==-1)
    {
      char answer;
      printf("Wrong datatype - should be c, uc, l, ul, f or d - usage: DIM_EST_PUB Filename Dimension_of_the_Data Maximal_Estimated_Dimension\n Press a Key");
      scanf("%s",&answer);
      exit(EXIT_FAILURE);
    }

    // third determine dimension of the data
    dim=atol(argv[3]);
    if(dim<=0)
    {
      char answer;
      printf("Dimension_of_the_Data could not be converted to a number - usage: DIM_EST_PUB Filename Dimension_of_the_Data Maximal_Estimated_Dimension\n Press a Key");
      scanf("%s",&answer);
      exit(EXIT_FAILURE);
    }

    // fourth determine maximal number of dimension
    max_dim=atol(argv[4]);
    if(max_dim<=0)
    {
      char answer;
      printf("Maximal_Estimated_Dimension could not be converted to a number - usage: DIM_EST_PUB Filename Dimension_of_the_Data Maximal_Estimated_Dimension\n Press a Key");
      scanf("%s",&answer);
      exit(EXIT_FAILURE);
    }

    // now determine the number of data points
    fseek(trainfile,0,SEEK_END);
    file_length=ftell(trainfile);

    switch(data_type)
    {
      case 0: image_c=new char[dim]; num=file_length/dim; break;
      case 1: image_uc=new unsigned char[dim]; num=file_length/dim; break;
      case 2: image_l=new long[dim];  num=file_length/(dim*4); break;
      case 3: image_ul=new unsigned long[dim]; num=file_length/(dim*4); break;
      case 4: image_f=new float[dim]; num=file_length/(dim*4); break;
      case 5: image_d=new double[dim]; num=file_length/(dim*8); break;
    }

    // go back to the beginning of the file
    num_read=fseek(trainfile,0,SEEK_SET);

    // allocate memory for the data
    data = new float*[num];
    num_read=0;
    for(i=0;i<num;i++)
    {
       data[i]=new float[dim];
       switch(data_type)
       {
         case 0: num_read+=fread(image_c,sizeof(char),dim,trainfile);  for(j=0;j<dim;j++) data[i][j]=image_c[j]; break; // read in the data from the file
         case 1: num_read+=fread(image_uc,sizeof(char),dim,trainfile); for(j=0;j<dim;j++) data[i][j]=image_uc[j]; break;
         case 2: num_read+=fread(image_l,sizeof(long),dim,trainfile);  for(j=0;j<dim;j++) data[i][j]=image_l[j];break;
         case 3: num_read+=fread(image_ul,sizeof(long),dim,trainfile); for(j=0;j<dim;j++) data[i][j]=image_ul[j]; break;
         case 4: num_read+=fread(image_f,sizeof(float),dim,trainfile); for(j=0;j<dim;j++) data[i][j]=image_f[j]; break;
         case 5: num_read+=fread(image_d,sizeof(double),dim,trainfile); for(j=0;j<dim;j++) data[i][j]=image_d[j]; break;
       }
    }
    if(num<10)
    {
      char answer;
      printf("Number of data points less than 10\n Press a Key");
      scanf("%s",&answer);
      exit(EXIT_FAILURE);
    }
    fclose(trainfile);
    printf("Read %d points from the file, should read %d\n",num_read/dim,num);
  }

  if(argc==5 && strlen(argv[1])<=3) // Second option - data is generated by the program
  {
     OPTION=1;
     // first determine the number of the dataset
    data_set=atol(argv[1]);
    // second determine the number of dimensions
    dim=atol(argv[2]);
    if(dim<=0)
    {
      char answer;
      printf("Dimension_of_the_Data could not be converted to a number - usage: DIM_EST_PUB Number_of_Dataset Dimension_of_Dataset Number_of_Dataset Number_of_Runs\n Press a Key");
      scanf("%s",&answer);
      exit(EXIT_FAILURE);
    }
    // third determine the number of datapoints
    num=atol(argv[3]);
    if(num<10)
    {
      char answer;
      printf("Number_of_Points could not be converted to a number or number of data points less than 10 - usage: DIM_EST_PUB Number_of_Dataset Dimension_of_Dataset Number_of_Dataset Number_of_Runs \n Press a Key");
      scanf("%s",&answer);
      exit(EXIT_FAILURE);
    }
    // third determine number of runs
    runs=atol(argv[4]);
    if(runs<=0)
    {
      char answer;
      printf("Number_of_Runs could not be converted to a number - usage: DIM_EST_PUB Number_of_Dataset Dimension_of_Dataset Number_of_Dataset Number_of_Runs\n Press a Key");
      scanf("%s",&answer);
      exit(EXIT_FAILURE);
    }
  }

  dim_est=new float[runs];
  dim_est_corr=new float[runs];
  dim_est_tak=new float[runs];

  long right=0,corr_right=0,tak_right=0;
  float CORR_DIM_MEAN=0,CORR_DIM_MEAN2=0;
  float TAK_DIM_MEAN=0, TAK_DIM_MEAN2=0;
  for(ctr=0;ctr<runs;ctr++)
  {
     // generate data
     if(OPTION==1)
      { data=generate_data(data_set);  max_dim=dim;}

     if(max_dim>dim)
      max_dim=dim;
     if(max_dim>500) max_dim=500;

     Calc_Dimension(data, dim, max_dim, num,TAKENS,CORR,SPARSE);

     dim_est[ctr]=DIM; // Intrinsic Dim.  Estimate

     if(OPTION==1)
     {
       if(DIM==dim_right) right++;
       if(fabs(CORR_DIM-dim_right)<0.5) corr_right++;
       if(fabs(TAK_DIM-dim_right)<0.5)  tak_right++;
     }
     dim_est_corr[ctr]=ceil(CORR_DIM-0.5);
     dim_est_tak[ctr]=ceil(TAK_DIM-0.5);

     CORR_DIM_MEAN+=CORR_DIM;
     TAK_DIM_MEAN+=TAK_DIM;

     printf("\n");
     printf("Direct Estimates:\n");
     printf("-----------------\n");
     printf("Intrinsic Dim. Estimate: %1.0f, Correlation Dimension: %1.4f, Takens Estimator: %1.4f\n",dim_est[ctr],CORR_DIM,TAK_DIM);
     printf("\n");
     printf("Rounded Estimates:\n");
     printf("-----------------\n");
     printf("Intrinsic Dim. Estimate: %1.0f, Correlation Dimension: %1.0f, Takens Estimator: %1.0f\n",dim_est[ctr],dim_est_corr[ctr],dim_est_tak[ctr]);
     printf("\n");

     if(ctr<runs-1)           // delete data in every run
     {
       for(j=0;j<num;j++)
        delete data[j];
       delete data;
     }

  }
  if(OPTION==1)
  {
    CORR_DIM_MEAN/=runs; TAK_DIM_MEAN/=runs;
    float mean_est=0;
    float var_est=0,CORR_VAR_EST, TAK_VAR_EST;

    for(i=0;i<runs;i++)
     { mean_est+=dim_est[i]; CORR_DIM_MEAN2+=dim_est_corr[i]; TAK_DIM_MEAN2+=dim_est_tak[i];}
    mean_est/=runs; CORR_DIM_MEAN2/=runs;    TAK_DIM_MEAN2/=runs;

    for(i=0;i<runs;i++)
    {
      var_est+=(dim_est[i]-mean_est)*(dim_est[i]-mean_est);
      CORR_VAR_EST+=(dim_est_corr[i]-CORR_DIM_MEAN2)*(dim_est_corr[i]-CORR_DIM_MEAN2);
      TAK_VAR_EST+=(dim_est_tak[i]-TAK_DIM_MEAN2)*(dim_est_tak[i]-TAK_DIM_MEAN2);
    }
    var_est/=runs; CORR_VAR_EST/=runs; TAK_VAR_EST/=runs;
    var_est=sqrt(var_est); CORR_VAR_EST=sqrt(CORR_VAR_EST); TAK_VAR_EST=sqrt(TAK_VAR_EST);
    printf("Dataset: %d, Dim: %d, Number of Points: %d\n",data_set,dim,num);
    printf("Correct Dimension of this Dataset: %d\n", dim_right);
    printf("\n");
    printf("Correct Intrinsic dimension estimates: %d out of %d\n",right,runs);
    printf("        Mean Est: %f, Var_Est: %f, \n", mean_est, var_est);
    if(CORR)
    {
      printf("Correct Correlation dimension estimates: %d out of %d\n",corr_right,runs);
      printf("        Mean Est Dir: %1.4f, Mean Rounded Est: %1.4f, Var: %1.4f\n",CORR_DIM_MEAN,CORR_DIM_MEAN2,CORR_VAR_EST);
    }
    if(TAKENS)
    {
      printf("Correct Takens dimension estimates: %d out of %d\n",tak_right,runs);
      printf("        Mean Est Dir: %1.4f, Mean Rounded Est: %1.4f, Var: %1.4f\n",TAK_DIM_MEAN,TAK_DIM_MEAN2,TAK_VAR_EST);
    }
  }

  for(j=0;j<num;j++)
    delete data[j];
  delete data;

  if(SPARSE)
  {
    for(i=0;i<num;i++) delete sp_index[i];
    delete sp_index;
    delete nr_elements;
  }
}


