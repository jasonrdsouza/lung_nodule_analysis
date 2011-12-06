/****************************************************************/
/* Example VisX4 program tcinterp                               */
/*  Compute nearest neighbor interpolation of the input 3d image*/
/* Syntax:                                                      */
/*        interp if=infile of=outfile dx= dy= dz= nx= ny= nz=   */
/****************************************************************/
/****************************************************************/
/*  This Software is Copyright (c) 1996, 2000 Anthony P. Reeves */
/*  All rights reserved.                                        */
/****************************************************************/

// THIS IS LEKIEN'S VERSION.  IT WORKS BUT IS INEFFICIENT

#include "VisXV4.h"          /* VisX structure include file     */
#include "Vutil.h"           /* VisX utility header files       */
#include "coeff.h"
#include <math.h>

VisXfile_t *VXin,            /* input file structure            */
	       *VXout;       /* output file structure           */
VisXelem_t *VXlist;          /* VisX data structure             */
VisXparam_t par[] =          /* command line structure          */
{
    "if=",    0,   /* input file       */
    "of=",    0,   /* output file      */
    "-v",     0,   /* visible flag     */
    "dx=",    0,   /* old voxel size   */
    "dy=",    0,   /* old voxel size   */
    "dz=",    0,   /* old voxel size   */
    "nx=",    0,   /* new voxel size   */
    "ny=",    0,   /* new voxel size   */
    "nz=",    0,   /* new voxel size   */
     0,       0    /* list termination */
};

#define  IVAL   par[0].val
#define  OVAL   par[1].val
#define  VFLAG  par[2].val
#define  OLDX   par[3].val
#define  OLDY   par[4].val
#define  OLDZ   par[5].val
#define  NEWX   par[6].val
#define  NEWY   par[7].val
#define  NEWZ   par[8].val

int ijk2n(int i, int j, int k) {
  return(i+4*j+16*k);
}

void tricubic_get_coeff_stacked(double a[64], double x[64]) {
  int i,j;
  for (i=0;i<64;i++) {
    a[i]=(double)(0.0);
    for (j=0;j<64;j++) {
      a[i]+=A[i][j]*x[j];
    }
  }
}

void tricubic_get_coeff(double a[64], double f[8], double dfdx[8], double dfdy[8], double dfdz[8], double d2fdxdy[8], double d2fdxdz[8], double d2fdydz[8], double d3fdxdydz[8]) {
  int i;
  double x[64];
  for (i=0;i<8;i++) {
    x[0+i]=f[i];
    x[8+i]=dfdx[i];
    x[16+i]=dfdy[i];
    x[24+i]=dfdz[i];
    x[32+i]=d2fdxdy[i];
    x[40+i]=d2fdxdz[i];
    x[48+i]=d2fdydz[i];
    x[56+i]=d3fdxdydz[i];
  }
  tricubic_get_coeff_stacked(a,x);
}

double tricubic_eval(double a[64], double x, double y, double z) {
  int i,j,k;
  double ret=(double)(0.0);
  /* TRICUBIC EVAL
     This is the short version of tricubic_eval. It is used to compute
     the value of the function at a given point (x,y,z). To compute
     partial derivatives of f, use the full version with the extra args.
  */
  for (i=0;i<4;i++) {
    for (j=0;j<4;j++) {
      for (k=0;k<4;k++) {
	ret+=a[ijk2n(i,j,k)]*pow(x,i)*pow(y,j)*pow(z,k);
      }
    }
  }
  return(ret);
}

main(argc, argv)
int argc;
char *argv[];
{

VisX3dim_t orig;                   /* input v3dim structure        */
VisX3dim_t tc_interp;              /* output v3dim structure       */
VisXelem_t *vptr;
int	   i,j,k,q;                  /* index counters               */
unsigned char t;
double      dx;
double      dy;
double      dz;
double      nx;
double      ny;
double      nz;
int	   nnx0, nnx1, nnx2, nnx3, xt;
int	   nny0, nny1, nny2, nny3, yt;
int	   nnz0, nnz1, nnz2, nnz3, zt;
double     f[8];
double     dfdx[8];
double     dfdy[8];
double     dfdz[8];
double     d2fdxdy[8];
double     d2fdxdz[8];
double     d2fdydz[8];
double     d3fdxdydz[8];
double     a[64];
double     result;

    Vparse(&argc, &argv, par);     /* parse the command line       */
    VXin  = VXopen(IVAL, 0);       /* open input file              */
    VXout = VXopen(OVAL, 1);       /* open the output file         */
    VXlist = VXread(VXin);         /* read file                    */
    dx= atof(OLDX);
    dy= atof(OLDY);
    dz= atof(OLDZ);
    nx= atof(NEWX);
    ny= atof(NEWY);
    nz= atof(NEWZ);

    if(VXNIL == (vptr = VXfind(VXlist, VX_PBYTE))){
       fprintf(stderr, "nninterp: no byte images found\n");
       exit(1);
    }

    VXset3dim(&orig, vptr, VXin);    /* initialize input structure   */

    float bbx[6]= {ceil(orig.bbx[0]*dx/nx), floor(orig.bbx[1]*dx/nx), ceil(orig.bbx[2]*dy/ny), floor(orig.bbx[3]*dy/ny), ceil(orig.bbx[4]*dz/nz), floor(orig.bbx[5]*dz/nz)};
    VXmake3dim(&tc_interp, VX_PBYTE, bbx, 1);

    if(VFLAG){
       fprintf(stderr,"bbx is %f %f %f %f %f %f\n", orig.bbx[0],
           orig.bbx[1],orig.bbx[2],orig.bbx[3],orig.bbx[4],orig.bbx[5]);
    }

    /* simple pixel computation -- notice the order of the loops */
    for (k = tc_interp.zlo + 2*dz/nz; k <= tc_interp.zhi - 2*dz/nz; k++){
      for (j = tc_interp.ylo + 2*dy/ny; j <= tc_interp.yhi - 2*dy/ny; j++){
        for (i = tc_interp.xlo + 2*dx/nx; i <= tc_interp.xhi - 2*dx/nx; i++){
            /* remap the pixels from the new voxel space to the old voxel space */
            xt= i;
            yt= j;
            zt= k;

            nnx1= (xt*nx/dx < orig.xlo ? ceil(xt*nx/dx) : floor(xt*nx/dx));
            nny1= (yt*ny/dy < orig.ylo ? ceil(yt*ny/dy) : floor(yt*ny/dy));
            nnz1= (zt*nz/dz < orig.zlo ? ceil(zt*nz/dz) : floor(zt*nz/dz));

            nnx2= (xt*nx/dx > orig.xhi ? floor(xt*nx/dx) : ceil(xt*nx/dx));
            nny2= (yt*ny/dy > orig.yhi ? floor(yt*ny/dy) : ceil(yt*ny/dy));
            nnz2= (zt*nz/dz > orig.zhi ? floor(zt*nz/dz) : ceil(zt*nz/dz));
            
            nnx0= nnx1 - 1;
            nny0= nny1 - 1;
            nnz0= nnz1 - 1;
           
            nnx3= nnx2 + 1;
            nny3= nny2 + 1;
            nnz3= nnz2 + 1;

            // 8 NN pixels
            f[0]= orig.u[nnz1][nny1][nnx1]; //f[8] (see tricubic_get_coeff)
            f[1]= orig.u[nnz1][nny1][nnx2];
            f[2]= orig.u[nnz1][nny2][nnx1];
            f[3]= orig.u[nnz1][nny2][nnx2];
 
            f[4]= orig.u[nnz2][nny1][nnx1];
            f[5]= orig.u[nnz2][nny1][nnx2];
            f[6]= orig.u[nnz2][nny2][nnx1];
            f[7]= orig.u[nnz2][nny2][nnx2];

            /* compute all the partial derivatives */

            //dfdx[8]
            //dfdx(i,j,k)= f(i+1,j,k)-f(i-1,j,k) / 2
            dfdx[0]= (orig.u[nnz1][nny1][nnx2] - orig.u[nnz1][nny1][nnx0]) / 2.0;
            dfdx[1]= (orig.u[nnz1][nny1][nnx3] - orig.u[nnz1][nny1][nnx1]) / 2.0;         
            dfdx[2]= (orig.u[nnz1][nny2][nnx2] - orig.u[nnz1][nny2][nnx0]) / 2.0;
            dfdx[3]= (orig.u[nnz1][nny2][nnx3] - orig.u[nnz1][nny2][nnx1]) / 2.0;
                                                                         
            dfdx[4]= (orig.u[nnz2][nny1][nnx2] - orig.u[nnz2][nny1][nnx0]) / 2.0;
            dfdx[5]= (orig.u[nnz2][nny1][nnx3] - orig.u[nnz2][nny1][nnx1]) / 2.0;
            dfdx[6]= (orig.u[nnz2][nny2][nnx2] - orig.u[nnz2][nny2][nnx0]) / 2.0;
            dfdx[7]= (orig.u[nnz2][nny2][nnx3] - orig.u[nnz2][nny2][nnx1]) / 2.0;
            
            //dfdy[8]
            dfdy[0]= (orig.u[nnz1][nny2][nnx1] - orig.u[nnz1][nny0][nnx1]) / 2.0;
            dfdy[1]= (orig.u[nnz1][nny2][nnx2] - orig.u[nnz1][nny0][nnx2]) / 2.0;
            dfdy[2]= (orig.u[nnz1][nny3][nnx1] - orig.u[nnz1][nny1][nnx1]) / 2.0;
            dfdy[3]= (orig.u[nnz1][nny3][nnx2] - orig.u[nnz1][nny1][nnx2]) / 2.0;
                                                                         
            dfdy[4]= (orig.u[nnz2][nny2][nnx1] - orig.u[nnz2][nny0][nnx1]) / 2.0;
            dfdy[5]= (orig.u[nnz2][nny2][nnx2] - orig.u[nnz2][nny0][nnx2]) / 2.0;
            dfdy[6]= (orig.u[nnz2][nny3][nnx1] - orig.u[nnz2][nny1][nnx1]) / 2.0;
            dfdy[7]= (orig.u[nnz2][nny3][nnx2] - orig.u[nnz2][nny1][nnx2]) / 2.0;
                                                                         
            //dfdz[8]                                                    
            dfdz[0]= (orig.u[nnz2][nny1][nnx1] - orig.u[nnz0][nny1][nnx1]) / 2.0;
            dfdz[1]= (orig.u[nnz2][nny1][nnx2] - orig.u[nnz0][nny1][nnx2]) / 2.0;
            dfdz[2]= (orig.u[nnz2][nny2][nnx1] - orig.u[nnz0][nny2][nnx1]) / 2.0;
            dfdz[3]= (orig.u[nnz2][nny2][nnx2] - orig.u[nnz0][nny2][nnx2]) / 2.0;
                                                                         
            dfdz[4]= (orig.u[nnz3][nny1][nnx1] - orig.u[nnz1][nny1][nnx1]) / 2.0;
            dfdz[5]= (orig.u[nnz3][nny1][nnx2] - orig.u[nnz1][nny1][nnx2]) / 2.0;
            dfdz[6]= (orig.u[nnz3][nny2][nnx1] - orig.u[nnz1][nny2][nnx1]) / 2.0;
            dfdz[7]= (orig.u[nnz3][nny2][nnx2] - orig.u[nnz1][nny2][nnx2]) / 2.0;

            //d2fdxdy
            //d2fdxdy(i,j,k)=(f(i+1,j+1,k)-f(i+1,j,k)-f(i,j+1,k)+f(i-1,j-1,k))/4
            d2fdxdy[0]= (orig.u[nnz1][nny2][nnx2] - orig.u[nnz1][nny1][nnx2] - orig.u[nnz1][nny2][nnx1] + orig.u[nnz1][nny0][nnx0]) / 4.0; 
            d2fdxdy[1]= (orig.u[nnz1][nny2][nnx3] - orig.u[nnz1][nny1][nnx3] - orig.u[nnz1][nny2][nnx2] + orig.u[nnz1][nny0][nnx1]) / 4.0;
            d2fdxdy[2]= (orig.u[nnz1][nny3][nnx2] - orig.u[nnz1][nny2][nnx2] - orig.u[nnz1][nny3][nnx1] + orig.u[nnz1][nny1][nnx0]) / 4.0;
            d2fdxdy[3]= (orig.u[nnz1][nny3][nnx3] - orig.u[nnz1][nny2][nnx3] - orig.u[nnz1][nny3][nnx2] + orig.u[nnz1][nny1][nnx1]) / 4.0;
                                                                                                                                  
            d2fdxdy[4]= (orig.u[nnz2][nny2][nnx2] - orig.u[nnz2][nny1][nnx2] - orig.u[nnz2][nny2][nnx1] + orig.u[nnz2][nny0][nnx0]) / 4.0;
            d2fdxdy[5]= (orig.u[nnz2][nny2][nnx3] - orig.u[nnz2][nny1][nnx3] - orig.u[nnz2][nny2][nnx2] + orig.u[nnz2][nny0][nnx1]) / 4.0;
            d2fdxdy[6]= (orig.u[nnz2][nny3][nnx2] - orig.u[nnz2][nny2][nnx2] - orig.u[nnz2][nny3][nnx1] + orig.u[nnz2][nny1][nnx0]) / 4.0;
            d2fdxdy[7]= (orig.u[nnz2][nny3][nnx3] - orig.u[nnz2][nny2][nnx3] - orig.u[nnz2][nny3][nnx2] + orig.u[nnz2][nny1][nnx1]) / 4.0;

            //d2fdxdz
            //d2fdxdz(i,j,k)=(f(i+1,j,k+1)-f(i+1,j,k)-f(i,j,k+1)+f(i-1,j,k-1))/4
            d2fdxdz[0]= (orig.u[nnz2][nny1][nnx2] - orig.u[nnz1][nny1][nnx2] - orig.u[nnz2][nny1][nnx1] + orig.u[nnz0][nny1][nnx0]) / 4.0;
            d2fdxdz[1]= (orig.u[nnz2][nny1][nnx3] - orig.u[nnz1][nny1][nnx3] - orig.u[nnz2][nny1][nnx2] + orig.u[nnz0][nny1][nnx1]) / 4.0;
            d2fdxdz[2]= (orig.u[nnz2][nny2][nnx2] - orig.u[nnz1][nny2][nnx2] - orig.u[nnz2][nny2][nnx1] + orig.u[nnz0][nny2][nnx0]) / 4.0;
            d2fdxdz[3]= (orig.u[nnz2][nny2][nnx3] - orig.u[nnz1][nny2][nnx3] - orig.u[nnz2][nny2][nnx2] + orig.u[nnz0][nny2][nnx1]) / 4.0;
                                                                                                                                  
            d2fdxdz[4]= (orig.u[nnz3][nny1][nnx2] - orig.u[nnz2][nny1][nnx2] - orig.u[nnz3][nny1][nnx1] + orig.u[nnz1][nny1][nnx0]) / 4.0;
            d2fdxdz[5]= (orig.u[nnz3][nny1][nnx3] - orig.u[nnz2][nny1][nnx3] - orig.u[nnz3][nny1][nnx2] + orig.u[nnz1][nny1][nnx1]) / 4.0;
            d2fdxdz[6]= (orig.u[nnz3][nny2][nnx2] - orig.u[nnz2][nny2][nnx2] - orig.u[nnz3][nny2][nnx1] + orig.u[nnz1][nny2][nnx0]) / 4.0;
            d2fdxdz[7]= (orig.u[nnz3][nny2][nnx3] - orig.u[nnz2][nny2][nnx3] - orig.u[nnz3][nny2][nnx2] + orig.u[nnz1][nny2][nnx1]) / 4.0;

            //d2fdydz
            //d2fdxdy(i,j,k)=(f(i,j+1,k+1)-f(i,j+1,k)-f(i,j,k+1)+f(i,j-1,k-1))/4
            d2fdydz[0]= (orig.u[nnz2][nny2][nnx1] - orig.u[nnz1][nny2][nnx1] - orig.u[nnz2][nny1][nnx1] + orig.u[nnz0][nny0][nnx1]) / 4.0;
            d2fdydz[1]= (orig.u[nnz2][nny2][nnx2] - orig.u[nnz1][nny2][nnx2] - orig.u[nnz2][nny1][nnx2] + orig.u[nnz0][nny0][nnx2]) / 4.0;
            d2fdydz[2]= (orig.u[nnz2][nny3][nnx1] - orig.u[nnz1][nny3][nnx1] - orig.u[nnz2][nny2][nnx1] + orig.u[nnz0][nny1][nnx1]) / 4.0;
            d2fdydz[3]= (orig.u[nnz2][nny3][nnx2] - orig.u[nnz1][nny3][nnx2] - orig.u[nnz2][nny2][nnx2] + orig.u[nnz0][nny1][nnx2]) / 4.0;
                                                                                                                                  
            d2fdydz[4]= (orig.u[nnz3][nny2][nnx1] - orig.u[nnz2][nny2][nnx1] - orig.u[nnz3][nny1][nnx1] + orig.u[nnz1][nny0][nnx1]) / 4.0;
            d2fdydz[5]= (orig.u[nnz3][nny2][nnx2] - orig.u[nnz2][nny2][nnx2] - orig.u[nnz3][nny1][nnx2] + orig.u[nnz1][nny0][nnx2]) / 4.0;
            d2fdydz[6]= (orig.u[nnz3][nny3][nnx1] - orig.u[nnz2][nny3][nnx1] - orig.u[nnz3][nny2][nnx1] + orig.u[nnz1][nny1][nnx1]) / 4.0;
            d2fdydz[7]= (orig.u[nnz3][nny3][nnx2] - orig.u[nnz2][nny3][nnx2] - orig.u[nnz3][nny2][nnx2] + orig.u[nnz1][nny1][nnx2]) / 4.0;

            //d3fdxdydz
            //d3fdxdydz(i,j,k)=      f(i+1,j+1,z+1) -             f(i+1,j,k+1) -             f(i,j+1,k+1) +           f(i-1,j-1,k+1) -           f(i+1,j+1,k-1) +             f(i+1,j,k-1) +             f(i,j+1,k-1) -            f(i-1,j-1,k-1)                     
            d3fdxdydz[0]= (orig.u[nnz2][nny2][nnx2] - orig.u[nnz2][nny1][nnx2] - orig.u[nnz2][nny2][nnx1] + orig.u[nnz2][nny0][nnx0] - orig.u[nnz0][nny2][nnx2] + orig.u[nnz0][nny1][nnx2] + orig.u[nnz0][nny2][nnx1] - orig.u[nnz0][nny0][nnx0]) / 8.0;     
            d3fdxdydz[1]= (orig.u[nnz2][nny2][nnx3] - orig.u[nnz2][nny1][nnx3] - orig.u[nnz2][nny2][nnx2] + orig.u[nnz2][nny0][nnx1] - orig.u[nnz0][nny2][nnx3] + orig.u[nnz0][nny1][nnx3] + orig.u[nnz0][nny2][nnx2] - orig.u[nnz0][nny0][nnx1]) / 8.0;      
            d3fdxdydz[2]= (orig.u[nnz2][nny3][nnx2] - orig.u[nnz2][nny2][nnx2] - orig.u[nnz2][nny3][nnx1] + orig.u[nnz2][nny1][nnx0] - orig.u[nnz0][nny3][nnx2] + orig.u[nnz0][nny2][nnx2] + orig.u[nnz0][nny3][nnx1] - orig.u[nnz0][nny1][nnx0]) / 8.0;      
            d3fdxdydz[3]= (orig.u[nnz2][nny3][nnx3] - orig.u[nnz2][nny2][nnx3] - orig.u[nnz2][nny3][nnx2] + orig.u[nnz2][nny1][nnx1] - orig.u[nnz0][nny3][nnx3] + orig.u[nnz0][nny2][nnx3] + orig.u[nnz0][nny3][nnx2] - orig.u[nnz0][nny1][nnx1]) / 8.0;      
                                                        
            d3fdxdydz[4]= (orig.u[nnz3][nny2][nnx2] - orig.u[nnz3][nny1][nnx2] - orig.u[nnz3][nny2][nnx1] + orig.u[nnz3][nny0][nnx0] - orig.u[nnz1][nny2][nnx2] + orig.u[nnz1][nny1][nnx2] + orig.u[nnz1][nny2][nnx1] - orig.u[nnz1][nny0][nnx0]) / 8.0;      
            d3fdxdydz[5]= (orig.u[nnz3][nny2][nnx3] - orig.u[nnz3][nny1][nnx3] - orig.u[nnz3][nny2][nnx2] + orig.u[nnz3][nny0][nnx1] - orig.u[nnz1][nny2][nnx3] + orig.u[nnz1][nny1][nnx3] + orig.u[nnz1][nny2][nnx2] - orig.u[nnz1][nny0][nnx1]) / 8.0;      
            d3fdxdydz[6]= (orig.u[nnz3][nny3][nnx2] - orig.u[nnz3][nny2][nnx2] - orig.u[nnz3][nny3][nnx1] + orig.u[nnz3][nny1][nnx0] - orig.u[nnz1][nny3][nnx2] + orig.u[nnz1][nny2][nnx2] + orig.u[nnz1][nny3][nnx1] - orig.u[nnz1][nny1][nnx0]) / 8.0;      
            d3fdxdydz[7]= (orig.u[nnz3][nny3][nnx3] - orig.u[nnz3][nny2][nnx3] - orig.u[nnz3][nny3][nnx2] + orig.u[nnz3][nny1][nnx1] - orig.u[nnz1][nny3][nnx3] + orig.u[nnz1][nny2][nnx3] + orig.u[nnz1][nny3][nnx2] - orig.u[nnz1][nny1][nnx1]) / 8.0;      
            
            //scale values
            for (q=0;q<8;q++) {
              f[q]*=1.0;
              dfdx[q]*=dx;
              dfdy[q]*=dy;
              dfdz[q]*=dz;
              d2fdxdy[q]*=dx*dy;
              d2fdxdz[q]*=dx*dz;
              d2fdydz[q]*=dy*dz;
              d3fdxdydz[q]*=dx*dy*dz;
            } 

            //perform tricubic interpolation
            tricubic_get_coeff(a, f, dfdx, dfdy, dfdz, d2fdxdy, d2fdxdz, d2fdydz, d3fdxdydz);
            result= tricubic_eval(a, (xt*nx - nnx1*dx), (yt*ny - nny1*dy), (zt*nz - nnz1*dz));
            tc_interp.u[k][j][i]= (unsigned char)result; 
	    }
      }
    }
	     
    VXwrite(VXout, tc_interp.list);  /* write data */
    VXclose(VXin);                   /* close files */
    VXclose(VXout);
    exit(0);
}
