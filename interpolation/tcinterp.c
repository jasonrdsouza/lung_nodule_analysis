/****************************************************************/
/* Example VisX4 program tcinterp                               */
/*  Compute tricubic interpolation of the input 3d image        */
/* Syntax:                                                      */
/*       tcinterp if=infile of=outfile dx= dy= dz= nx= ny= nz=  */
/****************************************************************/
/****************************************************************/
/*  This Software is Copyright (c) 1996, 2000 Anthony P. Reeves */
/*  All rights reserved.                                        */
/****************************************************************/

#include "VisXV4.h"          /* VisX structure include file     */
#include "Vutil.h"           /* VisX utility header files       */
#include <math.h>

VisXfile_t *VXin,            /* input file structure            */
           *VXout;           /* output file structure           */
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

//one-dimensional cubic interpolation
double CINT(double y0, double y1, double y2, double y3, double mu){
        double a0, a1, a2, a3, mu2;
        
        mu2= mu*mu;
        a0= y3 - y2 - y0 + y1;
        a1= y0 - y1 - a0;
        a2= y2 - y0;
        a3= y1;

        return a0*mu*mu2 + a1*mu2 + a2*mu + a3;
}

main(argc, argv)
int argc;
char *argv[];
{

VisX3dim_t orig;                   /* input v3dim structure        */
VisX3dim_t tc_interp;              /* output v3dim structure       */
VisXelem_t *vptr;
int	   i,j,k;                  /* index counters               */
unsigned char t;
float      dx;
float      dy;
float      dz;
float      nx;
float      ny;
float      nz;
int	   nnx0, nnx1, nnx2, nnx3, xt;
int	   nny0, nny1, nny2, nny3, yt;
int	   nnz0, nnz1, nnz2, nnz3, zt;
float      xd, yd, zd, p00, p01, p02, p03, p10, p11, p12, p13, p20, p21, p22, p23, p30, p31, p32, p33, w0, w1, w2, w3; //temp variables for linear interp

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
       fprintf(stderr, "linterp: no byte images found\n");
       exit(1);}

    VXset3dim(&orig, vptr, VXin);    /* initialize input structure   */

    float bbx[6]= {orig.bbx[0]*dx/nx, orig.bbx[1]*dx/nx, orig.bbx[2]*dy/ny, orig.bbx[3]*dy/ny, orig.bbx[4]*dz/nz, orig.bbx[5]*dz/nz};
    VXmake3dim(&tc_interp, VX_PBYTE, bbx, 1);

  
    if(VFLAG){
       fprintf(stderr,"bbx is %f %f %f %f %f %f\n", orig.bbx[0],
           orig.bbx[1],orig.bbx[2],orig.bbx[3],orig.bbx[4],orig.bbx[5]);
    }

    /* simple pixel computation -- notice the order of the loops */
    for (k = tc_interp.zlo+2*dz/nz; k <= tc_interp.zhi-2*dz/nz; k++){
      for (j = tc_interp.ylo+2*dy/ny; j <= tc_interp.yhi-2*dy/ny; j++){
        for (i = tc_interp.xlo+2*dx/nx; i <= tc_interp.xhi-2*dx/nx; i++){
            // remap from the new voxel space to the old voxel space
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
            orig.u[nnz1][nny1][nnx1];
            orig.u[nnz1][nny1][nnx2];
            orig.u[nnz1][nny2][nnx1];
            orig.u[nnz1][nny2][nnx2];

            orig.u[nnz2][nny1][nnx1];
            orig.u[nnz2][nny1][nnx2];
            orig.u[nnz2][nny2][nnx1];
            orig.u[nnz2][nny2][nnx2];

            // start interpolation calcuation
            // distance calculation - map lo and hi values
            xd= xt*nx/dx - nnx1;
            yd= yt*ny/dy - nny1;
            zd= zt*nz/dz - nnz1;
            if (xd < 0) xd *= -1;
            if (yd < 0) yd *= -1;
            if (zd < 0) zd *= -1;
	
            // interpolate in the Z direction
            p00= CINT(orig.u[nnz0][nny0][nnx0], orig.u[nnz1][nny0][nnx0], orig.u[nnz2][nny0][nnx0], orig.u[nnz3][nny0][nnx0], zd);
            p01= CINT(orig.u[nnz0][nny1][nnx0], orig.u[nnz1][nny1][nnx0], orig.u[nnz2][nny1][nnx0], orig.u[nnz3][nny1][nnx0], zd);
            p02= CINT(orig.u[nnz0][nny2][nnx0], orig.u[nnz1][nny2][nnx0], orig.u[nnz2][nny2][nnx0], orig.u[nnz3][nny2][nnx0], zd);
            p03= CINT(orig.u[nnz0][nny3][nnx0], orig.u[nnz1][nny3][nnx0], orig.u[nnz2][nny3][nnx0], orig.u[nnz3][nny3][nnx0], zd);
            p10= CINT(orig.u[nnz0][nny0][nnx1], orig.u[nnz1][nny0][nnx1], orig.u[nnz2][nny0][nnx1], orig.u[nnz3][nny0][nnx1], zd);
            p11= CINT(orig.u[nnz0][nny1][nnx1], orig.u[nnz1][nny1][nnx1], orig.u[nnz2][nny1][nnx1], orig.u[nnz3][nny1][nnx1], zd);
            p12= CINT(orig.u[nnz0][nny2][nnx1], orig.u[nnz1][nny2][nnx1], orig.u[nnz2][nny2][nnx1], orig.u[nnz3][nny2][nnx1], zd);
            p13= CINT(orig.u[nnz0][nny3][nnx1], orig.u[nnz1][nny3][nnx1], orig.u[nnz2][nny3][nnx1], orig.u[nnz3][nny3][nnx1], zd);
            p20= CINT(orig.u[nnz0][nny0][nnx2], orig.u[nnz1][nny0][nnx2], orig.u[nnz2][nny0][nnx2], orig.u[nnz3][nny0][nnx2], zd);
            p21= CINT(orig.u[nnz0][nny1][nnx2], orig.u[nnz1][nny1][nnx2], orig.u[nnz2][nny1][nnx2], orig.u[nnz3][nny1][nnx2], zd);
            p22= CINT(orig.u[nnz0][nny2][nnx2], orig.u[nnz1][nny2][nnx2], orig.u[nnz2][nny2][nnx2], orig.u[nnz3][nny2][nnx2], zd);
            p23= CINT(orig.u[nnz0][nny3][nnx2], orig.u[nnz1][nny3][nnx2], orig.u[nnz2][nny3][nnx2], orig.u[nnz3][nny3][nnx2], zd);
            p30= CINT(orig.u[nnz0][nny0][nnx3], orig.u[nnz1][nny0][nnx3], orig.u[nnz2][nny0][nnx3], orig.u[nnz3][nny0][nnx3], zd);
            p31= CINT(orig.u[nnz0][nny1][nnx3], orig.u[nnz1][nny1][nnx3], orig.u[nnz2][nny1][nnx3], orig.u[nnz3][nny1][nnx3], zd);
            p32= CINT(orig.u[nnz0][nny2][nnx3], orig.u[nnz1][nny2][nnx3], orig.u[nnz2][nny2][nnx3], orig.u[nnz3][nny2][nnx3], zd);
            p33= CINT(orig.u[nnz0][nny3][nnx3], orig.u[nnz1][nny3][nnx3], orig.u[nnz2][nny3][nnx3], orig.u[nnz3][nny3][nnx3], zd);
            
            // Y direction
            w0= CINT(p00, p01, p02, p03, yd);
            w1= CINT(p10, p11, p12, p13, yd);
            w2= CINT(p20, p21, p22, p23, yd);
            w3= CINT(p30, p31, p32, p33, yd);

            // X direction
            tc_interp.u[zt][yt][xt]= CINT(w0, w1, w2, w3, xd);
	    }
      }
    }
	     
    VXwrite(VXout, tc_interp.list);  // write data  
    VXclose(VXin);               	 // close files 
    VXclose(VXout);
    exit(0);
}
