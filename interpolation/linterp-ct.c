/****************************************************************/
/* Example VisX4 program linterp                                */
/*  Compute linear interpolation of the input 3d image          */
/* Syntax:                                                      */
/*       linterp if=infile of=outfile dx= dy= dz= nx= ny= nz=   */
/****************************************************************/
/****************************************************************/
/*  This Software is Copyright (c) 1996, 2000 Anthony P. Reeves */
/*  All rights reserved.                                        */
/****************************************************************/

// THIS VERSION WORKS FOR SHORT (CT) IMAGES

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

main(argc, argv)
int argc;
char *argv[];
{

VisX3dim_t orig;                   /* input v3dim structure        */
VisX3dim_t l_interp;              /* output v3dim structure       */
VisXelem_t *vptr;
int	   i,j,k;                  /* index counters               */
unsigned char t;
float      dx;
float      dy;
float      dz;
float      nx;
float      ny;
float      nz;
int	   nnxlo, nnxhi, xt;
int	   nnylo, nnyhi, yt;
int	   nnzlo, nnzhi, zt;
float      xd, yd, zd, x1, x2, y1, y2, w1, w2, p; //temp variables for linear interp

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

    if(VXNIL == (vptr = VXfind(VXlist, VX_PSHORT))){
       fprintf(stderr, "linterp: no byte images found\n");
       exit(1);}

    VXset3dim(&orig, vptr, VXin);    /* initialize input structure   */

    float bbx[6]= {ceil(orig.bbx[0]*dx/nx), floor(orig.bbx[1]*dx/nx), ceil(orig.bbx[2]*dy/ny), floor(orig.bbx[3]*dy/ny), ceil(orig.bbx[4]*dz/nz), floor(orig.bbx[5]*dz/nz)};
    VXmake3dim(&l_interp, VX_PSHORT, bbx, 1);

    if(VFLAG){
       fprintf(stderr,"bbx is %f %f %f %f %f %f\n", orig.bbx[0],
           orig.bbx[1],orig.bbx[2],orig.bbx[3],orig.bbx[4],orig.bbx[5]);
    }

    /* simple pixel computation -- notice the order of the loops */
    for (k = l_interp.zlo + dz/nz; k <= l_interp.zhi - dz/nz; k++){
      for (j = l_interp.ylo + dy/ny; j <= l_interp.yhi - dy/ny; j++){
        for (i = l_interp.xlo + dx/nx; i <= l_interp.xhi - dx/nx; i++){
	    // remap from the new voxel space to the old voxel space
            xt= i;
            yt= j;
            zt= k;

            nnxlo= (xt*nx/dx < orig.xlo ? ceil(xt*nx/dx) : floor(xt*nx/dx));
            nnylo= (yt*ny/dy < orig.ylo ? ceil(yt*ny/dy) : floor(yt*ny/dy));
            nnzlo= (zt*nz/dz < orig.zlo ? ceil(zt*nz/dz) : floor(zt*nz/dz));

            nnxhi= (xt*nx/dx > orig.xhi ? floor(xt*nx/dx) : ceil(xt*nx/dx));
            nnyhi= (yt*ny/dy > orig.yhi ? floor(yt*ny/dy) : ceil(yt*ny/dy));
            nnzhi= (zt*nz/dz > orig.zhi ? floor(zt*nz/dz) : ceil(zt*nz/dz));

            // 8 NN pixels
            orig.s[nnzlo][nnylo][nnxlo];
            orig.s[nnzlo][nnylo][nnxhi];
            orig.s[nnzlo][nnyhi][nnxlo];
            orig.s[nnzlo][nnyhi][nnxhi];

            orig.s[nnzhi][nnylo][nnxlo];
            orig.s[nnzhi][nnylo][nnxhi];
            orig.s[nnzhi][nnyhi][nnxlo];
            orig.s[nnzhi][nnyhi][nnxhi];

            // do interpolation cacluations
            // distance calculation - map lo and hi values
            xd= xt*nx/dx - nnxlo;
            yd= yt*ny/dy - nnylo;
            zd= zt*nz/dz - nnzlo;
            if (xd < 0) xd *= -1;
            if (yd < 0) yd *= -1;
            if (zd < 0) zd *= -1;

            // interpolate in the first direction
	    x1= orig.s[nnzlo][nnylo][nnxlo] * (1-zd) + orig.s[nnzhi][nnylo][nnxlo] * zd;
	    x2= orig.s[nnzlo][nnyhi][nnxlo] * (1-zd) + orig.s[nnzhi][nnyhi][nnxlo] * zd;
	    y1= orig.s[nnzlo][nnylo][nnxhi] * (1-zd) + orig.s[nnzhi][nnylo][nnxhi] * zd;
	    y2= orig.s[nnzlo][nnyhi][nnxhi] * (1-zd) + orig.s[nnzhi][nnyhi][nnxhi] * zd;

            // the second
            w1= x1 * (1-yd) + x2 * yd;
            w2= y1 * (1-yd) + y2 * yd;

            // final interpolation
            l_interp.s[zt][yt][xt]= w1 * (1-xd) + w2 * xd;
	}
      }
    }
	     
    VXwrite(VXout, l_interp.list);  // write data  
    VXclose(VXin);               	 // close files 
    VXclose(VXout);
    exit(0);
}
