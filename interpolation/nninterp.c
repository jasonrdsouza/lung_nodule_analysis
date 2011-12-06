/****************************************************************/
/* Example VisX4 program nninterp		                        */
/*  Compute nearest neighbor interpolation of the input 3d image*/
/* Syntax:                                                      */
/*        interp if=infile of=outfile dx= dy= dz= nx= ny= nz=   */
/****************************************************************/
/****************************************************************/
/*  This Software is Copyright (c) 1996, 2000 Anthony P. Reeves */
/*  All rights reserved.                                        */
/****************************************************************/

#include "VisXV4.h"          /* VisX structure include file     */
#include "Vutil.h"           /* VisX utility header files       */

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

main(argc, argv)
int argc;
char *argv[];
{

VisX3dim_t orig;                   /* input v3dim structure        */
VisX3dim_t nn_interp;              /* output v3dim structure       */
VisXelem_t *vptr;
int	   i,j,k;                  /* index counters               */
int xt, yt, zt;
unsigned char t;
float      dx;
float      dy;
float      dz;
float      nx;
float      ny;
float      nz;
int	   nnx;
int	   nny;
int	   nnz;

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

    float bbx[6]= {orig.bbx[0]*dx/nx, orig.bbx[1]*dx/nx, orig.bbx[2]*dy/ny, orig.bbx[3]*dy/ny, orig.bbx[4]*dz/nz, orig.bbx[5]*dz/nz};
    VXmake3dim(&nn_interp, VX_PBYTE, bbx, 1);

    if(VFLAG){
       fprintf(stderr,"bbx is %f %f %f %f %f %f\n", orig.bbx[0],
           orig.bbx[1],orig.bbx[2],orig.bbx[3],orig.bbx[4],orig.bbx[5]);
    }
    /* simple pixel computation -- notice the order of the loops */
    for (k = nn_interp.zlo; k <= nn_interp.zhi; k++){
      for (j = nn_interp.ylo; j <= nn_interp.yhi; j++){
        for (i = nn_interp.xlo; i <= nn_interp.xhi; i++){
	    // remap from the new voxel space to the old voxel space
            xt= i;
            yt= j;
            zt= k;

            // find the nearest pixels
            nnx= (xt*nx/dx < orig.xlo) ? ceil(xt*nx/dx)  :
                 (xt*nx/dx > orig.xhi) ? floor(xt*nx/dx) :
                                         round(xt*nx/dx);

            nny= (yt*ny/dy < orig.ylo) ? ceil(yt*ny/dy)  :
                 (yt*ny/dy > orig.yhi) ? floor(yt*ny/dy) :
                                         round(yt*ny/dy);

            nnz= (zt*nz/dz < orig.zlo) ? ceil(zt*nz/dz)  :
                 (zt*nz/dz > orig.zhi) ? floor(zt*nz/dz) :
                                         round(zt*nz/dz);
            // do the interpolation
            nn_interp.u[k][j][i]= orig.u[nnz][nny][nnx];
	}
      }
    }
	     
    VXwrite(VXout, nn_interp.list);       /* write data                   */
    VXclose(VXin);               	 /* close files                   */
    VXclose(VXout);
    exit(0);
}
