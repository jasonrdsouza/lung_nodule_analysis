/*********************************************************************/
/* vits:     Iterative threshold algorithm                           */
/*********************************************************************/

#include "VisXV4.h"          /* VisionX structure include file       */
#include "Vutil.h"           /* VisionX utility header files         */
VisXfile_t *VXin,            /* input file structure                 */
           *VXout;           /* output file structure                */
VisXelem_t *VXlist,*VXptr;   /* VisionX data structure               */
VXparam_t par[] =            /* command line structure               */
{
{  "if=",   0,   " input file, vtpeak: threshold between hgram peaks"},
{  "of=",   0,   " output file "},
{  "t=",    0,   " threshold value"},
{  "-v",    0,   "(verbose) print threshold information"},
{   0,      0,   0} /* list termination */
};
#define  IVAL   par[0].val
#define  OVAL   par[1].val
#define  TVAL   par[2].val
#define  VFLAG  par[3].val

main(argc, argv)
int argc;
char *argv[];
{

VisXimage_t im;                    /* input image structure          */
int        i,j;                    /* index counters                 */

    int hist[256];                 /* histogram bins                 */
    int thresh;                    /* threshold                      */

    int imgcount= 0;
    int imgavg= 0;                 /* average value of the image */
    int temp1= 0;                  /* temp value for calculating average values of R1*/
    int temp2= 0;                  /* temp value for calculating average values of R2*/
    int num_r1= 0;                 /* number of pixels in R1 */
    int num_r2= 0;                 /* number of pixels in R2 */
    int avg1= 0;                   /* average value of R1 (low) pixels */ 
    int avg2= 1;                   /* average value of R2 (low) pixels */
    int oldavg1= 1;                /* old average value of R1 */
    int oldavg2= 0;                /* old average value of R2 */
                                   // these need to be different so that we don't skip the while loop 			     
    VXparse(&argc, &argv, par);    /* parse the command line         */
    VXin  = VXopen(IVAL, 0);       /* open input file                */
    VXout = VXopen(OVAL, 1);       /* open the output file           */

    while((VXlist = VXptr = VXreadframe(VXin)) != VXNIL){ /* every frame */
        VXfupdate(VXout, VXin); /* update global constants */
	/* find next byte image */
        while (VXNIL != (VXptr = VXfind(VXptr, VX_PBYTE)))  { 
            VXsetimage(&im, VXptr, VXin); /* initialize input structure */

/***************** Application specific section **********************/

            /* clear the histogram */
            for (i = 0; i < 256; i++) hist[i] = 0;
 
            /* compute the histogram */
            for (i = im.ylo; i <= im.yhi; i++)
                for (j = im.xlo; j <= im.xhi; j++)
                    hist[im.u[i][j]]++;

                // compute average value of image
            imgavg= 0;
            imgcount= 0;
            for (i= 0; i < 256; i++) {
                imgavg += i*hist[i];
                imgcount += hist[i];
            }

            imgavg= imgavg / imgcount;

            // set threshold
            thresh= imgavg;
            if (TVAL) thresh= atoi(TVAL);
            if (thresh < 0 || thresh > 255) thresh= imgavg;

            if(VFLAG) fprintf(stderr, "Starting threshold = %d\n", thresh);

            /* compute the threshold */
            while ((avg1 != oldavg1) || (avg2 != oldavg2)) {
                temp1= 0; // clear values
                temp2= 0;
                num_r1= 0;
                num_r2= 0;

                oldavg1= avg1; // update old averages
                oldavg2= avg2; 

                for (i = 0; i < 256; i++) {
                    if (i > thresh) {
                        temp1 += i * hist[i];
                        num_r1 += hist[i];
                    } else if (i < thresh) {
                        temp2 += i * hist[i];
                        num_r2 += hist[i];
                    } else {}
                }

                avg1= (num_r1 == 0 ? thresh : temp1 / num_r1);
                avg2= (num_r2 == 0 ? thresh : temp2 / num_r2);
                thresh= (avg1 + avg2) / 2;

                if(VFLAG) fprintf(stderr, "Iteration threshold = %d\n", thresh);
            }
			
  
            if(VFLAG) fprintf(stderr, "End.\n");
  
            /* apply the threshold */
            for (i = im.ylo; i <= im.yhi; i++) {
                for (j = im.xlo; j <= im.xhi; j++) {
                    if (im.u[i][j] >= thresh) im.u[i][j] = 255;
                    else                      im.u[i][j] = 0;
                }
            }
 
/************** End of the Application specific section **************/

            VXresetimage(&im); /* free the im image structure  */
            VXptr = VXptr->next; /* move to the next image */
        } /* end of every image section */
        VXwriteframe(VXout,VXlist); /* write frame */
        VXdellist(VXlist);          /* delete the frame */
    } /* end of every frame section */
    VXclose(VXin);  /* close files */
    VXclose(VXout);
    exit(0);
}
