/*********************************************************************/
/* vtpeak:     Threshold image between two most sig. hgram peaks     */
/*********************************************************************/

#include "VisXV4.h"          /* VisionX structure include file       */
#include "Vutil.h"           /* VisionX utility header files         */
VisXfile_t *VXin,            /* input file structure                 */
           *VXout;           /* output file structure                */
VisXelem_t *VXlist,*vptr;   /* VisionX data structure               */
VXparam_t par[] =            /* command line structure               */
{
{  "if=",   0,   " input file, vtpeak: threshold between hgram peaks"},
{  "of=",   0,   " output file "},
{  "d=",    0,   " min dist between hgram peaks (default 10)"},
{  "-v",    0,   "(verbose) print threshold information"},
{  "b=",  0,  "number of bins"},
{   0,      0,   0} /* list termination */
};
#define  IVAL   par[0].val
#define  OVAL   par[1].val
#define  DVAL   par[2].val
#define  VFLAG  par[3].val
#define  BVAL  par[4].val

main(argc, argv)
int argc;
char *argv[];
{

//VisXimage_t im;                    /* input image structure          */
VisX3dim_t im;
int        i,j,k;                  /* index counters                 */

    int hist[5380];                 /* histogram bins                 */
    int thresh = 0;                /* threshold                      */
    int thresh_accum = 0;          /* accumulator for threshold      */
    int thresh_cnt = 0;
    long long avg1_accum = 0;
    long long avg2_accum = 0;
    int avg1_cnt = 0;
    int avg2_cnt = 0;
    int avg1 = 0;
    int avg2 = 0;
    int avg1_old = 0;
    int avg2_old = 0;
    int maxbin;                    /* maximum histogram bin          */
    int nxtbin;                    /* second maximum bin             */
    int minbin;                    /* minumim histogram bin          */
    int maxa, maxb;       /* second maximum bin above/below maxbin   */
    int dist;                   /* minimum distance between maxima   */
    int bins;
			     
    VXparse(&argc, &argv, par);    /* parse the command line         */
    VXin  = VXopen(IVAL, 0);       /* open input file                */
    VXout = VXopen(OVAL, 1);       /* open the output file           */
    VXlist = VXread(VXin);         /* read file                    */
    
/************ Parameter and initialization section *******************/

    dist = 10;                    /* default dist */
    if (DVAL) dist = atoi(DVAL);  /* if d= was specified, get value */
    if (dist < 0 || dist > 255) {
	fprintf(stderr, "d= must be between 0 and 255\nUsing d=10\n");
        dist = 10;
    }
    
    bins = 256;
    if (BVAL) bins = atoi(BVAL);
    
/************ End of Parameter and initialization section ************/

    if(VXNIL == (vptr = VXfind(VXlist, VX_PSHORT))){
       fprintf(stderr, "linterp: no byte images found\n");
       exit(1);}
       
    VXset3dim(&im, vptr, VXin); /* initialize input structure   */

/***************** Application specific section **********************/
      
            /* clear the histogram */
            for (i = 0; i < bins; i++) {
                hist[i] = 0;
            }

            /* compute the histogram */

            for(k = im.zlo; k <= im.zhi; k++) {
                for (i = im.ylo; i <= im.yhi; i++) {
                    for (j = im.xlo; j <= im.xhi; j++) {
                        if(im.s[k][i][j] >= bins) {
                            hist[bins-1]++;
                        }
                        else {
                            hist[im.s[k][i][j]]++;
                        }
                    }
                }
            }
            

             /* compute the average value */
             for (i = 0; i < bins; i++) {
                thresh_accum = thresh_accum + hist[i] * i;
                thresh_cnt = thresh_cnt + hist[i];
             }
             thresh = thresh_accum/thresh_cnt;

        
            for (k = 0; k < 100; k++) {
                /* determine the new threshold */
                for (i = 0; i < bins; i++) {
                    if (i >= thresh) {
                       avg1_accum = avg1_accum + hist[i] * i;
                       avg1_cnt = avg1_cnt + hist[i];
                    }
                    else {
                       avg2_accum = avg2_accum + hist[i] * i;
                       avg2_cnt = avg2_cnt + hist[i];
                    }
                }
                if (avg1_cnt == 0) {
                    avg1 = 0;
                }
                else {
                    avg1 = (avg1_accum/avg1_cnt);
                }
                if (avg2_cnt == 0) {
                    avg2 = 0;
                }
                else{
                    avg2 = (avg2_accum/avg2_cnt);
                }
                if (avg1 == avg1_old && avg2 == avg2_old) {
                    fprintf(stderr, "avgs are the same\n");
                    k = 100;
                    break;
                }
                thresh = (avg1 + avg2)/2;
                fprintf(stderr, "avg1: %d, avg2: %d \n", (avg1), (avg2));
                fprintf(stderr, "thresh:%d\n", thresh);
                avg1_accum = 0;
                avg2_accum = 0;
                avg1_cnt = 0;
                avg2_cnt = 0;
                avg1_old = avg1;
                avg2_old = avg2;
            }


            fprintf(stderr, "final thresh:%d\n", thresh);

            /* apply the threshold */
            for (k = im.zlo; k <= im.zhi; k++) {
                for (i = im.ylo; i <= im.yhi; i++) {
                    for (j = im.xlo; j <= im.xhi; j++) {
                        if (im.s[k][i][j] >= thresh) im.s[k][i][j] = 255;
                        else                           im.s[k][i][j] = 0;
                    }
                }
            }

/************** End of the Application specific section **************/

    VXwrite(VXout, im.list);  // write data  
    VXclose(VXin);               	 // close files 
    VXclose(VXout);
    exit(0);
}
