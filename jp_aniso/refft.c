/*
 *      correct the sign error in stanford fft.c program
 *      normalized for the inverse transform i.e. (sign < 0.)
 *      jeff g. s. pan          feb '85
 */
/*
 * AS WRITTEN THIS CALCULATES Y(f)= SUM X(t)exp(+i omega t) when asign=1
 */
#include "stdio.h"
#include "math.h"
char *alloc (size)
int size;
/*
 *   allocation with error detection
 */
   {
   char *ptr, *calloc();
   if ((ptr = calloc (size,1)) <= 0)
      err_("cant allocate %d bytes \n",size);
   return (ptr);
   }

refft_(x,ane,asign,amode)
register struct complex {float re, im;} *x;
int *ane, *asign, *amode;
/*
 *   radix 2 real <=> complex fast fourier transform
 *   x is the in place vector, n is the power of two transform length
 *      sign:
 *              1   forward transform
 *             -1   backward transform (will preserve parceval's thrm)
 *      mode:
 *      1   real to complex; real nyquist part in imaginary dc
 *      -1   complex to real; real nyquist part expected in imaginary dc
 *      2   real to complex; unpacked
 *      -2   complex to real; unpacked
 *   scaling for postive modes
 */
   {
   register struct complex *xp, *yp;
   int n, sign, mode, n2;
   float cn, sn, cd, sd, arg, scale;
   float aa, bb, ab, ba, real, imag;
   double sin(), atan2();

   n = *ane;
   sign = *asign;
   mode = *amode;
   n2 = n / 2;
   if (mode > 0) cefft_(x,&n2,&sign);

   /* pack */
   /* rfft butterfly */
   if (mode == -2) x->im = x[n/2].re;

   sn = 0.0;
   cn = 1.0;
   arg = 2. * atan2 (1.,0.) / n;
   aa = sin (arg);
   cd = 2.0 * aa * aa;
   sd = sin (arg+arg);
   if (sign < 0) sd = -sd; /* change sign for test */
   aa = x->re;
   bb = x->im;
   if (sign > 0)
      {
      x->re = aa + bb;
      x->im = aa - bb;
      }
   else
      {
      x->re = (aa + bb) * .5;
      x->im = (bb - aa) * .5;
      }
   for (xp=x+1, yp=x+n/2-1; xp<=yp; xp++, yp--)
      {
      aa = cd * cn + sd * sn;
      sn += sd * cn - cd * sn;
      cn -= aa;
      aa = (xp->re + yp->re) * .5;
      ab = (xp->re - yp->re) * .5;
      ba = (xp->im + yp->im) * .5;
      bb = (xp->im - yp->im) * .5;
      real = cn * ba + sn * ab;
      imag = sn * ba - cn * ab;
      yp->im = imag - bb;
      xp->im = imag + bb;
      yp->re = aa - real;
      xp->re = aa + real;
      }

   if (0 > mode)
      {
      cefft_(x,&n2,&sign);
      for (xp=x, yp=x+n/2; xp<yp; xp++) xp->im *= -1.;
      }

   /* unpack */
   if (mode == 2)
      {
      x[n/2].re = x[0].im;
      x[0].im = x[n/2].im = 0.;
      }
   }

cefft_(x,ane,asign)
register struct complex {float re, im;} *x;
int *ane, *asign;
/*
 *   radix 2 complex <=> complex fourier transform
 *   n is the number of complex points in x
 *   scaling on positive sign
 */
   {
   register struct complex *xp, *yp;
   struct complex *end;
   float cn, sn, cd, sd, real, imag, scale, *psintab;
   float temp;
   int i=0, j=0, step=0, n, sign;
   /* sines from pi/2 to pi/2097152 */
   static float sintab[] =
      {
      1.0000000000000000e-00,
      7.0710678118654747e-01,
      3.8268343236508974e-01,
      1.9509032201612825e-01,
      9.8017140329560596e-02,
      4.9067674327418010e-02,
      2.4541228522912286e-02,
      1.2271538285719925e-02,
      6.1358846491544749e-03,
      3.0679567629659760e-03,
      1.5339801862847655e-03,
      7.6699031874270447e-04,
      3.8349518757139556e-04,
      1.9174759731070329e-04,
      9.5873799095977337e-05,
      4.7936899603066881e-05,
      2.3968449808418217e-05,
      1.1984224905069705e-05,
      5.9921124526424274e-06,
      2.9960562263346605e-06

      };
   n = *ane;
   sign = *asign;

   /* bit reverse address swapping */
   for (xp=x, end=x+n, i=0; xp<end; xp++, i+=j)
      {
      if (xp < (yp=x+i))
         {
         temp = yp->re;
         yp->re = xp->re;
         xp->re = temp;
         temp = yp->im;
         yp->im = xp->im;
         xp->im = temp;
         }
      for (j=n>>1; j>=1 && i>=j;)
         {
         i -= j;
         j >>= 1;
         }
      }

   /* first butterfly */
   if (sign < 0) scale =  1. / n;
   for (xp=x, yp=x+1; xp<end; xp+=2, yp+=2)
      {
      if (sign < 0)
         {
         xp->re *= scale;
         xp->im *= scale;
         yp->re *= scale;
         yp->im *= scale;
         }
      temp = yp->re;
      yp->re = xp->re - temp;
      xp->re += temp;
      temp = yp->im;
      yp->im = xp->im - temp;
      xp->im += temp;
      }

   /* remaining butterflies */
   for (i=2, psintab=sintab; i<n; i=step)
      {
      step = i << 1;
      sd = *psintab++;
      temp = *psintab;
      cd = 2.0 * temp * temp;
      cn = 1.0;
      sn = 0.0;
      if (sign < 0) sd = -sd;
      for (j=0; j<i; j++)
         {
         for (xp=x+j; xp<end; xp+=step)
            {
            yp = xp + i;
            real = cn * yp->re - sn * yp->im;
            imag = sn * yp->re + cn * yp->im;
            yp->re = xp->re - real;
            yp->im = xp->im - imag;
            xp->re += real;
            xp->im += imag;
            }
         temp = cd * cn + sd * sn;
         sn += sd * cn - cd * sn;
         cn -=temp;
         }
      }
   }

err_(a,b,c,d,e,f,g,h)
char *a, *b, *c, *d, *e, *f, *g, *h;
/*
 *   error abortion subroutine
 *   a is a printf format string, while b-h are optional arguments
 *   I CANT FIND THIS ROUTINE SO COMMENT IT OUT -- JPARK 7/24/86
 */
   {
   printf("error \n");
/*   fprintf (stderr,a,b,c,d,e,f,g,h);*/
   exit (-1);
   }

