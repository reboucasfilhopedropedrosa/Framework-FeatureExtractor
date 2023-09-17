/******************************************************************************
Copyright 1994 by the Massachusetts Institute of Technology.  All
rights reserved.

Developed by Thomas P. Minka and Rosalind W. Picard at the Media
Laboratory, MIT, Cambridge, Massachusetts, with support from BT, PLC,
Hewlett-Packard, and NEC.

This distribution is approved by Nicholas Negroponte, Director of
the Media Laboratory, MIT.

Permission to use, copy, or modify this software and its documentation
for educational and research purposes only and without fee is hereby
granted, provided that this copyright notice and the original authors'
names appear on all copies and supporting documentation.  If individual
files are separated from this distribution directory structure, this
copyright notice must be included.  For any other uses of this software,
in original or modified form, including but not limited to distribution
in whole or in part, specific prior permission must be obtained from
MIT.  These programs shall not be used, rewritten, or adapted as the
basis of a commercial software or hardware product without first
obtaining appropriate licenses from MIT.  MIT. makes no representations
about the suitability of this software for any purpose.  It is provided
"as is" without express or implied warranty.
******************************************************************************/

/******************************************************************************
This program is a sample implementation of the texture features in

Tamura, H., Mori, S., and Yamawaki, T.
"Textural features corresponding to visual perception",
IEEE Trans. SMC, 8(6):460-473, 1978.

It is influenced by correspondence with William Equitz of the QBIC team.
It outputs a four-element feature vector. Each feature has a different dynamic
range, so similarity should be computed using Mahalanobis distance or
simplified Mahalanobis distance (i.e. with a diagonal inverse covariance 
matrix).

Compile with:

cc -o tamura tamura.c basic.c -lm

This assumes matrix.h is in the current directory.
******************************************************************************/
#include "tamura.h"
#include <math.h>
#include <cstring>

/* Coarseness */
Matrix nbr_avg(Matrix image, int factor);
Matrix nbr_diff(Matrix image, int factor);
//double Coarseness(Matrix image, int max_scale, double weight);

/* Contrast */
//double Contrast(Matrix image);

/* Directionality */
//void Directionality(Matrix image, int histo_bins, int histo_thres, double peak_thres, int *orientations, double *dir);
double *dir_histogram(Matrix image, int histo_bins, int threshold);
double *new_dir_histogram(Matrix image, int histo_bins);

//void main(int argc, char *argv[])
//{
//  Matrix image;
//  int orientations;
//  double crs, con, dir;
//  char *p,*q;
//  FILE *fp;
//
//  if(argc < 2) {
//    printf("Usage:\n%s <image_file>\n", argv[0]);
//    printf("To compute tamura features for a 128x128 character array\n");
//    printf("(coarseness, contrast, directionality, #orientations)\n");
//    exit(1);
//  }
//
//  image = MatrixCreate(128, 128);
//  fp = fopen(argv[1], "r");
//  if(!fp) { perror(""); exit(1); }
//  MatrixReadImage(image, fp);
//  fclose(fp);
//
//  crs = Coarseness(image, 5, 1);
//  con = Contrast(image);
//  Directionality(image, 16, 12, 2, &orientations, &dir);
//  
//  free(image);
//
//  for(p = q = argv[1];*q;q++) if(*q == '/') p=q+1;
//  printf("%s %lf %lf %lf %d\n", p, crs, con, dir, orientations);
//}

#define CHECK_EXTREMUM(a,b)            \
  if((diff[a] * diff[b] < 0) ||        \
     ((diff[a] * diff[b] == 0) &&      \
      (diff[a] + diff[b] != 0))) {     \
    extrema[n_extrema] = b;            \
    is_peak[n_extrema] = diff[b] < 0;  \
    n_extrema++;                       \
  }

int valid_peak(int *peak, double *histogram, double peak_thres)
{
  double max_valley, sharpness;

  max_valley = histogram[peak[0]];
  if(histogram[peak[2]] > max_valley) 
    max_valley = histogram[peak[2]];
  sharpness = histogram[peak[1]] / max_valley;
  max_valley = fabs(sharpness - peak_thres);
/*
  if(max_valley < 0.1) printf("close call\n");
*/
  return(sharpness > peak_thres);
}

double sd_histo(double *histogram, int histo_bins, int start, int stop)
{
  double *interval;
  double sum;
  int i, length, mp;

  length = stop - start;
  if(stop <= start) length += histo_bins;
  interval = Allocate(length, double);

  /* copy the sub-histogram into the interval */
  if(stop <= start) {
    for(i=start;i < histo_bins;i++) interval[i-start] = histogram[i];
    for(i=0;i<stop;i++) interval[histo_bins-start+i] = histogram[i];
  }
  else {
    for(i=start;i<stop;i++) interval[i-start] = histogram[i];
  }

  /* normalize to make a distribution */
  sum = 0.0;
  for(i=0;i<length;i++) sum += interval[i];
  for(i=0;i<length;i++) interval[i] /= sum;

  /* find mean */
  sum = 0.0;
  for(i=0;i<length;i++) {
    sum += (i+1)*interval[i];
  }
  mp = (int)(sum-0.5);

  /* compute variance */
  sum = 0.0;
  for(i=0;i<length;i++) {
    sum += (double)(i - mp) * (i - mp) * interval[i];
  }

  free(interval);
  return sqrt(sum);
}
  
void Directionality(Matrix image, int histo_bins, int histo_thres, 
		    double peak_thres,
		    int *orientations, double *dir)
{
  double *histogram, *diff;
  int *extrema, n_extrema, *is_peak;
  int num_peaks, i, j, k;
  int **pv;

  /* compute direction histogram */
  histogram = new_dir_histogram(image, histo_bins);
  //histogram = dir_histogram(image, histo_bins, histo_thres);

/*
  for(i=0;i<histo_bins;i++) printf("%lf\n", histogram[i]);
*/

  /* find extrema of histogram */
  extrema = Allocate(histo_bins,int);
  is_peak = Allocate(histo_bins,int);
  diff = Allocate(histo_bins,double);
  n_extrema = 0;

  /* compute difference vector */
  for(i=0;i < histo_bins-1;i++) 
    diff[i] = histogram[i+1] - histogram[i];
  diff[histo_bins-1] = histogram[0] - histogram[histo_bins-1];
  CHECK_EXTREMUM(histo_bins-1, 0);
  for(i=0;i < histo_bins-1;i++) {
    CHECK_EXTREMUM(i, i+1);
  }
  free(diff);
/*
  for(i=0;i<n_extrema;i++) printf("%d\n", extrema[i]);
*/

  /* extract salient peaks */
  pv = Allocate(n_extrema, int*);
  pv[0] = Allocate(3*n_extrema, int);
  for(i=1;i < n_extrema;i++) pv[i] = &pv[0][i * 3];
  
  j = 0;
  num_peaks = 0;
  for(i=0;i < n_extrema;i++) {
    if(is_peak[i]) {
      if(j == 0) {
	/* find valley to fill first slot */
        for(k=n_extrema-1;is_peak[k];k--);
	pv[num_peaks][0] = extrema[k];
	j++;
      }
      if(j == 1) {
	pv[num_peaks][1] = extrema[i];
	j++;
      }
    }
    else {
      if(j < 2) {
	pv[num_peaks][0] = extrema[i];
	j = 1;
      }
      else {
	pv[num_peaks][2] = extrema[i];
	if(valid_peak(pv[num_peaks], histogram, peak_thres)) {
	  num_peaks++;
	}
	pv[num_peaks][0] = extrema[i];
	j = 1;
      }
    }
  }
  if(j == 2) {
    for(k=0;is_peak[k];k++);
    pv[num_peaks][2] = extrema[k];
  }
  if((j < 2) || !valid_peak(pv[num_peaks], histogram, peak_thres)) 
    num_peaks--;
  num_peaks++;
/*
  for(i=0;i<num_peaks;i++) 
    printf("%d %d %d\n", pv[i][0], pv[i][1], pv[i][2]);
*/
  *orientations = num_peaks;

  /* compute 2nd moment about peaks */
  if(num_peaks == 0) {
    *dir = 1;
    return;
  }
  *dir = 0;
  for(i=0;i<num_peaks;i++) {
    *dir += sd_histo(histogram, histo_bins, pv[i][0], pv[i][2]);
  }
  *dir /= 5; /* normalize by maximum value */
}

double *dir_histogram(Matrix image, int histo_bins, int threshold)
{
  int y,x, bin;
  double delh, delv;
  double delG, theta;
  double sum;
  double *histo;
  double dont_care;

  histo = Allocate(histo_bins, double);
  for(y=0;y<histo_bins;y++) histo[y]=0;
  sum = 0;
  dont_care = 0;
  for(y=1;y<image->height-1;y++) {
    for(x=1;x<image->width-1;x++) {
      delh = (image->data[y-1][x+1] + 
	      image->data[y][x+1] + 
	      image->data[y+1][x+1]) -
	     (image->data[y-1][x-1] + 
	      image->data[y][x-1] + 
	      image->data[y+1][x-1]);
      delv = (image->data[y-1][x-1] + 
	      image->data[y-1][x] + 
	      image->data[y-1][x+1]) -
	     (image->data[y+1][x-1] + 
	      image->data[y+1][x] + 
	      image->data[y+1][x+1]);
      delG = (fabs(delh) + fabs(delv)) / 2;
      if(delG >= threshold) {
	theta = atan2(delv, delh);
	if(theta < 0) theta += M_PI;
	else if(theta >= M_PI) theta -= M_PI;
	bin = (int)(theta * histo_bins / M_PI + 0.5);
	if(bin == histo_bins) bin = 0;
	if(bin < 0 || bin >= histo_bins) {
	  printf("Bin error: %d\n", bin);
	  exit(1);
	}
	histo[bin]++;
	sum++;
      }
      else {
	dont_care++;
      }
    }
  }
  /* Equitz correction: add dont_care to every bin */
  sum += dont_care;
  dont_care /= histo_bins;
  if(sum) {
    for(y=0;y<histo_bins;y++) histo[y] = (histo[y] + dont_care) / sum;
  }
  return histo;
}

double *new_dir_histogram(Matrix image, int histo_bins)
{
  int y,x, bin;
  double delh, delv;
  double delG, theta;
  double sum;
  double *histo;

  histo = Allocate(histo_bins, double);
  for(y=0;y<histo_bins;y++) histo[y]=0;
  sum = 0;
  for(y=1;y<image->height-1;y++) {
    for(x=1;x<image->width-1;x++) {
      delh = (image->data[y-1][x+1] + 
	      image->data[y][x+1] + 
	      image->data[y+1][x+1]) -
	     (image->data[y-1][x-1] + 
	      image->data[y][x-1] + 
	      image->data[y+1][x-1]);
      delv = (image->data[y-1][x-1] + 
	      image->data[y-1][x] + 
	      image->data[y-1][x+1]) -
	     (image->data[y+1][x-1] + 
	      image->data[y+1][x] + 
	      image->data[y+1][x+1]);
      delG = (fabs(delh) + fabs(delv))/2;
      if(!delv && !delh) {
	theta = 0.0;
      }
      else {
	theta = atan2(delv, delh);
      }
      if(theta < 0) theta += M_PI;
      else if(theta >= M_PI) theta -= M_PI;
      bin = (int)(theta * histo_bins / M_PI + 0.5);
      if(bin == histo_bins) bin = 0;
      if(bin < 0 || bin >= histo_bins) {
	printf("Bin error: %d\n", bin);
	exit(1);
      }
      histo[bin] += delG;
      sum += delG;
    }
  }
  if(sum) {
    for(y=0;y<histo_bins;y++) histo[y] /= sum;
  }
  return histo;
}

double Contrast(Matrix image)
{
  double var, norm4, mean, diff;
  int y,x;

  mean = MatrixMean(image);
  var = 0;
  norm4 = 0;
  for(y=0;y<image->height;y++) {
    for(x=0;x<image->width;x++) {
      diff = image->data[y][x] - mean;
      var += diff*diff;
      norm4 += diff*diff*diff*diff;
    }
  }
  var /= (image->width * image->height);
  norm4 /= (image->width * image->height);
  if(norm4 == 0) return 0;
  return var / pow(norm4,0.25) / 128; /* normalize by maximum value */
}

double Coarseness(Matrix image, int max_scale, double weight)
{
  int scale, step;
  Matrix temp, diff, best, scales = nullptr;
  double crs;
  int y,x;
  //char str[100];

  best = NULL;
  image = MatrixCopy(image);
  for(scale = 1; scale <= max_scale; scale++) {
    temp = image;
    image = nbr_avg(image, scale);
    MatrixFree(temp);
    diff = nbr_diff(image, scale);
    MatrixScale(diff, pow(weight, scale-1));
    if(!best) {
      best = diff;
      scales = MatrixCreate(diff->width, diff->height);
      MatrixFill(scales, 1<<scale);
    }
    else {
      step = 1<<(scale-1);
      temp = best;
      best = MatrixSub(best, step, step, diff->height+step, diff->width+step);
      MatrixFree(temp);
      temp = scales;
      scales = MatrixSub(scales, step, step, 
			 diff->height+step, diff->width+step);
      MatrixFree(temp);

      for(y=0;y < diff->height;y++) {
	for(x=0;x < diff->width;x++) {
	  if(best->data[y][x] <= diff->data[y][x]) {
	    best->data[y][x] = diff->data[y][x];
	    scales->data[y][x] = 1<<scale;
	  }
	}
      }
      MatrixFree(diff);
    }
  }

  MatrixFree(image);
  MatrixFree(best);
  crs = MatrixMean(scales);
  MatrixFree(scales);
  return crs / (1<<max_scale); /* normalize by maximum value */
}

/* Assumes image->width == image->height */
Matrix nbr_diff(Matrix image, int factor)
{
  int step = 1<<factor;
  Matrix result;
  int y,x;
  double h_diff, v_diff;

  result = MatrixCreate(image->height-step, image->width-step);
  for(y=0; y < image->height-step; y++) {
    for(x=0; x < image->width-step; x++) {
      h_diff = image->data[y][x] - image->data[y][x+step];
      v_diff = image->data[y][x] - image->data[y+step][x];
      result->data[y][x] = fmax(fabs(h_diff), fabs(v_diff));
    }
  }
  return result;
}

/* Assumes image->width == image->height */
Matrix nbr_avg(Matrix image, int factor)
{
  int step;
  Matrix result;
  int y,x;

  step = 1<<(factor-1);
  result = MatrixCreate(image->height-step, image->width-step);
  for(y=0; y < image->height-step; y++) {
    for(x=0; x < image->width-step; x++) {
      result->data[y][x] = (image->data[y][x] +
			    image->data[y][x+step] +
			    image->data[y+step][x] +
			    image->data[y+step][x+step]) / 4;
    }
  }
  return result;
}



//basic.cpp

/* Matrix implementation. Simple operations.
 * A matrix is a structure of {height, width, data}.
 * The data is double precision numbers, addressable both as a
 * two-dimensional array (m->data[row][col]) or a row-packed vector
 * (m->data[0][index]).
 */

//#include "matrix.h"

/* Prototypes ****************************************************************/
double dot_product(double *a, double *b, int len);
Matrix MatrixFilter(Matrix matrix, Matrix filter);
void MatrixSubtractRow(Matrix matrix, Real *row);
Matrix MatrixColSum(Matrix matrix);
Matrix MatrixColMean(Matrix matrix);
void MatrixMinMax(Matrix matrix, Real *minp, Real *maxp);
void MatrixRescale(Matrix matrix, Real new_min, Real new_max);
void MatrixTruncate(Matrix matrix, Real new_min, Real new_max);
void MatrixTruncRescale(Matrix matrix, Real trunc_min, Real trunc_max,
                        Real new_min, Real new_max);
void MatrixApplyFunc(Matrix matrix, MathFunc *func);
Matrix MatrixSub(Matrix matrix, int y1, int x1, int y2, int x2);
Real MatrixSum(Matrix matrix);
Real MatrixMean(Matrix matrix);
void MatrixAddValue(Matrix matrix, Real value);
void MatrixFill(Matrix matrix, Real value);
void MatrixScale(Matrix matrix, Real factor);

void MatrixReadImage(Matrix matrix, FILE *fp);
void MatrixWriteImage(Matrix matrix, FILE *fp);
void MatrixWriteVector(Matrix matrix, FILE *fp);
void MatrixWrite(Matrix matrix, FILE *fp);

Matrix MatrixFromArray(int height, int width, double *data);
Matrix MatrixCreate(int height, int width);
void MatrixFree(Matrix matrix);
Matrix MatrixCopy(Matrix matrix);

/* Functions *****************************************************************/

/* Vector dot product. */
double dot_product(double *a, double *b, int len)
{
    double r = 0.0;
    while(len--) r += *a++ * *b++;
    return r;
}

/* Apply a filter to a matrix.
 * The middle value of the filter is at
 *   ((filter->height+1)/2, (filter->width+1)/2).
 * Not very efficient; could have implemented some special cases (3x3).
 */
Matrix MatrixFilter(Matrix matrix, Matrix filter)
{
    int i,j,k,l;
    int halfy = filter->height/2, halfx = filter->width/2;
    Matrix out = MatrixCreate(matrix->height, matrix->width);
    for(i=0;i<matrix->height;i++) {
        for(j=0;j<matrix->width;j++) {
            Real sum = 0.0;
            for(k=0;k<filter->height;k++) {
                for(l=0;l<filter->width;l++) {
                    int ii = (i+k-halfy+matrix->height)%matrix->height,
                    ij = (j+l-halfx+matrix->width)%matrix->width;
                    sum += filter->data[k][l] * matrix->data[ii][ij];
                }
            }
            out->data[i][j] = sum;
        }
    }
    return out;
}

/* Vector difference */
static void subtract(Real *a, Real *b, int len)
{
    while(len--) *a++ -= *b++;
}

/* Subtract a vector from each row of a matrix */
void MatrixSubtractRow(Matrix matrix, Real *row)
{
    int i;
    
    for(i=0;i<matrix->height;i++) {
        subtract(matrix->data[i], row, matrix->width);
    }
}

/* Compute a row vector of column sums */
Matrix MatrixColSum(Matrix matrix)
{
    Matrix row;
    int i,j;
    
    row = MatrixCreate(1, matrix->width);
    MatrixFill(row, 0.0);
    for(i=0;i<matrix->height;i++) {
        for(j=0;j<matrix->width;j++) {
            row->data[0][j] += matrix->data[i][j];
        }
    }
    return row;
}

/* Compute a row vector of column means */
Matrix MatrixColMean(Matrix matrix)
{
    Matrix row;
    
    row = MatrixColSum(matrix);
    MatrixScale(row, (Real)1/matrix->height);
    return row;
}

/* Compute the maximum and minimum values in a matrix */
void MatrixMinMax(Matrix matrix, Real *minp, Real *maxp)
{
    int y,x;
    
    *minp = *maxp = matrix->data[0][0];
    for(y=0;y < matrix->height;y++) {
        for(x=0;x < matrix->width;x++) {
            if(matrix->data[y][x] < *minp) *minp = matrix->data[y][x];
            else if(matrix->data[y][x] > *maxp) *maxp = matrix->data[y][x];
        }
    }
}

/* Rescale the values in a matrix to range from new_min to new_max */
void MatrixRescale(Matrix matrix, Real new_min, Real new_max)
{
    Real old_min, old_max;
    
    MatrixMinMax(matrix, &old_min, &old_max);
    MatrixAddValue(matrix, -old_min);
    if(old_max != old_min)
        MatrixScale(matrix, (new_max - new_min) / (old_max - old_min));
    MatrixAddValue(matrix, new_min);
}

/* Truncate the values in a matrix,
 * i.e. set anything below new_min to new_min, and same for new_max.
 */
void MatrixTruncate(Matrix matrix, Real new_min, Real new_max)
{
    int y,x;
    
    for(y=0;y<matrix->height;y++) {
        for(x=0;x<matrix->width;x++) {
            if(matrix->data[y][x] < new_min) matrix->data[y][x] = new_min;
            else if(matrix->data[y][x] > new_max) matrix->data[y][x] = new_max;
        }
    }
}

/* Truncate the values in a matrix to [trunc_min, trunc_max] and then
 * rescale to range from new_min to new_max.
 */
void MatrixTruncRescale(Matrix matrix, Real trunc_min, Real trunc_max,
                        Real new_min, Real new_max)
{
    int i;
    double factor;
    
    factor = (new_max - new_min) / (trunc_max - trunc_min);
    for(i=0;i<matrix->height*matrix->width;i++) {
        if(matrix->data[0][i] <= trunc_min) matrix->data[0][i] = new_min;
        else if(matrix->data[0][i] >= trunc_max) matrix->data[0][i] = new_max;
        else
            matrix->data[0][i] = (matrix->data[0][i] - trunc_min) * factor + new_min;
    }
}

/* Apply a unary math function to every element */
void MatrixApplyFunc(Matrix matrix, MathFunc *func)
{
    int y,x;
    
    for(y=0;y<matrix->height;y++) {
        for(x=0;x<matrix->width;x++) {
            matrix->data[y][x] = func(matrix->data[y][x]);
        }
    }
}

/* Extract the submatrix from (y1, x1) to (y2, x2), inclusive */
Matrix MatrixSub(Matrix matrix, int y1, int x1, int y2, int x2)
{
    Matrix result;
    int range_x, range_y, y, x;
    
    if((x1 < 0) || (y1 < 0) || (y2 < y1) || (x2 < x1) ||
       (x2 >= matrix->width) || (y2 >= matrix->height)) {
        //fprintf(stderr, "MatrixSub: bad region (%d, %d)-(%d, %d) on %d by %d matrix\n", y1, x1, y2, x2, matrix->height, matrix->width);
    }
    
    range_x = x2 - x1 + 1;
    range_y = y2 - y1 + 1;
    result = MatrixCreate(range_y, range_x);
    for(y=0;y < range_y;y++) {
        for(x=0;x < range_x;x++) {
            result->data[y][x] = matrix->data[y+y1][x+x1];
        }
    }
    return result;
}

/* Sum of all elements */
Real MatrixSum(Matrix matrix)
{
    Real sum;
    int y,x;
    
    sum = 0;
    for(y=0;y < matrix->height;y++) {
        for(x=0;x < matrix->width;x++) {
            sum += matrix->data[y][x];
        }
    }
    return sum;
}

/* Mean of all elements */
Real MatrixMean(Matrix matrix)
{
    return MatrixSum(matrix) / matrix->width / matrix->height;
}

/* Add a constant value to all elements */
void MatrixAddValue(Matrix matrix, Real value)
{
    int y,x;
    
    for(y=0;y<matrix->height;y++) {
        for(x=0;x<matrix->width;x++) {
            matrix->data[y][x] += value;
        }
    }
}

/* Copy a matrix */
Matrix MatrixCopy(Matrix matrix)
{
    Matrix result;
    
    result = MatrixCreate(matrix->height, matrix->width);
    memcpy(result->data[0], matrix->data[0],
           matrix->width*matrix->height*sizeof(Real));
    return result;
}

/* Set all elements to a value */
void MatrixFill(Matrix matrix, Real value)
{
    int y,x;
    
    for(y=0;y<matrix->height;y++) {
        for(x=0;x<matrix->width;x++) {
            matrix->data[y][x] = value;
        }
    }
}

/* Multiply all elements by a factor */
void MatrixScale(Matrix matrix, Real factor)
{
    int y,x;
    
    for(y=0;y<matrix->height;y++) {
        for(x=0;x<matrix->width;x++) {
            matrix->data[y][x] *= factor;
        }
    }
}

/* Read a character file as a matrix */
void MatrixReadImage(Matrix matrix, FILE *fp)
{
    unsigned char *buffer;
    int i,bufsize;
    
    bufsize = matrix->width * matrix->height;
    buffer = Allocate(bufsize, unsigned char);
    if(!buffer) {
        fprintf(stderr, "Cannot allocate %lu byte input buffer\n", bufsize*sizeof(unsigned char));
        exit(1);
    }
    fread(buffer, sizeof(unsigned char), bufsize, fp);
    
    for(i=0;i<bufsize;i++) {
        matrix->data[0][i] = (Real)buffer[i];
    }
    free(buffer);
}

/* Write a matrix as a character file by clipping */
void MatrixWriteImage(Matrix matrix, FILE *fp)
{
    unsigned char *buffer;
    int i,bufsize;
    Real v;
    
    bufsize = matrix->width * matrix->height;
    buffer = Allocate(bufsize, unsigned char);
    if(!buffer) {
        fprintf(stderr, "Cannot allocate %lu byte output buffer\n", bufsize*sizeof(unsigned char));
        exit(1);
    }
    for(i=0;i<bufsize;i++) {
        v = matrix->data[0][i];
        if(v > 255) v = 255;
        else if(v < 0) v = 0;
        buffer[i] = (unsigned char)v;
    }
    
    fwrite(buffer, sizeof(unsigned char), bufsize, fp);
    free(buffer);
}

/* Write a matrix as a vector of floating-point numbers, one per line. */
void MatrixWriteVector(Matrix matrix, FILE *fp)
{
    int i;
    
    for(i=0;i<matrix->height*matrix->width;i++) {
        fprintf(fp, "%lg\n", matrix->data[0][i]);
    }
}

/* Write a matrix as a rectangle of numbers, one row per line */
void MatrixWrite(Matrix matrix, FILE *fp)
{
    int i,j;
    for(i=0;i<matrix->height;i++) {
        for(j=0;j<matrix->width;j++) {
            fprintf(fp, "%g ", matrix->data[i][j]);
        }
        fprintf(fp, "\n");
    }
}

/* Creates a Matrix whose data is given. The array becomes property
 * of the matrix; it will be freed when the matrix is. The array must
 * be height*width elements long.
 */
Matrix MatrixFromArray(int height, int width, double *data)
{
    Matrix result;
    int i;
    
    result = Allocate(1, struct MatrixStruct);
    result->width = width;
    result->height = height;
    result->data = Allocate(height, Real*);
    result->data[0] = data;
    for(i=1;i<height;i++) {
        result->data[i] = &result->data[0][i*width];
    }
    return result;
}

/* Create a new, uninitialized matrix */
Matrix MatrixCreate(int height, int width)
{
    double *data = Allocate(width*height, Real);
    if(!data) {
        fprintf(stderr, "Cannot allocate %d by %d matrix\n", height, width);
        return NULL;
    }
    return MatrixFromArray(height, width, data);
}

/* Destroy a matrix */
void MatrixFree(Matrix matrix)
{
    free(matrix->data[0]);
    free(matrix->data);
    free(matrix);
}


