#ifndef FRMWRKHeader_h
#define FRMWRKHeader_h

//#ifndef round
//#define round(a) ((int)a+0.5> ? (int)a : (int)a+1)
//#endif

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>
#include "opencv2/opencv.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits>
#include <iostream>
#include <string.h>
#include <fstream>
#include <vector>
#include "dcmtk/config/osconfig.h"
#include "dcmtk/dcmdata/dctk.h"
#include "dcmtk/dcmimgle/dcmimage.h"

using namespace cv;
using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// atributos da MIDE/SIM/SCM
#define aCOR 0  // correlation [-1 1]
#define aMDI 1  // momento de diferença inverso [0 1]
#define aENT 2  // entropy
#define aCHI 3  // Chi quadrado da diagonal da matriz [0 1]
#define aCSD 3  // Chi quadrado da diagonal da matriz [0 1]
// atributos da SCM
#define aCSR 4  // distancia chi-quadrado entre quadrantes
#define aMDR 5  // Mean absolute difference ratio
#define aDKL 6  // Kullback and Leibler divergence
#define aCAD 7  // complementary absolute difference
//#define aW1  5  // Wasserstein metric 1
//#define aGINI 6 // Gini coefficient [0 1]
//
#define MAXATT 8
//
#define MAXLIN 4096 // quantidade de linhas
#define MAXCOL 4096 // quantidade de colunas
//
#define MAXNIV 8 // quantidade de niveis de cinza - quantizacao

// filtros
#define TAMFILT    3
//
#define fAVERAGE   0
#define fLAPLACIAN 1
#define fGAUSSIAN  2
#define fSOBELHV   3
#define fSOBELH    31
#define fSOBELV    32
//
#define MAXFILT    4

// particionamento
#define pQUANT    0
#define MAXPART   1

// atributos da linha de comando
#define pSCMGLOBAL 0
#define pSCMMAP    1

// declaracoes
struct Matrix2 {
    float *fDat;
    int nRows, nCols;
};
//
struct Filter {
    struct Matrix2 M;
    float fSum;
};
//
struct Signal {
    struct Matrix2 M;
    int nLevels;
    int minLevel;
    int maxLevel;
};
//
struct scM {
    struct Matrix2 M;
    int nLevels;
    float fAtt[MAXATT]; // atributos da MIDE
};


// Common Functions
long double factorial(int num);

void DicomToMat(Mat &output, int UHmax, int UHmin, int type);
Mat strelDisk(int Radius);
void loadImages(Mat &image, char *address);
void normalizationTo8UC1(Mat image, Mat &output, int max, int min);
void imgToStructSignal(Mat img, struct Signal *s, int type);

vector<double> mom_hu_extraction(Mat image);
vector<double> mom_central_extraction(Mat image);
vector<double> mom_statistical_extraction(Mat image);
vector<double> glcm_extraction(Mat image, int grayscale);
vector<double> lbp_extraction(Mat image);
vector<double> ahtd_lung_extraction(Mat image);
vector<double> ahtd_brain_extraction(Mat image);
vector<double> scm_extraction(struct Signal *Sinal, int mode, clock_t Tmp[2]);
vector<double> tamura_extraction(Mat image);
vector<double> mom_zernick_extraction(Mat image);
vector<double> fourier_extraction(Mat image, int N);

double timeCounter(clock_t TmpEnd, clock_t TmpBegin);

void initializeAdrsNames();
void initializeAddresses(bool SO, vector<int> methods, char *DesktopChar);

void saveFeatures(vector<int> methods, vector<vector<double>> data, int label, char *DesktopChar);
void saveTimes(vector<int> methods, vector<double> time, char *DesktopChar);

void initializeMethods(Mat image8U, Mat image16S, vector<int> methods, int label, char *DesktopChar);

#endif
