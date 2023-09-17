#ifndef glcm_h
#define glcm_h

#include "FRMWRKHeader.h"

void quantizerMatrix(Mat image, Mat &imageQuantized, int grayscale);

void getCoOccurrenceMatrix(Mat image, Mat &coOccurence, int grayscale);

void normalizeCoOccurrenceMatrix(Mat coOccurence, Mat &coOccurenceNormalized, Mat imageQuantized, int grayscale);

void getFeatures(vector<double> &glcm_features, Mat coOccurenceNormalized, int grayscale);

#endif
