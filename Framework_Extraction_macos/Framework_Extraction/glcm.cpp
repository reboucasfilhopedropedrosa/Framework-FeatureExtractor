#include "glcm.hpp"


void quantizerMatrix(Mat image, Mat &imageQuantized, int grayscale){
    
    imageQuantized=Mat::zeros(image.size(), image.type());
    
    for (int i=0; i<image.rows; i++) {
        for (int j=0; j<image.cols; j++) {
            imageQuantized.at<uchar>(i,j)=image.at<uchar>(i,j)/grayscale;
        }
    }
}

void getCoOccurrenceMatrix(Mat image, Mat &coOccurence, int grayscale){
    coOccurence=Mat::zeros(grayscale, grayscale, CV_32FC1);
    
    for (int i=0; i<image.rows; i++) {
        for (int j=0; j<image.cols-1; j++) {
            coOccurence.at<float>(image.at<uchar>(i,j),image.at<uchar>(i,j+1))+=1;
        }
    }
}

void normalizeCoOccurrenceMatrix(Mat coOccurence, Mat &coOccurenceNormalized, Mat imageQuantized, int grayscale){
    coOccurenceNormalized=Mat::zeros(grayscale, grayscale, CV_32FC1);
    //Mat=Mat+Mat.t();
    //Mat=Mat/sum(Mat)[0];
    for (int i=0; i<grayscale; i++){
        for (int j=0; j<grayscale; j++) {
            coOccurenceNormalized.at<float>(i,j)=coOccurence.at<float>(i,j)/(imageQuantized.rows*(imageQuantized.cols-1));
        }
    }
}

void getFeatures(vector<double> &glcm_features, Mat coOccurenceNormalized, int grayscale){
    /*
     >> glcm_features[0]  -> Contrast
     >> glcm_features[1]  -> Correlation
     >> glcm_features[2]  -> Energy
     >> glcm_features[3]  -> Homogeneity
     >> glcm_features[4]  -> Inverse Difference Moment
     >> glcm_features[5]  -> Entropy
     >> glcm_features[6]  -> Dissimilarity
     >> glcm_features[7]  -> Angular Second Moment
     >> glcm_features[8]  -> Variance
     */
    
    int i,j;
    for (int i=0; i<9; i++) {
        glcm_features.push_back(0);
    }
    
    double correlation=0,mean1=0,mean2=0,deviation=0,deviation1=0,deviation2=0;
    for(i=0; i<grayscale; i++){
        for(j=0; j<grayscale; j++){
            glcm_features[0]+=( (i-j) * (i-j) * (coOccurenceNormalized.at<float>(i,j)) );// Contrast
            
            glcm_features[2]+= ( ((coOccurenceNormalized.at<float>(i,j))) * ((coOccurenceNormalized.at<float>(i,j))) );// Energy
            
            glcm_features[3]+= ( (coOccurenceNormalized.at<float>(i,j))) / ( 1+abs(i-j) );// Homogeneity
            
            if(i!=j){
                glcm_features[4]+= (coOccurenceNormalized.at<float>(i,j) / (1 + (i-j)*(i-j))) ;// Inverse Difference Moment
                glcm_features[6]+= (coOccurenceNormalized.at<float>(i,j)*(abs(i-j)) );// Dissimilarity
            }
            if(coOccurenceNormalized.at<float>(i,j)!=0){
                glcm_features[5]+=coOccurenceNormalized.at<float>(i,j)*log10(coOccurenceNormalized.at<float>(i,j));// Entropy
            }
            
            glcm_features[7]+=(coOccurenceNormalized.at<float>(i,j)*coOccurenceNormalized.at<float>(i,j));// Angular Second Moment
            
            
            //Operacoes para Correlacao
            correlation+= ( (i*j) * ( (coOccurenceNormalized.at<float>(i,j)) ) );//posicoes da matriz glcm
            mean1+= (i * ( (coOccurenceNormalized.at<float>(i,j)) ));//Media x
            mean2+= (j * ( (coOccurenceNormalized.at<float>(i,j)) ));//Media y
            deviation1+= ( (i*i) * (coOccurenceNormalized.at<float>(i,j)) );//Desvio padrao x
            deviation2+= ( (j*j) * (coOccurenceNormalized.at<float>(i,j)) );//Desvio padrao y
        }
    }
    glcm_features[5]*=-1;//entropy
    
    //operações da correlacao
    deviation1-=mean1*mean1;
    deviation2-=mean2*mean2;
    deviation1=sqrt(deviation1);
    deviation2=sqrt(deviation2);
    deviation=(deviation1*deviation2);
    if (deviation==0) {
        deviation=1;
    }
    glcm_features[1]=(correlation-(mean1*mean2))/deviation;// Correlation
    
    //Atributos que necessitam de outros já calculados
    for(i=0; i<grayscale; i++){
        for(j=0; j<grayscale; j++){
            glcm_features[8]+=((i-(mean1+mean2)/2)*(i-(mean1+mean2)/2))*coOccurenceNormalized.at<float>(i,j);//Variance
            
        }
    }
}
