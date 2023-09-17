#include "FRMWRKHeader.h"
#include "scm.hpp"
#include "lbp.hpp"
#include "histogram.hpp"
#include "glcm.hpp"
#include "tamura.h"

FILE *file;
char namefile[500];
vector<string> adrsNames;

long double factorial(int num){
    if (num<0) {
        cout<<"num must be a non-negative double."<<endl;
        exit(-1);
    }
    else if(num==0) {
        return 1;
    }
    else{
        int i;
        long double fat=1;
        for (i=2; i<=num; i++) {
            fat=fat*i;
        }
        return fat;
    }
}

void DicomToMat(Mat &output, int UHmax, int UHmin, int type){
    //type 1 => 8UC1
    //type 2 => 16SC1
    
    DicomImage image_DCM(namefile);
    const DiPixel *dmp = NULL;
    Sint16 *pixelV = NULL;
    dmp = image_DCM.getInterData();
    pixelV = (Sint16 *)dmp->getData();
    
    if (type==1) {
        output = Mat::zeros((int)image_DCM.getWidth(), (int)image_DCM.getHeight(), CV_8UC1);
        int i, j, pix;
        
        for(i=0; i<output.rows; i++){
            for (j=0; j<output.cols; j++){
                pix = pixelV[i*output.rows + j];
                if(pix>=UHmax)
                    pix = 255;
                else if(pix<=UHmin)
                    pix = 0;
                else
                    pix = (int)(255*(pix-UHmin)/(UHmax-UHmin));
                output.at<uchar>(i,j) = pix;
            }
        }
    }
    else if(type==2){
        output = Mat::zeros((int)image_DCM.getWidth(), (int)image_DCM.getHeight(), CV_16SC1);
        int i, j;
        for(i=0; i<output.rows; i++){
            for (j=0; j<output.cols; j++){
                output.at<Sint16>(i,j) = pixelV[i*output.rows + j];
            }
        }
    }
    else{
        cout<<"Error: type must be 1 or 2"<<endl;
        exit(-1);
    }
}

Mat strelDisk(int Radius){
    // THIS RETURNS STREL('DISK',RADIUS) LIKE IN MATLAB FOR RADIUS = ODD NUMBER BETWEEN 3-->21
    Mat sel((2*Radius-1),(2*Radius-1),CV_8U,cv::Scalar(255));
    int borderWidth = 0;
    switch (Radius){
        case 1: borderWidth = 0; break;
        case 3: borderWidth = 0; break;
        case 5: borderWidth = 2; break;
        case 7: borderWidth = 2; break;
        case 9: borderWidth = 4; break;
        case 11: borderWidth = 6; break;
        case 13: borderWidth = 6; break;
        case 15: borderWidth = 8; break;
        case 17: borderWidth = 8; break;
        case 19: borderWidth = 10; break;
        case 21: borderWidth = 10; break;
    }
    for (int i=0; i<borderWidth; i++){
        for (int j=0; j<borderWidth; j++){
            if (i+j<8){
                sel.at<uchar>(i,j)=0;
                sel.at<uchar>(i,sel.cols-1-j)=0;
                sel.at<uchar>(sel.rows-1-i,j)=0;
                sel.at<uchar>(sel.rows-1-i,sel.cols-1-j)=0;
            }
        }
    }
    return sel;
}
void loadImages(Mat &image, char *address){
    char ff, *ffNum = nullptr;
    int Num;
    
    vector<vector<int>> lines;
    vector<int> line;
    vector<char> number;
    
    file=fopen(address, "r");
    
    while (!feof(file)) {
        ff=fgetc(file);
        if (ff=='\r'){
            lines.push_back(line);
            line.clear();
        }
        else{
            if (ff!=' ') {
                number.push_back(ff);
            }
            else{
                ffNum=(char*)calloc(sizeof(char*), number.size());
                for (int i=0; i<number.size(); i++) {
                    ffNum[i]=number[i];
                }
                Num=atoi(ffNum);
                line.push_back(Num);
                free(ffNum);
                number.clear();
            }
        }
    }
    line.clear();
    fclose(file);
    
    //Pass values to Mat
    image=Mat::zeros((int)lines.size(), (int)lines[0].size(), CV_32FC1);
    for (int i=0; i<image.rows; i++) {
        for (int j=0; j<image.cols; j++) {
            image.at<float>(j,i)=lines[i][j];
        }
    }
    lines.clear();
}
void normalizationTo8UC1(Mat image, Mat &output, int max, int min){
    Mat img;
    image.copyTo(img);
    float value;
    //Normalization
    //Lung application
    //max = 100;
    //min = -1000;
    for (int i=0; i<img.rows; i++) {
        for (int j=0; j<img.cols; j++) {
            img.at<float>(i,j)=(img.at<float>(i,j)-min)*255/(max-min);
            value=img.at<float>(i,j);
            if      (value>255)                    { img.at<float>(i,j) = 255; }
            else if (value<0)                      { img.at<float>(i,j) = 0;   }
            if      ( (value - (int)value) > 0.5  ){ img.at<float>(i,j)=(int)value+1; }
            else if ( (value - (int)value) <= 0.5 ){ img.at<float>(i,j)=(int)value;   }
        }
    }
    
    img.convertTo(output, CV_8UC1);
    img.release();
}
void imgToStructSignal(Mat img, struct Signal *s, int type){

    (*s).nLevels=0;
    (*s).minLevel=999;
    (*s).maxLevel=0;
    // carrega tamanho do sinal
    (*s).M.nRows=img.cols;
    (*s).M.nCols=img.rows;
    
    if ((*s).M.nRows>MAXLIN){(*s).M.nRows=MAXLIN;}
    if ((*s).M.nCols>MAXCOL){(*s).M.nCols=MAXCOL;}
    
    // aloca memoria do sinal
    (*s).M.fDat = (float*) malloc((*s).M.nRows*(*s).M.nCols*sizeof(float));
    
    // carrega sinal
    int sample = 0;
    for (int i=0; i<img.cols; i++) {
        for (int j=0; j<img.rows; j++) {
            
            if      (type==1)   { sample=img.at<uchar>(j,i); }
			else if (type==2)   { sample=img.at<float>(j,i); }
            //else if (type==2)   { sample=img.at<Sint16>(j,i); }
            
            (*s).M.fDat[i + j * (*s).M.nCols] = (float)sample;
            if (sample>(*s).nLevels)  { (*s).nLevels=sample;  }
            if (sample<(*s).minLevel) { (*s).minLevel=sample; }
            if (sample>(*s).maxLevel) { (*s).maxLevel=sample; }
        }
    }
    
    // quantizacao para a SCM (min-max)
    for (int i=0 ; i<(*s).M.nRows; i++){
        for (int j=0 ; j<(*s).M.nCols ; j++){
            (*s).M.fDat[i * (*s).M.nCols + j]= round(
                                                     ( (*s).M.fDat[i * (*s).M.nCols + j] - (float)(*s).minLevel )  *
                                                     (float)(MAXNIV-1) / ( (float)(*s).maxLevel - (float)(*s).minLevel )
                                                     );
        }
    }
    (*s).nLevels=MAXNIV;
}

vector<double> mom_hu_extraction(Mat image){
    vector<double> hu_features;
    
    double Hu[7];
    Moments momentos =  moments(image,false);
    HuMoments(momentos, Hu);
    
    for (int i=0; i<7; i++) {
        hu_features.push_back(1000000000000000*Hu[i]);
    }
    
    return hu_features;
}
vector<double> mom_central_extraction(Mat image){
    vector<double> central_features;
    
    double Central[7];
    Moments momentos =  moments(image, false);
    Central[0]=momentos.mu20;
    Central[1]=momentos.mu11;
    Central[2]=momentos.mu02;
    Central[3]=momentos.mu30;
    Central[4]=momentos.mu21;
    Central[5]=momentos.mu12;
    Central[6]=momentos.mu03;
    
    for (int i=0; i<7; i++) {
        central_features.push_back(Central[i]);
    }
    
    return central_features;
    
}
vector<double> mom_statistical_extraction(Mat image){
    vector<double> statistical_features;
    
    double Stat[10];
    Moments momentos =  moments(image, false);
    Stat[0]=momentos.m00;
    Stat[1]=momentos.m10;
    Stat[2]=momentos.m01;
    Stat[3]=momentos.m20;
    Stat[4]=momentos.m11;
    Stat[5]=momentos.m02;
    Stat[6]=momentos.m30;
    Stat[7]=momentos.m21;
    Stat[8]=momentos.m12;
    Stat[9]=momentos.m03;
    
    for (int i=0; i<10; i++) {
        statistical_features.push_back(Stat[i]);
    }
    
    return statistical_features;
    
}
vector<double> glcm_extraction(Mat image, int grayscale){
    vector<double> glcm_features;
    
    Mat imageQuantized, coOccurence, coOccurenceNormalized;
    image.copyTo(imageQuantized);
    
    /* Quantiza a imagem de entrada */
    //quantizerMatrix(image, imageQuantized, grayscale);
    
    /* Calcula a Matriz de Co-Ocorrencia */
    getCoOccurrenceMatrix(imageQuantized, coOccurence, grayscale);
    
    /* Normaliza a Matriz de Co-Ocorrencia */
    normalizeCoOccurrenceMatrix(coOccurence, coOccurenceNormalized, imageQuantized, grayscale);
    
    /* Calcula os atributos da GLCM */
    getFeatures(glcm_features, coOccurenceNormalized, grayscale);
    
    imageQuantized.release();
    coOccurence.release();
    coOccurenceNormalized.release();
    
    return glcm_features;
}
vector<double> lbp_extraction(Mat image){
    vector<double> lbp_features;
    Mat lbp_img = lbp::OLBP(image);
    int quant_tomCinza = 12;//Quantidade de Tons de Cinza
    normalize(lbp_img, lbp_img, 0, quant_tomCinza-1, NORM_MINMAX, CV_8UC1);
    //Parametros da Função Normalize
    //1 - imagem de entrada // 2 - imagem de saida
    //3, 4 - Variacao dos tons de cinza que a imagem vai ter. Neste caso -> 0 ~ 11
    //5 - A funcao NORM_MINMAX serve para normalizar entre os valores do parametros (3,4)
    //6 - CV_8UC1 -> Imagem 8 bits Uchar, 1 canal
    Mat spatial_hist = lbp::spatial_histogram(lbp_img, quant_tomCinza, 3, 3, 0);
    // Parametros da Funcao spatial_hist
    //1 - imagem de entrada // 2 - A variacao entre os parametros (3,4) da funcao Normalize no caso 12
    //3, 4 - Tamanho da mascara(vizinhanca)  3x3
    //
    
    for(int j = 0; j < spatial_hist.cols; ++j){
        if(spatial_hist.at<int>(0, j) != -1){
            lbp_features.push_back(spatial_hist.at<int>(0, j));
        }
    }
    return lbp_features;
}
vector<double> ahtd_lung_extraction(Mat image){
    vector<double> ahtd_lung_features;
    
    double osso=0;
    double enf=0;
    double normal=0;
    double poucoAerada=0;
    double naoAerada=0;
    
    for(int i = 0; i < image.cols; i++){
        for(int j = 0; j < image.rows; j++){
            
            //regiao considerada com enfisema
            //if(image.at<Sint16>(i, j) <-901 && image.at<Sint16>(i, j)>-1500)
			if(image.at<float>(i, j) <-901 && image.at<float>(i, j)>-1500)
                enf++;
            //regiao considerada normal
            //else if(image.at<Sint16>(i, j)>-900 && image.at<Sint16>(i, j)<-501)
			else if(image.at<float>(i, j)>-900 && image.at<float>(i, j)<-501)
                normal++;
            //regiao pouco aerada
            //else if(image.at<Sint16>(i, j)>-500 && image.at<Sint16>(i, j)<-101)
			else if(image.at<float>(i, j)>-500 && image.at<float>(i, j)<-101)
                poucoAerada++;
            //regiao considerada nao aerada(fora do pulmao)
            //else if(image.at<Sint16>(i, j)>-100 &&image.at<Sint16>(i, j)<+101)
			else if(image.at<float>(i, j)>-100 &&image.at<float>(i, j)<+101)
                naoAerada++;
            //else if(image.at<Sint16>(i, j)>600 && image.at<Sint16>(i, j)<2000)
			else if(image.at<float>(i, j)>600 && image.at<float>(i, j)<2000)
                osso++;
        }
    }
    
    double quant=enf+normal+poucoAerada+naoAerada+osso;
    ahtd_lung_features.push_back(enf/quant);
    ahtd_lung_features.push_back(normal/quant);
    ahtd_lung_features.push_back(poucoAerada/quant);
    ahtd_lung_features.push_back(naoAerada/quant);
    ahtd_lung_features.push_back(osso/quant);
    
    return ahtd_lung_features;
}
vector<double> ahtd_brain_extraction(Mat image){
    vector<double> ahtd_brain_features;
    
    double subBranca=0;
    double subCizenta=0;
    double liquido=0;
    double tecidoMassaOsso=0;
    double osso=0;
    
    for(int i = 0; i < image.cols; i++){
        for(int j = 0; j < image.rows; j++){
            
            //regiao Branca
            //if(image.at<Sint16>(i, j) > 20 && image.at<Sint16>(i, j) < 35)
			if(image.at<float>(i, j) > 20 && image.at<float>(i, j) < 35)
                subBranca++;
            //regiao Cizenta
            //else if(image.at<Sint16>(i, j)>34 && image.at<Sint16>(i, j)<50)
			else if(image.at<float>(i, j)>34 && image.at<float>(i, j)<50)
                subCizenta++;
            //regiao liquido
            //else if(image.at<Sint16>(i, j)>-100 && image.at<Sint16>(i, j)<20)
			else if(image.at<float>(i, j)>-100 && image.at<float>(i, j)<20)
                liquido++;
            //regiao tecidoMassaOsso
            //else if(image.at<Sint16>(i, j)>49 &&image.at<Sint16>(i, j)<201)
			else if(image.at<float>(i, j)>49 &&image.at<float>(i, j)<201)
                tecidoMassaOsso++;
            //regiao osso
            //else if(image.at<Sint16>(i, j)>200 && image.at<Sint16>(i, j)<2000)
			else if(image.at<float>(i, j)>200 && image.at<float>(i, j)<2000)
                osso++;
        }
    }
    
    double quant=subBranca+subCizenta+liquido+tecidoMassaOsso+osso;
    ahtd_brain_features.push_back(subBranca/quant);
    ahtd_brain_features.push_back(subCizenta/quant);
    ahtd_brain_features.push_back(liquido/quant);
    ahtd_brain_features.push_back(tecidoMassaOsso/quant);
    ahtd_brain_features.push_back(osso/quant);
    
    return ahtd_brain_features;
}
vector<double> scm_extraction(struct Signal *Sinal, int mode, clock_t Tmp[2]){
    
    vector<double> mide_features;
    int par_filtro=0;
    
    struct Signal *s; // original signal
    struct Signal *f; // reference signal
    struct Signal *f2 = nullptr; // filtered signal
    struct scM *scmglob; // SCM global
    struct Filter *filt;
    
    if      (mode == 0) { par_filtro = fAVERAGE;   }
    else if (mode == 1) { par_filtro = fLAPLACIAN; }
    else if (mode == 2) { par_filtro = fGAUSSIAN;  }
    else if (mode == 3) { par_filtro = fSOBELHV;   }
    else if (mode == 4) { par_filtro = fSOBELH;    }
    else if (mode == 5) { par_filtro = fSOBELV;    }
    else                { par_filtro = fAVERAGE;   }
    
    
    // aloca memoria dos sinais originais e filtrados
    s=(struct Signal*)malloc(1*sizeof(struct Signal));
    f=(struct Signal*)malloc(1*sizeof(struct Signal));
    
    //Copying Struct Signal Sinal to Struct Signal s
    gsignal_copy((struct Signal *)Sinal,(struct Signal *)s);
    
    // aloca memoria do filtro e do sinal filtrado
    filt=(struct Filter*)malloc(1*sizeof(struct Filter));
    (*filt).M.fDat = (float*) malloc(TAMFILT*TAMFILT*sizeof(float));
    (*filt).M.nRows=TAMFILT;
    (*filt).M.nCols=TAMFILT;
    if (par_filtro==fSOBELHV){
        // aloca memoria do segundo sinal filtrado para SOBELHV
        f2=(struct Signal*)malloc(1*sizeof(struct Signal));
        //gsignal_copy((struct Signal *)s,(struct Signal *)f);  // copia
        //gsignal_copy((struct Signal *)s,(struct Signal *)f2); // copia
        gcreate_filter((struct Filter*) filt, fSOBELH);       // sobel horizontal
        gsignal_filt((struct Signal *)s,(struct Signal *)f,(struct Filter *)filt);
        gcreate_filter((struct Filter*) filt, fSOBELV);       // sobel vertical
        gsignal_filt((struct Signal *)s,(struct Signal *)f2,(struct Filter *)filt);
        gsignal_max((struct Signal *)f,(struct Signal *)f2);  // sobel HV
    }else{
        gcreate_filter((struct Filter*) filt, par_filtro);
        gsignal_filt((struct Signal *)s,(struct Signal *)f,(struct Filter *)filt);
    }
    
    // calcula SCM global e atributos
    scmglob=(struct scM*)malloc(1*sizeof(struct scM));
    (*scmglob).M.fDat = (float*) malloc((*s).nLevels*(*s).nLevels*sizeof(float));
    
    Tmp[0] = clock();
    gscm((struct Signal *)s,(struct Signal *)f,(struct scM *)scmglob);
    Tmp[1] = clock();
    
    // Passando valores para vector mide_features
    mide_features.push_back((*scmglob).fAtt[aCOR]);
    mide_features.push_back((*scmglob).fAtt[aMDI]);
    mide_features.push_back((*scmglob).fAtt[aENT]);
    mide_features.push_back((*scmglob).fAtt[aCHI]);
    mide_features.push_back((*scmglob).fAtt[aCSR]);
    mide_features.push_back((*scmglob).fAtt[aMDR]);
    mide_features.push_back((*scmglob).fAtt[aDKL]);
    mide_features.push_back((*scmglob).fAtt[aCAD]);
    
    // libera memoria alocada
    free((*s).M.fDat);
    free(s);
    //
    free((*f).M.fDat);
    free(f);
    //
    free((*filt).M.fDat);
    free(filt);
    //
    if (par_filtro==fSOBELHV){
        free((*f2).M.fDat);
        free(f2);
    }
    
    //pSCMGLOBAL
    free((*scmglob).M.fDat);
    free(scmglob);
    
    return mide_features;
}
vector<double> tamura_extraction(Mat image){
    vector<double> tamura_features;
    
    Matrix M;
    M = MatrixCreate(image.rows, image.cols);
    
    int i,j;
    for (i=0; i<image.rows; i++) {
        for (j=0; j<image.cols; j++) {
            M->data[i][j] = (Real)image.at<uchar>(i,j);
        }
    }
    
    int orientations;
    double dir;
    tamura_features.push_back(Coarseness(M, 5, 1));
    tamura_features.push_back(Contrast(M));
    Directionality(M, 16, 12, 2, &orientations, &dir);
    tamura_features.push_back(dir);
    MatrixFree(M);
    
    return tamura_features;
}
vector<double> mom_zernick_extraction(Mat image){
    vector<double> zernick_features;
    int i, j;
    int i2, j2;
    int ord = 11, rep;
    int N = 0;
    Mat imageN;
    if (image.rows!=image.cols) {
        if (image.rows<image.cols) {
            N = (int)image.rows;
            resize(image, imageN, Size(N, N));
        }else{
            N = (int)image.cols;
            resize(image, imageN, Size(N, N));
        }
    }else{
        N = (int)image.rows;
        image.copyTo(imageN);
    }
    
    //*****Converter para imagem quadrada!
    
    //logical NOT on image
    threshold(imageN, imageN, 0, 255, CV_THRESH_BINARY_INV | CV_THRESH_OTSU);
    
    //get X and Y => [X,Y] = meshgrid(1:N,1:N)
    Mat X = Mat::zeros(N, N, CV_32FC1);
    Mat Y = Mat::zeros(N, N, CV_32FC1);
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            X.at<float>(i,j)= j+1;
            Y.at<float>(i,j)= i+1;
        }
    }
    //get R = sqrt((2.*X-N-1).^2+(2.*Y-N-1).^2)/N
    Mat R = Mat::zeros(N, N, CV_32FC1);
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            R.at<float>(i,j)= sqrt(pow((2*X.at<float>(i,j)-N-1), 2)+pow((2*Y.at<float>(i,j)-N-1), 2))/N;
        }
    }
    
    //get Theta = atan2((N-1-2.*Y+2),(2.*X-N+1-2))
    Mat Theta = Mat::zeros(N, N, CV_32FC1);
    
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            Theta.at<float>(i,j) =atan2( (N-1-2*Y.at<float>(i,j)+2) , (2*X.at<float>(i,j)-N+1-2) );
        }
    }
    //get R = (R<=1).*R
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            if (R.at<float>(i,j)>=1) {
                R.at<float>(i,j)=0;
            }
        }
    }
    //Zernike Moments
    Mat Rad;
    long double c=0.0;
    long int cnt=0;
    std::complex<double> c_i(0, -1);
    vector<std::complex<double>> Product;
    std::complex<double> Z(0, 0);
    
//    cout<<T.real()<<endl;
//    cout<<T.imag()<<endl;
    for (rep=1; rep<=ord; rep=rep+2) {
        //Radialpoly (Function)
        Rad = Mat::zeros(N, N, CV_32FC1);
        for (i=0; i<=(ord-abs(rep))/2; i++) {
            //do the operations
            //(n=ord and m=rep and s=i)
            //c = (-1)^s*factorial(n-s)/(factorial(s)*factorial((n+abs(m))/2-s)*...
            //    factorial((n-abs(m))/2-s));
            
            c = pow(-1, i)*factorial(ord-i);
            c = c / (factorial(i) * factorial(((ord+abs(rep))/2)-i) * factorial(((ord-abs(rep))/2)-i));
            
            //(rad=Rad and r=R and n=ord and s=i)
            //rad = rad + c*r.^(n-2*s);
            
            for (i2=0; i2<N; i2++) {
                for (j2=0; j2<N; j2++) {
                    Rad.at<float>(i2,j2) += c*pow(R.at<float>(i2,j2), ord-2*i);
                }
            }
        }
        //After, we calculate the moments
        //Product = p(x,y).*Rad.*exp(-1i*m*Theta);
        for (i=0; i<N; i++) {
            for (j=0; j<N; j++) {
                Product.push_back((double)(imageN.at<uchar>(i,j))*(double)Rad.at<float>(i,j)*exp(c_i*(double)rep*(double)Theta.at<float>(i,j)));
            }
        }
        
        //Z = sum(Product(:));        % calculate the moments
        Z=0;
        for (i=0; i<Product.size(); i++) {
            Z+=Product[i];
        }
        
        //cnt = nnz(R)+1;             % count the number of pixels inside the unit circle
        cnt=1;
        for (i=0; i<N; i++) {
            for (j=0; j<N; j++) {
                if (R.at<float>(i,j)!=0) {
                    cnt++;
                }
            }
        }
        
        //Z = (n+1)*Z/cnt;            % normalize the amplitude of moments
        Z = ((double)(ord+1)*Z)/(double)cnt;
        //A = abs(Z);                 % calculate the amplitude of the moment
        zernick_features.push_back(std::abs(Z));
        //Phi = angle(Z)*180/pi;      % calculate the phase of the mement (in degrees) (Nao precisa)
        
        
        Rad.release();
        Product.clear();
    }
    X.release();
    Y.release();
    R.release();
    Theta.release();
    imageN.release();
    
    return zernick_features;
}
vector<double> fourier_extraction(Mat image, int N){
    vector<double> fourier_features;
    
    int i,j;
    Mat I;
    //im2bw
    threshold(image, I, 127, 255, CV_THRESH_BINARY | CV_THRESH_OTSU);
    //Opening with the Disk structure element (5 = radius)
    Mat strel = strelDisk(5);
    erode(I, I, strel);
    dilate(I, I, strel);
    
    //get the contour
    Mat findI;
    I.copyTo(findI);
    
    vector<vector<Point>> contours;
    vector<Point> bestContour;
    findContours(findI, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE);
    
    if (contours.size()==0) {
        for (i=0; i<28; i++) {
            fourier_features.push_back(0);
        }
        return fourier_features;
    }
    
    int long max=0;
    int idc=0;
    for (i=0; i<contours.size(); i++) {
        if (max<(long int)contours[i].size()) {
            max=(long int)contours[i].size();
            idc=i;
        }
    }
    bestContour=contours[idc];
    contours.clear();
    
//    cout<<bestContour.size()<<endl;
//    Mat test=Mat::zeros(image.size(), CV_8UC1);
//    for (int j=0; j<bestContour.size(); j++) {
//        test.at<uchar>(bestContour[j])=255;
//    }
//    imshow("test", test); waitKey();
    
    vector<float> dt;
    float i2=0, deltaX=0, deltaY=0;
    for (i=0; i<bestContour.size(); i++) {
        i2 = abs(div(i-2, (int)bestContour.size()).rem + 1);
        deltaX = bestContour[i].x - bestContour[i2].x;
        deltaY = bestContour[i].y - bestContour[i2].y;
        dt.push_back(sqrt(pow(deltaX, 2)+pow(deltaY, 2)));
    }
    
    float T=0;
    for (i=0; i<dt.size(); i++) {
        T+=dt[i];
    }
    vector<float> ai;
    vector<float> bi;
    vector<float> ci;
    vector<float> di;
    int coef_i;
    for (j=1; j<=N; j++) {
        ai.push_back(0);
        bi.push_back(0);
        ci.push_back(0);
        di.push_back(0);
        //coef a
        deltaX=0;
        deltaY=0;
        for (coef_i=0; coef_i<bestContour.size(); coef_i++) {
            i2 = abs(div(i-2, (int)bestContour.size()).rem + 1);
            deltaX= bestContour[coef_i].x - bestContour[i2].x;
            deltaY+= dt[coef_i];//aux
            ai[j-1] += deltaX/dt[coef_i] * (cos(2*j*M_PI*deltaY/T) - cos(2*j*M_PI*(deltaY-dt[coef_i])/T));
        }
        ai[j-1] = ai[j-1]*T/(2*pow(j, 2)*pow(M_PI, 2));
        //coef b
        deltaX=0;
        deltaY=0;
        for (coef_i=0; coef_i<bestContour.size(); coef_i++) {
            i2 = abs(div(i-2, (int)bestContour.size()).rem + 1);
            deltaX= bestContour[coef_i].x - bestContour[i2].x;
            deltaY+= dt[coef_i];//aux
            bi[j-1] += deltaX/dt[coef_i] * (sin(2*j*M_PI*deltaY/T) - sin(2*j*M_PI*(deltaY-dt[coef_i])/T));
        }
        bi[j-1] = bi[j-1]*T/(2*pow(j, 2)*pow(M_PI, 2));
        //coef c
        deltaX=0;
        deltaY=0;
        for (coef_i=0; coef_i<bestContour.size(); coef_i++) {
            i2 = abs(div(i-2, (int)bestContour.size()).rem + 1);
            deltaY= bestContour[coef_i].y - bestContour[i2].y;
            deltaX+= dt[coef_i];//aux
            ci[j-1] += deltaX/dt[coef_i] * (cos(2*j*M_PI*deltaX/T) - cos(2*j*M_PI*(deltaX-dt[coef_i])/T));
        }
        ci[j-1] = ci[j-1]*T/(2*pow(j, 2)*pow(M_PI, 2));
        //coef d
        deltaX=0;
        deltaY=0;
        for (coef_i=0; coef_i<bestContour.size(); coef_i++) {
            i2 = abs(div(i-2, (int)bestContour.size()).rem + 1);
            deltaY= bestContour[coef_i].y - bestContour[i2].y;
            deltaX+= dt[coef_i];//aux
            di[j-1] += deltaX/dt[coef_i] * (sin(2*j*M_PI*deltaX/T) - sin(2*j*M_PI*(deltaX-dt[coef_i])/T));
        }
        di[j-1] = di[j-1]*T/(2*pow(j, 2)*pow(M_PI, 2));
        //
    }
    //Normalize the coefficients (Do it later)
//    % Normalize the coefficients
//    % --------------------------
//    
//    % 1. Phase invariance
//    % This is needed to receive the same fourier descriptor no matter which
//    % order of points you pick
//    [a,b,c,d] = normalizePhase(a,b,c,d,N);
//    
//    % Extract the scale-invariance to use it later
//    E = sqrt(a(1)^2+c(1)^2);
//    
//    % 2. Rotation invariance
//    if rotation
//        [a,b,c,d] = normalizeRotation(a,b,c,d,N);
//    end
//    
//    % 3. Scale invariance
//    a = a./E;
//    b = b./E;
//    c = c./E;
//    d = d./E;
    
    for (i=0; i<ai.size(); i++) {
        fourier_features.push_back(ai[i]);
    }
    for (i=0; i<bi.size(); i++) {
        fourier_features.push_back(bi[i]);
    }
    for (i=0; i<ci.size(); i++) {
        fourier_features.push_back(ci[i]);
    }
    for (i=0; i<di.size(); i++) {
        fourier_features.push_back(di[i]);
    }
    
    I.release();
    findI.release();
    contours.clear();
    bestContour.clear();
    dt.clear();
    ai.clear();
    bi.clear();
    ci.clear();
    di.clear();
    return fourier_features;
}

double timeCounter(clock_t TmpEnd, clock_t TmpBegin){
    double time;
    time = (TmpEnd-TmpBegin)*1000.0/CLOCKS_PER_SEC;
    return time;
}

void initializeAdrsNames(){
    adrsNames.push_back("MOMENTS_HU");
    adrsNames.push_back("MOMENTS_CENT");
    adrsNames.push_back("MOMENTS_STAT");
    adrsNames.push_back("GLCM");
    adrsNames.push_back("LBP");
    adrsNames.push_back("AHTD_LUNG");
    adrsNames.push_back("AHTD_BRAIN");
    adrsNames.push_back("SCM_AVERAGE");
    adrsNames.push_back("SCM_LAPLACIAN");
    adrsNames.push_back("SCM_GAUSSIAN");
    adrsNames.push_back("SCM_SOBELHV");
    adrsNames.push_back("SCM_SOBELH");
    adrsNames.push_back("SCM_SOBELV");
    adrsNames.push_back("TAMURA");
    adrsNames.push_back("ZERNICK");
    adrsNames.push_back("FOURIER");
}
void initializeAddresses(bool SO, vector<int> methods, char *DesktopChar){
    initializeAdrsNames();

    char namefile[500];
    char folders[500];

	//Create Main Folder
	if (SO==true)
		sprintf(folders, "mkdir %s", DesktopChar);
	else
		sprintf(folders, "md %s", DesktopChar);
    system(folders);

    for (int i=0; i<methods.size(); i++) {
        if (methods[i]==1 || methods[i]==2) {
			//Create Sub Folders
			if (SO==true)
				sprintf(folders, "mkdir %s/DB_%s",DesktopChar, adrsNames[i].c_str());
			else
				sprintf(folders, "md %s\\DB_%s",DesktopChar, adrsNames[i].c_str());
            system(folders);

			//Create DB archives
			if (SO==true)
				sprintf(namefile,"%s/DB_%s/DB_%s.txt", DesktopChar, adrsNames[i].c_str(), adrsNames[i].c_str());
			else
				sprintf(namefile,"%s\\DB_%s\\DB_%s.txt", DesktopChar, adrsNames[i].c_str(), adrsNames[i].c_str());
            file=fopen(namefile,"w");
            fclose(file);

			//Create Time archives
			if (SO==true)
				sprintf(namefile,"%s/Time_DB_%s.txt", DesktopChar, adrsNames[i].c_str());
			else
				sprintf(namefile,"%s\\Time_DB_%s.txt", DesktopChar, adrsNames[i].c_str());
            file=fopen(namefile,"w");
            fclose(file);
        }
    }
    
}

void saveFeatures(vector<int> methods, vector<vector<double>> data, int label, char *DesktopChar){
    int counter=0;
    for (int i=0; i<methods.size(); i++) {
        if (methods[i]==1 || methods[i]==2) {
            sprintf(namefile, "%s/DB_%s/DB_%s.txt", DesktopChar, adrsNames[i].c_str(), adrsNames[i].c_str());
            file=fopen(namefile, "a");
            for (int j=0; j<data[counter].size(); j++) {
                fprintf(file, "%lf,", data[counter][j]);
                //cout<<data[counter][j]<<", ";
            }
            //cout<<endl;
            fprintf(file, "%d\n",label);
            fclose(file);
            
            counter++;
        }
    }
    
    fclose(file);
}
void saveTimes(vector<int> methods, vector<double> time, char *DesktopChar){
    int counter=0;
    for (int i=0; i<methods.size(); i++) {
        if (methods[i]==1 || methods[i]==2) {
            sprintf(namefile,"%s/Time_DB_%s.txt", DesktopChar, adrsNames[i].c_str());
            file=fopen(namefile,"a");
            fprintf(file, "%lf\n",time[counter]);
            fclose(file);
            counter++;
        }
    }
    fclose(file);
}

void initializeMethods(Mat image8U, Mat image16S, vector<int> methods, int label, char *DesktopChar){
    Mat imageMethods;
    vector<vector<double>> allDatas;
    vector<double> data;
    vector<double> times;
    clock_t Tmp[2];
    
    for (int i=0; i<methods.size(); i++) {
        switch (i) {
            //*************** HU MOMENTS ***************
            case 0:
                if (methods[i]==1 || methods[i]==2) {
                    if      (methods[i]==1) { image8U.copyTo(imageMethods);  }
                    else if (methods[i]==2) { cout<<"The Hu moment method cannot use CV_16S image type"<<endl; exit(-1); }
                    //*** BEGIN
                    Tmp[0]=clock();
                    data = mom_hu_extraction(imageMethods);
                    Tmp[1]=clock();
                    //*** END
                    times.push_back(timeCounter(Tmp[1], Tmp[0]));
                    allDatas.push_back(data);
                    data.clear();
                    imageMethods.release();
                }else{ break; }
                break;
            //*************** CENTRAL MOMENTS ***************
            case 1:
                if (methods[i]==1 || methods[i]==2) {
                    if      (methods[i]==1) { image8U.copyTo(imageMethods);  }
                    else if (methods[i]==2) { cout<<"The Central moment method cannot use CV_16S image type"<<endl; exit(-1); }
                    //*** BEGIN
                    Tmp[0]=clock();
                    data = mom_central_extraction(imageMethods);
                    Tmp[1]=clock();
                    //*** END
                    times.push_back(timeCounter(Tmp[1], Tmp[0]));
                    allDatas.push_back(data);
                    data.clear();
                    imageMethods.release();
                }else{ break; }
                break;
            //*************** STATISTICAL MOMENTS ***************
            case 2:
                if (methods[i]==1 || methods[i]==2) {
                    if      (methods[i]==1) { image8U.copyTo(imageMethods);  }
                    else if (methods[i]==2) { cout<<"The Statistical moment method cannot use CV_16S image type"<<endl; exit(-1); }
                    //*** BEGIN
                    Tmp[0]=clock();
                    data = mom_statistical_extraction(imageMethods);
                    Tmp[1]=clock();
                    //*** END
                    times.push_back(timeCounter(Tmp[1], Tmp[0]));
                    allDatas.push_back(data);
                    data.clear();
                    imageMethods.release();
                }else{ break; }
                break;
            //***************  GLCM  ***************
            case 3:
                if (methods[i]==1 || methods[i]==2) {
                    if      (methods[i]==1) { image8U.copyTo(imageMethods);  }
                    else if (methods[i]==2) { cout<<"The GLCM method cannot use CV_16S image type"<<endl; exit(-1); }
                    //*** BEGIN
                    Tmp[0]=clock();
                    data = glcm_extraction(imageMethods, 256);
                    Tmp[1]=clock();
                    //*** END
                    times.push_back(timeCounter(Tmp[1], Tmp[0]));
                    allDatas.push_back(data);
                    data.clear();
                    imageMethods.release();
                }else{ break; }
                break;
            //***************  LBP  ***************
            case 4:
                if (methods[i]==1 || methods[i]==2) {
                    if      (methods[i]==1) { image8U.copyTo(imageMethods);  }
                    else if (methods[i]==2) { cout<<"The LBP method cannot use CV_16S image type"<<endl; exit(-1); }
                    //*** BEGIN
                    Tmp[0]=clock();
                    data = lbp_extraction(imageMethods);
                    Tmp[1]=clock();
                    //*** END
                    times.push_back(timeCounter(Tmp[1], Tmp[0]));
                    allDatas.push_back(data);
                    data.clear();
                    imageMethods.release();
                }else{ break; }
                break;
            //*************** AHTD LUNG ***************
            case 5:
                if (methods[i]==1 || methods[i]==2) {
                    if (methods[i]==1) { cout<<"The AHTD LUNG method cannot use CV_8U image type"<<endl; exit(-1); }
                    else if (methods[i]==2) { image16S.copyTo(imageMethods); }
                    //*** BEGIN
                    Tmp[0]=clock();
                    data = ahtd_lung_extraction(imageMethods);
                    Tmp[1]=clock();
                    //*** END
                    times.push_back(timeCounter(Tmp[1], Tmp[0]));
                    allDatas.push_back(data);
                    data.clear();
                    imageMethods.release();
                }else{ break; }
                break;
            //*************** AHTD BRAIN ***************
            case 6:
                if (methods[i]==1 || methods[i]==2) {
                    if (methods[i]==1) { cout<<"The AHTD BRAIN method cannot use CV_8U image type"<<endl; exit(-1); }
                    else if (methods[i]==2) { image16S.copyTo(imageMethods); }
                    //*** BEGIN
                    Tmp[0]=clock();
                    data = ahtd_brain_extraction(imageMethods);
                    Tmp[1]=clock();
                    //*** END
                    times.push_back(timeCounter(Tmp[1], Tmp[0]));
                    allDatas.push_back(data);
                    data.clear();
                    imageMethods.release();
                }else{ break; }
                break;
            //*************** SCM AVERAGE ***************
            case 7:
                if (methods[i]==1 || methods[i]==2) {
                    //Mat to Sinal
                    struct Signal *Sinal;
                    Sinal=(struct Signal*)malloc(1*sizeof(struct Signal));
                    //Convert Mat img to struct Signal
                    if      (methods[i]==1) { image8U.copyTo(imageMethods);  }
                    else if (methods[i]==2) { image16S.copyTo(imageMethods); }
                    imgToStructSignal(imageMethods, Sinal, methods[i]);
                    //*** BEGIN
                    data = scm_extraction(Sinal, 0, Tmp);
                    //*** END
                    times.push_back(timeCounter(Tmp[1], Tmp[0]));
                    allDatas.push_back(data);
                    data.clear();
                    free((*Sinal).M.fDat);
                    free(Sinal);
                    imageMethods.release();
                }else{ break; }
                break;
            //*************** SCM LAPLACIAN ***************
            case 8:
                if (methods[i]==1 || methods[i]==2) {
                    //Mat to Sinal
                    struct Signal *Sinal;
                    Sinal=(struct Signal*)malloc(1*sizeof(struct Signal));
                    //Convert Mat img to struct Signal
                    if      (methods[i]==1) { image8U.copyTo(imageMethods);  }
                    else if (methods[i]==2) { image16S.copyTo(imageMethods); }
                    imgToStructSignal(imageMethods, Sinal, methods[i]);
                    //*** BEGIN
                    data = scm_extraction(Sinal, 1, Tmp);
                    //*** END
                    times.push_back(timeCounter(Tmp[1], Tmp[0]));
                    allDatas.push_back(data);
                    data.clear();
                    free((*Sinal).M.fDat);
                    free(Sinal);
                    imageMethods.release();
                }else{ break; }
                break;
            //*************** SCM GAUSSIAN ***************
            case 9:
                if (methods[i]==1 || methods[i]==2) {
                    //Mat to Sinal
                    struct Signal *Sinal;
                    Sinal=(struct Signal*)malloc(1*sizeof(struct Signal));
                    //Convert Mat img to struct Signal
                    if      (methods[i]==1) { image8U.copyTo(imageMethods);  }
                    else if (methods[i]==2) { image16S.copyTo(imageMethods); }
                    imgToStructSignal(imageMethods, Sinal, methods[i]);
                    //*** BEGIN
                    data = scm_extraction(Sinal, 2, Tmp);
                    //*** END
                    times.push_back(timeCounter(Tmp[1], Tmp[0]));
                    allDatas.push_back(data);
                    data.clear();
                    free((*Sinal).M.fDat);
                    free(Sinal);
                    imageMethods.release();
                }else{ break; }
                break;
            //*************** SCM SOBEL HV ***************
            case 10:
                if (methods[i]==1 || methods[i]==2) {
                    //Mat to Sinal
                    struct Signal *Sinal;
                    Sinal=(struct Signal*)malloc(1*sizeof(struct Signal));
                    //Convert Mat img to struct Signal
                    if      (methods[i]==1) { image8U.copyTo(imageMethods);  }
                    else if (methods[i]==2) { image16S.copyTo(imageMethods); }
                    imgToStructSignal(imageMethods, Sinal, methods[i]);
                    //*** BEGIN
                    data = scm_extraction(Sinal, 3, Tmp);
                    //*** END
                    times.push_back(timeCounter(Tmp[1], Tmp[0]));
                    allDatas.push_back(data);
                    data.clear();
                    free((*Sinal).M.fDat);
                    free(Sinal);
                    imageMethods.release();
                }else{ break; }
                break;
            //*************** SCM SOBEL H ***************
            case 11:
                if (methods[i]==1 || methods[i]==2) {
                    //Mat to Sinal
                    struct Signal *Sinal;
                    Sinal=(struct Signal*)malloc(1*sizeof(struct Signal));
                    //Convert Mat img to struct Signal
                    if      (methods[i]==1) { image8U.copyTo(imageMethods);  }
                    else if (methods[i]==2) { image16S.copyTo(imageMethods); }
                    imgToStructSignal(imageMethods, Sinal, methods[i]);
                    //*** BEGIN
                    data = scm_extraction(Sinal, 4, Tmp);
                    //*** END
                    times.push_back(timeCounter(Tmp[1], Tmp[0]));
                    allDatas.push_back(data);
                    data.clear();
                    free((*Sinal).M.fDat);
                    free(Sinal);
                    imageMethods.release();
                }else{ break; }
                break;
            //*************** SCM SOBEL V ***************
            case 12:
                if (methods[i]==1 || methods[i]==2) {
                    //Mat to Sinal
                    struct Signal *Sinal;
                    Sinal=(struct Signal*)malloc(1*sizeof(struct Signal));
                    //Convert Mat img to struct Signal
                    if      (methods[i]==1) { image8U.copyTo(imageMethods);  }
                    else if (methods[i]==2) { image16S.copyTo(imageMethods); }
                    imgToStructSignal(imageMethods, Sinal, methods[i]);
                    //*** BEGIN
                    data = scm_extraction(Sinal, 5, Tmp);
                    //*** END
                    times.push_back(timeCounter(Tmp[1], Tmp[0]));
                    allDatas.push_back(data);
                    data.clear();
                    free((*Sinal).M.fDat);
                    free(Sinal);
                    imageMethods.release();
                }else{ break; }
                break;
            //****************** TAMURA ******************
            case 13:
                if (methods[i]==1 || methods[i]==2) {
                    if      (methods[i]==1) { image8U.copyTo(imageMethods);  }
                    else if (methods[i]==2) { cout<<"The Tamura's method cannot use CV_16S image type"<<endl; exit(-1); }
                    //*** BEGIN
                    Tmp[0]=clock();
                    //data = tamura_extraction(imageMethods);
                    Tmp[1]=clock();
                    //*** END
                    times.push_back(timeCounter(Tmp[1], Tmp[0]));
                    allDatas.push_back(data);
                    data.clear();
                    imageMethods.release();
                }else{ break; }
                break;
            //************** ZERNICK MOMENTS *************
            case 14:
                if (methods[i]==1 || methods[i]==2) {
                    if      (methods[i]==1) { image8U.copyTo(imageMethods);  }
                    else if (methods[i]==2) { cout<<"The Zernick's method cannot use CV_16S image type"<<endl; exit(-1); }
                    //*** BEGIN
                    Tmp[0]=clock();
                    data = mom_zernick_extraction(imageMethods);
                    Tmp[1]=clock();
                    //*** END
                    times.push_back(timeCounter(Tmp[1], Tmp[0]));
                    allDatas.push_back(data);
                    data.clear();
                    imageMethods.release();
                }else{ break; }
                break;
            //****************** FOURIER *****************
            case 15:
                if (methods[i]==1 || methods[i]==2) {
                    if      (methods[i]==1) { image8U.copyTo(imageMethods);  }
                    else if (methods[i]==2) { cout<<"The Fourier's method cannot use CV_16S image type"<<endl; exit(-1); }
                    //*** BEGIN
                    Tmp[0]=clock();
                    //segundo parametro: vai gerar 7*4 atributos (padrao)
                    data = fourier_extraction(imageMethods, 7);
                    Tmp[1]=clock();
                    //*** END
                    times.push_back(timeCounter(Tmp[1], Tmp[0]));
                    allDatas.push_back(data);
                    data.clear();
                    imageMethods.release();
                }else{ break; }
                break;
            //*************************************************
            default:
                cout<<endl<<"Something is wrong with the amount of selected methods"<<endl<<endl;
                exit(-1);
                break;
        }
    }
    
    //Save Features
    saveFeatures(methods, allDatas, label, DesktopChar);
    
    //Save Times
    saveTimes(methods, times, DesktopChar);
    
    //Clear
    allDatas.clear();
    times.clear();
    imageMethods.release();
    
}
