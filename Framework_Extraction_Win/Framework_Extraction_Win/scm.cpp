//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

//=======================================================
// SCM: Structural Cooccurrence Matrix
// COMPUTES SCM AND ITS ATTRIBUTES FOR 1D OR 2D SIGNALS
// OBS.: MIDE/SIM IS A PARTICULAR CASE OF SCM
//       MIDE: Matriz de Interdependência Espacial
//       SIM: Spatial Interdependence Matrix
//
// Geraldo Ramalho - jan/2016
//
// v1.00  - first version of SCM implemented in C - code based on MIDE - jan/2016
//
// please cite:
//    Ramalho et al. Rotation-invariant Feature Extraction using a Structural Co-occurrence Matrix.
//    Measurement, v.??, n.??, p.???-???, 2016.
//    doi:10.1016/j.measurement.2016.08.012
//    lapisco.ifce.edu.br/?page_id=191
//
// TODO LIST:
// - mudar a estrutura dos ponteiros para informar se ja foram alocados
// - otimizar código (1024x1024=5seg)
// - mapa da SCM com janela de tamanho variável
// - acrescentar matriz de pesos W
// - acrescentar funcao de particionamento P
//=======================================================
//
#include "scm.hpp"

//=======================================================
// show matrix
//=======================================================
void gmatrix_show(struct Matrix2* m){
    short i;
    short j;
    for (i=0 ; i<(*m).nRows; i++){
        printf("%d: ",i);
        for (j=0 ; j<(*m).nCols ; j++){
            printf("%.3f ",(*m).fDat[i * (*m).nCols + j]);
        }
        printf("\n");
    }
}

//=======================================================
// copy signal
//=======================================================
void gsignal_copy(struct Signal* s,struct Signal* f){
    short i;
    short j;
    
    // copia
    (*f).M.fDat = (float*) malloc((*s).M.nRows*(*s).M.nCols*sizeof(float));
    (*f).M.nRows=(*s).M.nRows;
    (*f).M.nCols=(*s).M.nCols;
    (*f).nLevels=(*s).nLevels;
    for (i=0 ; i<(*s).M.nRows; i++){
        for (j=0 ; j<(*s).M.nCols ; j++){
            (*f).M.fDat[i * (*s).M.nCols + j] = (*s).M.fDat[i * (*s).M.nCols + j];
        }
    }
}

//=======================================================
// S1 = maximum between signals S1, S2
//=======================================================
void gsignal_max(struct Signal* s1,struct Signal* s2){
    short i;
    short j;
    
    // substitui
    for (i=0 ; i<(*s1).M.nRows; i++){
        for (j=0 ; j<(*s1).M.nCols ; j++){
            if ((*s1).M.fDat[i * (*s1).M.nCols + j] < (*s2).M.fDat[i * (*s2).M.nCols + j])
            {
                (*s1).M.fDat[i * (*s1).M.nCols + j] = (*s2).M.fDat[i * (*s2).M.nCols + j];
            }
            else
            {
                (*s2).M.fDat[i * (*s1).M.nCols + j] = (*s1).M.fDat[i * (*s1).M.nCols + j];
            }
        }
    }
}


//=======================================================
// select filter
//=======================================================
void gcreate_filter(struct Filter* m, short f){
    short i;
    short j;
    float par;
    if (f==fAVERAGE){
        (*m).M.fDat[0 * (*m).M.nCols + 0]=0.1111;
        (*m).M.fDat[0 * (*m).M.nCols + 1]=0.1111;
        (*m).M.fDat[0 * (*m).M.nCols + 2]=0.1111;
        (*m).M.fDat[1 * (*m).M.nCols + 0]=0.1111;
        (*m).M.fDat[1 * (*m).M.nCols + 1]=0.1111;
        (*m).M.fDat[1 * (*m).M.nCols + 2]=0.1111;
        (*m).M.fDat[2 * (*m).M.nCols + 0]=0.1111;
        (*m).M.fDat[2 * (*m).M.nCols + 1]=0.1111;
        (*m).M.fDat[2 * (*m).M.nCols + 2]=0.1111;
        (*m).fSum=1;
    }
    if (f==fLAPLACIAN){
        (*m).M.fDat[0 * (*m).M.nCols + 0]=0.1667;
        (*m).M.fDat[0 * (*m).M.nCols + 1]=0.6667;
        (*m).M.fDat[0 * (*m).M.nCols + 2]=0.1667;
        (*m).M.fDat[1 * (*m).M.nCols + 0]=0.6667;
        (*m).M.fDat[1 * (*m).M.nCols + 1]=-3.3333;
        (*m).M.fDat[1 * (*m).M.nCols + 2]=0.6667;
        (*m).M.fDat[2 * (*m).M.nCols + 0]=0.1667;
        (*m).M.fDat[2 * (*m).M.nCols + 1]=0.6667;
        (*m).M.fDat[2 * (*m).M.nCols + 2]=0.1667;
        (*m).fSum=1;
    }
    if (f==fGAUSSIAN){
        (*m).M.fDat[0 * (*m).M.nCols + 0]=0.0113;
        (*m).M.fDat[0 * (*m).M.nCols + 1]=0.0838;
        (*m).M.fDat[0 * (*m).M.nCols + 2]=0.0113;
        (*m).M.fDat[1 * (*m).M.nCols + 0]=0.0838;
        (*m).M.fDat[1 * (*m).M.nCols + 1]=0.6193;
        (*m).M.fDat[1 * (*m).M.nCols + 2]=0.0838;
        (*m).M.fDat[2 * (*m).M.nCols + 0]=0.0113;
        (*m).M.fDat[2 * (*m).M.nCols + 1]=0.0838;
        (*m).M.fDat[2 * (*m).M.nCols + 2]=0.0133;
        (*m).fSum=1;
    }
    if (f==fSOBELH){
        (*m).M.fDat[0 * (*m).M.nCols + 0]=1;
        (*m).M.fDat[0 * (*m).M.nCols + 1]=2;
        (*m).M.fDat[0 * (*m).M.nCols + 2]=1;
        (*m).M.fDat[1 * (*m).M.nCols + 0]=0;
        (*m).M.fDat[1 * (*m).M.nCols + 1]=0;
        (*m).M.fDat[1 * (*m).M.nCols + 2]=0;
        (*m).M.fDat[2 * (*m).M.nCols + 0]=-1;
        (*m).M.fDat[2 * (*m).M.nCols + 1]=-2;
        (*m).M.fDat[2 * (*m).M.nCols + 2]=-1;
        (*m).fSum=1;
    }
    if (f==fSOBELV){
        (*m).M.fDat[0 * (*m).M.nCols + 0]=1;
        (*m).M.fDat[1 * (*m).M.nCols + 0]=2;
        (*m).M.fDat[2 * (*m).M.nCols + 0]=1;
        (*m).M.fDat[0 * (*m).M.nCols + 1]=0;
        (*m).M.fDat[1 * (*m).M.nCols + 1]=0;
        (*m).M.fDat[2 * (*m).M.nCols + 1]=0;
        (*m).M.fDat[0 * (*m).M.nCols + 2]=-1;
        (*m).M.fDat[1 * (*m).M.nCols + 2]=-2;
        (*m).M.fDat[2 * (*m).M.nCols + 2]=-1;
        (*m).fSum=1;
    }
    if (f==9){
        par=0.1;
        for (i=0 ; i<(*m).M.nRows; i++){
            for (j=0 ; j<(*m).M.nCols ; j++){
                (*m).M.fDat[i * (*m).M.nCols + j]=par;
            }
        }
        (*m).fSum=1;
    }
}

//=======================================================
// k: signal filtering
//=======================================================
void gsignal_filt(struct Signal* s,struct Signal* f, struct Filter* m){
    short i;
    short j;
    float *TMP; // memoria temporaria
    float tmp;
    
    // Signal temporario espelhado
    TMP = (float*) malloc(((*s).M.nRows+2)*((*s).M.nCols+2)*sizeof(float));
    
    // copia dados
    for (i=0 ; i<(*s).M.nRows; i++){
        for (j=0 ; j<(*s).M.nCols ; j++){
            TMP[(i+1) * ((*s).M.nCols+2) + (j+1)] = (*s).M.fDat[i * (*s).M.nCols + j];
        }
    }
    // espelhamento 1a. e ultima coluna
    for (i=0 ; i<(*s).M.nRows+2; i++){
        TMP[i * ((*s).M.nCols+2)					] = TMP[i * ((*s).M.nCols+2) + 1			 ];
        TMP[i * ((*s).M.nCols+2) + ((*s).M.nCols+1)] = TMP[i * ((*s).M.nCols+2) + ((*s).M.nCols)];
    }
    // espelhamento 1a. e ultima linha
    for (j=0 ; j<(*s).M.nCols+2; j++){
        TMP[									  j ] = TMP[				 ((*s).M.nCols+2) + j];
        TMP[((*s).M.nRows+1) * ((*s).M.nCols+2) + j] = TMP[((*s).M.nRows) * ((*s).M.nCols+2) + j];
    }
    // filtragem
    (*f).M.fDat = (float*) malloc((*s).M.nRows*(*s).M.nCols*sizeof(float));
    (*f).M.nRows=(*s).M.nRows;
    (*f).M.nCols=(*s).M.nCols;
    (*f).nLevels=(*s).nLevels;
    for (i=0 ; i<(*s).M.nRows; i++){
        for (j=0 ; j<(*s).M.nCols ; j++){
            tmp=((*m).M.fDat[	  	       0]*TMP[i	    * ((*s).M.nCols+2) + j]+(*m).M.fDat[				 1]*TMP[i	  * ((*s).M.nCols+2) + j+1]+(*m).M.fDat[				   2]*TMP[i		* ((*s).M.nCols+2) + j+2]+
                 (*m).M.fDat[	  (*m).M.nCols]*TMP[(i+1) * ((*s).M.nCols+2) + j]+(*m).M.fDat[	(*m).M.nCols + 1]*TMP[(i+1) * ((*s).M.nCols+2) + j+1]+(*m).M.fDat[	  (*m).M.nCols + 2]*TMP[(i+1) * ((*s).M.nCols+2) + j+2]+
                 (*m).M.fDat[2 * (*m).M.nCols]*TMP[(i+2) * ((*s).M.nCols+2) + j]+(*m).M.fDat[2 * (*m).M.nCols + 1]*TMP[(i+2) * ((*s).M.nCols+2) + j+1]+(*m).M.fDat[2 * (*m).M.nCols + 2]*TMP[(i+2) * ((*s).M.nCols+2) + j+2])
            /(*m).fSum;
            (*f).M.fDat[i * (*s).M.nCols + j] = round(tmp);
        }
    }
    
    // libera memoria temporaria
    free(TMP);
    TMP=NULL;
}


//=======================================================
// loads a matrix (or image in text file)
//=======================================================
void gsignal_read(char *filename,struct Signal* s){
    short i;
    short j;
    int p;
    
    FILE *fp;
    fp = fopen(filename, "r");
    if( fp == NULL ){
        printf("FILE NOT FOUND!");
    }
    else
    {
        (*s).nLevels=0;
        (*s).minLevel=999;
        (*s).maxLevel=0;
        // carrega tamanho do sinal
        fscanf ((FILE*)fp, "%d ", &(*s).M.nRows);
        fscanf ((FILE*)fp, "%d ", &(*s).M.nCols);
        //
        if ((*s).M.nRows>MAXLIN) (*s).M.nRows=MAXLIN;
        if ((*s).M.nCols>MAXCOL) (*s).M.nCols=MAXCOL;
        
        
        // aloca memoria do sinal
        (*s).M.fDat = (float*) malloc((*s).M.nRows*(*s).M.nCols*sizeof(float));
        // carrega sinal
        i=0; j=0;
        while ( !feof(fp ) &&	 i<(*s).M.nRows){
            fscanf ((FILE*)fp, "%d ", &p); // carrega amostra
            (*s).M.fDat[i * (*s).M.nCols + j]=(float)p;
            if (p>(*s).nLevels) (*s).nLevels=p;
            if (p<(*s).minLevel) (*s).minLevel=p;
            if (p>(*s).maxLevel) (*s).maxLevel=p;
            if (j++==(*s).M.nCols-1) {j=0; i++;}
        }
        fclose(fp);
        
        // quantizacao para a SCM (min-max)
        for (i=0 ; i<(*s).M.nRows; i++){
            for (j=0 ; j<(*s).M.nCols ; j++){
                (*s).M.fDat[i * (*s).M.nCols + j]= round(
                                                         ( (*s).M.fDat[i * (*s).M.nCols + j] - (float)(*s).minLevel )  *
                                                         (float)(MAXNIV-1) / ( (float)(*s).maxLevel - (float)(*s).minLevel )
                                                         );
            }
        }
        (*s).nLevels=MAXNIV;
    }
}

//=======================================================
// write a matrix in a text file
//=======================================================
void gsignal_write(char *filename,struct Signal* s){
    short i;
    short j;
    
    FILE *fp;
    fp = fopen(filename, "w");
    if( fp == NULL ){
        printf("FILE WRITE ERROR!");
    }
    else
    {
        // escreve tamanho do sinal
        fprintf ((FILE*)fp, "%d ", (*s).M.nRows);
        fprintf ((FILE*)fp, "%d\n", (*s).M.nCols);
        //
        
        // escreve sinal
        for (i=0 ; i<(*s).M.nRows; i++){
            for (j=0 ; j<(*s).M.nCols ; j++){
                fprintf ((FILE*)fp, "%f ", (*s).M.fDat[i * (*s).M.nCols + j]); // escreve amostra
            }
            fprintf ((FILE*)fp, "\n");
        }
        fclose(fp);
    }
}

//=======================================================
// write SCM attributes in text file
//=======================================================
void gatt_write(char *filename,struct scM* m){
    short i;
    
    FILE *fp;
    
    fp = fopen(filename, "w");
    if( fp == NULL ){
        printf("FILE WRITE ERROR!");
    }
    else
    {
        
        // escreve tamanho do sinal
        fprintf ((FILE*)fp, "%d ", 1);
        fprintf ((FILE*)fp, "%d\n", MAXATT);
        //
        
        // escreve sinal
        for (i=0 ; i<MAXATT; i++){
            fprintf ((FILE*)fp, "%f ", (*m).fAtt[i]); // escreve amostra
        }
        fprintf ((FILE*)fp, "\n");
        fclose(fp);
    }
}
//=======================================================
// computes SCM
//=======================================================
void gscm(struct Signal* s, struct Signal* f, struct scM *m){
    short i,j,n,i1,i2;
    float *pi; // memoria temporaria
    float *pj; // memoria temporaria
    float mi,mj,vi,vj,tmp;
    
    pi = (float*) malloc(((*s).nLevels)*((*s).nLevels)*sizeof(float)); // memoria temporaria
    pj = (float*) malloc(((*s).nLevels)*((*s).nLevels)*sizeof(float)); // memoria temporaria
    
    // inicializa matriz
    (*m).M.nRows=(*s).nLevels;
    (*m).M.nCols=(*s).nLevels;
    (*m).nLevels=(*s).nLevels;
    for (i=0 ; i<(*m).M.nRows; i++){
        for (j=0 ; j<(*m).M.nCols ; j++){
            (*m).M.fDat[i * (*m).M.nCols + j] = 0.0;
        }
    }
    
    // calcula frequencias
    float var;
    for (i=0 ; i<(*s).M.nRows; i++){
        for (j=0 ; j<(*s).M.nCols; j++){
            i1=(*s).M.fDat[i * (*s).M.nCols + j];
            i2=(*f).M.fDat[i * (*s).M.nCols + j];
            //(*m).M.fDat[(int)(i1 * (*m).nLevels + i2)] += 1.0;
            var=(i1 * (*m).nLevels + i2);
            if (var>=((*m).nLevels)*((*m).nLevels)) {
                var=(*m).nLevels*(*m).nLevels-1;
            }
            if (var<0) {
                var=0;
            }
            (*m).M.fDat[(int)var] += 1.0;
        }
    }
    /*
     W= ones(N,N);
     for i=1:N
     for j=1:N
     W(i,j)= abs(i-j);
     end
     end
     W=(W+1)/N; % diagonal-distance weigths
     end
     // multiplica pelos pesos Wij=abs(i-j)/(MAXNIV-1)
     for (i=0 ; i<(*m).nLevels; i++){
     for (j=0 ; j<(*m).nLevels ; j++){
     (*m).M.fDat[(int)(i * (*m).nLevels + j)]*=fabsf((float)i-(float)j)/(MAXNIV-1);
     }
     }
     */

    // normaliza a matriz
    tmp=(*s).M.nRows*(*s).M.nCols;
    for (i=0 ; i<(*m).nLevels; i++){
        for (j=0 ; j<(*m).nLevels ; j++){
            (*m).M.fDat[(int)(i * (*m).nLevels + j)]/=tmp;
        }
    }
    
    // -------- atributos da SCM assimetrica
    //
    // calcula probabilidades Oi e Ei da SCM assimetrica
    for (i=0 ; i<(*m).nLevels; i++){
        pi[i]=0.0;
        pj[i]=0.0;
        for (j=0 ; j<(*m).nLevels; j++){
            if (i==j) pi[i]+=(*m).M.fDat[(i * (*m).nLevels + j)]; // Oi - diagonal
            if (j>=i) pj[i]+=(*m).M.fDat[(i * (*m).nLevels + j)]; // Ei - triang. sup
        }
    }
    // CHI ou CSD  Chi-square distance between expected and observed diagonals
    (*m).fAtt[aCHI]=0.0; // Chi quadrado
    for (i=0 ; i<(*m).nLevels; i++){
        if (pj[i]>0) (*m).fAtt[aCHI] += (pi[i]-pj[i])*(pi[i]-pj[i])/pj[i];
    }
    
    // calcula probabilidades marginais da SCM assimetrica
    tmp=0;
    for (i=0 ; i<(*m).nLevels; i++){
        pi[i]=0.0; // Pr
        pj[i]=0.0; // Pc
        for (j=0 ; j<(*m).nLevels ; j++){
            pi[i]+=(*m).M.fDat[i * (*m).nLevels + j]; // prob.marginal linhas Pr
            pj[i]+=(*m).M.fDat[j * (*m).nLevels + i]; // prob.marginal colunas Pc
        }
    }
    
    // DKL, CAD
    (*m).fAtt[aDKL]=0.0; // Kullback and Leibler divergence
    (*m).fAtt[aCAD]=0.0; // complementary absolute difference
    for (i=0 ; i<(*m).nLevels; i++){
        if (pj[i]>0 && pi[i]>0) (*m).fAtt[aDKL] += log(pi[i]/pj[i]) * (pj[i]);
        else (*m).fAtt[aDKL]=0;
        (*m).fAtt[aCAD] += fabsf(pi[i] - pj[i]);
    }
    (*m).fAtt[aCAD] = 1-(*m).fAtt[aCAD];
    
    // MDR
    // calcula valor esperado das diferenças absolutas
    mi=0.0; // MDPr
    mj=0.0; // MDPc
    for (i=0 ; i<(*m).nLevels; i++){
        for (j=0 ; j<(*m).nLevels ; j++){
            if (pi[i]>0 && pj[j]>0)	mi+=pi[i]*pi[j] * fabsf((float)i-(float)j);
        }
    }
    tmp=mi;
    if (mj<mi) tmp=mj; // MD1=min(MDPc,MDPr)
    if (mi>mj) mj=mi; // MD2=max(MDPc,MDPr)
    (*m).fAtt[aMDR]=tmp/mj; // mean absolute difference ratio = MD1/MD2
    
    
    // CSR
    // calcula probabilidades de QI e QII da SCM assimetrica
    // prob. QI
    n=round((*m).nLevels/2);
    for (i=0 ; i<n; i++){
        pi[i]=0.0; // PI
        for (j=0 ; j<n ; j++){
            pi[i]+=(*m).M.fDat[i * (*m).nLevels + j];
        }
    }
    for (i=n ; i<(*m).nLevels; i++){
        pj[i]=0.0; // Pm=(PI+PII)/2
        for (j=n ; j<(*m).nLevels ; j++){
            pj[j]+=(pi[j]+(*m).M.fDat[(i-n) * (*m).nLevels + (j-n)])/2;
        }
    }
    (*m).fAtt[aCSR]=0.0; // ratio between QI and QIII of matrix coefficientes - distancia chi-quadrado
    for (j=0 ; j<n ; j++){
        if (pj[i]>0) (*m).fAtt[aCSR]+=(pi[j]-pj[j]) * (pi[j]-pj[j]) / pj[j];
    }
    
    
    // calcula a matriz simetrica
    for (i=0 ; i<(*m).nLevels; i++){
        for (j=i ; j<(*m).nLevels ; j++){
            tmp=((*m).M.fDat[i * (*m).nLevels + j]+(*m).M.fDat[j * (*m).nLevels + i])/2;
            (*m).M.fDat[i * (*m).nLevels + j]=tmp;
            (*m).M.fDat[j * (*m).nLevels + i]=tmp;
        }
    }
    
    // -------- atributos da SCM simetrica
    //
    // calcula probabilidades marginais
    for (i=0 ; i<(*m).nLevels; i++){
        pi[i]=0.0; // Pr
        pj[i]=0.0; // Pc
        for (j=0 ; j<(*m).nLevels ; j++){
            pi[i]+=(*m).M.fDat[i * (*m).nLevels + j]; // prob.marginal linhas
            pj[i]+=(*m).M.fDat[j * (*m).nLevels + i]; // prob.marginal colunas
        }
    }
    // calcula media e variancia
    mi=0.0;
    mj=0.0;
    vi=0.0;
    vj=0.0;
    for (i=0 ; i<(*m).nLevels ; i++){
        mi+=pi[i]*(float)i;
        mj+=pj[i]*(float)i;
    }
    for (i=0 ; i<(*m).nLevels ; i++){
        vi+=((float)i-mi)*((float)i-mi)*pi[i];
        vj+=((float)i-mj)*((float)i-mj)*pj[i];
    }
    
    // COR e MDI e ENT
    (*m).fAtt[aCOR]=0.0; // correlacao
    (*m).fAtt[aMDI]=0.0; // momento de diferencas inverso
    (*m).fAtt[aENT]=0.0; // entropia
    tmp=(mi*mj);
    for (i=0 ; i<(*m).nLevels; i++){
        for (j=0 ; j<(*m).nLevels ; j++){
            (*m).fAtt[aCOR] += ( (*m).M.fDat[i * (*m).nLevels + j] * ((float)i - mi) * ((float)j - mj ) );
            (*m).fAtt[aMDI] += (*m).M.fDat[i * (*m).nLevels + j] / (1 + (float)abs(i -j));
            if ((*m).M.fDat[i * (*m).nLevels + j]>0)
                (*m).fAtt[aENT] += (*m).M.fDat[i * (*m).nLevels + j]*log((*m).M.fDat[i * (*m).nLevels + j]);
        }
    }
    if (sqrt( vi * vj ) >0 )
        (*m).fAtt[aCOR]/=sqrt( vi * vj );
    else
        (*m).fAtt[aCOR]=1.0;
    (*m).fAtt[aENT]=-(*m).fAtt[aENT];
    
    
    // libera memoria temporaria alocada
    free(pi); pi=NULL;
    free(pj); pj=NULL;
    
}

//=======================================================
// computes SCM map (windowing SCM)
//=======================================================
void gscm_map(struct Signal* s,struct Signal* f,struct Signal* map, short att){
# define TAMJAN 3
    short i,wi;
    short j,wj;
    float *TMPS; // memoria temporaria
    float *TMPF; // memoria temporaria
    float tmp,wss,wsf,wsmax,wsmin,wfmax,wfmin;
    struct Signal *WS;
    struct Signal *WF;
    struct scM *mw;	// SCM da janela NxN
    
    // Signal original e filtrado espelhado
    TMPS = (float*) malloc(((*s).M.nRows+2)*((*s).M.nCols+2)*sizeof(float));
    TMPF = (float*) malloc(((*s).M.nRows+2)*((*s).M.nCols+2)*sizeof(float));
    // copia dados
    for (i=0 ; i<(*s).M.nRows; i++){
        for (j=0 ; j<(*s).M.nCols ; j++){
            TMPS[(i+1) * ((*s).M.nCols+2) + (j+1)] = (*s).M.fDat[i * (*s).M.nCols + j];
            TMPF[(i+1) * ((*f).M.nCols+2) + (j+1)] = (*f).M.fDat[i * (*f).M.nCols + j];
        }
    }
    // espelhamento 1a. e ultima coluna
    for (i=0 ; i<(*s).M.nRows+2; i++){
        TMPS[i * ((*s).M.nCols+2)					 ] = TMPS[i * ((*s).M.nCols+2) + 1			   ];
        TMPS[i * ((*s).M.nCols+2) + ((*s).M.nCols+1)] = TMPS[i * ((*s).M.nCols+2) + ((*s).M.nCols)];
        TMPF[i * ((*f).M.nCols+2)					 ] = TMPF[i * ((*f).M.nCols+2) + 1			   ];
        TMPF[i * ((*f).M.nCols+2) + ((*f).M.nCols+1)] = TMPF[i * ((*f).M.nCols+2) + ((*f).M.nCols)];
    }
    // espelhamento 1a. e ultima linha
    for (j=0 ; j<(*s).M.nCols+2; j++){
        TMPS[										j] = TMPS[				   ((*s).M.nCols+2) + j];
        TMPS[((*s).M.nRows+1) * ((*s).M.nCols+2) + j] = TMPF[((*s).M.nRows) * ((*s).M.nCols+2) + j];
        TMPF[										j] = TMPF[				   ((*f).M.nCols+2) + j];
        TMPF[((*f).M.nRows+1) * ((*f).M.nCols+2) + j] = TMPF[((*f).M.nRows) * ((*f).M.nCols+2) + j];
    }
    //
    WS=(struct Signal*)malloc(1*sizeof(struct Signal));
    (*WS).M.fDat = (float*) malloc(TAMJAN*TAMJAN*sizeof(float));
    (*WS).M.nRows=TAMJAN;
    (*WS).M.nCols=TAMJAN;
    (*WS).nLevels=(*s).nLevels;
    WF=(struct Signal*)malloc(1*sizeof(struct Signal));
    (*WF).M.fDat = (float*) malloc(TAMJAN*TAMJAN*sizeof(float));
    (*WF).M.nRows=TAMJAN;
    (*WF).M.nCols=TAMJAN;
    (*WF).nLevels=(*f).nLevels;
    //
    mw=(struct scM*)malloc(1*sizeof(struct scM));
    (*mw).M.fDat = (float*) malloc((*s).nLevels*(*s).nLevels*sizeof(float));
    
    // mapa da SCM TAMJANxTAMJAN
    (*map).M.fDat = (float*) malloc((*s).M.nRows*(*s).M.nCols*sizeof(float));
    (*map).M.nRows=(*s).M.nRows;
    (*map).M.nCols=(*s).M.nCols;
    (*map).nLevels=(*s).nLevels;
    //
    for (i=0 ; i<(*s).M.nRows; i++){
        for (j=0 ; j<(*s).M.nCols ; j++){
            //
            wsmax=-99999; wsmin=99999;
            wfmax=-99999; wfmin=99999;
            for (wi=0 ; wi<TAMJAN ; wi++){
                for (wj=0 ; wj<TAMJAN ; wj++){
                    (*WS).M.fDat[wi * (*WS).M.nCols + wj] = (TMPS[(i+wi) * ((*s).M.nCols+2) + j+wj]);
                    (*WF).M.fDat[wi * (*WF).M.nCols + wj] = (TMPF[(i+wi) * ((*f).M.nCols+2) + j+wj]);
                    if ((*WS).M.fDat[wi * (*WS).M.nCols + wj]>wsmax) wsmax=(*WS).M.fDat[wi * (*WS).M.nCols + wj];
                    if ((*WS).M.fDat[wi * (*WS).M.nCols + wj]<wsmin) wsmin=(*WS).M.fDat[wi * (*WS).M.nCols + wj];
                    if ((*WF).M.fDat[wi * (*WF).M.nCols + wj]>wfmax) wfmax=(*WF).M.fDat[wi * (*WF).M.nCols + wj];
                    if ((*WF).M.fDat[wi * (*WF).M.nCols + wj]<wfmin) wfmin=(*WF).M.fDat[wi * (*WF).M.nCols + wj];
                }
            }
            wss=0; wsf=0;
            if ((wsmax-wsmin)>0) wss=(MAXNIV-1)/(wsmax-wsmin);
            if ((wfmax-wfmin)>0) wsf=(MAXNIV-1)/(wfmax-wfmin);
            
            // renormaliza e requantiza a janela antes de calcular a SCM
            for (wi=0 ; wi<TAMJAN ; wi++){
                for (wj=0 ; wj<TAMJAN ; wj++){
                    (*WS).M.fDat[wi * (*WS).M.nCols + wj] = round(((*WS).M.fDat[wi * (*WS).M.nCols + wj]-wsmin)*wss);
                    (*WF).M.fDat[wi * (*WF).M.nCols + wj] = round(((*WF).M.fDat[wi * (*WF).M.nCols + wj]-wfmin)*wsf);
                }
            }
            
            
            // calcula a SCM da janela
            gscm((struct Signal*) WS,(struct Signal*) WF,(struct scM*) mw);
            tmp=(*mw).fAtt[att];
            (*map).M.fDat[i * (*map).M.nCols + j] = tmp;
        }
    }
    
    // libera memoria temporaria
    free((*WS).M.fDat);	(*WS).M.fDat=NULL;
    free(WS);	WS = NULL;
    free((*WF).M.fDat);	(*WF).M.fDat=NULL;
    free(WF);	WF = NULL;
    free((*mw).M.fDat);	(*mw).M.fDat=NULL;
    free(mw);	mw = NULL;
    free(TMPS);	TMPS=NULL;
    free(TMPF);	TMPF=NULL;
}
