#include "FRMWRKHeader.h"

int main(){
    //***********************************************************************************************************//
    //******************************* FRAMEWORK DE EXTRACAO DE ATRIBUTOS - LAPISCO ******************************//
    //***********************************************************************************************************//
    //
    // MAIS INFORMACOES SOBRE FRAMEWORK
    // Email: wellmend0@gmail.com
    // Whatsapp: +55 85 9 8915 7995
    //
    // ============================================================================================================
    // ==========================================___ CONFIGURANDO LIBS ___=========================================
    // ============================================================================================================
    //
    //       ==========================================___ MACOS ___=========================================
    //
    // OBS.: Todas as libs serao instaladas com auxilio do gerenciador de pacotes "Homebrew".
    //
    // >>>>> [ homebrew ] <<<<<
    //
    // # Instale o brew com o passo a passo nesse site: https://brew.sh
    //
    // # Insira esses comandos no terminal do macos:
    // - /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
    // - brew install wget
    //
    // >>>>> [ opencv ] <<<<<
    //
    // # Insira esses comandos no terminal do macos:
    // - brew tap homebrew/science
    // - brew install opencv@2
    // - brew link -f opencv@2
    //
    // # As .dylib estaram localizadas na pasta /usr/local/lib
    //
    // >>>>> [ dcmtk ] <<<<<
    //
    // # Insira esses comandos no terminal do macos:
    // - brew install cmake
    // - brew install dcmtk --with-libiconv
    // - brew link -f libiconv
    // - echo 'export PATH="/usr/local/opt/libiconv/bin:$PATH"' >> ~/.bash_profile
    // - brew install zlib
    // - brew link -f zlib
    // - brew link -f dcmtk
    //
    // # As .dylib estaram localizadas na pasta /usr/local/lib
    //
    //
    //       =========================================___ WINDOWS ___========================================
    //
    // Requisitos:
    // >> opencv2
    // >> dcmtk
    //
    // ============================================================================================================
    // ==========================================___     UTILIZACAO    ___=========================================
    //
    // > MacOS: Va em FRMWRKHeader.h e comente as linhas 4,5,6
    // > Windows: Va em FRMWRKHeader.h e descomente as linhas 4,5,6
    //
    // OBS.:
    // - Separe todas as imagens em pastas para cada uma das classes.
    // - A pasta deve ser nomeada com o numero da classe.
    // - As imagens devem ser nomeadas da seguinte forma: image_NUMERO
    // Exemplo 1: image_0
    // Exemplo 2: image_523
    //
    // - As seguintes instrucoes estao comentadas abaixo
    //
    // OBS2.:
    // Caso haja imagens medicas use a funcao DicomToMat descomentando.
    // DicomToMat(Mat &output, UHMax, UHMin, type)
    // > UHMax e UHMin devem ser definidos pelo usuario.
    //   Sao as faixas de valores que serao consideradas em UH.
    // > type
    //    1 - Imagem 8UC1  : Grayscale - Faixa padrão 0-255
    //    2 - Imagem 16SC1 : UH        - Unidade de Hounsfield
    //
    //
    // >>>>> SELECIONE OS METODOS PARA EXTRACAO <<<<<
    // 0 - Metodo NAO escolhido
    // 1 - Metodo escolhido utilizando imagem NORMALIZADA     (CV_8UC1 )(Faixa padrão 0-255        )
    // 2 - Metodo escolhido utilizando imagem NAO NORMALIZADA (CV_16SC1)(Diversas faixas de valores)
    //
    
    // >>>>> DEFINICAO DAS VARIAVEIS <<<<<
    vector<int> methods;
    char desktopDataAdrs[300];
    char imageFileAdrs[500];
    vector<int> samplesAmount;
    Mat image8U;
    Mat image16S;
    int counter=0;
    int label,sample;
    int i;
    
    // >>>>> SELECAO DOS METODOS <<<<<
    methods.push_back(0);//Hu Moments                  // Aceita somente 1
    methods.push_back(0);//Central Moments             // Aceita somente 1
    methods.push_back(0);//Statistical Moments         // Aceita somente 1
    methods.push_back(0);//GLCM                        // Aceita somente 1
    methods.push_back(0);//LBP                         // Aceita somente 1
    methods.push_back(0);//AHTD Lung                   // Aceita somente 2
    methods.push_back(0);//AHTD Brain                  // Aceita somente 2
    methods.push_back(0);//SCM AVERAGE                 // Aceita 1 ou 2
    methods.push_back(0);//SCM LAPLACIAN               // Aceita 1 ou 2
    methods.push_back(0);//SCM GAUSSIAN                // Aceita 1 ou 2
    methods.push_back(0);//SCM SOBEL HV                // Aceita 1 ou 2
    methods.push_back(0);//SCM SOBEL H                 // Aceita 1 ou 2
    methods.push_back(0);//SCM SOBEL V                 // Aceita 1 ou 2
    methods.push_back(0);//Tamura                      // Aceita somente 1
    methods.push_back(0);//Zernick Moments             // Aceita somente 1
    methods.push_back(0);//Fourier                     // Aceita somente 1 (FASE DE TESTE)
    
    // >>>>> AJUSTAGEM DAS VARIAVEIS <<<<<
    sprintf(desktopDataAdrs,"/Users/imacpedrosa/Desktop/DB");// Local onde sera salvo os bancos de dados gerados
    initializeAddresses(1, methods, desktopDataAdrs);// 1 - MacOS | 0 - Windows
    
    for (i=0; i<15; i++) {
        samplesAmount.push_back(50); // Quantidade de amostras da classe i
    }
    //Caso haja mais alguma classe, realize o mesmo processo : samplesAmount.push_back(quantidade_de_amostras);
    
    // >>>>> EXECUCAO DA EXTRACAO DE ATRIBUTOS DAS IMAGENS <<<<<
    for(label=0; label<samplesAmount.size(); label++){
        for(sample=0; sample<samplesAmount[label]; sample++){
            cout<<"image "<<counter+1<<" - "<<sum(samplesAmount)[0]<<endl;
            counter++;
            
            /* Carregamento de imagens */
            //sprintf(imageFileAdrs, "C:\\Users\\Wellington\\Desktop\\omni_virtual\\classe %.2d\\Omn%.2d_%.2d.jpg",label+1,label+1,sample+1);
            sprintf(imageFileAdrs, "/Users/imacpedrosa/Desktop/omni_real/classe %.2d/Omn%.2d_%.2d.jpg",label+1,label+1,sample+1);
            image8U = imread(imageFileAdrs, CV_LOAD_IMAGE_GRAYSCALE);
            
            /* DicomImage to 8U,16S */
            //DicomToMat(image8U, 120, 10, 1);
            //DicomToMat(image16S, 0, 0, 2);
            
            /* Inicialização dos Metodos */
            //Obs: Caso nao sejam usadas images DICOM passe a image8U nos dois primeiros parametros dessa funcao!
            initializeMethods(image8U, image8U, methods, label, desktopDataAdrs);
            
            /* Free */
            image8U.release();
            image16S.release();
        }
    }
    
    methods.clear();
    
    return 0;
}
