#ifndef scm_hpp
#define scm_hpp

//=======================================================
// CALCULO DA SCM E DOS SEUS ATRIBUTOS PARA SINAIS 1D OU 2D
// MIDE/SIM É UMA VERSÃO PARTICULAR DA SCM
// ESTE CÓDIGO ESTÁ BASEADO NO CÓDIGO DA MIDE
//
// DECLARACOES DE VARIÁVEIS
// Geraldo Ramalho - jan/2016
//=======================================================

#include "FRMWRKHeader.h"

// operacao com arquivo
void gsignal_read(char *filename, struct Signal* s);
void gsignal_write(char *filename,struct Signal* s);
void gmatrix_show(struct Matrix2* m);
void gatt_write(char *filename,struct scM* m);
// dll
void gcreate_filter(struct Filter* m, short f);
void gsignal_filt(struct Signal* s,struct Signal* f, struct Filter* m);
void gsignal_copy(struct Signal* s,struct Signal* f);
void gsignal_max(struct Signal* s1,struct Signal* s2);
void gscm(struct Signal* s,struct Signal* f,struct scM* m);
void gscm_map(struct Signal* s,struct Signal* f,struct Signal* map, short att);



#endif /* scm_hpp */
