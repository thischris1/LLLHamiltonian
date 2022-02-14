#ifndef SORT_EIGST_HPP
#define SORT_EIGST_HPP

#include<eigSt.hpp>
#include<fstream>
#include<iostream>
#include<myaccessories.hpp>

class sort_eigSt{
private:
  typedef struct list_eigst {eigSt* p_st; struct list_eigst *p_next;} LIST_EIGST;
  
  static const int FILENAME_MAX_LEN=30;

  int sysType;
  int Ne,Nm;
  int nrOutChnl;
  eigSt *res;  
  LIST_EIGST *tmp;
  LIST_EIGST **sr;  
  int *srEnd;
  int s,itotS,itotSz,printLowestNStates,n_sorted,eigsToFind,reqSz;
  FILE **outF_En,**outF_vec,**outF_bas;
  FILE *basisfile;
  char filename[FILENAME_MAX_LEN];

  void ctorBody(int new_Ne,int new_Nm,int new_sysType,int new_eigsToFind,int new_printLowestNStates,int NrParSteps);

public:
  sort_eigSt(int new_Ne,
	     int new_Nm,
	     int new_sysType,
	     int new_eigsToFind,
	     int new_printLowestNStates);
  
  sort_eigSt(int new_Ne,
	     int new_Nm,
	     int new_sysType,
	     int new_eigsToFind,
	     int new_printLowestNStates,
	     int NrParSteps);

  ~sort_eigSt();

  int sortAndWrite(int res_ind,int write_vectors,double aspect);
  int sortAndWrite(int res_ind,int write_vectors,double aspect,int NrParStep);
  FILE* getBasisFile(int reqSz,int reqJ);
  FILE* getBasisFile(int reqJ);
  eigSt *getRes();
};


#endif
