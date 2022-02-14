#include<sort_eigSt.hpp>
#include <utils/logger.hpp>
#include <ERRORS.h>
sort_eigSt::sort_eigSt(int new_Ne,int new_Nm,int new_sysType,int new_eigsToFind,int new_printLowestNStates)
{
  if(new_sysType!=0 && new_sysType!=1) 
    { std::cerr << "\n sort_eigSt::sort_eigSt() Dont use this ctor for other sysType than 0,1.\n";
    exit(1);
    }
  ctorBody(new_Ne,new_Nm,new_sysType,new_eigsToFind,
	   new_printLowestNStates,0);
}

sort_eigSt::sort_eigSt(int new_Ne,int new_Nm,int new_sysType,int new_eigsToFind,int new_printLowestNStates,int NrParSteps)
{
  ctorBody(new_Ne,new_Nm,new_sysType,new_eigsToFind,
	   new_printLowestNStates,NrParSteps);
}

void sort_eigSt::ctorBody(int new_Ne,int new_Nm,int new_sysType,int new_eigsToFind,int new_printLowestNStates,int NrParSteps)
{
  Ne=new_Ne;
  Nm=new_Nm;
  sysType=new_sysType;
  eigsToFind=new_eigsToFind;
  printLowestNStates=new_printLowestNStates;
  switch(sysType){
  case 0: // homogeneous system
    nrOutChnl=((Ne/2+1)*(Ne/2+2))/2;
    break;
  case 1: // magnetic impurity
    nrOutChnl=Ne/2+1;
    break;
  case 2: // sweep_B_IMI
    nrOutChnl=1;
    break;
  default:
    std::cerr << "\nsort_eigSt::sort_eigSt. Unknown sysType.\n";exit(1);break;
  }
  res=new eigSt[Nm*eigsToFind*nrOutChnl];
  sr=new LIST_EIGST*[nrOutChnl];
  srEnd=new int[nrOutChnl];
  outF_En=new FILE*[nrOutChnl];
  switch(sysType){
  case 0:
  case 1:
    outF_bas=new FILE*[(Ne/2+1)*Nm];
    outF_vec=new FILE*[nrOutChnl];
    break;
  case 2:
    outF_bas=new FILE*[Nm];
    outF_vec=new FILE*[NrParSteps];
    break;
  }
  s=0;


// Initialize arrays for sorting the eigSt and open output files for energies and eigvecs
  switch(sysType){
  case 0: // homogeneous system
    for(int totS=(Ne%2);totS<=Ne;totS=totS+2)
      for(reqSz=(Ne%2);reqSz<=totS;reqSz=reqSz+2)
	{
	  sr[s]=new LIST_EIGST[Nm*eigsToFind+1];
	  sr[s][0].p_next=NULL;
	  sr[s][0].p_st=NULL;
	  srEnd[s]=1;
	  sprintf(filename,"out/vec%1d%1d.dat",s/10,s%10);
	  outF_vec[s]=testfopenW(filename);
	  sprintf(filename,"out/En%1d%1d.dat",s/10,s%10);
	  outF_En[s]=testfopenW(filename);
	  fprintf(outF_En[s],"# Aspect ratios vs. Energies per particle; totS=%3.1f, totSz=%3.1f\n",totS/2.,reqSz/2.);
	  fprintf(outF_vec[s],"# Energies (as in %s) and corresponding eigvecs (only for the first aspect ratio)\n",filename);
	  s++;
	}
    break;

  case 1: // magnetic impurity
    for(reqSz=(Ne%2);reqSz<=Ne;reqSz=reqSz+2)
	{
	  sr[s]=new LIST_EIGST[Nm*eigsToFind+1];
	  sr[s][0].p_next=NULL;
	  sr[s][0].p_st=NULL;
	  srEnd[s]=1;
	  sprintf(filename,"out/vec%1d%1d.dat",s/10,s%10);
	  outF_vec[s]=testfopenW(filename);
	  sprintf(filename,"out/En%1d%1d.dat",s/10,s%10);
	  outF_En[s]=testfopenW(filename);
	  fprintf(outF_En[s],"# Aspect ratios vs. Energies per particle; totSz=%3.1f\n",reqSz/2.);
	  fprintf(outF_vec[s],"# Energies (as in %s) and corresponding eigvecs (only for the first aspect ratio)\n",filename);
	  s++;
	}
    break;
  case 2: // sweep_B_IMI
    sr[0]=new LIST_EIGST[Nm*eigsToFind+1];
    sr[0][0].p_next=NULL;
    sr[0][0].p_st=NULL;
    srEnd[0]=1;
    sprintf(filename,"out/En00.dat");
    outF_En[0]=testfopenW(filename);
    fprintf(outF_En[0],"# Mag. fields vs. Energies per particle\n");
    for(int i=0;i<NrParSteps;i++)
      {
	sprintf(filename,"out/vec00_step%1d%1d%1d.dat",i/100,(i/10)%10,i%10);
	outF_vec[i]=testfopenW(filename);
	sprintf(filename,"out/En00.dat");	
	fprintf(outF_vec[i],"# Energies (as in %s, row %d) and corresponding eigvecs\n",filename,i);
      }
    s++;
    break;
  }
  if(s!=nrOutChnl) {printf("\n Problem: Count all the S,Sz states again.\n");exit(2);}


  // Open files for bases 
  switch(sysType){
  case 0:
  case 1:
    for(int i=0;i<(Ne/2)+1;i++)
      for(int reqJ=0;reqJ<Nm;reqJ++)
	{
	  sprintf(filename,"out/basSz%1d%1d_J%1d%1d.dat",i/10,i%10,reqJ/10,reqJ%10);
	  outF_bas[i*Nm+reqJ]=testfopenW(filename);	
	  fprintf(outF_bas[i*Nm+reqJ],"# Basis for the J=%2d eigvecs with 2*Sz=%2d (the same for all values of the parameter)\n",reqJ,2*i+(Ne%2));
	}
    break;
  case 2:
    for(int reqJ=0;reqJ<Nm;reqJ++)
      {
	sprintf(filename,"out/basJ%1d%1d.dat",reqJ/10,reqJ%10);
	outF_bas[reqJ]=testfopenW(filename);	
	fprintf(outF_bas[reqJ],"# Basis for the J=%2d eigvecs (the same for all values of the parameter)\n",reqJ);
      }
    break;
  }    
}

sort_eigSt::~sort_eigSt()
{  
  for(int j=0;j<nrOutChnl;j++)
    {delete[] sr[j];fclose(outF_En[j]);fclose(outF_vec[j]);}
  switch(sysType){
  case 0:
  case 1:
    for(int i=0;i<((Ne/2)+1)*Nm;i++)
      fclose(outF_bas[i]);
    break;
  case 2:
    for(int i=0;i<Nm;i++)
      fclose(outF_bas[i]);
    break;
  }          
  delete[] sr;delete[] srEnd;delete[] outF_En;
  delete[] outF_vec;delete[] outF_bas;
  delete[] res;
}


int sort_eigSt::sortAndWrite(int res_ind,int write_vectors,double par)
{
  return sortAndWrite(res_ind,write_vectors,par,0);
}

int sort_eigSt::sortAndWrite(int res_ind,int write_vectors,double par,int NrParStep)
{
  // Sorting results; e.g. S=0(Sz=0),S=1(Sz=0,1),S=2(Sz=0,1,2),etc.
  //                      s=      0         1 2         3 4 5  ...
  for(int k=0;k<res_ind;k++)
    {
      switch(sysType){
      case 0: // homogeneous system
	itotS=(int)rint(2*res[k].totS)-2;itotSz=(int)rint(2*res[k].Sz);
	if(itotS<0) s=itotSz/2;
	else
	  s=((itotS/2+1)*(itotS/2+2))/2+itotSz/2;
	break;
      case 1: // magnetic impurity
	itotSz=(int)rint(2*res[k].Sz);
	s=itotSz/2;
	break;
      case 2: // sweep_B_IMI
	s=0;
	break;
      }
      //if(srEnd[s]>printLowestNStates) continue;
      tmp=sr[s];
      n_sorted=0;
      while(tmp->p_next!=NULL)
	{
	  if(((res+k)->En)<(((tmp->p_next)->p_st)->En)) break;
	  tmp=tmp->p_next;
	  //if(n_sorted>printLowestNStates) break;
	  //n_sorted++;
	}
      //if(n_sorted>printLowestNStates) continue;
      if(srEnd[s]>Nm*eigsToFind) 
	{
	  printf("\n Problem while sorting states: srEnd>=Nm*eigsToFind\n");
	  exit(1);
	}
      sr[s][srEnd[s]].p_next=tmp->p_next;
      sr[s][srEnd[s]].p_st=res+k;
      tmp->p_next=&(sr[s][srEnd[s]]);
      srEnd[s]++;
    }


  // Flush energies and eigvecs out to files
  int ind_vecFile=0;   // File to write the eigvecs to
  switch(sysType){
  case 0:
  case 1:
    ind_vecFile=s;break;
  case 2:
    ind_vecFile=NrParStep;
    fseek(outF_vec[ind_vecFile],-1,SEEK_CUR);
    fprintf(outF_vec[ind_vecFile],"; B=%f T\n",par);
  }
  for(s=0;s<nrOutChnl;s++) // Go through all channels
    {
      // Write value of the parameter
      fprintf(outF_En[s],"%f ",par);
      tmp=sr[s]->p_next;
      // In each channel go through all states
      for(int j=1;j<srEnd[s];j++)
	{
	  // Write energy
	  fprintf(outF_En[s],"%f ",tmp->p_st->En/Ne);
	  // Write vectors
	  if(write_vectors)
	    {
	      fprintf(outF_vec[ind_vecFile],"%f ",tmp->p_st->En/Ne);
	      for(int l=0;l<tmp->p_st->dimension();l++)
		fprintf(outF_vec[ind_vecFile],"%f ",tmp->p_st->outCoef(l));
	      fprintf(outF_vec[ind_vecFile],"\n");
	    }
	  tmp=tmp->p_next;
	  if(j>=printLowestNStates && printLowestNStates>0) break;
	}
      // Write J's
      fprintf(outF_En[s],"\n");
      fprintf(outF_En[s],"# totJ's  ");
      tmp=sr[s]->p_next;
      for(int j=1;j<srEnd[s];j++)
	{
	  fprintf(outF_En[s],"%5.2f     ",tmp->p_st->totJ);
	  tmp=tmp->p_next;
	  if(j>=printLowestNStates && printLowestNStates>0) break;
	}
      fprintf(outF_En[s],"\n");
      // Where relevant (magnetic impurity, sweep_B_IMI), write totS's
      if(sysType==1 || sysType==2) 
	{
	  fprintf(outF_En[s],"# totS's  ");
	  tmp=sr[s]->p_next;
	  for(int j=1;j<srEnd[s];j++)
	    {
	      fprintf(outF_En[s],"%8.5f  ",tmp->p_st->totS);
	      tmp=tmp->p_next;
	      if(j>=printLowestNStates && printLowestNStates>0) break;
	    }
	  fprintf(outF_En[s],"\n");
	}
      // Where relevant (sweep_B_IMI), write totSz's
      if(sysType==2) 
	{
	  fprintf(outF_En[s],"# totSz's ");
	  tmp=sr[s]->p_next;
	  for(int j=1;j<srEnd[s];j++)
	    {
	      fprintf(outF_En[s],"%8.5f  ",tmp->p_st->Sz);
	      //fprintf(outF_En[s],"unknown   ");	      
	      tmp=tmp->p_next;
	      if(j>=printLowestNStates && printLowestNStates>0) break;
	    }
	  fprintf(outF_En[s],"\n");
	}
    }


  // Reinitialize the arrays used for sorted states
  for(s=0;s<nrOutChnl;s++)
    {
      sr[s][0].p_next=NULL;
      sr[s][0].p_st=NULL;
      srEnd[s]=1;
      fflush(outF_En[s]);
    }
  return 0;
}

FILE* sort_eigSt::getBasisFile(int reqSz,int reqJ)
{
  if(sysType!=0 && sysType!=1)
    {
      std::cerr << "\nsort_eigSt::getBasisFile(int,int) Use this method for sysType=0,1 only.\n";
      exit(1);
    }
  return outF_bas[(reqSz/2)*Nm+reqJ];
}

FILE* sort_eigSt::getBasisFile(int reqJ)
{
  if(sysType!=2)
    {
      std::cerr << "\nsort_eigSt::getBasisFile(int) Use this method for sysType=2 only.\n";
      exit(1);
    }
  return outF_bas[reqJ];
}

eigSt *sort_eigSt::getRes()
{return res;}
