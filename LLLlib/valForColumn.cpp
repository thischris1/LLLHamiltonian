#include <valForColumn.hpp>
#include <eigSt.hpp>

valForColumn::valForColumn(int new_whatToCompute,
			   yosBasis *new_yosbas,
			   int new_matType,
			   eigSt &stateToEval,
			   double aTob,
			   double alpha1,
			   double alpha2)
{
  //eigSt state(new_yosbas);

  whatToCompute=new_whatToCompute;
  yosbas=new_yosbas;
  matType=new_matType;
  
  // initialize all (not needed) objects to null
  //densOp = NULL;
  pDensOp = NULL;
  pJxOp = NULL;
  pJyOp = NULL;
  spinDensOp = NULL;
  spinDownDensOp = NULL;
  spinZDensOp = NULL;
  spinXDensOp = NULL;


  switch(whatToCompute){
  case 0:
    strcpy(descr,"constant");
    printf("   - none\n");
    break;

  case 1:
    strcpy(descr,"particle density");
    printf("   - %s (<n(r)>)\n",descr);
    pDensOp = new DensOperator (*yosbas, stateToEval, matType);
    break;

  case 2: // To calculate spin density
    strcpy(descr,"spin density");
    printf("   - %s (<S^2.n(r)>, ** i.e. local value of S(S+1) **)\n",descr);
    spinDensOp=new spinDensOperator(*yosbas, stateToEval, matType);  
    break;

  case 3:  // To calculate spin down density  
    strcpy(descr,"spin-down density");
    printf("   - %s (<n^-(x)>)\n",descr);
    spinDownDensOp=new spinDownDensOperator(*yosbas, stateToEval, matType);
    break;

  case 4:  // To calculate Sz density
    strcpy(descr,"spin z-component density");
    printf("   - %s (<S_z.n(r)>)\n",descr);
    spinZDensOp=new spinZDensOperator(*yosbas, stateToEval, matType);
    break;
    
  case 5: // To calculate spin density / particle density
    strcpy(descr,"spin density/particle density");
    printf("   - %s (<S^2.n(r)>/<n(r)>, ** S(S+1) recalc-ed to S **)\n",descr);
    spinDensOp=new spinDensOperator(*yosbas, stateToEval, matType);  
    pDensOp=new DensOperator(*yosbas, stateToEval, matType);
    break;

  case 6:  // To calculate spin down density / particle density  
    strcpy(descr,"spin-down density/particle density");
    printf("   - %s (<n^-(r)>/<n(r)>)\n",descr);
    spinDownDensOp=new spinDownDensOperator(*yosbas, stateToEval, matType);
    pDensOp=new DensOperator(*yosbas, stateToEval, matType);
    break;

  case 7:  // To calculate Sz density / particle density
    strcpy(descr,"spin z-component density/particle density");
    printf("   - %s (<S_z.n(r)>/<n(r)>)\n",descr);
    spinZDensOp=new spinZDensOperator(*yosbas, stateToEval, matType);
    pDensOp=new DensOperator(*yosbas, stateToEval, matType);
    break;

  case 8:
    strcpy (descr, "current density x-component");
    pJxOp = new jxOp (*yosbas, stateToEval, matType);
    pJxOp->setSolenoidFluxes (alpha1, alpha2);
    break;

  case 9:
    strcpy (descr, "current density y-component");
    pJyOp = new jyOp (*yosbas, stateToEval, aTob,  matType);
    pJyOp->setSolenoidFluxes (alpha1, alpha2);
    break;

  case 10:  // To calculate Sx density
    strcpy(descr,"spin x-component density");
    printf("   - %s (<S_x.n(r)>)\n",descr);
    spinXDensOp=new spinXDensOperator(*yosbas, stateToEval, matType);
    break;
    
  case 11:  // To calculate Sx density / particle density
    strcpy(descr,"spin x-component density/particle density");
    printf("   - %s (<S_x.n(r)>/<n(r)>)\n",descr);
    spinXDensOp=new spinXDensOperator(*yosbas, stateToEval, matType);
    pDensOp=new DensOperator(*yosbas, stateToEval, matType);
    break;

  default:
    std::cerr << "valForColumn::valForColumn Unimplemented value of the quantity to be computed.\n";
    exit(1);
  };
}

valForColumn::~valForColumn() 
{
  //  if (densOp != NULL) delete densOp;
  if (pDensOp != NULL) delete pDensOp;
  if (pJxOp != NULL) delete pJxOp;
  if (pJyOp != NULL) delete pJyOp;
  if (spinDensOp != NULL) delete spinDensOp;
  if (spinZDensOp != NULL) delete spinZDensOp;
  if (spinXDensOp != NULL) delete spinXDensOp;
  if (spinDownDensOp != NULL) 
    {
      /*  std::ofstream fsLandau;
      int Nm = yosbas->getNm();
      fsLandau.open ("spinDownLandau.dat");
      for (int i=0; i<Nm; i++)
	fsLandau << spinDownDensOp->allLandauMatricesElement(0,0,i,i) << "\n";
      fsLandau.close();
      */
      delete spinDownDensOp;
    }
}

double valForColumn::getValForColumn(double x,double y)
{
  double dens,spSqDens;

  //std::cout << "getValForColumn (" << whatToCompute << ")\n";

  switch(whatToCompute){
  case 0:
  default:
    break;
    
  case 1: // To calculate particle density
    return pDensOp->evalAt(x,y);
    break;

  case 2: // To calculate spin density
    //spSqDens=spinDensOp->DensAtXY(wf,yosbas->dimension(),x,y);
    spSqDens=spinDensOp->evalAt(x,y);
    return spSqDens;  
    break;

  case 3:  // To calculate spin down density  
    return spinDownDensOp->evalAt(x,y);
    break;

  case 4:  // To calculate Sz density
    return spinZDensOp->evalAt(x,y);
    break;

  case 5: // To calculate spin density/particle density
    if((dens=pDensOp->evalAt(x,y))==0)
      return 0;
    spSqDens=spinDensOp->evalAt(x,y);
    return (sqrt(1+4*spSqDens/dens)-1)/2.0;  // compute S from S(S+1)
    break;

  case 6:  // To calculate spin down density/particle density
    if((dens=pDensOp->evalAt(x,y))==0)
      return 0;
    return spinDownDensOp->evalAt(x,y)/dens;
    break;

  case 7:  // To calculate Sz density/particle density
    if((dens=pDensOp->evalAt(x,y))==0)
      return 0;
    return spinZDensOp->evalAt(x,y)/dens;
    break;

  case 8:
    return pJxOp->evalAt(x,y);
    break;

  case 9:
    return pJyOp->evalAt(x,y);
    break;

  case 10:  // To calculate Sx density
    return spinXDensOp->evalAt(x,y);
    break;

  case 11:  // To calculate Sx density/particle density
    if((dens=pDensOp->evalAt(x,y))==0)
      return 0;
    return spinXDensOp->evalAt(x,y)/dens;
    break;

  }
  return 1.0;
}

const char *valForColumn::getDescrOfColumn()
{
  return descr;
}

void valForColumn::writeLandauDiag (char *strFile)
{
 
  if (pDensOp != NULL)
    {
      TOCCNUMLIST occnum;
      if ( pDensOp->getOccupationNumbers (occnum) == 0)
	{
	  std::ofstream fsLandau (strFile);

	  TOCCNUMLIST::iterator run;
	  for (run = occnum.begin(); run != occnum.end(); run++)
	    fsLandau << run->maxj << ' ' << run->diag << ' ' << run->occ << '\n';

	  fsLandau.close();
	}
    }
}

























