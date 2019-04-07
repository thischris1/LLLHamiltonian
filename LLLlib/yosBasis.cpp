/*! \file yosBasis.cpp
 */

/*!
  \class yosBasis : public Basis

  \brief Basis consisting of yosStates (see yosState ).

*/

// Implementation of the class yosBasis

#include <yosBasis.hpp>
#include <ERRORS.h>
/*! \fn yosBasis::yosBasis(int new_Ne,int new_Nm,int new_spinYes)
   \brief   Stadard constructor (some of the basis-generating routine should
   be called after that. 
   
   \param int new_Ne         Nr. of electrons
   \param int new_Nm         Nr. of flux quanta
   \param int new_spinYes    spinYes=1 -> els. with spin (i.e. not necessarily
                             spin polarized states); spinYes=0 -> els. fully
			     polarized, i.e. effectivelly without spin

*/
yosBasis::yosBasis(int new_Ne,int new_Nm,int new_spinYes) : Basis(NULL,NULL,0),
Ne_(new_Ne),
Nm_(new_Nm),
capacity(0),
jsUpSignature  (0),
jsDoSignature  (0),
logfile(0),
m_WF_x(-1.0),
m_WF_y(-1.0),
m_pWaveFct (0),
  m_pWaveFctdx (0),
  m_pWaveFctdy (0),
state_(0)
{ 
  
  obs_totJ.knownInAdvance = false;
  obs_totS.knownInAdvance = false;
  obs_totSz.knownInAdvance = false;
  if(new_spinYes) spinYes=1; else spinYes=0;
    

}


/*! \fn yosBasis::yosBasis(FILE *in,FILE *logfile)
  \brief Get the basis from a file _with header_ (as produced e.g. by
  writeBasisHeaderToFile() + writeBasisToFile() ).
*/
yosBasis::yosBasis(FILE *in,FILE *logfile) : Basis(NULL,NULL,0),
spinYes(1),
capacity(0),
logfile(0),
state_ (0)
{
   setLogFile(logfile);
  findStartLine(in);
  Ne_=getdec(in);
  Nm_=getdec(in);
  dim=getdec(in);
  findStartLine(in);
  getBasisFromFile(in,dim);

  m_WF_x = -1.0;
  m_WF_y = -1.0;
  m_pWaveFct = NULL;
  m_pWaveFctdx = NULL;
  m_pWaveFctdy = NULL;

  obs_totJ.knownInAdvance = false;
  obs_totS.knownInAdvance = false;
  obs_totSz.knownInAdvance = false;

  jsDoSignature = NULL; jsUpSignature = NULL;
}


yosBasis::yosBasis(const yosBasis &rhs):
        Basis(rhs),
        Ne_(rhs.getNe()),
        Nm_(rhs.getNm()),
        spinYes ( rhs.getSpinYes()),
        jsUpSignature(0),
        jsDoSignature(0),
	logfile(0),
	m_pWaveFct(0),
        m_pWaveFctdx(0),
        m_pWaveFctdy(0),
        state_(0)
{
        // Do copying here?

	// copy states (if available)

}
yosBasis::yosBasis(const std::string & fileName)
:
Basis(0,0,0),
Ne_(0),
Nm_(0),
jsUpSignature(0),
jsDoSignature(0),
logfile(0),
m_pWaveFct(0),
m_pWaveFctdx(0),
m_pWaveFctdy(0),
state_(0)
{
FILE *inFile = fopen(fileName.c_str(),"r");
readFromFile(inFile);

}


yosBasis yosBasis::operator =(const yosBasis &rhs)
{
        if (this != &rhs)
        {
                Ne_=rhs.getNe();
                Nm_=rhs.getNm();
                spinYes = rhs.getSpinYes();



        }
        return (*this);
}

yosBasis::~yosBasis() 
{
  //std::cout << "!!! Destructor yosBasis !!! \n";
	if (state_) 
	{
  		for(int i=0;i<dim;i++) 
  		{ 
  			delete state_[i];
  		}
  		delete[] state_;
	}

  if (m_pWaveFct != NULL) delete[] m_pWaveFct;
  if (m_pWaveFctdx != NULL) delete[] m_pWaveFctdx;
  if (m_pWaveFctdy != NULL) delete[] m_pWaveFctdy;
  if (jsDoSignature != NULL) 
  	{
  		delete [] jsDoSignature;
  		jsDoSignature = 0;
  	}
  if (jsUpSignature) 
  {
  		delete [] jsUpSignature;
  		jsUpSignature =0;
  }
  
}

yosState *yosBasis::operator[](long int i) const 
{
		return state_[i];
}



yosState * yosBasis::getState(long index) const
{
	 return state_[index];
	
}

int yosBasis::getNe() const 
{
	return Ne_;
}

int yosBasis::getNm() const 
{
	return Nm_;
}

int yosBasis::getSpinYes() const 
{
	return spinYes;
}

int yosBasis::increaseNm(int increment) 
{
  if(increment<0) return -1;
  Nm_ += increment;
  for(int i=0;i<dim;i++)
   // state_[i]->increaseNm(increment);
  glLogger.warning("new Nm_:%d",Nm_);
  return 0;
}

double yosBasis::norm(int i_st) {return 1.0;}


/* Set this aspect ratio for all states of the basis.
   See also yosState::setAspect. 
 */
int yosBasis::broadcastAspectToAllStates(double new_aToB)
{
  int retVal=0;
  for(int i=0;i<dim;i++)
    retVal+=state_[i]->setAspect(new_aToB);
  if(retVal) return 1; /* There has been some problem while 
			  setting the aspect ratios. */
  return 0;
}




/*********************************************	 
 Output routines 
 *********************************************/

void yosBasis::setLogFile(FILE *new_logfile)
{
  logfile=new_logfile;
}


void yosBasis::writeBasisToFile(FILE *out)
{
  for(int i=0;i<dim;i++)
    {
      for(int l=0;l<Ne_;l++)
	{
	  fprintf(out,"%d",state_[i]->j(l));
	  if(state_[i]->spin(l)) fprintf(out,"+");
	  else fprintf(out,"-");
	}
      fprintf(out,"\n");fflush(out);
    }
}

void yosBasis::writeBasisHeaderToFile(FILE *out)
{
  fprintf(out,"%3d                 # Ne\n",Ne_);
  fprintf(out,"%3d                 # Nm\n",Nm_);
  fprintf(out,"%10ld          # dim\n",dim);
  fprintf(out,"# States of the basis are listed below:\n");
}


/* Get the basis from a file _without header_. */

int yosBasis::getBasisFromFile(FILE *in,int new_dim)
{
  int *j_l,*spin_l;
  char c;

  //if(logfile) fprintf(logfile,"Got state ");
  dim=new_dim;
  j_l=allocK<int>(Ne_);
  spin_l=allocK<int>(Ne_);
  state_=new yosState*[dim];
  capacity = dim;
  for(int i=0;i<dim;i++)
    {
      for(int l=0;l<Ne_;l++)
	{
	  switch(getJSCharByChar(in,j_l+l,spin_l+l)){
	  case -1 : ;
	  case  1 :throw CxErrors("Wrong type of the input file");
	  default:;
	  }
	}
      if((c=fgetc(in))!='\n')throw CxErrors("Wrong type of the input file");       
      state_[i]=new yosState(Ne_,Nm_,j_l,spin_l);
   //   if(logfile) {state_[i]->fprint(logfile);fprintf(logfile,", ");}
    }
  return 0;
}

double yosBasis::getAspectRatio()const
{
  return (state_[0]->getAspect());
}
int yosBasis::mergeBasis(yosBasis &basToAdd)
{  
  int new_dim;
  if(Ne_!=basToAdd.getNe() || Nm_!=basToAdd.getNm())
    {glLogger.error("yosBasis::mergeBasis Bases incompatible (Ne or Nm are not equal).\n");exit(1);}
  new_dim=dimension()+basToAdd.dimension();
  yosState **new_state=new yosState*[new_dim];
  int i;
  for(i=0;i<dimension();i++) // The new basis consists of the old one...
    new_state[i]=state_[i];
  for(;i<new_dim;i++) // ...plus the new one.
    new_state[i]= new yosState(*(basToAdd[i-dimension()])); // deep copies of the basToAdd states
  state_=new_state;
  dim=new_dim;
  capacity = new_dim;
  return 0;
}



int yosBasis::amITheSameAs(FILE *in)
{
  int j_l,spin_l;
  char c;
  for(int i=0;i<dim;i++)
    {
      for(int l=0;l<Ne_;l++)
	{
	  switch(getJSCharByChar(in,&j_l,&spin_l)){
	  case  1 : 
	  	throw CxErrors("TheSameAs has been reading from a file"
			   "which is of different type than a file"
			   "produced by yosBasis::writeBasisToFile",
			   __FILE__, __LINE__); 
	  	case -1 : return 0; // The bases are different; dimensions are not the same
	  default:;
	  }
	  if(j_l-state_[i]->j(l)) return 0;
	  // The bases are different; one pair of j's is not the same
	  if(spin_l-state_[i]->spin(l)) return 0;
	  // The bases are different; one pair of spins is not the same
	}
      if((c=fgetc(in))!='\n') 
	{
	  ungetc(c,in);
	  if(getJSCharByChar(in,&j_l,&spin_l)) return 0; 
	  // The bases are not the same, one of the states is 'longer' (more particles)
	  
	  throw CxErrors("TheSameAs has been reading from a file"
			 "which is of different type than a file"
			 "produced by yosBasis::writeBasisToFile",
			 __FILE__, __LINE__); 
	  
	}
    }
  if(feof(in)) return 1; // The bases are the same
  if(getJSCharByChar(in,&j_l,&spin_l)) 
  {
  	
  // The bases are not the same, basis in the file has larger dimension
  	throw CxErrors("TheSameAs has been reading from a file"
		 "which is of different type than a file"
		 "produced by yosBasis::writeBasisToFile",
		 __FILE__, __LINE__); 
  	
  }
  return 0;
}

  // Private auxiliary routine for amITheSameAs()
  
  int yosBasis::getJSCharByChar(FILE *in,int *j_l,int *spin)
  {
    char c=0;
    *j_l=0;
    while(!feof(in))
      {
	c=fgetc(in);
	if(c<'0' || c>'9') break;
	*j_l=(*j_l)*10+(c-'0');
      }
    if(c=='+') {*spin=1;return 0;}
    if(c=='-') {*spin=0;return 0;}
    if(feof(in)) return -1; // premature EOF => basis probably too short (dim too low)
    return 1; // unexpected character encountered.
  }




/********************************
Routines for generation of various bases
 ********************************/

// genBasisForConst_J_noSpin(...) and related methods

/*! \fn int yosBasis::genBasisForConst_J_noSpin(int reqJ,long int guessDim)
  \brief Construct a basis composed of all (fully polarized) states
  with totJ (j_1+..+j_Ne mod Nm = totJ) equal to reqJ. 

  Dimension of the basis must not exceed guessDim.
*/
int yosBasis::genBasisForConst_J_noSpin(int reqJ,long int guessDim)
{
  int nr,next;
  int *js;

  obs_totJ.knownInAdvance = true; obs_totJ.iVal=reqJ;
  obs_totS.knownInAdvance = false;
  obs_totSz.knownInAdvance = false;

  if(spinYes)  throw CxErrors("yosBasis::genBasisForConst_J_Sz cannot be used to"
		   "initialize a spinpolarized basis",
		   __FILE__, __LINE__);
  
  glLogger.info(" - Spin polarized, all %d-electron states with total J=%d.\n",Ne_,reqJ);
  dim=0;
  js=new int[Ne_];
  for(int i=0;i<Ne_;i++) js[i]=0;
  setCapacity(guessDim);
  //state_= new yosState*[guessDim];
  nr=0;next=0;  
  while(!next) 
    {
      if(totJ(js)==reqJ)
	if(checkPauli_withoutSpin(js)) 
	  {
	    if(dim>=guessDim) return 1;
	    state_[dim]=new yosState(Ne_,Nm_,js);
	    if(logfile!=NULL)
	      {
		fprintf(logfile,"%ld: ",dim);state_[dim]->fprint(logfile);
		fprintf(logfile,"\n");fflush(logfile);
	      }
	    dim++;	      
	  }
      next=nextSt(js);
    }
  if(!dim) {
  	 throw CxErrors("yosBasis::genBasisForConst_J_noSpin found no vector matching the conditions.", __FILE__, __LINE__);
  }
  adjustCapacityToSize();
  return 0;
}

/*
int yosBasis::genBasis_allJ_noSpin()
{
  int nr,next;
  int *js;

  int new_dim;

  // dimension = # of states(Nm) over # of fermions(Ne)
  double tmp_dim = 1.0;
  for (int i=0; i<Ne_; i++) tmp_dim = tmp_dim * (Nm_-i)/(Ne_-i);
  new_dim = (int)rint(tmp_dim);

  state_ = new (yosState*)[new_dim];

  if(spinYes) err_msg(2);
  printf(" - Spin polarized, all %d-electron states \n",Ne_);
  if(logfile!=NULL)
    fprintf(logfile," - Spin polarized, all %d-electron states\n",Ne_);

  dim=0;
  js=new int[Ne_];
  for(int i=0;i<Ne_;i++) js[i]=0;
  
  nr=0;next=0;  
  while(!next) 
    {
      if(checkPauli_withoutSpin(js)) 
	{	  
	  yosState *pNew_state = new yosState(Ne_,Nm_,js);
	  state_[dim] = pNew_state;

	  if(logfile!=NULL)
	    {
	      fprintf(logfile,"%ld: ",dim);pNew_state->fprint(logfile);
	      fprintf(logfile,"\n");fflush(logfile);
	    }
	  dim++;	      
	}
      next=nextSt(js);
    }
  if(!dim) err_msg(3);

  assert (dim == new_dim);

  return 0;
}
*/

/*! \fn int yosBasis::genBasis_allJ_noSpin()
  \brief Construct basis from all possible Ne-electron fully polarized
  states.
*/
int yosBasis::genBasis_allJ_noSpin()
{
  int nr,next;
  int *js;

  int new_dim;

  obs_totJ.knownInAdvance = false;
  obs_totS.knownInAdvance = false;
  obs_totSz.knownInAdvance = false;


  // dimension = # of states(Nm) over # of fermions(Ne)
  double tmp_dim = 1.0;
  for (int i=0; i<Ne_; i++) tmp_dim = tmp_dim * (Nm_-i)/(Ne_-i);
  new_dim = (int)rint(tmp_dim);

  setCapacity(new_dim);
  //state_ = new yosState*[new_dim];

  if(spinYes) 
   	throw CxErrors ("yosBasis::genBasisForConst_J_noSpin can be used to initialize only a spin polarized basis.", __FILE__, __LINE__);

  glLogger.info(" - Spin polarized, all %d-electron states \n",Ne_);
  dim=0;
  js=new int[Ne_];
  for(int i=0;i<Ne_;i++) js[i]=i;
  
  nr=0;
  do 
    {
      yosState *pNew_state = new yosState(Ne_,Nm_,js);
      state_[dim] = pNew_state;

      if(logfile!=NULL)
	{
	  fprintf(logfile,"%ld: ",dim);pNew_state->fprint(logfile);
	  fprintf(logfile,"\n");fflush(logfile);
	}
      dim++;
      next = nextSt_nospin(js, Ne_-1, Nm_);

    }
  while (!next);

  if(!dim) {
  	throw CxErrors("yosBasis::genBasisForConst_J_noSpin found no vector matching the conditions.", __FILE__, __LINE__);
  }

  assert (dim == new_dim);

  adjustCapacityToSize();

  return 0;
}

/*! \fn int yosBasis::genBasis_allJ_noSpin()
  \brief Construct basis from all possible Ne-electron fully polarized
  states.
*/
int yosBasis::genBasis_allSz_allJ()
{
  int nr,next;

  int new_dim;

  obs_totJ.knownInAdvance = false;
  obs_totS.knownInAdvance = false;
  obs_totSz.knownInAdvance = false;


  // dimension = (2 * # of states(Nm)) over # of fermions(Ne) (2* for spin)
  double tmp_dim = 1.0;
  for (int i=0; i<Ne_; i++) tmp_dim = tmp_dim * (2*Nm_-i)/(Ne_-i);
  new_dim = (int)rint(tmp_dim);

  setCapacity(new_dim);
  //state_ = new yosState*[new_dim];

  if(!spinYes) 
  {  
    glLogger.error("yosBasis::genBasis_allSz_allJ Attempted to use a spin polarized basis.\n");
    exit(1);  
  }
  glLogger.info(" - All %d-electron states, all J's, all Sz's \n",Ne_);


  int *js=new int[Ne_];
  int *spin=new int[Ne_];
  int maxS=1;
  int totSz,itmp;
  for(int i=0;i<Ne_;i++) {js[i]=0;maxS*=2;}
  nr=0;next=0;  
  while(!next) 
    {
      for(int i=0;i<maxS;i++)
      {
	totSz=0;
	itmp=i;
	for(int j=0;j<Ne_;j++)
	{spin[j]=itmp%2;totSz+=itmp%2;itmp=itmp/2;}
	totSz=2*totSz-Ne_;  // actually not needed here
	if(!checkPauli_withSpin(js,spin)) continue;
	if(dim>=new_dim) 
	{
	  glLogger.error(
	    "yosBasis::genBasis_allSz_allJ Learn me to count the basis states correctly!\n");
	  exit(1);
	}
	state_[dim]=new yosState(Ne_,Nm_,js,spin);
	if(logfile!=NULL)
	{
	  fprintf(logfile,"%ld: ",dim);state_[dim]->fprint(logfile);
	  fprintf(logfile,"\n");fflush(logfile);
	}
	dim++;	      
      }
      next=nextSt(js);
    }
  if(!dim) 
    {
      glLogger.warning("yosBasis::genBasis_allSz_allJ"); 
      glLogger.warning("found no vector matching the conditions.");
      glLogger.warning(
       "Note that totSz is 2*required tot. Sz, i.e. for 3 particles, totSz ought to be odd.");
    }
  adjustCapacityToSize();

  return 0;
}



int yosBasis::checkPauli_withoutSpin(int *spat)
  /* returns: 0=excluded by Pauli principle; 1=allowed;
     it's simple here: only states with two equal j's are excluded
  */
{ 
  for(int i=0;i<Ne_-1;i++)
    if(spat[i]==spat[i+1]) return 0;
  return(1);                            
}

// End of methods related to genBasisForConst_J_noSpin()


/* genBasis_allSz_constJ(...) (some methods dedicated primarily to
   genBasisForConst_J_Sz() are used; see below)
*/

/*! \fn int yosBasis::genBasis_allSz_constJ(int reqJ,long int guessDim)
  \brief Construct a basis composed of all states (with spin freedom)
  with totJ (j_1+..+j_Ne mod Nm = totJ) equal to reqJ. 
  

  Dimension of the basis must not exceed guessDim.
 */
int yosBasis::genBasis_allSz_constJ(int reqJ,long int guessDim)
{
  int nr,totSz,maxS,itmp,next;
  int *js,*spin;

  obs_totJ.knownInAdvance = true; obs_totJ.iVal=reqJ;
  obs_totS.knownInAdvance = false;
  obs_totSz.knownInAdvance = false;

  if(!spinYes) {
  	
  	 throw CxErrors("yosBasis::genBasisForConst_J_Sz cannot be used to"
		    "initialize a spinpolarized basis",
		    __FILE__, __LINE__);
  }
  glLogger.info(" - Spin unpolarized, all %d-electron states with total J=%d and all possible spin orientations.\n",Ne_,reqJ);fflush(stdout);
  dim=0;
  js=new int[Ne_];
  spin=new int[Ne_];
  maxS=1;
  for(int i=0;i<Ne_;i++) {js[i]=0;maxS*=2;}
  
  setCapacity(guessDim);
  //state_=new yosState*[guessDim];
  nr=0;next=0;  
  while(!next) 
    {
      if(totJ(js)==reqJ)
	{
	  for(int i=0;i<maxS;i++)
	    {
	      totSz=0;
	      itmp=i;
	      for(int j=0;j<Ne_;j++)
		{spin[j]=itmp%2;totSz+=itmp%2;itmp=itmp/2;}
	      totSz=2*totSz-Ne_;  // actually not needed here
  	      if(!checkPauli_withSpin(js,spin)) continue;
	      if(dim>=guessDim) return 1;
	      state_[dim]=new yosState(Ne_,Nm_,js,spin);
	      if(logfile!=NULL)
		{
		  fprintf(logfile,"%ld: ",dim);state_[dim]->fprint(logfile);
		  fprintf(logfile,"\n");fflush(logfile);
		}
	      dim++;	      
	    }
	}
      next=nextSt(js);
    }

  if(dim==0) 
    {
      glLogger.warning("yosBasis::genBasis_allSz_constJ");
      glLogger.warning("found no vector matching the conditions. ");
      glLogger.warning(
       "Note that totSz is 2*required tot. Sz, i.e. for 3 particles, totSz ought to be odd.");
    }
  adjustCapacityToSize();
  return 0;
}



// genBasisForConst_J_Sz(...) and related methods
/*! \fn int yosBasis::genBasisForConst_J_Sz(int reqJ,int reqSz,long int guessDim)
  \brief Construct a basis composed of all states
  with totJ (j_1+..+j_Ne mod Nm = totJ) equal to reqJ with z-component of
  the total spin equal to reqSz. 

  Dimension of the basis must not exceed guessDim.
 */
int yosBasis::genBasisForConst_J_Sz(int reqJ,int reqSz,long int guessDim)
{


  obs_totJ.knownInAdvance = true; obs_totJ.iVal=reqJ;
  obs_totS.knownInAdvance = false;
  obs_totSz.knownInAdvance = true; 
  obs_totSz.dVal = reqSz*1./2.;obs_totSz.iVal = reqSz;

  if(!spinYes) 
    {
      glLogger.error("yosBasis::genBasisForConst_J_Sz"); 
      glLogger.error("cannot be used to initialize a spinpolarized basis.");
      exit(1);
    }
  glLogger.info(" - Spin unpolarized, all %d-electron states with total J=%d and total spin 2*totSz=%d.\n",Ne_,reqJ,reqSz);

  dim = createAll_oneSz_oneJ(Ne_,Nm_,reqJ,reqSz,guessDim);
  
  return 0;
}

/*! \fn int yosBasis::createAll_oneSz_oneJ(int Ne,int Nm,int reqJ,int reqSz,long int guessDim,yosState **output)
				  
  \brief The core of genBasisForConst_J_Sz(...) (protected method). Puts the
  generated states into output.
*/
int yosBasis::createAll_oneSz_oneJ(int Ne,int Nm,
				   int reqJ,int reqSz,long int guessDim)
{
  vector<yosState*> tmpVec(guessDim);
  const int tmpStr_len = 5*Ne+3;
  char *tmpStr=new char[tmpStr_len];
  int nr,totSz=0,maxS,itmp,next;
  int *js,*spin;

  if(Ne != Ne_) 
    {
      glLogger.warning("yosBasis::createAll_oneSz_oneJ");
      glLogger.warning("called with Ne != Ne_");
    }

  glLogger.debug("((%d/%d) states, reqSz=%d, reqJ=%d)",Ne,Nm,reqSz,reqJ);

  int dim_intern=0;
  js=new int[Ne];
  spin=new int[Ne];
  maxS=1;
  for(int i=0;i<Ne;i++) {js[i]=0;maxS*=2;}


  //state=new yosState[guessDim](Ne,Nm,js);
  nr=0;next=0;  
  while(!next) 
    {
      if(totJ(js,Nm)==reqJ)
	{
	  for(int i=0;i<maxS;i++)
	    {
	      totSz=0;
	      itmp=i;
	      for(int j=0;j<Ne;j++)
		{spin[j]=itmp%2;totSz+=itmp%2;itmp=itmp/2;}
	      totSz=2*totSz-Ne;
	      if(totSz!=reqSz) continue;
  	      if(!checkPauli_withSpin(js,spin)) continue;
	      if(dim_intern>=guessDim) 
		{
		  glLogger.warning("yosBasis::createAll_oneSz_oneJ");
		  glLogger.warning("guessDim exceeded. Increasing it.");
		  guessDim*=2;
		  tmpVec.resize(guessDim);
		}
	      tmpVec[dim_intern]=new yosState(Ne,Nm,js,spin);
	      tmpVec[dim_intern]->sPrint_num(tmpStr,tmpStr_len); 
	      glLogger.error("%ld: %s\n",dim_intern,tmpStr); 
	      dim_intern++;	      
	    }
	}
      next=nextSt(js,Nm);
    }
  if(dim_intern==0) 
    {
      glLogger.warning("yosBasis::createAll_oneSz_oneJ");
      glLogger.warning("found no vector matching the conditions. Note that totSz is 2*required tot. Sz, i.e. for 3 particles, totSz ought to be odd.");
    }
  //state_=new yosState*[dim_intern];
  setCapacity(dim_intern);
  for(int i=0;i<dim_intern;i++)
    state_[i]=tmpVec[i];
  delete [] tmpStr;

  dim = dim_intern;
  return dim_intern;
}

int yosBasis::nextSt(int *my_js)
{
  return nextSt(my_js,Nm_);
}

int yosBasis::nextSt(int *my_js,int Nm)
{
  int j;
  my_js[Ne_-1]++;
  if(my_js[Ne_-1]<Nm) return(0);
  j=Ne_-2;
  while(my_js[j]>=Nm-1 && j>=0) j--;
  if(j<0) {for(j=0;j<Ne_;my_js[j]=0,j++);return(1);}
  my_js[j]++;
  for(int i=j+1;i<Ne_;i++) 
    my_js[i]=my_js[j];
  return(0);
}

int yosBasis::totJ(int *my_js) 
{
  return totJ(my_js, Nm_);
}

int yosBasis::totJ(int *my_js,int Nm) 
{
  int sum=0;
  for(int j=0;j<Ne_;sum+=my_js[j],j++);
  return(sum % Nm);
}

inline int yosBasis::checkSymmetry(int *spat,int *spin,int i1,int i2)
{
  if(spat[i1]!=spat[i2]) return(0);
  if(spin[i1]!=spin[i2]) return(0);
  return(1);
}

int yosBasis::checkPauli_withSpin(int *spat,int *spin)
//returns: 0=excluded by Pauli principle; 1=allowed
{ 
  for(int i=0;i<Ne_-1;i++)
    {
      // states |...2+2-...> and |...2-2+...> lead to the same antisym. state
      if(spat[i]==spat[i+1]) if(spin[i]>=spin[i+1]) return(0);
      for(int j=i+1;j<Ne_;j++)
	if(checkSymmetry(spat,spin,i,j)) return(0);    
    }
  return(1);                            
}

int yosBasis::nextSt_nospin(int *my_js, int tmp_Pos, int tmp_Nm)
{
  my_js[tmp_Pos]++;
  if (my_js[tmp_Pos] < tmp_Nm) return(0);

  if (tmp_Pos == 0) return 1;

  // recursive call
  int notOk = nextSt_nospin (my_js, tmp_Pos-1, tmp_Nm-1); 

  int next = my_js[tmp_Pos-1] + 1;

  if (!notOk && next < tmp_Nm) 
    {
      my_js[tmp_Pos] = next;
      return 0;
    }
  else 
    return 1;
}

// End of methods related to genBasisForConst_J_Sz()

// Methods for inserting one state
/*! \fn int yosBasis::setCapacity(int new_cap)
  \brief Increase capacity of basis (not its dimension!). When inserting
  a new state to basis, its capacity has to be at least by 1 greater
  than its current dimension.
 */
int yosBasis::setCapacity(int new_cap)
{
  if(new_cap<dim) 
    {
      glLogger.error("yosBasis::setCapacity");
      glLogger.error("capacity<dim");
      exit(1);
    }
  if(new_cap==dim) return 0; // nothing to do

  yosState **p_old = state_;
  state_ = new yosState* [new_cap];
  for(int i=0;i<dim;i++)
    state_[i] = p_old[i];
  if(dim>0 && p_old != NULL) delete [] p_old;
  capacity = new_cap;
  return 0;
}



int yosBasis::adjustCapacityToSize()
{
  if(capacity==dim) return 0; // nothing to do

  yosState **p_old = state_;
  state_ = new yosState* [dim];
  for(int i=0;i<dim;i++)
    state_[i] = p_old[i];
  if(dim>0 && p_old != NULL) delete [] p_old;
  capacity = dim;
  return 0;
}

int yosBasis::addOneState(int *new_js,int *new_spins)
{
  if(dim == capacity)
    {
      glLogger.error("yosBasis::addOneState");
      glLogger.error("capacity == dim; increase capacity first");
      exit(1);
    }
  if( !spinYes )
    {
      glLogger.error("yosBasis::addOneState");
      glLogger.error("not implemented for spinYes==0");
      exit(1);
    }
  state_[dim] = new yosState(Ne_,Nm_,new_js,new_spins);
  dim++;
  return 0;
}



// methods for wavefunction-cache

void yosBasis::initWFCache (void)
{
  // std::cout << "initWFcache Nm=" << Nm_ << "\n";
  if (m_pWaveFct == NULL)
    m_pWaveFct = allocK< std::complex<double> >(Nm_);

  if (m_pWaveFctdx == NULL)
    m_pWaveFctdx = allocK< std::complex<double> >(Nm_);

  if (m_pWaveFctdy == NULL)
    m_pWaveFctdy = allocK< std::complex<double> >(Nm_);
}

void yosBasis::setWFfluxes (double alpha1, double alpha2)
{
  // tell all states
  for (int i=0; i<dimension(); i++)
    state_[i]->setSolenoidFluxes (alpha1, alpha2);
}

void yosBasis::setWFCacheForXY (double x, double y)
{
  assert (m_pWaveFct != NULL);   // if one of these assertions fails, you probably didn't call
  assert (m_pWaveFctdx != NULL); // initWFcache() before
  assert (m_pWaveFctdy != NULL);
  
  if (x != m_WF_x || y != m_WF_y)
    {
      //std::cout << " - calculating wavefunctions for (x=" << x << ", y=" << y << ")\n";
          
      for (int i=0; i<Nm_; i++)
	{
	  m_pWaveFct[i] = state_[0]->SPwaveFct (i, x, y);
	  m_pWaveFctdx[i] = state_[0]->SPwaveFctdx (i, x, y);
	  m_pWaveFctdy[i] = state_[0]->SPwaveFctdy (i, x, y);
	}

      m_WF_x = x;
      m_WF_y = y;
    }  
}

std::complex<double> yosBasis::getCachedWF (int j) const
{
  assert (j >= 0 && j < Nm_);
  return m_pWaveFct[j];
}

std::complex<double> yosBasis::getCachedWFdx (int j) const
{
  assert (j >= 0 && j < Nm_);
  return m_pWaveFctdx[j]; 
}

std::complex<double> yosBasis::getCachedWFdy (int j) const
{
  assert (j >= 0 && j < Nm_);
  return m_pWaveFctdy[j];
}

/*! \fn int yosBasis::shiftRight_allStates(int shift) 
  \brief Change (j_1,..,j_Ne) -> (j_1+shift,..,j_Ne+shift) for
  all states.
 */
int yosBasis::shiftRight_allStates(int shift)   
{
  for(int i=0;i<dim;i++)
    state_[i]->shiftRight_State(shift);
  return 0;
}



/*! \fn int test_totJ()
  \brief Tests whether all states of the basis have the same totJ and sets
  obs_totJ correspondingly.
  \return totJ if all the states have the same totJ, -1 otherwise.
*/
int yosBasis::test_totJ()
{
  int totJ=0,
	all_totJ=0;

  for(int i=0;i<dim;i++)
    {
      totJ = 0;
      for(int l=0;l<Ne_;l++)
	totJ += state_[i]->j(l);
      totJ = (totJ % Nm_);
      if(i==0) all_totJ = totJ;
      if(totJ != all_totJ)
	{
	  obs_totJ.knownInAdvance = false;
	  return -1;
	}
    }
  obs_totJ.knownInAdvance = true;
  obs_totJ.iVal = totJ;
  obs_totJ.dVal = totJ * 1.0;
  return totJ;
}

  /*! \fn int test_totSz()
    \brief Tests whether all states of the basis have the same totSz and sets
    obs_totSz correspondingly.
    \return totSz (twice the value, i.e. always integer) 
    if all the states have the same totSz, otherwise totSz of the last state tested.
   */
int yosBasis::test_totSz()
{
  int totSz=0,
all_totSz=0;

  for(int i=0;i<dim;i++)
    {
      totSz = 0;
      for(int l=0;l<Ne_;l++)
	totSz += state_[i]->spin(l);
      totSz = 2*totSz - Ne_;
      if(i==0) all_totSz = totSz;
      if(totSz != all_totSz)
	{
	  obs_totSz.knownInAdvance = false;
	  return totSz;
	}
    }
  obs_totSz.knownInAdvance = true;
  obs_totSz.iVal = totSz;
  obs_totSz.dVal = totSz * .5;
  return totSz;
}


int yosBasis::createSignatures()
{
  if(Nm_>32) 
    {
      glLogger.error("yosBasis::createSignature(): Nm>32 and int has only 32 bits; try long int");
      exit(1);
    }

  if (jsDoSignature != NULL) delete [] jsDoSignature;
  if (jsUpSignature != NULL) delete [] jsUpSignature;
  
  jsDoSignature = new int[dim];
  jsUpSignature = new int[dim];
  
  int tmp_j;

  for(int i=0;i<dim;i++)
    {
      jsDoSignature[i] = 0; jsUpSignature[i] = 0; 
      for(int l=0;l<Ne_;l++)
	{
	  tmp_j = 1 << (state_[i])->j(l); 
	  if((state_[i])->spin(l))
	    jsUpSignature[i] = jsUpSignature[i] | tmp_j;
	  else
	    jsDoSignature[i] = jsDoSignature[i] | tmp_j;
	}
    }
  return 0;
}


int yosBasis::shuffleSpins(int i,int k_from,int k_to)
{
  if(shuffled) return 0;
  if(k_from==k_to) return 0;
  int *new_spins=new int[Ne_];
  int *js=new int[Ne_];
  if(k_from < k_to)
    {
      int j=0;
      for(j=0;j<k_from;j++)
	new_spins[j]=(state_[i])->spin(j);
      for(int j=k_from;j<k_to;j++) 
	new_spins[j]=(state_[i])->spin(j+1);
      new_spins[k_to]=(state_[i])->spin(k_from);
      for(int j=k_to+1;j<Ne_;j++) 
	new_spins[j]=(state_[i])->spin(j);
    }
  else
    {
      int j=0;
      for(j=0;j<k_to;j++)
	new_spins[j]=(state_[i])->spin(j);
      new_spins[k_to]=(state_[i])->spin(k_from);
      for(int j=k_to+1;j<k_from+1;j++) 
	new_spins[j]=(state_[i])->spin(j-1);
      for(int j=k_from+1;j<Ne_;j++) 
	new_spins[j]=(state_[i])->spin(j);
    }
  for(int j=0;j<Ne_;j++)
    js[j]=(state_[i])->j(j);
  st_origUnshuffled=state_[i];
  state_[i]=new yosState(Ne_,Nm_,js,new_spins);
  shuffled=1;
  exit(1);
  delete []js;delete []new_spins;
  return 1;
}


int yosBasis::restoreSpins(int i)
{
  if(!shuffled) return 0;
  //state_[i]->~yosState();
  delete state_[i];
  state_[i]=st_origUnshuffled;
  shuffled=0;
  return 1;
}




bool yosBasis::readFromFile(FILE * in)
{
if (in) {
  findStartLine(in);
  Ne_=getdec(in);
  Nm_=getdec(in);
  dim=getdec(in);
  findStartLine(in);
  spinYes=1;
  myBasis=NULL;
  basis_type=YOS_BASIS;
  getBasisFromFile(in,dim);

  m_WF_x = -1.0;
  m_WF_y = -1.0;
  m_pWaveFct = NULL;
  m_pWaveFctdx = NULL;
  m_pWaveFctdy = NULL;
        return (true);
}
return (false);




}









