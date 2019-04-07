/* Some useful hand-made routines. */
/* Last updated 27.09.02.          */

#include<myaccessories.hpp>

void SeekLF(FILE* file)
{
  char c;
  while((c=getc(file))!='\n');
}

int findeof(FILE *file)
{
  char c;
  if(feof(file)) return(1);
  while((c=fgetc(file))<33 && !feof(file));
  if(feof(file)) return(1);
  ungetc(c,file);
  return(0);
}

int getdec(FILE *file)
{
  int n;
  fscanf(file,"%d",&n);     
  SeekLF(file);
  return(n);
}

double getdbl(FILE *file)
{
  double x;
  fscanf(file,"%lf",&x);     
  SeekLF(file);
  return(x);
}




/* Reads a string delimited by " from the file. From the actual position in
file searches first for the first ", then reads the string until the second
" is reached. */
void myFscanf(FILE *inpt,char *str)
{
  int n;
  char ch;
  str[0]=0;
  n=0;
  ch=fgetc(inpt);                 // 34 = "
  while((ch!=34)||feof(inpt)) ch=fgetc(inpt);  
  ch=fgetc(inpt);
  while(ch!=34)
    {      
      if(feof(inpt)) exit(1); 
      str[n]=ch;
      ch=fgetc(inpt);
      n++;
    }  
  str[n]=0;
  SeekLF(inpt);
}

/* x^y also for x=0 and/or y=0. */
double mypow(double x,double y)
{
  if(x == 0)
    {
    if(y==0)
      {
	return(1.0);
      }
      else 
	{
	  return (0.0);
	}
    }
  return(pow(x,y));
}



FILE *testfopen(char *name)
{
  FILE *file;
  file=fopen(name,"r");
  if(file==NULL)
    {
      printf("No file ");printf(name);printf(" so quitting.\n");
      exit(1);
    }         
  return(file);
}

FILE *testfopenW(char *name)
{
  FILE *file;
  file=fopen(name,"w");
  if(file==NULL)
    {
      printf("Could not open file ");
      printf(name);
      printf(" so quitting.\n");
      exit(1);
    }         
  return(file);
}


void findStartLine(FILE *inpt)
{
  char ch;
  ch=fgetc(inpt);                 // 35 = #
  while(ch==35) 
    {
      while((ch!='\n') && !feof(inpt)) ch=fgetc(inpt);
      ch=fgetc(inpt);
    }
  if(feof(inpt)) 
    { 
      printf("No uncommented lines in one of the input files. Quitting.\n ");
      exit(1);
    }  
  ungetc(ch,inpt);
}



/* The same as myFscanf(), but tests for max. string length maxn. */
void myMFscanf(FILE *inpt,char *str,int maxn)
{
  int n;
  char ch;
  str[0]=0;
  n=0;
  ch=fgetc(inpt);                 // 34 = "
  while((ch!=34)||feof(inpt)) ch=fgetc(inpt);  
  ch=fgetc(inpt);
  while(ch!=34)
    {      
      if(feof(inpt)) exit(1); 
      str[n]=ch;
      ch=fgetc(inpt);
      n++;
      if(n>=maxn) {printf("String too long. Exiting.\n");exit(1);}
    }  
  str[n]=0;
  SeekLF(inpt);
}

void Nbytes(long int n,char *str)
{
  double x;
  if(n<10000) sprintf(str,"%ld B",n);
  else{
    x=n/1024.0;
    if(x<1000) sprintf(str,"%5.1f K",x);
    else{
      x=x/1024;
      if(x<1000) sprintf(str,"%5.1f M",x);
      else{
	x=x/1024;
	if(x<1024) sprintf(str,"%5.1f G",x);
	else
	  sprintf(str,">1 T");
      }
    }
  }
}

void Nbytes(double x,char *str)
{
  long int xInt = (long int)(fabs(x));
  return Nbytes(xInt, str);
  /*
  if(x<10000) sprintf(str,"%4.0f B",x);
  else{
    x=x/1024.0;
    if(x<1000) sprintf(str,"%5.1f K",x);
    else{
      x=x/1024;
      if(x<1000) sprintf(str,"%5.1f M",x);
      else{
	x=x/1024;
	if(x<1024) sprintf(str,"%5.1f G",x);
	else
	  sprintf(str,">1 T");
      }
    }
  }
  */
}

void Ntime(double x,char *str)
{
  if(x<60) sprintf(str,"%2.0f s",x);
  else{
    x=x/60;
    if(x<60) sprintf(str,"%4.1f min",x);
    else{
      x=x/60;
      if(x<100) sprintf(str,"%4.1f h",x);
      else{
	x=x/24;
	sprintf(str,"%4.1f d",x);
      }
    }
  }
}

void Npercent(double x,FILE *f)
{
  static const int LENGTH=9;
  if(x<0.) {for(int i=0;i<LENGTH;i++) fprintf(f," ");fflush(f);return;}
  for(int i=0;i<LENGTH;i++) fprintf(f,"%c",8);
  if(x<=1.) {fprintf(f," %7.3f%%",100.*x);}
  fflush(f);
}






/*!
  \todo there is certainlay a mathlib function somewhere
*/

long factorial(long in)
{
  
  if (in < 0)
    {
      return -1;
    }
  if (in == 0)
    {
      return (1);
    }
  else
    {
      return (in * factorial(in-1));
    }
}


