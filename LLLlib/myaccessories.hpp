/* Some useful routines...
 */

#ifndef MYACCESSORIES_HPP
#define MYACCESSORIES_HPP


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include<new>

void SeekLF(FILE* file);  // go to the EOL in file
int findeof(FILE *file);  // something???
int getdec(FILE *file);   /* read an integer from file and 
			     go to next line (GTNL)*/
double getdbl(FILE *file);// read a double from file and GTNL
void myFscanf(FILE *inpt,char *str); /* read a string from file 
					(see .cpp file) and GTNL */
double mypow(double x,double y); // x^y, 0^0=1, 0^x=0.
FILE *testfopen(char *name);     /* open file for reading and tell
				    if the file does not exist */
FILE *testfopenW(char *name);    /* open file for writing and tell
				    if the file cannot be opened (e.g.
				    when a directory doesnot exist). */
void findStartLine(FILE *inpt);  /* find a new line which is not 
				    a comment line (starting with #)*/
void myMFscanf(FILE *inpt,char *str,int maxn);
                          /* as myFscanf but watches for the string 
			     length */
void Nbytes(long int n,char *str); /* n=number of bytes, convert it
	     into appropriate unit (K,M,G) and put it into string */
void Nbytes(double x,char *str); // if n is too large, turn it to dbl
void Ntime(double x,char *str);  // like Nbytes but with time in seconds
void Npercent(double x,FILE *f); 

long factorial (long);

// Templates defined right away...

/*!
  Allocate an array of 'size' objects of the type Num_Type. If there's
   not enough memory, say it. 
   Thats something ypou really dont need. new does throw an error anyways.
*/
template<class Num_Type> 
Num_Type *allocK(long int size)
{
  Num_Type *tmp;
  char stmp[10];

  if((tmp=new(std::nothrow) Num_Type[size]))
    {
      return tmp;
    }
  else 
    {
      Nbytes((1.*size)*sizeof(Num_Type),stmp);
      printf("Not enough memory to allocate an array, %s would be needed.\n",stmp);
      exit(1);
    }
  return 0;
}

#endif
