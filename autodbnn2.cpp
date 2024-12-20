/***************************************************************************
                          main.cpp  -  description
                             -------------------
    begin                : Sat Oct 13 11:58:21 IST 2001
    copyright            : (C) 2001 by Ninan Sajeeth Philip
    email                : nsp@stthom.ernet.in
			 : nspp@iucaa.ernet.in
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

/*
November 2010 Version 6
Most part rewritten. Algorithm changed from compensation to actual bin separation in feature 
space. Expected to be run with par files.
Jan 2010; Version 5.1
mask_dist added for each vector along with min max masks. This gives autodbnn the ability to 
handle chessboard type datasets where max,min masks may not work.Factor Analysis removed.
Replaced variable declarations with vector so that dynamic memory allocations are possible.
Nov. 2003 : Version 4.2
Factor Analysis included. The significance of each factor is saved in a 
file <filename>.factor
June 2003 : Version 4.1
DBNN can now be run in automated mode using par files. This is useful for automated
selection of training vectors from a sample of human classified objects. See related
publication.
April 2002 : Version 4.01
DBNN can now determine the number of bins on its own. Array structures changed. You can have
more number of feature vectors and output classes.
Nov. 2001 : Version 2.2
Sept.20th 2001
This is version 2.1 of the DBNN software. The major change here is that a network topology is
considered better over the other by computing the total probability of success when there in no
increase in the number of success counts. rslt orslt prslt changed to floats from integer mode.
A new variable pcnt and pocnt now holds the success counts. The second change is in initializing
the arrays - new c++ in RedHat 7.x require this to avoid segmentation error!

Program to generate anticipation inputs for the training/test set(c)Sajith Philip 1999.
This program is Public Domain. You have all rights to modify and use it. But please acknowledge
the original author and the cite that contains the details of this network.
You may download a technical paper on this network from:
                    
                        http://www.iucaa.ernet.in/~nspp

To run the program, you need an information file and a training and test data set. The format
of these files are explained later in this file.

Enjoy!

Sajith Philip
*/
//**********************************************************
// Command to compile:   c++  -lm -O3  main.cpp -o autodbnn
//**********************************************************

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
using namespace std;
#include <stdlib.h>
#include<sys/times.h> // times() fun. is here.
#include <time.h>
#include <vector>
using std::vector;
/**************************
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

#include "vsipc.h"
***************************/
#define max_resol 1600
#define features 100
#define classes 500           // If you get segmentation errors, reduce these values to fit your
//#define oneround 100       //  Memory size.
//#define fst_gain 1.0 Moved to be a floating variable
static double bgain,gain,dmyclass[classes],classval[classes],cmax,c2max,tmp2_wts,totprob,oldj;
static double LoC=0.65;
static double nLoC=0.0;
static int resol=100,nresol=0,nLoCcnt=0,skpchk=0;
static double omax,omin,rslt,rslt2,orslt,orslt2,prslt,nrslt,fst_gain;
clock_t start,stop;
static int argfnd, oneround=100,kmax,k2max,ans1,rnn,rnd,i,j,k,l,n,nd1,m,c1cnt,c2cnt,pcnt,pocnt,invcnt,innodes=100,outnodes=100;
char fln[256],fltmp[256],urchoice,urchoicex,bcchoice,savedpar;
FILE *fl1,*fl2,*fl3,*fl4,*fl5,*fl6,*fl7,*fl8,*fl9,*fl10;
int main(int argv, char *argp[256])
/*
 Important note for revision in ver. 4.02
 Compute the bin center of gravities and the Kernel fit that holds the probability to
 find them. The slope of it should be used as bgain for each bin.
*/
/*
 Revision in June 2003
 You can now run dbnn in automated mode by specifying the parameters in 0.par and 1.par
 files. Also dbnn can now use the bin values from the saved apf file.
*/
{
   if(argv > 3)
   {
    argfnd=1;
    cout << "The selected option is " << *argp[3] <<"\n";    
    switch(*argp[3])
    {
    case '0':
    
    ans1=0;
    if((fl2=fopen("0.par","r"))!=NULL)
     {
      fscanf(fl2,"%c\n",&bcchoice);
      fscanf(fl2,"%c\n",&urchoice);
      fscanf(fl2,"%c\n",&savedpar);
      fscanf(fl2,"%c\n",&urchoicex);
      fclose(fl2);			
      }
      else
     { cout << "No Parameter File... existing..";
     exit(1);
     }
    break;
    case '1':

    
    ans1=1;
    if((fl2=fopen("1.par","r"))!=NULL)
     {
      fscanf(fl2,"%lf",&gain);
      fscanf(fl2,"%d",&oneround);
      fclose(fl2);			
     }
     else
     { cout << "No Parameter File... existing..";
     exit(1);
     }
    break;
    case '2':
    
    ans1=2;
    break;

    case '3':
    
    ans1=3;
    break;
    
    default:

     cout << "Create the APF file(0) or Create the Weights file (1) or Classify Data(2,3) ?";
     cin >> ans1;
     break;
    }
    }
    else
    {
     argfnd=0;
     cout << "Create the APF file(0) or Create the Weights file (1) or Classify Data(2,3) ?";
     cin >> ans1;
    }
    if(ans1 == 2)
     {
     if(argfnd==1)
     bgain=0.0;
     else
     {
     cout << "Allowed relaxation on the boundary (in % use 0 for default from training data) :";
     cin >> bgain;
//     if (bgain < 0) bgain = 0;
     bgain=bgain*1.0;
     }
     }
     else
     bgain= 0;  // During training we are strict on boundary constraints.

    if(argv < 3)
     {
     cout << "Enter the name of the input file without extension (dat) :";
     cin >> fln;
     //scanf ("%s",&fln);
     }
     else
     {
     strcpy(fln,argp[1]);
     }
     strcpy(fltmp,fln);
     strcat(fltmp,".dat");
/*
  The structure of the data file is:
  Feature1 Feature2 Feature3 ....(etc upto innodes) ActualClass
  Feature1 Feature2 Feature3 ....(etc upto innodes) ActualClass
  Feature1 Feature2 Feature3 ....(etc upto innodes) ActualClass
  The delimiters are spaces and not tabs!!
*/
    if((fl1=fopen(fltmp,"r"))!=NULL)
     {
	    strcpy(fltmp,fln);
	    strcat(fltmp,".inf");
/*
  The format of the info file is: (in each line enter)
  innodes
  outnodes
  margin   <- This addition is required for regression problems.
  1.0       <- You can give any real positive value here. It is just a label.
  2.0
  ... (etc. upto no of classes)
  0.65 <- The Margin or Line of Control for marginal values.
  100 <- By default, the maximum bin size is set to 100. You can change this if required.
*/	
     	if((fl2=fopen(fltmp,"r"))!=NULL)
	    {
		 i=0;
		 fscanf(fl2,"%d",&innodes);
		 fscanf(fl2,"%d",&outnodes);
		 for (i=0;i<=outnodes;i++) // dmyclass[0] contains margin others are expected values.
		 fscanf(fl2,"%lf",&dmyclass[i]);
		 fscanf(fl2,"%lf",&LoC);   // New parameter to specify the Line Of Control
		 fscanf(fl2,"%d",&nresol);
		 cout <<"You have "<< innodes << " input nodes and " << outnodes <<" Output nodes with " << "margin set to " << LoC << "\n";
		 cout << "The target outputs are\n";
		 for (i=0;i<=outnodes;i++) cout << dmyclass[i] <<"\n";
		 if(nresol >0) 
		 {
		 resol=nresol;cout << "The maximum binsize is: " << resol <<"\n";
		 }	
		 else
		 {cout << "The maximum binsize is: " << resol<<"\n";
		 }
		 fst_gain=1.0/outnodes;
	    }
	    else
	    {
		 cout << "Unable to find the Info file. Exiting !!";
		 exit(1);
	    }
  
     } // program ends.
     else   // data file read error.
     {
      cout << "Unable to open the data file";
      exit(1);
      }
      gain=2.0;
      cout << "Going to initialise the arrays\n";
 /**************** Let us Define the Network Structure *********************************/
 double(* TabLookUp);
 TabLookUp = new double[int(max_resol)+4];
double mask_disp_maxres; // Space to save max resol for normalisation of mask_dist
double vects[innodes+outnodes+2],vectso[innodes+outnodes+2],tmpv,max[innodes+2],min[innodes+2];
if(TabLookUp==NULL){cout << "Out of Memory to Run Code at TabLookUp.. Exiting\n";exit(1);}
vector<vector<vector<vector<vector<int> > > > > anti_net(innodes+2,vector<vector<vector<vector<int> > > >(resol+2, vector<vector<vector<int> > >(innodes+2,vector<vector<int> >(resol+2,vector<int>(outnodes+2,0)))));
vector<vector<vector<vector<vector<double> > > > > anti_wts(innodes+2,vector<vector<vector<vector<double> > > >(resol+2, vector<vector<vector<double> > >(innodes+2,vector<vector<double> >(resol+2,vector<double>(outnodes+2,0)))));
vector<vector<vector<vector<vector<double> > > > > antit_wts(innodes+2,vector<vector<vector<vector<double> > > >(resol+2, vector<vector<vector<double> > >(innodes+2,vector<vector<double> >(resol+2,vector<double>(outnodes+2,0)))));
vector<vector<vector<vector<vector<double> > > > > antip_wts(innodes+2,vector<vector<vector<vector<double> > > >(resol+2, vector<vector<vector<double> > >(innodes+2,vector<vector<double> >(resol+2,vector<double>(outnodes+2,0)))));
 int resolution[innodes+8];
 double classtot[innodes+2][resol+2];           // Total Prob. computed
 if(classtot==NULL){cout << "Out of Memory to Run Code at classtot.. Exiting\n";exit(1);}
 double binloc[innodes+2][resol+8];
 if(binloc==NULL){cout << "Out of Memory to Run Code at binloc.. Exiting\n";exit(1);}

cout << "Going to format the allocated space\n";
  /***************************Let us put up the Network***********************************/
//    Start the counter for case 2 here.................
  start = times(NULL);
   if (ans1==0)
   {
	      n=0;
	      omax=-400;
	      omin=400;
	      while (!feof(fl1))
	      {
		 for(i=1;i<=innodes;i++)
		 if (n==0)
		 {
		   fscanf(fl1,"%lf",&vects[i]);
		   min[i]=vects[i];
		   max[i]=vects[i];
		 }
		 else
		 {
		   fscanf(fl1,"%lf",&vects[i]);
		   if( vects[i]> max[i]) max[i]=vects[i];
		   if (min[i] > vects[i]) min[i]=vects[i];
		 }
		 fscanf(fl1,"%lf\n",&tmpv);
		 if(tmpv>omax) omax = tmpv;
                 if(tmpv<omin) omin =tmpv;
		 k=1;
		 j=1;
		 n++;
	      }
              
	      cout << "No of vectors =" << n <<" and i/n is= " << 1.0/n << "\n";

        if(argfnd==0)
              {
	      cout <<"Do you want to use the saved parameters (Y/N)? ";
	      cin >>savedpar;
              }
	      if (savedpar == 'y') savedpar='Y';
	      else
	      if(savedpar == 'n') savedpar='N';
	      if((savedpar == 'Y') || (savedpar=='y'))
	      {
	       strcpy(fltmp,fln);
	       strcat(fltmp,".apf");
	       fl2=NULL;
	       if((fl2=fopen(fltmp,"r"))!=NULL)
	        {
	        cout << "Reading from the saved information\n";
	   	for (i=1;i<=innodes;i++)
		{
		fscanf(fl2,"%d",&resolution[i]);
                for(j=0;j<=resolution[i];j++) binloc[i][j+1]=j*1.0;
                
//		cout<<". ";
	        }
		cout << innodes << " itemes read from " << fltmp <<"\n";
		}
		else
		{
		cout << "ERROR: File " << fltmp << " not found" << "\n";
		exit(1);
		}
	//	fclose(fl2);
	      }
	      else
	      for(i=1;i<=innodes;i++)
	      {
                cout << "Enter the resolution required for node " << i << "[" << resol << "] (1 to " << max[i]-min[i] << "): ";
                cin >> resolution[i];
                for(j=0;j<=resolution[i];j++) binloc[i][j+1]=j*1.0;
	      }
              for(k=1;k<=outnodes;k++)
	      for(i=1;i<=innodes;i++)
	      for(j=0;j<=resolution[i];j++)
	      for(l=1;l<=innodes;l++)
	      for(m=0;m<=resolution[l];m++)
	      {
//	      cout << "[i,j,k]= [" << i << ","<< j <<","<< k <<"]\n"; 
	      anti_net[i][j][l][m][k]=1;
	      anti_wts[i][j][l][m][k]=(double)(1.0);
	      }
	      // Start the counter now...............
	      start = times(NULL);
           
	     rewind(fl1);
//	     cout << "I am here \n";
	      while (!feof(fl1))
	      {
	        for (i=1;i<=innodes;i++) fscanf(fl1,"%lf",&vects[i]);
	 	fscanf(fl1,"%lf\n",&tmpv);
	 	for(i=1;i<=innodes;i++)
                {
                vectso[i]=vects[i];
	 	vects[i]=round((vects[i]-min[i])/(max[i]-min[i])*resolution[i]);
                 }
	 	for (i=1;i<=innodes;i++)
	 	for (l=1;l<=innodes;l++)
	 	{
		  j=0;
		  m=0;
	         oldj=(double)2*resolution[i];
                  while (fabs(vects[i]-binloc[i][j+1]) < oldj)
                  {
                  oldj=fabs(vects[i]-binloc[i][j+1]);
                  j++;
                  }
                  if(j >0)j--;
		  oldj=(double)2*resolution[l];
                  while (fabs(vects[l]-binloc[l][m+1]) < oldj)
                  {
                  oldj=fabs(vects[l]-binloc[l][m+1]);
                  m++;
                  }
                  if(m >0)m--;
		  k=1;
		  while ((fabs(tmpv - dmyclass[k])) > dmyclass[0]) k++;
		   (anti_net[i][j][l][m][k])++;
		   (anti_net[i][j][l][m][0])++;
        	  }
	      }
	      fclose(fl1);
	      fclose(fl2);
   
         /*
            The conditional Probability,
	    P(A|B) = P(A intersection B)/P(B) is the
	    probability for the occurance of A(k) if B(ij) has happened =
	    Share of B(ij) that is held by A(k) / Probability of total B(ij)
	    in that particular feature i with resolution j.

                      */
              strcpy(fltmp,fln);
	      strcat(fltmp,".awf");      // This file holds the weights
	      fl6=fopen(fltmp,"w+");
	      strcpy(fltmp,fln);
	      strcat(fltmp,".apf");     // This file holds the estimated probability
	      if((fl1=fopen(fltmp,"w+"))!=NULL)
	      {
		 for(i=1;i<=innodes;i++) fprintf(fl1,"%d ",resolution[i]);
		 fprintf(fl1,"\n%lf %lf \n",omax,omin);
		 for(i=1;i<=innodes;i++) fprintf(fl1,"%lf ",max[i]);
		 fprintf(fl1,"\n");
		 for(i=1;i<=innodes;i++) fprintf(fl1,"%lf ",min[i]);
		 fprintf(fl1,"\n");
		 for(k=1;k<=outnodes;k++)
                 {
		   for(i=1;i<=innodes;i++)
		   for(j=0;j<=resolution[i];j++)
                   {
		    for(l=1;l<=innodes;l++)
		    for(m=0;m<=resolution[l];m++)
                    {
//cout << (double)anti_wts[i][j][l][m][k] << " ";
                      fprintf(fl1,"%d ",anti_net[i][j][l][m][k]);
                      fprintf(fl6,"%lf ",(double)anti_wts[i][j][l][m][k]);
  		    }
		   fprintf(fl6,"\n");
		   fprintf(fl1,"\n");
		   }
		   fprintf(fl6,"\n");
	           fprintf(fl1,"\n");
                 }
                 fprintf(fl6,"\n");
	         fprintf(fl1,"\n");
              }
	      else
	      {
		 cout << "Unable to create file for output\n";
		 exit(1);
	      }
              for(i=1;i<=innodes;i++)
	      for(j=1;j<=resolution[i];j++)
	      fprintf(fl6,"%lf\n", (double)binloc[i][j]);                 /// Let us print the bins.
	      fclose(fl1);
	      fclose(fl6);
	      fflush(NULL);
	      cout << "Creating the Anticipated Weights data file\n";
}
/**********************************End of Case 0 ******************************/
if(ans1==1)
{
    start = times(NULL);
    pcnt=0;
    pocnt=0;
    rslt=0.0;
    rslt2=0.0;
    orslt=rslt;
    orslt2=rslt2;
    cout << "The programe will now modify the compensatory weights\n";
    if(argfnd==0)
    {
    cout << "Please enter the gain:";
    cin >> gain;
    cout << "Please enter the number of training epochs:";
    cin >> oneround;
    }
    // Start the counter in this round here...................
    start = times(NULL);
    strcpy(fltmp,fln);
    strcat(fltmp,".awf");
    if((fl6=fopen(fltmp,"r"))!=NULL)
    {
//     fl6=fopen(fltmp,"r");
     strcpy(fltmp,fln);
     strcat(fltmp,".apf");
     fl2=NULL;
     if((fl2=fopen(fltmp,"r"))!=NULL)
       {
	for (i=1;i<=innodes;i++) fscanf(fl2,"%d",&resolution[i]);
	fscanf(fl2,"\n%lf",&omax);
	fscanf(fl2,"%lf",&omin);
	fscanf(fl2,"\n");
	for(i=1;i<=innodes;i++) fscanf(fl2,"%lf",&max[i]);
	fscanf(fl2,"\n");
	for(i=1;i<=innodes;i++) fscanf(fl2,"%lf",&min[i]);
	fscanf(fl2,"\n");
        for(k=1;k<=outnodes;k++)
        {
         for(i=1;i<=innodes;i++)
	 for(j=0;j<=resolution[i];j++)
	 {
	  for(l=1;l<=innodes;l++)
	  for(m=0;m<=resolution[l];m++)
	  {
	   fscanf(fl2,"%d",&anti_net[i][j][l][m][k]);
	   anti_net[i][j][l][m][0]+=anti_net[i][j][l][m][k];
	   fscanf(fl6,"%lf",&anti_wts[i][j][l][m][k]);
	   antit_wts[i][j][l][m][k]=anti_wts[i][j][l][m][k];
	   antip_wts[i][j][l][m][k]=anti_wts[i][j][l][m][k];
	  }
	   fscanf(fl2,"\n");
	   fscanf(fl6,"\n");
	 }
	  fscanf(fl2,"\n");
	  fscanf(fl6,"\n");
	 }
  	 for(i=1;i<=innodes;i++)
	 for(j=1;j<=resolution[i];j++)
	 fscanf(fl6,"%lf\n", &binloc[i][j]);                 /// Let us print the bins.
        }
        else
	{
 	cout << "Unable to Open the APF information file\n";
	exit(1);
	}
 	fclose(fl2);
      }
      else
      {
 	cout << "Unable to Open the AWF information file\n";
	exit(1);
      }
      fclose(fl6);
      
      for(rnd=0;rnd<=oneround;rnd++)     // Training round starts here....
      {
       if((n==pocnt)&& (n>0))  break;
       strcpy(fltmp,fln);
       strcat(fltmp,".dat");
       fl1=fopen(fltmp,"r");
  	   n=0;
       rslt=0.0;
       rslt2=0.0;
       pcnt=0;
     	   while (!feof(fl1))
	   {
              for(k=1;k<=outnodes;k++) classval[k]=1.0;
 	      n++;
	      if(ans1==3)
	      {
	      for (i=1;i<innodes;i++) fscanf(fl1,"%lf",&vects[i]);
	      fscanf(fl1,"%lf\n",&vects[innodes]);
	      }
	      else
	      {
	      for (i=1;i<=innodes;i++) fscanf(fl1,"%lf",&vects[i]);
	      fscanf(fl1,"%lf\n",&tmpv);
	      }
	      for(i=1;i<=innodes;i++)
              {
              vectso[i]=vects[i];
	      vects[i]=round((vects[i]-min[i])/(max[i]-min[i])*resolution[i]);
              }
	      for (i=1;i<=innodes;i++)
	      {
		j=0;
	        k=1;
	        oldj=(double)2*resolution[i];
                while (fabs(vects[i]-binloc[i][j+1]) < oldj)
                {
                 oldj=fabs(vects[i]-binloc[i][j+1]);
                 j++;
                }
                if(j>0)j--;
	      for (l=1;l<=innodes;l++)
	      {
		m=0;
	        k=1;
	        oldj=(double)2*resolution[l];
                while (fabs(vects[l]-binloc[l][m+1]) < oldj)
                {
                 oldj=fabs(vects[l]-binloc[l][m+1]);
                 m++;
                }
                if(m>0)m--;
	        for (k=1;k<=outnodes;k++)
		{
                 if(anti_net[i][j][l][m][0] > 0)
                 tmp2_wts=(double)anti_net[i][j][l][m][k]/(anti_net[i][j][l][m][0]);
                 else
 		 tmp2_wts=1/outnodes;
	         classval[k]*=(double)tmp2_wts*anti_wts[i][j][l][m][k];
	         }
	        } 

	       }
	       kmax=1;
	       cmax=0;
	       for (k=1;k<=outnodes;k++)
	       {
	        if (classval[k] > cmax)
	        {
	         cmax=classval[k];
	         kmax=k;
	        }
	       }
	       if ((fabs(dmyclass[kmax]-tmpv) > dmyclass[0]) && (rnd >0))
	       {
	        for (i=1;i<=innodes;i++)
	        {
	         j=0;
	         k=1;
                 oldj=(double)2*resolution[i];
                 while (fabs(vects[i]-binloc[i][j+1]) < oldj)
                 {
                  oldj=fabs(vects[i]-binloc[i][j+1]);
                  j++;
                 }
                 if(j>0)j--;
	         for (l=1;l<=innodes;l++)
	         {
		  m=0;
	          k=1;
	          oldj=(double)2*resolution[l];
                  while (fabs(vects[l]-binloc[l][m+1]) < oldj)
                 {
                  oldj=fabs(vects[l]-binloc[l][m+1]);
                  m++;
                 }
                 if(m>0)m--;
	         while (fabs(dmyclass[k]-tmpv) > dmyclass[0]) k++;
//	         {
                  if((classval[(int)kmax] >0)&&(classval[k]<classval[(int)kmax]))
	          anti_wts[i][j][l][m][k]+=(double)gain*(1.0-(classval[k]/classval[(int)kmax]));
      		  if(anti_wts[i][j][l][m][k] <= 0.0)
		  cout << k << " "<< tmpv << "[" << dmyclass[1] << "]" << dmyclass[outnodes] << "\n";
//                 }
		 }	
	        }
	       } // kmax che
	  } // while not eof check
		 // Now save the wieights
	 fclose(fl1);
	 strcpy(fltmp,fln);
         strcat(fltmp,".dat");
         fl1=fopen(fltmp,"r");
         m=n;
  	 n=0;
         rslt=0.0;
         rslt2=0.0;
         pcnt=0;
         while (!feof(fl1))                    // Test round...
	 {
         for(k=1;k<=outnodes;k++) classval[k]=1.0;
 	   n++;
	   kmax=1;
	   cmax=0;
	   if(ans1==3)
	   {
	    for (i=1;i<innodes;i++) fscanf(fl1,"%lf",&vects[i]);
	    fscanf(fl1,"%lf\n",&vects[innodes]);
	   }
	   else
	   {
	    for (i=1;i<=innodes;i++) fscanf(fl1,"%lf",&vects[i]);
	    fscanf(fl1,"%lf\n",&tmpv);
	   }
	   for(i=1;i<=innodes;i++)
           {
	    vects[i]=round((vects[i]-min[i])/(max[i]-min[i])*resolution[i]);
            if (vects[i] < 0) vects[i]=0;             // let us be bounded. #Oct 2001.
           }
	   for (i=1;i<=innodes;i++)
	   {
	     j=0;
	     k=1;
             oldj=(double)2*resolution[i];
             while (fabs(vects[i]-binloc[i][j+1]) < oldj)
             {
              oldj=fabs(vects[i]-binloc[i][j+1]);
              j++;
             }
             if(j>0)j--;
	     for (l=1;l<=innodes;l++)
	     {
		m=0;
	        k=1;
	        oldj=(double)2*resolution[l];
                while (fabs(vects[l]-binloc[l][m+1]) < oldj)
                {
                 oldj=fabs(vects[l]-binloc[l][m+1]);
                 m++;
                }
                if(m>0)m--;
	        for (k=1;k<=outnodes;k++)
		{
                 if(anti_net[i][j][l][m][0] > 0)
                 tmp2_wts=(double)anti_net[i][j][l][m][k]/(anti_net[i][j][l][m][0]);
                 else
 		 tmp2_wts=1/outnodes;
	         classval[k]*=(double)tmp2_wts*anti_wts[i][j][l][m][k];
	        }
	     }
	    
	   }
	   for (k=1;k<=outnodes;k++)
	   {
	    if (classval[k] > cmax)
	    {
	     cmax=classval[k];
	     kmax=k;
	    }
           }
           if (fabs(dmyclass[kmax]-tmpv) <= dmyclass[0])
	   {
	    rslt2+=cmax;
	    pcnt++;
	   }
	   else
	   {
	    k=1;
	    while (fabs(dmyclass[k]-tmpv) > dmyclass[0]) k++;
	    rslt+=cmax-classval[k];
	   }
	 } // while not eof check
	 // Now save the wieights
	 fclose(fl1);
	 kmax=1;
         if(orslt2==0) orslt2=rslt2;
         if(orslt==0) orslt=rslt;
	 prslt=(rslt2-orslt2);
         if(rslt > 0)
	 nrslt=(orslt/rslt);
	 if(pcnt>pocnt) 
	 {
	   rnn=rnd;
	   pocnt=pcnt;   // The best result is now saved in pocnt
	   for(k=1;k<=outnodes;k++)
           for(i=1;i<=innodes;i++)
	   for(j=0;j<=resolution[i];j++)
	   for(l=1;l<=innodes;l++)
	   for(m=0;m<=resolution[l];m++)
           {
	     antip_wts[i][j][l][m][k]=antit_wts[i][j][l][m][k];
	     antit_wts[i][j][l][m][k]=anti_wts[i][j][l][m][k];
  	   }
	   cout << "Round:" << rnn << "| TProb["<<prslt<<"," <<nrslt<<"] | Passed count:" << pocnt << endl;
	   if(orslt2 <rslt2) orslt2=rslt2;
	   if(rslt < orslt) orslt=rslt;
	 }
	 n=m;
         // prslt=rslt+rslt2;
	}  //rnd inc.
	strcpy(fltmp,fln);
	strcat(fltmp,".awf");
	fl6=fopen(fltmp,"w+");
	kmax=1;
        for(k=1;k<=outnodes;k++)
        {
         for(i=1;i<=innodes;i++)
         for(j=0;j<=resolution[i];j++)
         {
          for(l=1;l<=innodes;l++)
          for(m=0;m<=resolution[l];m++)
          {
           fprintf(fl6,"%lf ",antit_wts[i][j][l][m][k]);
          }
          fprintf(fl6,"\n");
         }
         fprintf(fl6,"\n");
        }
        fprintf(fl6,"\n");
//cout <<"I am here...";
        for(i=1;i<=innodes;i++)
        for(j=1;j<=resolution[i];j++)
        fprintf(fl6,"%lf\n", binloc[i][j]);                 /// Let us print the bins.
 	fflush(fl6);
	fclose(fl6);
	fl6=NULL;
	cout << "Best result at round " << rnn<< endl;
      }  // ans <> 1
/***********************************End of Case 1*******************************/	
      strcpy(fltmp,fln);
      strcat(fltmp,".dat");
      fl1=fopen(fltmp,"r");
      strcpy(fltmp,fln);
      strcat(fltmp,".awf");
      fl6=NULL;
      fl6=fopen(fltmp,"r");
      strcpy(fltmp,fln);
      strcat(fltmp,".apf");
      fl2=NULL;
      if((fl2=fopen(fltmp,"r"))!=NULL)
      {
        cout << "Creating the Anticipated Network outputs\n";
	for (i=1;i<=innodes;i++) fscanf(fl2,"%d\n",&resolution[i]);
	fscanf(fl2,"%lf",&omax);
	fscanf(fl2,"%lf",&omin);
	fscanf(fl2,"\n");
        for(i=1;i<=innodes;i++) fscanf(fl2,"%lf",&max[i]);
        fscanf(fl2,"\n");
    	for(i=1;i<=innodes;i++) fscanf(fl2,"%lf",&min[i]);
 	fscanf(fl2,"\n");
        for(k=1;k<=outnodes;k++)
        {
         for(i=1;i<=innodes;i++)
	 for(j=0;j<=resolution[i];j++)
	 {
	  for(l=1;l<=innodes;l++)
	  for(m=0;m<=resolution[l];m++)
	  {
	   fscanf(fl2,"%d",&anti_net[i][j][l][m][k]);
	//   anti_net[i][j][l][m][0]+=anti_net[i][j][l][m][k];
	   fscanf(fl6,"%lf",&anti_wts[i][j][l][m][k]);
	   antit_wts[i][j][l][m][k]=anti_wts[i][j][l][m][k];
	   antip_wts[i][j][l][m][k]=anti_wts[i][j][l][m][k];
	  }
	   fscanf(fl2,"\n");
	   fscanf(fl6,"\n");
	 }
	  fscanf(fl2,"\n");
	  fscanf(fl6,"\n");
	 }
       }
       else
       {
  	cout << "Unable to Open the APF information file";
	exit(1);
       }
       for(i=1;i<=innodes;i++)
       for(j=1;j<=resolution[i];j++)
       {
        fscanf(fl6,"%lf\n",&binloc[i][j]);                 /// Let us print the bins.
       }
       fclose(fl6);
       fl4=fopen("output.dat","w+");  // Network Output values
       cout << "Read all input parameters\n";
// *********** case 3 ***********************************************
      if (ans1 !=3)
      {
       fl5=fopen("actual.dat","w+");  // Expected Output Values
       strcpy(fltmp,fln);
       strcat(fltmp,argp[2]);
       strcpy(fltmp,fln);
       strcat(fltmp,argp[2]);
       strcat(fltmp,".cmp");         // Lets see how well the classification went.
       fl7=fopen(fltmp,"w+");
       fprintf(fl7,"Sample     Predicted   Predicted    Actual            Prediction \n");
       fprintf(fl7," No.       Ist Choice  2nd Choice   item              Confidence\n");
       c1cnt=0;
       c2cnt=0;
       invcnt=0;
       n=0;
      }
 // Create classtot values ***********************
      for(i=1;i<=innodes;i++)for(j=0;j<=resolution[i];j++) for(l=1;l<=innodes;l++)for(m=0;m<=resolution[l];m++) anti_net[i][j][l][m][0] =0;
      for(k=1;k<=outnodes;k++)
      for(i=1;i<=innodes;i++)
      for(j=0;j<=resolution[i];j++)
      for(l=1;l<=innodes;l++)
      for(m=0;m<=resolution[l];m++)
      {
      anti_net[i][j][l][m][0]+=(double)(anti_net[i][j][l][m][k]);
//      cout << "[" << (anti_net[i][j][l][m][k]) <<"]";
      }
     while (!feof(fl1))
	 {
          for(k=1;k<=outnodes;k++) classval[k]=1.0;
          n++;
	  cmax= 0.0;
	  c2max=0.0;
	  kmax=0;
	  k2max=0;
          classval[0]=0.0;
	      if(ans1==3)
	      {
               
	      for (i=1;i<innodes;i++) fscanf(fl1,"%lf",&vects[i]);
	      fscanf(fl1,"%lf\n",&vects[innodes]);
	      }
	      else
	      {
	      for (i=1;i<=innodes;i++) fscanf(fl1,"%lf",&vects[i]);
	      fscanf(fl1,"%lf\n",&tmpv);
	      }
              skpchk=0;
	      for(i=1;i<=innodes;i++)
               {
               vectso[i]=vects[i]; 
               if((max[i]-min[i]) >0) 
               vects[i]=round(((vects[i]-min[i])/(max[i]-min[i]))*resolution[i]);
               else
               skpchk=1;
               if ( resolution[i] < vects[i]) 
               {
//               cout << "Data on line " << n << " has vect[" << i << "]= "<< vects[i] << " while resolution[i] = " << resolution[i] <<" is problematic \n";
	       skpchk=1;
               }
               }
               if(skpchk==0)
               {
		 for (i=1;i<=innodes;i++)
		 {
         	     j=0;
	             k=1;
	             oldj=(double)2*resolution[i];
                     while ((fabs(vects[i]-binloc[i][j+1]) < oldj)&& (j<= resolution[i]+1))
                     {
                       oldj=fabs(vects[i]-binloc[i][j+1]);
                       j++;
                     }
                     if(j >0)j--;
	             for (l=1;l<=innodes;l++)
	             {
		      m=0;
	              k=1;
	              oldj=(double)2*resolution[l];
                      while ((fabs(vects[l]-binloc[l][m+1]) < oldj)&& (m<= resolution[l]+1))
                      {       
                        oldj=fabs(vects[l]-binloc[l][m+1]);
                        m++;
                      }
                      if(m>0)m--;
	              for (k=1;k<=outnodes;k++)
		      {
                        if(anti_net[i][j][l][m][0] > 0)
                        tmp2_wts=(double)anti_net[i][j][l][m][k]*1.0/(anti_net[i][j][l][m][0]);
                        else
 		        tmp2_wts=1/outnodes;
	 //               if(l==1) classval[k]=(double)tmp2_wts*anti_wts[i][j][l][m][k];
	   //             else
	                classval[k]*=(double)tmp2_wts*anti_wts[i][j][l][m][k];
//                    cout << "[" << classval[k] <<"]";
 	              }
	              }
	             }

	             totprob=0.0;
	             for (k=1;k<=outnodes;k++)
	             {
//					cout  << "classval[" << k << "] is" << classval[k];
	               if (classval[k] > cmax)
	               {
                        c2max=classval[kmax];
	                k2max=kmax;
	                cmax=classval[k];
	                kmax=k;
                       }
                       else
                       if (classval[k]>c2max)
                       {
                        c2max=classval[k];
	                k2max=k;
                       }
	               totprob += (double)classval[k];
                     }
	          }
                  else
                  {
                   kmax=0;
                   k2max=0;
                   totprob=0.0;
                   classval[kmax]=0.0;
                   classval[k2max]=0.0;
	          }
                  if(ans1 ==3)
		  {
                   fprintf(fl4,"%d  %lf %lf %lf %lf",n, dmyclass[(int)kmax],100.0*((classval[kmax])/totprob),dmyclass[(int)k2max],100.0*((classval[k2max])/totprob));
		   if((fabs(classval[kmax]-classval[k2max]))<0.01*classval[kmax]) //classval[kmax])
		   {
		    nLoC+=classval[kmax]/totprob;
	            nLoCcnt++;
		    if(classval[kmax]>totprob*LoC)    //LoC)
		    {
		     fprintf(fl4, " <-- Either of it"); 
		    }
		    else
		    {
                      fprintf(fl4, " <-- Rejected");
		    }
		   }
		   else
		   {
		    if(classval[kmax]>totprob*LoC)    //LoC)
		    {
		      fprintf(fl4, " <-- confident");
		    }
		    else
		    {
		     fprintf(fl4, " <-- Rejected");
		    }
		   }
		   fprintf(fl4,"\n");
		 
		 }
		 if(ans1 !=3)
		 {
		  fprintf(fl4,"%d  %lf\n",n, dmyclass[(int)kmax]);
		  fprintf(fl7, "%-8d    %lf   %lf     %lf    ",n,dmyclass[(int)kmax],dmyclass[(int)k2max],tmpv);
		  if(fabs(dmyclass[kmax]-tmpv) > dmyclass[0])
		  {
		   if (classval[kmax]==0.0)
		   {
		    invcnt++;
                    fprintf(fl7, "%-7.4f %% <-Out of range %-7.4f %% \n",100.0*((classval[kmax])/totprob),100.0*((classval[k2max])/totprob));
		   }
		   else
		   {
                    if (fabs(dmyclass[k2max]-tmpv) <=dmyclass[0])
		    {
// The Next line defines the margin level required to pick up less confident examples!
//		      if ((classval[kmax]-classval[k2max]) > (0.2*totprob/outnodes))
//		      if (fabs(classval[kmax]-classval[k2max]) < 0.2*totprob/outnodes)
		     if((fabs(classval[kmax]-classval[k2max]))<0.01*classval[k2max]) //classval[kmax])
		     {
		       nLoC+=classval[kmax]/totprob;
		       nLoCcnt++;
		       if (classval[kmax]>totprob*LoC) // LoC)
		       {
		         c2cnt++;  // No more differences. NSP (OCT 2001)
		         fprintf(fl7, "%-7.4f %% <-F(1)P(2) %-7.4f %% \n",100.0*((classval[kmax])/totprob),100.0*((classval[k2max])/totprob));
		       }
		       else
		       {
		        fprintf(fl7, "%-7.4f %%  <-FMC %-7.4f %% \n",100.0*((classval[kmax])/totprob),100.0*((classval[k2max])/totprob));
		        invcnt++;
		       }
		     }
		     else
		     {
		      if (classval[kmax]>totprob*LoC) // LoC)
		      {
		       fprintf(fl7, "%-7.4f %% <-Failed %-7.4f %% \n",100.0*((classval[kmax])/totprob),100.0*((classval[k2max])/totprob));
		      }
		      else
		      {
		       fprintf(fl7, "%-7.4f %% <-FMC %-7.4f %% \n",100.0*((classval[kmax])/totprob),100.0*((classval[k2max])/totprob));
		       invcnt++;
		      }
		     }
		    }
		    else
		    {
		      if (classval[kmax]>totprob*LoC) // LoC)
		       {
		        fprintf(fl7, "%-7.4f %% <-Failed %-7.4f %% \n",100.0*((classval[kmax])/totprob),100.0*((classval[k2max])/totprob));
		       }
		       else
		       {
		        fprintf(fl7, "%-7.4f %% <-FMC %-7.4f %% \n",100.0*((classval[kmax])/totprob),100.0*((classval[k2max])/totprob));
		        invcnt++;
		       }
		    }
	          }
		}
		else
		{
		  if((fabs(classval[kmax]-classval[k2max]))<0.01*classval[kmax])
		  {
	           nLoC+=classval[kmax]/totprob;
		   nLoCcnt++;
                   if (classval[kmax]>totprob*LoC) // LoC)
		   {
		    fprintf(fl7, "%-7.4f %% <-P(1)F(2) %-7.4f %% \n",100.0*((classval[kmax])/totprob),100.0*((classval[k2max])/totprob));
		    c1cnt++;
		   }
		   else
		   {
		    invcnt++;
               	    fprintf(fl7, "%-7.4f %% <-PMC %-7.4f %% \n",100.0*((classval[kmax])/totprob),100.0*((classval[k2max])/totprob));
		   }
		  }
		  else
		 { 
                  if (classval[kmax]>totprob*LoC) // LoC)
		  {
		    fprintf(fl7, "%-7.4f %% <-Passed %-7.4f %% \n",100.0*((classval[kmax])/totprob),100.0*((classval[k2max])/totprob));
		    c1cnt++;
		  }
		  else
		  {
		    invcnt++;
		    fprintf(fl7, "%-7.4f %% <-PMC %-7.4f %% \n",100.0*((classval[kmax])/totprob),100.0*((classval[k2max])/totprob));
		  }
		 }
	       }
	       fprintf(fl5,"%d %e \n",n,(double) tmpv);
	  } // ans1 != 3 ends here ******************
	 }
	 cout << "The suggested LoC is " << nLoC/nLoCcnt << "\n";
	 fclose(fl1);
	 fclose(fl2);
	 fclose(fl4);
	 if(ans1 < 3)
	 {
	   strcpy(fltmp,fln);
	   tmp2_wts=0.0;
	   fclose(fl5);
	   fprintf(fl7,"*________________________________________________________________________\n");
	   fprintf(fl7,"*Total    Success in   Success in   Non classified   Real success in    \n");
	   cout << "*________________________________________________________________________\n";
	   cout << "*Total    Success in   Success in   Non classified   Real success in    \n";
	   if (outnodes > 3)
	    {
     		 fprintf(fl7,"* No.    Ist Choice  2nd Choice     items           two chances    \n");
     		 fprintf(fl7,"* %d       %d          %d           %d             %-7.4f %% \n",n,c1cnt,c2cnt,invcnt,(double)100.0*(c1cnt+c2cnt)/(n-invcnt));
     		 cout << "* No.    Ist Choice  2nd Choice     items           two chances    \n";
     		 printf("* %d       %d          %d           %d             %-7.4f %% \n",n,c1cnt,c2cnt,invcnt,(double)100.0*(c1cnt+c2cnt)/(n-invcnt));
     	    }
     	    else
     	    {
     		 fprintf(fl7,"* No.    Ist Choice  2nd Choice     items           First chance    \n");
     		 fprintf(fl7,"* %d       %d          %d           %d             %-7.4f %% \n",n,c1cnt,c2cnt,invcnt,(double)100.0*(c1cnt)/(n-invcnt));
     		 cout << "* No.    Ist Choice  2nd Choice     items           First chance    \n";
     		 printf("* %d       %d          %d           %d             %-7.4f %% \n",n,c1cnt,c2cnt,invcnt,(double)100.0*(c1cnt)/(n-invcnt));
     	    }
	  	 fprintf(fl7,"*________________________________________________________________________\n");
		 printf("*________________________________________________________________________\n");
		 fclose(fl7);
		 } // ******** ans1!=3 ends here *************
                 delete[] TabLookUp;
		 cout << "Done.\n";
		 stop = times(NULL);
		 cout << "The computation took " << fabs(start - stop)*10000/(CLOCKS_PER_SEC) << " Secs.\n";

 } //end main



