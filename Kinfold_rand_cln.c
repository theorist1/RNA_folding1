// TO RUN, THIS PROGRAM REQUIRES THAT ""MersenneTwister.h", WHICH IS R.J. WAGNER'S 2003 C++ IMPLEMENTATION OF MATSUMOTO AND NISHIMURA'S 1998 MERSENNE TWISTER 
// RANDOM NUMBER GENERATOR (RNG), BE INSTALLED.
// [Matsumoto M, Nishimura T (1998) Mersenne twister: A 623-dimensionally equidistributed uniform pseudo-random number generator. ACM Trans Model Comput Simulat 8:3–30.]
// A MODIFIED VERSION OF THIS CODE CAN BE FOUND HERE:  http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/C-LANG/MersenneTwister.h
// ALTERNATELY, IT CAN BE RUN USING C++'s BUILT-IN RNG.

#include <iostream>
#include <cstdlib>
#include <ctime>
#include "MersenneTwister.h"
using namespace std;

int main()
{

	int n=3000; // number of bases in each sequence
	int nseq=500; // number of random sequences generated
	int g=3; // min. hairpin size
	int tmp_id = 0;
	
	cout<<endl<<"sequence length = "<<n<<"  min. hairpin length = "<<g<<"  number sequences generated = "<<nseq<<endl<<endl;
	
	char baselt;
	int kmax; // maximum duplex length

	// MTRand mtrand1; //initialize Mersenne random number generator
	MTRand mtrand1 (n);//initialize Mersenne random number generator with a seed whose value (for convenience, should we wish to repeat the analysis) equals the sequence length

//***************ALLOCATE MEMORY FOR ARRAYS base AND index******************************

	size_t width = n+1;  
	size_t height = n+1;  
	char** base = new char* [height];
	char* block = new char [width * height];  
	for (size_t x = 0; x < height; ++x)  
	base[x] = block + width * x; 
   	
	int** index = new int* [height];
	int* block2 = new int [width * height];  
	for (size_t x = 0; x < height; ++x)  
	index[x] = block2 + width * x; 
	
	int pair[n+1]; // keeps track of the development of pairing, and gives the final pairing

	char s[n+1]; // stores the base sequence

	int maxladder[nseq+1]; // stores the max ladder distance for each of the y sequences
	int percentpair[nseq+1]; // stores the percent of bases paired for each of the y sequences
	int degeneracy[nseq+1]; // stores the "degeneracy" of each of the y sequences

//************************************************************	

	for(int y=1;y<=nseq;y++) // FOR EACH EACH SEQUENCE . . .
	{
		for (int j=0;j<=n;j++) pair[j]=0; // RE-INITIALIZE  PAIR[] FOR EACH NEW SEQUENCE
		int degeneracy=0;
		// GENERATE RANDOM SEQUENCE
		for (int d=1; d<=n; d++)
		{	
		// tmp_id = (rand() % 4); //use gcc-supplied algorithm to randomly generate integers 0 through 3
			tmp_id = mtrand1.randInt( 3 ); // use Mersenne to randomly generate integers 0 through 3
				if (tmp_id==0) {baselt='A';}
				if (tmp_id==1) {baselt='C';}
				if (tmp_id==2) {baselt='G';}
				if (tmp_id==3) {baselt='U';}
				s[d] = baselt;
				// cout<<s[d];
		}
// INSTEAD OF GENERATING THE SEQUENCE RANDOMLY, YOU CAN INSTEAD ENTER A SPECIFIC SEQUENCE AS FOLLOWS.  
// NOTE THAT SEQUENCE STARTS WITH S[1], WHICH IS WHY A DUMMY BASE, X, IS USED FOR S[0],
// AND THAT YOU NEED TO SET NSEQ TO 1, AND N TO THE NUMBER OF BASES IN YOUR SEQUENCE (NOT INCL. THE DUMMY BASE).
// char s[]={'X','C','C','A','A','G','G','A','G','G','A','A','A','A','C','C','C','C','A','G','G','G','A','G','G','G','G','A','A','A','A','C','C','C','C','A','A','A','A','C'};

//************************************************************
		int counter=0; // KEEPS TRACK OF NUMBER OF ITERATIONS
		int count=1;
		while (count!=0) // WHILE WE CAN STILL FIND NEW PAIRINGS.... (STOP WHEN WE CAN'T)
		{
		//******RE-INITIALIZE BASE[][] AND INDEX[][] FOR EACH NEW CONFIGURATION*********
		for (int j=0;j<=n;j++)  
		{  
			for (int q=0;q<=n;q++)  
			{  
			base[j][q]='X';
			}  
		}

		for (int j=0;j<=n;j++)  
		{  
			for (int q=0;q<=n;q++)  
			{  
			index[j][q]=0;
			}  
		}		
			
			count=0;
			counter++;

			/*DETERMINE LOCATION OF FIRST UNPAIRED BASE, WHICH WILL BE THE STARTING
			POINT OF LOOP #1; THE ARRAY st[] HOLDS THE STARTING POINT OF EACH LOOP*/
	
			int st[n+1];
			for (int e=0; e<=n; e++) st[e]=0; // intialize start array, assuming max of 1000 loops
			// note that this array needs to be re-initialized for every run-through of the looping, since
			// what's, say, loop 2 at one stage may become loop 4 at another.
			
			for (int a=n+1; a>0; a--)
			{
				if (pair[a]==0) 
				{
					st[1]=a; // set a equal to the starting number of the first base in the first loop
				}
			}

//*******************************************************

			int li[n+1]; // array to keep track of the current value of the index number within a loop l
			for (int f=0; f<=n; f++) li[f]=1;
			int lmax[n]; // array that holds the number of bases in each loop l; i.e., l is the number of loops 
			for (int f=0; f<=n; f++) lmax[f]=0;
			int l=1; // l is number of loops; start with one loop
			int i;
			int nl=1;
			for (i=st[l]; i<=n; i++) // start with base i = st[l] assume max of n loops
			{
			// procedure for determining h{st[l],j} (ladder distance between bases
			// st[l] and j)
				if (pair[i]==0) // if unpaired
				{					
					int check=0; // reset check to 0 (do this for each new base i)
					for (l=1; l<=nl; l++) // check to see if base i is in any of the existing loops
					// and if not, create a new loop for it 
					{
						int ladder = 0; // reset ladder to 0 (do this for each new loop l)
						for (int base_number = st[l]; base_number<=i; base_number++) // check to see if it is in same loop as st[1]
						{
							if (pair[base_number] != 0)
							{
								if (pair[base_number]>=st[l] && pair[base_number]<=base_number)
								{
									ladder=ladder-1;
								}
								else 
								{
									ladder=ladder+1;
								}
							}
						}
						if (ladder==0) // if it is in same loop
						{	
							base[l][li[l]]=s[i]; // enters, into array base, the base associated with the ith element of loop l  
							index[l][li[l]]=i;
							lmax[l]=li[l];
							li[l]=li[l]+1;
							l=nl+10; // ends the loop if ladder==0
						}
						else
						{
							check=check+1;
							if(l==nl)
							{
								nl=nl+1; // know that there is an additional loop in the structure, so increase nl (the number of loops) by 1
								base[nl][li[nl]]=s[i]; // enters, into array base, the base associated with the ith element of loop l 
								index[nl][li[nl]]=i;
								lmax[nl]=li[nl];
								li[nl]=li[nl]+1;
								st[nl]=i;
								l=nl+10; // ends the loop
							}
						}
					}
				}	
			}
			
//**********ROTATION, DONE FOR EACH LOOP*************************	
		
			for(int t=1;t<=nl;t++) // go through each loop (t is the loop number
			{		
				int lenmax=0;
				int kmax=0;
				int rmax=0;
				int rot=0;
				if ((lmax[t])%2!=0) // if the loop has an ODD number of bases, do this rotation
				{			
					for (int k=0;k<lmax[t];k++)  // rotate to a new config
					{  
						int len=0;
						// reset duplex length to zero
						for (int r=1;r<=(lmax[t]/2);r++) // go through the existing config
						{  
							if ( (abs(index[t][((k+r-1)%lmax[t])+1] - index[t][((lmax[t]+k-r)%lmax[t])+1])<=g) || (! ( ((base[t][((k+r-1)%lmax[t])+1]=='A') && (base[t][((lmax[t]+k-r)%lmax[t])+1]=='U')) || ((base[t][((k+r-1)%lmax[t])+1]=='U') && (base[t][((lmax[t]+k-r)%lmax[t])+1]=='A'))  || ((base[t][((k+r-1)%lmax[t])+1]=='C') && (base[t][((lmax[t]+k-r)%lmax[t])+1]=='G'))  || ((base[t][((k+r-1)%lmax[t])+1]=='G') && (base[t][((lmax[t]+k-r)%lmax[t])+1]=='C')) )  ) )
							// DETERMINE MAX DUPLEX
							{ // if there's no pairing, or the hairping restriction is violated:
								if (len>lenmax) // if, in addition, we are at the end of a duplex, BANK IT IF IT'S LONGER
								{
									lenmax=len; 
									kmax=k;
									rmax=r-1; // here we're actually just past the end of a duplex, so we need to subtract 1
									len=0;
								}
								len=0; // if we are at the end of a duplex, regardless of length, reset len to 0
							}
							else // otherwise, if there is pairing, 
							{
								len=len+1; // we record the increase in length
								count=count+1; // and note that we have a new pairing										
								if(r==lmax[t]/2 && len>lenmax) // and if, in addition, we're at the end of the iteration, BANK IT IF IT'S LONGER
								{
									lenmax=len; 
									kmax=k;
									rmax=r;
									len=0;
								}
								if (r==lmax[t]/2) len=0; // if we are at the end of an iteration, regardless of length, reset len to 0	
							}
						}									
					}						
				}

//*****************************************************************************
					
				if ((lmax[t])%2==0) // if the loop has an EVEN number of bases, do these two rotations				{
				{	
					for (int k=0;k<lmax[t]/2;k++)  // rotate to a new config
					{  
						int len=0;
						// reset duplex length to zero
						for (int r=1;r<=(lmax[t]/2);r++)  // go through the existing config
						{  
							if ( (abs(index[t][((k+r-1)%lmax[t])+1] - index[t][((lmax[t]+k-r)%lmax[t])+1])<=g) || (! ( ((base[t][((k+r-1)%lmax[t])+1]=='A') && (base[t][((lmax[t]+k-r)%lmax[t])+1]=='U')) || ((base[t][((k+r-1)%lmax[t])+1]=='U') && (base[t][((lmax[t]+k-r)%lmax[t])+1]=='A'))  || ((base[t][((k+r-1)%lmax[t])+1]=='C') && (base[t][((lmax[t]+k-r)%lmax[t])+1]=='G'))  || ((base[t][((k+r-1)%lmax[t])+1]=='G') && (base[t][((lmax[t]+k-r)%lmax[t])+1]=='C')) )  ) )
							// DETERMINE MAX DUPLEX
							{ // if there's no pairing, or the hairping restriction is violated:
								if (len>lenmax) // if, in addition, we are at the end of a duplex, BANK IT IF IT'S LONGER
								{
									lenmax=len; 
									kmax=k;
									rmax=r-1; // here we're actually just past the end of a duplex, so we need to subtract 1
									len=0;
								}
								len=0; // if we are at the end of a duplex, regardless of length, reset len to 0
							}
							else // otherwise, if there is pairing, 
							{
								len=len+1; // we record the increase in length
								count=count+1; // and note that we have a new pairing								
								if(r==lmax[t]/2 && len>lenmax) // and if, in addition, we're at the end of the iteration, BANK IT IF IT'S LONGER
								{
									lenmax=len; 
									kmax=k;
									rmax=r;
									len=0;
								}
								if (r==lmax[t]/2) len=0; // if we are at the end of an iteration, regardless of length, reset len to 0	
							}
						}									
					}
					for (int k=0;k<lmax[t];k++)  // rotate to a new config
					// need to go thru entire rotation, not just top, to completely fill pair[]
					{  
						int len=0;			
						// reset duplex length to zero
						for (int r=1;r<=((lmax[t]/2)-1);r++)  // go through the existing config
						{  
							if ( (abs(index[t][((k+r-1)%lmax[t])+1] - index[t][((lmax[t]+k-r-1)%lmax[t])+1])<=g) || (! ( ((base[t][((k+r-1)%lmax[t])+1]=='A') && (base[t][((lmax[t]+k-r-1)%lmax[t])+1]=='U')) || ((base[t][((k+r-1)%lmax[t])+1]=='U') && (base[t][((lmax[t]+k-r-1)%lmax[t])+1]=='A'))  || ((base[t][((k+r-1)%lmax[t])+1]=='C') && (base[t][((lmax[t]+k-r-1)%lmax[t])+1]=='G'))  || ((base[t][((k+r-1)%lmax[t])+1]=='G') && (base[t][((lmax[t]+k-r-1)%lmax[t])+1]=='C')) )  ) )
							// DETERMINE MAX DUPLEX
							{ // if there's no pairing, or the hairping restriction is violated:
								if (len>lenmax) // if, in addition, we are at the end of a duplex, BANK IT IF IT'S LONGER
								{
									lenmax=len; 
									kmax=k;
									rmax=r-1; // here we're actually just past the end of a duplex, so we need to subtract 1
									len=0;
								}
								len=0; // if we are at the end of a duplex, regardless of length, reset len to 0
							}
							else // otherwise, if there is pairing, 
							{
								len=len+1; // we record the increase in length
								count=count+1; // and note that we have a new pairing									
								if(r==lmax[t]/2-1 && len>lenmax) // and if, in addition, we're at the end of the iteration, BANK IT IF IT'S LONGER
								{
									lenmax=len; 
									kmax=k;
									rmax=r;
									len=0;
								}
								if (r==lmax[t]/2-1) len=0; // if we are at the end of an iteration, regardless of length, reset len to 0	
							}
						}									
					}						
				}
				// NOW THAT WE'VE IDENTIFIED THE LONGEST DUPLEX FOR LOOP "t," WE NEED TO FIND OUT WHERE IT IS LOCATED RELATIVE TO THE ORIGINAL
				// INDEXING, AS RECORDED IN PAIR[], AND ADJUST PAIR[] ACCORDINGLY:
				if(rot==0)
				{
					for (int r=rmax-lenmax+1;r<=rmax;r++)  // go from the beginning to the end of the longest duplex
					{
						pair[ index[t][((kmax+r-1)%lmax[t])+1] ] = index[t][((lmax[t]+kmax-r)%lmax[t])+1];
						pair[ index[t][((lmax[t]+kmax-r)%lmax[t])+1] ] = index[t][((kmax+r-1)%lmax[t])+1];
					}
				}	
				if(rot==1)
				{
					for (int r=rmax-lenmax+1;r<=rmax;r++) // go from the beginning to the end of the longest duplex
					{
						pair[ index[t][((kmax+r-1)%lmax[t])+1] ] = index[t][((lmax[t]+kmax-r-1)%lmax[t])+1];
						pair[ index[t][((lmax[t]+kmax-r-1)%lmax[t])+1] ] = index[t][((kmax+r-1)%lmax[t])+1];
					}
				}	
			}			
		} //********END OF "WHILE" LOOP, WHICH MEANS THAT WE HAVE COMPLETED ALL OF THE FOLDING FOR A SINGLE SEQUENCE**************

//*******************CALCULATE MAX LADDER DISTANCE FOR THE SEQUENCE WE JUST FINISHED FOLDING***************
		int maxladder=0;
		int templadder=0;
		int base_num = 0;
		for (int start=1; start<n; start++)
		{ 
			if (maxladder > templadder) 
			{ 
				templadder = maxladder; // reset temp 
			}
			maxladder = 0; // reset ladder		
			for (int end=start+1; end<=n; end++)
			{
				if (maxladder > templadder) 
				{ 
					templadder = maxladder; // reset temp 					
				}
				maxladder = 0; // reset maxladder
				for (base_num = start; base_num<=end; base_num++)
				{				
					if (pair[base_num] != 0)
					{
						if (pair[base_num]>=start && pair[base_num]<=base_num)
							{
								maxladder=maxladder-1;
							}
						else {maxladder=maxladder+1;}
					}
				}
			}		
		}

//************CALCULATE PERCENT BASE PAIRING FOR THE SEQUENCE WE JUST FINISHED FOLDING*************************		
		double pbp = 0.0;
		double bp=0.0;
		for (int j=0;j<=n;j++) 
		{
			if (pair[j]!=0) bp=bp+1;
		}	
		cout<<"number paired bases = "<<bp<<"  number iterations = "<<counter-1<<endl;
		pbp=(bp*100/n);
		
//***************OUTPUT RESULTS FOR THAT SEQUENCE!!!****************************************
		cout<<"maxladder = "<<templadder<<"  percentage of bases in pairs = "<<pbp<<endl;
		cout<<"{"<<templadder<<","<<pbp<<","<<counter-1<<"},"<<endl<<endl;	
						
	} //**********************************RETURN TO LINE 47 AND GENERATE A NEW SEQUENCE!!!
	return 0;
}