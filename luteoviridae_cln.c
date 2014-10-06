// THIS PROGRAM IS DESIGNED TO ANALYZE MULTIPLE VIRAL SEQUENCES, FROM THE LUETEOVIRIDAE FAMILY, USING VIENNA (specifically, RNAfold, which generates a Boltzmann-weighted sampling of 1000 sequences; doing a straight average
// over these gives a thermal average over the ensemble).  YOU FIRST NEED TO CREAT A FILE WITH THE NCBI ACCESSION NUMBERS OF THE SEQUENCES.  THIS PROGRAM WILL AUTOMATICALLY ACCESS NCBI, 
// DOWNLOAD EACH SEQUENCE IN TURN, AND ANALYZE IT.  
// TO RUN, THIS PROGRAM REQUIRES THAT THE VIENNA FOLDING PROGRAM BE INSTALLED (including the b2ct utility), AND IT ALSO REQUIRES A SHELL SCRIPT, "tmp3b2ct.sh".  HERE IS A COPY OF THAT SCRIPT, 
// WHICH SHOULD BE COPIED INTO A FILE AND SAVED UNDER THAT NAME (WITH THE COMMENT MARKS, "/*" AND "*/", REMOVED, OF COURSE):

/*
#!/bin/sh
rm tmp3.ct; rm tmp3.tmp; rm tmp3.ok
touch tmp3.ct
i="1"
while [ $i -lt 1001 ]
do
j=$i
# sed -n 'j,j p' < tmp3.b | tr -d '\012' > tmp3.tmp
head -n$i tmp3.b | tail -n1 | tr -d '\012' > tmp3.tmp
cat tmp3.seq tmp3.tmp end > tmp3.ok
cat tmp3.ok | b2ct >> tmp3.ct
i=$[$i+1]
# i='expr $i+1'
done
exit
*/

// NOTE: ALL OF THE FOLLOWING MUST BE PROPERLY SET BEFORE RUNNING PROGRAM:  !!!!!!!!!
// 1.  NAME OF FAMILY OF VIRAL SEQUENCES IS SET ON ~LINE 36 ("seqfamily")
// 2. YOU MUST CREATE FOLDERS HAVING THE NAME "seqfamily" IN: FOLDING_RESULTS, RAW_FOLDING_RESULTS, SEQUENCE_BANK, AND RNA_PICTURES.
// 3.  EACH OF THESE MUST BE RUN FROM THE DIRECTORY /Users/Aron/vienna/"seqfamiy"

// SUMMARY RESULTS ARE AUTOMATICALLY SENT TO:
// "/Users/Aron/DOCUMENTS/FOLDING_RESULTS/" << seqfamily<<"_RAW_SUMMARY_VIENNA.csv";
// [THE CSV FILE IS AUTOMATICALLY OPENED IN EXCEL.]

// Raw results are automatically sent to: (~LINE 456)
// "/Users/Aron/DOCUMENTS/RAW_FOLDING_DATA/"<<seqfamily<<"/"<<seqfamily<<"_results_VIENNA" << seqnum << ".txt";
///////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include <sstream>
#include "/Users/Aron/vienna/MersenneTwister.h"
using namespace std;

typedef int* intptr;

int main()
{

//////////////////////// create the starting sequence (with a specified base composition)

string seqfamily="Vienna_Luteoviridae"; //"Rand_Uniform_N200_C20000"
	

	int start = 1; // run from SPECIFIED START to number of sequences (the number of sequences is determined automatically from the number of accession numbers in the accession number file).
		// Watch out! In C++ an array declared as array[10] actually has elements
	// array[0], array[1],..,array[9].
	string 	name1,fname1,fnamep,fname1m, fnamesum, fname3, fname5, fnamemma, fnameRNAfold,  fseqname, accnum[1001], seqname, fname6;
	ifstream	*in1,*in2, *inp, *inRNAfold, *in5, in5a, *inc, *ine;
	ofstream	*out2, *out3, *out1m, *outsum, *out5, *outmma;
	int		n,k,l,st,num,tmp,calc,n1,n2,n3,b1,b2,lstr,i,j,m,imxLD,jmxLD,nG,nC,nA,nU,z,bin1[201],sbin,tb1,tb2;
	int	   iMALD,jMALD,number1,basenum,linewidth, seqnum, rand, temp, position, acclinenum, numbersequences;
	int		counter_a, counter_c, counter_g, counter_u, basen, loopsize, distss, distds;
	char    next, next1, word1[10],word2[10],is,name2[20],base,base1, seqname1[500];
	double	conv,AMLD,Ztot,AALD,MALD,number2,AMLDsq,sqAMLD,stdevAMLD,nfloat,temp2, temp3, dup1, dup2;
	double	stdevMALD, stdevGC, stdevAU, stdevGU, perGC, perGCsq, perAU, perAUsq, perGU, perGUsq, pG, pC, pA, pU, pGC;
	double	stdevperGU, stdevperGC, stdevperAU, NhalfALD, stdevNhalfALD, perpair, perpairsq, stdevperpair,mxLD;
	double  Adistss, Adistds, distssfloat, distdsfloat, Adistsssq, Adistdssq, sqAdistss, sqAdistds, stdevAdistss, stdevAdistds;	
	double	mfe, ensemble_fe, Abpdist, freqmfe, ensemble_diversity,RNAfold_perpair, RNAfold_MALD, RNAfold_AALD;
	double	Adup1, Adup2, Adup1sq, sqAdup1, stdevAdup1, Adup2sq, sqAdup2, stdevAdup2, numss, lenss, Alenss, Alensssq, sqAlenss, stdevAlenss; 
	//don't need to initialize branch1 ..., Abranch1 ..., etc., because these are stored using the arrays branch[] and Abranch.
	double	Abranch1sq, Abranch2sq, Abranch3sq, Abranch4sq, Abranch5sq, Abranch6sq, sqAbranch1, sqAbranch2, sqAbranch3, sqAbranch4, sqAbranch5, sqAbranch6;
	double	stdevAbranch1, stdevAbranch2, stdevAbranch3, stdevAbranch4, stdevAbranch5, stdevAbranch6;
	double	branch3p, branch4p, branch5p, branch6p, branch15p, Abranch3p, Abranch4p, Abranch5p, Abranch6p, Abranch15p;
	double  sqAbranch3p, sqAbranch4p, sqAbranch5p, sqAbranch6p, sqAbranch15p, Abranch3psq, Abranch4psq, Abranch5psq, Abranch6psq, Abranch15psq;
	double	stdevAbranch3p, stdevAbranch4p, stdevAbranch5p, stdevAbranch6p, stdevAbranch15p, maxbranch;
	double	maxLD[1001],imaxLD[1001],jmaxLD[1001],GC[1001],AU[1001],GU[1001],branch[101],Abranch[101];
	char pm;
	double branch6[1001], branch3[1001];
	pm='\xB1'; // Assign "pm" ASCII character value of 241, which gives +/- symbol.
	
	conv=1.62359; // conversion from KCal/mol to kT



	  for (k=0; k<=1000; k++)
   {
		accnum[k]="";
   }

	
	stringstream accnumfilestream; // create array containing accession numbers
	accnumfilestream << "/Users/Aron/SEQUENCE_BANK/"<< seqfamily <<"/"<< seqfamily << "_ACCESSION_NUMBERS";
	string accnumfilestring = accnumfilestream.str();
	in5 = new ifstream(accnumfilestring.c_str());
	
	if (in5->fail())
	{
		cout << "Accession number file opening failed.\n";
		exit(1);
	}
	acclinenum=1;
		while (! in5->eof())
	{
		*in5 >> accnum[acclinenum];
		acclinenum++;
	}
	numbersequences = acclinenum-2; // ASSUMING THERE IS A LINE RETURN AT THE END OF THE ACCESSION NUMBER FILE!! (OTHERWISE, SUBTRACT 1 INSTEAD OF 2)
	
	in5->close();
	delete in5;
	for(k=1;k<=numbersequences;k++)
	
	cout <<"SEQUENCE NUMBER:" <<'\t' << k<<'\t' <<"ACCESSION NUMBER:" <<'\t' <<accnum[k]<<endl;

	cout <<endl<<endl<< "NUMBER OF SEQUENCES:" <<'\t' << numbersequences<< endl<<endl;
	
	stringstream filenamesumstream;
	filenamesumstream << "/Users/Aron/DOCUMENTS/FOLDING_RESULTS/" << seqfamily<<"/"<<seqfamily<<"_SUM_"<<start<<".csv";
	fnamesum = filenamesumstream.str();
	outsum=new ofstream(fnamesum.c_str());
	
		if (outsum->fail())
	{
		cout << "Output file for summary failed to open.\n";

		exit(1);
	}


	stringstream filenamemmastream;
	filenamemmastream << "/Users/Aron/DOCUMENTS/FOLDING_RESULTS/" << seqfamily<<"/"<<seqfamily<<"_MMA_"<<start<<".csv";
	fnamemma = filenamemmastream.str();
	outmma=new ofstream(fnamemma.c_str());
	
		if (outmma->fail())
	{
		cout << "Output file for mathematica-compatible summary failed to open.\n";

		exit(1);
	}


	
*outsum<<"NCBI Accession No."<<","<<"Sequence"<<","<<"Genome Type"<<","<<"Sequence Length"<<","<<"%GC Content"<<","<<"% G Content"<<","<<"% C Content"<<","<<"% A Content"<<","<<"% U Content"<<","<<"No. Folds Generated"<<","<< "MFE/kT"<<","<< "Ensemb. FE/kT"<<","<<"Prob. MFE Struct."<<","<<"Ensemb. Divers."<<","<<"Avg. BP Dist."<<","<< "AMLD"<<","<<"Disp.(AMLD)"<<","<<"MALD"<<","<<"Disp.(MALD)"<<","<<"RNAfold MALD"<<","<<"N/2ALD"<<","<<"Disp.(N/2ALD)"<<","<<"AALD"<<","<<"RNAfold AALD"<<","<< "% of bases in pairs" <<","<< "Disp.(% in pairs)"<<","<<"RNAfold % of bases paired"<<","<< "% bases in GC pairs" <<","<< "Disp.(% in GC pairs)"<<","<<"% bases in AU pairs" <<","<<"Disp.(% in AU pairs)"<<","<<"% bases in GU pairs" <<","<<"Disp.(% in GU pairs)"<<","<<"Avg. Len. Pure Duplexes"<<","<<"Disp.(Avg. Len. Pure Dup.)"<<","<<"Avg. Len. Dup. w Bulges"<<","<<"Disp.(Avg. Len. Dup. w Bul.)"<<","<<"ss 5'-3' dist."<<","<<"Disp.(ss 5'-3' dist.)"<<","<<"ds 5'-3' dist. & deg. ext. loop"<<","<<"Disp(ds 5'-3' dist.)"<<","<<"Ends"<<","<<"Disp.(Ends)"<<","<<"Juncts."<<","<<"Disp.(Juncts.)"<<","<<"V4+"<<","<<"Disp.(V4+)"<<","<<"V5+"<<","<<"Disp.(V5+)"<<","<<"V6+"<<","<<"Disp.(V6+)"<<","<<"V15+"<<","<<"Disp.(V15+)"<<","<<"V2"<<","<<"Disp.(V2)"<<","<<"V3"<<","<<"Disp.(V3)"<<","<<"V4"<<","<<"Disp.(V4)"<<","<<"V5"<<","<<"Disp.(V5)"<<","<<"V6"<<","<<"Disp.(V6)"<<","<<"V7"<<","<<"V8"<<","<<"V9"<<","<<"V10"<<","<<"V11"<<","<<"V12"<<","<<"V13"<<","<<"V14"<<","<<"V15"<<","<<"Max. Branch"<<","<<"Centroid %Pairs"<<","<<"Centroid %GC Pairs"<<","<<"Centroid %AU Pairs"<<","<<"Centroid %GU Pairs"<<","<<"Centroid FE"<<","<< "Centroid MLD" <<endl;

for (seqnum=start; seqnum<=numbersequences; seqnum++)
 // loop through all the sequences

	{
		maxbranch = 0;
		cout<<endl<<endl<<"SEQUENCE NUMBER:" <<'\t' << seqnum<<endl<<endl;
	////////////////////////////
	// Initialize arrays and variables.
	AMLD=0;
	Ztot=1000;// WE ARE DOING A STRAIGHT AVERAGE OVER 1000 SEQUENCES
	AMLDsq=0;
	sqAMLD=0;
	pG=0, pC=0, pU=0, pA=0; 
	perGU=0, perGUsq=0, perGC=0, perGCsq=0, perAU=0, perAUsq=0,basenum=0;
	perpair=0,	perpairsq=0, Adup1=0, Adup2=0, Adistss=0, Adistds=0, Abranch3p=0, Abranch6p=0, Abranch15p=0;
	Adistsssq=0, Adistdssq=0, sqAdistss=0, sqAdistds=0, stdevAdistss=0, stdevAdistds=0;
	Adup1sq=0, sqAdup1=0, stdevAdup1=0, Adup2sq=0, sqAdup2=0, stdevAdup2=0; 
	Alenss=0, Alensssq=0, sqAlenss=0, stdevAlenss=0, Abranch1sq=0, Abranch2sq=0, Abranch3sq=0; 
	Abranch4sq=0, Abranch5sq=0, Abranch6sq=0, sqAbranch1=0, sqAbranch2=0, sqAbranch3=0, sqAbranch4=0, sqAbranch5=0, sqAbranch6=0;
	stdevAbranch1=0, stdevAbranch2=0, stdevAbranch3=0, stdevAbranch4=0, stdevAbranch5=0, stdevAbranch6=0;
	Abranch3p=0, Abranch4p=0, Abranch5p=0, Abranch6p=0, Abranch15p=0;
	sqAbranch3p=0, sqAbranch4p=0, sqAbranch5p=0, sqAbranch6p=0, sqAbranch15p=0, Abranch3psq=0, Abranch4psq=0, Abranch5psq=0, Abranch6psq=0, Abranch15psq=0;
	stdevAbranch3p=0, stdevAbranch4p=0, stdevAbranch5p=0, stdevAbranch6p=0, stdevAbranch15p=0;
	
   
   for (k=0; k<=1000; k++)
   {
		maxLD[k]=0;
		imaxLD[k]=0;
		jmaxLD[k]=0;
		GC[k]=0;
		AU[k]=0;
		GU[k]=0;
		branch3[k]=0;
		branch6[k]=0;

   }
	
	
		
	for (k=0; k<=100; k++)
	{
		Abranch[k]=0;
	}

	
	////////////////////////////////	

	stringstream curlstream;
	curlstream << "curl \"http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=genome&id=NC_"<<accnum[seqnum]<<"&retmode=text&rettype=fasta\" > tmp3.out && cat tmp3.out | sed -e 1d >tmp3.seq && cat tmp3.out | grep \"complete\" | cut -d \'|\' -f5 > tmp3a.name && awk \'{sub(/^[ \\t]+/,\"\"); print}\' tmp3a.name > tmp3.name";
	string curlstring = curlstream.str();
	
	
	if (std::system(0)) // If a command processor is available.
	{	
		
		std::system(&curlstring[0]);
	}
	else
	{
		std::cerr << "A command processor is not available.\n";
	}
	
//////////////////
// get name of sequence; it's just a temporary file, so we store it locally
	in5a.open("tmp3.name");
	
	if (in5a.fail())
	{
		cout << "first RNA name file opening failed.\n";
		exit(1);
	}
int control = 0;	
while ((! in5a.eof())&& (control==0)) // Use control variable because only want first line.
{
	getline(in5a,seqname,'\n');
	control++;
}
	in5a.close();

	in5a.open("tmp3.name");
	
	if (in5a.fail())
	{
		cout << "second RNA name file opening failed.\n";
		exit(1);
	}

z=0;
next1=' ';
while (next1 == ' ')
{
	in5a.get(next1);
	if (next1 != ' ')
	{
		seqname1[z]=next1;
		z++;
	}
}
		
while (next1 !=',')
{
in5a.get(next1);
if (next1==' ') {next1='_';}
if (next1==',') 
	break;
seqname1[z] = next1;
z++;
}	
z++;
seqname1[z]='\0';

in5a.close();	
	
cout<<endl<<endl<<"SEQUENCE NUMBER:" <<'\t' << seqnum<<'\t' <<"SEQUENCE NAME:" <<'\t' <<seqname<<endl<<endl;
		
/////////////////////	

	
	in2 = new ifstream("tmp3.seq");
	
	if (in2->fail())
	{
		cout << "Input fasta sequence file opening failed.\n";
		exit(1);
	}


	stringstream filename5stream;
	filename5stream << "/Users/Aron/SEQUENCE_BANK/"<<seqfamily<<"/"<<seqfamily<<"_"<<seqnum<<".seq";
	fname5 = filename5stream.str();
	out5 = new ofstream(fname5.c_str());
	if (out5->fail())
	{
		cout << "Output file for sequence banks failed to open.\n";

		exit(1);
	}
	
	while (! in2->eof())
	{
		in2->get(next);
		if (in2->eof()) {break;}
		*out5 << next;
	}			
	out5->close();
	delete out5;
	in2->close();
	delete in2;
			
////////////////////////////////////////////////////////////////////
        
if (std::system(0)) // If a command processor is available.
	{	// CALL Vienna in traceback mode
		std::system("cat tmp3.seq |  tr -d \'\012\' > tmp3a.seq; cat /Users/Aron/vienna/return >> tmp3a.seq; mv tmp3a.seq tmp3.seq");
	}
	else
	{
		std::cerr << "A command processor is not available.\n";
	}


////////////////////////////////////////////////////////////////CALL VIENNA IN TRACEBACK MODE
	
	if (std::system(0)) // If a command processor is available.
	{	// CALL Vienna in traceback mode
		std::system("/Users/Aron/vienna/Progs/RNAsubopt -d2 -p 1000 < tmp3.seq > tmp3.b");// RUN_TYPE=html
	}
	else
	{
		std::cerr << "A command processor is not available.\n";
	}


//////////////////////////////////////////////////////////////// CALL VIENNA IN PARTITION FUNCTION MODE
	
	if (std::system(0)) // If a command processor is available.
	{	// CALL Vienna in partition function mode
		std::system("/Users/Aron/vienna/Progs/RNAfold -d2 -p -S 0.9 < tmp3.seq > partition.txt");// RUN_TYPE=html
	}
	else
	{
		std::cerr << "A command processor is not available.\n";
	}

//////////////////////////////////////////////////////////////// CONVERT THE tmp3.b FILE TO A tmp3.ct file
	
	
	if (std::system(0)) // If a command processor is available.
	{	
		std::system("/Users/Aron/vienna/tmp3b2ct.sh");// RUN_TYPE=html
	}
	else
	{
		std::cerr << "A command processor is not available.\n";
	}


///////////////////


	fname1="tmp3.ct";
	in1 = new ifstream(fname1.c_str());

	if (in1->fail())
	{
		cout << "Input .ct file opening failed.\n";
		exit(1);
	}
	
	
   // Our third use of the sequence file is to create an array to put sequence in and fill it.
   // Determine number of bases.
	*in1 >> n;

   in2 = new ifstream("tmp3.seq");
	if (in2->fail())
	{
		cout << "Input .seq file opening failed.\n";
		exit(1);
	}
	
	nG=0;
	nC=0;
	nA=0;
	nU=0;
	

  	intptr sen;
   sen=new int[n+1];
   k = 1;
	while (! in2->eof())
	{
		in2->get(next);
		if (in2->eof()) {break;}
		if (next=='A'||next=='a')
		{
			sen[k]=3;
			nA++;
			k++;
		}
		if (next=='C'||next=='c')
		{
			sen[k]=4;
			nC++;
			k++;
		}
		if (next=='G'||next=='g')
		{
			sen[k]=1;
			nG++;
			k++;
		}
		if (next=='U'||next=='T'||next=='u'||next=='t')
		{
			sen[k]=0;
			nU++;
			k++;
		}
	}
	in2->close();
	delete in2;
	
	
	
   // Now create all the matrices and initialize them.ALDsq[i][j]
	size_t width = n+1;
	size_t height = n+1;
	
	double** LD = new double* [height];
	double* blockLD = new double [width * height];
	for (size_t x = 0; x < height; ++x)
	LD[x] = blockLD + width * x;
	
	double** ALD = new double* [height];
	double* blockALD = new double [width * height];
	for (size_t x = 0; x < height; ++x)
	ALD[x] = blockALD + width * x;
	
	double** ALDsq = new double* [height];
	double* blockALDsq = new double [width * height];
	for (size_t x = 0; x < height; ++x)
	ALDsq[x] = blockALDsq + width * x;
		
	for (k=0; k<=n; k++)
	{
		for (l=0; l<=n; l++)
		{
			LD[k][l]=0;
			ALD[k][l]=0;
			ALDsq[k][l]=0;
		}
	}
	
	// Create the matrices to put the structure in and initialize them.
	intptr structure1,structure2,ct1,ct2;
	structure1=new int[n+1];
	structure2=new int[n+1];
	ct1=new int[n+1];
	ct2=new int[n+1];
	for (k=0; k<=n; k++)
	{
		structure1[k]=0;
		structure2[k]=0;
		ct1[k]=0;
		ct2[k]=0;
	}


	st=1;
	num=1;

	// Read the first line.
	*in1 >> word1 >> is >> number2 >> word2;
	in1->get(next);
	while (! (next=='\n'))
	{
		in1->get(next);
	}

	// Read the rest.
	

	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// THIS IS THE START OF BIG "WHILE" LOOP, IN WHICH WE GO THROUGH THE CT FILE, ONE STRUCTURE AT A TIME:
	
	while (st)
	{
		lstr=0;
		for (k=1; k<=n; k++) //this reads in ONE structure, starting from one line beyond the end of the last structure read in
		{
			*in1 >> b1 >> base >> n1 >> n2 >> b2 >> n3;
			
			ct1[k]=b1;//these two lines read in the full .ct pairing info., including unpaired bases
			ct2[k]=b2;//these two arrays will be used to calculate the degree of branching
			
			if (b2!=0)
         	{
            	if (b1<b2)
				{
					lstr++;
					structure1[lstr]=b1;
					structure2[lstr]=b2;
				}
			}
			
		}
		// Now do the calculations.
		// First the ladder distances.
		for (k=0; k<=n; k++)
		{
			for (l=0; l<=n; l++)
			{
				LD[k][l]=0;
			}
		}
		for (k=1; k<=lstr; k++)
		{
			i=structure1[k];
			j=structure2[k];
			for (l=1; l<=i-1; l++)
			{
				for (m=i+1; m<=j-1; m++)
				{
					LD[l][m]++;
				}
			}
			for (l=i+1; l<=j-1; l++)
			{
				for (m=j+1; m<=n; m++)
				{
					LD[l][m]++;
				}
			}
		}
		// Determine MLD;
		mxLD=0;
		imxLD=0;
		jmxLD=0;
		for (i=1; i<=n-1; i++)
		{
			for (j=i+1; j<=n; j++)
			{
				if (LD[i][j]>mxLD)
				{
					mxLD=LD[i][j];
					imxLD=i;
					jmxLD=j;
				}
			}
		}
		maxLD[num]=mxLD;
		imaxLD[num]=imxLD;
		jmaxLD[num]=jmxLD;
		// boltz=exp(-(energy[num]-energy[1]));
		// Ztot=Ztot+boltz;
        
		// Collect data for MALD.
		for (i=1; i<=n-1; i++)
		{
			for (j=i+1; j<=n; j++)
			{
				ALD[i][j]=ALD[i][j]+LD[i][j];
				ALDsq[i][j]=ALDsq[i][j]+LD[i][j]*LD[i][j];
			}
		}
		// Collect data for AMLD.
		AMLD=AMLD+mxLD;
		AMLDsq=AMLDsq+mxLD*mxLD;
		
		//////////////////////////////////////////////////////////////////////////   CALCULATE DEGREE OF BRANCHING!!!!!
		for (k=0; k<=100; k++)
		{
		branch[k]=0;
		}
		basen=1;
		
		while (basen<=n)
		{
			loopsize=1; // this is the minimum loopsize (1 = hairpin)
			if (ct2[basen]==0 || ct1[basen]>ct2[basen]) {basen++;} // if the base is unpaired, or paired to a base with a lower no., proceed to the next base
			else // if the base is paired to a base with a higher no.
			{
				i=ct1[basen]; // no. of the base
				j=ct2[basen]; // no. of its pair
				basen++; // go to the next base
				while (basen!=j) // keep doing the following until you come around to the pairing base
				{
					if (ct2[basen]==0 || ct1[basen]>ct2[basen]) {basen++;} // if the base is unpaired, or paired to a base with a lower no., proceed to the next base
					else // if it's paired to a base with a higher no., we have a duplex coming off of the loop
					{
						loopsize++; // count the no. of duplexes coming off of the loop
						basen=ct2[basen]; // jump across to the pairing base
						if (basen!=j) {basen++;} // if we're not at the end of the loop (if we are, we keep basen=j, thus causing us to exit from the while loop
					}	
				}		
				branch[loopsize]++; // increase by one the running count of the number of branches of degree "loopsize"
				basen=i+1; // go to the next possible loop
			}	
		}

		// CALCULATE DEGREE OF BRANCHING FOR EXTERIOR LOOP, AS WELL AS 5' - 3' DISTANCE!!
			numss=0, lenss=0;
			distss=0;
			distds=0;
			loopsize=0;
			branch3p=0, branch4p=0, branch5p=0, branch6p=0, branch15p=0;
			
			basen=1;
			
		while (basen<n) // don't want to do this operation for the last base, since this will incorrectly increase distss by 1
		{
			if (ct2[basen]==0 || ct1[basen]>ct2[basen]) // except for the last base, should never have a case in which  ct1[basen]>ct2[basen]			
			{
				basen++;
				distss++;
			} 
			else 
			{
				basen=ct2[basen];
				distds++;
				loopsize++;			
			}
		}
	
			branch[loopsize]++;
		
			for (k=1; k<=100; k++)
		{
		
				Abranch[k]=Abranch[k]+branch[k];

		}	
	
		Abranch1sq=Abranch1sq+branch[1]*branch[1];	
		Abranch2sq=Abranch2sq+branch[2]*branch[2];	
		Abranch3sq=Abranch3sq+branch[3]*branch[3];	
		Abranch4sq=Abranch4sq+branch[4]*branch[4];	
		Abranch5sq=Abranch5sq+branch[5]*branch[5];	
		Abranch6sq=Abranch6sq+branch[6]*branch[6];
		
			
		branch3[num]=branch[3];				
		branch6[num]=branch[6];

						
						/////////CONSOLIDATE BRANCHING RESULTS
								
		for (k=15; k<=100; k++)
		{
		branch15p=branch15p+branch[k];
		}
		
		branch6p=branch15p;
		for (k=6; k<15; k++)
		{
		branch6p=branch6p+branch[k];
		}

		branch5p=branch6p+branch[5];
		branch4p=branch5p+branch[4];
		branch3p=branch4p+branch[3];
		
		Abranch3p=Abranch3p+branch3p;
		Abranch3psq=Abranch3psq+branch3p*branch3p;
		Abranch4p=Abranch4p+branch4p;
		Abranch4psq=Abranch4psq+branch4p*branch4p;
		Abranch5p=Abranch5p+branch5p;
		Abranch5psq=Abranch5psq+branch5p*branch5p;
		Abranch6p=Abranch6p+branch6p;
		Abranch6psq=Abranch6psq+branch6p*branch6p;
		Abranch15p=Abranch15p+branch15p;
		Abranch15psq=Abranch15psq+branch15p*branch15p;
			
		for (k=1; k<=100; k++) // determine loop with largest braching among ALL the structures in the ensemble 
		// (which is why we initialize OUTSIDE of the large while loop)
		{	
			if (branch[k]!=0 && k>maxbranch)
			{
				maxbranch=k;
			}
		}	
		
		// determine number of ss sections in the exterior loop
		numss=distds+1; // the number of ss sections is equal to the number of ds "dividers" plus one, unless:
		if (ct1[1]!=0 || ct2[n]!=0) {numss=numss-1;} // if at least one of the ends is paired
		if (ct1[1]!=0 && ct2[n]!=0) {numss=1;} // if both ends are paired, there are no single-stranded sections;
		// hence we set numss=1 just because we don't want to divide by 0.
		// determine the avg. length of ss sections in the exterior loop
		lenss=distss/numss;		

		
	
		distssfloat=static_cast<double>(distss);
		distdsfloat=static_cast<double>(distds);
		Adistss=Adistss+distssfloat;
		Adistds=Adistds+distdsfloat; // the double-stranded 5'-3' distance is the same as the degree of branching of the exterior loop
		Adistsssq=Adistsssq+distss*distss;
		Adistdssq=Adistdssq+distds*distds;
					
		Alenss=Alenss+lenss;
		Alensssq=Alensssq+lenss*lenss;


		
	//////////////////////////////////////////////////////////////////////////////  CALCULATE AVERAGE DUPLEX LENGTH
	


	// Calculate duplex size, not allowing bulges of size 1.

temp2=0;
temp3=0;
dup1=0;
   for (i=1; i<=200; i++)
   {
    	bin1[i]=0;
   }

	b1=structure1[1];
	b2=structure2[1];
   sbin=1;
   tb1=b1;
   tb2=b2;
	for (i=2; i<=lstr; i++)
   {
	b1=structure1[i];
	b2=structure2[i];

      if ((b1-tb1==1)&&(tb2-b2==1))
      {
         sbin++;
      }
      else
      {
        	bin1[sbin]++;
       		   sbin=1;
		}
      tb1=b1;
      tb2=b2;
   }
   bin1[sbin]++;
   
   for (i=1; i<=200; i++)
   {
    	temp2=temp2+bin1[i];
		temp3=temp3+bin1[i]*i;
   }
	dup1=temp3/temp2;
	Adup1=Adup1+dup1;


	// Calculate duplex size, allowing bulges of size 1.
	temp2=0;
	temp3=0;
	dup2=0;
	
   for (i=1; i<=200; i++)
   {
    	bin1[i]=0;
   }

	b1=structure1[1];
	b2=structure2[1];
   sbin=1;
   tb1=b1;
   tb2=b2;
	for (i=2; i<=lstr; i++)
   {
	b1=structure1[i];
	b2=structure2[i];

     if (((b1-tb1==1)&&((tb2-b2==1)||(tb2-b2==2)))||((b1-tb1==2)&&(tb2-b2==1)))
      {
         sbin++;
      }
      else
      {
        	bin1[sbin]++;
       		   sbin=1;
		}
      tb1=b1;
      tb2=b2;
   }
   bin1[sbin]++;
   
      for (i=1; i<=200; i++)
   {
    	temp2=temp2+bin1[i];
		temp3=temp3+bin1[i]*i;
   }
	dup2=temp3/temp2;
	Adup2=Adup2+dup2;
	
	Adup1sq=Adup1sq+dup1*dup1;
	Adup2sq=Adup2sq+dup2*dup2;
	
	


		///////////////////////////////////////////////////////////////////////// IF USE RAG, IT GOES HERE!
		
		////////////////////////////////////////////////////////////////////////////
		
		for (k=1; k<=lstr; k++)
		{
			i=structure1[k];
			j=structure2[k];
			b1=sen[i];
			b2=sen[j];
			if (((b1==1)&&(b2==4))||((b1==4)&&(b2==1)))
			{
				GC[num]++;
			}
			if (((b1==0)&&(b2==3))||((b1==3)&&(b2==0)))
			{
				AU[num]++;
			}
			if (((b1==0)&&(b2==1))||((b1==1)&&(b2==0)))
			{
			GU[num]++;
			}
		}
		nfloat=static_cast<double>(n);
		perGC=perGC+GC[num]*200.0/nfloat;//2 bases per pair
		perGCsq=perGCsq+GC[num]*GC[num]*40000.0/(nfloat*nfloat);
		perAU=perAU+AU[num]*200.0/nfloat;
		perAUsq=perAUsq+AU[num]*AU[num]*40000.0/(nfloat*nfloat);
		perGU=perGU+GU[num]*200.0/nfloat;
		perGUsq=perGUsq+GU[num]*GU[num]*40000.0/(nfloat*nfloat);
		perpair=perpair+(GU[num]+GC[num]+AU[num])*200.0/nfloat;//2 bases per pair
		perpairsq=perpairsq+(GU[num]+GC[num]+AU[num])*40000.0*(GU[num]+GC[num]+AU[num])/(nfloat*nfloat);
	
		//cout << num << " " << maxLD[num] << "	V1:	"<<branch[1]<<"	V3:	"<<branch[3]<<"	V4:	"<<branch[4]<<"	SSDIST:	"<<distss<<"	DSDIST:		"<<distds<<"	LOOPSZ:		"<<loopsize<<endl;
		

				
		//cout << num << " " << energy[num] << " " << maxLD[num] << " "
		//<< imaxLD[num] << " " << jmaxLD[num] << " " << GC[num] << " "
		//<< AU[num] << " " << GU[num] <<" "<<dup2<<" "<<branch[1]<<" "<<branch[2]<<" "<<branch[3]<<" "<<branch[4]<<" "<<distss<<" "<<distds<<" "<<loopsize<<endl;
		
		cout << num << " "<< maxLD[num] << " "
		<< imaxLD[num] << " " << jmaxLD[num] << " " << GC[num] << " "
		<< AU[num] << " " << GU[num] << '\n';
      // You have to read an extra character to get to the end of the file if
      // that's where you are. Otherwise it doesn't matter.
		in1->get(next);
		in1->get(next);
		if (in1->eof())
		{
			st=0;
		}
		else
		{
			num++;
			*in1 >> number1 >> word1 >> is >> number2 >> word2;
			in1->get(next);
			while (! (next=='\n'))
			{
				in1->get(next);
			}
		}
   }   // END OF LOOP BEGINNING "while (st)" @ LINE 481!!!!!!!!!!!
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   in1->close();
   delete in1;

   // Calculate MALD, AALD and AMLD.

   AALD=0;
   AMLD=AMLD/Ztot;
   AMLDsq=AMLDsq/Ztot;
   sqAMLD=AMLD*AMLD;
   stdevAMLD=sqrt(fabs(AMLDsq-sqAMLD));
	perGU=perGU/Ztot;
	perGUsq=perGUsq/Ztot;
	perAU=perAU/Ztot;
	perAUsq=perAUsq/Ztot;
	perGC=perGC/Ztot;
	perGCsq=perGCsq/Ztot;
   perpair=perpair/Ztot;
   perpairsq=perpairsq/Ztot;
   stdevperGU=sqrt(fabs(perGUsq-perGU*perGU));
   stdevperGC=sqrt(fabs(perGCsq-perGC*perGC));
   stdevperAU=sqrt(fabs(perAUsq-perAU*perAU));
   stdevperpair=sqrt(fabs(perpairsq-perpair*perpair));

	Adistss=Adistss/Ztot;
	Adistds=Adistds/Ztot;
	Adistsssq=Adistsssq/Ztot;
	Adistdssq=Adistdssq/Ztot;
	sqAdistss=Adistss*Adistss;
	sqAdistds=Adistds*Adistds;
	stdevAdistss=sqrt(fabs(Adistsssq-sqAdistss));
	stdevAdistds=sqrt(fabs(Adistdssq-sqAdistds));
	Alenss=Alenss/Ztot;
	Alensssq=Alensssq/Ztot;
	sqAlenss=Alenss*Alenss;
	stdevAlenss=sqrt(fabs(Alensssq-sqAlenss));
		
	Adup1=Adup1/Ztot;
	Adup1sq=Adup1sq/Ztot;
	sqAdup1=Adup1*Adup1;
	stdevAdup1=sqrt(fabs(Adup1sq-sqAdup1));
	
	Adup2=Adup2/Ztot;
	Adup2sq=Adup2sq/Ztot;
	sqAdup2=Adup2*Adup2;
	stdevAdup2=sqrt(fabs(Adup2sq-sqAdup2));
	
	
		for (k=1; k<=15; k++) // determine no. of loop sizes 1 through 15
		{
		
				Abranch[k]=Abranch[k]/Ztot;

		}	
		
		// determine dispersions for number of loops, sizes 3 through 6

		Abranch1sq=Abranch1sq/Ztot;
		Abranch2sq=Abranch2sq/Ztot;
		Abranch3sq=Abranch3sq/Ztot;
		Abranch4sq=Abranch4sq/Ztot;
		Abranch5sq=Abranch5sq/Ztot;
		Abranch6sq=Abranch6sq/Ztot;
		
		sqAbranch1=Abranch[1]*Abranch[1];
		sqAbranch2=Abranch[2]*Abranch[2];
		sqAbranch3=Abranch[3]*Abranch[3];
		sqAbranch4=Abranch[4]*Abranch[4];
		sqAbranch5=Abranch[5]*Abranch[5];
		sqAbranch6=Abranch[6]*Abranch[6];
	
		stdevAbranch1=sqrt(fabs(Abranch1sq-sqAbranch1));
		stdevAbranch2=sqrt(fabs(Abranch2sq-sqAbranch2));
		stdevAbranch3=sqrt(fabs(Abranch3sq-sqAbranch3));
		stdevAbranch4=sqrt(fabs(Abranch4sq-sqAbranch4));
		stdevAbranch5=sqrt(fabs(Abranch5sq-sqAbranch5));
		stdevAbranch6=sqrt(fabs(Abranch6sq-sqAbranch6));
			
		// determine no. of loops, 3+ through 6+, and 15+, and the dispersions of these values		
								
		Abranch3p=Abranch3p/Ztot;
		Abranch3psq=Abranch3psq/Ztot;
		sqAbranch3p=Abranch3p*Abranch3p;
		stdevAbranch3p=sqrt(fabs(Abranch3psq-sqAbranch3p));	

		Abranch4p=Abranch4p/Ztot;
		Abranch4psq=Abranch4psq/Ztot;
		sqAbranch4p=Abranch4p*Abranch4p;
		stdevAbranch4p=sqrt(fabs(Abranch4psq-sqAbranch4p));	
		
		Abranch5p=Abranch5p/Ztot;
		Abranch5psq=Abranch5psq/Ztot;
		sqAbranch5p=Abranch5p*Abranch5p;
		stdevAbranch5p=sqrt(fabs(Abranch5psq-sqAbranch5p));	
		
		Abranch6p=Abranch6p/Ztot;
		Abranch6psq=Abranch6psq/Ztot;
		sqAbranch6p=Abranch6p*Abranch6p;
		stdevAbranch6p=sqrt(fabs(Abranch6psq-sqAbranch6p));	
		
		Abranch15p=Abranch15p/Ztot;
		Abranch15psq=Abranch15psq/Ztot;
		sqAbranch15p=Abranch15p*Abranch15p;
		stdevAbranch15p=sqrt(fabs(Abranch15psq-sqAbranch15p));	

   
   /////////////////// Collect screen output of RNAfold -p
	mfe=0, ensemble_fe=0, Abpdist=0, freqmfe=0, ensemble_diversity=0;
	
	
	
	
if (std::system(0)) //If a command processor is available.
{	
		std::system("grep - partition.txt | awk '{print $NF}' | head -n1 | tail -n1 |  sed 's/[()]//g' > partition.out && grep - partition.txt | awk '{print $NF}' | head -n2 | tail -n1 |  sed 's/\\[//g' | sed 's/\\]//g' >> partition.out && grep d partition.txt | head -n1 |  awk '{ print $NF }' | sed '$s/.$//' | sed 's/d=//g' >> partition.out && grep frequency partition.txt | awk '{print $7}'  | sed 's/;//' >> partition.out && grep frequency partition.txt | awk '{print $10}'  >> partition.out");
	}
	else
	{
		std::cerr << "A command processor is not available.\n";
	}


	fnamep="partition.out";
	inp = new ifstream(fnamep.c_str());

	if (inp->fail())
	{
		cout << "Input partition.out data file opening failed.\n";
		exit(1);
	}
	mfe=0, ensemble_fe=0, Abpdist=0, freqmfe=0, ensemble_diversity=0;
	

	*inp >> mfe>> ensemble_fe>> Abpdist>> freqmfe>> ensemble_diversity;  	
	mfe=mfe*conv, ensemble_fe=ensemble_fe*conv; //convert mfe and ensemble_fe from KCal/mol to kT

////////////////////////////// Convert dot.ps to pairprob3.txt, so that petermald.exe can read it
	
	if (std::system(0)) // If a command processor is available.
	{	// CALL Vienna in traceback mode
		std::system("grep ubox dot.ps | sed \'1,2d\' > pairprob3.txt");
	}
	else
	{
		std::cerr << "A command processor is not available.\n";
	}


//////////////////////////// Call petermald.exe and capture output
if (std::system(0)) // If a command processor is available.
	{	// CALL Vienna in traceback mode
		std::system("/Users/Aron/vienna/petermald.exe tmp3");
	}
	else
	{
		std::cerr << "A command processor is not available.\n";
	}

	fnameRNAfold="tmp3.dat";
	inRNAfold = new ifstream(fnameRNAfold.c_str());

	if (inRNAfold->fail())
	{
		cout << "Input tmp3.dat function data file opening failed.\n";
		exit(1);
	}
	

	*inRNAfold >> RNAfold_perpair >> RNAfold_MALD >> RNAfold_AALD;  	

inRNAfold->close();
delete inRNAfold;
 
//////////////////////////
//////////////////////////  DETERMINE CENTROID
//////////////////////////////////////////////
double centroid_energy=0; double nfloat2=0; double MLDcent=0;
int GCcent = 0; int AUcent=0; int GUcent = 0; double perGCcent=0; double perAUcent=0;double perGUcent=0; double perpaircent=0; 
	

if (std::system(0)) // If a command processor is available.
	{	// CALL Vienna in traceback mode
		std::system("sed -n '1,1 p' partition.txt > test && cat partition.txt | sed -n '4,4 p' | sed s/\\{/\\(/g | sed s/\\}/\\)/g >> test && cat test | /home/aron/b2ct.exe > centroid.ct");
		std::system("head -n1 centroid.ct | awk '{print $4}' > centroid_energy");
		std::system("sed '1d' centroid.ct > centroid2.ct");

	}
	else
	{
		std::cerr << "A command processor is not available.\n";
	}

	ine = new ifstream("centroid_energy");
	if (ine->fail())
	{
		cout << "Input centroid2.ct file opening failed.\n";
		exit(1);
	}
	
	*ine >> centroid_energy;
	ine->close();
	delete ine;
	
	cout << "CENTROID ENERGY:	" << centroid_energy << endl;





	inc = new ifstream("centroid2.ct");
	if (inc->fail())
	{
		cout << "Input centroid2.ct file opening failed.\n";
		exit(1);
	}
	

	
	lstr=0;
		for (k=1; k<=n; k++) // this reads in ONE structure, starting from one line beyond the end of the last structure read in
		{
			*inc >> b1 >> base >> n1 >> n2 >> b2 >> n3;
			
			ct1[k]=b1; // these two lines read in the full .ct pairing info., including unpaired bases
			ct2[k]=b2; // these two arrays will be used to calculate the degree of branching
			
			if (b2!=0)
         	{
            	if (b1<b2)
				{
					lstr++;
					structure1[lstr]=b1;
					structure2[lstr]=b2;
				}
			}
			
		}
		// Now do the calculations.
		// First the ladder distances.
		for (k=0; k<=n; k++)
		{
			for (l=0; l<=n; l++)
			{
				LD[k][l]=0;
			}
		}
		for (k=1; k<=lstr; k++)
		{
			i=structure1[k];
			j=structure2[k];
			for (l=1; l<=i-1; l++)
			{
				for (m=i+1; m<=j-1; m++)
				{
					LD[l][m]++;
				}
			}
			for (l=i+1; l<=j-1; l++)
			{
				for (m=j+1; m<=n; m++)
				{
					LD[l][m]++;
				}
			}
		}
		// Determine MLD;
		mxLD=0;
		imxLD=0;
		jmxLD=0;
		for (i=1; i<=n-1; i++)
		{
			for (j=i+1; j<=n; j++)
			{
				if (LD[i][j]>mxLD)
				{
					mxLD=LD[i][j];
					imxLD=i;
					jmxLD=j;
				}
			}
		}
		
	MLDcent = mxLD;		
		
		for (k=1; k<=lstr; k++)
		{
			i=structure1[k];
			j=structure2[k];
			b1=sen[i];
			b2=sen[j];
			if (((b1==1)&&(b2==4))||((b1==4)&&(b2==1)))
			{
				GCcent++;
			}
			if (((b1==0)&&(b2==3))||((b1==3)&&(b2==0)))
			{
				AUcent++;
			}
			if (((b1==0)&&(b2==1))||((b1==1)&&(b2==0)))
			{
			GUcent++;
			}
		}

		nfloat2=static_cast<double>(n);
		perGCcent=GCcent*200.0/nfloat2;//2 bases per pair
		perAUcent=AUcent*200.0/nfloat2;
		perGUcent=GUcent*200.0/nfloat2;
		perpaircent=(GUcent+GCcent+AUcent)*200.0/nfloat2;//2 bases per pair
	
inc->close();
delete inc;

///////////////////////////
//////////////////////////

   
	for (i=1; i<=n-1; i++)
	{
		for (j=i+1; j<=n; j++)
		{
			ALD[i][j]=ALD[i][j]/Ztot;
		}
	}
	MALD=0;
	iMALD=0;
	jMALD=0;
	for (i=1; i<=n-1; i++)
	{
		for (j=i+1; j<=n; j++)
		{
			AALD=AALD+ALD[i][j];
			if (ALD[i][j]>MALD)
			{
				MALD=ALD[i][j];
				iMALD=i;
				jMALD=j;
			}
		}
	}
	
	int midpoint = static_cast<int>(floor(nfloat/2+1));
	AALD=AALD/((n*(n-1))/2);
	NhalfALD=ALD[1][midpoint];
	
	stdevNhalfALD=sqrt(fabs(ALDsq[1][midpoint]/Ztot-(NhalfALD*NhalfALD)));
	
	stdevMALD=sqrt(fabs(ALDsq[iMALD][jMALD]/Ztot-(MALD*MALD)));
	
	pGC = (static_cast<double>(nG)+static_cast<double>(nC))*100.0/nfloat;
	pG = static_cast<double>(nG)*100.0/nfloat;
	pC = static_cast<double>(nC)*100.0/nfloat;
	pA = static_cast<double>(nA)*100.0/nfloat;
	pU = static_cast<double>(nU)*100.0/nfloat;
	
	// Write to file.
	// Open the .mlds file.
	stringstream filename3stream;
	filename3stream << "/Users/Aron/DOCUMENTS/RAW_FOLDING_DATA/"<<seqfamily<<"/"<<seqfamily<<"_"<<seqnum<<"_results_VIENNA.txt";
	fname3 = filename3stream.str();
	out3 = new ofstream(fname3.c_str());
	if (out3->fail())
	{
		cout << "Output file for all data failed to open.\n";

		exit(1);
	}
	*out3 <<" RNA"<<'\t'<<seqfamily<<"_"<<seqnum<<endl;
	*out3 << "Sequence Length"<<'\t'<< n << '\n';
	*out3 << "%GC Content"<<'\t'<< pGC << '\n';
	*out3 << "% G Content" <<'\t'<< pG << '\n';
	*out3 << "% C Content" <<'\t'<< pC << '\n';
	*out3 << "% A Content" <<'\t'<< pA << '\n';
	*out3 << "% U Content" <<'\t'<< pU << '\n';
	*out3 << "No. Folds Generated"<<'\t'<< num<< '\n';
	*out3 << "MFE/kT" <<'\t'<<mfe << '\n';
	*out3 << "AMLD" <<'\t'<< AMLD <<" "<<pm<<" "<< stdevAMLD << '\n';
	*out3 << "MALD" <<'\t'<< MALD<<" "<<pm<<" "<< stdevMALD << '\n';
	*out3 << "Index of first base:  " << iMALD << '\n';
	*out3 << "Index of second base: " << jMALD << '\n';
	*out3 << "N/2ALD" <<'\t'<< NhalfALD<<" "<<pm<<" "<< stdevNhalfALD << '\n';
	*out3 << "AALD" <<'\t'<< AALD << '\n';
	*out3 << "% of bases in pairs" <<'\t'<<  perpair<<" "<<pm<<" "<< stdevperpair << '\n';
	*out3 << "% bases in GC pairs" <<'\t'<<  perGC<<" "<<pm<<" "<< stdevperGC << '\n';
	*out3 << "% bases in AU pairs" <<'\t'<<  perAU<<" "<<pm<<" "<< stdevperAU << '\n';
	*out3 << "% bases in GU pairs" <<'\t'<<  perGU<<" "<<pm<<" "<< stdevperGU << '\n';
	*out3 << "avg. len. pure duplexes" << '\t' << Adup1 <<endl;
	*out3 << "avg. len. duplexes incl. bulges" << '\t' << Adup2 <<endl;
	*out3 << "Index of first base:  " << iMALD << '\n';
	*out3 << "Index of second base: " << jMALD << '\n';
	*out3 << '\n';
	*out3 << "Number, Energy (kT), MLD, i, j, nGC, nAU, nGU: \n";

	*outsum << accnum[seqnum] <<","<< seqname <<","<<  n <<","<<  pGC <<","<<  pG <<","<<  pC <<","<<  pA <<","<<   pU <<","<<  num <<","<<mfe<<","<< ensemble_fe<<","<<freqmfe<<","<<ensemble_diversity<<","<<Abpdist<<","<<AMLD<<","<< stdevAMLD<<","<< MALD<<","<< stdevMALD<<","<<RNAfold_MALD<<","<<NhalfALD<<","<<stdevNhalfALD<<","<< AALD<<","<<RNAfold_AALD<<","<<perpair<<","<< stdevperpair<<","<<RNAfold_perpair<<","<<perGC<<","<<stdevperGC<<","<<perAU<<","<<stdevperAU<<","<<perGU<<","<<stdevperGU<<","<<Adup1<<","<<stdevAdup1<<","<<Adup2<<","<<stdevAdup2<<","<<Adistss<<","<<stdevAdistss<<","<<Adistds<<","<<stdevAdistds<<","<<Abranch[1]<<","<<stdevAbranch1<<","<<Abranch3p<<","<<stdevAbranch3p<<","<<Abranch4p<<","<<stdevAbranch4p<<","<<Abranch5p<<","<<stdevAbranch5p<<","<<Abranch6p<<","<<stdevAbranch6p<<","<<Abranch15p<<","<<stdevAbranch15p<<","<<Abranch[2]<<","<<stdevAbranch2<<","<<Abranch[3]<<","<<stdevAbranch3<<","<<Abranch[4]<<","<<stdevAbranch4<<","<<Abranch[5]<<","<<stdevAbranch5<<","<<Abranch[6]<<","<<stdevAbranch6<<","<<Abranch[7]<<","<<Abranch[8]<<","<<Abranch[9]<<","<<Abranch[10]<<","<<Abranch[11]<<","<<Abranch[12]<<","<<Abranch[13]<<","<<Abranch[14]<<","<<Abranch[15]<<","<<maxbranch<<","<<perpaircent<<","<<perGCcent<<","<<perAUcent<<","<<perGUcent<<","<<centroid_energy<<","<<MLDcent<<endl;  	*outmma << accnum[seqnum] <<","<<  n <<","<<  pGC <<","<<  pG <<","<<  pC <<","<<  pA <<","<<   pU <<","<<  num <<","<<mfe<<","<< ensemble_fe<<","<<freqmfe<<","<<ensemble_diversity<<","<<Abpdist<<","<<AMLD<<","<< stdevAMLD<<","<< MALD<<","<< stdevMALD<<","<<RNAfold_MALD<<","<<NhalfALD<<","<<stdevNhalfALD<<","<< AALD<<","<<RNAfold_AALD<<","<<perpair<<","<< stdevperpair<<","<<RNAfold_perpair<<","<<perGC<<","<<stdevperGC<<","<<perAU<<","<<stdevperAU<<","<<perGU<<","<<stdevperGU<<","<<Adup1<<","<<stdevAdup1<<","<<Adup2<<","<<stdevAdup2<<","<<Adistss<<","<<stdevAdistss<<","<<Adistds<<","<<stdevAdistds<<","<<Abranch[1]<<","<<stdevAbranch1<<","<<Abranch3p<<","<<stdevAbranch3p<<","<<Abranch4p<<","<<stdevAbranch4p<<","<<Abranch5p<<","<<stdevAbranch5p<<","<<Abranch6p<<","<<stdevAbranch6p<<","<<Abranch15p<<","<<stdevAbranch15p<<","<<Abranch[2]<<","<<stdevAbranch2<<","<<Abranch[3]<<","<<stdevAbranch3<<","<<Abranch[4]<<","<<stdevAbranch4<<","<<Abranch[5]<<","<<stdevAbranch5<<","<<Abranch[6]<<","<<stdevAbranch6<<","<<Abranch[7]<<","<<Abranch[8]<<","<<Abranch[9]<<","<<Abranch[10]<<","<<Abranch[11]<<","<<Abranch[12]<<","<<Abranch[13]<<","<<Abranch[14]<<","<<Abranch[15]<<","<<maxbranch<<","<<perpaircent<<","<<perGCcent<<","<<perAUcent<<","<<perGUcent<<","<<centroid_energy<<","<<MLDcent<<endl;  

	for (k=1; k<=num; k++)
	{
		*out3 << k << "  " << mfe << "  "  << maxLD[k] << "  " << imaxLD[k]
   	  	  << "  " << jmaxLD[k] << "  " << GC[k] << "  " << AU[k] << "  "
           << GU[k] << "\n";
	}
	out3->close();
	delete out3;	
	
//////////////////////////////
	stringstream pairprobstream;
	pairprobstream << "cp pairprob3.txt /Users/Aron/Documents/RAW_FOLDING_DATA/"<<seqfamily<<"/"<<seqfamily<<"_PairProbs_"<<seqnum<<".txt";
	string pairprobstring = pairprobstream.str();
		
	if (std::system(0)) // If a command processor is available.
	{	
		
		std::system(&pairprobstring[0]);
	}
	else
	{
		std::cerr << "A command processor is not available.\n";
	}

//////////////////////////////
	stringstream bracketstream;
	bracketstream << "cp tmp3.b /Users/Aron/Documents/RAW_FOLDING_DATA/"<<seqfamily<<"/"<<seqfamily<<"_Structures_"<<seqnum<<".b";
	string bracketstring = bracketstream.str();
		
	if (std::system(0)) // If a command processor is available.
	{	
		
		std::system(&bracketstring[0]);
	}
	else
	{
		std::cerr << "A command processor is not available.\n";
	}
//////////////////////////////
	stringstream centroidstream;
	centroidstream << "cp partition.txt /Users/Aron/Documents/RAW_FOLDING_DATA/"<<seqfamily<<"/"<<seqfamily<<"_Centroid_"<<seqnum<<".b";
	string centroidstring = centroidstream.str();
		
	if (std::system(0)) // If a command processor is available.
	{	
		
		std::system(&centroidstring[0]);
	}
	else
	{
		std::cerr << "A command processor is not available.\n";
	}
////////////////////////////////////////////	
	
	
// destroy dynamically allocated arrays
	delete[] *LD;
	delete[] LD;
	delete[] *ALD;
	delete[] ALD;
	delete[] *ALDsq;
	delete[] ALDsq;
	delete[] sen;
	delete[] structure1;
	delete[] structure2;	
	

// std::system("rm tmp3*");///WARNING: IF WANT TO DO THIS, NEED TO CHANGE NAME OF TMP3B2CT.SH; OTHERWISE, THIS CRITICAL EXECUTABLE WILL BE REMOVED AND THE PROGRAM WILL ONLY RUN ONCE!	
} // END BIG "FOR" LOOP!!!!

outsum->close();
delete outsum;
outmma->close();
delete outmma;
	return 0;
}