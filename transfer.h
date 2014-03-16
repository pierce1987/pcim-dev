#include "pcim.h"

using namespace std;

pcim* CTBNtransfer(istream& is);

//on the CTBN side. use the following code after a ctbndyn of type CTBNDyn
/* 

cout<<ctbndyn.NumofNodes()<<endl;
for(int p = 0; p<ctbndyn.NumofNodes(); p++){

	const Dynamics *bnode = ctbndyn.Node(p);
	int var = bnode->Domain().VarList()[0];
	cout<<var<<endl;

	int ds = ((MarkovDyn*)bnode)->Domain().Size();
	int cds = ((MarkovDyn*)bnode)->CondDomain().Size();
	Context contexts = ((MarkovDyn*)bnode)->CondDomain();
	cout<<ds<<endl;
	cout<<contexts.VarList().size()<<endl;

	for(int i=0; i<contexts.VarList().size(); i++){
		int id = contexts.VarList()[i];
		int card = contexts.Cardinality(id);
		cout<<id<<endl;
		cout<<card<<endl;
	}
	
	//cout<<"ds: "<<ds<<endl;
	//cout<<"cds: "<<cds<<endl;
	int i = 0;
	for (int i=0; i<cds; i++){
		//((MarkovDyn*)bnode)-> CondDomain().Index(i)

		matrix &m = ((MarkovDyn*)bnode)->operator()(((MarkovDyn*)bnode)-> CondDomain().Index(i))->Intensity();
 		
 		for(int j=0; j<ds; j++){
 			for(int k=0; k<ds; k++){
 				cout<< m[j][k] <<endl ;
 			}
 		//cout<<endl;	
 		}
	//cout<<endl<<endl;
	}		
	
}
*/
/*example:
var0 and var1 are both binary with var0 <-> var1
for var0:
when var1=0, Q = [-1 1 ; 2 -2]
when var1=1, Q = [-10 10 ; 20 -20]

for var1:
when var0=0, Q = [-5 5; 6 -6]
when var0=1, Q = [-7 7; 8 -8]

Transferred PCIM model:
if var == 0 [0,0]
   if most recent 1 == 0 [0,0]
      if most recent 0 == 0 [0,0]
         if event == (0,0) [0,0]
            rate = 0[0,0]
         else
            rate = 1[0,0]
      else
         if event == (0,0) [0,0]
            rate = 2[0,0]
         else
            rate = 0[0,0]
   else
      if most recent 0 == 0 [0,0]
         if event == (0,0) [0,0]
            rate = 0[0,0]
         else
            rate = 10[0,0]
      else
         if event == (0,0) [0,0]
            rate = 20[0,0]
         else
            rate = 0[0,0]
else
   if most recent 0 == 0 [0,0]
      if most recent 1 == 0 [0,0]
         if event == (1,0) [0,0]
            rate = 0[0,0]
         else
            rate = 5[0,0]
      else
         if event == (1,0) [0,0]
            rate = 6[0,0]
         else
            rate = 0[0,0]
   else
      if most recent 1 == 0 [0,0]
         if event == (1,0) [0,0]
            rate = 0[0,0]
         else
            rate = 7[0,0]
      else
         if event == (1,0) [0,0]
            rate = 8[0,0]
         else
            rate = 0[0,0]
*/



