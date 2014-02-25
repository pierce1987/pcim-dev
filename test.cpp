#include "pcim.h"
#include "transfer.h"
#include <iostream>
#include <vector>
#include <random>
#include <sstream>

using namespace std;

#define NPROC 1
//#define NPROC 4

//LD_LIBRARY_PATH=/home/rlair/lib:/usr/csshare/pkgs/gcc-4.8.1/lib:/usr/csshare/pkgs/gcc-4.8.1/lib64

//export LD_LIBRARY_PATH

int main(int argc, char **argv) {
	int nsamp = argc>1 ? atoi(argv[1]) : 10;
	//cout << nsamp << endl;
/*
	pcim truemodel(new counttest(2,0,0.0,2.0),
					new pcim(new counttest(2,-1,2.0,4.0),
							new pcim(2.0,-1.0,2.0),
							new pcim(2.0,-2.0,0.1)),
					new pcim(new meantest(-0.5,1,0.0,2.0),
						new pcim(0.5,0.0,1.0),
						new pcim(new vartest(1),
							new pcim(3.0,0.0,5.0),
							new pcim(1.0,0.0,5.0))));

	pcim truemodel(new timetest(1,4,5),
				new pcim(new eventcounttest(1, 1, 0, 1, 2),new pcim(20.0),new pcim(1.15)),
				new pcim(new varcounttest(1,0),new pcim(0.01),new pcim(1)));

*/
/*	pcim truemodel(new eventcounttest(1, 1, 0, 1, 2),
				new pcim(20.0),
				new pcim(0.01));
	truemodel.print(cout); cout << endl;
	random_device rd;

	int nvar = 3;

	ctbn::Context contexts;
	contexts.AddVar(0, 3);
	contexts.AddVar(1, 3);
	contexts.AddVar(2, 3);

	unsigned int seed = rd();
	cout << "seed = " << seed << endl;
	mt19937 randgen(seed);
	vector<ctbn::Trajectory> data;

	for(int i=0;i<nsamp;i++) {
		ctbn::Trajectory tr = truemodel.sample(100.0,nvar,randgen,contexts);
		//printtr(cout,tr);
		data.push_back(tr);
	}
	//for(auto &x : data) printtr(cout,x);
	cout << "done sampling" << endl;

	vector<shptr<pcimtest>> tests;
	for(int v=-1;v<nvar;v++) {
		tests.emplace_back(new vartest(v));
		for(double t0=0.0;t0<=4.0;t0+=2.0)
			for(double t1=0.0;t1<t0;t1+=2.0) {
				for(int i=1;i<3;i++)
					tests.emplace_back(new varcounttest(i,v,t0,t1));
			}
	}
	tests.emplace_back(new timetest(1,4,5));
	tests.emplace_back(new timetest(0,3,10));
	tests.emplace_back(new timetest(2,4,10));
	tests.emplace_back(new timetest(5,8,10));
	tests.emplace_back(new timetest(15,30,60));
	tests.emplace_back(new timetest(0,15,60));
	tests.emplace_back(new timetest(0,45,60));
	tests.emplace_back(new eventcounttest(1, 1, 0, 1, 2));
	tests.emplace_back(new varcounttest(1,0));

	pcim::pcimparams p(1,1,1.0/tests.size(),0,NPROC);

	pcim model(data,tests,p,contexts);
	model.print(cout);
	cout << endl;
*/

	////////Given a ctbndyn as input, generate a pcim
	int NumofNodes;
	cin>>NumofNodes;

	pcim * PCIMfromCTBN = transfer(NumofNodes, 1);

	PCIMfromCTBN->print(cout);

}
