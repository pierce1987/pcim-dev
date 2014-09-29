#include "pcim.h"
#include "transfer.h"
#include "gibbssampler.h"
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
	nsamp = 1;
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


	pcim truemodel(new timetest(1,4,5),
				new pcim(2.0),
				new pcim(0.01));

	pcim truemodel(new lasttest(-1, 1),
				new pcim(new lasttest(1, 1),new pcim(10.0),new pcim(2)),
				new pcim(new eventtest(-1, 1),new pcim(3),new pcim(4)));
*/

//	pcim truemodel(new varcounttest(1, 0, 0.0, 5.0),				
//				new pcim(0.1),				
//				new pcim(new varcounttest(1, 1, 0.0, 1.0),new pcim(5),new pcim(0.3)));

	pcim truemodel(new varcounttest(1, 1, 0.0, 1.0),				
				new pcim(5),				
				new pcim(new varcounttest(1, 0, 0.0, 5.0),new pcim(0.5),new pcim(0.3)));

	truemodel.print(cout); cout << endl;
	random_device rd;

	int nvar = 2;

	ctbn::Context contexts;
	contexts.AddVar(0, 1);
	contexts.AddVar(1, 1);
	//contexts.AddVar(2, 1);

	unsigned int seed = rd();
	cout << "seed = " << seed << endl;
	mt19937 randgen(seed);
/*	vector<ctbn::Trajectory> data;

	for(int i=0;i<nsamp;i++) {
		ctbn::Trajectory tr = truemodel.sample(15.0,nvar,randgen,contexts);
		printtr(cout,tr,3);
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

//testing sampler
	//let data[0] be the evidence, model be the pcim model
	data[0].AddTransition(0, 2.0, -1);//stopped observing
	data[0].AddTransition(0, 4.0, -2);//starts observing again
	data[0].AddTransition(0, 10.0, -1);//stopped observing
	data[0].AddTransition(0, 12.0, -2);//starts observing again
	data[0].AddTransition(1, 2.1, -1);//stopped observing
	data[0].AddTransition(1, 8.1, -2);//starts observing again
*/

	ctbn::Trajectory tr = ctbn::Trajectory();
	/*tr.AddTransition(0, 1, 0);
	tr.AddTransition(0, 3, 1);
	tr.AddTransition(0, 4, -1);
	tr.AddTransition(0, 6.9, -2);
	tr.AddTransition(0, 7.1, -1);
	tr.AddTransition(0, 8, -2);
	//tr.AddTransition(0, 0.7, 0);
	tr.AddTransition(1, 2, 0);
	tr.AddTransition(1, 6.01, 1);
	tr.AddTransition(1, 4, -1);
	tr.AddTransition(1, 5, -2);*/

	tr.AddTransition(0, 1, -2);
	tr.AddTransition(0, 10, -3);
	//tr.AddTransition(0, 1, 0);
	tr.AddTransition(1, 5, 0);
	//tr.AddTransition(0, 1, 0);
	//tr.AddTransition(0, 2, 0);
	/*tr.AddTransition(1, 7.0, 0);
	tr.AddTransition(1, 7.1, 0);
	tr.AddTransition(1, 7.2, 0);
	tr.AddTransition(1, 7.3, 0);
	tr.AddTransition(1, 7.4, 0);
	tr.AddTransition(1, 7.5, 0);
tr.AddTransition(1, 7.6, 0);
tr.AddTransition(1, 7.7, 0);
tr.AddTransition(1, 7.8, 0);
tr.AddTransition(1, 7.9, 0);*/
        //tr.AddTransition(1, 10, 0);

	vector<ctbn::Trajectory> t;//vector of sampled trajectories
	vector<double> w;
	GibbsAuxSampler sampler(&truemodel, &tr, &contexts, 100); //tr is evidence, last param is burn-in round
	sampler.SampleTrajectories(t,w,10,randgen);//3rd param: # of samples wanted


	//printtr(cout,sampler.tr,3);



	////////Given a ctbndyn as input, generate a pcim
/*	pcim * PCIMfromCTBN = CTBNtransfer(cin);

	PCIMfromCTBN->print(cout);

//or read from a file
	ifstream fs;
	fs.open("1.txt");
	pcim * PCIMfromCTBN = CTBNtransfer(fs);

	PCIMfromCTBN->print(cout);

	fs.close();
*/
}
