#include "pcim.h"
#include <iostream>
#include <vector>
#include <random>
#include <sstream>
#include <numeric>

using namespace std;

#define NPROC 1
//#define NPROC 4

//LD_LIBRARY_PATH=/home/rlair/lib:/usr/csshare/pkgs/gcc-4.8.1/lib:/usr/csshare/pkgs/gcc-4.8.1/lib64

//export LD_LIBRARY_PATH

int main(int argc, char **argv) {
	int nsamp = argc>1 ? atoi(argv[1]) : 10;
	cout << nsamp << endl;
/*
	pcim truemodel(new counttest(2,0,0.0,2.0),
					new pcim(new counttest(2,-1,2.0,4.0),
							new pcim(2.0,-1.0,2.0),
							new pcim(2.0,-2.0,0.1)),
					new pcim(new meantest(-0.5,1,0.0,2.0),
						new pcim(0.5,0.0,1.0),
						new pcim(new eventtypetest(1),
							new pcim(3.0,0.0,5.0),
							new pcim(1.0,0.0,5.0))));
*/
	pcim truemodel(new timetest(2,4,10),
				new pcim(20.0,0.0,1.0),
				new pcim(0.01,10.0,1.0));
	truemodel.print(cout); cout << endl;
	random_device rd;

	int nvar = 3;
	int states[] = {2,1,1};
	vector<int> statenum (states, states + sizeof(states) /sizeof(int));
	if(statenum.size() != nvar)
		cout<<"error"<<endl;
	int nevent = accumulate( statenum.begin(), statenum.end(), 0);

	map<pair<int, int>, int> EventIndexMap;
	int counter = 0;
	for(int i=0; i<nvar; i++){
		for(int j=0; j<statenum[i]; j++){
			EventIndexMap.insert(map<pair<int, int>, int>::value_type(make_pair(i, j), counter));
			counter++;
		}
	}

	unsigned int seed = rd();
	cout << "seed = " << seed << endl;
	mt19937 randgen(seed);
	vector<traj> data;

	for(int i=0;i<nsamp;i++) {
		traj tr = truemodel.sample(100.0,nevent,randgen);
		//printtr(cout,tr);
		data.push_back(tr);
	}
	//for(auto &x : data) printtr(cout,x);

	cout << "done sampling" << endl;

	vector<shptr<pcimtest>> tests;
	for(int v=-1;v<nevent;v++) {
		tests.emplace_back(new eventtypetest(v));
		for(double t0=0.0;t0<=4.0;t0+=2.0)
			for(double t1=0.0;t1<t0;t1+=2.0) {
				for(int i=1;i<3;i++)
					tests.emplace_back(new counttest(i,v,t0,t1));
				for(double t=-0.5;t<=0.5;t+=0.5) 
					tests.emplace_back(new meantest(t,v,t0,t1));
			}
	}
	tests.emplace_back(new timetest(1,4,5));
	tests.emplace_back(new timetest(0,3,10));
	tests.emplace_back(new timetest(2,4,10));
	tests.emplace_back(new timetest(5,8,10));
	tests.emplace_back(new timetest(15,30,60));
	tests.emplace_back(new timetest(0,15,60));
	tests.emplace_back(new timetest(0,45,60));

	pcim::pcimparams p(1,1,1.0/tests.size(),0.0,1.0,2.0,2.0,0,NPROC);

	pcim model(data,tests,p);
	model.print(cout);
	cout << endl;
}
