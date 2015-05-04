#include "pcim.h"
#include "transfer.h"
#include "gibbssampler.h"
#include "em.h"
#include <iostream>
#include <vector>
#include <random>
#include <sstream>
#include <chrono>

using namespace std;
using namespace std::chrono;

#define NPROC 1
//#define NPROC 4

//LD_LIBRARY_PATH=/home/rlair/lib:/usr/csshare/pkgs/gcc-4.8.1/lib:/usr/csshare/pkgs/gcc-4.8.1/lib64

//export LD_LIBRARY_PATH


// main for ising model
int main(int argc, char **argv) {
	int nsamp = argc>1 ? atoi(argv[1]) : 10;
	nsamp = 20;
/*
	pcim truemodel(new varcounttest(1, -1, 0, 1),
				new pcim(new varcounttest(1,-1,1,2),
					new pcim(0.01),
					new pcim(2.0)),
				new pcim(new vartest(0),
					new pcim(0.5),
					new pcim(new varcounttest(1,0,0,5),
						new pcim(0.1),
						new pcim(0.01))));

	truemodel.print(cout);
	cout<<endl;
*/


	//ifstream fs;
	//fs.open("ising.txt");
	//pcim *truemodel = CTBNtransfer(fs);
	//truemodel->print(cout); cout<<endl;
	//fs.close();

set<int> othervar0;
	othervar0.insert(0);

set<int> othervar1;
	othervar1.insert(1);

	pcim truemodel(new vartest(0),
				new pcim(new varcounttest(1,0,0,0.5),
					new pcim(10.0),	
					new pcim(new varcounttest(1,0,0,1),
							new pcim(new lastvartest(1,othervar0),
								new pcim(5.0),
								new pcim(1.0)),
							new pcim(1.0))),
				new pcim(new lastvartest(1,othervar0),
					new pcim(new varcounttest(1,0,0,0.5),
							new pcim(1.0),
							new pcim(5.0)),
					new pcim(new timetest(0,0.5,1),
							new pcim(1.0),
							new pcim(10.0))));

	truemodel.print(cout);
	cout<<endl;


/*	set<int> othervar1;
	othervar1.insert(1);
	pcim truemodel(new lastvartest(0, othervar1),
				new pcim(0.5),
				new pcim(1));	

	truemodel.print(cout);
	cout<<endl;
*/
	random_device rd;

	int nvar = 2;

	ctbn::Context contexts;
	contexts.AddVar(0, 1);
	contexts.AddVar(1, 1);



	ctbn::Trajectory tr = ctbn::Trajectory();
	tr.SetEndTime(5.0);
	tr.AddTransition(0, 0.4, 0);	
	tr.AddTransition(0, 0.6, 0);
	tr.AddTransition(0, 1.8, 0);
	tr.AddTransition(0, 4.7, 0);

	tr.AddTransition(1, 0.1, 0);
	tr.AddTransition(1, 0.2, 0);
	tr.AddTransition(1, 3.4, 0);
	tr.AddTransition(1, 3.6, 0);
	tr.AddTransition(1, 3.7, 0);

	tr.AddTransition(0, 2.0, -2);
	tr.AddTransition(0, 4.0, -3);
	tr.AddTransition(1, 1.0, -2);
	tr.AddTransition(1, 3.0, -3);

	int N = 50;
	cerr<<"N:"<<N<<endl;
	for (int round = 1; round<=100; round++) {
	unsigned int seed = rd();
	//cout << "seed = " << seed << endl;
	mt19937 randgen(seed);
	vector<ctbn::Trajectory> t;//vector of sampled trajectories
	vector<double> w;

	GibbsAuxSampler sampler(&truemodel, &tr, &contexts, N/10); //tr is evidence, last param is burn-in round
	sampler.SampleTrajectories(t,w,N,randgen);//3rd param: # of samples wanted

	int countA = 0;
	int countB = 0;

	for (ctbn::Trajectory traj : t) {

			countA += traj.GetVarTraj(0).size();
			countB += traj.GetVarTraj(1).size();

	}

	cerr<<((double)countA/N)<<" "<<((double)countB/N)<<endl;

	}






/*

	vector<ctbn::Trajectory> data;
	vector<ctbn::Trajectory> holdout;
	for (int i = 0; i < nsamp; ++i) {
		ctbn::Trajectory tr = truemodel.sample(10.0, nvar, randgen, contexts);
		//printtr(cout,tr,2);
		data.push_back(tr);
	}
	for (int i = 0; i < 100; ++i) {
		ctbn::Trajectory tr = truemodel.sample(10.0, nvar, randgen, contexts);
		//printtr(cout,tr,2);
		holdout.push_back(tr);
	}
	//cout<<truemodel.calcDataLikelihood(data, contexts)<<endl;

	pcim::pcimparams p(1,1,0.1,0,NPROC);
	//ParaEstimateEM(truemodel, p, contexts, data, 0.4, randgen, holdout);


*/

	// for structre EM
/*	vector<shptr<pcimtest>> tests;
	for (int v=-1; v<nvar; v++) {
		tests.emplace_back(new vartest(v));
		for(double t0 = 0.0; t0<=1.0; t0+=1.0) {
			for (double t1 = 1.0; t1<= 6.0; t1+=1.0) {
				for (int i = 1; i < 3; i++) {
					tests.emplace_back(new varcounttest(i,v,t0,t1));
				}
			}
		}
	}*/
	/*tests.emplace_back(new varcounttest(1,0,0,5));
	tests.emplace_back(new varcounttest(1,-1,1,2));
	tests.emplace_back(new varcounttest(1,-1,0,1));
	tests.emplace_back(new vartest(0));
	tests.emplace_back(new vartest(1));
	*/
	//StructureEM(p, tests, contexts, data, 0.6, randgen, holdout);


	//gamma_distribution<double> dis(1.0,1.0);
	//for (int i = 0; i < 100; ++i) {
	//	cerr<<" "<<dis(randgen);
	//}
	











/*
	//ctbn::Trajectory tr = ctbn::Trajectory();
	//tr.SetEndTime(1.2);

	vector<ctbn::Trajectory> t;//vector of sampled trajectories
	vector<double> w;

	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	GibbsAuxSampler sampler(truemodel, &tr, &contexts, 0); //tr is evidence, last param is burn-in round
	sampler.SampleTrajectories(t,w,5,randgen);//3rd param: # of samples wanted


	// get result for each var
	vector<int> countzero(9,0);
	vector<int> countone(9,0);
	//cerr<<"size: "<<t.size()<<endl;
	for (ctbn::Trajectory traj : t) {
		for (int i = 0; i <= 8; i++) {
			auto it = traj.GetVarTraj(i).upper_bound(0.6);
			it--;
			if (it->second == 0) {
				countzero[i]++;
			} 
			if (it->second == 1) {
				countone[i]++;
			} 
		}
	}

	cerr<<"countzero:"<<endl;
	for (int i : countzero) {
		cerr<<i<<" ";
	}
	cerr<<endl;

	cerr<<"countone:"<<endl;
	for (int i : countone) {
		cerr<<i<<" ";
	}
	cerr<<endl;

	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::microseconds>(t2 - t1).count();
	cout<<duration<<endl;
*/

}
