#ifndef AKIPCIM_H
#define AKIPCIM_H

#include "load.h"
#include "pcim.h"

std::vector<shptr<pcimtest>> akitests(const datainfo &info) {
     const std::set<double> lags{0.5, 1.0, 6.0, 12.0, 24.0};
     const std::set<double> threshs{-1.0,-0.5,0.0,0.5,1.0};
     const std::set<int> counts{1,5,10};

	std::vector<shptr<pcimtest>> tests;
	for(int v=-1;v<(int)(info.dvarnames.size());v++) {
		tests.emplace_back(new vartest(v));
		for(double thresh : threshs) {
			tests.emplace_back(new lasttest(v,thresh));
			for(double lag : lags) 
				tests.emplace_back(new meantest(thresh,v,0.0,lag));
		}

		for(int c : counts)
			for(double lag : lags)
				tests.emplace_back(new counttest(c,v,0.0,lag));
	}
	auto stidloc = info.svarid.find(std::string("age"));
	if (stidloc!=info.svarid.end()) {
		int id = stidloc->second;
		tests.emplace_back(new staticgreqtest(1,id));
		tests.emplace_back(new staticgreqtest(6,id));
		tests.emplace_back(new staticgreqtest(12,id));
		tests.emplace_back(new staticgreqtest(12*2,id));
		tests.emplace_back(new staticgreqtest(12*4,id));
		tests.emplace_back(new staticgreqtest(12*8,id));
		tests.emplace_back(new staticgreqtest(12*12,id));
	}
	stidloc = info.svarid.find(std::string("weight"));
	if (stidloc!=info.svarid.end()) {
		int id = stidloc->second;
		tests.emplace_back(new staticgreqtest(-1.0,id));
		tests.emplace_back(new staticgreqtest(-0.0,id));
		tests.emplace_back(new staticgreqtest(1.0,id));
	}
	stidloc = info.svarid.find(std::string("admitcategory"));
	if (stidloc!=info.svarid.end()) {
		int id = stidloc->second;
		tests.emplace_back(new staticeqtest(0,id));
		tests.emplace_back(new staticeqtest(1,id));
		tests.emplace_back(new staticeqtest(2,id));
		tests.emplace_back(new staticeqtest(3,id));
		tests.emplace_back(new staticeqtest(4,id));
		tests.emplace_back(new staticeqtest(5,id));
		tests.emplace_back(new staticeqtest(6,id));
	}
	stidloc = info.svarid.find(std::string("race"));
	if (stidloc!=info.svarid.end()) {
		int id = stidloc->second;
		tests.emplace_back(new staticeqtest(0,id));
		tests.emplace_back(new staticeqtest(1,id));
		tests.emplace_back(new staticeqtest(2,id));
		tests.emplace_back(new staticeqtest(3,id));
		tests.emplace_back(new staticeqtest(4,id));
		tests.emplace_back(new staticeqtest(5,id));
		tests.emplace_back(new staticeqtest(6,id));
		tests.emplace_back(new staticeqtest(7,id));
		tests.emplace_back(new staticeqtest(8,id));
	}
	
	stidloc = info.svarid.find(std::string("starttod"));
	if (stidloc!=info.svarid.end()) {
		int id = stidloc->second;
		tests.emplace_back(new timetest(55.0/60.0,05.0/60.0,1,id));
		tests.emplace_back(new timetest(05.0/60.0,25.0/60.0,1,id));
		tests.emplace_back(new timetest(25.0/60.0,35.0/60.0,1,id));
		tests.emplace_back(new timetest(35.0/60.0,55.0/60.0,1,id));
	}
	return tests;
}


pcim::pcimparams akipcimparams(int ntests, int nproc=1) {
	// parameter order: alpha, beta (for duration gamma prior)
	// 				kappa (for exponential structure size prior)
	//                  mu, kappa, alpha, beta (for normal-gamma value prior)
	//                  minimum number of events in leaves
	//                  num of threads to use
	//pcim::pcimparams p(1,1,1.0/tests.size(),0.0,1.0,2.0,2.0,0,nproc);
#ifdef USEPERSIST
	pcim::wtT mu;
	pcim::wtvarT kappa;
	double alpha,beta;
	mu << 0.0, 1.0;
	kappa << 1.0, 0.0, 0.0, 1.0; // identity (why not?)
	alpha = 2.0;
	beta = 1.0;
	pcim::pcimparams p(1,1,1.0/ntests,
			mu,kappa,alpha,beta,
			1000,nproc);
#else
	pcim::pcimparams p(1,1,1.0/ntests,0.0,1.0,1.5,0.5,1000,nproc);
#endif
	return p;
}

#endif
