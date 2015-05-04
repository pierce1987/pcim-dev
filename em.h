#ifndef PCIM_EM_H
#define PCIM_EM_H

#include "pcim.h"
#include "gibbssampler.h"

using namespace std;


void ClearedTraj(vector<ctbn::Trajectory> &data,  const ctbn::Context &contexts, vector<ctbn::Trajectory> &cleardata) {
	decltype(data[0].GetVarTraj(0).begin()) it, tmpend;
	for (ctbn::Trajectory &tr : data) { 
		ctbn::Trajectory newdata = ctbn::Trajectory(contexts.VarList().size());
		newdata.SetBeginTime(tr.TimeBegin());
		newdata.SetEndTime(tr.TimeEnd());
		bool shouldsave = true;
		for(int i=0;i<contexts.VarList().size();i++){
			int varid = contexts.VarList()[i];
			if (tr.GetVarTraj(varid).empty()) 
				continue;		
			it = tr.GetVarTraj(varid).begin();
			tmpend = tr.GetVarTraj(varid).end();
			while(it!=tmpend){
				if(shouldsave && it->second != -2)		
					newdata.AddTransition(varid, it->first, it->second);
				if(it->second == -2)
					shouldsave = false;
				if(it->second == -3)
					shouldsave = true;
				it++;
			}
		}
		cleardata.push_back(newdata);
	}	
}

template<typename R>
void ParaEstimateEM(pcim &model, const pcim::pcimparams &params, const ctbn::Context &contexts, vector<ctbn::Trajectory> &data, double prop, R &rand, vector<ctbn::Trajectory> &holdout) {
	
	cerr<<"true model: "<<endl;
	cerr<<"LL: "<<model.calcDataLikelihood(holdout, contexts)<<endl;

	uniform_real_distribution<> unifdist(0.0,1.0);

	// generate unobserved intervals		
	for (ctbn::Trajectory &tr : data) { 
		double duration = tr.TimeEnd() * prop;	
		//cerr<<"ending time1: "<<tr.TimeEnd()<<endl;
		for (int i = 0; i < contexts.VarList().size(); ++i) {
			double p = unifdist(rand);
			double start = (tr.TimeEnd() - duration) * p;
			double ending = start + duration;	
			tr.AddTransition(i,start,-2);
			tr.AddTransition(i,ending,-3);		
		}

		//printtr(cout,tr,2);		
	}
	
	// generate cleared trajs for initial model
	vector<ctbn::Trajectory> cleardata;
	ClearedTraj(data, contexts, cleardata);
	model.updateparams(cleardata, contexts, params);
	cerr<<"model without filling in: "<<endl;
	cerr<<"LL: "<<model.calcDataLikelihood(holdout, contexts)<<endl;
	//model.print(cout);

	// intialize model paramters
	//model.resetparams(params, rand);
	//cerr<<"after reset parameters"<<endl;
	//cerr<<"LL: "<<model.calcDataLikelihood(holdout, contexts)<<endl;
	//model.print(cout);

	vector<ctbn::Trajectory> t;
	vector<double> w;

	for (int itr = 0; itr < 10; itr++) {
		// Given model, generate samples 
		t.clear();
		w.clear();
		for (ctbn::Trajectory &tr : data) {
			//printtr(cout,tr,2);

			GibbsAuxSampler sampler(&model, &tr, &contexts, 10);
			sampler.SampleTrajectories(t,w,2,rand);
		}

		// Given trajs, estimate model parameters
		model.updateparams(t, contexts, params);
		cerr<<"LL: "<<model.calcDataLikelihood(holdout, contexts)<<endl;
		//model.print(cout);
	
	}
	cerr<<"after EM:"<<endl;
	model.print(cout);
	return;

}


// Structure EM
template<typename R>
void StructureEM(const pcim::pcimparams &params, vector<shptr<pcimtest>> &tests, const ctbn::Context &contexts, vector<ctbn::Trajectory> &data, double prop, R &rand, vector<ctbn::Trajectory> &holdout) {
	uniform_real_distribution<> unifdist(0.0,1.0);
		
	pcim model2(data, tests, params, contexts);
	cerr<<"model with complete data: "<<endl;
	cerr<<"LL: "<<model2.calcDataLikelihood(holdout, contexts)<<endl;
	model2.print(cout);

	// generate unobserved intervals		
	for (ctbn::Trajectory &tr : data) { 
		double duration = tr.TimeEnd() * prop;	
		for (int i = 0; i < contexts.VarList().size(); ++i) {
			double p = unifdist(rand);
			double start = (tr.TimeEnd() - duration) * p;
			double ending = start + duration;	
			tr.AddTransition(i,start,-2);
			tr.AddTransition(i,ending,-3);		
		}

		//printtr(cout,tr,2);		
	}
	// generate cleared trajs for initial model
	vector<ctbn::Trajectory> cleardata;
	ClearedTraj(data, contexts, cleardata);
	//for (int i = 0; i < 5; i++) {
	//	printtr(cout, cleardata[i], 2);
	//}
	// Get initial model
	pcim model(cleardata, tests, params, contexts);
	cerr<<"initial model: "<<endl;
	cerr<<"LL: "<<model.calcDataLikelihood(holdout, contexts)<<endl;
	model.print(cout);

	//model.learnNewModel(data, tests, params, contexts);
	//cerr<<"test: "<<endl;
	//model.print(cout);

	vector<ctbn::Trajectory> t;
	vector<double> w;	
	for (int itr = 0; itr < 10; itr++) {
		cerr<<"iteration...."<<endl;
		// Given model, generate samples 
		t.clear();
		w.clear();
		for (ctbn::Trajectory &tr : data) {
			//printtr(cout,tr,2);

			GibbsAuxSampler sampler(&model, &tr, &contexts, 5); 
			sampler.SampleTrajectories(t,w,1,rand);
		}

		// Given trajs, estimate model parameters
		model.learnNewModel(t, tests, params, contexts);
		cerr<<"LL: "<<model.calcDataLikelihood(holdout, contexts)<<endl;
		model.print(cout);
	
	}

	cerr<<"after EM:"<<endl;
	
	model.print(cout);
}


#endif
