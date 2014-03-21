#ifndef PCIM_GIBBSSAMPLER_H
#define PCIM_GIBBSSAMPLER_H

#include "pcim.h"

using namespace std;

bool IsInUnobserved(std::vector<double> &starts, std::vector<double> &ends, double t);

class GibbsAuxSampler{

public:
	GibbsAuxSampler(const pcim *model, const ctbn::Trajectory *evidence, const ctbn::Context *contexts, int burnin);

	~GibbsAuxSampler();

	void SetTrajectory(const ctbn::Trajectory *traj);

	template<typename R>
	void SampleTrajectories(std::vector<ctbn::Trajectory> &traj, std::vector<double> &w,
					int numsamples, R &rand){	
		BurnIn(rand);
		cout<<"after burn in: "<<endl;
		tr.AddTransition(0, 7, 0);//sampled from previous iteration
		printtr(cout,tr,3);
		for (int i=0; i<numsamples; ++i) {
			traj.push_back(Get());
			w.push_back(0.0);  // log weight
			Next(rand);
			cout<<"after sample: "<<endl;
			printtr(cout,tr,3);
		}
	}

mutable ctbn::Trajectory tr,oldtr;

protected:
	template<typename R>
	void BurnIn(R &rand) const{
		if (burntin) return;
		if (!init_traj) {
			SampleInitialTrajectory(); //now only clean up evid and set it to tr
		} else {
			tr = *init_traj;
		}
		Next(rand, numBurninIter);
		burntin = true;
	}

	// one step in the markov chain state transition; 
	template<typename R>
	void Next(R &rand, int num_iter = 1) const{
		for (int i=0; i<num_iter; ++i) {
			for (size_t var=0; var!=own_var_list.size(); ++var)
				SampleVariable(own_var_list[var], rand);
		}			
	}

	const ctbn::Trajectory &Get() const { return tr; }

	const ctbn::Trajectory &GetInitTraj() const { return *init_traj; }

	// Sample initial trajectory that agrees with evidence *evid.
	void SampleInitialTrajectory() const;
	void GetUnobservedIntervals(int varid) const;
	void GetAuxRates(int varid, int card) const;
	void ClearInitTraj();

	// Resample the entire trajectory of v given all the other variables' full trajectory.
	template<typename R>
	void SampleVariable(int var, R &rand) const{
		ctbn::Context varcontext;
		varcontext.AddVar(var, context->Cardinality(var));
		GetUnobservedIntervals(var);
		GetAuxRates(var, context->Cardinality(var));
		cout<<"auxrates:"<<endl;
		for(int i = 0; i<auxstarts.size(); i++)
		{
			cerr<<auxrates[i]<<" in ( "<<auxstarts[i]<<","<<auxends[i]<<")"<<endl;
		}
		oldtr = tr;
		std::exponential_distribution<> expdist(1.0);
		std::uniform_real_distribution<> unifdist(0.0,1.0);
		std::normal_distribution<> normdist(0.0,1.0);		
	
		for(int i = 0; i < starts.size(); i++){

			//samplecomplete
			//double t = evid->TimeBegin();
			//double T = evid->TimeEnd();
			double t = starts[i];
			double T = ends[i];	
			double lastt=t;
			while((t = m->geteventaux(tr,lastt,expdist(rand),unifdist(rand),normdist(rand),var,T,varcontext,auxstarts,auxends,auxrates))<T) {
				cerr<<"sampled: "<<"var: "<<var<<" t: "<<t<<endl;
				//should add to a container - to do		
				lastt = t; //proceed no matter event kept or not
			}
			//thinning tr - to do
		}//for	
	}

	int numBurninIter;
	const ctbn::Context *context;//contexts of all vars
	std::vector<int> own_var_list;
	mutable std::vector<double> starts; //times when unobserved intervals start, for the current variable
	mutable std::vector<double> ends; //times when unobserved intervals ends, for the current variable
	const pcim *m;
	const ctbn::Trajectory *evid;
	ctbn::Trajectory *init_traj; 
	//mutable ctbn::Trajectory tr,oldtr; 
	double begintime;
	double endtime;
	mutable bool burntin;

	mutable vector<double> auxstarts;
	mutable vector<double> auxends;
	mutable vector<double> auxrates;

private:
};


#endif
