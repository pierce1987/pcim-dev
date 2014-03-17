#ifndef PCIM_GIBBSSAMPLER_H
#define PCIM_GIBBSSAMPLER_H

#include "pcim.h"

class GibbsAuxSampler{

public:
	GibbsAuxSampler(const pcim *model, const ctbn::Trajectory *evidence, const ctbn::Context *contexts, int burnin);

	~GibbsAuxSampler();

	void SetTrajectory(const ctbn::Trajectory *traj);

	template<typename R>
	void SampleTrajectories(std::vector<ctbn::Trajectory> &traj, std::vector<double> &w,
					int numsamples, R &rand){	
		BurnIn(rand);
		for (int i=0; i<numsamples; ++i) {
			traj.push_back(Get());
			w.push_back(0.0);  // log weight
			//Next();
		}
	}

//mutable ctbn::Trajectory tr,oldtr;

protected:
	template<typename R>
	void BurnIn(R &rand) const{
		if (burntin) return;
		if (!init_traj) {
			SampleInitialTrajectory(); //now only clean up evid and set it to tr
		} else {
			tr = *init_traj;
		}
		//Next(numBurninIter, rand);
		burntin = true;
	}

	//void Next(int num_iter = 1, Random &rand = randomizer) const;

	const ctbn::Trajectory &Get() const { return tr; }

	const ctbn::Trajectory &GetInitTraj() const { return *init_traj; }

	// Sample initial trajectory that agrees with evidence *evid.
	void SampleInitialTrajectory() const;

	void ClearInitTraj();

	int numBurninIter;
	const ctbn::Context *context;//contexts of all vars
	const pcim *m;
	const ctbn::Trajectory *evid;
	ctbn::Trajectory *init_traj; 
	mutable ctbn::Trajectory tr,oldtr; 
	double begintime;
	double endtime;
	mutable bool burntin;

private:
};


#endif
