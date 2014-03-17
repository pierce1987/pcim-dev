#include "gibbssampler.h"

using namespace std;

GibbsAuxSampler::GibbsAuxSampler(const pcim *model, const ctbn::Trajectory *evidence, const ctbn::Context *contexts, int burnin) {

	m = model;
	init_traj = NULL;
	evid = evidence;
	numBurninIter = burnin;
	context = contexts;
	burntin = false;
	//Initialize();
}

GibbsAuxSampler::~GibbsAuxSampler() {
 	if (init_traj) delete init_traj;
}

void GibbsAuxSampler::SetTrajectory(const ctbn::Trajectory *evidence) { 
	evid = evidence;
	begintime = evid->TimeBegin();
	endtime = evid->TimeEnd();
	burntin=false;
}

//now only set init_traj as the cleaned up evid (nothing between -1 and -2)
void GibbsAuxSampler::SampleInitialTrajectory() const {
	tr = ctbn::Trajectory();
	tr.SetBeginTime(evid->TimeBegin());
	tr.SetEndTime(evid->TimeEnd());
	decltype(evid->GetVarTraj(0).begin()) it, tmpend;
	bool shouldsave = true;

	for(int i=0;i<context->VarList().size();i++){
		int varid = context->VarList()[i];
		if (evid->GetVarTraj(varid).empty()) 
			continue;		
		it = evid->GetVarTraj(varid).begin();
		tmpend = evid->GetVarTraj(varid).end();
		while(it!=tmpend){
			if(shouldsave)		
				tr.AddTransition(varid, it->first, it->second);
			if(it->second == -1)
				shouldsave = false;
			if(it->second == -2){
				shouldsave = true;
				tr.AddTransition(varid, it->first, -2);
			}
			it++;
		}
	}
}


void GibbsAuxSampler::ClearInitTraj() {
	if (init_traj) {
		delete init_traj;
		init_traj = NULL;
	}
}










