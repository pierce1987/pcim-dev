#include "gibbssampler.h"

using namespace std;

bool IsInUnobserved(std::vector<double> &starts, std::vector<double> &ends, double t){
	if(starts.empty())
		return false;
	for(int i = 0; i<starts.size(); i++){
		if(t > starts[i] && t < ends[i])
			return true;
	}
	return false;
}


double GibbsAuxSampler::getnextevent(double t0, int &event) const{
	double min_time = 10000000000;
	double old_min = min_time;
	map<double, int>::const_iterator tempitr;
	for (size_t var=0; var!=own_var_list.size(); ++var){
		tempitr = oldtr.GetVarTraj(own_var_list[var]).upper_bound(t0);
		if( tempitr != oldtr.GetVarTraj(own_var_list[var]).end()){
			if(tempitr->first < min_time){
				min_time = tempitr->first;
				event = own_var_list[var];
			}
		}
	}

	if(min_time == old_min)
		return -1.0;

	else
		return min_time;
}

bool GibbsAuxSampler::IsVirtual(double t0, int event, int varid) const{
	if(event != varid)
		return false;
	for(int i=0; i<starts.size(); i++){
		if(t0>=starts[i] && t0<=ends[i])
			return true;
	}
	return false;
}

double GibbsAuxSampler::Getkeepprob(double rate, double t0) const{
	for(int i = 0; i<auxstarts.size(); i++){
		if(t0 > auxstarts[i] && t0 < auxends[i])
			return rate/auxrates[i];
	}
	cerr<<"error!!!!"<<endl;

}

GibbsAuxSampler::GibbsAuxSampler(const pcim *model, const ctbn::Trajectory *evidence, const ctbn::Context *contexts, int burnin) {

	m = model;
	init_traj = NULL;
	evid = evidence;
	numBurninIter = burnin;
	context = contexts;
	burntin = false;
	own_var_list = contexts->VarList();
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
			if(shouldsave && it->second != -1)		
				tr.AddTransition(varid, it->first, it->second);
			if(it->second == -1)
				shouldsave = false;
			if(it->second == -2)
				shouldsave = true;
			it++;
		}
	}
	//should sample in the [-1,-2] intervals, but not necessary - to do
}


void GibbsAuxSampler::ClearInitTraj() {
	if (init_traj) {
		delete init_traj;
		init_traj = NULL;
	}
}

void GibbsAuxSampler::GetUnobservedIntervals(int varid) const{

	if (evid->GetVarTraj(varid).empty()) 
		return;
	starts.clear();
	ends.clear();
	decltype(evid->GetVarTraj(varid).begin()) it, tmpend;
	it = evid->GetVarTraj(varid).begin();
	tmpend = evid->GetVarTraj(varid).end();
	while(it!=tmpend){
		if(it->second == -1)
			starts.push_back(it->first);
		if(it->second == -2)
			ends.push_back(it->first);
		it++;
	}
}

//fill auxstarts auxends and auxrates;
void GibbsAuxSampler::GetAuxRates(int varid, int card) const{

	auxstarts.clear();
	auxends.clear();
	auxrates.clear();

	for(int i = 0; i < starts.size(); i++){
		double t = starts[i];
		double T = ends[i];
		double lastt=t;
		double until;
		double r;

		while((t = m->getauxrates(tr,lastt,card,until,r,varid))<T) {
			if(until >= T){ until = T; break;}
			//cerr<<"here: "<<"r: "<<r<<" t: "<<t<<" until: "<<until<<endl;
			//merge intervals when we can
			if(!auxrates.empty() && abs(2*r - auxrates[auxrates.size()-1])<0.00001 && abs(t - auxends[auxends.size()-1])<0.00001)
				auxends[auxends.size()-1] = until;
			else{
				auxstarts.push_back(t);
				auxends.push_back(until);
				auxrates.push_back(2*r);
			}
			lastt = until; //proceed no matter event kept or not
		}
		//cerr<<"here: "<<"r: "<<r<<" t: "<<t<<" until: "<<until<<endl;
		if(!auxrates.empty() && abs(2*r - auxrates[auxrates.size()-1])<0.00001 && abs(t - auxends[auxends.size()-1])<0.00001)
			auxends[auxends.size()-1] = until;
		else{
			auxstarts.push_back(t);
			auxends.push_back(until);
			auxrates.push_back(2*r);
		}
	}
	
}

//thinning, performed on oldtr







