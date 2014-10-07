#include "gibbssampler.h"

using namespace std;

// This function takes log p1, log p2, and returns log(p1+p2)
double log_add(double a, double b) {
  double maximum = max(a,b);
  double minimum = min(a,b);
  if (fabs(maximum - minimum) > 30) {
    return maximum;
  }
  return maximum + log(1 + exp(minimum - maximum));
}

// Samples an event from log probs. The sum of probs can be < 1.
int sample_unnorm(vector<double> &input, double r) {
  double max_log_prob = *(max_element(input.begin(), input.end()));
  vector<double> CDF;
  CDF.reserve(input.size());
  CDF.push_back(input[0]);
  for (int i = 1; i < input.size(); ++i) {
    CDF.push_back(log_add(CDF[CDF.size() - 1], input[i]));
  }
  r += CDF.back();
  
  return lower_bound(CDF.begin(), CDF.end(),r) - CDF.begin();
  
}

bool IsInUnobserved(std::vector<double> &starts, std::vector<double> &ends, double t){
	if(starts.empty())
		return false;
	for(int i = 0; i<starts.size(); i++){
		if(t > starts[i] && t < ends[i])
			return true;
	}
	return false;
}

bool GetPreviousState(map<vector<shptr<generic_state> >, vector<pair<vector<shptr<generic_state> >, pair<double,bool> > >, ssumpcomp> &transmap, vector<shptr<generic_state> > &jointstate, double p){
	//cerr<<"starting..."<<endl;
	bool keep = false;
        vector<double> logprobs;
	auto iter = transmap.find(jointstate);
	logprobs.reserve(iter->second.size());
	for(auto iter1 = iter->second.begin(); iter1!= iter->second.end(); iter1++){
		logprobs.push_back(iter1->second.first);
	}

	auto tempit = iter->second.begin();
	advance(tempit,sample_unnorm(logprobs, p));
	jointstate = tempit ->first;	
	keep = tempit->second.second;

	return keep;
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
	cerr<<"error1!!!!"<<endl;

}

GibbsAuxSampler::GibbsAuxSampler(const pcim *model, const ctbn::Trajectory *evidence, const ctbn::Context *contexts, int burnin) {

	m = model;
	init_traj = NULL;
	evid = evidence;
	numBurninIter = burnin;
	context = contexts;
	burntin = false;
	own_var_list = contexts->VarList();
	testindexes.resize(m->counttest());		
	model->Makeindex(testindexes,0);
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

//now only set init_traj as the cleaned up evid (nothing between -2 and -3)
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
			if(shouldsave && it->second != -2)		
				tr.AddTransition(varid, it->first, it->second);
			if(it->second == -2)
				shouldsave = false;
			if(it->second == -3)
				shouldsave = true;
			it++;
		}
	}
	//should sample in the [-2,-3] intervals, but not necessary - to do
}

void GibbsAuxSampler::Clearcurrentvar(int varid) const{

	tr.SetUnknown(varid,true);
	decltype(evid->GetVarTraj(0).begin()) it, tmpend;
	bool shouldsave = true;

	if (evid->GetVarTraj(varid).empty()) 
		return;		
	it = evid->GetVarTraj(varid).begin();
	tmpend = evid->GetVarTraj(varid).end();
	while(it!=tmpend){
		if(shouldsave && it->second != -2)		
			tr.AddTransition(varid, it->first, it->second);
		if(it->second == -2)
			shouldsave = false;
		if(it->second == -3)
			shouldsave = true;
		it++;
	}

}


void GibbsAuxSampler::ClearInitTraj() {
	if (init_traj) {
		delete init_traj;
		init_traj = NULL;
	}
}

//TODO calculate only once?
void GibbsAuxSampler::GetUnobservedIntervals(int varid) const{

	if (evid->GetVarTraj(varid).empty()) 
		return;
	starts.clear();
	ends.clear();
	decltype(evid->GetVarTraj(varid).begin()) it, tmpend;
	it = evid->GetVarTraj(varid).begin();
	tmpend = evid->GetVarTraj(varid).end();
	while(it!=tmpend){
		if(it->second == -2)
			starts.push_back(it->first);
		if(it->second == -3)
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









