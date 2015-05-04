#include "pcim.h"
#include <algorithm>
#include <cmath>
#include <mutex>
#include <future>
#include "serial.h"

using namespace std;

double getactualrate(double r, double t, double until, vector<double> &auxstarts, vector<double> &auxends, vector<double> &auxrates){

	for(int i = 0; i<auxstarts.size(); i++){
		if(t>=auxstarts[i] && t<auxends[i]){
			//if(until>auxstarts[i] && until<=auxends[i] || until == numeric_limits<double>::infinity())
				return auxrates[i] - r;
		}

	}
	return 0.0;
	cout<<"problem!!!!!!!!!!!"<<endl;
	
	
}

void GetObservedIntervals(double t_previous, double t0, const vector<double> &starts, const vector<double> &ends, vector<double> &obstarts, vector<double> &obends) {
	if (starts.empty()) {
		obstarts.push_back(t_previous);
		obends.push_back(t0);
		return;
	}
	vector<double> rs;//removed start
	vector<double> re;
	
	for (int i = 0; i < starts.size(); ++i) {
		//cerr<<"starts: "<<starts[i]<<endl;
		//cerr<<"ends: "<<ends[i]<<endl;
		if (ends[i] <= t_previous || starts[i] >= t0) {
			continue;
		}
		if (starts[i] <= t_previous && ends[i] >= t0) {
			return;
		} else if (starts[i] <= t_previous && ends[i] <= t0) {
			rs.push_back(t_previous);
			re.push_back(ends[i]);
		} else if (starts[i] >= t_previous && ends[i] <= t0) {
			rs.push_back(starts[i]);
			re.push_back(ends[i]);
		} else if (starts[i] >= t_previous && ends[i] >= t0) {
			rs.push_back(starts[i]);
			re.push_back(t0);
		}
	}
	
	if (rs.size() == 0) {
		obstarts.push_back(t_previous);
		obends.push_back(t0);
		return;
	}

	if (t_previous < rs[0]) {
		obstarts.push_back(t_previous);
		obends.push_back(rs[0]);
	}

	for (int i = 0; i < rs.size()-1; ++i) {
		obstarts.push_back(re[i]);
		obends.push_back(rs[i+1]);
	}

	if (t0 > re.back()) {
		obstarts.push_back(re.back());
		obends.push_back(t0);
	}	
	return;
}

inline vector<vartrajrange> torange(const vector<ctbn::Trajectory> &data, const ctbn::Context &contexts) {
	vector<vartrajrange> ret;
	if (data.empty()) return ret;
	//int nv = data[0].GetTraj().size();
	for(auto &x : data) for(auto &p : contexts) for(int s=0;s<p.second;s++)
		ret.emplace_back(&x,eventtype(p.first,s));
	return ret;
}

pcim::pcim(const vector<ctbn::Trajectory> &data,
		const vector<shptr<pcimtest>> &tests,
		const pcimparams &params,
		const ctbn::Context &contexts) {
	const vector<vartrajrange> &d = torange(data, contexts);
	ss s = suffstats(d);
	globalm = data.size();
	build(d,s,tests,score(s,params),params);
}

void pcim::learnNewModel(const std::vector<ctbn::Trajectory> &data, const std::vector<shptr<pcimtest>> &tests, const pcimparams &params, const ctbn::Context &contexts) {
	const vector<vartrajrange> &d = torange(data, contexts);
	ss s = suffstats(d);
	globalm = data.size();
	build(d,s,tests,score(s,params),params);
}

pcim::ss pcim::suffstats(const std::vector<vartrajrange> &data) {
	ss ret;
	ret.n=0.0;
	ret.t=0.0;

	for(const auto &x : data) {
		ret.t += x.range.second-x.range.first;
		const ctbn::VarTrajectory &vtr = (*(x.tr)).GetVarTraj(x.event.var);
		auto i0 = vtr.upper_bound(x.range.first);
		auto i1 = vtr.upper_bound(x.range.second);
		//ret.n += distance(i0,i1);		
		for(auto i = i0;i!=i1;++i) {
			if(i->second == x.event.state && i->second != -2 && i->second != -3)
				ret.n++;
		}
	}
	return ret;
}


double pcim::score(const ss &d, const pcimparams &p) {
	double a_n = p.a+d.n;
	double b_n = p.b+d.t;

	return p.lk
		+ lgamma(a_n) - p.lga
		+ p.alb - a_n*log(b_n);
}

void pcim::calcleaf(const ss &d, const pcimparams &p) {
	rate = (p.a+d.n)/(p.b+d.t);
}

pcim::pcim(const vector<vartrajrange> &data, const pcim::ss &s,
		const vector<shptr<pcimtest>> &tests,
		double basescore, const pcimparams &params, int m) {
	globalm = m;
	build(data,s,tests,basescore,params);
}

pcim::testpick pcim::picktest(const vector<vartrajrange> &data,
               const vector<shptr<pcimtest>> &tests,
			const pcimparams &params,
               int procnum, int nproc) const {
	double sc = -numeric_limits<double>::infinity();
	testpick ret;
	ret.testnum = -1; ret.s1 = ret.s2 = sc;
	for(int i=procnum;i<tests.size();i+=nproc) {
		auto &t = tests[i];
		vector<vartrajrange> td1,td2;
		for(auto &x : data) t->chop(x,td1,td2);//chop, passed along the tree
		if (td1.empty() || td2.empty()) continue;
		ss tss1 = suffstats(td1), tss2 = suffstats(td2);
		if (tss1.n<params.mne || tss2.n<params.mne) continue;
		double rs1 = score(tss1,params);
		double rs2 = score(tss2,params);
		if (rs1+rs2>sc) {
			ret.testnum = i;
			ret.s1 = rs1; ret.s2=rs2; sc = rs1+rs2;
			ret.ss1 = move(tss1); ret.ss2 = move(tss2);
			ret.test = t;
			ret.d1 = move(td1);
			ret.d2 = move(td2);
		}
	}
	return ret;
}

void pcim::build(const vector<vartrajrange> &data, const ss &s,
		const vector<shptr<pcimtest>> &tests,
		double basescore, const pcimparams &params) {
	assert(!data.empty());
	vector<vartrajrange> d1,d2;
	test.reset();
	double sc = basescore;
	testpick pick;
	if (params.nproc<=1) {
		pick = picktest(data,tests,params);
		if (pick.s1+pick.s2>sc) sc = pick.s1+pick.s2;
	} else {
		vector<future<testpick>> futs(params.nproc);
		for(int i=0;i<params.nproc;i++)
			futs[i] = async(launch::async,
							&pcim::picktest,this,data,tests,params,i,params.nproc);
		int picki = tests.size();
		for(auto &f : futs) {
			testpick p = f.get();
			double newsc = p.s1+p.s2;
			if (newsc>sc || (newsc==sc && p.testnum<picki)) {
				pick = move(p);
				picki = p.testnum;
				sc = newsc;
			}
		}
	}
	if (sc > basescore) {
		test = pick.test;

		ttree = shptr<pcim>(new pcim(pick.d1,pick.ss1,
									tests,pick.s1,params,globalm));
		ftree = shptr<pcim>(new pcim(pick.d2,pick.ss2,
									tests,pick.s2,params,globalm));

	} else {
		ttree.reset();
		ftree.reset();
	}	
	calcleaf(s,params);
	stats = s;
}

double pcim::getevent(const ctbn::Trajectory &tr, double &t, double expsamp, double unisamp,
		double normsamp, int &var, int &state, double maxt, const ctbn::Context &contexts) const {
	//cout << "-----" << endl;
	//tr.print(cout); cout << endl;
	//cout << t << " w/ " << expsamp << endl;
	double until;
	map<eventtype, const pcim *, eventcomp> leaves;
	double r = getrate(tr,t,until,leaves,contexts); //know which event map to which leaf
	while(expsamp>(until-t)*r) {
		expsamp -= (until-t)*r;
		//cout << r << " until " << until << " (" << expsamp << ")" << endl;
		if (until>maxt) return maxt;
		t = until;
		r = getrate(tr,t,until,leaves,contexts);
	}
	//cout << r << " through " << t+expsamp/r << " [" << until << "]" << endl;
	//var = leaves.size()-1;//not necessary??????
	for(map<eventtype, const pcim *>::iterator it = leaves.begin();it!=leaves.end();it++) {
		unisamp -= it->second->rate/r;
		if (unisamp<=0) { var = it->first.var; state = it->first.state; break; } //var and state of sampled event!!
	}
	return t+expsamp/r;//time of sampled event
}

// Samples auxilary events
double pcim::geteventaux(const ctbn::Trajectory &tr, double &t, double expsamp, double unisamp,
		double normsamp, int &var, double maxt, const ctbn::Context &contexts, vector<double> &auxstarts, vector<double> &auxends, vector<double> &auxrates) const {

	double until;
	map<eventtype, const pcim *, eventcomp> leaves;
	double r = getrate(tr,t,until,leaves,contexts);
	//cerr<<"no1: "<<"r: "<<r<<" t: "<<t<<" until: "<<until<<endl;
	r = getactualrate(r, t, until, auxstarts, auxends, auxrates);
	//cout<<"actual rate1: "<<r<<endl;
	while(expsamp>(until-t)*r) {
		expsamp -= (until-t)*r;
		if (until>maxt) return maxt;
		t = until;
		r = getrate(tr,t,until,leaves,contexts);
		//cerr<<"no2: "<<"r: "<<r<<" t: "<<t<<" until: "<<until<<endl;
		r = getactualrate(r, t, until, auxstarts, auxends, auxrates);
		//cout<<"actual rate2: "<<r<<endl;
	}
	//no need to iterate through leaves, since we only sample the one and only var (not care about state yet)
	return t+expsamp/r;//time of sampled event
}


double pcim::getrate(const ctbn::Trajectory &tr, double t, double &until,
			map<eventtype, const pcim *, eventcomp> &ret, const ctbn::Context &contexts) const {
	until = numeric_limits<double>::infinity();
	double r = 0.0;
	for(int i=0;i<contexts.VarList().size();i++){
		int varid = contexts.VarList()[i]; //fixed 3/17/14 
		for(int s=0; s<contexts.Cardinality(varid); s++){
			r += getratevar(tr,varid,s,t,until,ret[eventtype(varid,s)]);
		}
	}
	return r;
}

//new
double pcim::getauxrates(const ctbn::Trajectory &tr, double &t, int card, double &until, double &r, double varid) const {
	until = numeric_limits<double>::infinity();
	r = 0.0;	
	for(int s = 0; s < card; s++){
		double tmp = getratevaraux(tr,varid,s,t,until);
		//cerr<<"tmp: "<<tmp<<endl;
		r += tmp;
	}
	//cerr<<"r: "<<r<<endl;
	//only sample one variable in the gibbs sampler, same for all states so just multiply - problem for eventtest
	//r = card * getratevaraux(tr,varid,0,t,until);
	return t;
}

double pcim::getratevar(const ctbn::Trajectory &tr, int var, int state, double t, double &until,
			const pcim *&leaf) const {
	if (!test) { leaf = this; return rate; }//reached leaf
	double til;
	bool dir = test->eval(tr,eventtype(var, state),t,til);
	if (til<until) until = til;
	return (dir ? ttree : ftree)->getratevar(tr,var,state,t,until,leaf);
}

// getratevar that considers state.
double pcim::getratevar_state(const ctbn::Trajectory &tr, std::vector<shptr<generic_state> > &jointstate, int varid, eventtype testevent, double t, double &until, const std::vector<int> &testindexes, int index) const {
	if (!test) {return rate;}//reached leaf
	double til;
	//cerr<<"111111111"<<endl;
	//cerr<<"size: "<<jointstate.size()<<endl;
	//cerr<<"index: "<<index<<endl;
	//for(int i = 0; i<jointstate.size(); i++){
		//cerr<<"content!!!!!!: "<<i;
		//jointstate[i]->print();
		//cerr<<endl;
	//}
	bool dir = test->neweval(tr, jointstate[index], varid, testevent,t,til);
	//cerr<<"true or false? "<<dir<<endl;
	if (til<until) until = til;
	return dir ? ttree -> getratevar_state(tr,jointstate, varid, testevent,t,until, testindexes, index + 1) : ftree -> getratevar_state(tr,jointstate, varid, testevent,t,until, testindexes, index+testindexes[index]);
}

// Get aux rate. If the current test result may depend on the current var, take the maximum of both branches.
double pcim::getratevaraux(const ctbn::Trajectory &tr, int varid, int state, double t, double &until) const {
	if (!test) { return rate; }
	double til;
	if(test->getauxv() == -1 || test->getauxv() == varid){ //use maximum!
		double until1 = until;
		double until2 = until;
		double rate1 = ttree -> getratevaraux(tr,varid,state,t,until1);
		double rate2 = ftree -> getratevaraux(tr,varid,state,t,until2);	
		until = until1 < until2? until1 : until2;
		return rate1 > rate2? rate1 : rate2;				
	}
	else{
		bool dir = test->eval(tr,eventtype(varid, state),t,til);
		if (til < until) until = til;
		return (dir ? ttree : ftree)->getratevaraux(tr,varid,state,t,until);			
	}
}

// First argument should be the sampled varid, second argument is event. 
// When calculating the likelihood, if not the sampled var, the
// procedure is the same (in the else statement). If the sampled var, need to calculate likelihood in all
// observed areas between t_previous and t0, no matter evidence or not (rate for evidence handled out of the
// function).
double pcim::Getlikelihood(int varid, eventtype event, ctbn::Trajectory &tr, std::vector<shptr<generic_state> > &jointstate, const std::vector<int> &testindexes, const std::vector<int> &own_var_list, double t_previous, double t0, std::vector<double> &rates, const ctbn::Context &contexts, const std::vector<double> &starts, const std::vector<double> &ends) const{
	//cerr<<"starting...."<<endl;
	rates.reserve(contexts.Cardinality(event.var));
	//printtr(cout,tr,3);
	double P = 0.0;//log
	for(int i = 0; i < own_var_list.size(); i++){
		//cerr<<"t_previous: "<<t_previous<<endl;
		//cerr<<"t0: "<<t0<<endl;
		//cerr<<"current var: "<<own_var_list[i]<<endl;
		//cerr<<"varid: "<<varid<<endl;
		
		// for the sampled var only
		if(varid == own_var_list[i]) {
			//cerr<<"in here?"<<endl;
			//cerr<<"t0: "<<t0<<endl;
			// first get observed intervals between t_previous and t0, then propograte each interval
			// and get likelihood
			vector<double> obstarts;
			vector<double> obends;
			GetObservedIntervals(t_previous, t0, starts, ends, obstarts, obends);
			for (int state = 0; state < contexts.Cardinality(own_var_list[i]); ++state) {
				for (int j = 0; j < obstarts.size(); ++j) {
					// This two mimics t and t0, but only for one observed interval
					double t_small = obstarts[j];
					double t0_small = obends[j];
					while(t_small < t0_small){

						double until = numeric_limits<double>::infinity();

						double temprate = getratevar_state(tr, jointstate, varid, eventtype(own_var_list[i],state), t_small, until, testindexes, 0);
						//cerr<<"rate: "<<temprate<<endl;
						//cerr<<"until: "<<until<<endl;
						if(until < t0_small){
							//cerr<<"until-t: "<<until-t<<endl; 
							P += -1*temprate*(until-t_small);
						}
						else{
							P += -1*temprate*(t0_small-t_small);
						}
						t_small = until;
					}
				}
			}
		}
		//do this for all event. If the sampled event, do not add likelihood (taken care above.)
		for (int state = 0; state < contexts.Cardinality(own_var_list[i]); ++state) {
			double t = t_previous;
			//cerr<<"current state: "<<state<<endl;
			while(t < t0){
				double until = numeric_limits<double>::infinity();
				double temprate = getratevar_state(tr, jointstate, varid, eventtype(own_var_list[i],state), t, until, testindexes, 0);
				//cerr<<"rate: "<<temprate<<endl;
				//cerr<<"until: "<<until<<endl;
				if(until < t0){
					//cerr<<"until-t: "<<until-t<<endl; 
					if (varid != own_var_list[i]) {
						P += -1*temprate*(until-t);
					}
				}
				else{
					if (varid != own_var_list[i]) {
						P += -1*temprate*(t0-t);
					}
					//cerr<<"var: "<<event.var<<" "<<own_var_list[i]<<endl;
					if (event.var == own_var_list[i]) {
						rates.push_back(temprate);
					}
				}
				t = until;
			}
		}
	}
	//cerr<<"rate size: "<<rates.size()<<endl;
	if (rates.empty()) {
		cerr<<"Error, did not get correct rate"<<endl;	
	}
	//cerr<<"finished a loop"<<endl;
	return P;
}


int pcim::Makeindex(vector<int> &indexes, int i) const{
	if(!test) {return 0;}
	int nleft = ttree->Makeindex(indexes, i+1);
	indexes[i] = nleft+1;
	return ftree->Makeindex(indexes, i+indexes[i]) + indexes[i];
}

// Get initilized state vector in order.
void pcim::StateInit(std::vector<shptr<generic_state> > &jointstate) const{
	if(!test) {return;}	
	jointstate.push_back(test->getteststate()->initialize());
	ttree -> StateInit(jointstate);
	ftree -> StateInit(jointstate);	
}

int pcim::counttest() const{
	if(!test) return 0;
	return 1 + ttree -> counttest() + ftree->counttest();
}

void pcim::getnewstates(std::vector<shptr<generic_state> > &jointstate, const std::vector<int> &testindexes, eventtype event, double t0, int index, int varid) const{
	if(!test) {return;}
	jointstate[index] = test->stateupdate(jointstate[index], event, t0, varid);
	ttree -> getnewstates(jointstate, testindexes, event, t0, index+1, varid);
	ftree -> getnewstates(jointstate, testindexes, event, t0, index+testindexes[index], varid);	
}

double pcim::calcDataLikelihood(const vector<ctbn::Trajectory> &data, const ctbn::Context &contexts) {
	double ret = 0;
	const vector<vartrajrange> &d = torange(data, contexts);
	ret = passDownforLL(d);
	return ret;
}

void pcim::updateparams(const std::vector<ctbn::Trajectory> &data, const ctbn::Context &contexts, const pcimparams &params) {
	const vector<vartrajrange> &d = torange(data, contexts);
	passDownforParam(d, params);
	return;	
}

void pcim::passDownforParam(const vector<vartrajrange> &data, const pcimparams &params) {
	if (!test) { //reached leaf
		ss stats = suffstats(data);
		rate = (params.a + stats.n)/(params.b + stats.t);
		return;
	}
	vector<vartrajrange> td1,td2;
	for (auto &x:data) test->chop(x,td1,td2);
	ttree->passDownforParam(td1,params);
	ftree->passDownforParam(td2,params);
	return;
}

double pcim::passDownforLL(const vector<vartrajrange> &data) {
	if (!test) { //reached leaf
		double r = rate;
		ss stats = suffstats(data);
		cerr<<"rate: "<<r<<" n: "<<stats.n<<" t: "<<stats.t<<endl;
		return stats.n*log(r) - r*stats.t;
	}
	vector<vartrajrange> td1,td2;
	for (auto &x:data) test->chop(x,td1,td2);
	return ttree->passDownforLL(td1) + ftree->passDownforLL(td2);
}

void pcim::print(ostream &os) const {
	printhelp(os,0);
}
void pcim::print(ostream &os, const datainfo &info) const {
	printhelp(os,0,&info);
}

void pcim::todot(ostream &os, const datainfo &info) const {
	os << "digraph {" << endl;
	int nn = 0;
	todothelp(os,-1,false,nn,info);
	os << "}" << endl;
}

void pcim::todothelp(ostream &os, int par, bool istrue, int &nn, const datainfo &info) const {
	int mynode = nn++;
	os << "\tNODE" << mynode << " [label=\"";
	if (!ttree) os << "rate = " << rate;
	else test->print(os,info);
	os << "\"];" << endl;
	if (par>=0) os << "\tNODE" << par << " -> NODE" << mynode << " [label=\"" <<
			(istrue ? "Y: " : "N: ") << stats.n << ',' << stats.t << "\"];" << endl;
	if (ttree) {
		ttree->todothelp(os,mynode,true,nn,info);
		ftree->todothelp(os,mynode,false,nn,info);
	}
}

void pcim::printhelp(ostream &os, int lvl, const datainfo *info) const {
	for(int i=0;i<lvl;i++) os << "   ";
	if (!ttree)
		os << "rate = " << rate << "[" << stats.n << ',' << stats.t << "]" << endl;
	else {
		os << "if ";
		if (info!=nullptr) test->print(os,*info);
		else test->print(os);
		os << " [" << stats.n << ',' << stats.t << "]" << endl;
		ttree->printhelp(os,lvl+1,info);
		for(int i=0;i<lvl;i++) os << "   ";
		os << "else" << endl;
		ftree->printhelp(os,lvl+1,info);
	}
}

BOOST_CLASS_EXPORT_IMPLEMENT(pcimtest)
BOOST_CLASS_EXPORT_IMPLEMENT(lasttest)
BOOST_CLASS_EXPORT_IMPLEMENT(timetest)
BOOST_CLASS_EXPORT_IMPLEMENT(varcounttest)
BOOST_CLASS_EXPORT_IMPLEMENT(varstattest<varcounttest>)
BOOST_CLASS_EXPORT_IMPLEMENT(eventcounttest)
BOOST_CLASS_EXPORT_IMPLEMENT(varstattest<eventcounttest>)
BOOST_CLASS_EXPORT_IMPLEMENT(vartest)
BOOST_CLASS_EXPORT_IMPLEMENT(eventtest)
BOOST_CLASS_EXPORT_IMPLEMENT(pcim)
