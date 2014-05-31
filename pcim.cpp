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
			//cout<<"error!!!!!!!!!!!!!!!!"<<endl;
		}

	}
	cout<<"problem!!!!!!!!!!!"<<endl;
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
			if(i->second == x.event.state)
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
	return t+expsamp/r;//time of sampled event!!
}

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
	return t+expsamp/r;//time of sampled event!!
}

//new

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
	//cerr<<"begin: "<<"card: "<<card<<" t: "<<t<<" until: "<<until<<" r: "<<r<<" id:"<<varid<<endl;
	until = numeric_limits<double>::infinity();
	r = 0.0;
	//only sample one variable in the gibbs sampler
	for(int s=0; s<card; s++){
		r += getratevaraux(tr,varid,s,t,until);
	}
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

double pcim::getratevaraux(const ctbn::Trajectory &tr, int varid, int state, double t, double &until) const {
	if (!test) { return rate; }
	double til;
	if(test->getauxv() == -1 || test->getauxv() == varid){ //use maximum!
		//test->eval(tr,eventtype(varid, state),t,til);
		//cout<<"til1: "<<til<<endl;
		//if (til<until) until = til;
		double until1 = until;
		double until2 = until;
		double rate1 = ttree -> getratevaraux(tr,varid,state,t,until1);
		double rate2 = ftree -> getratevaraux(tr,varid,state,t,until2);	
		until = until1<until2? until1 : until2;
		//cout<<"rate1: "<<rate1<<" rate2: "<<rate2<<endl;
		return rate1 > rate2? rate1 : rate2;				
	}
	else{
		bool dir = test->eval(tr,eventtype(varid, state),t,til);
		//cout<<"til2: "<<til<<endl;
		if (til<until) until = til;
		return (dir ? ttree : ftree)->getratevaraux(tr,varid,state,t,until);			
	}
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
