#include "pcim.h"
#include <algorithm>
#include <cmath>
#include <mutex>
#include <future>
#include "serial.h"

using namespace std;

inline vector<vartrajrange> torange(const vector<Trajectory> &data, const vector<int> states) {
	vector<vartrajrange> ret;
	if (data.empty()) return ret;
	int nv = data[0].size();
	for(auto &x : data) for(int v=0;v<nv;v++) for(int s=0;s<states[v];s++)
		ret.emplace_back(&x,eventtype(v,s));
	return ret;
}

pcim::pcim(const vector<Trajectory> &data,
		const vector<shptr<pcimtest>> &tests,
		const pcimparams &params,
		const vector<int> &states) {
	const vector<vartrajrange> &d = torange(data, states);
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
		const vartraj &vtr = (*(x.tr)).find(x.event.var)->second;
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

constexpr double log2pi() { return std::log(8.0*std::atan(1)); }


double pcim::score(const ss &d, const pcimparams &p) {
	double a_n = p.a+d.n;
	double b_n = p.b+d.t;
	//double hn = d.n/2.0;

	return p.lk

		+ lgamma(a_n) - p.lga
		+ p.alb - a_n*log(b_n);

		//- hn*log2pi();//?
}

void pcim::calcxxinvsqrt(const ss &d) {
}

void pcim::calcleaf(const ss &d, const pcimparams &p) {
	rate = (p.a+d.n)/(p.b+d.t);
	calcxxinvsqrt(d);
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

//vector<shptr<pcimtest>> ancestors;

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
/*
		//assert(pick.ss1.n+pick.ss2.n==s.n);
		if (pick.ss1.n+pick.ss2.n != s.n) {
			test->print(cout); cout << " " << pick.ss1.n << ' ' << pick.ss2.n << endl;
			bool dbool = false;
			for(auto & a : ancestors) dbool |= (a==test);
			if (dbool) {
				cout << "here!" << endl;
				vector<vartrajrange> td1,td2;
				for(auto &x : data) test->chop(x,td1,td2);
				cout << "done!" << endl;
			} else {
				for(auto &x : data) {
					vector<vartrajrange> td1,td2;
					vector<vartrajrange> xx {x};
					test->chop(x,td1,td2);
					ss tss1 = suffstats(td1), tss2 = suffstats(td2);
					ss press = suffstats(xx);
					if (press.n != tss1.n + tss2.n) {
						cout << "ah ha!" << endl;
						td1.clear(); td2.clear();
						test->chop(x,td1,td2);
						assert(0);
					} else cout << "nope" << endl;
				}
			}
			assert(0);
		}
		ancestors.push_back(test);
*/
		ttree = shptr<pcim>(new pcim(pick.d1,pick.ss1,
									tests,pick.s1,params,globalm));
		ftree = shptr<pcim>(new pcim(pick.d2,pick.ss2,
									tests,pick.s2,params,globalm));
		//ancestors.pop_back();
	} else {
		ttree.reset();
		ftree.reset();
	}
	calcleaf(s,params);
	stats = s;
}

double pcim::getevent(const Trajectory &tr, double &t, double expsamp, double unisamp,
		double normsamp, int &var, int &state, double maxt, const vector<int> &states) const {
	//cout << "-----" << endl;
	//tr.print(cout); cout << endl;
	//cout << t << " w/ " << expsamp << endl;
	double until;
	map<eventtype, const pcim *, comparator> leaves;
	double r = getrate(tr,t,until,leaves,states);
	while(expsamp>(until-t)*r) {
		expsamp -= (until-t)*r;
		//cout << r << " until " << until << " (" << expsamp << ")" << endl;
		if (until>maxt) return maxt;
		t = until;
		r = getrate(tr,t,until,leaves,states);
	}
	//cout << r << " through " << t+expsamp/r << " [" << until << "]" << endl;
	var = leaves.size()-1;//?
	for(map<eventtype, const pcim *>::iterator it = leaves.begin();it!=leaves.end();it++) {
		unisamp -= it->second->rate/r;
		if (unisamp<=0) { var = it->first.var; state = it->first.state; break; } //var and state of sampled event!!
	}
	return t+expsamp/r;//time of sampled event!!
}

double pcim::getrate(const Trajectory &tr, double t, double &until,
			map<eventtype, const pcim *, comparator> &ret, const vector<int> &states) const {
	until = numeric_limits<double>::infinity();
	double r = 0.0;
	for(int i=0;i<tr.size();i++)
		for(int s=0; s<states[i]; s++){
			r += getratevar(tr,i,s,t,until,ret[eventtype(i,s)]);
		}
	return r;
}

double pcim::getratevar(const Trajectory &tr, int var, int state, double t, double &until,
			const pcim *&leaf) const {
	if (!test) { leaf = this; return rate; }//reached leaf
	double til;
	bool dir = test->eval(tr,eventtype(var, state),t,til);
	if (til<until) until = til;
	return (dir ? ttree : ftree)->getratevar(tr,var,state,t,until,leaf);
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

void pcim::getleaffeature(const vector<vartrajrange> &tr, array<double,nleaffeat> &f) const {
	ss d = suffstats(tr);
//leaf features
	f[0] = d.n - rate*d.t;

/*
	double N = 0.0, D = 0.0, E = 0.0, T = 0.0;
	for(auto &x : tr) {
		T += x.range.second-x.range.first;
		const vartraj &vtr = (*(x.tr))[x.var];
		auto i0 = vtr.upper_bound(x.range.first);
		auto i1 = vtr.upper_bound(x.range.second);
		for(auto i = i0;i!=i1;++i) {
			D += i->second.v - mu;
			E += (i->second.v - mu)*(i->second.v - mu);
			N++;
		}
	}
	f[0] = N-rate*T;
	f[1] = D/sigma;
	f[2] = E/(sigma*sigma) - N;
*/
}

void pcim::featurenames(vector<string> &ret, string prefix) const {
	if (!ttree) {
		array<string,nleaffeat> ln = getleaffeaturenames();
		for(auto &n : ln) ret.push_back(prefix+n);
	} else {
		ttree->featurenames(ret,prefix+"T");
		ftree->featurenames(ret,prefix+"F");
	}
}

array<string,pcim::nleaffeat> pcim::getleaffeaturenames() const {

	return {"rate"};

}



void pcim::trajtofeatures(const vector<vartrajrange> &tr,
				vector<double> &f) const {
	constexpr double sqrt2 = std::sqrt(2.0);
	if (!ttree) {
		double w = sqrt(globalm/stats.n);
		array<double,nleaffeat> lf;
		getleaffeature(tr,lf);
		//std::cout << "features for leaf (" << w << ',' << globalm << ',' << stats.n << "): " << lf[0] << ',' << lf[1] << ',' << lf[2] << std::endl;
		for(int i=0;i<nleaffeat-1;i++)
			f.push_back(w*lf[i]);
		f.push_back(w/sqrt2*lf[nleaffeat-1]);
	} else {
		vector<vartrajrange> ttr,ftr;
		for(auto &x : tr) test->chop(x,ttr,ftr);
		ttree->trajtofeatures(ttr,f);
		ftree->trajtofeatures(ftr,f);
	}
}

double pcim::similarity(const vector<vartrajrange> &tr1,
			const vector<vartrajrange> &tr2) const {
	if (!ttree) {
		double w = globalm/stats.n;
		array<double,nleaffeat> lf1,lf2;
		getleaffeature(tr1,lf1);
		getleaffeature(tr2,lf2);
		double ret = 0.0;
		for(int i=0;i<nleaffeat-1;i++)
			ret += lf1[i]*lf2[i];
		ret += (lf1[nleaffeat-1]*lf2[nleaffeat-1])/2.0;
		return ret*w;
	} else {
		vector<vartrajrange> ttr1,ftr1;
		for(auto &x : tr1) test->chop(x,ttr1,ftr1);
		vector<vartrajrange> ttr2,ftr2;
		for(auto &x : tr2) test->chop(x,ttr2,ftr2);
		return ttree->similarity(ttr1,ttr2) + ftree->similarity(ftr1,ftr2);
	}
}


BOOST_CLASS_EXPORT_IMPLEMENT(pcimtest)
BOOST_CLASS_EXPORT_IMPLEMENT(lasttest)
BOOST_CLASS_EXPORT_IMPLEMENT(timetest)
BOOST_CLASS_EXPORT_IMPLEMENT(counttest)
BOOST_CLASS_EXPORT_IMPLEMENT(varstattest<counttest>)
BOOST_CLASS_EXPORT_IMPLEMENT(counteventtest)
BOOST_CLASS_EXPORT_IMPLEMENT(eventstattest<counteventtest>)
BOOST_CLASS_EXPORT_IMPLEMENT(vartest)
BOOST_CLASS_EXPORT_IMPLEMENT(staticgreqtest)
BOOST_CLASS_EXPORT_IMPLEMENT(staticeqtest)
BOOST_CLASS_EXPORT_IMPLEMENT(pcim)
