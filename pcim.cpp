#include "pcim.h"
#include <algorithm>
#include <cmath>
#include <mutex>
#include <future>
#include "serial.h"

using namespace std;

inline vector<vartrajrange> torange(const vector<traj> &data) {
	vector<vartrajrange> ret;
	if (data.empty()) return ret;
	int nv = data[0].size();
	for(auto &x : data) for(int v=0;v<nv;v++)
		ret.emplace_back(&x,v);
	return ret;
}

pcim::pcim(const vector<traj> &data,
		const vector<shptr<pcimtest>> &tests,
		const pcimparams &params) {
	const vector<vartrajrange> &d = torange(data);
	ss s = suffstats(d);
	globalm = data.size();
	build(d,s,tests,score(s,params),params);
}


pcim::ss pcim::suffstats(const std::vector<vartrajrange> &data) {
	ss ret;
	ret.n=0.0;
	ret.t=0.0;
	ret.sum = zerowt;
	ret.invar = zerowtvar;
	ret.sum2 = 0.0;
	for(const auto &x : data) {
		ret.t += x.range.second-x.range.first;
		const vartraj &vtr = (*(x.tr))[x.eventtype];
		auto i0 = vtr.upper_bound(x.range.first);
		auto i1 = vtr.upper_bound(x.range.second);
#ifdef USEPERSIST
		double lasty = 0.0;
		if (i0!=vtr.begin() && !vtr.empty()) lasty = std::prev(i0)->second.v;
#endif
		//ret.n += distance(i0,i1);
		for(auto i = i0;i!=i1;++i) {
			double y = i->second.v;
#ifdef USEPERSIST
			wtT x = predfeat(lasty);
			ret.sum += x*y;
			ret.invar += x*x.transpose();
			lasty = y;
#else
			ret.sum += y;
			ret.invar += 1.0;
#endif
			ret.sum2 += y*y;
			ret.n++;
		}
	}
	//ret.sum2 = ret.var/ret.n - ret.sum*ret.sum/ret.n/ret.n; // maybe not numerically
				// stable enough?
	return ret;
}

constexpr double log2pi() { return std::log(8.0*std::atan(1)); }


double pcim::score(const ss &d, const pcimparams &p) {
	double a_n = p.a+d.n;
	double b_n = p.b+d.t;

	double hn = d.n/2.0;
	double va_n = p.va+hn;
	wtvarT vk_n = p.vk + d.invar;
	wtT    mu_n = div(p.vk*p.vm + d.sum , vk_n);
	//double meandiff = d.mean-p.vm;
	//double vb_n = p.vb + hn*d.var + hn*p.vk*meandiff*meandiff/vk_n;
	double vb_n = p.vb + (p.m2k - quad(mu_n,vk_n) + d.sum2)/2.0;

	return p.lk

		+ lgamma(a_n) - p.lga
		+ p.alb - a_n*log(b_n)

		+ lgamma(va_n) - p.lgva
// check if below should be reversed! (I think it's okay, but not sure)
		+ p.valvb - va_n*log(vb_n)
#ifdef USEPERSIST
		+ p.lvk/2.0 - log(vk_n.determinant())/2.0
#else
		+ p.lvk/2.0 - log(vk_n)/2.0
#endif
		- hn*log2pi();
}

void pcim::calcxxinvsqrt(const ss &d) {
#ifdef USEPERSIST
	wtvarT ivtemp = d.invar;
	for(int i=0;i<npredfeat;i++) ivtemp(i,i) += invarreg;
	Eigen::SelfAdjointEigenSolver<wtvarT> es(ivtemp/(d.n+invarreg));
	xxinvsqrt = es.operatorInverseSqrt();
#endif
}

void pcim::calcleaf(const ss &d, const pcimparams &p) {
	rate = (p.a+d.n)/(p.b+d.t);
	mu = div(p.vk*p.vm + d.sum, p.vk+d.invar);
	double hn = d.n/2.0;
	double va_n = p.va+hn;
	wtvarT vk_n = p.vk + d.invar;
	//double meandiff = d.mean-p.vm;
	//double vb_n = p.vb + hn*d.var + hn*p.vk*meandiff*meandiff/vk_n;
	double vb_n = p.vb + (p.m2k - quad(mu,vk_n) + d.sum2)/2.0;
	// below is "Bayes equivalent"
	//sigma = sqrt(vb_n*(vk_n+1)/(vk_n*(va_n-1)));
	// instead, use MAP:
	sigma = sqrt(vb_n/(va_n-0.5));
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

double pcim::getevent(const traj &tr, double &t, double expsamp, double unisamp,
		double normsamp, int &var, double &val, double maxt) const {
	//cout << "-----" << endl;
	//tr.print(cout); cout << endl;
	//cout << t << " w/ " << expsamp << endl;
	double until;
	vector<const pcim *> leaves;
	double r = getrate(tr,t,until,leaves);
	while(expsamp>(until-t)*r) {
		expsamp -= (until-t)*r;
		//cout << r << " until " << until << " (" << expsamp << ")" << endl;
		if (until>maxt) return maxt;
		t = until;
		r = getrate(tr,t,until,leaves);
	}
	//cout << r << " through " << t+expsamp/r << " [" << until << "]" << endl;
	var = leaves.size()-1;
	for(int i=0;i<leaves.size();i++) {
		unisamp -= leaves[i]->rate/r;
		if (unisamp<=0) { var = i; break; } //var of sampled event!!
	}
#ifdef USEPERSIST
	double x = 0.0;
	if (!tr[var].empty()) x = tr[var].rbegin()->second.v;
	val = predfeat(x).dot(leaves[var]->mu);
#else
	val = leaves[var]->mu;
#endif
	val += normsamp*leaves[var]->sigma;//value of sampled event!!
	return t+expsamp/r;//time of sampled event!!
}

double pcim::getrate(const traj &tr, double t, double &until,
			vector<const pcim *> &ret) const {
	until = numeric_limits<double>::infinity();
	ret.resize(tr.size());
	double r = 0.0;
	for(int i=0;i<tr.size();i++)
		r += getratevar(tr,i,t,until,ret[i]);
	return r;
}

double pcim::getratevar(const traj &tr, int var, double t, double &until,
			const pcim *&leaf) const {
	if (!test) { leaf = this; return rate; }//reached leaf
	double til;
	bool dir = test->eval(tr,var,t,til);
	if (til<until) until = til;
	return (dir ? ttree : ftree)->getratevar(tr,var,t,until,leaf);
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
	if (!ttree) os << "rate = " << rate << " (" << mu << ',' << sigma << ")";
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
		os << "rate = " << rate << " (" << mu << ',' << sigma << ") [" << stats.n << ',' << stats.t << ',' << stats.sum << ',' << stats.sum2 << "]" << endl;
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
#ifdef USEPERSIST
	f[0] = d.n - rate*d.t;
	wtT dff = xxinvsqrt*(d.sum - d.invar*mu);
	for(int i=0;i<npredfeat;i++)
		f[i+1] = dff(i)/sigma;
	f[npredfeat+1] = (d.sum2 - 2*mu.dot(d.sum) + mu.transpose()*d.invar*mu)
					/ (sigma*sigma) - d.n;
	bool problem = false;
	for(auto &x : f) if (std::isnan(x)) problem=true;
	if (problem) {
		cout << "problem:" << endl;
		for(auto &x : f) cout << x << ' ';
		cout << endl;
		cout << d.n << ' ' << d.t << ' ' << d.sum2 << ' ' << d.sum << ' ' << d.invar << endl;
		cout << rate << ' ' << sigma << ' ' << mu << endl;
		exit(1);
	}
#else    //leaf features
	f[0] = d.n - rate*d.t;
	f[1] = (d.sum - mu*d.n)/sigma;
	f[2] = (d.sum2 - 2*mu*d.sum + mu*mu*d.n)/(sigma*sigma) - d.n;
#endif
/*
	double N = 0.0, D = 0.0, E = 0.0, T = 0.0;
	for(auto &x : tr) {
		T += x.range.second-x.range.first;
		const vartraj &vtr = (*(x.tr))[x.eventtype];
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
#ifdef USEPERSIST
	array<string,nleaffeat> ret;
	ret[0] = "rate";
	for(int i=0;i<npredfeat;i++)
		ret[i+1] = string("mean")+to_string(i);
	ret[npredfeat+1] = "var";
	return ret;
#else
	return {"rate","mean","var"};
#endif
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

#ifdef USEPERSIST
const pcim::wtT pcim::zerowt = pcim::wtT::Zero(pcim::npredfeat,1);
const pcim::wtvarT pcim::zerowtvar = pcim::wtvarT::Zero(pcim::npredfeat,pcim::npredfeat);
#endif

BOOST_CLASS_EXPORT_IMPLEMENT(pcimtest)
BOOST_CLASS_EXPORT_IMPLEMENT(timetest)
BOOST_CLASS_EXPORT_IMPLEMENT(counttest)
BOOST_CLASS_EXPORT_IMPLEMENT(varstattest<counttest>)
BOOST_CLASS_EXPORT_IMPLEMENT(meantest)
BOOST_CLASS_EXPORT_IMPLEMENT(varstattest<meantest>)
BOOST_CLASS_EXPORT_IMPLEMENT(vartest)
BOOST_CLASS_EXPORT_IMPLEMENT(staticgreqtest)
BOOST_CLASS_EXPORT_IMPLEMENT(staticeqtest)
BOOST_CLASS_EXPORT_IMPLEMENT(pcim)
