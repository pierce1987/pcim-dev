#include "pcim.h"
#include "akipcim.h"
#include <vector>
#include "oldmethod.h"

using namespace std;
using namespace boost::archive;

const int maxn = 100000000;
const double predthresh = 1.5;
const string akidir{"/data/aki"};
const int nproc = 8;
//const double predtime = 48.0; // in hours
//const double mintime = 48.0+predtime;
const double predtime = 24.0; // in hours
const double mintime = 24.0+predtime;
const string predname{"futurecreatinine"};
const string basevarname{"baselinecreatinine"};
const string dynname{"Creatinine"};

pair<pcim,datainfo> loadmodel(string filename) {
     ifstream fs(filename.c_str());
     xml_iarchive ar(fs);
     datainfo info;
     pcim model;
     ar >> BOOST_SERIALIZATION_NVP(info) >> BOOST_SERIALIZATION_NVP(model);
     return make_pair(model,info);
}

void savemodel(string filename, const pcim &model, const datainfo &info) {
     ofstream fs(filename.c_str());
     xml_oarchive ar(fs);
     ar << BOOST_SERIALIZATION_NVP(info) << BOOST_SERIALIZATION_NVP(model);
}

template<typename X>
int get(const std::map<X,int> &m, const X &x) {
	auto l = m.find(x);
	if (l==m.end()) return -1;
	return l->second;
}

template<typename Y>
shptr<Y> get(const std::vector<shptr<Y>> &m, const int &x) {
	if (x<0 || x>=m.size()) return shptr<Y>{};
	return m[x];
}

void topred(const dataset &base, dataset &pred) {
	int outid = get(pred.info.svarid,predname);
	int ageid = get(pred.info.svarid,string("startage"));
	int sexid = get(pred.info.svarid,string("sex"));
	int baseid = get(base.info.svarid,basevarname);
	int dynid = get(base.info.dvarid,dynname);
	shptr<normmethod> dynnorm = get(base.info.dnorm,dynid);
	assert(ageid>=0 && sexid>=0 && dynnorm);
	// don't go past when one of these is recorded
	set<int> invalidval{get(pred.info.dvarid,string("Peritoneal Dialysis")),
				get(pred.info.dvarid,string("Hemofiltration"))};
	int i=0;
	for(auto &tr : base.ds) {
		double baseval = baseid>=0 ? tr.sx[baseid] : 1.0;
		double tlimit = numeric_limits<double>::infinity();
		for(auto &i : invalidval) 
			for(auto &x : tr[i])
				if (x.first < tlimit) tlimit = x.first;
		for(auto &x : tr[dynid]) {
			if (x.first >= tlimit) break;
			bool predval = (dynnorm->unnormalize
					(x.second,tr.sx[ageid],
							tr.sx[sexid])/baseval)>=predthresh;
			if (x.first >= mintime) {
				pred.ds.emplace_back(tr,x.first-predtime);
				if (pred.ds.back().sx.size() <= outid)
					pred.ds.back().sx.resize(outid+1);
				pred.ds.back().sx[outid] = predval;
			}
			// don't predict past first positive
			if (predval) break;
		}
		i++;
	}
}

void predictlast(const dataset &ds, ostream &os) {
	int outid = get(ds.info.svarid,predname);
	int ageid = get(ds.info.svarid,string("startage"));
	int sexid = get(ds.info.svarid,string("sex"));
	int baseid = get(ds.info.svarid,basevarname);
	int dynid = get(ds.info.dvarid,dynname);
	shptr<normmethod> dynnorm = get(ds.info.dnorm,dynid);
	assert(ageid>=0 && sexid>=0 && dynnorm);
	for(auto &tr : ds.ds) {
		double baseval = baseid>=0 ? tr.sx[baseid] : 1.0;
		double age = ageid>=0 ? tr.sx[ageid] : 0.0;
		double sex = sexid>=0 ? tr.sx[sexid] : 0.0;
		os << tr.sx[outid] << ' ';
		if (tr[dynid].empty()) os << 1.0 << endl;
		else
			os << dynnorm->unnormalize(
				tr[dynid].rbegin()->second,age,sex)
					/baseval << endl;
	}
}

template<typename R>
void predict(const dataset &ds, const pcim &model, ostream &os, R &rand) {
	int outid = get(ds.info.svarid,predname);
	int ageid = get(ds.info.svarid,string("startage"));
	int sexid = get(ds.info.svarid,string("sex"));
	int baseid = get(ds.info.svarid,basevarname);
	int dynid = get(ds.info.dvarid,dynname);
	shptr<normmethod> dynnorm = get(ds.info.dnorm,dynid);
	assert(ageid>=0 && sexid>=0 && dynnorm);
	for(auto &tr : ds.ds) {
		double baseval = baseid>=0 ? tr.sx[baseid] : 1.0;
		double age = ageid>=0 ? tr.sx[ageid] : 0.0;
		double sex = sexid>=0 ? tr.sx[sexid] : 0.0;
		os << tr.sx[outid] << ' ';
/*
		auto chk = [baseval,dynnorm,age,sex](double v) {
					return dynnorm->unnormalize(v,age,sex)/baseval
								> predthresh; };
		os << model.eventprob(tr,tr[0].endtime()+predtime,
				dynid,chk,rand,nproc) << "   ";
*/
		double thresh = dynnorm->normalize(baseval*predthresh,age,sex);
		os << model.eventprobthresh(tr,tr[0].endtime()+predtime,
				dynid,thresh,rand,nproc) << endl;
	}
}

void dumpfeat(const dataset &data, const pcim &model, const string &fn) {
	ofstream fs(fn.c_str());
	int outid = data.info.svarid.at(predname);

	fs << predname;
	for(auto &s : data.info.svarnames)
		if (s!=predname) fs << ',' << s;
	auto fnames = model.featurenames();
	for(auto &s : fnames)
		fs << ',' << s;
	fs << endl;
		
	for(auto &tr : data.ds) {
		auto f = model.trajtofeatures(tr);
		fs << tr.sx[outid];
		for(int i=0;i<tr.sx.size();i++)
			if (i!=outid) fs << ',' << tr.sx[i];
		for(auto &x : f) fs << ',' << x;
		fs << endl;
	}
}

int main(int argc, char **argv) {
	int n = min(akicount(akidir),maxn);
	int ntrain = n/2.0;
	if (argc>1) {
		string ntstr(argv[1]);
		if (ntstr.back()=='%') {
			double frac = stof(ntstr)/100.0;
			ntrain = (int)floor(n*frac);
		} else ntrain = stoi(ntstr);
		if (ntrain>n) ntrain = n;
	}

	dataset train,test;
	pcim model;

	if (argc>2) {
		tie(model,train.info) = loadmodel(argv[2]);
		akiloadgiveninfo(akidir,train,ntrain,0);
	} else {
		akiload(akidir,train,ntrain,0);
		auto tests = akitests(train.info);
		auto params = akipcimparams(tests.size(),nproc);
		model = pcim{train.ds,tests,params};
		savemodel("akimodel.pcim",model,train.info);
	}
	test.info = train.info;
	akiloadgiveninfo(akidir,test,n-ntrain,ntrain);

	datainfo info=train.info;
	info.svarid[predname] = info.svarnames.size();
	info.snorm.emplace_back(new normmethod());
	info.svarnames.emplace_back(predname);
	info.sisdiscrete.push_back(false);

	dataset trainpred, testpred;
	trainpred.info = info;
	testpred.info = info;
	topred(train,trainpred);
	topred(test,testpred);

	ofstream trainpredlast("akitrainlast.pred");
	predictlast(trainpred,trainpredlast);
	ofstream testpredlast("akitestlast.pred");
	predictlast(testpred,testpredlast);

	dumpfeat(trainpred,model,"akitrain.csv");
	dumpfeat(testpred,model,"akitest.csv");
	//savemodel("akimodel2.pcim",model,train.info);
	//exit(1);

	dumpoldfeat(trainpred,"akioldtrain.csv",predname,mintime-predtime);
	dumpoldfeat(testpred,"akioldtest.csv",predname,mintime-predtime);

	random_device rd;
	int seed = rd();
	mt19937 randgen(seed);
	ofstream trainpredof("akitrain.pred");
	predict(trainpred,model,trainpredof,randgen);
	ofstream testpredof("akitest.pred");
	predict(testpred,model,testpredof,randgen);
}
	
