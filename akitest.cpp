#include "pcim.h"
#include "akipcim.h"
#include <vector>
//#include "oldmethod.h"

using namespace std;
using namespace boost::archive;

const double predthresh = 1.5;
const string akidir{"/data/aki"};
const double predtime = 24.0; // in hours
const double mintime = 24.0+predtime;
const string predname{"futurecreatinine"};
const string basevarname{"baselinecreatinine"};
const string dynname{"Creatinine"};


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

map<int,int> getepnummap(void) {
	ifstream f("/data/aki/eplist.txt");
	string l;
	map<int,int> ret;
	while(!f.eof()) {
		getline(f,l);
		if (f.eof()) return ret;
		auto pi = l.find(",");
		if (pi!=string::npos && l[0]>='0' && l[0]<='9') {
			int epnum = stoi(l.substr(0,pi));
			int ind = stoi(l.substr(l.size()-8,4));
			ret[ind] = epnum;
		}
	}
	return ret;
}


int main(int argc, char **argv) {
	dataset ds;
	akiload(akidir,ds,1000000,0);

	int ageid = get(ds.info.svarid,string("startage"));
	int sexid = get(ds.info.svarid,string("sex"));
	int baseid = get(ds.info.svarid,basevarname);
	int dynid = get(ds.info.dvarid,dynname);
	shptr<normmethod> dynnorm = get(ds.info.dnorm,dynid);
	assert(ageid>=0 && sexid>=0 && dynnorm);

	auto epnum = getepnummap();

	double big = 100000000;
	int n0=0,n1=0,n2=0,n3=0;
	ofstream f1("list1.csv"), f2("list2.csv"), f3("list3.csv"), f4("list4.csv");
	ofstream fnn("listnn.csv"), fnr("listnr.csv"), fni("listni.csv"), fnf("listnf.csv");
	ofstream frr("listrr.csv"), fri("listri.csv"), frf("listrf.csv");
	ofstream fii("listii.csv"), fif("listif.csv");
	ofstream fff("listff.csv"), ffp("listfp.csv");
	
	int i=0;
	for(auto &tr : ds.ds) {
		double baseval = baseid>=0 ? tr.sx[baseid] : 1.0;
		double age = ageid>=0 ? tr.sx[ageid] : 0.0;
		double sex = sexid>=0 ? tr.sx[sexid] : 0.0;
		double tht = big;
		double maxval = 0,thval = 0, maxt=0, fval, ft=big;
		for(auto &pt : tr[dynid]) {
			double val = dynnorm->unnormalize(pt.second,age,sex)/baseval+1e-6;
			if (val >=predthresh && tht==big) {
				tht = pt.first;
				thval = val;
			}
			if (val>maxval) {
				maxt = pt.first;
				maxval = val;
			}
			if (ft==big) {
				fval = val;
				ft = pt.first;
			} 
		}
		stringstream ss;
		ss << i << ' ' << epnum[i] << ' ' << ft << ' ' << fval << ' ' << tht << ' ' << thval << ' ' << maxt << ' ' << maxval << ' ' << baseval << endl;
		if (tht<predtime) {
			n0++;
			f1 << ss.str();
		} else if (tht<mintime) {
			n1++;
			f2 << ss.str();
		} else if (tht<big) {
			n2++;
			f3 << ss.str();
		} else {
			n3++;
			f4 << ss.str();
		}
		if (fval<1.5) {
			if (maxval<1.5) fnn << ss.str();
			else if (maxval<2.0) fnr << ss.str();
			else if (maxval<3.0) fni << ss.str();
			else fnf << ss.str();
		} else if (fval<2.0) {
			if (maxval<2.0) frr << ss.str();
			else if (maxval<3.0) fri << ss.str();
			else frf << ss.str();
		} else if (fval<3.0) {
			if (maxval<3.0) fii << ss.str();
			else fif << ss.str();
		} else {
			if (maxval<fval*1.5) fff << ss.str();
			else ffp << ss.str();
		}
		i++;
	}
	cout << n0 << ',' << n1 << ',' << n2 << ',' << n3 << endl;
}
