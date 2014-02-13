#include "traj.h"
#include <sstream>
#include <iomanip>

using namespace std;

void printtime(ostream &os, double t, int colw, bool ishrs) {
	if (!ishrs) os << setw(colw) << t;
	else {
		int sec = floor(t*60*60);
		int min = sec/60;
		sec %= 60;
		int hr = min/60;
		min %= 60;
		int day = hr/24;
		hr %= 24;
		stringstream ss;
		ss << setfill('0');
		ss << day << "d " << setw(2) << hr << ':' << setw(2) << min;
		for(int i=ss.str().size();i<colw;i++) os << ' ';
		os << ss.str();
	}
}

void printdynamic(ostream &os, const Trajectory &tr, bool incolumns, bool ishrs) {
	vector<decltype(tr.find(0)->second.begin())> it;
	vector<int> varid;
	int ndone = 0;
	for(int i=0;i<tr.size();i++) {
		if (tr.find(i)->second.empty()) continue;
		it.push_back(tr.find(i)->second.begin());
		varid.push_back(i);
	}
	int colw = 10;
	os << "events:" << endl;
	auto osf = os.setf(ios_base::fixed);
	if (incolumns) {
		os << "time      ";
		for(int i=0;i<varid.size();i++)
			os << "|" << setw(colw) << varid[i];
		os << endl;
	}
	while(ndone<it.size()) {
		int i=-1;
		double t = numeric_limits<double>::infinity();
		for(int j=0;j<it.size();j++)
			if (it[j]!=tr.find(varid[j])->second.end() && it[j]->first<t) {
				t = it[j]->first;
				i = j;
			}
		if (incolumns) {
			printtime(os,t,colw,ishrs);
			for(int j=0;j<i;j++) os << "|          ";
			os << "|" << setw(colw) << it[i]->second;//state
			for(int j=i+1;j<varid.size();j++) os << "|          ";
			os << endl;
		} else {
			os << setw(4) << i << " @ ";
			printtime(os,t,colw,ishrs);
			os << " " << setw(colw) << it[i]->second << endl;
		}
		if (++it[i] == tr.find(varid[i])->second.end()) ndone++;
	}
	os.setf(osf);
}

void printtr(ostream &os, const Trajectory &tr, bool incolumns, bool ishrs) {
	printdynamic(os,tr,incolumns,ishrs);
}

