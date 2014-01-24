#include "pcim.h"
#include <iostream>
#include <vector>
#include <random>
#include <sstream>
#include "load.h"
#include <set>
#include "akipcim.h"

#define NPROC 8

using namespace std;
using namespace boost::archive;

pair<pcim,datainfo> buildmodel(string filename, int npts, int nproc) {
	dataset data;
	akiload(filename,data,npts);
	auto tests = akitests(data.info);
	return make_pair(pcim{data.ds,tests,akipcimparams(tests.size(),nproc)},
					data.info);
}

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


int main(int argc, char **argv) {
	if (argc<2) exit(0);
	datainfo info;
	pcim model;
	if (argv[1][0] < '0' || argv[1][0] > '9') {
		tie(model,info) = loadmodel(string(argv[1]));
	} else {
		double nproc = argc>2 ? atoi(argv[2]) : NPROC;
		tie(model,info) = buildmodel("/data/aki",atoi(argv[1]),nproc);
		savemodel("aki.pcim",model,info);
		ofstream dotf("aki.dot");
		model.todot(dotf,info);
	}
	model.print(cout,info);
	cout << endl;
}
