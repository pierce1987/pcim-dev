#ifndef VARTRAJ_H
#define VARTRAJ_H

#include <map>
#include <boost/range/adaptor/reversed.hpp>
#include <vector>
#include <iostream>
#include <queue>

class vartraj {
public:
	struct rec {
		rec(double val) { v = val; }
		double v;
		operator double() const { return v; }
	};
	typedef std::map<double,rec> trajT;
	//typedef trajT::iterator iterator;
	typedef trajT::const_iterator const_iterator;

	vartraj() = default;

	vartraj(std::initializer_list<std::pair<double,double>> data) {
		for(auto &x : data)
			traj.insert(std::make_pair(x.first,rec(x.second)));
	}
	trajT::iterator begin() { return traj.begin(); }
	trajT::const_iterator begin() const { return traj.cbegin(); }
	trajT::const_iterator cbegin() const { return traj.cbegin(); }
	trajT::iterator end() { return traj.end(); }
	trajT::const_iterator end() const { return traj.cend(); }
	trajT::const_iterator cend() const { return traj.cend(); }

	trajT::reverse_iterator rbegin() { return traj.rbegin(); }
	trajT::const_reverse_iterator rbegin() const { return traj.crbegin(); }
	trajT::const_reverse_iterator crbegin() const { return traj.crbegin(); }
	trajT::reverse_iterator rend() { return traj.rend(); }
	trajT::const_reverse_iterator rend() const { return traj.crend(); }
	trajT::const_reverse_iterator crend() const { return traj.crend(); }

	bool empty() const { return traj.empty(); }

	trajT::const_iterator find(double t) const { return traj.find(t); }
	trajT::const_iterator lower_bound(double t) const { return traj.lower_bound(t); }
	trajT::const_iterator upper_bound(double t) const { return traj.upper_bound(t); }

	vartraj(const vartraj &t) : traj(t.traj) {
		tstart = t.tstart;
		tend = t.tend;
	}

	vartraj(const vartraj &t, double endt) {
		tstart = t.tstart;
		tend = endt;
		for(auto &x : t.traj) {
			if (x.first >= endt) break;
			traj.insert(x);
		}
	}
		

	vartraj(vartraj &&t) : traj(std::move(t.traj)) {
		tstart = t.tstart;
		tend = t.tend;
	}

	vartraj &operator=(const vartraj &t) {
		if (this==&t) return *this;
		traj = t.traj;
		tstart = t.tstart;
		tend = t.tend;
		return *this;
	}

	vartraj &operator=(vartraj &&t) {
		traj = std::move(t.traj);
		tstart = t.tstart;
		tend = t.tend;
		return *this;
	}

	void insert(double t, double v) {
		traj.insert(std::make_pair(t,rec(v)));
		if (tstart>t) tstart = t;
		if (tend<t) tend = t;
	}

	void reindex() {
	}

	double &starttime() { return tstart; }
	const double &starttime() const { return tstart; }
	double &endtime() { return tend; }
	const double &endtime() const { return tend; }

	void print(std::ostream &os) const {
		for(auto &x : *this) os << x.first << ": " << x.second << std::endl;
	}


private:

	double tstart,tend;
	trajT traj;
};


class traj : public std::vector<vartraj> {
public:
	// only initializes non-static vars
	traj(std::initializer_list<vartraj> data) {
		for(auto &x : data)
			push_back(x);
	}
	traj() = default;
	traj(int nevent) : std::vector<vartraj>(nevent) {
	}

	traj(const traj &tr, double endt) : sx(tr.sx) {
		for(auto &vtr : tr) emplace_back(vtr,endt);
	}
	
	std::vector<double> sx; // static attributes
};

void printtr(std::ostream &os, const traj &tr, bool incolumns=true, bool ishrs=false);

#endif
