#ifndef VARTRAJ_H
#define VARTRAJ_H

#include <map>
#include <boost/range/adaptor/reversed.hpp>
#include <vector>
#include <iostream>
#include <queue>

class vartraj {
public:

	typedef std::map<double,int> trajT;
	//typedef trajT::iterator iterator;
	typedef trajT::const_iterator const_iterator;

	vartraj() = default;

	vartraj(std::initializer_list<std::pair<double,int>> data) {
		for(auto &x : data)
			traj.insert(std::make_pair(x.first,x.second));
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

	void insert(double t, int v) {
		traj.insert(std::make_pair(t,v));
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


class Trajectory : public std::map<int, vartraj> {
public:
	int counter = 0;
	Trajectory(std::initializer_list<vartraj> data) {
		for(auto &x : data){
			insert(std::pair<int, vartraj>(counter,x));
			counter++;
		}
	}
	Trajectory() = default;
	//Trajectory(int nvar) : std::vector<vartraj>(nvar) {
	//}

	Trajectory(const Trajectory &tr, double endt){
		for(auto &vtr : tr) emplace(std::make_pair(vtr.first,vartraj(vtr.second,endt)));//??
	}
};

void printtr(std::ostream &os, const Trajectory &tr, bool incolumns=true, bool ishrs=false);

#endif
