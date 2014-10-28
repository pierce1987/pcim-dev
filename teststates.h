#ifndef TESTSTATES_H
#define TESTSTATES_H

#include "trajectory.h"
#include <memory>
#include <iostream>
#include <utility>
#include <cmath>
#include <queue>
#include <random>
#include <string>
#include <iterator>
#include "serial.h"
#include "load.h"
#include <vector>
#include <array>
#include <future>
#include <boost/make_shared.hpp>

class generic_state{ //for test state
public:
	virtual shptr<generic_state> getnewstate(shptr<generic_state>, double t, int var) const{};
	virtual shptr<generic_state> initialize() {};
	virtual std::string getsig() {return "empty";};
	virtual void print() const{std::cerr<<"nothing"<<std::endl;}
	virtual bool isequal(shptr<generic_state> rhs) const{};
	virtual bool islessthan(shptr<generic_state> rhs) const{};
};

// debugging only
class state_double : public generic_state{
public:
	state_double(){lasttime=-100;}
	state_double(double time) {lasttime = time;}
	virtual shptr<generic_state> initialize()  { return boost::make_shared<state_double>(); 
	}	
	virtual std::string getsig() {return "ssum_double";}	
	virtual void print() const{std::cerr<<"double: "<<lasttime<<std::endl;}
	virtual bool isequal(shptr<generic_state> rhs) const{
		return fabs(this->lasttime - boost::dynamic_pointer_cast<state_double>(rhs)->lasttime) < 0.00001;
	}
	virtual bool islessthan(shptr<generic_state> rhs) const{
		return this->lasttime < boost::dynamic_pointer_cast<state_double>(rhs)->lasttime - 0.00001;
	}	
	double lasttime;
};

// This is used to maintain the state of varcounttest. A queue is used to maintain events from current - maxlag
// to current time. 
class varcount_state : public generic_state{
public:
	varcount_state():times() {}
	varcount_state(std::queue<double> q) : times(q) {}
	virtual shptr<generic_state> initialize() {
		return boost::make_shared<varcount_state>(); 
	}
	// Debug info.
	virtual std::string getsig() {
		return "varcount_state";
	}
	virtual void print() const {
		std::cerr<<"time queue: ";
		std::queue<double> temp = times;
		while (!temp.empty()) {
			std::cerr<<temp.front()<<" ";
			temp.pop();		
		}
		std::cerr<<std::endl;
	}
	virtual bool isequal(shptr<generic_state> rhs) const{
		return (this->times == boost::dynamic_pointer_cast<varcount_state>(rhs)->times);
	}
	virtual bool islessthan(shptr<generic_state> rhs) const{
		return (this->times < boost::dynamic_pointer_cast<varcount_state>(rhs)->times);
	}
	// TODO maybe a queue of pair<double, int> for eventtypes? NO? 
	std::queue<double> times;
};

// This is used to maintain the state of eventcounttest. A queue is used to maintain events from current - maxlag
// to current time. 
class eventcount_state : public generic_state{
public:
	eventcount_state():times() {}
	eventcount_state(double time):times() {
		times.push(time);
	}
	eventcount_state(std::queue<double> q) : times(q) {}
	virtual shptr<generic_state> initialize() {
		return boost::make_shared<eventcount_state>(); 
	}
	// Debug info.
	virtual std::string getsig() {
		return "eventcount_state";
	}
	virtual void print() const {
		std::cerr<<"eventcount queue: ";
		std::queue<double> temp = times;
		while (!temp.empty()) {
			std::cerr<<temp.front()<<" ";
			temp.pop();		
		}
		std::cerr<<std::endl;
	}
	virtual bool isequal(shptr<generic_state> rhs) const{
		return (this->times == boost::dynamic_pointer_cast<eventcount_state>(rhs)->times);
	}
	virtual bool islessthan(shptr<generic_state> rhs) const{
		return (this->times < boost::dynamic_pointer_cast<eventcount_state>(rhs)->times);
	}

	std::queue<double> times;
};


#endif
