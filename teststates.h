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

// This is used to maintain the state of lasttest. An int is used to maintain the last state of testvar
// default value is 0 (even no last state, consider it to be 0) 
class last_state : public generic_state{
public:
	last_state() { laststate = 0; }
	last_state(int i) { laststate = i; }
	virtual shptr<generic_state> initialize() {
		return boost::make_shared<last_state>(); 
	}
	// Debug info.
	virtual std::string getsig() {
		return "last_state";
	}
	virtual void print() const {
		std::cerr<<"last state: ";
		std::cerr<<laststate<<std::endl;
	}
	virtual bool isequal(shptr<generic_state> rhs) const{
		return (this->laststate == boost::dynamic_pointer_cast<last_state>(rhs)->laststate);
	}
	virtual bool islessthan(shptr<generic_state> rhs) const{
		return (this->laststate < boost::dynamic_pointer_cast<last_state>(rhs)->laststate);
	}

	int laststate;
};

// This is used to maintain the state of timetest. No data structure is needed.
class time_state : public generic_state{
public:
	time_state() {}
	virtual shptr<generic_state> initialize() {
		return boost::make_shared<time_state>(); 
	}
	// Debug info.
	virtual std::string getsig() {
		return "time_state";
	}
	virtual void print() const {
		std::cerr<<"time state: ";
	}
	virtual bool isequal(shptr<generic_state> rhs) const{
		return true;
	}
	virtual bool islessthan(shptr<generic_state> rhs) const{
		return false;
	}
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

	std::queue<double> times;
};

// This is used to maintain the state of eventcounttest. A queue is used to maintain events from current - maxlag
// to current time. 
class eventcount_state : public generic_state{
public:
	eventcount_state():times() {}
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

// This is used to maintain the state of vartest. No data structure is needed.
class var_state : public generic_state{
public:
	var_state() {}
	virtual shptr<generic_state> initialize() {
		return boost::make_shared<var_state>(); 
	}
	// Debug info.
	virtual std::string getsig() {
		return "var_state";
	}
	virtual void print() const {
		std::cerr<<"var state: ";
	}
	virtual bool isequal(shptr<generic_state> rhs) const{
		return true;
	}
	virtual bool islessthan(shptr<generic_state> rhs) const{
		return false;
	}
};

// This is used to maintain the state of eventtest. No data structure is needed.
class event_state : public generic_state{
public:
	event_state() {}
	virtual shptr<generic_state> initialize() {
		return boost::make_shared<event_state>(); 
	}
	// Debug info.
	virtual std::string getsig() {
		return "event_state";
	}
	virtual void print() const {
		std::cerr<<"event state: ";
	}
	virtual bool isequal(shptr<generic_state> rhs) const{
		return true;
	}
	virtual bool islessthan(shptr<generic_state> rhs) const{
		return false;
	}
};

#endif
