#ifndef TESTSTATES_H
#define TESTSTATES_H

#include "trajectory.h"
#include <memory>
#include <iostream>
#include <utility>
#include <cmath>
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

// This is used to maintain the state of varcounttest(1, 0, 3.0, 5.0)
class state_double1 : public generic_state{
public:
	state_double1(){lasttime=-100;}
	state_double1(double time) {lasttime = time;}
	virtual shptr<generic_state> initialize()  { return boost::make_shared<state_double1>(); 
	}
	virtual std::string getsig() {return "ssum_double1";}	
	virtual void print() const{std::cerr<<"double1: "<<lasttime<<std::endl;}
	virtual bool isequal(shptr<generic_state> rhs) const{
		return fabs(this->lasttime - boost::dynamic_pointer_cast<state_double1>(rhs)->lasttime) < 0.00001;
	}
	virtual bool islessthan(shptr<generic_state> rhs) const{
		return this->lasttime < boost::dynamic_pointer_cast<state_double1>(rhs)->lasttime - 0.00001;
	}
	double lasttime;
};


#endif
