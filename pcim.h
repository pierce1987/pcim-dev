#ifndef PCIM_H
#define PCIM_H

#include "trajectory.h"
#include "teststates.h"
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
#include <queue>
#include <boost/make_shared.hpp>

namespace boost { namespace serialization { 
	class access;
}}

typedef std::pair<double,double> timerange;

struct eventtype{
	eventtype(int v = 0, int s = 0) { var = v; state = s;}
	int var;
	int state;
};

struct eventcomp{
	bool operator() (const eventtype& lhs, const eventtype& rhs) const{
		return std::make_pair(lhs.var, lhs.state) < std::make_pair(rhs.var, rhs.state);		
	}
};

struct vartrajrange {
	vartrajrange(const ctbn::Trajectory *traject, eventtype e, double t0, double t1) : range(t0,t1), tr(traject) {
		event = e;
	}
	vartrajrange(const vartrajrange &vtr, double t0, double t1) : range(t0,t1) {
		tr = vtr.tr;
		event = vtr.event;
	}
	vartrajrange(const ctbn::Trajectory *traject, int v){
		tr = traject;
		event.var = v;
		range.first = (*traject).TimeBegin(),
		range.second = (*traject).TimeEnd();
		
	}
	vartrajrange(const ctbn::Trajectory *traject, eventtype e){
		tr = traject;
		event = e;
		range.first = (*traject).TimeBegin(),
		range.second = (*traject).TimeEnd();
	}
	const ctbn::Trajectory *tr;
	eventtype event;
	timerange range;
};

class pcimtest {
public:
	virtual ~pcimtest() {} ;
	virtual void print(std::ostream &os) const = 0;
	virtual void print(std::ostream &os, const datainfo &info) const = 0;
	virtual int getauxv() const = 0;
	// adds to (does not replace) outtrue and outfalse
	virtual void chop(const vartrajrange &in,
			std::vector<vartrajrange> &outtrue,
			std::vector<vartrajrange> &outfalse) const = 0;
	virtual bool eval(const ctbn::Trajectory &tr, eventtype event, double t) const {
		double toss; return eval(tr,event,t,toss);
	}
	virtual bool eval(const ctbn::Trajectory &tr, eventtype event, double t, double &until) const = 0;
	// This assumes that the queue has been updated correctly, for example, for eventcounttest, the queue should only contain times of the queried event. Also, if auxv = -1, the queue should only have events of the sampled var. 
	virtual bool neweval(const ctbn::Trajectory &tr, shptr<generic_state> &teststate, int varid, eventtype event, double t, double &until) const {};
	virtual shptr<generic_state> stateupdate(shptr<generic_state> &teststate, eventtype event, double t0, int varid) const{};
	virtual generic_state* getteststate() = 0;
	//virtual void updatetraj(shptr<generic_state> teststate, ctbn::Trajectory &temptr, int varid) {};

private:
	
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
	}
};

// tests if last state of testvar == teststate (if no last state, s taken to be 0)
// testvar can be -1.
class lasttest : public pcimtest {
public:
	lasttest(int testvar=0, int teststate=0) {
		v = testvar;
		s = teststate;
		auxv = v;
	}
	virtual ~lasttest() {}
	virtual void print(std::ostream &os) const {
		os << "most recent " << v << " == " << s;
	}
	virtual void print(std::ostream &os, const datainfo &info) const {
		os << "most recent state for " << info.dvarname(v) << " == " << s;
	}
	virtual int getauxv() const {
		return auxv;
	}
	virtual void chop(const vartrajrange &in,
			std::vector<vartrajrange> &outtrue,
			std::vector<vartrajrange> &outfalse) const {
		const ctbn::VarTrajectory &tr = (*(in.tr)).GetVarTraj(v==-1?in.event.var:v);
		const auto e = tr.cend();
		double t0 = in.range.first;
		double tend = in.range.second;
		auto i0 = tr.upper_bound(t0);
		if (i0!=tr.begin()) --i0;
		auto i1 = tr.lower_bound(tend);
		bool currval = (i0==e ? 0 : i0->second) == s;
		double t1 = i0->first;
		while(t1<tend && i0!=e) {
			++i0;
			while (i0!=e && i0->first<=t1) ++i0;
			if (i0==e || i0->first>=tend) t1 = tend;
			else t1 = i0->first;
			bool newval = (i0==e ? 0 : i0->second) == s;
			if (!newval && currval) {
				outtrue.emplace_back(in,t0,t1);
				t0 = t1;
				currval = false;
			} else if (newval && !currval) {
				outfalse.emplace_back(in,t0,t1);
				t0 = t1;
				currval = true;
			}
		}
		if (t0<tend) {
			if (currval) outtrue.emplace_back(in,t0,tend);
			else outfalse.emplace_back(in,t0,tend);
		}
	}
	virtual bool eval(const ctbn::Trajectory &tr, eventtype event, double t) const {
		const ctbn::VarTrajectory &vtr = tr.GetVarTraj(v==-1?event.var:v);
		if (vtr.empty()) return 0 == s;
		auto i0 = vtr.lower_bound(t);
		if (i0==vtr.begin()) { 
			if (t==i0->first) {
				return s == i0->second;
			}
			return 0 == s;
		} // try to fix the next line. 
		if (i0==vtr.cend() || i0->first>t) --i0; //if i0 points to the first event, --i0 would do nothing, in the next line i0->second may not be 0
		return (i0==vtr.cend() ? 0 : i0->second) == s;
	}

	virtual bool eval(const ctbn::Trajectory &tr, eventtype event, double t, double &until) const {
		const ctbn::VarTrajectory &vtr = tr.GetVarTraj(v==-1?event.var:v);
		if (vtr.empty()) {
			until = std::numeric_limits<double>::infinity();
			return 0 == s;
		}
		auto i0 = vtr.lower_bound(t);
		auto e = vtr.cend();
		auto i1 = i0;
		if (i0 == vtr.begin()) {
			if (t==i0->first) {
				++i1;
				until = i1!=e ? i1->first : std::numeric_limits<double>::infinity();
				assert(until>t);
				return i0->second == s;
			}
			until = i0->first;
			assert(until>t);
			return 0 == s;
		}
		if (i0==e || i0->first>t) --i0;
		else ++i1;
		until = i1!=e ? i1->first : std::numeric_limits<double>::infinity();
		assert(until>t);
		return (i0==e ? 0 : i0->second) == s;
	}

	virtual generic_state* getteststate() {return &teststate;}

	// The int needs to maintain last state of testvar.
	// auxv = -1 means test current var. event = -1 means do not keep the virtual event.
	virtual shptr<generic_state> stateupdate(shptr<generic_state> &state, eventtype event, double t0, int varid) const{
		
		if ((auxv != -1 && auxv != varid) || event.var != varid) {
			return state;
		}
		return boost::make_shared<last_state>(event.state); 	
	}

	virtual bool neweval(const ctbn::Trajectory &tr, shptr<generic_state> &state, int varid, eventtype event, double t, double &until) const {
		if (auxv != varid && (auxv != -1 || event.var != varid)) {
			return eval(tr, event, t, until);
		}
		// use state
		until = std::numeric_limits<double>::infinity();
		int lstate = boost::dynamic_pointer_cast<last_state>(state)->laststate;	
		return lstate == s;	
	}

protected:
	int s;
	int v;
	int auxv; //used to decide aux rate, -2 means does not care 
	last_state teststate;
private:
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(pcimtest);
		ar & BOOST_SERIALIZATION_NVP(v) & BOOST_SERIALIZATION_NVP(s) & BOOST_SERIALIZATION_NVP(auxv);
	}
};

class timetest : public pcimtest {
public:
	timetest(double tstart=0, double tend=1, double mod=1) {
		t0 = tstart; t1 = tend; m = mod; auxv = -2;
	}
	virtual ~timetest() {}
	virtual void print(std::ostream &os) const {
		os << "time (%" << m << ") in (" << t0 << ',' << t1 << ')';
	}
	virtual void print(std::ostream &os, const datainfo &info) const {
		os << "time (%" << m << ") in (" << t0 << ',' << t1 << ')';
	}
	virtual int getauxv() const {
		return auxv;
	}
	inline double remerge(double base, double inc) const {
		return base+inc;
	}
        
	inline double breakup(double v, double &n) const {
		double ret = std::remainder(v,m);
		if (ret<0) ret += m;
		n = v - ret;
		return ret;
	}
			
	virtual void chop(const vartrajrange &in,
			std::vector<vartrajrange> &outtrue,
			std::vector<vartrajrange> &outfalse) const {
		int igt = outtrue.size();
		int igf = outfalse.size();
		double temp;
		double myt0 = breakup(t0,temp);
		double myt1 = breakup(t1,temp);
		double startbase,endbase;
		double start = breakup(in.range.first,startbase);
		double end = breakup(in.range.second,endbase);	// This assumes that the queue has been updated correctly, for example, for eventcounttest, the queue should only contain times of the queried event. Also, if auxv = -1, the queue should only have events of the sampled var. 
		double tmin,tmax;
		std::tie(tmin,tmax) = std::minmax(myt0,myt1);
		std::vector<vartrajrange> &out0 = myt0<myt1 ? outfalse : outfalse;
		std::vector<vartrajrange> &out1 = myt0<myt1 ? outtrue: outtrue;
		auto startpt = in.range.first; //remerge(startbase,start);
		while(startbase<endbase) {
			if (start<tmin) {
				auto nextpt = remerge(startbase,tmin);
				if (startpt<nextpt)
					out0.emplace_back(in,startpt,nextpt);
				startpt = nextpt;
				start = tmin;
			}
			if (start<tmax) {
				auto nextpt = remerge(startbase,tmax);
				if (startpt<nextpt)
					out1.emplace_back(in,startpt,nextpt);
				startpt = nextpt;
				start = tmax;
			}
			start -= m;
			startbase += m;
		}
		if (start<tmin) {
			if (end<tmin) {
				auto nextpt = in.range.second;
				if (startpt<nextpt)
					out0.emplace_back(in,startpt,nextpt);
				//dump(in,outtrue,outfalse,igt,igf);
				return;
			} else {
				auto nextpt = remerge(startbase,tmin);
				if (startpt<nextpt)
					out0.emplace_back(in,startpt,nextpt);
				start = tmin;
				startpt = nextpt;
			}
		}
		if (start<tmax) {
			if (end<tmax) {
				auto nextpt = in.range.second;
				if (startpt<nextpt)
					out1.emplace_back(in,startpt,nextpt);
				//dump(in,outtrue,outfalse,igt,igf);
				return;
			} else {
				auto nextpt = remerge(startbase,tmax);
				if (startpt<nextpt)
					out1.emplace_back(in,startpt,nextpt);
				start = tmax;
				startpt = nextpt;
			}
		}
		if (start<end) {
			auto nextpt = in.range.second;
			if (startpt<nextpt)
				out0.emplace_back(in,startpt,nextpt);
		}
		//dump(in,outtrue,outfalse,igt,igf);
	}

	void dump(const vartrajrange &in,
			std::vector<vartrajrange> &outtrue,
			std::vector<vartrajrange> &outfalse,
			int igt, int igf) const {
		std::cout << "%" << m << " in (" << t0 <<','<< t1 <<")" << std::endl;
		std::cout << "in: " << in.range.first << ' ' << in.range.second << std::endl;
		std::cout << "outtrue: " << std::endl;
		for(int i=igt;i<outtrue.size();i++) {
			auto &x = outtrue[i];
			std::cout << "\t" << x.range.first << ' ' << x.range.second << std::endl;
		}
		std::cout << "outfalse: " << std::endl;
		for(int i=igf;i<outfalse.size();i++) {
			auto &x = outfalse[i];
			std::cout << "\t" << x.range.first << ' ' << x.range.second << std::endl;
		}
	}
		
	virtual bool eval(const ctbn::Trajectory &tr, eventtype event, double t) const {
		double temp;
		double myt0 = breakup(t0,temp);
		double myt1 = breakup(t1,temp);
		double tbase;
		double tmod = breakup(t,tbase);
		double tmin,tmax;
		std::tie(tmin,tmax) = std::minmax(myt0,myt1);
		if (tmod<tmin || tmod>tmax) return myt0>myt1;
		return myt1>myt0;
	}
	virtual bool eval(const ctbn::Trajectory &tr, eventtype event, double t, double &until) const {
		double temp;
		double myt0 = breakup(t0,temp);
		double myt1 = breakup(t1,temp);
		double tbase;
		double tmod = breakup(t,tbase);
		double tmin,tmax;
		std::tie(tmin,tmax) = std::minmax(myt0,myt1);
		if (tmod<tmin) {
			until = remerge(tbase,tmin);
			assert(until>t);
			return myt0>myt1;
		}
		if (tmod<tmax) {
			until = remerge(tbase,tmax);
			assert(until>t);
			return myt1>myt0;
		}
		until = remerge(tbase,tmin)+m;
		assert(until>t);
		return myt0>myt1;
	}

	virtual generic_state* getteststate() {return &teststate;}

	// The state is always the same
	virtual shptr<generic_state> stateupdate(shptr<generic_state> &state, eventtype event, double t0, int varid) const{
		return state;	
	}

	virtual bool neweval(const ctbn::Trajectory &tr, shptr<generic_state> &state, int varid, eventtype event, double t, double &until) const {
		return eval(tr, event, t, until);	
	}
private:
	double t0,t1,m;
	int auxv; // -2, do not relate to event
	time_state teststate;
private:
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(pcimtest);
		ar & BOOST_SERIALIZATION_NVP(t0) & BOOST_SERIALIZATION_NVP(t1) & BOOST_SERIALIZATION_NVP(m) & BOOST_SERIALIZATION_NVP(auxv);
	}
};

//abstract class
template<typename D>
class varstattest : public pcimtest {
public:
	varstattest(int testvar=0, double lag0=0, double lag1=1, int teststate=-1) {//if teststate omitted (== -1), we only care about var, not eventtypes.
		v = testvar;
		s = teststate;
		maxlag = std::max(lag0,lag1);
		minlag = std::min(lag0,lag1);
	}
	virtual ~varstattest() {}
	virtual void print(std::ostream &os) const = 0;
	virtual void print(std::ostream &os, const datainfo &info) const =0;
	virtual int getauxv() const = 0;
	virtual void chop(const vartrajrange &in,
			std::vector<vartrajrange> &outtrue,
			std::vector<vartrajrange> &outfalse) const {
		const ctbn::VarTrajectory &tr = (*(in.tr)).GetVarTraj(v==-1?in.event.var:v);
		const auto &e = tr.cend();
		double t0 = in.range.first;
		double tend = in.range.second;
		auto i0 = tr.upper_bound(t0-maxlag);
		auto i1 = tr.upper_bound(t0-minlag);
		typename D::statT stat;
		for(auto i=i0;i!=i1;i++) if(s == -1 || i->second == s) stat.add(i->first,i->second);

		double t1 = t0;
		bool currval = static_cast<const D *>(this)->evalstat(stat);
		while(t1<tend && i0!=e) {
			t1 = std::min(i0==e ?
				std::numeric_limits<double>::infinity() : i0->first+maxlag,
						i1==e ?
				std::numeric_limits<double>::infinity() : i1->first+minlag);

			while (i0!=e && i0->first+maxlag<=t1) {
				if(s == -1 || i0->second == s) stat.del(i0->first,i0->second);
				++i0;
			}
			while (i1!=e && i1->first+minlag<=t1) {
				if(s == -1 || i1->second == s) stat.add(i1->first,i1->second);
				++i1;
			}
			if (t1>=tend) t1 = tend;
			bool newval = static_cast<const D *>(this)->evalstat(stat);
			if (t0<t1) {
				if (!newval && currval) {
					outtrue.emplace_back(in,t0,t1);
					t0 = t1;
					currval = false;
				} else if (newval && !currval) {
					outfalse.emplace_back(in,t0,t1);
					t0 = t1;
					currval = true;
				}
			}
		}
		if (t0<tend) {
			if (currval) outtrue.emplace_back(in,t0,tend);
			else outfalse.emplace_back(in,t0,tend);
		}

	}
	virtual bool eval(const ctbn::Trajectory &tr, eventtype event, double t) const {
		const ctbn::VarTrajectory &vtr = tr.GetVarTraj(v==-1?event.var:v);
		double tnext
			= std::nextafter(tnext,std::numeric_limits<double>::infinity());
		double t0 = tnext-maxlag;
		if (t0==t-maxlag)
			t0 = std::nextafter(t0,std::numeric_limits<double>::infinity());
		double t1 = tnext-minlag;
		if (t1==t-minlag)
			t1 = std::nextafter(t1,std::numeric_limits<double>::infinity());
		auto i0 = vtr.lower_bound(t0);
		auto i1 = vtr.lower_bound(t1);

		typename D::statT stat;
		for(auto i=i0;i!=i1;i++) if(s == -1 || i->second == s) stat.add(i->first,i->second);
		return static_cast<const D *>(this)->evalstat(stat);
	}

	virtual bool eval(const ctbn::Trajectory &tr, eventtype event, double t, double &until) const {
		const ctbn::VarTrajectory &vtr = tr.GetVarTraj(v==-1?event.var:v);
		const auto &e = vtr.cend();
		double tnext
			= std::nextafter(t,std::numeric_limits<double>::infinity());
		double t0 = tnext-maxlag;
		if (t0==t-maxlag)
			t0 = std::nextafter(t0,std::numeric_limits<double>::infinity());
		double t1 = tnext-minlag;
		if (t1==t-minlag)
			t1 = std::nextafter(t1,std::numeric_limits<double>::infinity());
		auto i0 = vtr.lower_bound(t0);
		auto i1 = vtr.lower_bound(t1);
		typename D::statT stat;
		for(auto i=i0;i!=i1;i++) if((s == -1 || i->second == s)) {
			//std::cerr<<"adding.."<<std::endl;
			stat.add(i->first,i->second);
		}
		until = std::min(i0!=e ? i0->first+maxlag
						: std::numeric_limits<double>::infinity(),
				i1!=e ? i1->first+minlag
						: std::numeric_limits<double>::infinity());
		assert(until>t);
		return static_cast<const D *>(this)->evalstat(stat);
	}

protected:
	double minlag,maxlag;
	int v;
	int s;
private:
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(pcimtest);
		ar & BOOST_SERIALIZATION_NVP(minlag) & BOOST_SERIALIZATION_NVP(maxlag) & BOOST_SERIALIZATION_NVP(v)& BOOST_SERIALIZATION_NVP(s);
	}
};
	

// test if count of number of events of var testvar (-1 == currvar) from t-lag0 to t-lag1
// is >= thresh. If we care teststate (not -1), should use eventcounttest
class varcounttest : public varstattest<varcounttest> {
public:
	varcounttest(int thresh=0, int testvar=0, double lag0=0, double lag1=1, int teststate = -1)
			: varstattest<varcounttest>(testvar,lag0,lag1,teststate) { theta=thresh; auxv = testvar;}
	virtual ~varcounttest() {}
	virtual void print(std::ostream &os) const {
		os << "# " << v << " in [" << maxlag << ',' << minlag << ") >= "
				<< theta;
	}
	virtual void print(std::ostream &os, const datainfo &info) const {
		os << "# " << info.dvarname(v) << " measurements in [" 
			<< maxlag << ',' << minlag << ") >= " << theta;
	}
	virtual int getauxv() const{
		return auxv;
	}
	struct statT {
		statT() { n=0; }
		int n;
		void add(double,int) { n++; }
		void del(double,int) { n--; }
	};

	bool evalstat(const statT &s) const {
		return s.n >= theta;
	}

	double getmaxlag() const{
		return maxlag;
	}

	virtual generic_state* getteststate() {return &teststate;}		

	// The queue needs to maintain every instance from t0 - minlag until t0, between 
	// t0 - maxlag and t0 - minlag, only at most theta instances are needed.
	// auxv = -1 means test current var in the test. event = -1 means do not keep the virtual event.
	virtual shptr<generic_state> stateupdate(shptr<generic_state> &state, eventtype event, double t0, int varid) const{
		
		if (auxv != -1 && auxv != varid) {
			return state;
		}
		std::queue<double> eventtimes = boost::dynamic_pointer_cast<varcount_state>(state)->times;
		if (event.var == varid) { // takes care of event.var == -1
			eventtimes.push(t0);
		}		
		while (!eventtimes.empty() && eventtimes.front() < t0 - maxlag) {
			eventtimes.pop();	
		}
		// In many cases minlag == 0, can save some time
		if (minlag == 0) {
			while (eventtimes.size() > theta) {
				eventtimes.pop();
			}
			return boost::make_shared<varcount_state>(eventtimes); 			
		}
		std::queue<double> result;
		while (!eventtimes.empty() && eventtimes.front() >= t0 - maxlag && eventtimes.front() <= t0 - minlag) {		
			result.push(eventtimes.front());
			eventtimes.pop();
			if (result.size() > theta) {
				result.pop();
			}
		}
		while (!eventtimes.empty()) {		
			result.push(eventtimes.front());
			eventtimes.pop();
		}
		return boost::make_shared<varcount_state>(result); 	
	}

	virtual bool neweval(const ctbn::Trajectory &tr, shptr<generic_state> &state, int varid, eventtype event, double t, double &until) const {
		if (auxv != varid && (auxv != -1 || event.var != varid)) {
			return eval(tr, event, t, until);
		}
		// use state
		std::queue<double> eventtimes = boost::dynamic_pointer_cast<varcount_state>(state)->times;	
		double i0_t = -1.0;
		double i1_t = -1.0;
		int count = 0;
		while (!eventtimes.empty()) {
			if (eventtimes.front() >= t - maxlag && eventtimes.front() <= t - minlag) {
				count++;
				if (i0_t == -1.0) {
					i0_t = eventtimes.front();
				}			
			}
			if (eventtimes.front() > t - minlag) {
				if (i0_t == -1.0) {
					i0_t = eventtimes.front();
				}
				if (i1_t == -1.0) {
					i1_t = eventtimes.front();
				}
			}
			eventtimes.pop();
		}

		until = std::min(i0_t != -1.0 ? std::nextafter(i0_t + maxlag, std::numeric_limits<double>::infinity())
				: std::numeric_limits<double>::infinity(),
				i1_t != -1.0 ? std::nextafter(i1_t + minlag, std::numeric_limits<double>::infinity())
				: std::numeric_limits<double>::infinity());


		//std::cerr<<"until: "<<until<<" t: "<<t<<std::endl;
		assert(until>t);
		return count >= theta;	
	}
	
protected:
	int theta;
	int auxv;
	varcount_state teststate;
private:
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
		//ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(varstattest<counttest>);
		ar & boost::serialization::make_nvp("varstattest",
			boost::serialization::base_object<varstattest<varcounttest>>(*this));
		ar & BOOST_SERIALIZATION_NVP(theta) & BOOST_SERIALIZATION_NVP(auxv);
	}
};

// test if count of number of events of var testvar (-1 == currvar) at state teststate
// from t-lag0 to t-lag1 is >= thresh
class eventcounttest : public varstattest<eventcounttest> {
public:
	eventcounttest(int thresh=0, int testvar=0, double lag0=0, double lag1=1, int teststate = 0)
			: varstattest<eventcounttest>(testvar,lag0,lag1,teststate) { theta=thresh; auxv = testvar;}
	virtual ~eventcounttest() {}
	virtual void print(std::ostream &os) const {
		os << "# (" << v <<","<<s<<")"<< " in [" << maxlag << ',' << minlag << ") >= "
				<< theta;
	}
	virtual void print(std::ostream &os, const datainfo &info) const {
		os << "# " << info.dvarname(v) << " measurements in [" 
			<< maxlag << ',' << minlag << ") >= " << theta;
	}
	virtual int getauxv() const{
		return auxv;
	}
	struct statT {
		statT() { n=0; }
		int n;
		void add(double t,int s) { n++; }
		void del(double t,int s) { n--; }
	};

	bool evalstat(const statT &s) const {
		return s.n >= theta;
	}
	virtual generic_state* getteststate() {return &teststate;}
	// The queue needs to maintain every instance from t0 - minlag until t0, between 
	// t0 - maxlag and t0 - minlag, only at most theta instances are needed.
	// auxv = -1 means test current var in the test. event = -1 means do not keep the virtual event.
	virtual shptr<generic_state> stateupdate(shptr<generic_state> &state, eventtype event, double t0, int varid) const{
		
		if (auxv != -1 && auxv != varid) { //this condition is not generic
			return state;
		}
		std::queue<double> eventtimes = boost::dynamic_pointer_cast<eventcount_state>(state)->times;
		if (event.var == varid && event.state == s) { // takes care of event.var == -1
			eventtimes.push(t0);
		}		
		while (!eventtimes.empty() && eventtimes.front() < t0 - maxlag) {
			eventtimes.pop();	
		}
		// In many cases minlag == 0, can save some time
		if (minlag == 0) {
			while (eventtimes.size() > theta) {
				eventtimes.pop();
			}
			return boost::make_shared<eventcount_state>(eventtimes); 			
		}
		std::queue<double> result;
		while (!eventtimes.empty() && eventtimes.front() >= t0 - maxlag && eventtimes.front() <= t0 - minlag) {		
			result.push(eventtimes.front());
			eventtimes.pop();
			if (result.size() > theta) {
				result.pop();
			}
		}
		while (!eventtimes.empty()) {		
			result.push(eventtimes.front());
			eventtimes.pop();
		}
		return boost::make_shared<eventcount_state>(result); 	
	}

	virtual bool neweval(const ctbn::Trajectory &tr, shptr<generic_state> &state, int varid, eventtype event, double t, double &until) const {
		// use state when ask about the current sampled var. That would be auxv==varid, or auxv==-1 && event.var=varid
		if (auxv != varid && (auxv != -1 || event.var != varid)) {
			return eval(tr, event, t, until);
		}
		// use state
		std::queue<double> eventtimes = boost::dynamic_pointer_cast<eventcount_state>(state)->times;	
		double i0_t = -1.0;
		double i1_t = -1.0;
		int count = 0;
		while (!eventtimes.empty()) {
			if (eventtimes.front() >= t - maxlag && eventtimes.front() <= t - minlag) {
				count++;
				if (i0_t == -1.0) {
					i0_t = eventtimes.front();
				}			
			}
			if (eventtimes.front() > t - minlag) {
				if (i0_t == -1.0) {
					i0_t = eventtimes.front();
				}
				if (i1_t == -1.0) {
					i1_t = eventtimes.front();
				}
			}
			eventtimes.pop();
		}

		until = std::min(i0_t != -1.0 ? std::nextafter(i0_t + maxlag, std::numeric_limits<double>::infinity())
				: std::numeric_limits<double>::infinity(),
				i1_t != -1.0 ? std::nextafter(i1_t + minlag, std::numeric_limits<double>::infinity())
				: std::numeric_limits<double>::infinity());


		//std::cerr<<"until: "<<until<<" t: "<<t<<std::endl;
		assert(until>t);
		return count >= theta;	
	}	
protected:
	int theta;
	int auxv;
	eventcount_state teststate;
private:
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
		ar & boost::serialization::make_nvp("varstattest",
			boost::serialization::base_object<varstattest<eventcounttest>>(*this));
		ar & BOOST_SERIALIZATION_NVP(theta)& BOOST_SERIALIZATION_NVP(auxv);
	}
};

// test if current variable == testvar
class vartest : public pcimtest {
public:
	vartest(int testvar=0) : pcimtest() { v = testvar; auxv = -2;};
	virtual ~vartest() {} ;
	virtual void print(std::ostream &os) const { os << "var == " << v; }
	virtual void print(std::ostream &os, const datainfo &info) const {
		os << "X == " << info.dvarname(v);
	}
	virtual int getauxv() const{
		return auxv;
	}
	virtual void chop(const vartrajrange &in,
			std::vector<vartrajrange> &outtrue,
			std::vector<vartrajrange> &outfalse) const {
		if (in.event.var==v) outtrue.emplace_back(in);
		else outfalse.emplace_back(in);
	}
	virtual bool eval(const ctbn::Trajectory &tr, eventtype event, double t) const {
		return event.var==v;
	}
	virtual bool eval(const ctbn::Trajectory &tr, eventtype event, double t, double &until) const {
		until = std::numeric_limits<double>::infinity();
		return event.var==v;
	}
	virtual generic_state* getteststate() {return &teststate;}

	// The state is always the same
	virtual shptr<generic_state> stateupdate(shptr<generic_state> &state, eventtype event, double t0, int varid) const{
		return state;	
	}
	virtual bool neweval(const ctbn::Trajectory &tr, shptr<generic_state> &state, int varid, eventtype event, double t, double &until) const {
		return eval(tr, event, t, until);	
	}
private:
	int v;
	int auxv;
	var_state teststate;
private:
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(pcimtest);
		ar & BOOST_SERIALIZATION_NVP(v)& BOOST_SERIALIZATION_NVP(auxv);
	}
};

// test if current event == testevent, this supports testvar to be -1, thus only checks state.
class eventtest : public pcimtest {
public:
	eventtest(int testvar=0, int teststate=0) : pcimtest() { v = testvar; s = teststate; auxv = -2;};
	virtual ~eventtest() {} ;
	virtual void print(std::ostream &os) const { os << "event == (" << v <<"," << s <<")"; }
	virtual void print(std::ostream &os, const datainfo &info) const {
		os << "E == (" << info.dvarname(v)<<","<<info.dvarname(s)<<")";
	}
	virtual int getauxv() const{
		return auxv;
	}
	virtual void chop(const vartrajrange &in,
			std::vector<vartrajrange> &outtrue,
			std::vector<vartrajrange> &outfalse) const {
		if (in.event.var == v && in.event.state==s) outtrue.emplace_back(in);
		else outfalse.emplace_back(in);
	}
	virtual bool eval(const ctbn::Trajectory &tr, eventtype event, double t) const {
		if(v!=-1)return (event.var == v && event.state==s);
		return (event.state==s);
	}
	virtual bool eval(const ctbn::Trajectory &tr, eventtype event, double t, double &until) const {
		until = std::numeric_limits<double>::infinity();
		if(v!=-1)return (event.var == v && event.state==s);
		return (event.state==s);
	}
	virtual generic_state* getteststate() {return &teststate;}
	// The state is always the same
	virtual shptr<generic_state> stateupdate(shptr<generic_state> &state, eventtype event, double t0, int varid) const{
		return state;	
	}
	virtual bool neweval(const ctbn::Trajectory &tr, shptr<generic_state> &state, int varid, eventtype event, double t, double &until) const {
		return eval(tr, event, t, until);	
	}
private:
	int v;
	int s;
	int auxv;
	event_state teststate;
private:
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(pcimtest);
		ar & BOOST_SERIALIZATION_NVP(v);
		ar & BOOST_SERIALIZATION_NVP(s);
		ar & BOOST_SERIALIZATION_NVP(auxv);
	}
};


class pcim {
public:
	typedef double wtT;
	typedef double wtvarT;
	static constexpr double zerowt = 0.0;
	static constexpr double zerowtvar = 0.0;
	inline static double quad(const wtT &x, const wtvarT &M) {
		return x*M*x;
	}
	inline static wtT div(const wtT &x, const wtvarT &M) {
		return x/M;
	}

	class pcimparams {
	friend pcim;
	public:
		// alpha -- pseudo event counts
		// beta -- pseudo durations
		// kappa -- kappa^n is prior on tree of size n (kappa<=1!)
		pcimparams(double alpha, double beta, double kappa, int minnumevents=0, int numproc=1) {
			a = alpha; b = beta; k = kappa;
			alb = alpha*::log(beta);
			lga = std::lgamma(alpha);
			lk = ::log(kappa);
			nproc = numproc;
			mne = minnumevents;
		}
	private:
		double a,b,k;
		double alb; // alpha*log(beta)
		double lga; // log(gamma(alpha))
		double lk; // log(kappa)
		int nproc,mne;
	};

	virtual ~pcim() {
	}

	pcim(const std::vector<ctbn::Trajectory> &data, const std::vector<shptr<pcimtest>> &tests,
		const pcimparams &params, const ctbn::Context &contexts);

	pcim(shptr<pcimtest> tst, 
		shptr<pcim> truebranch, shptr<pcim> falsebranch)
			: test(tst), ttree(truebranch), ftree(falsebranch) {
		rate = 1.0; 
	}
	pcim(pcimtest *tst, pcim *truebranch, pcim *falsebranch)
			: test(tst), ttree(truebranch), ftree(falsebranch) {
		rate = 1.0; 
	}

	pcim(double lambda=1.0): ttree(), ftree(), test() {
		rate =lambda;
		globalm = 0;
	}

	template<typename R>
	double samplecomplete(ctbn::Trajectory &ret, double T, R &rand, ctbn::Context &contexts) const {
		double t=0.0;
		int var;
		int state;
		std::exponential_distribution<> expdist(1.0);
		std::uniform_real_distribution<> unifdist(0.0,1.0);
		std::normal_distribution<> normdist(0.0,1.0);
		double lastt=t;
		while((t = getevent(ret,lastt,expdist(rand),unifdist(rand),normdist(rand),var,state,T,contexts))<T) {
			ret.AddTransition(var, t, state);
			lastt = t;
		}
		return lastt;
	}

	template<typename R>
	ctbn::Trajectory sample(double T,int nvar,R &rand, ctbn::Context &contexts) const {
		ctbn::Trajectory ret(nvar);
		if (T<0.0) return ret;

		ret.SetBeginTime(0.0);
		ret.SetEndTime(T);
		samplecomplete(ret,T,rand,contexts);
		return ret;
	}

	// returns relevant leaves in ret and sum as return value
	double getrate(const ctbn::Trajectory &tr, double t, double &until, std::map<eventtype, const pcim *, eventcomp> &ret, const ctbn::Context &contexts) const;
	double getauxrates(const ctbn::Trajectory &tr, double &t, int card, double &until, double &r, double varid) const;
	// returns new time and sets var and val to the variable and its value
	double getevent(const ctbn::Trajectory &tr, double &t, double expsamp, double unisamp, double normsamp,
					int &var, int &state, double maxt, const ctbn::Context &contexts) const;
	double geteventaux(const ctbn::Trajectory &tr, double &t, double expsamp, double unisamp, double normsamp,
					int &var, double maxt, const ctbn::Context &contexts, std::vector<double> &auxstarts, std::vector<double> &auxends, std::vector<double> &auxrates) const;
	double Getlikelihood(int varid, eventtype event, ctbn::Trajectory &tr, std::vector<shptr<generic_state> > &jointstate, const std::vector<int> &testindexes, const std::vector<int> &own_var_list, double t_previous, double t0, std::vector<double> &rates, const ctbn::Context &contexts, const std::vector<double> &starts, const std::vector<double> &ends) const;
	void StateInit(std::vector<shptr<generic_state> > &jointstate) const;
	int Makeindex(std::vector<int> &indexes, int i) const;
	void getnewstates(std::vector<shptr<generic_state> > &jointstate, const std::vector<int> &testindexes, eventtype event, double t0, int index, int varid) const;
	int counttest() const;
	void print(std::ostream &os) const;
	void print(std::ostream &os, const datainfo &info) const;
	void todot(std::ostream &os, const datainfo &info) const;

	void save(std::ostream &os) const;
	void load(std::ostream &os);


private:
	class ss {
	public:
		double n; // not int, in case need expected value
		double t; // total time
	private:
		friend class boost::serialization::access;
		template<typename Ar>
		void serialize(Ar &ar, const unsigned int ver) {
			if (ver==0) {
				ar & BOOST_SERIALIZATION_NVP(n) & BOOST_SERIALIZATION_NVP(t);
			} else {
				ar & BOOST_SERIALIZATION_NVP(n) & BOOST_SERIALIZATION_NVP(t);
			}
		}
	};
	static constexpr size_t nleaffeat=3;

	pcim(const std::vector<vartrajrange> &data, const pcim::ss &s,
		const std::vector<shptr<pcimtest>> &tests,
		double basescore, const pcimparams &params, int m);


	void build(const std::vector<vartrajrange> &data, const ss &s,
		const std::vector<shptr<pcimtest>> &tests,
		double basescore, const pcimparams &params);
	
	static ss suffstats(const std::vector<vartrajrange> &data);
	static double score(const ss &d, const pcimparams &p);
	void calcleaf(const ss &d, const pcimparams &p);

	double getratevar(const ctbn::Trajectory &tr, int var, int state, double t, double &until, const pcim *&leaf) const;
	double getratevar_state(const ctbn::Trajectory &tr, std::vector<shptr<generic_state> > &jointstate, int varid, eventtype testevent, double t, double &until, const std::vector<int> &testindexes, int index) const;
	double getratevaraux(const ctbn::Trajectory &tr, int varid, int state, double t, double &until) const;


	void printhelp(std::ostream &os, int lvl, const datainfo *info=nullptr) const;
	void todothelp(std::ostream &os, int par, bool istrue, int &nn, const datainfo &info) const;

	struct testpick {
		int testnum;
		shptr<pcimtest> test;
		std::vector<vartrajrange> d1,d2;
		double s1,s2;
		ss ss1,ss2;
	};

	testpick picktest(const std::vector<vartrajrange> &data,
			const std::vector<shptr<pcimtest>> &tests,
			const pcimparams &params,
			int procnum=0, int nproc=1) const;

	shptr<pcim> ftree,ttree;
	shptr<pcimtest> test;
	double rate;
	ss stats;
	int globalm; // number of traj it was built from

private:

	friend class boost::serialization::access;

	template<typename Ar>
	void save(Ar &ar, const unsigned int ver) const {
		ar & BOOST_SERIALIZATION_NVP(globalm);
		ar & BOOST_SERIALIZATION_NVP(rate);
		ar & BOOST_SERIALIZATION_NVP(stats);
		ar & BOOST_SERIALIZATION_NVP(test);
		ar & BOOST_SERIALIZATION_NVP(ttree);
		ar & BOOST_SERIALIZATION_NVP(ftree);
	}
	template<typename Ar>
	void load(Ar &ar, const unsigned int ver) {
		ar & BOOST_SERIALIZATION_NVP(globalm);
		ar & BOOST_SERIALIZATION_NVP(rate);
		ar & BOOST_SERIALIZATION_NVP(stats);
		ar & BOOST_SERIALIZATION_NVP(test);
		ar & BOOST_SERIALIZATION_NVP(ttree);
		ar & BOOST_SERIALIZATION_NVP(ftree);
	}
	BOOST_SERIALIZATION_SPLIT_MEMBER()
};

BOOST_CLASS_EXPORT_KEY(pcimtest)
BOOST_CLASS_EXPORT_KEY(lasttest)
BOOST_CLASS_EXPORT_KEY(timetest)
BOOST_CLASS_EXPORT_KEY(varcounttest)
BOOST_CLASS_EXPORT_KEY(varstattest<varcounttest>)
BOOST_CLASS_EXPORT_KEY(eventcounttest)
BOOST_CLASS_EXPORT_KEY(varstattest<eventcounttest>)
BOOST_CLASS_EXPORT_KEY(vartest)
BOOST_CLASS_EXPORT_KEY(eventtest)
BOOST_CLASS_EXPORT_KEY(pcim)
BOOST_SERIALIZATION_SHARED_PTR(pcimtest)
BOOST_SERIALIZATION_SHARED_PTR(pcim)

BOOST_CLASS_VERSION(pcim::ss,1)


#endif
