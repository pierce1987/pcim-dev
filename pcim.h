#ifndef PCIM_H
#define PCIM_H

#include "traj.h"
#include <memory>
#include <iostream>
#include <utility>
#include <cmath>
#include <random>
#include <string>
#include "serial.h"
#include "load.h"
#include <vector>
#include <array>
#include <future>

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
		return std::make_pair(lhs.var, lhs.state)<std::make_pair(rhs.var, rhs.state);		
	}
};

struct vartrajrange {
	vartrajrange(const Trajectory *traject, eventtype e, double t0, double t1) : range(t0,t1), tr(traject) {
		event = e;
	}
	vartrajrange(const vartrajrange &vtr, double t0, double t1) : range(t0,t1) {
		tr = vtr.tr;
		event = vtr.event;
	}
	vartrajrange(const Trajectory *traject, int v){
		tr = traject;
		event.var = v;
		range.first = (*traject).find(v)->second.starttime(),
		range.second = (*traject).find(v)->second.endtime();
		
	}
	vartrajrange(const Trajectory *traject, eventtype e){
		tr = traject;
		event = e;
		range.first = (*traject).find(e.var)->second.starttime(),
		range.second = (*traject).find(e.var)->second.endtime();
	}
	const Trajectory *tr;
	eventtype event;
	timerange range;
};

class pcimtest {
public:
	virtual ~pcimtest() {} ;
	virtual void print(std::ostream &os) const = 0;
	virtual void print(std::ostream &os, const datainfo &info) const = 0;
	// adds to (does not replace) outtrue and outfalse
	virtual void chop(const vartrajrange &in,
			std::vector<vartrajrange> &outtrue,
			std::vector<vartrajrange> &outfalse) const = 0;
	virtual bool eval(const Trajectory &tr, eventtype event, double t) const {
		double toss; return eval(tr,event,t,toss);
	}
	virtual bool eval(const Trajectory &tr, eventtype event, double t, double &until) const = 0;
private:
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
	}
};

// tests if last state of testvar == teststate (if no last state, state taken to be 0)
class lasttest : public pcimtest {
public:
	lasttest(int testvar=0, int teststate=0) {
		v = testvar;
		state = teststate;
	}
	virtual ~lasttest() {}
	virtual void print(std::ostream &os) const {
		os << "most recent " << v << " == " << state;
	}
	virtual void print(std::ostream &os, const datainfo &info) const {
		os << "most recent state for " << info.dvarname(v) << " == " << state;
	}
	virtual void chop(vartrajrange &in,
			std::vector<vartrajrange> &outtrue,
			std::vector<vartrajrange> &outfalse) const {
		const vartraj &tr = (*(in.tr)).find(v==-1?in.event.var:v)->second;
		const auto e = tr.cend();
		double t0 = in.range.first;
		double tend = in.range.second;
		auto i0 = tr.upper_bound(t0);
		if (i0!=tr.begin()) --i0;
		auto i1 = tr.lower_bound(tend);
		bool currval = (i0==e ? 0 : i0->second) == state;
		double t1 = i0->first;
		while(t1<tend && i0!=e) {
			++i0;
			while (i0!=e && i0->first<=t1) ++i0;
			if (i0==e || i0->first>=tend) t1 = tend;
			else t1 = i0->first;
			bool newval = (i0==e ? 0 : i0->second) == state;
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
	virtual bool eval(const Trajectory &tr, eventtype event, double t) const {
		const vartraj &vtr = tr.find(v==-1?event.var:v)->second;
		if (vtr.empty()) return 0 == state;
		auto i0 = vtr.lower_bound(t);
		if (i0==vtr.cend() || i0->first>t) --i0;
		return (i0==vtr.cend() ? 0 : i0->second) == state;
	}

	virtual bool eval(const Trajectory &tr, eventtype event, double t, double &until) const {
		const vartraj &vtr = tr.find(v==-1?event.var:v)->second;
		if (vtr.empty()) {
			until = std::numeric_limits<double>::infinity();
			return 0 == state;
		}
		auto i0 = vtr.lower_bound(t);
		auto e = vtr.cend();
		auto i1 = i0;
		if (i0==e || i0->first>t) --i0;
		else ++i1;
		until = i1!=e ? i1->first : std::numeric_limits<double>::infinity();
		assert(until>t);
		return (i0==e ? 0 : i0->second) == state;
	}

protected:
	int state;
	int v;
private:
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(pcimtest);
		ar & BOOST_SERIALIZATION_NVP(v) & BOOST_SERIALIZATION_NVP(state);
	}
};

class timetest : public pcimtest {
public:
	timetest(double tstart=0, double tend=1, double mod=1) {
		t0 = tstart; t1 = tend; m = mod;
	}
	virtual ~timetest() {}
	virtual void print(std::ostream &os) const {
		os << "time (%" << m << ") in (" << t0 << ',' << t1 << ')';
	}
	virtual void print(std::ostream &os, const datainfo &info) const {
		os << "time (%" << m << ") in (" << t0 << ',' << t1 << ')';
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
		double end = breakup(in.range.second,endbase);
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
		
	virtual bool eval(const Trajectory &tr, eventtype event, double t) const {
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
	virtual bool eval(const Trajectory &tr, eventtype event, double t, double &until) const {
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

private:
	double t0,t1,m;
private:
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(pcimtest);
		ar & BOOST_SERIALIZATION_NVP(t0) & BOOST_SERIALIZATION_NVP(t1) & BOOST_SERIALIZATION_NVP(m);
	}
};


template<typename D>
class varstattest : public pcimtest {
public:
	varstattest(int testvar=0, double lag0=0, double lag1=1, int teststate=0) {
		v = testvar;
		s = teststate;
		maxlag = std::max(lag0,lag1);
		minlag = std::min(lag0,lag1);
	}
	virtual ~varstattest() {}
	virtual void print(std::ostream &os) const = 0;
	virtual void print(std::ostream &os, const datainfo &info) const =0;
	virtual void chop(const vartrajrange &in,
			std::vector<vartrajrange> &outtrue,
			std::vector<vartrajrange> &outfalse) const {
		const vartraj &tr = (*(in.tr)).find(v==-1?in.event.var:v)->second;
		const auto &e = tr.cend();
		double t0 = in.range.first;
		double tend = in.range.second;
		auto i0 = tr.upper_bound(t0-maxlag);
		auto i1 = tr.upper_bound(t0-minlag);
		typename D::statT stat;
		for(auto i=i0;i!=i1;i++) stat.add(i->first,i->second);
		double t1 = t0;
		bool currval = static_cast<const D *>(this)->evalstat(stat);
		while(t1<tend && i0!=e) {
			t1 = std::min(i0==e ?
				std::numeric_limits<double>::infinity() : i0->first+maxlag,
						i1==e ?
				std::numeric_limits<double>::infinity() : i1->first+minlag);
			while (i0!=e && i0->first+maxlag<=t1) {
				stat.del(i0->first,i0->second);
				++i0;
			}
			while (i1!=e && i1->first+minlag<=t1) {
				stat.add(i1->first,i1->second);
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
	virtual bool eval(const Trajectory &tr, eventtype event, double t) const {
		const vartraj &vtr = tr.find(v==-1?event.var:v)->second;
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
		for(auto i=i0;i!=i1;i++) stat.add(i->first,i->second);
		return static_cast<const D *>(this)->evalstat(stat);
	}

	virtual bool eval(const Trajectory &tr, eventtype event, double t, double &until) const {
		const vartraj &vtr = tr.find(v==-1?event.var:v)->second;
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
		for(auto i=i0;i!=i1;i++) stat.add(i->first,i->second);
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
	

// test if count of number of events of var testvar (-1 == currvar)
// from t-lag0 to t-lag1 is greater than thresh
class counttest : public varstattest<counttest> {
public:
	counttest(int thresh=0, int testvar=0, double lag0=0, double lag1=1, int teststate = 0)
			: varstattest<counttest>(testvar,lag0,lag1, teststate) { theta=thresh; }
	virtual ~counttest() {}
	virtual void print(std::ostream &os) const {
		os << "# " << v << " in [" << maxlag << ',' << minlag << ") >= "
				<< theta;
	}
	virtual void print(std::ostream &os, const datainfo &info) const {
		os << "# " << info.dvarname(v) << " measurements in [" 
			<< maxlag << ',' << minlag << ") >= " << theta;
	}

	struct statT {
		statT() { n=0; }
		int n;
		void add(double,int) { n++; }
		void del(double,int) { n--; }
	};

	bool evalstat(const statT &s) const {
		return s.n>=theta;
	}
	
protected:
	int theta;
private:
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
		//ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(varstattest<counttest>);
		ar & boost::serialization::make_nvp("varstattest",
			boost::serialization::base_object<varstattest<counttest>>(*this));
		ar & BOOST_SERIALIZATION_NVP(theta);
	}
};

//new, associate to eventtype
template<typename D>
class eventstattest : public pcimtest {
public:
	eventstattest(int testvar=0,int teststate=0, double lag0=0, double lag1=1) {
		v = testvar;
		s = teststate;
		maxlag = std::max(lag0,lag1);
		minlag = std::min(lag0,lag1);
	}
	virtual ~eventstattest() {}
	virtual void print(std::ostream &os) const = 0;
	virtual void print(std::ostream &os, const datainfo &info) const =0;
	virtual void chop(const vartrajrange &in,
			std::vector<vartrajrange> &outtrue,
			std::vector<vartrajrange> &outfalse) const {
		const vartraj &tr = (*(in.tr)).find(v==-1?in.event.var:v)->second;
		const auto &e = tr.cend();
		double t0 = in.range.first;
		double tend = in.range.second;
		auto i0 = tr.upper_bound(t0-maxlag);
		auto i1 = tr.upper_bound(t0-minlag);
		typename D::statE stat;
		for(auto i=i0;i!=i1;i++) stat.add(i->first,i->second,s);
		double t1 = t0;
		bool currval = static_cast<const D *>(this)->evalstat(stat);
		while(t1<tend && i0!=e) {
			t1 = std::min(i0==e ?
				std::numeric_limits<double>::infinity() : i0->first+maxlag,
						i1==e ?
				std::numeric_limits<double>::infinity() : i1->first+minlag);
			while (i0!=e && i0->first+maxlag<=t1) {
				stat.del(i0->first,i0->second,s);
				++i0;
			}
			while (i1!=e && i1->first+minlag<=t1) {
				stat.add(i1->first,i1->second,s);
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
	virtual bool eval(const Trajectory &tr, eventtype event, double t) const {
		const vartraj &vtr = tr.find(v==-1?event.var:v)->second;
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

		typename D::statE stat;
		for(auto i=i0;i!=i1;i++) stat.add(i->first,i->second,s);
		return static_cast<const D *>(this)->evalstat(stat);
	}

	virtual bool eval(const Trajectory &tr, eventtype event, double t, double &until) const {
		const vartraj &vtr = tr.find(v==-1?event.var:v)->second;
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
		typename D::statE stat;
		for(auto i=i0;i!=i1;i++) stat.add(i->first,i->second,s);
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
	


// test if count of number of events of var testvar (-1 == currvar) at state teststate
// from t-lag0 to t-lag1 is greater than thresh
class counteventtest : public eventstattest<counteventtest> {
public:
	counteventtest(int thresh=0, int teststate = 0, int testvar=0, double lag0=0, double lag1=1)
			: eventstattest<counteventtest>(testvar,teststate,lag0,lag1) { theta=thresh;}
	virtual ~counteventtest() {}
	virtual void print(std::ostream &os) const {
		os << "# " << v << " in [" << maxlag << ',' << minlag << ") >= "
				<< theta;
	}
	virtual void print(std::ostream &os, const datainfo &info) const {
		os << "# " << info.dvarname(v) << " measurements in [" 
			<< maxlag << ',' << minlag << ") >= " << theta;
	}

	struct statE {
		statE() { n=0; }
		int n;
		void add(double t,int s,int teststate) { if(teststate == s) n++; }
		void del(double t,int s,int teststate) { if(teststate == s) n--; }
	};

	bool evalstat(const statE &s) const {
		return s.n>=theta;
	}
	
protected:
	int theta;
private:
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
		ar & boost::serialization::make_nvp("eventstattest",
			boost::serialization::base_object<eventstattest<counteventtest>>(*this));
		ar & BOOST_SERIALIZATION_NVP(theta);
	}
};

// test if current variable == testvar
class vartest : public pcimtest {
public:
	vartest(int testvar=0) : pcimtest() { v = testvar; };
	virtual ~vartest() {} ;
	virtual void print(std::ostream &os) const { os << "var == " << v; }
	virtual void print(std::ostream &os, const datainfo &info) const {
		os << "X == " << info.dvarname(v);
	}
	virtual void chop(const vartrajrange &in,
			std::vector<vartrajrange> &outtrue,
			std::vector<vartrajrange> &outfalse) const {
		if (in.event.var==v) outtrue.emplace_back(in);
		else outfalse.emplace_back(in);
	}
	virtual bool eval(const Trajectory &tr, eventtype event, double t) const {
		return event.var==v;
	}
	virtual bool eval(const Trajectory &tr, eventtype event, double t, double &until) const {
		until = std::numeric_limits<double>::infinity();
		return event.var==v;
	}
private:
	int v;
private:
	friend class boost::serialization::access;
	template<typename Ar>
	void serialize(Ar &ar, const unsigned int ver) {
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(pcimtest);
		ar & BOOST_SERIALIZATION_NVP(v);
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

	pcim(const std::vector<Trajectory> &data, const std::vector<shptr<pcimtest>> &tests,
		const pcimparams &params, const std::vector<int> &states);

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
	double samplecomplete(Trajectory &ret, double T, R &rand, std::vector<int> &states) const {
		double t=0.0;
		int var;
		int state;
		std::exponential_distribution<> expdist(1.0);
		std::uniform_real_distribution<> unifdist(0.0,1.0);
		std::normal_distribution<> normdist(0.0,1.0);
		double lastt=t;
		while((t = getevent(ret,lastt,expdist(rand),unifdist(rand),normdist(rand),var,state,T,states))<T) {
			ret[var].insert(t,state);
			lastt = t;
		}
		return lastt;
	}

	template<typename R>
	Trajectory sample(double T,int nvar,R &rand, std::vector<int> &states) const {
		Trajectory ret;
		if (T<0.0) return ret;
		for(int i=0;i<nvar;i++) {
			ret[i].starttime() = 0.0;
			ret[i].endtime() = T;
		}
		samplecomplete(ret,T,rand,states);
		return ret;
	}

	// returns relevant leaves in ret and sum as return value
	double getrate(const Trajectory &tr, double t, double &until, std::map<eventtype, const pcim *, eventcomp> &ret, const std::vector<int> &states) const;
	// returns new time and sets var and val to the variable and its value
	double getevent(const Trajectory &tr, double &t, double expsamp, double unisamp, double normsamp,
					int &var, int &state, double maxt, const std::vector<int> &states) const;

	void print(std::ostream &os) const;
	void print(std::ostream &os, const datainfo &info) const;
	void todot(std::ostream &os, const datainfo &info) const;

	void save(std::ostream &os) const;
	void load(std::ostream &os);

	std::vector<std::string> featurenames() const {
		std::vector<std::string> ret;
		featurenames(ret,"");
		return ret;
	}
	std::vector<double> trajtofeatures(const Trajectory &tr) const {
		std::vector<double> ret;
		std::vector<vartrajrange> vtr;
		for(int v=0;v<tr.size();v++)
			vtr.emplace_back(&tr,v);
		trajtofeatures(std::vector<vartrajrange>{vtr},ret);
		return ret;
	}
	double similarity(const Trajectory &tr1, const Trajectory &tr2) const {//measure similarity based on features
		std::vector<vartrajrange> vtr1;
		for(int v=0;v<tr1.size();v++)
			vtr1.emplace_back(&tr1,v);
		std::vector<vartrajrange> vtr2;
		for(int v=0;v<tr2.size();v++)
			vtr2.emplace_back(&tr2,v);
		return similarity(std::vector<vartrajrange>{vtr1},
					std::vector<vartrajrange>{vtr2});
	}
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

	double getratevar(const Trajectory &tr, int var, int state, double t, double &until, const pcim *&leaf) const;

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



	void trajtofeatures(const std::vector<vartrajrange> &tr,
				std::vector<double> &f) const;


	double similarity(const std::vector<vartrajrange> &tr1,
				const std::vector<vartrajrange> &tr2) const;

	shptr<pcim> ftree,ttree;
	shptr<pcimtest> test;
	double rate;
	ss stats;
	int globalm; // number of traj it was built from

private:
	void getleaffeature(const std::vector<vartrajrange> &tr,
			std::array<double,nleaffeat> &f) const;

	void featurenames(std::vector<std::string> &ret, std::string prefix) const;
	std::array<std::string,nleaffeat> getleaffeaturenames() const;

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
BOOST_CLASS_EXPORT_KEY(counttest)
BOOST_CLASS_EXPORT_KEY(varstattest<counttest>)
BOOST_CLASS_EXPORT_KEY(counteventtest)
BOOST_CLASS_EXPORT_KEY(eventstattest<counteventtest>)
BOOST_CLASS_EXPORT_KEY(vartest)
BOOST_CLASS_EXPORT_KEY(pcim)
BOOST_SERIALIZATION_SHARED_PTR(pcimtest)
BOOST_SERIALIZATION_SHARED_PTR(pcim)

BOOST_CLASS_VERSION(pcim::ss,1)


#endif
