#ifndef PCIM_GIBBSSAMPLER_H
#define PCIM_GIBBSSAMPLER_H

#include "pcim.h"

using namespace std;

typedef vector<shptr<generic_state> > js;

struct ssumpcomp{
	bool operator() (const js &lhs, const js &rhs){
		for(int i = 0; i < lhs.size(); i++){
			if(lhs[i]->isequal(rhs[i]))
				continue;
			else
				return lhs[i]->islessthan(rhs[i]);
		}
		return false;
	}
};

double log_add(double a, double b);
int sample_unnorm(vector<double> &input, double r); 
int GetPreviousState(map<js, vector<pair<js, pair<double,int>>>, ssumpcomp> &transmap, js &jointstate, double prob);

class GibbsAuxSampler{
public:
	GibbsAuxSampler(const pcim *model, const ctbn::Trajectory *evidence, const ctbn::Context *contexts, int burnin);

	~GibbsAuxSampler();

	void SetTrajectory(const ctbn::Trajectory *traj);

	template<typename R>
	void SampleTrajectories(std::vector<ctbn::Trajectory> &traj, std::vector<double> &w,
					int numsamples, R &rand){	
		BurnIn(rand);
		cerr<<"Burn In finished!"<<endl;
		printtr(cout,tr,3);
		for (int i=0; i<numsamples; ++i) {
			traj.push_back(Get());
			w.push_back(0.0);  // log weight
			Next(rand);
			printtr(cout,tr,3);
		}
	}

// tr is the output,which takes the evidence from the previous example, clears unobserved
// intervals, then add in event via thinning. oldtr is used in thinning process.
mutable ctbn::Trajectory tr,oldtr;

protected:
	template<typename R>
	void BurnIn(R &rand) const{
		if (burntin) return;
		if (!init_traj) {
			SampleInitialTrajectory(); //only clean up evid and set it to tr
		} else {
			tr = *init_traj;
		}
		Next(rand, numBurninIter);
		burntin = true;
	}

	// one step in the markov chain state transition; 
	template<typename R>
	void Next(R &rand, int num_iter = 1) const{
		for (int i=0; i<num_iter; ++i) {
			for (size_t var=0; var!=own_var_list.size(); ++var)
				SampleVariable(own_var_list[var], rand);
		}			
	}

	const ctbn::Trajectory &Get() const { return tr; }

	const ctbn::Trajectory &GetInitTraj() const { return *init_traj; }

	// Sample initial trajectory that agrees with evidence *evid.
	void SampleInitialTrajectory() const;
	void Clearcurrentvar(int var) const;
	void GetAuxRates(int varid, int card) const;
	void ClearInitTraj();
	double getnextevent(double t0, eventtype &event) const;
	bool IsVirtual(double t0, eventtype event, int varid) const;
	vector<double> Getkeepprob(vector<double> rates, double t0) const;

	// Assumes that the trajectory ends with an event
	template<typename R>
	void Thinning(int varid, R &rand) const{
		std::uniform_real_distribution<> unifdist(0.0,1.0);
		//forward pass
		int T_event = 0; //event count
		js jointstate, jointstate1;//temp container of one joint state
		m->StateInit(jointstate);

		for(int i = 0; i<jointstate.size(); i++){
			//cerr<<"content: ";
			//jointstate[i]->print();
			//cerr<<endl;
		}
		
		vector<double> times;
		//TODO: use typedef
		map<js, double, ssumpcomp> timestate; //state table at one time,temp
		map<js, vector<pair<js, pair<double,int>>>, ssumpcomp> transmap;//state transition at one time
		vector<map<js, double, ssumpcomp>> allstates;
		vector<map<js, vector<pair<js, pair<double,int>>>, ssumpcomp>> alltrans;
		timestate.insert(pair<js, double>(jointstate, 0.0));//initial state, with log prob
		allstates.push_back(timestate);
		//cerr<<"size: "<<jointstate.size()<<endl;
		//cerr<<"value:"<<timestate.begin()->second<<endl;

		//jointstate.clear();

		//m->StateInit(jointstate);

		//auto it = timestate.find(jointstate); 
		//if(it != timestate.end()){
		//	it->second = 1.0;
		//} 


		//cerr<<"size: "<<jointstate.size()<<endl;
		//cerr<<"value:"<<timestate.begin()->second<<endl;
		//cerr<<"MAPSIZE:"<<timestate.size()<<endl;

		double t0 = 0;//starts[0]-0.001; should not use starts[0] if doing tr insertion
		double t_previous = 0;
		eventtype event;
		for(;;) { //until the last event
		
			timestate.clear();
			transmap.clear();
			t0 = getnextevent(t0, event); //event gets the next event label
			//cerr<<"next event is "<<event<<" at "<<t0<<endl;
			if(t0 == -1.0) break;
			times.push_back(t0);
			T_event++; //push time (event count) forward
			bool isvirtual = IsVirtual(t0, event, varid);
			if(isvirtual){
				//get actual rate
				//cerr<<"event: "<<event<<"t0:"<<t0<<endl;
				
				//for each exsiting state
				//cerr<<"SIZE: "<<allstates[T_event-1].size()<<endl;
			/*for(auto iter = allstates[T_event-1].begin(); iter!=allstates[T_event-1].end();iter++){
				for(int i = 0; i<iter->first.size(); i++){
					cerr<<"content: ";
					iter->first[i]->print();
					cerr<<endl;
				}
				double Prob = iter->second;
				cerr<<"with log P: "<<Prob<<endl<<endl;
				//cerr<<"with P log: "<<log(Prob)<<endl<<endl;
			}*/
				for(auto iter = allstates[T_event-1].begin(); iter!=allstates[T_event-1].end();iter++){
					jointstate = iter->first;
					double p_previous = iter->second;
					//temptr =  ctbn::Trajectory();
					//temptr.SetUnknown(varid,true);
					vector<double> rates;
					p_previous += m->Getlikelihood(varid, event, tr, jointstate, testindexes, own_var_list, t_previous, t0, rates, *context, allstarts[varid], allends[varid]); //log p
					//cerr<<"p_previous: "<<p_previous<<endl;
					//cerr<<"actual rate: "<<rate<<endl;
					vector<double> p_keep = Getkeepprob(rates, t0); //not log prob
					//double p_keep = Getkeepprob(m->getrate_test(event, t0, testindexes), t0);
					//cerr<<"p_keep: "<<log(p_keep)<<endl;
					jointstate1 = jointstate;//jointstate1: the previous state
					//Case: keep event,consider all states
					for (int s = 0; s < p_keep.size(); ++s) {
						jointstate = jointstate1;
						m->getnewstates(jointstate, testindexes, eventtype(event.var, s), t0, 0, varid);

						auto it = timestate.find(jointstate); 
						if(it != timestate.end())
							it->second = log_add(it->second, (p_previous + log(p_keep[s])));
						else
							timestate.insert(pair<js, double>(jointstate, p_previous + log(p_keep[s])));		
					

						auto transit = transmap.find(jointstate); 
						if(transit != transmap.end()){
							(transit->second).push_back(pair<js, pair<double,int>>(jointstate1, make_pair(p_previous + log(p_keep[s]),s)));
						}
						
						else{
							vector<pair<js, pair<double,int>>> tempvec;
							tempvec.push_back(pair<js, pair<double,int>>(jointstate1, make_pair(p_previous + log(p_keep[s]),s)));
							transmap.insert(pair<js, vector<pair<js, pair<double,int>>>>(jointstate, tempvec));		
						}
					}
					//case: do not keep event
					jointstate = jointstate1;
					m->getnewstates(jointstate, testindexes, eventtype(-1,0), t0, 0, varid);
					auto it = timestate.find(jointstate); 
					double p_sum = accumulate(p_keep.begin(), p_keep.end(), 0);
					if(it != timestate.end())
						it->second = log_add(it->second, (p_previous + log(1 - p_sum)));
					else
						timestate.insert(pair<js, double>(jointstate, p_previous + log(1 - p_sum)));		
					
					auto transit = transmap.find(jointstate); 
					if(transit != transmap.end()) {
						(transit->second).push_back(pair<js, pair<double,int>>(jointstate1, make_pair(p_previous + log(1 - p_sum),-1)));
					}
						
					else{
						vector<pair<js, pair<double,int>>> tempvec;
						tempvec.push_back(pair<js, pair<double,int>>(jointstate1, make_pair(p_previous + log(1 - p_sum),-1)));
						transmap.insert(pair<js, vector<pair<js, pair<double,int>>>>(jointstate, tempvec));		
					}

				}	
			}
			else{//evidence
				for(auto iter = allstates[T_event-1].begin(); iter!=allstates[T_event-1].end();iter++){
					jointstate = iter->first;
/*				for(int i = 0; i<jointstate.size(); i++){
					cerr<<"content!!!!!!: ";
					jointstate[i]->print();
					cerr<<endl;
				}
*/
					jointstate1 = jointstate;
					double p_previous = iter->second;

					//temptr =  ctbn::Trajectory();
					//temptr.SetUnknown(varid,true);
					vector<double> rates;
					//cerr<<"previous1: "<<p_previous<<endl;
					p_previous += m->Getlikelihood(varid, event, tr, jointstate, testindexes, own_var_list, t_previous, t0, rates, *context, allstarts[varid], allends[varid]);
					//cerr<<"previous2: "<<p_previous<<endl;
					p_previous += log(rates[event.state]);//evidence needs this
					//cerr<<"previous3: "<<p_previous<<endl;


			 		m->getnewstates(jointstate, testindexes, event, t0, 0, varid);

					auto it = timestate.find(jointstate); 
					if(it != timestate.end())
						it->second = log_add(it->second, p_previous);
					else
						timestate.insert(pair<js, double>(jointstate, p_previous));	

					auto transit = transmap.find(jointstate); 
					if(transit != transmap.end()) {
						(transit->second).push_back(pair<js, pair<double,int>>(jointstate1, make_pair(p_previous,-1)));
					}
						
					else{
						vector<pair<js, pair<double,int>>> tempvec;
						tempvec.push_back(pair<js, pair<double,int>>(jointstate1, make_pair(p_previous,-1)));
						transmap.insert(pair<js, vector<pair<js, pair<double,int>>>>(jointstate, tempvec));		
					}	
					
				}
			}
			//cerr<<"time: "<<t0<<" event: "<<event<<endl;
			allstates.push_back(timestate);
			alltrans.push_back(transmap);
			t_previous = t0;

		}//end of forward pass

/*
		for(int i=0; i<allstates.size(); i++){
			cerr<<"states at time instance"<<i<<endl;
			for(auto iter = allstates[i].begin(); iter!=allstates[i].end();iter++){
				for(int i = 0; i<iter->first.size(); i++){
					cerr<<"content: ";
					iter->first[i]->print();
					cerr<<endl;
				}
				double Prob = iter->second;
				cerr<<"with log P: "<<Prob<<endl<<endl;
				//cerr<<"with P log: "<<log(Prob)<<endl<<endl;
			}
		}
		
		for(int i=0; i<alltrans.size(); i++){
			cerr<<"transition: states at time "<<i<<endl;
			for(auto iter = alltrans[i].begin(); iter!=alltrans[i].end();iter++){
				for (auto it2 = iter->second.begin(); it2 != iter->second.end(); it2++){
					//cerr<<"trans to "<<endl;

					for(int i = 0; i<iter->first.size(); i++){
						cerr<<"content: ";
						iter->first[i]->print();
						cerr<<endl;
					}

					cerr<<"from"<<endl;
					//it2->first->print();
					for(int i = 0; i<it2->first.size(); i++){
						cerr<<"content: ";
						it2->first[i]->print();
						cerr<<endl;
					}

					cerr<<"with P: "<<(it2->second.first)<<endl;
					cerr<<"keep?: "<<(it2->second.second)<<endl;
				}

			}
			cerr<<endl;
		}

		cerr<<T_event<<endl;
		for(int i=0; i<times.size(); i++){
			cerr<<times[i]<<" ";
		}
*/
	//backward pass
	//get final state	

	double p = log(unifdist(rand));//0.99;
/*
	double sump_final = 0.;
	for(auto iter = allstates[T_event].begin(); iter != allstates[T_event].end(); iter++){
		sump_final += iter->second;
	}
	for(auto iter = allstates[T_event].begin(); iter != allstates[T_event].end(); iter++){
		if(p <= iter->second/sump_final){
			jointstate = iter->first;
			break;
		}
		else
			p -= iter->second/sump_final;
	}
*/

        vector<double> logprobs;
	logprobs.reserve(allstates[T_event].size());
	for(auto iter = allstates[T_event].begin(); iter != allstates[T_event].end(); iter++){
		logprobs.push_back(iter->second);
	}

	auto tempit = allstates[T_event].begin();
	//cerr<<allstates[T_event].size()<<endl;
	//cerr<<logprobs.size()<<endl;
	advance(tempit,sample_unnorm(logprobs, p));

	jointstate = tempit->first;	
        
	/*cerr<<"sampled final state:"<<endl;
	for(int i = 0; i<jointstate.size(); i++){
		cerr<<"content: ";
		jointstate[i]->print();
		cerr<<endl;
	}
	cerr<<T_event<<endl;*/
	T_event--;

	for(;T_event >= 0; T_event--){
	p = log(unifdist(rand));
	int keep = GetPreviousState(alltrans[T_event], jointstate, p);	
	//cerr<<"time: "<<times[T_event]<<endl;
	if(keep == -2) {
		assert(1==0);
	}
	if(keep != -1)
		tr.AddTransition(varid, times[T_event], keep);
/*	cerr<<"previous state:"<<endl;
	for(int i = 0; i<jointstate.size(); i++){
		cerr<<"content: ";
		jointstate[i]->print();
		cerr<<endl;
	}
	cerr<<"keep? "<<keep<<endl;*/
		

	}
	

	}//end of Thinning

	// Resample the entire trajectory of v given all the other variables' full trajectory.
	template<typename R>
	void SampleVariable(int var, R &rand) const{
		if (allstarts[var].empty()) {
			return;
		}
		//only use Context for the current variable
		ctbn::Context varcontext;
		varcontext.AddVar(var, context->Cardinality(var));
		//get omega intervals (info in auxstarts, auxends, and auxrates)
		GetAuxRates(var, context->Cardinality(var));
		cerr<<"auxrates:"<<endl;
		for(int i = 0; i<auxstarts.size(); i++)
		{
			cerr<<auxrates[i]<<" in ( "<<auxstarts[i]<<","<<auxends[i]<<")"<<endl;
		}
		oldtr = tr;//use oldtr to maintain the previous sample

		std::exponential_distribution<> expdist(1.0);
		std::uniform_real_distribution<> unifdist(0.0,1.0);
		std::normal_distribution<> normdist(0.0,1.0);		
		for(int i = 0; i < allstarts[var].size(); i++) { // for each unobserved intervals
			double t = allstarts[var][i];
			double T = allends[var][i];	
			double lastt=t;
			while((t = m->geteventaux(tr,lastt,expdist(rand),unifdist(rand),normdist(rand),var,T,varcontext,auxstarts,auxends,auxrates)) < T) {				
				oldtr.AddTransition(var, t, 0);	
				cerr<<"sampled: "<<"var: "<<var<<" t: "<<t<<endl;
				lastt = t; //proceed no matter event kept or not
			}
		}

		Clearcurrentvar(var);//clear events in unobserved areas for *tr*. Use thinning result to add events
		Thinning(var, rand); //thinning oldtr, in backward pass add events to tr.
	}

	int numBurninIter;
	const ctbn::Context *context;//contexts of all vars
	std::vector<int> testindexes;
	std::vector<int> own_var_list;
	std::vector<std::vector<double>> allstarts; //times when unobserved intervals start, for all variables, initialized in constructor
	std::vector<std::vector<double>> allends; //times when unobserved intervals ends, for all variables
	const pcim *m;
	const ctbn::Trajectory *evid;
	ctbn::Trajectory *init_traj; 
	double begintime;
	double endtime;
	mutable bool burntin;
	mutable vector<double> auxstarts;
	mutable vector<double> auxends;
	mutable vector<double> auxrates;

private:
};


#endif
