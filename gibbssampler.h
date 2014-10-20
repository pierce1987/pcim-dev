#ifndef PCIM_GIBBSSAMPLER_H
#define PCIM_GIBBSSAMPLER_H

#include "pcim.h"

using namespace std;

typedef vector<shptr<generic_state> > js;

struct ssumpcomp{
	bool operator() (const js &lhs, const js &rhs){
		for(int i = 0; i<lhs.size(); i++){
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
//Not used for now.
bool IsInUnobserved(std::vector<double> &starts, std::vector<double> &ends, double t);

bool GetPreviousState(map<js, vector<pair<js, pair<double,bool> > >, ssumpcomp> &transmap, js &jointstate, double prob);

class GibbsAuxSampler{

public:
	GibbsAuxSampler(const pcim *model, const ctbn::Trajectory *evidence, const ctbn::Context *contexts, int burnin);

	~GibbsAuxSampler();

	void SetTrajectory(const ctbn::Trajectory *traj);

	template<typename R>
	void SampleTrajectories(std::vector<ctbn::Trajectory> &traj, std::vector<double> &w,
					int numsamples, R &rand){	
		BurnIn(rand);
		cout<<"after burn in: "<<endl;
		//tr.AddTransition(0, 3.0, 0);//sampled from previous iteration
		printtr(cout,tr,3);
		for (int i=0; i<numsamples; ++i) {
			traj.push_back(Get());
			w.push_back(0.0);  // log weight
			Next(rand);
			//cout<<"after sample: "<<endl;
			printtr(cout,tr,3);
		}
	}

// tr is the output,which takes the evidence from the previous example, clears unobserved
// intervals, then add in event via thinning. oldtr is used in thinning process. temptr is
// used to take a jointstate and calculate rates. Initiliazed as tr (or oldtr), clear the whole
// sampled var traj, then insert events according to the state.
mutable ctbn::Trajectory tr,oldtr;

protected:
	template<typename R>
	void BurnIn(R &rand) const{
		if (burntin) return;
		if (!init_traj) {
			SampleInitialTrajectory(); //now only clean up evid and set it to tr
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
	void GetUnobservedIntervals(int varid) const;
	void GetAuxRates(int varid, int card) const;
	void ClearInitTraj();
	double getnextevent(double t0, int &event) const;
	bool IsVirtual(double t0, int event, int varid) const;
	double Getkeepprob(double rate, double t0) const;

	template<typename R>
	void Thinning(int varid, R &rand) const{
		if(starts.empty())
			return;
		std::uniform_real_distribution<> unifdist(0.0,1.0);
		//forward pass
		int T_event = 0; //event count

		/////
		js jointstate, jointstate1;//temp container of one joint state

		m->StateInit(jointstate);

		//cerr<<"testindex:"<<endl;
		for(int i=0; i<testindexes.size(); i++){
			//cerr<<testindexes[i]<<" ";
		}
		//cerr<<endl;

		for(int i = 0; i<jointstate.size(); i++){
			//cerr<<"content: ";
			//jointstate[i]->print();
			//cerr<<endl;
		}
		
		vector<double> times;
		//TODO: use typedef
		map<js, double, ssumpcomp> timestate; //state table at one time,temp
		map<js, vector<pair<js, pair<double,bool> > >, ssumpcomp> transmap;//state transition at one time
		vector<map<js, double, ssumpcomp> > allstates;
		vector<map<js, vector<pair<js, pair<double,bool> > >, ssumpcomp> > alltrans;
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
		int event = varid;
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
					double rate = -1.0;
					p_previous += m->Getlikelihood(varid, event, tr, jointstate, testindexes, own_var_list, t_previous, t0, rate, starts, ends); //log p
					//cerr<<"p_previous: "<<p_previous<<endl;
					//cerr<<"actual rate: "<<rate<<endl;
					double p_keep = Getkeepprob(rate, t0); //not log prob
					//double p_keep = Getkeepprob(m->getrate_test(event, t0, testindexes), t0);
					//cerr<<"p_keep: "<<log(p_keep)<<endl;
	
					//case: keep event
					jointstate1 = jointstate;//jointstate1: the previous state
					m->getnewstates(jointstate, testindexes, event, t0, 0, varid);

					auto it = timestate.find(jointstate); 
					if(it != timestate.end())
						//it->second = it->second + (p_previous*p_keep);
						it->second = log_add(it->second, (p_previous + log(p_keep)));
					else
						timestate.insert(pair<js, double>(jointstate, p_previous + log(p_keep)));		
					

					auto transit = transmap.find(jointstate); 
					if(transit != transmap.end()){
						(transit->second).push_back(pair<js, pair<double,bool> >(jointstate1, make_pair(p_previous + log(p_keep),true)));
					}
						
					else{
						vector<pair<js, pair<double,bool> > > tempvec;
						tempvec.push_back(pair<js, pair<double,bool> >(jointstate1, make_pair(p_previous + log(p_keep),true)));
						transmap.insert(pair<js, vector<pair<js, pair<double,bool> > > >(jointstate, tempvec));		
					}
					
					//case: do not keep event
					jointstate = jointstate1;
					m->getnewstates(jointstate, testindexes, -1, t0, 0, varid);
					it = timestate.find(jointstate); 
					if(it != timestate.end())
						it->second = log_add(it->second, (p_previous + log(1-p_keep)));
					else
						timestate.insert(pair<js, double>(jointstate, p_previous + log(1-p_keep)));		
					
					transit = transmap.find(jointstate); 
					if(transit != transmap.end()){
						(transit->second).push_back(pair<js, pair<double,bool> >(jointstate1, make_pair(p_previous + log(1 - p_keep),false)));
					}
						
					else{
						vector<pair<js, pair<double,bool> > > tempvec;
						tempvec.push_back(pair<js, pair<double,bool> >(jointstate1, make_pair(p_previous + log(1 - p_keep),false)));
						transmap.insert(pair<js, vector<pair<js, pair<double,bool> > > >(jointstate, tempvec));		
					}

				}	
			}
			else{//evidence
				for(auto iter = allstates[T_event-1].begin(); iter!=allstates[T_event-1].end();iter++){
					//cerr<<"!!!!!!!!!!!!!!"<<endl;
					jointstate = iter->first;
/*
				for(int i = 0; i<jointstate.size(); i++){
					cerr<<"content!!!!!!: ";
					jointstate[i]->print();
					cerr<<endl;
				}
*/

					jointstate1 = jointstate;
					double p_previous = iter->second;

					//temptr =  ctbn::Trajectory();
					//temptr.SetUnknown(varid,true);
					double rate = -1.0;
					//cerr<<"previous1: "<<p_previous<<endl;
					p_previous += m->Getlikelihood(varid, event, tr, jointstate, testindexes, own_var_list, t_previous, t0, rate, starts, ends);
					//cerr<<"previous2: "<<p_previous<<endl;
					p_previous += log(rate);//evidence needs this
					//cerr<<"previous3: "<<p_previous<<endl;


			 		m->getnewstates(jointstate, testindexes, event, t0, 0, varid);



					auto it = timestate.find(jointstate); 
					if(it != timestate.end())
						it->second = log_add(it->second, p_previous);
					else
						timestate.insert(pair<js, double>(jointstate, p_previous));	


					auto transit = transmap.find(jointstate); 
					if(transit != transmap.end()){
						(transit->second).push_back(pair<js, pair<double,bool> >(jointstate1, make_pair(p_previous,false)));
					}
						
					else{
						vector<pair<js, pair<double,bool> > > tempvec;
						tempvec.push_back(pair<js, pair<double,bool> >(jointstate1, make_pair(p_previous,false)));
						transmap.insert(pair<js, vector<pair<js, pair<double,bool> > > >(jointstate, tempvec));		
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
	T_event--;//4

	//p = unifdist(rand);//0.7, position should be adjusted?;move into the next loop
	



	for(;T_event>=0; T_event--){
	p = log(unifdist(rand));
	bool keep = GetPreviousState(alltrans[T_event], jointstate, p);	
	//cerr<<"time: "<<times[T_event]<<endl;
	if(keep)
		tr.AddTransition(varid, times[T_event], 0);
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
		//only use Context for the current variable
		ctbn::Context varcontext;
		varcontext.AddVar(var, context->Cardinality(var));
		//get unobserved intervals for the current variable (info in starts and ends)
		GetUnobservedIntervals(var);
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
		//cerr<<"BEFORE"<<endl;printtr(cout,oldtr,3);
		for(int i = 0; i < starts.size(); i++){ // for each unobserved intervals
			//from samplecomplete
			double t = starts[i];
			double T = ends[i];	
			double lastt=t;

			// tr->oldtr?? no, should use tr, which is fixed during sampling virtual events
			while((t = m->geteventaux(tr,lastt,expdist(rand),unifdist(rand),normdist(rand),var,T,varcontext,auxstarts,auxends,auxrates)) < T) {
				//cerr<<"sampled: "<<"var: "<<var<<" t: "<<t<<endl;
				oldtr.AddTransition(var, t, 0);	
				lastt = t; //proceed no matter event kept or not
			}

		}
		//cerr<<"AFTER"<<endl;printtr(cout,oldtr,3);


		//virtual events as (for testing)
		//oldtr.AddTransition(0, 3.5, 0);
		//oldtr.AddTransition(0, 4.0, 0);
		//oldtr.AddTransition(0, 4, 0);
		//oldtr.AddTransition(0, 8, 0);
		Clearcurrentvar(var);//clear events in unobserved areas for *tr*! Use thinning result to add events
		Thinning(var, rand); //thinning oldtr, in backward pass add events to tr.
	}

	int numBurninIter;
	const ctbn::Context *context;//contexts of all vars
	std::vector<int> testindexes;
	std::vector<int> own_var_list;
	mutable std::vector<double> starts; //times when unobserved intervals start, for the current variable
	mutable std::vector<double> ends; //times when unobserved intervals ends, for the current variable
	const pcim *m;
	const ctbn::Trajectory *evid;
	ctbn::Trajectory *init_traj; 
	//mutable ctbn::Trajectory tr,oldtr; 
	double begintime;
	double endtime;
	mutable bool burntin;

	mutable vector<double> auxstarts;
	mutable vector<double> auxends;
	mutable vector<double> auxrates;

private:
};


#endif
