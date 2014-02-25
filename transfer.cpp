#include "transfer.h"

using namespace std;

pcim* rwogenerator(int var, int state, int ds){
	double rate;
	cin>>rate;
	if(state == ds-2){
		double last;
		cin>>last;
		pcim *temp = new pcim(new eventtest(var, state), new pcim(rate<0? 0:rate), new pcim(last<0? 0: last));
		return temp;
	}
	else{
		pcim *temp = new pcim(new eventtest(var, state), new pcim(rate<0? 0:rate), rwogenerator(var, state+1, ds));
		return temp;
	}		

}

pcim* matrixgenerator(int var, int state, int ds){
		
	if(state == ds-2){
		pcim* temp = new pcim(new lasttest(var, state), rwogenerator(var, 0, ds), rwogenerator(var, 0, ds));
		return temp;
	}	
	else{	
		pcim* temp = new pcim(new lasttest(var, state), rwogenerator(var, 0, ds), matrixgenerator(var, state+1, ds));
		return temp;
	}
}

pcim* nodegenerator(int var, int parentindex, int state, ctbn::Context &contexts, int ds){//only when var has parents
	int id = contexts.VarList()[parentindex];
	if(parentindex == contexts.VarList().size()-1 && (state == contexts.Cardinality(id)-1)){
		pcim *temp = matrixgenerator(var, 0, ds);
		return temp;
	}	

	if(parentindex == contexts.VarList().size()-1){
		pcim *temp = new pcim(new lasttest(id, state), matrixgenerator(var, 0, ds), nodegenerator(var, parentindex, state+1, contexts, ds));
		return temp;
	}
		
	pcim *temp = new pcim(new lasttest(id, state), nodegenerator(var, parentindex+1, 0, contexts, ds), nodegenerator(var, parentindex, state+1, contexts, ds));
	return temp;
	
}

pcim* transfer(int NumofNodes, int count){
	
	ctbn::Context contexts; //parent context

	int var;
	cin>>var;
	int ds;
	cin>>ds;
	int num;
	cin>>num;	

	for(int i=0; i<num; i++){
		int id;
		int card;
		cin>>id;
		cin>>card;
		contexts.AddVar(id, card);
	}

	if(count == NumofNodes){
		if(num == 0){//no parent
			pcim* temp = matrixgenerator(var, 0, ds);	
			return temp;		
		}	
		else{
			pcim* temp = nodegenerator(var, 0, 0, contexts, ds);
			return temp;
		}		
	}

	else{
		if(num == 0){//no parent
			pcim* temp = new pcim(new vartest(var), matrixgenerator(var, 0, ds),transfer(NumofNodes, count+1));	
			return temp;
}	
		else{
			pcim* temp = new pcim(new vartest(var), nodegenerator(var, 0, 0, contexts, ds), transfer(NumofNodes, count+1));
			return temp;
		}
	}
}


