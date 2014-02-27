#include "transfer.h"

using namespace std;

pcim* rwogenerator(istream& is, int var, int state, int ds){
	double rate;
	is>>rate;
	if(state == ds-2){
		double last;
		is>>last;
		pcim *temp = new pcim(new eventtest(var, state), new pcim(rate<0? 0:rate), new pcim(last<0? 0: last));
		return temp;
	}
	else{
		pcim *temp = new pcim(new eventtest(var, state), new pcim(rate<0? 0:rate), rwogenerator(is, var, state+1, ds));
		return temp;
	}		

}

pcim* matrixgenerator(istream& is, int var, int state, int ds){
		
	if(state == ds-2){
		pcim* temp = new pcim(new lasttest(var, state), rwogenerator(is, var, 0, ds), rwogenerator(is, var, 0, ds));
		return temp;
	}	
	else{	
		pcim* temp = new pcim(new lasttest(var, state), rwogenerator(is, var, 0, ds), matrixgenerator(is, var, state+1, ds));
		return temp;
	}
}

pcim* nodegenerator(istream& is, int var, int parentindex, int state, ctbn::Context &contexts, int ds){//only when var has parents
	int id = contexts.VarList()[parentindex];
	if(parentindex == contexts.VarList().size()-1 && (state == contexts.Cardinality(id)-1)){
		pcim *temp = matrixgenerator(is, var, 0, ds);
		return temp;
	}	

	if(parentindex == contexts.VarList().size()-1){
		pcim *temp = new pcim(new lasttest(id, state), matrixgenerator(is, var, 0, ds), nodegenerator(is, var, parentindex, state+1, contexts, ds));
		return temp;
	}
		
	pcim *temp = new pcim(new lasttest(id, state), nodegenerator(is, var, parentindex+1, 0, contexts, ds), nodegenerator(is, var, parentindex, state+1, contexts, ds));
	return temp;
	
}

pcim* graphgenerator(istream& is, int NumofNodes, int count = 1){
	
	ctbn::Context contexts; //parent context

	int var;
	is>>var;
	int ds;
	is>>ds;
	int num;
	is>>num;	

	for(int i=0; i<num; i++){
		int id;
		int card;
		is>>id;
		is>>card;
		contexts.AddVar(id, card);
	}

	if(count == NumofNodes){
		if(num == 0){//no parent
			pcim* temp = matrixgenerator(is, var, 0, ds);	
			return temp;		
		}	
		else{
			pcim* temp = nodegenerator(is, var, 0, 0, contexts, ds);
			return temp;
		}		
	}

	else{
		if(num == 0){//no parent
			pcim* temp = new pcim(new vartest(var), matrixgenerator(is, var, 0, ds),graphgenerator(is, NumofNodes, count+1));	
			return temp;
}	
		else{
			pcim* temp = new pcim(new vartest(var), nodegenerator(is, var, 0, 0, contexts, ds), graphgenerator(is, NumofNodes, count+1));
			return temp;
		}
	}
}

pcim* CTBNtransfer(istream& is){
	int NumofNodes;
	is>>NumofNodes;
	return graphgenerator(is, NumofNodes);
}


