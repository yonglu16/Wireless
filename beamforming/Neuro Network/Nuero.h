class inputNode:public Nueron
{
public:
	double AoA;
	double AoD;
	double dis;
	input_triple():AoA(0),AoD(0),dis(0){}
	input_triple(double a,double b, double c):AoA(a),AoD(b),dis(c){}
	~input_triple(){}
	virtual activate();
};
class hiddenNode:public Neuron
{
	
};
class outputNode:public Nueron
{
	
};
class Edge
{
public:
	Edge():target(0),weight(0);
	~Edge();
	tail *Nueron;
	double weight;
}£»
class Nueron
{
public:
	Nueron(){}
	~Nueron(){}
	virtual activate();
	vector<double> ;
};
class Network
{
public:
	Network();
	~Network)();
	vector<Nueron> input_layer;
	vector<Nueron> hidden_layer;
	vector<Nueron> output_layer;
};