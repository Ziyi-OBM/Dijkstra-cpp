/*

06/28/2019
Bug Fixed:
*Fixed the problem that the shortest path algorithm only works when starting at node 0


06/12/2019
Summary:
This C++ programm implements the Dijkstra's algorithm, and applies it on a set of randomly generated maps.

Output:
	Randomly generating graph of density: 0.2 and distance range: 1 ~ 10
	Randomly generating graph of density: 0.4 and distance range: 1 ~ 10
	Average distance of a map with density 20% and distance range 1~10 is 7.03116
	Average distance of a map with density 40% and distance range 1~10 is 4.74452

The program used 3 classes and 2 helper class to implement the algorithem.
Graph class: To represent a graph in matrix representation. Also provides a method for the random generation of a graph.
priorityQueue class: Used minHeap to implement a priority queue, which will be used in the Dijkstra implementation.
shortestPath class: Takes a pointer to a Graph object as input and provides a method that finds the shortest path between two points using Dijkstra's algorithm.

Since the distance between nodes are 1~10, the results suggest that it usually takes only 1~2 steps to reach any points in the graph.
This is expected because such randomly generated map has a high level of symmetry that each node is essentially the same and has 20%/40% chance of directly
connecting to a arbitary node in the map. 
Thus, it is clear that a graph of density 20%~40% with 50 nodes is already a reasonable dense map. 

*/



#include <iostream>
#include <vector>
#include <algorithm>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */


	
//A class used to wrap the information of a edge, i.e. starting node, ending node and cost.		
class edgeQ{
	public:
	
	edgeQ():cost(0.0),start(0),ends(0){}
	edgeQ(double cost, int start, int ends):cost(cost),start(start),ends(ends){}
	
	double cost=0.0;
	int start=0;
	int ends=0;
	
	~edgeQ(){}
};


//A class used to wrap the information of a shortest path solution, i.e. total cost and the predecessors
class pathData{
	public:
	
	double cost=0.0;
	std::vector<int> route;
	
	pathData():cost(-1.0),route(){}

	~pathData(){}
};

std::ostream& operator<< (std::ostream& out, const pathData& path){
	if (path.cost!=-1){
		std::cout << "The cost of the shortest path is: "<< path.cost <<std::endl;
		
		std::cout << "The route is: "  << std::endl;
		std::vector<int> toPrint = path.route;
		
		//https://stackoverflow.com/questions/18446946/printing-vector-in-reverse-order/54448537
		for(auto it = toPrint.crbegin(); it != toPrint.crend(); it++){
			std::cout<<" -> "<<*it;
		}
		std::cout << std::endl;
	}else if(path.cost==-1)
		std::cout << "No path exists."  << std::endl;
	return out;
}


std::ostream& operator<< (std::ostream& out, const edgeQ& edge){
	out << "Edge Object:\nCost: "<<edge.cost<<"; Start at: "<< edge.start<<"; Ends at: "<<edge.ends;
	return out;}

class Graph{
		

	inline double prob(){return (rand() / ((double)RAND_MAX+1));}
	template<typename T>
	inline void disp(T vec2D) const {	//Print vector borrowed from stack overflow
		for(auto vec : vec2D){
			for(auto x : vec)
				std::cout<<x<<" , ";
			std::cout << std::endl;
		}
	}
	template<typename T>
	inline int count(T vec2D) const {
		int sum = 0;
		for(auto vec : vec2D)
		{
			for(auto x : vec)
				if(x!=0.0)++sum;
		}
		return sum;
	}

		
	public:
	//Parameters:
	const int V;	//# of vertices
	
	//Default constructor
	Graph():V(1){std::cout << "Default constructor called" << std::endl;}

	//Parameterized constructor
	Graph(int ver):V(ver){}

	friend std::ostream& operator<< (std::ostream& out, const Graph& g);
	
	//Function to randomly initialize the graph 
	void generateRanomGraph(double edgeDensity, double distRangeLb, double distRangeUb) { 
		std::cout << "Randomly generating graph of density: " << edgeDensity << " and distance range: " << distRangeLb <<" ~ "<< distRangeUb <<std::endl; 
		srand(clock());   
		// srand(20190612);  Setting a constant seed for debugging
		if (distRangeLb>distRangeUb){
			std::cout << "Invalidd Distance Bounds" << std::endl;
			return;
		}
		double range = distRangeUb - distRangeLb;
		  for (int i=0;i<V;++i){
			for (int j=i;j<V;++j){
					if(i==j)G[i][j]=false;
						else{
							G[i][j]=G[j][i]=(prob()<edgeDensity);
							//if(G[i][j])E[i][j]=E[j][i]=(distRange*(1-prob()));
							if(G[i][j])E[i][j]=E[j][i]=(range*prob()+distRangeLb);
						}
			}
		  }	  
		
		//Print the resulting array to debug
		//disp(G);
		//std::cout << "\n" << std::endl;
		//disp(E);	   
	} 

	//Interface methods
	void printNodes() const {disp(G);}	//Print the Nodes as a matrix

	void printEdges() const {disp(E);}	//Print the Edges as a matrix

	int vertices() const {return V;} //Returns the number of vertices in the graph

	int edges() const {return count(E);}	//Returns the number of edges in the graph

	bool adjacent(int i, int j) const {return G[i][j];}	//Tests whether there is an edge from node i to node j.

	std::vector<int> neighbors(int i) const {	//Print all nodes j such that there is an edge from i to j.
		std::vector<int> neighborNodes;
		for (int j=0;j<V;++j){
			if(G[i][j]){
				neighborNodes.push_back(j);
			}
		}
		return neighborNodes;
	}

	void addEdge (int i, int j, double cost){	//adds to G the edge from i to j with cost x, if it is not there.
			if(!G[i][j]){
				G[i][j]=G[j][i]=true;
				E[i][j]=E[j][i]=cost;
			}else{std::cout<<"Error: Edge exists. The value is unchanged. Use setEdgeValue to assign new value"<<std::endl;}
	}

	void deleteEdge (int i, int j){	//removes the edge from i to j, if it exists.
			if(G[i][j]){
				G[i][j]=G[j][i]=false;
				E[i][j]=E[j][i]=0.0;
			}else{std::cout<<"Warning: Edge does not exists. No deletion peformed"<<std::endl;}
	}

	bool getNodeValue(int i, int j) const {return G[i][j];}	//returns the value associated to the node (i,j) (True/False=there Is/Not a node).

	// set_node_value is not necessary because its the same as the "addEdge" method

	double getEdgeValue(int i, int j) const {return E[i][j];}	//returns the value associated to the edge (i,j).

	double setEdgeValue(int i, int j, double cost){	//sets the value associated to the edge (i,j) to v.
		//If the two nodes are not connected, remind the user to use "addEdge" to add a new connection between the two nodes
		if(!G[i][j])std::cout<<"Error: No edge exists between given nodes. please use the 'addEdge' method"<<std::endl;
						else{E[i][j]=E[j][i]=cost;}
	}


	//Destructor
	~Graph() {} 


	//Private variables
	private:


	//Initialize a VxV bool matrix to represent the nodes
	std::vector<bool> rowg = std::vector<bool>(V);
	std::vector<std::vector<bool> >   G = std::vector< std::vector<bool> >(V,rowg);

	//Initialize a VxV double matrix to store the cost of the edges
	std::vector<double> rowe = std::vector<double>(V,0.0);
	std::vector<std::vector<double> >  E = std::vector< std::vector<double> >(V,rowe);

		
};

std::ostream& operator<< (std::ostream& out, const Graph& g){
	std::cout << "Graph Object: "  << std::endl;
	std::cout << "Node Matrix: "  << std::endl;
	for(auto vec : g.G){
			for(auto x : vec)
				std::cout<<x<<" , ";
			std::cout << std::endl;
		}
		
	std::cout << "Edge Matrix: "  << std::endl;
	for(auto vec : g.E){
			for(auto x : vec)
				std::cout<<x<<" , ";
			std::cout << std::endl;
		}	
		
	return out;
}
	
	
	
	
class priorityQueue{
	//Implement priority queue using min heap	
	//Reference: https://www.geeksforgeeks.org/binary-heap/
	template<typename T>
	inline void swap(T *a, T *b){	//Print vector borrowed from stack overflow
		T temp = *a;
		*a=*b;
		*b=temp;
	}	

	inline int parentId(int i) { return (i-1)/2; } 
	inline int LchildId(int i) { return (2*i)+1; } 
	inline int RchildId(int i) { return (2*i)+2; } 
	inline int smaller(int i,int j) const {
		if (minHeap[i].cost<minHeap[j].cost){
			return i;
		}else{
			return j;
		}
		
	}

	//Defining member variables
	int capacity;
	edgeQ edgeInit=edgeQ();
	std::vector<edgeQ> minHeap;
	int heapEnd=0;	
	
	public: 	
	//Default constructor: Create an empty queue of max size 5000
	priorityQueue():capacity(5000), minHeap(capacity,edgeInit) {}
	//Constructor: Create an empty queue of max size cap
	priorityQueue(int cap):capacity(cap), minHeap(capacity,edgeInit) {}
	
	void getCapacity() const { // Print capacity
				std::cout << "The capacity of the Queue is: " << capacity <<std::endl; 
	}
	
	//Test if the queue contains a path to a certain destination
	//Returns the index of the edge in the queue, or return -1 if it does not contain such path
	int contains(int destination) const{
		int itIndex=0;
		for (auto it : minHeap){
			if (it.ends == destination){
				return itIndex;
			}
			itIndex++;
		}
		return -1;
	}
	
	//Increase the priority(decrease the cost) of a edge in the queue which has index: nodeIndex in the minHeap vector, 
	//also update the start point of this new edge/path 
	void increasePriority (int nodeIndex, int newStart, double newCost){
		if (nodeIndex==-1){
			return;
		}
		if (minHeap[nodeIndex].cost<newCost){
			//std::cout<<"The updated route does not improve from the existing one."<<std::endl;
			return;
		}		
		minHeap[nodeIndex].start=newStart;
		minHeap[nodeIndex].cost=newCost;		
		moveUp(nodeIndex);
		return;
	}
		
	
	void insert(edgeQ edgeIn){ // insert queue_element into queue
		if (heapEnd == capacity){ 
				std::cout << "\nInsertion failed: Queue reached max capacity\n" << std::endl; 
				return;
			} 
		minHeap[heapEnd]=edgeIn;	

		moveUp(heapEnd);
		heapEnd++; 
	}
	
	edgeQ extract(){ // Return and delete the root element
		if (heapEnd == 0){ 
				std::cout << "\nExtract failed: Queue is empty\n" << std::endl; 
				return edgeQ();
			} 
			
		if (heapEnd == 1){ 
				heapEnd--;
				return minHeap[0];
			} 

		edgeQ root = minHeap[0];

		
		heapEnd--;
		minHeap[0]=minHeap[heapEnd];

		
		moveDown(0);
		
		return root;
	}

	edgeQ top() const { // Return the root element
		if (heapEnd == 0){ 
				std::cout << "\nExtract failed: Queue is empty\n" << std::endl; 
				return edgeQ();
			} 
		
		return minHeap[0];
	}
	
	
	int sizeQ() const { // Return the root element
		std::cout << "The size is: " << heapEnd << "\n" << std::endl; 
		return heapEnd;
	}
	
	bool isEmpty() const { // Return the root element
		if (heapEnd == 0){ 
			return true;
		}else{
			return false;
		}
	}
	
	void clearAll(){
		heapEnd=0;		
	}

	//Destructor
	~priorityQueue() {} 


	private:

	//Compare a node with its parents and move it upward to maintain heap property
	void moveUp(int cursor){
	int parent = parentId(cursor);
		while(cursor!=0 & minHeap[parent].cost>minHeap[cursor].cost){
			
			swap(&minHeap[parent],&minHeap[cursor]);
			
			cursor = parent;
			parent = parentId(cursor);

		}
	}
	
	//Compare a node with its childs and move it downward to maintain heap property
	void moveDown(int cursor){
		
		if (LchildId(cursor)>=heapEnd){
			return;
		}
		
		int smallerChild= LchildId(cursor);
		
		if (RchildId(cursor) < heapEnd & minHeap[smallerChild].cost > minHeap[RchildId(cursor)].cost){
			smallerChild = RchildId(cursor);
		}
		
		if ( minHeap[cursor].cost > minHeap[smallerChild].cost){
			swap(&minHeap[cursor],&minHeap[smallerChild]);	
			
			cursor=smallerChild;
			moveDown(cursor);
		}
		
	}

};


//Class to find the shortest path between two verties, given a graph
class shortestPath{
	
	Graph const *G;
	
	int vertices;
	public:
	
	//Initialize the object
	shortestPath(Graph const *graphIn):G(graphIn),vertices(G->vertices()){
		//Debugging: Printing the input information
		//std::cout << "Number of Nodes in the Graph: " << vertices <<std::endl;
		//std::cout << "Number of Edges in the Graph: " << G->edges() <<std::endl;
		//Q.getCapacity();
		//std::cout << *G <<std::endl;
	}
	
	//Method to find the shortest path
	pathData findPath(int u, int w){
		//Create a priority queue that serves as the openset discussed in the algorithm
		priorityQueue Q(G->edges());
		
		
		std::vector<int> neighborNodes = G->neighbors(u);
		std::vector<double> dist(vertices,-1);	//A vector that holds the shortest distance from the start point to each node
		std::vector<int> prev(vertices,-1);	//A vector that holds the previous nodes on the shortest path
		
		int closedSetSize = 0; 	//The close set is the visited nodes that we've found the shorted distance
								//It is the number of non-"-1" element in the dist and prev vector
		dist[u]=0;
		prev[u]=u;	
		closedSetSize++;
		
		//Add the neibours of the starting node to the openset
		for (auto x : neighborNodes){
			edgeQ edgeIn=edgeQ(G->getEdgeValue(u,x),u,x);
			if (dist[edgeIn.ends]==-1)
				Q.insert(edgeIn);
		//Debugging: Printing elements that get added to the openset.
		//std::cout << "Inserting: " << std::endl;
		//std::cout << edgeIn << std::endl;			
			
		}
		
		//While we have not reach destination and the openset is not exhausted, carry on the algorithm
		while (dist[w] == -1 & !Q.isEmpty()){
			
			//Extract the highest priority element in the openset and move it to the closed set
			edgeQ currentEdge=Q.extract();
			//Debugging: Printing elements that get added to the close set.
			//std::cout << "Evaluating: " << std::endl;
			//std::cout << currentEdge << std::endl;
		
			//Save the distance and path to the dist and prev vector
			dist[currentEdge.ends]=currentEdge.cost;
			closedSetSize++;
			prev[currentEdge.ends]=currentEdge.start;
			
			//Explore the neibours of the point that was added to the closed set
			neighborNodes = G->neighbors(currentEdge.ends);
			for (auto x : neighborNodes){
				//Gather the information of the edges to its neibour
				edgeQ edgeIn(G->getEdgeValue(currentEdge.ends,x)+dist[currentEdge.ends],currentEdge.ends,x);
				
				//Update the openset base on these neibours
				int isInQ = Q.contains(edgeIn.ends);
				if (isInQ==-1 & dist[edgeIn.ends]==-1) Q.insert(edgeIn);			
				//This method will test whether it improves over the previous path, and update it when it does			
				else if (isInQ!=-1 & dist[edgeIn.ends]==-1) Q.increasePriority (isInQ, edgeIn.start, edgeIn.cost);
				
			}			
		}
		
		std::vector<int> route;	//A variable to hold the vertices along the shortest path
		
		//After the growing process
		//First test if we have reached the destination
		if (dist[w]!=-1){
			//If true, wrap the result to a helper class
			int cursor = w;
			route.push_back(w);
			while (cursor != u){
				route.push_back(prev[cursor]);
				cursor=prev[cursor];
			}
			double totalCost = dist[w];	
			pathData pathOut;
			pathOut.route = route;
			pathOut.cost = totalCost;
			return pathOut;
			//If destination is not reached, return a empty path with cost -1
		}else if (Q.isEmpty()) return pathData();
		
	}	
	
	//Destructor
	~shortestPath() {} 
		
	private:

	
	
};


int main(){
		
	// //Testing the Graph class	
	// std::cout << "\n\n Testing the Graph class" << std::endl;
	// Graph testMap(5);
	// testMap.generateRanomGraph(0.5, 0.0, 10.0);
	// std::cout << "Printing the generated map: " << std::endl;	
	// std::cout << testMap << std::endl;	


	// std::cout << "Number of vertices = " << testMap.vertices() << std::endl;	
	// std::cout << "Number of edges = " << testMap.edges() << std::endl;	
	// std::cout << "Node 3 is connected to node " << std::endl;	
	// //std::cout << testMap.neighbors(3) << std::endl;

	// for(auto  x : testMap.neighbors(3)){
		// std::cout << x <<" , ";
	// }
	// std::cout << std::endl;		

	// testMap.addEdge (4, 4, 5.5);
	// std::cout << "The connectivity matrix" << std::endl;	
	// testMap.printNodes();

	// std::cout << "The edge cost" << std::endl;
	// testMap.printEdges();

	// std::cout << "The edge cost at (4,4) is " << testMap.getEdgeValue(4, 4) << std::endl;	

	// //Test the priorityQueue class
	// std::cout << "\n\n Testing the priorityQueue class" << std::endl;
	// double A[5]={5,4,3,7,1};

	// priorityQueue testQ(100);
	// for (int i=0; i<5; i++){
		// std::cout << "\nInserting " << A[i] << ": " << std::endl;
		// edgeQ edgeInsert = edgeQ(A[i],0,i);
		// testQ.insert(edgeInsert);	
		// edgeQ printTop=testQ.top();
		// std::cout << "Top element is: " << std::endl; 	
		// std::cout << printTop << std::endl;
		// //std::cout << printTop.cost <<"   "<< printTop.start <<"   "<< printTop.ends <<std::endl;
	// }

	// testQ.sizeQ();

	// while (!testQ.isEmpty()){
		// edgeQ root=testQ.extract();
		// std::cout << "\nExtraction: " << std::endl;
		// std::cout << root << std::endl;
		// //std::cout << root.cost <<"   "<< root.start <<"   "<< root.ends <<std::endl;
	// }

	
	// //Testing shortest pass
	// //https://www.geeksforgeeks.org/dijkstras-shortest-path-algorithm-greedy-algo-7/
	// Graph testDijk(9);
	// testDijk.addEdge(0,1,4);
	// testDijk.addEdge(0,7,8);
	// testDijk.addEdge(1,7,11);
	// testDijk.addEdge(1,2,8);
	// testDijk.addEdge(2,3,7);
	// testDijk.addEdge(2,5,4);
	// testDijk.addEdge(2,8,2);
	// testDijk.addEdge(3,4,9);
	// testDijk.addEdge(3,5,14);
	// testDijk.addEdge(4,5,10);
	// testDijk.addEdge(5,6,2);
	// testDijk.addEdge(6,7,1);
	// testDijk.addEdge(6,8,6);
	// testDijk.addEdge(7,8,7);
	
	// std::cout << testDijk << std::endl;	
	
	// shortestPath finder(&testDijk);
	
	// pathData solution = finder.findPath(0,4);
	// std::cout<<solution<<std::endl;
	
	
	
	
	
	
	//Run random graph simulation
	
	//First trial: Density = 20%
	std::vector<double> trialLengh(49);

	Graph map1(50);
	map1.generateRanomGraph(0.2,1,10);
	shortestPath pathMap1(&map1);
	for (int j=1;j<50;j++){
		pathData solutionMap1 = pathMap1.findPath(0,j);
		trialLengh[j-1]=solutionMap1.cost;
	}
	//accumulate( trialLengh.begin(), trialLengh.end(), 0.0)/trialLengh.size(); 	Not good because it induces error when there is disconnection
	double avgLength=0;
	double runningSum=0;
	int count=0;
	for (int i=0;i<50;i++){
		if (trialLengh[i]!=-1){
			runningSum+=trialLengh[i];
			count++;
		}
	}
	avgLength=runningSum/count;
	
	
	
	//Second trial: Density = 40%
	std::vector<double> trialLengh2(49);
	
	Graph map2(50);
	map2.generateRanomGraph(0.4,1,10);
	shortestPath pathMap2(&map2);
	for (int j=1;j<50;j++){
		pathData solutionMap1 = pathMap2.findPath(0,j);
		trialLengh2[j-1]=solutionMap1.cost;
	}
	double avgLength2=0;
	double runningSum2=0;
	int count2=0;
	for (int i=0;i<50;i++){
		if (trialLengh2[i]!=-1){
			runningSum2+=trialLengh2[i];
			count2++;
		}
	}
	

	avgLength2=runningSum2/count2;
	
	
	std::cout<<"Average distance of a map with density 20% and distance range 1~10 is "<< avgLength <<std::endl;
	//Debugging: Printing the shortest path length
	//for (auto x:trialLengh){
	//	std::cout<<x<<", ";
	// }
	
	std::cout<<"Average distance of a map with density 40% and distance range 1~10 is "<< avgLength2 <<std::endl;

	//Debugging: Printing the shortest path length
	//for (auto x:trialLengh2){
	//	std::cout<<x<<", ";
	// }
	
	
	

	return 0;	
}