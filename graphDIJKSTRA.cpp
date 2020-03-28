//              **HOMEWORK 2: DIJKSTRA'S ALGORITHM IMPLEMENTATION**
//                                by ADAM WOLF
//                 Note: Please compile in a C++11 enabled compiler.

#include <exception>
#include <iostream>
#include <limits>
#include <queue>
#include <vector>
#include <random>
#include <tuple>
#include <ctime>

using namespace std;

class Graph
{
public:
    ///Initialize a graph with v vertices, and use a 1D array with v*v locations for the adjacency matrix.
    ///@param v number of vertices:
    Graph(unsigned int v) : v_(v), e_(0), adj_(v * v, 0) {}

    ///@return the number of vertices.
    unsigned int V() const { return v_; }

    ///@return the number of edges in the graph.
    unsigned int E() const { return e_; }

    ///@return true IF there exists an edge from node-x to node-y.
    bool adjacent(unsigned int x, unsigned int y) const
    {
        return adj_[index_for(x, y)] > 0;
    }

    ///@return list of all nodes y that are connected to node-x.
    vector<unsigned int> neighbors(unsigned int x) const
    {
        vector<unsigned int> result;
        for (int i = 0; i < v_; ++i)
            if (adjacent(x, i))
                result.push_back(i);

        return result;
    }

    ///Add edge from node-x to node-y.
    void add_edge(unsigned int x, unsigned int y, double weight)
    {
        if (weight < 0)
            throw invalid_argument("Edge cost cannot be negative!");
        if (adj_[index_for(x, y)] == 0)
            ++e_;

        adj_[index_for(x, y)] = weight;
        adj_[index_for(y, x)] = weight;
    }

    ///Remove edge from node-x to node-y.
    void remove_edge(unsigned int x, unsigned int y)
    {
        adj_[index_for(x, y)] = 0;
        adj_[index_for(y, x)] = 0;
        --e_;
    }

    ///@return the edge values from x -> y.
    double get_edge_value(unsigned int x, unsigned int y)
    {
        return adj_[index_for(x, y)];
    }

protected:
    ///@return Convert the index of 2D-array to a 1D-array.
    unsigned int index_for(unsigned int x, unsigned int y) const
    {
        if (x >= v_ || y >= v_)
            throw invalid_argument("Index is out of bounds!");

        return x * v_ + y;
    }

private:
    const unsigned int v_;
    unsigned int e_;
    vector<double> adj_;
};

class ShortestPath
{
public:
    ///Here we compute the shortest path from the node source for the graph "g" to all remaining nodes.
    ///@param g the graph @param source the starting node
    ShortestPath(Graph g, unsigned int source = 0) : g_(g),
                                                     dist_(g.V(), numeric_limits<double>::max()),
                                                     source_(source)
    {
        compute();
    }

    ///@return Here we obtain the distance from the source node to the node "idx"
    double operator[](int idx) const
    {
        return dist_[idx];
    }

protected:
    ///Now we can compute the shortest path from the given start node to all remaining nodes with the Dijkstra's algorithm.
    void compute()
    {
        ///Queue for next node to visit each element is a "tuple" of the distance to the node and the node's ID. 
        ///This priority queue is sorted as minimum distance at the top.
        priority_queue<DistNode, vector<DistNode>, less<DistNode>> nodes;

        ///Start the queue with the source node, which by definition has a distance to itself of 0.
        nodes.push(make_tuple(0, source_));

        while (!nodes.empty())
        {
            ///Next, we extract the next node and its distance from the queue.
            double d;
            unsigned int n;
            tie(d, n) = nodes.top();
            nodes.pop();

            ///Here we update the distance array to determine the next nodes to hit.
            if (d < dist_[n])
                dist_[n] = d;

            for (const auto &w : g_.neighbors(n))
            {
                auto new_dist = dist_[n] + g_.get_edge_value(n, w);
                if (new_dist < dist_[w])
                {
                    dist_[w] = new_dist;
                    nodes.push(make_tuple(new_dist, w));
                }
            }
        }
    }

private:
    Graph g_;
    ///Final distance from the source node.
    vector<double> dist_;
    unsigned int source_;

    ///Tuple type for the distance from the source node to a node and its corresponding ID.
    using DistNode = tuple<double, unsigned int>;
};

class Simulator
{
public:
    ///Create a simulator with the given parameters:
    ///@param num gives number of vertices in a graph.
    ///@param density of the generated graph.
    ///@param distance_min
    ///@param distance_max
    ///@param times for number of times to run the simulation.
    Simulator(int num = 50, double density = 0.2,
              double distance_min = 1, double distance_max = 10,
              int times = 50) : num_(num),
                                den_(density),
                                random_generator_(time(NULL)),
                                distance_distribution_(distance_min, distance_max),
                                existence_distribution_(0.0, 1.0)
    {
        double sum_avg_dist = 0;
        for (int i = 0; i < times; ++i)
        {
            sum_avg_dist += run_simulation();
        }
        avg_ = sum_avg_dist / times;
    }

    ///@return the average distance for the shortest path.
    double average_distance() const
    {
        return avg_;
    }

protected:
    ///Generate a graph and run the simulation. 
    ///@return the average shortest path of the generated graph
    double run_simulation()
    {
        Graph g = generate_graph();
        ShortestPath sp(g, 0);
        int count = 0;
        double sum = 0;

        for (int i = 1; i < g.V(); ++i)
            if (sp[i] < numeric_limits<double>::max())
            {
                ++count;
                sum += sp[i];
            }

        return sum / count;
    }
    
    ///@return graph generated with the given parameter.
    Graph generate_graph() 
    {
        Graph g(num_);
        for (int i = 0; i < num_ - 1; ++i)
            for (int j = i + 1; j < num_; ++j)
                if (existence_distribution_(random_generator_) < den_)
                    g.add_edge(i, j, distance_distribution_(random_generator_));

        return g;
    }

private:
    double avg_;
    int num_;
    double den_;

    ///Mersenne Twister 19937 generator: a pseudo-random generator of 32-bit numbers with a state size of 19937 bits.
    mt19937 random_generator_; 

    ///Random number distribution that produces floating-point values according to a uniform distribution. 
    ///This distribution (also know as rectangular distribution) produces random numbers in a range [a,b) 
    ///where all intervals of the same length within it are equally probable. The distribution parameters, 
    ///a and b, are set on construction. To produce a random value following this distribution, call its member function operator().*/
    uniform_real_distribution<double> distance_distribution_;
    uniform_real_distribution<double> existence_distribution_;
};

int main()
{
    cout << "Executing..." << endl;

    Simulator sim20;
    cout << "20% Density Graph: average distance is "
         << sim20.average_distance() << endl;

    Simulator sim40(50, 0.4);
    cout << "40% Density Graph: average distance is "
         << sim40.average_distance() << endl;

    return 0;
}
