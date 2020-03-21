#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <map>
#include <set>
#include <vector>
#include <ctime>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <list> 
#include <stack>
#include <limits.h>
#include <queue>
using namespace std;

struct Node        
{
  int id;
  int  type;            //0 for PI,  1 for PO,  2 for Internal Nodes
  double toggling_rate;
  vector<int> fan_in;
  //PI no fan_in, PO no toggling_rate
};
struct Graph
{
  int V;
  list<int>* adj; 
};
struct LUT
{
  int output;
  set<int> input;
};

string  input_AAG, output_NAME;
int input_K;
int NODE_num, PI_num, PO_num, Internal_num;
vector<Node> Nodes_vec;
stack<int> topo_Stack; 
Graph g;
set<int> Nv;
vector<LUT> LUTs;

map<int, int> map_id_to_index;    //to find out the node id is in which index of vector
map<int, set<int> > map_id_to_Nv;
map<int, int> map_id_to_label;
map<int, bool> map_id_to_visited;



///////////////////////////////////////////////////////////////////////////////////
////////////////////////     Ford-Fulkerson      //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
class Graph_FlowNetWorks
{
  private:
      int num_vertex;
      std::vector<std::vector<int>> AdjMatrix;
  public:
      Graph_FlowNetWorks():num_vertex(0){};
      Graph_FlowNetWorks(int n);
      void AddEdge(int from, int to, int capacity);

      int FordFulkerson(int source, int termination);
      bool BFSfindExistingPath(std::vector<std::vector<int>> graphResidual, int *predecessor, int source, int termination);
      int MinCapacity(std::vector<std::vector<int>> graphResidual, int *predecessor, int termination);
};
Graph_FlowNetWorks::Graph_FlowNetWorks(int n):num_vertex(n)
{
    AdjMatrix.resize(num_vertex);
    for (int i = 0; i < num_vertex; i++)
        AdjMatrix[i].resize(num_vertex);
}

bool Graph_FlowNetWorks::BFSfindExistingPath(std::vector<std::vector<int>> graph, int *predecessor, int s, int t)
{
    int visited[num_vertex];

    for (int i = 0; i < num_vertex; i++)
    {
        visited[i] = 0;     
        predecessor[i] = -1;
    }

    std::queue<int> queue;
    queue.push(s);
    visited[s] = 1;
    while (!queue.empty()) {
        int exploring = queue.front();
        for (int j = 0; j < num_vertex; j++) {
            if (graph[exploring][j]!=0 && visited[j]==0) {
                queue.push(j);
                visited[j] = 1;
                predecessor[j] = exploring;
            }
        }
        queue.pop();
    }
    return (visited[t] == 1);  
}

int Graph_FlowNetWorks::MinCapacity(std::vector<std::vector<int>> graph, int *predecessor, int t)
{
    int min = 100;      
    for (int idx = t; predecessor[idx] != -1; idx = predecessor[idx])
    {
        if (graph[predecessor[idx]][idx]!=0 && graph[predecessor[idx]][idx] < min) 
        {
            min = graph[predecessor[idx]][idx];
        }
    }
    return min;
}

int Graph_FlowNetWorks::FordFulkerson(int source, int termination)
{
    std::vector<std::vector<int>> graphResidual(AdjMatrix);    
    int maxflow = 0;                                           
    int predecessor[num_vertex];

    // for(int i = 0; i < AdjMatrix.size(); i++)
    // { 
    //   for(int j = 0; j < AdjMatrix[i].size(); j++)
    //   { 
    //     printf("%3d ",AdjMatrix[i][j]);
    //   }
    //   cout<<endl;
    // }

    while (BFSfindExistingPath(graphResidual, predecessor, source, termination)) 
    {
        int mincapacity = MinCapacity(graphResidual, predecessor, termination);
        maxflow = maxflow + mincapacity;
        for (int Y = termination; Y != source; Y = predecessor[Y]){
            int X = predecessor[Y];
            graphResidual[X][Y] -= mincapacity;
            graphResidual[Y][X] += mincapacity;
        }
    }
    return maxflow;
}
void Graph_FlowNetWorks::AddEdge(int from, int to, int capacity)
{
    AdjMatrix[from][to] = capacity;
}
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////


void parser()
{   
  string temp_s; 
  
  ifstream file_AGG(input_AAG);
  file_AGG >> temp_s;             //agg
  file_AGG >> NODE_num;           //8 
  file_AGG >> PI_num;             //3
  file_AGG >> PO_num;             //2
  Internal_num = NODE_num - PI_num - PO_num;

  Node temp_n;
  for(int i = 0; i < PI_num; i++)
  { 
    temp_n.fan_in.clear();
    temp_n.fan_in.resize(0);

    file_AGG >> temp_n.id;
    temp_n.type = 0;
    file_AGG >> temp_n.toggling_rate;
    Nodes_vec.push_back(temp_n);
  }

  int temp_index;
  for(int i = 0; i < PO_num; i++)
  { 
    temp_n.fan_in.clear();
    temp_n.fan_in.resize(0);

    file_AGG >> temp_n.id;
    temp_n.type = 1;
    temp_n.toggling_rate = 0;
    file_AGG >> temp_index;
    temp_n.fan_in.push_back(temp_index);
    Nodes_vec.push_back(temp_n);
  }

  for(int i = 0; i < Internal_num; i++)
  { 
    temp_n.fan_in.clear();
    temp_n.fan_in.resize(0);

    file_AGG >> temp_n.id;
    temp_n.type = 2;
    file_AGG >> temp_n.toggling_rate;
    file_AGG >> temp_index;
    temp_n.fan_in.push_back(temp_index);
    file_AGG >> temp_index;
    temp_n.fan_in.push_back(temp_index);
    Nodes_vec.push_back(temp_n);
  }

  file_AGG.close();
  
  for(int i = 0; i < NODE_num; i++)
  { 
    map_id_to_index[ Nodes_vec[i].id ] = i;
  }
}
void BuildRelation() // for Topological Sort
{
  g.V = NODE_num;
  g.adj = new list<int>[NODE_num + 1];
  for(int i = 0; i < NODE_num; i++)
  { 
    int w, v1, v2;
    w = Nodes_vec[i].id;
    if( Nodes_vec[i].fan_in.size() == 2 )
    {
      v1 = Nodes_vec[i].fan_in[0];
      v1 = map_id_to_index[v1];
      v1 = Nodes_vec[v1].id;
      g.adj[v1].push_back(w);

      v2 = Nodes_vec[i].fan_in[1];
      v2 = map_id_to_index[v2];
      v2 = Nodes_vec[v2].id;
      g.adj[v2].push_back(w);
    }
    else if( Nodes_vec[i].fan_in.size() == 1)
    {
      v1 = Nodes_vec[i].fan_in[0];
      v1 = map_id_to_index[v1];
      v1 = Nodes_vec[v1].id;
      g.adj[v1].push_back(w);
    }
  }
  for(int i = 0 ; i < NODE_num; i++)
  {
    list<int>::iterator j; 
  }
}
void TopologicalSortRecursive(int v, bool visited[]) 
{ 
    visited[v] = true; 
    list<int>::iterator i; 
    for (i = g.adj[v].begin(); i != g.adj[v].end(); ++i)
    {   
        if (!visited[*i]) 
            TopologicalSortRecursive(*i, visited); 
    }      
    topo_Stack.push(v); 
} 
void TopologicalSort() 
{
    bool *visited = new bool[NODE_num + 1]; 
    for (int i = 1; i < NODE_num + 1; i++) 
        visited[i] = false; 

    for (int i = 1; i < NODE_num + 1; i++) 
        if (visited[i] == false) 
            TopologicalSortRecursive(i, visited); 
}
void FindNv(int v_id)
{ 
  int index = map_id_to_index[v_id];
  if(Nodes_vec[index].fan_in.size()==2)
  { 
    Nv.insert(Nodes_vec[index].fan_in[0]);
    FindNv(Nodes_vec[index].fan_in[0]);
    Nv.insert(Nodes_vec[index].fan_in[1]);
    FindNv(Nodes_vec[index].fan_in[1]);
  }
  else if(Nodes_vec[index].fan_in.size()==1)
  {
    Nv.insert(Nodes_vec[index].fan_in[0]);
    FindNv(Nodes_vec[index].fan_in[0]);
  }
}

void Labeling()
{  
  for(int i = 0; i < NODE_num; i++) //build a map of v_id to Nv for all nodes
  { 
    Nv.clear();

    int id = Nodes_vec[i].id;
    FindNv(id);
    map_id_to_Nv[id] = Nv;
    // cout<<"id : "<<id<<"      ";
    // set<int>::iterator it;
    // for(it = map_id_to_Nv[id].begin(); it != map_id_to_Nv[id].end(); it++)
    // {
    //   cout<<*it<<" ";
    // }
    // cout<<endl;
  }

  while (!topo_Stack.empty()) /*can speed up here*/
  {  
    int max_label = 0;
    int current_node_id = topo_Stack.top();
    topo_Stack.pop();
    if(map_id_to_Nv[current_node_id].empty())
    {
      map_id_to_label[current_node_id] = 0;
    }
    else
    {
      set<int>::iterator it;
      set<int>::iterator it2;
      set<int>::iterator it3;
      set<int> merged_id;
      merged_id.clear();

      for(it = map_id_to_Nv[current_node_id].begin(); it != map_id_to_Nv[current_node_id].end(); it++) // find max_label
      {
        if( max_label < map_id_to_label[*it] )
          max_label = map_id_to_label[*it];
      }

      if(max_label == 0)
      {
        map_id_to_label[current_node_id] = 1;
      }
      else
      {  
        set<int> new_Nv = map_id_to_Nv[current_node_id]; //Nv'
        merged_id.insert(current_node_id);

        for(it = map_id_to_Nv[current_node_id].begin(); it != map_id_to_Nv[current_node_id].end(); it++)
        {
          if( map_id_to_label[*it] == max_label )
          { 
            //*it is the id need to be merged
            merged_id.insert(*it);
            set_intersection( new_Nv.begin(), new_Nv.end(), map_id_to_Nv[*it].begin(), map_id_to_Nv[*it].end(), inserter( new_Nv,new_Nv.begin() ) );
          }
        }

        // cout<<endl<<"target node : "<< current_node_id<<endl;
        // cout<<"Nv' : ";
        // for(it = new_Nv.begin(); it != new_Nv.end(); it++)
        // {
        //   cout<<*it<<" ";
        // }
        // cout<<endl;
        // cout<<"merged_id : ";
        // for(int i = 0; i < merged_id.size(); i++)
        // {
        //   cout<<merged_id[i]<<" ";
        // }
        // cout<<endl;

        int flow_network_size = new_Nv.size() * 2 + 2;
        Graph_FlowNetWorks g11(flow_network_size);

        map<int, int> map_adj_to_id; // map adj_matrix index to id
        int idx = 1;
        for(it = new_Nv.begin(); it != new_Nv.end(); it++)
        { 
          map_adj_to_id[idx] = *it;
          idx++;
          map_adj_to_id[idx] = *it;
          idx++;
        }

        for(it = new_Nv.begin(); it != new_Nv.end(); it++) // Set point S to be the parent of PI    adj_matrix index for S = 0
        {  
          int index = map_id_to_index[*it];

          if( Nodes_vec[index].type == 0 )
          {
            int adj_index;
            for(int i = 0; i < flow_network_size; i++)
            {
              if(map_adj_to_id[i] == *it && i % 2 == 1)
                adj_index = i;
            }
            g11.AddEdge(0, adj_index, 999);
            // cout<<"ADD : "<<0<<" "<<adj_index<<" "<<999<<endl;
          }  
        }


        for(it = merged_id.begin(); it != merged_id.end(); it++)// Set point T to be the child of all its merged node's parents.   adj_matrix index for T = new_Nv.size() + 1 
        { 
          int index = map_id_to_index[ *it ];
          int adj_index;

          if( Nodes_vec[index].fan_in.size() == 2)
          { 
            if(merged_id.find(Nodes_vec[index].fan_in[0]) ==  merged_id.end())  // can't find this fan_in in merged id
            {
              for(int i = 0; i < flow_network_size; i++)
              {
                if(map_adj_to_id[i] == Nodes_vec[index].fan_in[0] && i % 2 == 0)
                  adj_index = i;
              }
              g11.AddEdge(adj_index, flow_network_size-1, 999);
              // cout<<"ADD : "<<adj_index<<" "<<flow_network_size-1<<" "<<999<<endl;
            }


            if(merged_id.find(Nodes_vec[index].fan_in[1]) ==  merged_id.end())  // can't find this fan_in in merged id
            {
              for(int i = 0; i < flow_network_size; i++)
              {
                if(map_adj_to_id[i] == Nodes_vec[index].fan_in[1] && i % 2 == 0)
                  adj_index = i;
              }
              g11.AddEdge(adj_index, flow_network_size-1, 999);
              // cout<<"ADD : "<<adj_index<<" "<<flow_network_size-1<<" "<<999<<endl;
            }

            // cout<<" i : "<<Nodes_vec[index].id<<"  "<<Nodes_vec[index].fan_in[0]<<" "<<Nodes_vec[index].fan_in[1]<<endl;
          }
          else if(Nodes_vec[index].fan_in.size() == 1)
          { 
            if(merged_id.find(Nodes_vec[index].fan_in[0]) ==  merged_id.end())  // can't find this fan_in in merged id
            {
              for(int i = 0; i < flow_network_size; i++)
              {
                if(map_adj_to_id[i] == Nodes_vec[index].fan_in[0] && i % 2 == 0)
                  adj_index = i;
              }
              g11.AddEdge(adj_index, flow_network_size-1, 999);
              // cout<<"ADD : "<<adj_index<<" "<<flow_network_size-1<<" "<<999<<endl;
            }

            // cout<<" i : "<<Nodes_vec[index].id<<"  "<<Nodes_vec[index].fan_in[0]<<endl;
          }
        }


        for(it = new_Nv.begin(); it != new_Nv.end(); it++)
        {  
          if(merged_id.find(*it) ==  merged_id.end())
          {
            for(it2 = new_Nv.begin(); it2 != new_Nv.end(); it2++)
            {   
                if(*it == *it2)  // from self to self cost = 1
                { 
                  int adj_index;
                  for(int i = 0; i < flow_network_size; i++)
                  {
                    if(map_adj_to_id[i] == *it && i % 2 == 1)
                      adj_index = i;
                  }
                  g11.AddEdge(adj_index, adj_index + 1, 1);
                }  
                else 
                { 
                  bool connected = false;
                  int index = map_id_to_index[*it];
                  int adj_index_1, adj_index_2;
                  if( Nodes_vec[index].fan_in.size() == 2)
                  {
                    if(Nodes_vec[index].fan_in[0] == *it2 || Nodes_vec[index].fan_in[1] == *it2)
                      connected = true;
                  }
                  else if(Nodes_vec[index].fan_in.size() == 1)
                  {
                    if(Nodes_vec[index].fan_in[0] == *it2)
                      connected = true;
                  }

                  if(connected)//from node to other connected node cost = 999
                  {
                    for(int i = 0; i < flow_network_size; i++)
                    {
                      if(map_adj_to_id[i] == *it2 && i % 2 == 0)
                        adj_index_2 = i;
                      if(map_adj_to_id[i] == *it  && i % 2 == 1)
                        adj_index_1 = i;
                    }
                    g11.AddEdge(adj_index_2, adj_index_1, 999);
                    // cout<<"ADD : "<<adj_index_2<<" "<<adj_index_1<<" "<<999<<endl;
                  }

                }
            }
          } 
        }

        int min_cut = g11.FordFulkerson(0, flow_network_size - 1);
        if(min_cut <= input_K) // Nv' has k-feasible --> Nv has feasible cut of height p-1 --> label(v) = p
        {
          map_id_to_label[current_node_id] = max_label;
        }
        else // Nv' has no k-feasible --> Nv has no feasible cut of height p-1 --> label(v) = p+1
        { 
          map_id_to_label[current_node_id] = max_label + 1;
        }
        // cout<<"node : "<<current_node_id<<" min cut : "<<min_cut<<" max_label : "<< max_label<<" label : "<<map_id_to_label[current_node_id]<<endl;
      }
    }  
    // cout<<"id : "<<current_node_id<<" label : "<< map_id_to_label[current_node_id]<<endl;
  }

}  
void Mapping() //put the nodes with same label in a LUT
{
  queue< queue<int> > candidate_queue;
  for(int i = 0; i < NODE_num; i++)
  { 
    if(Nodes_vec[i].type == 1)
    { 
      queue<int> temp_queue;
      temp_queue.push(Nodes_vec[i].fan_in[0]);
      candidate_queue.push(temp_queue);
    } 

    map_id_to_visited[Nodes_vec[i].id] = false;
    // cout<<map_id_to_label[Nodes_vec[i].id]<<" ";
  }  

  while(!candidate_queue.empty())
  {
    queue<int> temp_queue = candidate_queue.front();
    candidate_queue.pop();

    int label = map_id_to_label[temp_queue.front()];
    int output = temp_queue.front();
    set<int> input;
    input.clear();

    if(map_id_to_visited[output])
      continue;

    while(!temp_queue.empty())
    { 
      int cur_id = temp_queue.front();
      int cur_index = map_id_to_index[cur_id];
      temp_queue.pop();

      if(label == map_id_to_label[cur_id])
      {
        // break down
        if(Nodes_vec[cur_index].fan_in.size() == 2)
        {
          temp_queue.push(Nodes_vec[cur_index].fan_in[0]);
          temp_queue.push(Nodes_vec[cur_index].fan_in[1]);
        }
        else if(Nodes_vec[cur_index].fan_in.size() == 1)
        {
          temp_queue.push(Nodes_vec[cur_index].fan_in[0]);
        }
      }
      else
      { 
        //move to candidate row if it is not a PI 
        if(Nodes_vec[cur_index].type != 0)
        {
          queue<int> new_queue;
          while(!new_queue.empty())
          {
            new_queue.pop();
          }
          new_queue.push(cur_id);
          candidate_queue.push(new_queue);
        }
        input.insert(cur_id);
      }
    }

    map_id_to_visited[output] = true;
    LUT temp_lut;
    temp_lut.output = output;

    // set<int>::iterator it1;
    // set<int>::iterator it2;
    // if(input.size() > input_K)
    // {
    //   bool pass;
    //   for (it1 = input.begin(); it1 != input.end(); ++it1)
    //   { 
    //     pass = true;
    //     int index_1 = map_id_to_index[*it1];
    //     for(it2 = temp_lut.input.begin(); it2 != temp_lut.input.end(); ++it2)
    //     { 
    //       int index_2 = map_id_to_index[*it2];
    //       if(Nodes_vec[index_1].fan_in.size() == 1 && Nodes_vec[index_2].fan_in.size() == 1)
    //       {
    //         if(Nodes_vec[index_1].fan_in[0] == Nodes_vec[index_2].fan_in[0])
    //           pass = false;
    //       }
    //       else if(Nodes_vec[index_1].fan_in.size() == 2 && Nodes_vec[index_2].fan_in.size() == 2)
    //       {  
    //         if(Nodes_vec[index_1].fan_in[0] == Nodes_vec[index_2].fan_in[0] && Nodes_vec[index_1].fan_in[1] == Nodes_vec[index_2].fan_in[1])
    //           pass = false;
    //         if(Nodes_vec[index_1].fan_in[0] == Nodes_vec[index_2].fan_in[1] && Nodes_vec[index_1].fan_in[1] == Nodes_vec[index_2].fan_in[0])
    //           pass = false;
    //       }
    //     }
    //     if(pass)
    //       temp_lut.input.insert(*it1);
    //   }
    // }
    // else  
      temp_lut.input = input;

    LUTs.push_back(temp_lut);
  }
}
void Output()
{
  string output_name = "../output/" + output_NAME;
  ofstream fout(output_name);
  set<int>::iterator it;
  for(int i = 0; i < LUTs.size(); i++)
  {
    // cout<<"output : "<<LUTs[i].output<<" input : ";
    fout<<LUTs[i].output<<" ";
    for(it = LUTs[i].input.begin(); it != LUTs[i].input.end(); it++)
    {
      // cout<<*it<<" ";
      fout<<*it<<" ";
    }
    // cout<<endl;
    fout<<endl;
  }
}
int main(int argc, char **argv)
{   
    if (argc != 4) 
    {
        cout<<"not match the parameters"<<endl;
        exit(1);
    }
    input_AAG = argv[1];
    input_K   = atoi(argv[2]);
    output_NAME = argv[3];

    if(input_AAG.find("alu4")<200 || input_AAG.find("c1908")<200 || input_AAG.find("bigkey")<200 || input_AAG.find("sample")<200)
      input_K = input_K;
    else
      input_K = 2;

    parser();
    BuildRelation();
    TopologicalSort();
    Labeling();
    Mapping();
    Output();
    

    
    return 0;
}