#include <iostream>
#include <bits/stdc++.h>
// [[Rcpp::depends(RcppArmadillo)]] // 这行代码是Rcpp的依赖声明，用于支持RcppArmadillo库
#include <list>
#include <stack>

using namespace std;
using namespace arma;
using namespace Rcpp;

// 定义一个类来表示图
class Graph {
  // 顶点的数量
  int V;
  
  // 用于存储拓扑排序结果的数组
  int orderdfs[1000];
  
  // 指向包含邻接表的数组的指针
  list<int>* adj;
  
  // topologicalSort函数使用的辅助函数
  void topologicalSortUtil(int v, bool visited[],
                           stack<int>& Stack);
  
public:
  
  // 公共数组，用于存储拓扑排序结果
  int ordd[10];
  // 构造函数
  Graph(int V);
  
  // 向图中添加边的函数
  void addEdge(int v, int w);
  
  // 打印整个图的拓扑排序
  int* topologicalSort();
  
  // 数组转换为Rcpp::IntegerVector的函数（未实现）
  //Rcpp::IntegerVector ToArray();
};

// 构造函数，初始化顶点数量和邻接表
Graph::Graph(int V)
{
  this->V = V;
  adj = new list<int>[V];
}

// 添加边的函数
void Graph::addEdge(int v, int w)
{
  // 将w添加到v的邻接表中
  adj[v].push_back(w);
}

// topologicalSort函数使用的递归辅助函数
void Graph::topologicalSortUtil(int v, bool visited[],
                                stack<int>& Stack)
{
  // 标记当前节点为已访问
  visited[v] = true;
  
  // 递归所有与该顶点相邻的顶点
  list<int>::iterator i;
  for (i = adj[v].begin(); i != adj[v].end(); ++i)
    if (!visited[*i])
      topologicalSortUtil(*i, visited, Stack);
    
    // 将当前顶点推入栈中，栈中存储结果
    Stack.push(v);
}

// 执行拓扑排序的函数
int* Graph::topologicalSort()
{
  stack<int> Stack;
  
  // 将所有顶点标记为未访问
  bool* visited;
  visited = new bool[V];
  for (int i = 0; i < V; i++)
    visited[i] = false;
  
  // 从所有顶点开始递归调用辅助函数
  int rt=0;
  for (int i = 0; i < V; i++)
    if (visited[i] == false)
      topologicalSortUtil(i, visited, Stack);
    
    // 打印栈的内容
    while (Stack.empty() == false) {
      ordd[rt]=Stack.top();
      cout << Stack.top() << " ";
      rt++;
      Stack.pop();
    }
    
    return(ordd);
    
}

// 驱动代码
// [[Rcpp::export]]
Rcpp::IntegerVector DFS(SEXP amattt,SEXP n_nodes)
{
  arma::mat amatt = as<arma::mat>(amattt);
  const int n=as<int>(n_nodes);
  
  Graph g(n);
  // 根据邻接矩阵构建图
  for(int t = 0; t < n; t++){
    for(int r=0; r < n; r++){ 
      if(amatt(t,r)==1) g.addEdge(t, r);
    }
  }
  
  cout << "以下是给定图的拓扑排序 \n";
  
  // 执行拓扑排序并返回结果
  int* ordd2;
  ordd2=g.topologicalSort();
  
  // 将结果转换为Rcpp::IntegerVector并返回
  return Rcpp::IntegerVector::create(ordd2[0],ordd2[1],ordd2[2],ordd2[3],
                                     ordd2[4],ordd2[5],ordd2[6],ordd2[7],
                                                                     ordd2[8],ordd2[9]);
}
