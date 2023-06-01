/*
*
* NGraph++ : a simple graph library
*
* v. 2.5  03/15/2010
*
* Mathematical and Computational Sciences Division
* National Institute of Technology,
* Gaithersburg, MD USA
*
*
* This software was developed at the National Institute of Standards and
* Technology (NIST) by employees of the Federal Government in the course
* of their official duties. Pursuant to title 17 Section 105 of the
* United States Code, this software is not subject to copyright protection
* and is in the public domain. NIST assumes no responsibility whatsoever for
* its use by other parties, and makes no guarantees, expressed or implied,
* about its quality, reliability, or any other characteristic.
*
*/





#ifndef NGRAPH_H_
#define NGRAPH_H_


// DIRECTED GRAPH (with in-out adjacency list)
//

#include <iostream>
#include <set>
#include <map>
#include <utility>      // for std::pair
#include <iterator>     // for inserter()
#include <vector>       // for exporting edge_list
#include "set_ops.hpp"



//
// T is the vertex type
//
//  An adjacency graph format lists for each vertex, a set of neighbors
//  the it connects to (outlinks) and optionally another set of neighbors
//  which connect to it.  (The second set provides a quick way to find
//  out who is pointing to you.)
//
//     vertex  {out-neighbors}  {in-neighbors}
//

/**

    @brief A mathematical graph object: a simple, directed, connected graph, 
    where nodes are of arbitrary type (colores, cities, names, etc.)  
    Operations for adding and removing edges and vertices, together with
    functions for finding neighbors, subgraphs, and other properties are
    included.

   Example:


<pre>
    enum color {blue, green, red, yellow, pink, black};
    tGraph<color> A;

    A.insert_edge(blue, red);
    A.insert_edge(yellow, blue);
    A.insert_edge(blue, black);

    tGraph<color> B(A), S=A.subgraph(blue);;

</pre>


*/
namespace NGraph
{

template <typename T>
class tGraph
{

  public:

    typedef T vertex;
    typedef T value_type;
    typedef std::pair<vertex,vertex> edge;
    typedef std::set<vertex> vertex_set;
    typedef std::set<edge> edge_set;
    typedef std::pair<vertex_set, vertex_set> in_out_edge_sets;
    typedef std::map<vertex, in_out_edge_sets>  adj_graph;

    typedef typename edge_set::iterator edge_iterator;
    typedef typename edge_set::const_iterator const_edge_iterator;

    typedef typename vertex_set::iterator vertex_iterator;
    typedef typename vertex_set::const_iterator const_vertex_iterator;

  private:
  
    adj_graph G_;
    unsigned int num_edges_;


  public:

    // 
    //

    /**
      tGraph::iterator refers to pair<const vertex, in_out_edge_sets> .

<pre>
      tGraph::const_iterator p = G.begin();
      
      tGraph::vertex a = node(p);
      tGraph::vertex_set  &in_edges = in_neighbors(p);
      tGraph::vertex_set  &out_edges = out_neighbors(p);

</pre>
    */

    typedef typename adj_graph::iterator iterator;
    typedef typename adj_graph::const_iterator const_iterator;

    typedef iterator node_iterator;
    typedef const_iterator const_node_iterator;

    unsigned int num_vertices() const { return G_.size(); }
    unsigned int num_nodes() const { return G_.size(); }
    unsigned int num_edges() const { return num_edges_; }

    iterator begin() { return G_.begin(); }
    const_iterator begin() const { return G_.begin(); }
    iterator end() { return G_.end(); }
    const_iterator end() const { return G_.end(); }
    
    
    tGraph(): G_(), num_edges_(0){};
    tGraph(const tGraph &B) : G_(B.G_), num_edges_(B.num_edges_){}


    const vertex_set &in_neighbors(const vertex &a) const
              { return (G_.find(a))->second.first;}
          vertex_set &in_neighbors(const vertex &a)      
              { return (G_[a]).first; }

    const vertex_set &out_neighbors(const vertex &a)const
                {return G_.find(a)->second.second;}
          vertex_set &out_neighbors(const vertex &a)      
                {return G_[a].second; }


     iterator find(const vertex &a)
     {
        return G_.find(a);
     }

     const_iterator find(const vertex  &a) const
     {
        return G_.find(a);
     }


     unsigned int in_degree(const vertex &a) const
     {
        return in_neighbors(a).size();
     }

     unsigned int out_degree(const vertex &a) const
     {
        return out_neighbors(a).size();
     }

     unsigned int degree(const vertex &a) const
     {
        return out_neighbors(a).size();
     }


    void insert_vertex(const vertex &a)
    {
      G_[a];
    }

    void insert_new_vertex_inout_list(const vertex &a, const vertex_set &IN,
              const vertex_set &OUT)
    {
      typename adj_graph::iterator p = G_.find(a);

      // exit if "a" already in graph
      if (p != G_.end())
      {
            // remove the old number of OUT vertices
            num_edges_ -= p->second.second.size();
      }

      G_[a] = make_pair(IN, OUT);
      num_edges_ += OUT.size();
    }



    void insert_edge(iterator pa, iterator pb)
    {
      if (pa==pb) return;
      vertex a = node(pa);
      vertex b = node(pb);

      unsigned int old_size = out_neighbors(pa).size();

      out_neighbors(pa).insert(b);
      in_neighbors(pb).insert(a);

      unsigned int new_size = out_neighbors(a).size();
      if (new_size > old_size)
      {
          num_edges_++;
      }

    }


    void insert_edge(const vertex &a, const vertex& b)
    {
      if (a==b) return;

      iterator pa = find(a);
      if (pa == G_.end())
      {
          insert_vertex(a);
          pa = find(a);
      }

      iterator pb = find(b);
      if (pb == G_.end())
      {
          insert_vertex(b);
          pb = find(b);
      }

      
      insert_edge( pa, pb );
    }

   void insert_undirected_edge(const vertex &a, const vertex &b)
   {
      (a < b ) ?  insert_edge(a,b) : insert_edge(b,a);
   }


    void insert_edge(const edge &E)
    {
        insert_edge(E.first, E.second);
    }

    void insert_undirected_edge(const edge &E)
    {
      insert_undirected_edge(E.first, E.second);
    }

    tGraph(const edge_set &E)
    {
      for (typename edge_set::const_iterator p = E.begin(); p != E.end(); p++)
        insert_edge(*p);
    }

   

    void remove_edge(iterator pa, iterator pb)
    {
      unsigned int old_size = out_neighbors(pa).size();
      out_neighbors(pa).erase(node(pb));
      in_neighbors(pb).erase(node(pa));
      if (out_neighbors(pa).size() < old_size)
        num_edges_ --;

    }
   

    void remove_edge(const vertex &a, const vertex& b)
    {
      remove_edge( find(a), find(b) );
    }

    void remove_edge(const edge &E)
    {
        remove_edge(E.first, E.second);
    }

    void remove_undirected_edge(const vertex &a, const vertex& b)
    {
      (a < b) ? remove_edge(a,b) : remove_edge(b,a);
    }

    void remove_undirected_edge(const edge &e)
    {
      remove_undirected_edge(e.first, e.second);
    }

    void remove_vertex(iterator pa)
    {
      vertex_set & out_edges = out_neighbors(pa);
      vertex_set & in_edges =  in_neighbors(pa);

      // remove out-going edges
      for (typename vertex_set::iterator p = out_edges.begin(); 
                  p!=out_edges.end(); p++)
      {
          remove_edge(pa, find(*p));
      }


      // remove in-coming edges
      for (typename vertex_set::iterator p = in_edges.begin(); 
                  p!=in_edges.end(); p++)
      {
          remove_edge(find(*p), pa);
      }


      G_.erase(node(pa));
    }


    void remove_vertex(const vertex &a)
    {

          remove_vertex( find(a));

    }

  /**
        Is vertex 'a' included in graph?

        @return true, if vertex a is present; false, otherwise.
  */
   bool includes_vertex(const vertex &a) const
   {
        return  (G_.find(a) != G_.end());
   }


   /**
        Is edge (a,b) included in graph?

        @return true, if edge is present; false, otherwise.

   */
   bool includes_edge(const vertex &a, const vertex &b) const
   {
      return (includes_vertex(a)) ? 
          includes_elm(out_neighbors(a),b): false;
      
      //const vertex_set &out = out_neighbors(a);
      // return ( out.find(b) != out.end());

   }

  
  bool includes_edge(const edge& e) const
  {
    return includes_edge(e.first, e.second);
  }

    // convert to a simple edge list for exporting
    //
    /**
        Create a new representation of graph as a list
        of vertex pairs (a,b).

        @return std::vector<edge> a vector (list) of vertex pairs
    */
    std::vector<edge> edge_list() const;


/**
    @param A vertex set of nodes (subset of G)
    @return a new subgraph containing all nodes of A
*/
    tGraph subgraph(const vertex_set &A) const  
    {
        tGraph G;

        for (typename vertex_set::const_iterator p = A.begin(); p!=A.end(); p++)
        {
           
            vertex_set new_in =  (A * in_neighbors(*p));
            vertex_set new_out = (A * out_neighbors(*p));

            G.insert_new_vertex_inout_list(*p, new_in, new_out);
        }
        return G;
    }

    

  void print() const;


/* tGraph iterator methods */

/**
    @param p tGraph::const_iterator
*/
static const vertex &node(const_iterator p) 
{
    return p->first;
}

static const vertex &node(iterator p) 
{
    return p->first;
}


static const vertex_set & in_neighbors(const_iterator p)
    { return (p->second).first; }

static    vertex_set & in_neighbors(iterator p)      
      { return (p->second).first; }

static const vertex_set& out_neighbors(const_iterator p) 
      { return (p->second).second; }

static     vertex_set& out_neighbors(iterator p) 
      { return (p->second).second; }

static  unsigned int num_edges(const_iterator p)  
  {
     return out_neighbors(p).size();
  }
    
static  unsigned int degree(const_iterator p)  
  {
     return out_neighbors(p).size();
  }
    
static  unsigned int num_edges(iterator p)  
  {
     return out_neighbors(p).size();
  }
    
/**
    @param p tGraph::const_iterator
    @return number of edges going out (out-degree) at node pointed to by p.
*/
static  unsigned int out_degree(const_iterator p) 
  {
    return (p->second).second.size();
  }

static  unsigned int out_degree(iterator p) 
  {
    return (p->second).second.size();
  }



/**
      abosrb(a,b):   'a' absorbs 'b', b gets removed from graph
*/
void absorb(iterator pa, iterator pb)
{
    if (pa == pb)
      return;

    
    // first remove edge (a,b) to avoid self-loops
    remove_edge(pa, pb);

    // chnage edges (b,i) to a(i,j)
    //
    {
    vertex_set b_out = out_neighbors(pb);
    for (typename vertex_set::iterator p = b_out.begin(); 
              p!=b_out.end(); p++)
    {
      iterator pi = find(*p);
      remove_edge(pb, pi);
      //std::cout<<"\t remove_edge("<<node(pb)<< ", " << node(pi) <<")\n";
      insert_edge(pa, pi);
      //std::cout<<"\t insert_edge("<<node(pa)<< ", " << node(pi) <<")\n";
    }
    }

    // change edges (i,b) to (i,a)
    {
    vertex_set b_in = in_neighbors(pb);
    for (typename vertex_set::iterator p = b_in.begin(); 
              p!=b_in.end(); p++)
    {
      iterator pi = find(*p);
      remove_edge(pi, pb);
      //std::cout<<"\t remove_edge("<<node(pi)<< ", " << node(pb) <<")\n";
      insert_edge(pi, pa);
      //std::cout<<"\t insert_edge("<<node(pi)<< ", " << node(pa) <<")\n";
    }
    }


    //std::cout<<"\t about to remove vertex: "<< node(pb) << "\n";
    remove_vertex(pb);
    //std::cout<<"\t removed_vertex.\n";


}


void absorb(vertex a, vertex b)
{
    if (a == b)
      return ;


    absorb( find(a), find(b) );


}


};
// end tGraph<T>


typedef tGraph<unsigned int> Graph;
typedef tGraph<int> iGraph;
typedef tGraph<std::string> sGraph;


template <class T>
std::vector<typename tGraph<T>::edge> tGraph<T>::edge_list() const
    {
        //std::vector<tGraph::edge> E(num_edges());
        std::vector<typename tGraph<T>::edge> E;

        for (typename tGraph::const_iterator p = begin(); p!=end(); p++)
        {
            const vertex &a = tGraph::node(p);
            const vertex_set &out = tGraph::out_neighbors(p);
            for (typename vertex_set::const_iterator t = out.begin(); 
                        t != out.end(); t++)
            {
                E.push_back( edge(a, *t));
            }
        }
        return E;
    }

// global functions
//

template <typename T>
std::istream & operator>>(std::istream &s, tGraph<T> &G)
{

    while (!s.eof())
    {
      T v1, v2;
      s >> v1 >> v2;
      G.insert_edge(v1, v2);
    }

    return s;

}

template <typename T>
std::ostream & operator<<(std::ostream &s, const tGraph<T> &G)
{
  for (typename tGraph<T>::const_node_iterator p=G.begin(); p != G.end(); p++)
  {
    const typename tGraph<T>::vertex_set &out = tGraph<T>::out_neighbors(p);
    for ( typename tGraph<T>::vertex_set::const_iterator q=out.begin(); 
                q!=out.end(); q++)
       s << p->first << " " << *q << "\n";
  }
  return s;
}


template <typename T>
void tGraph<T>::print() const 
    {

       std::cerr << "# vertices: " <<  num_vertices()  << "\n";
       std::cerr << "# edges:    " <<  num_edges()  << "\n";

        for (const_iterator p=G_.begin(); 
              p != G_.end(); p++)
        {
          const vertex_set   &out =  out_neighbors(p);

          for (typename vertex_set::const_iterator q=out.begin(); 
                          q!=out.end(); q++)
              std::cerr << p->first << "  -->  " << *q << "\n";
        }
        std::cerr << std::endl;

    }

}
// namespace NGraph



#endif
// NGRAPH_H_
