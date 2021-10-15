#ifndef __UNDIRECTED_GRAPH_T_HPP__
#define __UNDIRECTED_GRAPH_T_HPP__

/*
 *
 *  This file contains the complete source code of the undirected_graph_t class of the AEGIS
 *    package to analyze the structure of undirected, unweighted and simple graph. While having
 *    one single file is not the most clear and organized choice for source codes, this format has
 *    been chosen to faciliate portability of the code.
 *
 *  Compilation example: g++ -O3 my_code.cpp
 *
 *  Author:  Antoine Allard
 *  WWW:     antoineallard.info
 *  Date:    July 2019
 *
 *  Version: 1.0.8
 *
 *
 *  Copyright (C) 2019 Antoine Allard
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as published
 *  by the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 */

// Standard Template Library
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>


class undirected_graph_t
{
  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  // Internal storage of the graph object.
  public:
    // Name to ID conversion. Every vertex is assigned a numerical ID in [0, |V|).
    std::map<std::string, int> Name2ID;
    // Edgelist.
    std::set< std::pair<int, int> > edgelist;
    // Adjacency list.
    std::vector<std::vector<int> > adjacency_list;
  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  // Typedefs and accessors.
  public:
    typedef std::set< std::pair<int, int> >::iterator edgelist_iterator;
    edgelist_iterator edgelist_begin() { return edgelist.begin(); }
    edgelist_iterator edgelist_end()   { return edgelist.end();   }
  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  // Objects/functions related to outputs.
  public:
    // Custom types to indicate how to identify vertices in outputs.
    enum vID_t { vID_none, vID_num, vID_name };
    static const bool header_true = true;
    static const bool header_false = false;
    // ID to name conversion.
    std::vector< std::string > ID2Name;
    // Build the ID2Name vector.
    void build_ID2Name();
  private:
    // Default width of columns.
    static const int default_column_width = 15;
    // Available vertex properties.
    std::set< std::string > available_vertex_prop;
    // Available vertex integer properties.
    std::set< std::string > available_vertex_integer_prop;
    // Headers associated with vertex properties.
    std::map< std::string, std::string > v_prop_header;
  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  // Properties.
  public:
    // Properties of the graph.
    std::map< std::string, double> g_prop;
    // Properties of vertices.
    std::map< std::string, std::vector<double> > v_prop;
    // Statistics of vertices.
    std::map< std::string, std::map<int, double> > v_stat;
    // List of all triangles.
    std::vector< std::vector<int> > triangles;
    // Onion spectrum.
    std::map< int, std::vector<double> > onion_spect;
    // Shortest path lengths.
    int shortest_path_length(int v1, int v2);
    std::vector<double> average_shortest_path_lengths;
    std::vector< std::vector<int> > shortest_path_lengths;
    std::vector< std::vector<double> > shortest_path_length_distribution;
    // Diameters of the components.
    std::vector<double> diameters;
    // Sizes of the connected components.
    std::vector<int> connected_components_size;
    // Connected components ordered by size.
    std::set< std::pair<int, int> > ordered_connected_components;
    // Unique integer properties of vertices.
    std::map< std::string, std::set<int> > v_class_1p;
    std::map< std::pair<std::string, std::string>, std::set< std::pair<int, int> > > v_class_2p;
    // Number of vertices in each unique integer properties of vertices.
    std::map< std::string, std::map<int, int> > v_class_count_1p;
    std::map< std::pair<std::string, std::string>, std::map<std::pair<int, int>, int> > v_class_count_2p;
    // Average value of a property binned in function of another integer property.
    std::map< std::pair<std::string, std::string>, std::map<int, double> > v_class_average_1p;
    // Number of edges in each unique integer properties of edges.
    std::map< std::string, std::map< std::pair<int, int>, int> > e_class_count_1p;
    std::map< std::pair< std::string, std::string >, std::map< std::multiset< std::pair<int, int> >, int> > e_class_count_2p;
  private:
    // Function verifying whether a vertex property exists.
    void is_vertex_property(std::string prop);
    // Function verifying wether a vertex property qualifies as a "integer" property.
    void is_vertex_integer_property(std::string prop);
  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  // Functions loading/saving edgelists and vertices/edges properties.
  private:
    void add_edge(std::string name1_str, std::string name2_str, int& nb_vertices);
  public:
    // Loads the graph structure from an edgelist in a file.
    void load_graph_from_edgelist_file(std::string edgelist_filename);
    // Loads the graph structure from a edgelist object.
    void load_graph_from_edgelist_object(std::vector< std::pair<std::string, std::string> > edgelist_object);
    void load_graph_from_edgelist_object(std::vector< std::pair<int, int> > edgelist_object);
    // Saves the edgelist in a file.
    void save_edgelist(std::string edgelist_filename, std::string id = "names", bool header = header_true, int base = 0);
    void save_edgelist_using_names(std::string edgelist_filename, bool header = header_true);
    void save_edgelist_using_num(std::string edgelist_filename, bool header = true, int _base = 0);
    // Loads vertices properties from a file.
    void load_vertices_properties(std::string prop_filename, int usecol, std::string prop, std::string prop_header = "");
    // // Loads edges properties from a file.
    // void load_edges_properties(std::string prop_filename, int usecol, std::string prop, std::string prop_header = "");
    // Outputing the vertices properties (the last 3 inputs can be omitted and/or put in any order).
    void save_vertices_properties(std::string filename, std::vector<std::string> props_id, vID_t vID = vID_name, int width = default_column_width, bool header = header_true);
    void save_vertices_properties(std::string filename, std::vector<std::string> props_id, vID_t vID,            bool header,                      int width = default_column_width)                { save_vertices_properties(filename, props_id, vID, width, header); };
    void save_vertices_properties(std::string filename, std::vector<std::string> props_id, bool header,          vID_t vID = vID_name,             int width = default_column_width)                { save_vertices_properties(filename, props_id, vID, width, header); };
    void save_vertices_properties(std::string filename, std::vector<std::string> props_id, bool header,          int width,                        vID_t vID = vID_name)                            { save_vertices_properties(filename, props_id, vID, width, header); };
    void save_vertices_properties(std::string filename, std::vector<std::string> props_id, int width,            vID_t vID = vID_name,             bool header = header_true)                       { save_vertices_properties(filename, props_id, vID, width, header); };
    void save_vertices_properties(std::string filename, std::vector<std::string> props_id, int width,            bool header,                      vID_t vID = vID_name)                            { save_vertices_properties(filename, props_id, vID, width, header); };
    // Outputing the number of vertices in each vertex class.
    void save_number_vertices_in_classes(std::string filename, std::string prop, bool header = header_true, int width = default_column_width);
    void save_number_vertices_in_classes(std::string filename, std::string prop, int width,                 bool header = header_true)         { save_number_vertices_in_classes(filename, prop, header, width); };
    void save_number_vertices_in_classes(std::string filename, std::pair<std::string, std::string> props,   bool header = header_true, int width = default_column_width);
    void save_number_vertices_in_classes(std::string filename, std::pair<std::string, std::string> props,   int width,                 bool header = header_true)         { save_number_vertices_in_classes(filename, props, header, width); };
    // Outputing the number of edges in each vertex class.
    void save_number_edges_in_classes(std::string filename, std::string prop, bool header = header_true, int width = default_column_width);
    void save_number_edges_in_classes(std::string filename, std::string prop, int width, bool header = header_true)                         { save_number_edges_in_classes(filename, prop, header, width); };
    void save_number_edges_in_classes(std::string filename, std::pair<std::string, std::string> props, bool header = header_true, int width = default_column_width);
    void save_number_edges_in_classes(std::string filename, std::pair<std::string, std::string> props, int width, bool header = header_true)                         { save_number_edges_in_classes(filename, props, header, width); };
  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  // Functions building secondary objects related to the graph.
  public:
    // Builds the adjacency list.
    void build_adjacency_list();
  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  // Functions extracting properties of graph/vertices/edges.
  public:
    // Computes the average degree.
    void average_degree();
    // Computes the average degree of neighbors. It can be normalized or not.
    void average_degree_of_neighbors(bool normalized = false);
    // Computes the average degree of neighbors spectrum. It can be normalized or not.
    void average_degree_of_neighbors_spectrum(bool normalized = false);
    // Computes the average local clustering coefficient.
    void average_local_clustering_coefficient();
    // Computes the average value of a vertex property.
    void average_vertex_prop(std::string prop);
    // Computes the average value of of a vertex property in function of the class of vertex.
    void average_vertex_prop_by_classes_of_vertices(std::string prop1, std::string prop2);
    // Computes the average value of a vertex property of neighbors.
    void average_vertex_prop_of_neighbors(std::string prop);
    // Computes the closeness centrality.
    void closeness_centrality(bool correct_for_component_size = true);
    // Computes the degree correlation coefficient.
    void degree_correlation_coefficient();
    // Extracts the degree of vertices.
    void degrees();
    // Degree distribution.
    void degree_distribution();
    // Clustering spectrum.
    void clustering_spectrum();
    // Onion_spectrum.
    void onion_spectrum();
    // Erases vertex (integer) properties.
    void erase_vertex_integer_property(std::string prop);
    void erase_vertex_property(std::string prop);
    // Computes the global clustering coefficient.
    void global_clustering_coefficient();
    // Computes the harmonic centrality.
    void harmonic_centrality();
    // Extracts the k-core decomposition.
    void kcore_decomposition();
    // Extracts the local clustering coefficients.
    void local_clustering_coefficients();
    // Computes the modularity given the memberships provided via its name prop.
    void modularity(std::string prop, std::string g_modularity_prop_name = "modularity");
    // Adding new vertex (integer) properties.
    void new_vertex_integer_property(std::string prop, std::string prop_header = "");
    void new_vertex_property(std::string prop, std::string prop_header = "");
    // // Adding a new edge property.
    // void new_edge_property(std::string prop, std::string prop_header = "");
    // Counts the number of triangles around every vertex.
    void number_of_triangles();
    // Extracts the onion decomposition (includes the k-core decomposition as a by-product).
    void onion_decomposition();
    // Identifies the connected components to which vertices belong.
    void survey_connected_components();
    // Finds the length of all shortest paths (equal to the number of vertices if no path exists).
    void survey_shortest_path_lengths();
    // Computes the length distribution of the shortest paths for each connected component.
    void survey_shortest_path_length_distribution();
    // Compiles a list of all triangles in the graph.
    void survey_triangles();
    // Computes the average shortest path lengths and the diameters.
    void topological_dimensions();
    // Identifies unique integer properties of vertices (ex. degree or degree-layer).
    void identify_classes_of_vertices(std::string prop);
    void identify_classes_of_vertices(std::pair<std::string, std::string> props);
    // Counts the number of vertices in each class of vertices.
    void count_vertices_in_classes(std::string prop);
    void count_vertices_in_classes(std::pair<std::string, std::string> props);
    // Counts the number of edges in each class of edges.
    void count_edges_in_classes(std::string prop);
    void count_edges_in_classes(std::pair<std::string, std::string> props);
  private:
    // Counts the number of triangles around a given vertex.
    int count_triangles_around_vertex(int v1);
    // Function associated to the extraction of the components.
    int get_root(int i, std::vector<int> &clust_id);
    void merge_clusters(std::vector<int> &size, std::vector<int> &clust_id, int nb_vertices);
  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  // Constructors (and related functions).
  private:
    // Function setting default values (to avoid requiring the C++11 standard).
    void initialization();
  public:
    // Empty constructor.
    undirected_graph_t() { initialization(); };
    // Constructor with edgelist.
    undirected_graph_t(std::string edgelist_filename) { initialization(); load_graph_from_edgelist_file(edgelist_filename); };
    undirected_graph_t(std::vector< std::pair<std::string, std::string> > edgelist_object) { initialization(); load_graph_from_edgelist_object(edgelist_object); };
    undirected_graph_t(std::vector< std::pair<int,         int        > > edgelist_object) { initialization(); load_graph_from_edgelist_object(edgelist_object); };
  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  // Convert to python.
  public:
    // Graph properties.
    int get_nb_vertices();
    int get_nb_edges();
    int get_total_nb_triangles();
    int get_nb_connected_components();
    int get_shortest_path(std::string v1, std::string v2);
    double get_average_degree();
    double get_average_local_clustering_coefficient();
    double get_degree_correlation_coefficient();
    double get_global_clustering_coefficient();
    std::vector<int> get_diameter();
    std::vector<int> get_connected_components_size();
    std::vector<double> get_average_shortest_path();
    std::vector< std::string > get_vertices();
    std::vector< std::vector< double > > get_shortest_path_length_distribution();
    std::vector< std::vector< std::string > > get_triangles();
    std::vector< std::vector< std::string > > get_connected_components();
    std::vector< std::pair< std::string, std::string > > get_edges();
    // Vertex statistics.
    // std::map<int, int> get_count_vertices_in_classes(std::string prop);
    std::map<int, double> get_average_degree_of_neighbors_spectrum(bool normalized = false);
    std::map<int, double> get_clustering_spectrum();
    std::map<int, double> get_degree_distribution();
    std::map<int, std::vector<double> > get_onion_spectrum();
    // Vertex properties.
    std::map<std::string, double> get_average_degree_of_neighbors(bool normalized = false);
    std::map<std::string, double> get_closeness_centrality(bool correct_for_component_size = true);
    std::map<std::string, double> get_degrees();
    std::map<std::string, double> get_k_shells();
    std::map<std::string, double> get_harmonic_centrality();
    std::map<std::string, double> get_local_clustering_coefficients();
    std::map<std::string, double> get_nb_triangles();
    std::map<std::string, double> get_od_layers();
};










// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// Structure containing the desired method to compared strings (put shorter ones before longer ones).
struct compare_names
{
bool operator()(const std::pair<std::string, int>& lhs, const std::pair<std::string, int>& rhs) const
{
  if(lhs.first.size() == rhs.first.size())
  {
    if(lhs.first == rhs.first)
    {
      return lhs.second < rhs.second;
    }
    else
    {
      return lhs.first < rhs.first;
    }
  }
  else
  {
    return lhs.first.size() < rhs.first.size();
  }
}
};


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::initialization()
{
  // Available vertex properties.
  available_vertex_prop.insert("degrees");
  available_vertex_prop.insert("local_clust");
  available_vertex_prop.insert("kcore");
  available_vertex_prop.insert("od_layer");
  available_vertex_prop.insert("nb_triangles");
  available_vertex_prop.insert("conn_comp");
  available_vertex_prop.insert("close_cent");
  available_vertex_prop.insert("corr_close_cent");
  available_vertex_prop.insert("harmo_cent");
  available_vertex_prop.insert("avg_neigh_degrees");
  available_vertex_prop.insert("avg_neigh_degrees_norm");
  available_vertex_prop.insert("avg_neigh_local_clust");
  available_vertex_prop.insert("avg_neigh_kcore");
  available_vertex_prop.insert("avg_neigh_od_layer");
  available_vertex_prop.insert("avg_neigh_nb_triangles");
  // Available vertex integer properties.
  available_vertex_integer_prop.insert("degrees");
  available_vertex_integer_prop.insert("kcore");
  available_vertex_integer_prop.insert("od_layer");
  available_vertex_integer_prop.insert("nb_triangles");
  available_vertex_integer_prop.insert("conn_comp");
  // Headers for vertex properties.
  v_prop_header["degrees"]                = "Degree";
  v_prop_header["local_clust"]            = "Clust.";
  v_prop_header["kcore"]                  = "Coreness";
  v_prop_header["od_layer"]               = "ODlayer";
  v_prop_header["nb_triangles"]           = "NbTriang.";
  v_prop_header["conn_comp"]              = "Conn.Comp.";
  v_prop_header["close_cent"]             = "Close.Cent.";
  v_prop_header["corr_close_cent"]        = "CClose.Cent.";
  v_prop_header["harmo_cent"]             = "Harmo.Cent.";
  v_prop_header["avg_neigh_degrees"]      = "ANDegree";
  v_prop_header["avg_neigh_degrees_norm"] = "ANDegreeNorm";
  v_prop_header["avg_neigh_local_clust"]  = "ANClust.";
  v_prop_header["avg_neigh_kcore"]        = "ANCoreness";
  v_prop_header["avg_neigh_od_layer"]     = "ANODlayer";
  v_prop_header["avg_neigh_nb_triangles"] = "ANNbTriang.";
  // // Sets some default values.
  // header_true = true;
  // header_false = false;
  // default_column_width = 15;
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::count_edges_in_classes(std::string prop)
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  std::vector<double>& Vertex2Prop = v_prop[prop];

  e_class_count_1p[prop].clear();
  std::map< std::pair<int, int>, int>& PropCount = e_class_count_1p[prop];
  // ===============================================================================================

  // Loops over every edges and counts the number of edges in every unique class of edges.
  int p1, p2;
  std::pair<int, int> p;
  std::set< std::pair<int, int> >::iterator it = edgelist.begin();
  std::set< std::pair<int, int> >::iterator end = edgelist.end();
  for(; it!=end; ++it)
  {
    p1 = Vertex2Prop[it->first];
    p2 = Vertex2Prop[it->second];
    if(p1 < p2)
    {
      p = std::make_pair(p1, p2);
    }
    else
    {
      p = std::make_pair(p2, p1);
    }
    if(PropCount.find(p) == PropCount.end())
    {
      PropCount[p] = 1;
    }
    else
    {
      PropCount[p] += 1;
    }
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::count_edges_in_classes(std::pair<std::string, std::string> props)
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  std::vector<double>& Vertex2Prop1 = v_prop[props.first];
  std::vector<double>& Vertex2Prop2 = v_prop[props.second];

  e_class_count_2p[props].clear();
  std::map< std::multiset< std::pair<int, int> >, int>& PropCount = e_class_count_2p[props];
  // ===============================================================================================

  // Loops over every edges and counts the number of edges in every unique class of edges.
  int p1, p2;
  std::multiset< std::pair<int, int> > p;
  std::multiset< std::pair<int, int> >::iterator it = edgelist.begin();
  std::multiset< std::pair<int, int> >::iterator end = edgelist.end();
  for(; it!=end; ++it)
  {
    p.clear();
    p1 = Vertex2Prop1[it->first];
    p2 = Vertex2Prop2[it->first];
    p.insert(std::make_pair(p1, p2));
    p1 = Vertex2Prop1[it->second];
    p2 = Vertex2Prop2[it->second];
    p.insert(std::make_pair(p1, p2));
    if(PropCount.find(p) == PropCount.end())
    {
      PropCount[p] = 1;
    }
    else
    {
      PropCount[p] += 1;
    }
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::count_vertices_in_classes(std::string prop)
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];

  std::vector<double>& Vertex2Prop = v_prop[prop];

  std::set<int>& PropClass = v_class_1p[prop];

  v_class_count_1p[prop].clear();
  std::map<int, int>& PropCount = v_class_count_1p[prop];
  // ===============================================================================================

  // Ensures that the relevant classes of vertices have been identified.
  if(PropClass.size() == 0)
  {
    identify_classes_of_vertices(prop);
  }

  // Initializes the count structure.
  std::set<int>::iterator it = PropClass.begin();
  std::set<int>::iterator end = PropClass.end();
  for(; it!=end; ++it)
  {
    PropCount[*it] = 0;
  }

  // Loops over every vertices and counts the number of vertices in every unique class of vertex.
  int p;
  for(int v(0); v<nb_vertices; ++v)
  {
    p = Vertex2Prop[v];
    PropCount[p] += 1;
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::count_vertices_in_classes(std::pair<std::string, std::string> props)
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];

  std::vector<double>& Vertex2Prop1 = v_prop[props.first];
  std::vector<double>& Vertex2Prop2 = v_prop[props.second];

  v_class_count_2p[props].clear();
  std::map<std::pair<int, int>, int>& PropCount = v_class_count_2p[props];

  std::set< std::pair<int, int> >& PropClass = v_class_2p[props];
  // ===============================================================================================

  // Ensures that the relevant classes of vertices have been identified.
  if(PropClass.size() == 0)
  {
    identify_classes_of_vertices(props);
  }

  // Initializes the count structure.
  std::set< std::pair<int, int> >::iterator it = PropClass.begin();
  std::set< std::pair<int, int> >::iterator end = PropClass.end();
  for(; it!=end; ++it)
  {
    PropCount[*it] = 0;
  }

  // Loops over every vertices and counts the number of vertices in every unique class of vertex.
  int p1, p2;
  for(int v(0); v<nb_vertices; ++v)
  {
    p1 = Vertex2Prop1[v];
    p2 = Vertex2Prop2[v];
    PropCount[std::make_pair(p1, p2)] += 1;
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::erase_vertex_integer_property(std::string prop)
{
  is_vertex_integer_property(prop);
  available_vertex_integer_prop.erase( available_vertex_integer_prop.find(prop) );
  erase_vertex_property(prop);
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::erase_vertex_property(std::string prop)
{
  is_vertex_property(prop);
  available_vertex_prop.erase( available_vertex_prop.find(prop) );
  v_prop_header.erase( v_prop_header.find(prop) );
  v_prop.erase( v_prop.find(prop) );
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::identify_classes_of_vertices(std::string prop)
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];

  std::vector<double>& Vertex2Prop = v_prop[prop];

  v_class_1p[prop].clear();
  std::set<int>& PropClass = v_class_1p[prop];
  // ===============================================================================================

  is_vertex_property(prop);
  is_vertex_integer_property(prop);
  if(v_prop[prop].size() != nb_vertices)
  {
    std::cerr << "ERROR: The property " << prop << " has not been extracted/computed." << std::endl;
    std::terminate();
  }

  // Loops over every vertices and identifies all unique class of vertex.
  int p;
  for(int v(0); v<nb_vertices; ++v)
  {
    p = Vertex2Prop[v];
    PropClass.insert(p);
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::identify_classes_of_vertices(std::pair<std::string, std::string> props)
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];

  std::vector<double>& Vertex2Prop1 = v_prop[props.first];
  std::vector<double>& Vertex2Prop2 = v_prop[props.second];

  v_class_2p[props].clear();
  std::set< std::pair<int, int> >& PropClass = v_class_2p[props];
  // ===============================================================================================

  is_vertex_property(props.first);
  is_vertex_integer_property(props.first);
  if(v_prop[props.first].size() != nb_vertices)
  {
    std::cerr << "ERROR: The property " << props.first << " has not been extracted/computed." << std::endl;
    std::terminate();
  }

  is_vertex_property(props.second);
  is_vertex_integer_property(props.second);
  if(v_prop[props.second].size() != nb_vertices)
  {
    std::cerr << "ERROR: The property " << props.second << " has not been extracted/computed." << std::endl;
    std::terminate();
  }

  // Loops over every vertices and identifies all unique class of vertex.
  int p1, p2;
  for(int v(0); v<nb_vertices; ++v)
  {
    p1 = Vertex2Prop1[v];
    p2 = Vertex2Prop2[v];
    PropClass.insert(std::make_pair(p1, p2));
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::is_vertex_integer_property(std::string prop)
{
  if(available_vertex_integer_prop.find(prop) == available_vertex_integer_prop.end())
  {
    std::cerr << "ERROR: " << prop << " is not a valid vertex property." << std::endl;
    std::terminate();
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::is_vertex_property(std::string prop)
{
  if(available_vertex_prop.find(prop) == available_vertex_prop.end())
  {
    std::cerr << "ERROR: " << prop << " is not a valid vertex property." << std::endl;
    std::terminate();
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::load_vertices_properties(std::string prop_filename, int usecol, std::string prop, std::string prop_header)
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  new_vertex_property(prop, prop_header);
  std::vector<double>& Vertex2Prop = v_prop[prop];
  // ===============================================================================================
  // Stream object.
  std::stringstream one_line;
  // String objects.
  std::string full_line, name_str, value_str;
  // Variable.
  int idx;
  double tmp_value, dummy;
  bool is_integer = true;
  // Opens the stream and terminates if the operation did not succeed.
  std::fstream prop_file(prop_filename.c_str(), std::fstream::in);
  if( !prop_file.is_open() )
  {
    std::cerr << "Could not open file: " << prop_filename << "." << std::endl;
    std::terminate();
  }
  // Reads the propperty file line by line.
  while( !prop_file.eof() )
  {
    // Reads a line of the file.
    std::getline(prop_file, full_line);
    prop_file >> std::ws;
    one_line.str(full_line);
    one_line >> std::ws;
    one_line >> name_str >> std::ws;
    // Skips lines of comment.
    // if(name_str == "#")
    if(name_str.compare(0, 1, "#") == 0)
    {
      one_line.clear();
      continue;
    }
    // Identifies the vertex.
    if(Name2ID.find(name_str) == Name2ID.end())
    {
      one_line.clear();
      continue;
    }
    idx = Name2ID[name_str];
    // Inputs the property.
    for(int r(1); r<usecol; ++r)
    {
      one_line >> value_str >> std::ws;
    }
    one_line >> value_str >> std::ws;
    tmp_value = std::atof(value_str.c_str());
    Vertex2Prop[idx] = tmp_value;
    if( is_integer && (std::modf(tmp_value, &dummy) != 0) )
    {
      is_integer = false;
    }
    // Resets the stringstream.
    one_line.clear();
  }
  // Closes the stream.
  prop_file.close();
  // ===============================================================================================
  // Initializes relevant objects of the class.
  if(is_integer)
  {
    available_vertex_integer_prop.insert(prop);
  }
  // ===============================================================================================
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::modularity(std::string prop, std::string g_modularity_prop_name)
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_edges = g_prop["nb_edges"];
  int nb_vertices = g_prop["nb_vertices"];
  std::vector<double>& Vertex2Cluster = v_prop[prop];
  // ===============================================================================================

  // Checks if the property is an integer.
  is_vertex_integer_property(prop);

  // Containers counting edges/stubs in and within communities.
  std::map<int, double> a;
  std::map<int, double> e;
  // Covering the edgelist.
  edgelist_iterator it = edgelist_begin();
  edgelist_iterator end = edgelist_end();
  for(int v1, v2, c1, c2, d1, d2; it!=end; ++it)
  {
    // Gets information about the vertices.
    v1 = it->first;
    c1 = Vertex2Cluster[v1];
    v2 = it->second;
    c2 = Vertex2Cluster[v2];
    // Counts intra-community edges if they belong to the same community.
    if(c1 == c2)
    {
      if( e.find(c1) == e.end() )
      {
        e[c1] = 0;
      }
      e[c1] += 1;
    }
    // Counts the stubs in each community.
    if( a.find(c1) == a.end() )
    {
      a[c1] = 0;
    }
    a[c1] += 1;
    if( a.find(c2) == a.end() )
    {
      a[c2] = 0;
    }
    a[c2] += 1;
  }
  // Computes the modularity.
  double modularity = 0;
  std::map<int, double>::iterator it1 = e.begin();
  std::map<int, double>::iterator end1 = e.end();
  for(; it1!=end1; ++it1)
  {
    modularity += it1->second / nb_edges;
  }
  it1 = a.begin();
  end1 = a.end();
  for(; it1!=end1; ++it1)
  {
    modularity -= std::pow( it1->second / (2 * nb_edges), 2);
  }
  // ===============================================================================================
  // Initializes relevant objects of the class.
  // if(g_modularity_prop_name == "")
  // {
  //   g_modularity_prop_name = "modularity";
  // }
  g_prop[g_modularity_prop_name] = modularity;
  // ===============================================================================================
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::new_vertex_integer_property(std::string prop, std::string prop_header)
{
  new_vertex_property(prop, prop_header);
  if(available_vertex_integer_prop.find(prop) != available_vertex_integer_prop.end())
  {
    std::cerr << "ERROR: Vertex integer property " << prop << " already exists." << std::endl;
    std::terminate();
  }
  available_vertex_integer_prop.insert(prop);
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::new_vertex_property(std::string prop, std::string prop_header)
{
  if(prop_header == "")
  {
    prop_header = prop;
  }
  if(available_vertex_prop.find(prop) != available_vertex_prop.end())
  {
    std::cerr << "ERROR: Vertex property " << prop << " already exists." << std::endl;
    std::terminate();
  }
  available_vertex_prop.insert(prop);
  int nb_vertices = g_prop["nb_vertices"];
  v_prop_header[prop] = prop_header;
  v_prop[prop] = std::vector<double>(nb_vertices, 0);
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::save_number_edges_in_classes(std::string filename, std::string prop, bool header, int width)
{
  // Stream objects.
  std::fstream output_file;
  // Opens the stream and terminates if the operation did not succeed.
  output_file.open(filename.c_str(), std::fstream::out);
  if( !output_file.is_open() )
  {
    std::cerr << "ERROR: Could not open file: " << filename << "." << std::endl;
    std::terminate();
  }


  // ===============================================================================================
  // Initializes relevant objects of the class.
  std::map< std::pair<int, int>, int>& PropCount = e_class_count_1p[prop];
  // ===============================================================================================


  if(PropCount.size() == 0)
  {
    std::cerr << "ERROR: the number of vertices in classes <" << prop << "> is not available." << std::endl;
    std::terminate();
  }

  // Prints the header (if required).
  if(header)
  {
    output_file << "#" << std::setw(width - 2) << v_prop_header[prop] << "1 ";
    output_file << std::setw(width - 1) << v_prop_header[prop] << "2 ";
    output_file << std::setw(width) << "Nb. edges" << " ";
    output_file << std::endl;
  }

  // Prints the properties.
  std::map< std::pair<int, int>, int>::iterator it = PropCount.begin();
  std::map< std::pair<int, int>, int>::iterator end = PropCount.end();
  for(; it!=end; ++it)
  {
    output_file << std::setw(width) << it->first.first << " ";
    output_file << std::setw(width) << it->first.second << " ";
    output_file << std::setw(width) << it->second << " ";
    output_file << std::endl;
  }

  // Closes the stream.
  output_file.close();
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::save_number_edges_in_classes(std::string filename, std::pair<std::string, std::string> props, bool header, int width)
{
  // Stream objects.
  std::fstream output_file;
  // Opens the stream and terminates if the operation did not succeed.
  output_file.open(filename.c_str(), std::fstream::out);
  if( !output_file.is_open() )
  {
    std::cerr << "ERROR: Could not open file: " << filename << "." << std::endl;
    std::terminate();
  }


  // ===============================================================================================
  // Initializes relevant objects of the class.
  std::map< std::multiset< std::pair<int, int> >, int>& PropCount = e_class_count_2p[props];
  // ===============================================================================================


  if(PropCount.size() == 0)
  {
    std::cerr << "ERROR: the number of vertices in classes <" << props.first << "," << props.second << "> is not available." << std::endl;
    std::terminate();
  }

  // Prints the header (if required).
  if(header)
  {
    output_file << "#" << std::setw(width - 2) << v_prop_header[props.first] << "1 ";
    output_file << std::setw(width - 1) << v_prop_header[props.second] << "1 ";
    output_file << std::setw(width - 1) << v_prop_header[props.first] << "2 ";
    output_file << std::setw(width - 1) << v_prop_header[props.second] << "2 ";
    output_file << std::setw(width) << "Nb. edges" << " ";
    output_file << std::endl;
  }

  // Prints the properties.
  std::map< std::multiset< std::pair< int, int > >, int>::iterator it = PropCount.begin();
  std::map< std::multiset< std::pair< int, int > >, int>::iterator end = PropCount.end();
  for(; it!=end; ++it)
  {
    output_file << std::setw(width) << (it->first).begin()->first << " ";
    output_file << std::setw(width) << (it->first).begin()->second << " ";
    output_file << std::setw(width) << (++(it->first.begin()))->first << " ";
    output_file << std::setw(width) << (++(it->first.begin()))->second << " ";
    output_file << std::setw(width) << it->second << " ";
    output_file << std::endl;
  }

  // Closes the stream.
  output_file.close();
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::save_number_vertices_in_classes(std::string filename, std::string prop, bool header, int width)
{
  // Stream objects.
  std::fstream output_file;
  // Opens the stream and terminates if the operation did not succeed.
  output_file.open(filename.c_str(), std::fstream::out);
  if( !output_file.is_open() )
  {
    std::cerr << "ERROR: Could not open file: " << filename << "." << std::endl;
    std::terminate();
  }


  // ===============================================================================================
  // Initializes relevant objects of the class.
  std::map<int, int>& PropCount = v_class_count_1p[prop];
  // ===============================================================================================


  if(PropCount.size() == 0)
  {
    std::cerr << "ERROR: the number of vertices in classes <" << prop << "> is not available." << std::endl;
    std::terminate();
  }

  // Prints the header (if required).
  if(header)
  {
    output_file << "#" << std::setw(width - 1) << v_prop_header[prop] << " ";
    output_file << std::setw(width) << "Nb. vertices" << " ";
    output_file << std::endl;
  }

  // Prints the properties.
  std::map<int, int>::iterator it = PropCount.begin();
  std::map<int, int>::iterator end = PropCount.end();
  for(; it!=end; ++it)
  {
    output_file << std::setw(width) << it->first << " ";
    output_file << std::setw(width) << it->second << " ";
    output_file << std::endl;
  }

  // Closes the stream.
  output_file.close();
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::save_number_vertices_in_classes(std::string filename, std::pair<std::string, std::string> props, bool header, int width)
{
  // Stream objects.
  std::fstream output_file;
  // Opens the stream and terminates if the operation did not succeed.
  output_file.open(filename.c_str(), std::fstream::out);
  if( !output_file.is_open() )
  {
    std::cerr << "ERROR: Could not open file: " << filename << "." << std::endl;
    std::terminate();
  }


  // ===============================================================================================
  // Initializes relevant objects of the class.
  std::map< std::pair<int, int>, int>& PropCount = v_class_count_2p[props];
  // ===============================================================================================


  if(PropCount.size() == 0)
  {
    std::cerr << "ERROR: the number of vertices in classes <" << props.first << "," << props.second << "> is not available." << std::endl;
    std::terminate();
  }

  // Prints the header (if required).
  if(header)
  {
    output_file << "#" << std::setw(width - 1) << v_prop_header[props.first] << " ";
    output_file << std::setw(width) << v_prop_header[props.second] << " ";
    output_file << std::setw(width) << "Nb. vertices" << " ";
    output_file << std::endl;
  }

  // Prints the properties.
  std::map< std::pair<int, int>, int>::iterator it = PropCount.begin();
  std::map< std::pair<int, int>, int>::iterator end = PropCount.end();
  for(; it!=end; ++it)
  {
    output_file << std::setw(width) << it->first.first << " ";
    output_file << std::setw(width) << it->first.second << " ";
    output_file << std::setw(width) << it->second << " ";
    output_file << std::endl;
  }

  // Closes the stream.
  output_file.close();
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::save_vertices_properties(std::string filename, std::vector<std::string> props_id, vID_t vID, int width, bool header)
{
  // Stream objects.
  std::fstream output_file;
  // Opens the stream and terminates if the operation did not succeed.
  output_file.open(filename.c_str(), std::fstream::out);
  if( !output_file.is_open() )
  {
    std::cerr << "ERROR: Could not open file: " << filename << "." << std::endl;
    std::terminate();
  }


  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  // ===============================================================================================


  // Number of properties.
  int nb_props = props_id.size();

  // Checks whether all properties have been extracted/computed.
  for(int i(0); i<nb_props; ++i)
  {
    is_vertex_property(props_id[i]);
    if(v_prop[props_id[i]].size() != nb_vertices)
    {
      std::cerr << "ERROR: The property " << props_id[i] << " has not been extracted/computed." << std::endl;
      std::terminate();
    }
  }

  // Ensures that the ID2Name has been built. Otherwise, build it.
  if(ID2Name.size() != nb_vertices)
  {
    build_ID2Name();
  }


  // Prints the header (if required).
  if(header)
  {
    int i = 0;
    if(vID == vID_name || vID == vID_num)
    {
      output_file << "#" << std::setw(width - 1) << "Vertex" << " ";
      output_file << std::setw(width) << v_prop_header[props_id[i]] << " ";
      ++i;
    }
    else if(vID == vID_none)
    {
      output_file << "#" << std::setw(width - 1) << v_prop_header[props_id[i]] << " ";
      ++i;
    }
    else
    {
      std::cerr << "ERROR: Unknown vertex identifier type." << std::endl;
      std::terminate();
    }
    for(; i<nb_props; ++i)
    {
      output_file << std::setw(width) << v_prop_header[props_id[i]] << " ";
    }
    output_file << std::endl;
  }

  // Prints the properties.
  std::set< std::pair<std::string, int>, compare_names > ordered_names;
  for(int v(0); v<nb_vertices; ++v)
  {
    ordered_names.insert(std::make_pair(ID2Name[v], v));
  }
  std::set< std::pair<std::string, int> >::iterator it  = ordered_names.begin();
  std::set< std::pair<std::string, int> >::iterator end = ordered_names.end();
  for(int v; it!=end; ++it)
  {
    v = it->second;
    if(vID == vID_name)
    {
      output_file << std::setw(width) << ID2Name[v] << " ";
    }
    else if(vID == vID_num)
    {
      output_file << std::setw(width) << v << " ";
    }
    else if(vID == vID_none) { }
    else
    {
      std::cerr << "ERROR: Unknown vertex identifier type." << std::endl;
      std::terminate();
    }
    for(int i(0); i<nb_props; ++i)
    {
      output_file << std::setw(width) << v_prop[props_id[i]][v] << " ";
    }
    output_file << std::endl;
  }

  // Closes the stream.
  output_file.close();
}










// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// 0. INTERNAL MANAGMENT: GENERATING SECONDARY OBJECTS
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::build_adjacency_list()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  adjacency_list.clear();
  adjacency_list.resize(nb_vertices);
  // ===============================================================================================

  // Loops over all edges.
  int v1, v2;
  std::set< std::pair<int, int> >::iterator it = edgelist.begin();
  std::set< std::pair<int, int> >::iterator end = edgelist.end();
  for(; it!=end; ++it)
  {
    // Identifies the vertices.
    v1 = it->first;
    v2 = it->second;
    // Adds the percolated edge.
    adjacency_list[v1].push_back(v2);
    adjacency_list[v2].push_back(v1);
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::build_ID2Name()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  ID2Name.clear();
  ID2Name.resize(nb_vertices);
  // ===============================================================================================

  // Loops over all names.
  std::map<std::string, int>::iterator it = Name2ID.begin();
  std::map<std::string, int>::iterator end = Name2ID.end();
  for(; it!=end; ++it)
  {
    ID2Name[it->second] = it->first;
  }
}










// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// 1. INTERNAL MANAGMENT: AVERAGE VALUES OF VERTEX PROPERTIES
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::average_vertex_prop(std::string prop)
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  std::vector<double>& Vertex2Prop = v_prop[prop];
  // ===============================================================================================

  // Checks if the property has been extracted/computed already.
  if(v_prop[prop].size() != nb_vertices)
  {
    std::cerr << "ERROR: The property " << prop << " has not been extracted/computed." << std::endl;
    std::terminate();
  }

  // Computes the average degree of vertices.
  double sum_of_values = 0;
  for(int v(0); v<nb_vertices; ++v)
  {
    sum_of_values += Vertex2Prop[v];
  }

  // ===============================================================================================
  // Updates the properties of the graph.
  std::string avg_v_prop = "avg_" + prop;
  g_prop[avg_v_prop] = sum_of_values / nb_vertices;
  // ===============================================================================================
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::average_vertex_prop_by_classes_of_vertices(std::string prop1, std::string prop2)
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];

  std::vector<double>& Vertex2Prop1 = v_prop[prop1];
  std::vector<double>& Vertex2Prop2 = v_prop[prop2];

  std::set<int>& PropClass = v_class_1p[prop1];

  v_class_count_1p[prop1].clear();
  std::map<int, int>& PropCount = v_class_count_1p[prop1];

  v_class_average_1p[std::make_pair(prop1, prop2)].clear();
  std::map<int, double>& PropAverage = v_class_average_1p[std::make_pair(prop1, prop2)];
  // ===============================================================================================

  // Ensures that the relevant classes of vertices have been identified.
  if(PropClass.size() == 0)
  {
    identify_classes_of_vertices(prop1);
  }

  // Ensures that the relevant classes of vertices have been identified.
  if(PropCount.size() == 0)
  {
    count_vertices_in_classes(prop1);
  }

  // Initializes the count structure.
  std::set<int>::iterator it = PropClass.begin();
  std::set<int>::iterator end = PropClass.end();
  for(; it!=end; ++it)
  {
    PropAverage[*it] = 0;
  }

  // Loops over every vertices and sums the vertex property in the proper class.
  int p;
  for(int v(0); v<nb_vertices; ++v)
  {
    p = Vertex2Prop1[v];
    PropAverage[p] += Vertex2Prop2[v];
  }
  // Normalises each average.
  it = PropClass.begin();
  end = PropClass.end();
  for(; it!=end; ++it)
  {
    PropAverage[*it] /= PropCount[*it];
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::average_vertex_prop_of_neighbors(std::string prop)
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];

  std::vector<double>& Vertex2Prop = v_prop[prop];

  std::string avg_neigh_prop = "avg_neigh_" + prop;
  v_prop[avg_neigh_prop].clear();
  std::vector<double>& Vertex2ANProp = v_prop[avg_neigh_prop];
  Vertex2ANProp.resize(nb_vertices, 0);
  // ===============================================================================================

  // Ensures that the adjacency list has been built.
  if(adjacency_list.size() != nb_vertices)
  {
    build_adjacency_list();
  }

  // Computes the average vertex property of neighbors.
  int d1, v2;
  double value;
  for(int v1(0); v1<nb_vertices; ++v1)
  {
    value = 0;
    d1 = adjacency_list[v1].size();
    for(int s(0); s<d1; ++s)
    {
      v2 = adjacency_list[v1][s];
      value += Vertex2Prop[v2];
    }
    Vertex2ANProp[v1] = value / d1;
  }
}










// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// 2. INPUT
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::add_edge(std::string name1_str, std::string name2_str, int& nb_vertices)
{
  // Variables.
  int v1, v2;
  // Iterators.
  std::map< std::string, int >::iterator name_it;

  // Ignores self-loops.
  if(name1_str != name2_str)
  {

    // Is name1 new?
    name_it = Name2ID.find(name1_str);
    if(name_it == Name2ID.end())
    {
      Name2ID[name1_str] = nb_vertices;
      v1 = nb_vertices;
      ++nb_vertices;
    }
    else
    {
      v1 = name_it->second;
    }

    // Is name2 new?
    name_it = Name2ID.find(name2_str);
    if(name_it == Name2ID.end())
    {
      Name2ID[name2_str] = nb_vertices;
      v2 = nb_vertices;
      ++nb_vertices;
    }
    else
    {
      v2 = name_it->second;
    }

    // Adds the edge (multiedges are ignored).
    if(v1 > v2)
    {
      std::swap(v1, v2);
    }
    edgelist.insert(std::make_pair(v1, v2));
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::load_graph_from_edgelist_file(std::string edgelist_filename)
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  Name2ID.clear();
  edgelist.clear();
  // ===============================================================================================


  // Stream objects.
  std::fstream edgelist_file;
  std::stringstream one_line;
  // Initializes the number of vertices.
  int nb_vertices = 0;
  // String objects.
  std::string full_line, name1_str, name2_str;


  // Opens the stream and terminates if the operation did not succeed.
  edgelist_file.open(edgelist_filename.c_str(), std::fstream::in);
  if( !edgelist_file.is_open() )
  {
    std::cerr << "ERROR: Could not open file: " << edgelist_filename << "." << std::endl;
    std::terminate();
  }
  else
  {
    // Reads the file line by line.
    while( !edgelist_file.eof() )
    {
      // Reads a line of the file.
      std::getline(edgelist_file, full_line);
      edgelist_file >> std::ws;
      one_line.str(full_line);
      one_line >> std::ws;
      one_line >> name1_str >> std::ws;
      // Skips a line of comment.
      if(name1_str == "#")
      {
        one_line.clear();
        continue;
      }
      one_line >> name2_str >> std::ws;
      one_line.clear();

      // Adds the edge.
      add_edge(name1_str, name2_str, nb_vertices);
    }
  }
  // Closes the stream.
  edgelist_file.close();

  // ===============================================================================================
  // Updates the properties of the graph.
  g_prop["nb_vertices"] = nb_vertices;
  g_prop["nb_edges"] = edgelist.size();
  // ===============================================================================================
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::load_graph_from_edgelist_object(std::vector< std::pair<std::string, std::string> > edgelist_object)
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  Name2ID.clear();
  edgelist.clear();
  // ===============================================================================================


  // Initializes the number of vertices.
  int nb_vertices = 0;
  // String objects.
  std::string name1_str, name2_str;


  // Loops over the edgelist.
  for(int e(0), ee(edgelist_object.size()); e<ee; ++e)
  {
    // Reads an edge.
    name1_str = edgelist_object[e].first;
    name2_str = edgelist_object[e].second;

    // Adds the edge.
    add_edge(name1_str, name2_str, nb_vertices);
  }

  // ===============================================================================================
  // Updates the properties of the graph.
  g_prop["nb_vertices"] = nb_vertices;
  g_prop["nb_edges"] = edgelist.size();
  // ===============================================================================================
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::load_graph_from_edgelist_object(std::vector< std::pair<int, int> > edgelist_object)
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  Name2ID.clear();
  edgelist.clear();
  // ===============================================================================================


  // Stream objects.
  std::stringstream value;
  // Initializes the number of vertices.
  int nb_vertices = 0;
  // String objects.
  std::string name1_str, name2_str;


  // Loops over the edgelist.
  for(int e(0), ee(edgelist_object.size()); e<ee; ++e)
  {
    // Reads an edge.
    value.str("");
    value << edgelist_object[e].first;
    name1_str = value.str();
    value.str("");
    value << edgelist_object[e].second;
    name2_str = value.str();

    // Adds the edge.
    add_edge(name1_str, name2_str, nb_vertices);
  }

  // ===============================================================================================
  // Updates the properties of the graph.
  g_prop["nb_vertices"] = nb_vertices;
  g_prop["nb_edges"] = edgelist.size();
  // ===============================================================================================
}










// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// 3. OUTPUT
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::save_edgelist(std::string edgelist_filename, std::string id, bool header, int base)
{
  if(id == "names")
  {
    save_edgelist_using_names(edgelist_filename, header);
  }
  else if(id == "int")
  {
    save_edgelist_using_num(edgelist_filename, header, base);
  }
  else
  {
    std::cerr << "ERROR: Unknown vertex identifier type." << std::endl;
    std::terminate();
  }
}

// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::save_edgelist_using_names(std::string edgelist_filename, bool header)
{
  // Stream objects.
  std::fstream output_file;
  // Opens the stream and terminates if the operation did not succeed.
  output_file.open(edgelist_filename.c_str(), std::fstream::out);
  if( !output_file.is_open() )
  {
    std::cerr << "ERROR: Could not open file: " << edgelist_filename << "." << std::endl;
    std::terminate();
  }

  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  // ===============================================================================================

  // Ensures that the ID2Name has been built. Otherwise, build it.
  if(ID2Name.size() != nb_vertices)
  {
    build_ID2Name();
  }

  // Sets the width of columns.
  int width = 7;
  for(int v(0), l; v<nb_vertices; ++v)
  {
    l = ID2Name[v].length();
    if(l > width)
    {
      width = l;
    }
  }
  width += 2;

  // Variables.
  int v1, v2;
  edgelist_iterator it  = edgelist_begin();
  edgelist_iterator end = edgelist_end();
  // Writes the shuffled edgelist into the file.
  if(header)
  {
    output_file << "#" << std::setw(width - 1) << "Vertex1" << " ";
    output_file        << std::setw(width)     << "Vertex2" << " ";
    output_file << std::endl;
  }
  for(; it!=end; ++it)
  {
    v1 = it->first;
    v2 = it->second;
    if(v1 > v2)
    {
      std::swap(v1, v2);
    }
    output_file << std::setw(width) << ID2Name[v1] << " ";
    output_file << std::setw(width) << ID2Name[v2] << " ";
    output_file << std::endl;
  }
  // Closes the stream.
  output_file.close();
}

// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::save_edgelist_using_num(std::string edgelist_filename, bool header, int _base)
{
  // Stream objects.
  std::fstream output_file;
  // Opens the stream and terminates if the operation did not succeed.
  output_file.open(edgelist_filename.c_str(), std::fstream::out);
  if( !output_file.is_open() )
  {
    std::cerr << "ERROR: Could not open file: " << edgelist_filename << "." << std::endl;
    std::terminate();
  }

  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  // ===============================================================================================

  // Sets the width of columns.
  int width = std::log10(nb_vertices + _base);
  width = std::max( width, 7) + 2;

  // Variables.
  int v1, v2;
  edgelist_iterator it  = edgelist_begin();
  edgelist_iterator end = edgelist_end();
  // Writes the shuffled edgelist into the file.
  if(header)
  {
    output_file << "#" << std::setw(width - 1) << "Vertex1" << " ";
    output_file        << std::setw(width)     << "Vertex2" << " ";
    output_file << std::endl;
  }
  for(; it!=end; ++it)
  {
    v1 = it->first  + _base;
    v2 = it->second + _base;
    if(v1 > v2)
    {
      std::swap(v1, v2);
    }
    output_file << std::setw(width) << v1 << " ";
    output_file << std::setw(width) << v2 << " ";
    output_file << std::endl;
  }
  // Closes the stream.
  output_file.close();
}










// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// 4. DEGREE
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::average_degree()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  std::vector<double>& Vertex2Prop = v_prop["degrees"];
  // ===============================================================================================
  // Extracts the degrees of the vertices.
  if(Vertex2Prop.size() != nb_vertices)
  {
    degrees();
  }
  // Checks if the value of the degree correlation coefficient, if not already calculated.
  average_vertex_prop("degrees");
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::average_degree_of_neighbors(bool normalized)
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  std::vector<double>& Vertex2Degree = v_prop["degrees"];
  // ===============================================================================================

  // Extracts the degrees of the vertices.
  if(Vertex2Degree.size() != nb_vertices)
  {
    degrees();
  }

  // Computes the average degree of neighbors.
  average_vertex_prop_of_neighbors("degrees");

  if(normalized)
  {
    // =============================================================================================
    // Initializes relevant objects of the class.
    v_prop["avg_neigh_degrees_norm"].clear();
    std::vector<double>& Vertex2ANDegree = v_prop["avg_neigh_degrees_norm"];
    Vertex2ANDegree.resize(nb_vertices, 0);

    std::vector<double>& Vertex2ADegree = v_prop["avg_neigh_degrees"];
    // =============================================================================================

    // Computes the two first "moments" of the degree distribution.
    double m1 = 0;
    double m2 = 0;
    double d;
    for(int v(0); v<nb_vertices; ++v)
    {
      d = Vertex2Degree[v];
      m1 += d;
      m2 += d * d;
    }

    // Compute the normalized values of the average degree of neighbors.
    for(int v(0); v<nb_vertices; ++v)
    {
      Vertex2ANDegree[v] = Vertex2ADegree[v] * m1 / m2;
    }
  }
}


// =================================================================================================
// =================================================================================================
void undirected_graph_t::average_degree_of_neighbors_spectrum(bool normalized)
{
  std::string prop = "avg_neigh_degrees";
  if(normalized)
  {
    prop = "avg_neigh_degrees_norm";
  }
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  v_class_average_1p[std::make_pair("degrees", prop)].clear();
  // ===============================================================================================

  // Extracts the degrees of the vertices, if not already done.
  if(v_prop["degrees"].size() != nb_vertices)
  {
    degrees();
  }

  // Extracts the local clustering coefficient of the vertices, if not already done.
  if(v_prop[prop].size() != nb_vertices)
  {
    average_degree_of_neighbors(normalized);
  }

  // Counts the number of vertices in each degree class.
  average_vertex_prop_by_classes_of_vertices("degrees", prop);
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::degree_correlation_coefficient()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  int nb_edges = g_prop["nb_edges"];
  std::set<int>& DegreeClasses = v_class_1p["degrees"];
  // ===============================================================================================

  // Extracts the degrees of the vertices.
  if(v_prop["degrees"].size() != nb_vertices)
  {
    degrees();
  }

  // Computes the average degree of vertices.
  average_vertex_prop("degrees");

  // Compiles the degree distribution.
  identify_classes_of_vertices("degrees");
  count_vertices_in_classes("degrees");

  // Computes the excess degree distribution and its variance.
  int k;
  double avg_degree = g_prop["avg_degrees"];
  double var_q, value;
  double m1 = 0;
  double m2 = 0;
  std::map<int, double> q;
  // std::map<int, double> e;
  std::map<int,int>::iterator it1 = v_class_count_1p["degrees"].begin();
  std::map<int,int>::iterator end1 = v_class_count_1p["degrees"].end();
  for(; it1!=end1; ++it1)
  {
    k = it1->first;
    value = k * it1->second / avg_degree / nb_vertices;
    q[k - 1] = value;
    // e[k - 1] = 0;
    m2 += (k - 1) * (k - 1) * value;
    m1 += (k - 1) * value;
  }
  var_q = m2 - (m1 * m1);

  // Extracts the degree-degree correlation matrix.
  count_edges_in_classes("degrees");

  // Computes the degree correlation coefficient.
  double k1, k2;
  double corr_coef = 0;
  std::map< std::pair<int, int>, int>::iterator it2 = e_class_count_1p["degrees"].begin();
  std::map< std::pair<int, int>, int>::iterator end2 = e_class_count_1p["degrees"].end();
  for(; it2!=end2; ++it2)
  {
    k1 = it2->first.first;
    k2 = it2->first.second;
    corr_coef += (k1 - 1) * (k2 - 1) * it2->second / nb_edges;
  }
  std::map<int, double>::iterator it3 = q.begin();
  std::map<int, double>::iterator end3 = q.end();
  std::map<int, double>::iterator it4;
  std::map<int, double>::iterator end4 = q.end();
  for(; it3!=end3; ++it3)
  {
    it4 = it3;
    corr_coef -= it3->first * it4->first * it3->second * it4->second;
    for(++it4; it4!=end4; ++it4)
    {
      corr_coef -= 2 * it3->first * it4->first * it3->second * it4->second;
    }
  }

  // ===============================================================================================
  // Updates the properties of the graph.
  g_prop["degree_correlation_coefficient"] = corr_coef / var_q;
  // ===============================================================================================
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::degree_distribution()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];

  std::vector<double>& Vertex2Prop = v_prop["degrees"];

  v_stat["degree_dist"].clear();
  std::map<int, double>& VertStat = v_stat["degree_dist"];
  // ===============================================================================================

  // Extracts the degrees of the vertices, if not already done.
  if(Vertex2Prop.size() != nb_vertices)
  {
    degrees();
  }

  // Counts the number of vertices in each degree class.
  count_vertices_in_classes("degrees");
  // Normalizes the histogram.
  std::map<int, int>::iterator it  = v_class_count_1p["degrees"].begin();
  std::map<int, int>::iterator end = v_class_count_1p["degrees"].end();
  for(; it!=end; ++it)
  {
    VertStat[it->first] = (double) it->second / nb_vertices;
  }
  v_class_count_1p["degrees"].clear();
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::degrees()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  v_prop["degrees"].clear();
  std::vector<double>& Vertex2Degree = v_prop["degrees"];
  Vertex2Degree.resize(nb_vertices, 0);
  // ===============================================================================================


  // Extracts the degrees from the adjacency list (if already generated).
  if(adjacency_list.size() == nb_vertices)
  {
    // Loops over all vertices.
    for(int v(0); v<nb_vertices; ++v)
    {
      Vertex2Degree[v] = adjacency_list[v].size();
    }
  }
  // Otherwise extracts the degrees from the edgelist.
  else
  {
    // Loops over all edges.
    int v1, v2;
    std::set< std::pair<int, int> >::iterator it = edgelist.begin();
    std::set< std::pair<int, int> >::iterator end = edgelist.end();
    for(; it!=end; ++it)
    {
      // Identifies the vertices.
      v1 = it->first;
      v2 = it->second;
      // Adds the contribution to the degree of vertices.
      Vertex2Degree[v1] += 1;
      Vertex2Degree[v2] += 1;
    }
  }
}










// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// 5. CLUSTERING / TRIANGLES
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::average_local_clustering_coefficient()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  std::vector<double>& Vertex2Prop = v_prop["local_clust"];
  // ===============================================================================================
  // Extracts the degrees of the vertices.
  if(Vertex2Prop.size() != nb_vertices)
  {
    local_clustering_coefficients();
  }
  // Checks if the value of the degree correlation coefficient, if not already calculated.
  average_vertex_prop("local_clust");
}


// =================================================================================================
// =================================================================================================
void undirected_graph_t::clustering_spectrum()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  v_class_average_1p[std::make_pair("degrees", "local_clust")].clear();
  // ===============================================================================================

  // Extracts the degrees of the vertices, if not already done.
  if(v_prop["degrees"].size() != nb_vertices)
  {
    degrees();
  }

  // Extracts the local clustering coefficient of the vertices, if not already done.
  if(v_prop["local_clust"].size() != nb_vertices)
  {
    local_clustering_coefficients();
  }

  // Counts the number of vertices in each degree class.
  average_vertex_prop_by_classes_of_vertices("degrees", "local_clust");
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
int undirected_graph_t::count_triangles_around_vertex(int v1)
{
  // Variables.
  int v2, v3, d1, d2;
  int nb_triangles = 0;
  // Vector objects.
  std::vector<int> intersection;
  // Set objects.
  std::set<int> neighbours_v1, neighbours_v2;
  // Iterator objects.
  std::vector<int>::iterator it;
  // Degree of vertex v1.
  d1 = adjacency_list[v1].size();
  // Performs the calculation only if d1>1.
  if( d1 > 1 )
  {
    // Builds an ordered list of the neighbourhood of v1
    neighbours_v1.clear();
    neighbours_v1.insert(adjacency_list[v1].begin(), adjacency_list[v1].end());
    // Loops over the neighbours of vertex v1.
    for(int n1(0); n1<d1; ++n1)
    {
      // Identity and degree of vertex 2.
      v2 = adjacency_list[v1][n1];
      d2 = adjacency_list[v2].size();
      // Performs the calculation only if d2>1.
      if( d2 > 1 )
      {
        // Builds an ordered list of the neighbourhood of v2
        neighbours_v2.clear();
        for(int n2(0); n2<d2; ++n2)
        {
          // Identifies the neighbors.
          v3 = adjacency_list[v2][n2];
          if(v2 < v3) // Ensures that triangles will be counted only once.
          {
            neighbours_v2.insert(v3);
          }
        }
        // Identifies the triangles.
        d2 = neighbours_v2.size();
        intersection.clear();
        intersection.resize(std::min(d1, d2));
        it = std::set_intersection(neighbours_v1.begin(), neighbours_v1.end(), neighbours_v2.begin(), neighbours_v2.end(), intersection.begin());
        // intersection.resize(it-intersection.begin());
        // Computes the local clustering coeficient.
        // nb_triangles += intersection.size();
        nb_triangles += it-intersection.begin();
      }
    }
    // Returns the number of triangles.
    return nb_triangles;
  }
  else
  {
    return 0;
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::global_clustering_coefficient()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];

  std::vector<double>& Vertex2Degree = v_prop["degrees"];

  std::vector<double>& Vertex2NbTriangles = v_prop["nb_triangles"];
  // ===============================================================================================

  // Extracts the degree (if not already generated).
  if(Vertex2Degree.size() != nb_vertices)
  {
    degrees();
  }

  // Extracts the number of triangles to which each vertex contributes (if not already generated).
  if(Vertex2NbTriangles.size() != nb_vertices)
  {
    number_of_triangles();
  }

  // Computes the global clustering coefficient.
  double d;
  double global_number_of_triangles = 0;
  double global_number_of_wedges = 0;
  for(int v(0); v<nb_vertices; ++v)
  {
    d = Vertex2Degree[v];
    if(d > 1)
    {
      global_number_of_wedges += d * (d - 1) / 2;
    }
    global_number_of_triangles += Vertex2NbTriangles[v];
  }

  // ===============================================================================================
  // Updates the properties of the graph.
  g_prop["global_clustering_coefficient"] = global_number_of_triangles / global_number_of_wedges;
  // ===============================================================================================
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::local_clustering_coefficients()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];

  v_prop["local_clust"].clear();
  std::vector<double>& Vertex2Prop = v_prop["local_clust"];
  Vertex2Prop.resize(nb_vertices, 0);
  // ===============================================================================================

  // Ensures that the adjacency list has been built.
  if(adjacency_list.size() != nb_vertices)
  {
    build_adjacency_list();
  }

  // Computes the intersection for the in- and out- neighbourhoods of each node.
  int d1;
  double nb_triangles, tmp_value;
  double tmp_average_local_clustering_coefficient = 0;
  for(int v1(0); v1<nb_vertices; ++v1)
  {
    // The degree of the vertex.
    d1 = adjacency_list[v1].size();
    // Counts the number of triangle.
    nb_triangles = count_triangles_around_vertex(v1);
    tmp_value = 0;
    // Computes the coefficient of clustering.
    if(d1 > 1)
    {
      tmp_value = 2.0 * nb_triangles / d1 / (d1 - 1);
      Vertex2Prop[v1] = tmp_value;
      tmp_average_local_clustering_coefficient += tmp_value;
    }
    else
    {
      Vertex2Prop[v1] = 0;
    }
  }
  // ===============================================================================================
  // Updates the properties of the graph.
  g_prop["avg_local_clust"] = tmp_average_local_clustering_coefficient / g_prop["nb_vertices"];
  // ===============================================================================================
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::number_of_triangles()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  v_prop["nb_triangles"].clear();
  std::vector<double>& Vertex2Prop = v_prop["nb_triangles"];
  Vertex2Prop.resize(nb_vertices, 0);
  // ===============================================================================================


  // Ensures that the adjacency list has been built.
  if(adjacency_list.size() != nb_vertices)
  {
    build_adjacency_list();
  }

  // Counts the number of triangles to which each vertex participates.
  int tmp_value;
  int nb_triangles = 0;
  for(int v1(0); v1<nb_vertices; ++v1)
  {
    tmp_value = count_triangles_around_vertex(v1);
    Vertex2Prop[v1] = tmp_value;
    nb_triangles += tmp_value;
  }

  // ===============================================================================================
  // Updates the properties of the graph.
  g_prop["nb_triangles"] = nb_triangles / 3;
  // ===============================================================================================
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::survey_triangles()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  triangles.clear();
  // ===============================================================================================

  // Variables.
  int v2, v3, d1, d2;
  // Vector objects.
  std::vector<int> intersection;
  // Set objects.
  std::set<int> neighbours_v1, neighbours_v2;
  // Iterator objects.
  std::vector<int>::iterator it;
  // Newly found triangle.
  std::vector<int> triangle(3);

  // Ensures that the adjacency list has been built.
  if(adjacency_list.size() != nb_vertices)
  {
    build_adjacency_list();
  }

  // Finds all the triangles.
  for(int v1(0); v1<nb_vertices; ++v1)
  {
    // Degree of vertex v1.
    d1 = adjacency_list[v1].size();
    // Performs the calculation only if d1>1.
    if( d1 > 1 )
    {
      // Builds an ordered list of the neighbourhood of v1
      neighbours_v1.clear();
      neighbours_v1.insert(adjacency_list[v1].begin(), adjacency_list[v1].end());
      // Loops over the neighbours of vertex v1.
      for(int n1(0); n1<d1; ++n1)
      {
        // Identity and degree of vertex 2.
        v2 = adjacency_list[v1][n1];
        d2 = adjacency_list[v2].size();
        // Performs the calculation only if d2>1 and if v2>v1 (ensures that each triangle is counted once).
        if( v1 < v2 && d2 > 1 )
        {
          // Builds an ordered list of the neighbourhood of v2
          neighbours_v2.clear();
          for(int n2(0); n2<d2; ++n2)
          {
            // Identifies the neighbors.
            v3 = adjacency_list[v2][n2];
            if(v2 < v3) // Ensures that triangles will be counted only once.
            {
              neighbours_v2.insert(v3);
            }
          }
          // Identifies the triangles.
          d2 = neighbours_v2.size();
          intersection.clear();
          intersection.resize(std::min(d1, d2));
          it = std::set_intersection(neighbours_v1.begin(), neighbours_v1.end(), neighbours_v2.begin(), neighbours_v2.end(), intersection.begin());
          intersection.resize(it-intersection.begin());
          // Loops over the common neighbours of vertices v1 and v2.
          for(int n(0), nn(intersection.size()); n<nn; ++n)
          {
            // Adds the triangle to the list.
            triangle[0] = v1;
            triangle[1] = v2;
            triangle[2] = intersection[n];
            triangles.push_back(triangle);
          }
        }
      }
    }
  }

  // ===============================================================================================
  // Updates the properties of the graph.
  g_prop["nb_triangles"] = triangles.size();
  // ===============================================================================================
}










// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// 6. CENTRALITY
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::closeness_centrality(bool correct_for_component_size)
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  std::string v_prop_name;
  if(correct_for_component_size)
  {
    v_prop_name = "corr_close_cent";
  }
  else
  {
    v_prop_name = "close_cent";
  }
  std::vector<double>& Vertex2Cent = v_prop[v_prop_name];
  Vertex2Cent.clear();
  Vertex2Cent.resize(nb_vertices, 0);
  std::vector<double>& Vertex2Comp = v_prop["conn_comp"];
  // ===============================================================================================

  // Ensures that the shortest path lengths have been computed.
  if(shortest_path_lengths.size() != nb_vertices)
  {
    survey_shortest_path_lengths();
  }

  // Ensures that the connected components have been surveyed.
  if(connected_components_size.size() != nb_vertices)
  {
    survey_connected_components();
  }

  // Compiles the closeness centrality of every vertex (ignores paths that does not exist).
  int length, cnt;
  double value;
  for(int v1(0); v1<nb_vertices; ++v1)
  {
    cnt = connected_components_size[ Vertex2Comp[v1] ];
    value = 0;
    for(int v2(0); v2<nb_vertices; ++v2)
    {
      length = shortest_path_length(v1, v2);
      if(length != nb_vertices)
      {
        value += length;
        // ++cnt;
      }
    }
    if(correct_for_component_size)
    {
      // There is also the Wasserman and Faust correction.
      Vertex2Cent[v1] = ((cnt - 1) / ((double) (nb_vertices - 1))) * ((cnt - 1) / value);
    }
    else
    {
      Vertex2Cent[v1] = (cnt - 1) / value;
    }
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::harmonic_centrality()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  std::vector<double>& Vertex2Cent = v_prop["harmo_cent"];
  Vertex2Cent.clear();
  Vertex2Cent.resize(nb_vertices, 0);
  // ===============================================================================================


  // Ensures that the shortest path lengths have been computed.
  if(shortest_path_lengths.size() != nb_vertices)
  {
    survey_shortest_path_lengths();
  }

  // Compiles the harmonic centrality of every vertex.
  int length;
  double value;
  for(int v1(0); v1<nb_vertices; ++v1)
  {
    value = 0;
    for(int v2(0); v2<nb_vertices; ++v2)
    {
      length = shortest_path_length(v1, v2);
      if(length != nb_vertices && length != 0)
      {
        value += 1 / (double) length;
      }
    }
    Vertex2Cent[v1] = value / (nb_vertices - 1);
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::kcore_decomposition()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];

  v_prop["kcore"].clear();
  std::vector<double>& Vertex2Coreness = v_prop["kcore"];
  Vertex2Coreness.resize(nb_vertices, 0);

  std::vector<double>& Vertex2Degree = v_prop["degrees"];
  // ===============================================================================================

  // Ensures that the adjacency list has been built.
  if(adjacency_list.size() != nb_vertices)
  {
    build_adjacency_list();
  }

  // Extracts the degree (if not already generated).
  if(Vertex2Degree.size() != nb_vertices)
  {
    degrees();
  }

  // Builds two lists (std::vector, std::set) of the degree of the vertices.
  std::vector<int> DegreeVec(nb_vertices);
  std::set<std::pair<int, int> > DegreeSet;
  for(int v(0); v<nb_vertices; ++v)
  {
    DegreeSet.insert(std::make_pair(Vertex2Degree[v], v));
    DegreeVec[v] = Vertex2Degree[v];
  }

  // Determines the coreness and the layer based on the algorithm of Batagelj and Zaversnik.
  int v1, v2, d1, d2;
  std::set< std::pair<int, int> >::iterator m_it;
  while(DegreeSet.size() > 0)
  {
    // Sets the coreness of the first vertex in the list and then removes it.
    m_it = DegreeSet.begin();
    d1 = m_it->first;
    v1 = m_it->second;
    Vertex2Coreness[v1] = d1;
    DegreeSet.erase(m_it);
    // Reduces the "effective" degree of its neighbours.
    for(int n(0), nn(adjacency_list[v1].size()); n<nn; ++n)
    {
      // Identifies the neighbor.
      v2 = adjacency_list[v1][n];
      d2 = DegreeVec[v2];
      // Finds the neighbor in the list "effective" degrees.
      m_it = DegreeSet.find(std::make_pair(d2, v2));
      if(m_it != DegreeSet.end())
      {
        if(d2 > d1)
        {
          DegreeVec[v2] = d2 - 1;
          DegreeSet.erase(m_it);
          DegreeSet.insert(std::make_pair(d2 - 1, v2));
        }
      }
    }
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::onion_decomposition()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];

  v_prop["kcore"].clear();
  std::vector<double>& Vertex2Coreness = v_prop["kcore"];
  Vertex2Coreness.resize(nb_vertices, 0);

  v_prop["od_layer"].clear();
  std::vector<double>& Vertex2Layer = v_prop["od_layer"];
  Vertex2Layer.resize(nb_vertices, 0);

  std::vector<double>& Vertex2Degree = v_prop["degrees"];
  // ===============================================================================================

  // Ensures that the adjacency list has been built.
  if(adjacency_list.size() != nb_vertices)
  {
    build_adjacency_list();
  }

  // Extracts the degree (if not already generated).
  if(Vertex2Degree.size() != nb_vertices)
  {
    degrees();
  }

  // Builds two lists (std::vector, std::set) of the degree of the vertices.
  std::vector<int> DegreeVec(nb_vertices);
  std::set<std::pair<int, int> > DegreeSet;
  for(int v(0); v<nb_vertices; ++v)
  {
    DegreeSet.insert(std::make_pair(Vertex2Degree[v], v));
    DegreeVec[v] = Vertex2Degree[v];
  }


  // Determines the coreness and the layer based on the modified algorithm of Batagelj and
  //   Zaversnik by Hébert-Dufresne, Grochow and Allard.
  int v1, v2, d1, d2;
  int layer = 0;
  std::set< std::pair<int, int> > LayerSet;
  std::set< std::pair<int, int> >::iterator m_it;
  while(DegreeSet.size() > 0)
  {
    // Increases the layer id.
    layer += 1;
    // Populates the set containing the vertices belonging to the same layer.
    m_it = DegreeSet.begin();
    d1 = m_it->first;
    // Sets the coreness and the layer the vertices with the same degree.
    while(m_it->first == d1 && m_it != DegreeSet.end())
    {
      // Sets the coreness and the layer.
      v1 = m_it->second;
      Vertex2Coreness[v1] = d1;
      Vertex2Layer[v1] = layer;
      // Looks at the next vertex.
      ++m_it;
    }
    // Adds the vertices of the layer to the set.
    LayerSet.insert(DegreeSet.begin(), m_it);
    // Removes the vertices of the current layer.
    DegreeSet.erase(DegreeSet.begin(), m_it);
    // Modifies the "effective" degree of the neighbors of the vertices in the layer.
    while(LayerSet.size() > 0)
    {
      // Gets information about the next vertex of the layer.
      v1 = LayerSet.begin()->second;
      // Reduces the "effective" degree of its neighbours.
      for(int n(0), nn(adjacency_list[v1].size()); n<nn; ++n)
      {
        // Identifies the neighbor.
        v2 = adjacency_list[v1][n];
        d2 = DegreeVec[v2];
        // Finds the neighbor in the list "effective" degrees.
        m_it = DegreeSet.find(std::make_pair(d2, v2));
        if(m_it != DegreeSet.end())
        {
          if(d2 > d1)
          {
            DegreeVec[v2] = d2 - 1;
            DegreeSet.erase(m_it);
            DegreeSet.insert(std::make_pair(d2 - 1, v2));
          }
        }
      }
      // Removes the vertices from the LayerSet.
      LayerSet.erase(LayerSet.begin());
    }
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::onion_spectrum()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];

  std::vector<double>& Vertex2Onion = v_prop["od_layer"];
  std::vector<double>& Vertex2kcore = v_prop["kcore"];

  onion_spect.clear();
  // ===============================================================================================

  // Extracts the degrees of the vertices, if not already done.
  if(Vertex2Onion.size() != nb_vertices)
  {
    onion_decomposition();
  }

  // Number of layers.
  int nb_layers = *std::max_element(Vertex2Onion.begin(), Vertex2Onion.end());

  // Mapping onion layers to k-shells and building the unsegregated onion spectrum.
  std::vector<int> od_layer2k_shell(nb_layers + 1);
  std::vector<double> unsegregated_onion_spectrum(nb_layers + 1, 0);
  for(int v(0), o; v<nb_vertices; ++v)
  {
    o = Vertex2Onion[v];
    od_layer2k_shell[o] = Vertex2kcore[v];
    unsegregated_onion_spectrum[o] += 1;
  }

  // Segregates the onion spectrum.
  for(int o(1); o<=nb_layers; ++o)
  {
    onion_spect[ od_layer2k_shell[o] ].push_back(unsegregated_onion_spectrum[o] / nb_vertices);
  }
}










// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// 7. CONNECTIVITY
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
int undirected_graph_t::get_root(int i, std::vector<int> &clust_id)
{
  while(i != clust_id[i])
  {
    clust_id[i] = clust_id[clust_id[i]];
    i = clust_id[i];
  }
  return i;
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::merge_clusters(std::vector<int> &size, std::vector<int> &clust_id, int nb_vertices)
{
  // Variables.
  int v1, v2, v3, v4;
  // Loops over the vertices.
  for(int i(0); i<nb_vertices; ++i)
  {
    // Loops over the neighbors.
    for(int j(0), jj(adjacency_list[i].size()); j<jj; ++j)
    {
      if(get_root(i, clust_id) != get_root(adjacency_list[i][j], clust_id))
      {
        // Adjust the root of vertices.
        v1 = i;
        v2 = adjacency_list[i][j];
        if(size[v2] > size[v1])
          std::swap(v1, v2);
        v3 = get_root(v1, clust_id);
        v4 = get_root(v2, clust_id);
        clust_id[v4] = v3;
        size[v3] += size[v4];
      }
    }
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::survey_connected_components()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];

  connected_components_size.clear();
  ordered_connected_components.clear();

  v_prop["conn_comp"].clear();
  std::vector<double>& Vertex2Prop = v_prop["conn_comp"];
  Vertex2Prop.resize(nb_vertices, 0);
  // ===============================================================================================

  // Ensures that the adjacency list has been built.
  if(adjacency_list.size() != nb_vertices)
  {
    build_adjacency_list();
  }

  // Starts with every vertex as an isolated cluster.
  std::vector<int> clust_id(nb_vertices);
  std::vector<int> clust_size(nb_vertices, 1);
  for(int v(0); v<nb_vertices; ++v)
  {
    clust_id[v] = v;
  }
  // Merges clusters until the minimal set is obtained.
  merge_clusters(clust_size, clust_id, nb_vertices);
  clust_size.clear();
  // Identifies the connected component to which each vertex belongs.
  int nb_conn_comp = 0;
  int comp_id;
  std::map<int, int> CompID;
  for(int v(0); v<nb_vertices; ++v)
  {
    comp_id = get_root(v, clust_id);
    if(CompID.find(comp_id) == CompID.end())
    {
      CompID[comp_id] = nb_conn_comp;
      connected_components_size.push_back(0);
      ++nb_conn_comp;
    }
    Vertex2Prop[v] = CompID[comp_id];
    connected_components_size[CompID[comp_id]] += 1;
  }

  // Orders the size of the components.
  for(int c(0); c<nb_conn_comp; ++c)
  {
    ordered_connected_components.insert( std::make_pair(connected_components_size[c], c) );
  }

  // ===============================================================================================
  // Updates the properties of the graph.
  g_prop["nb_conn_comp"] = nb_conn_comp;
  g_prop["lcc_id"] = (--ordered_connected_components.end())->second;
  g_prop["lcc_size"] = (--ordered_connected_components.end())->first;
  g_prop["second_lcc_id"] = (----ordered_connected_components.end())->second;
  g_prop["second_lcc_size"] = (----ordered_connected_components.end())->first;
  // ===============================================================================================
}










// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// 8. SHORTEST PATHS
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
int undirected_graph_t::shortest_path_length(int v1, int v2)
{
  // Ensures that the shortest path lengths have been computed.
  if(shortest_path_lengths.size() != g_prop["nb_vertices"])
  {
    survey_shortest_path_lengths();
  }
  if(v2 > v1)
  {
    return shortest_path_lengths[v2][v1];
  }
  else
  {
    return shortest_path_lengths[v1][v2];
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::survey_shortest_path_length_distribution()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  int nb_conn_comp = g_prop["nb_conn_comp"];
  std::vector<double>& Vertex2Comp = v_prop["conn_comp"];
  // ===============================================================================================


  // Ensures that the shortest path lengths have been computed.
  if(shortest_path_lengths.size() != nb_vertices)
  {
    survey_shortest_path_lengths();
  }

  // Ensures that the connected components have been identified.
  if(Vertex2Comp.size() != nb_vertices)
  {
    survey_connected_components();
    nb_conn_comp = g_prop["nb_conn_comp"];
  }

  // Ensures that the average shortest path lengths have been computed.
  if(diameters.size() != nb_conn_comp)
  {
    topological_dimensions();
  }

  // Resets the container.
  shortest_path_length_distribution.clear();
  shortest_path_length_distribution.resize(nb_conn_comp);
  for(int c(0); c<nb_conn_comp; ++c)
  {
    shortest_path_length_distribution[c].resize(diameters[c] + 1, 0);
  }

  // Computes the shortest path length distribution for each connected component.
  int c;
  double l;
  std::vector<int> comp_size(nb_conn_comp, 0);
  for(int i(0), ii(shortest_path_lengths.size()); i<ii; ++i)
  {
    c = Vertex2Comp[i];
    for(int j(0), jj(shortest_path_lengths[i].size()-1); j<jj; ++j)
    {
      l = shortest_path_lengths[i][j];
      if(l < nb_vertices)
      {
        comp_size[c] += 1;
        shortest_path_length_distribution[c][l] += 1;
      }
    }
  }
  for(int c(0); c<nb_conn_comp; ++c)
  {
    for(int l(0), ll(shortest_path_length_distribution[c].size()); l<ll; ++l)
    {
      shortest_path_length_distribution[c][l] /= comp_size[c];
    }
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::survey_shortest_path_lengths()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  std::vector<double>& Vertex2Degree = v_prop["degrees"];
  shortest_path_lengths.clear();
  // ===============================================================================================


  // Ensures that the adjacency list has been built.
  if(adjacency_list.size() != nb_vertices)
  {
    build_adjacency_list();
  }

  // Extracts the degree (if not already generated).
  if(Vertex2Degree.size() != nb_vertices)
  {
    degrees();
  }

  // Resizes the container for the shortest path lengths.
  shortest_path_lengths.resize(nb_vertices);
  for(int v(0); v<nb_vertices; ++v)
  {
    shortest_path_lengths[v].resize(v+1, nb_vertices);
  }

  // Finds the length of shortest paths via a breadth first search.
  int v1, v2;
  std::vector<int> to_visit;
  to_visit.resize(nb_vertices);
  int next_current, last_current, last_next, current_distance;
  bool keep_going;
  std::vector<bool> is_new(nb_vertices);
  for(int v(0); v<nb_vertices; ++v)
  {
    // No vertex has been visited yet, except "v".
    std::fill(is_new.begin(), is_new.end(), true);
    is_new[v] = false;
    to_visit[0] = v;

    next_current = 0;
    last_current = 0;
    last_next = 0;
    current_distance = 0;
    keep_going = true;
    while(keep_going)
    {
      // Sets the shortest path distance lengths for vertices in the current layer and populates the
      //   next layer.
      for(; next_current <= last_current; ++next_current)
      {
        // Identifies the vertex that has just been reached.
        v1 = to_visit[next_current];

        // Sets the shortest path length.
        if(v > v1)
          shortest_path_lengths[v][v1] = current_distance;
        else
          shortest_path_lengths[v1][v] = current_distance;

        // Loops over the neighbors of the vertex that has just been reached.
        for(int s(0), ss(Vertex2Degree[v1]); s<ss; ++s)
        {
          // Identifies the neighbor.
          v2 = adjacency_list[v1][s];
          // Adds the vertex to the next layer if it has not been reached yet.
          if(is_new[v2])
          {
            is_new[v2] = false;
            ++last_next;
            to_visit[last_next] = v2;
          }
        }
      }

      // Checks if all vertices have been reached.
      if(last_current == last_next)
      {
        keep_going = false;
      }
      else
      {
        last_current = last_next;
        ++current_distance;
      }
    }
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void undirected_graph_t::topological_dimensions()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  int nb_conn_comp = g_prop["nb_conn_comp"];
  std::vector<double>& Vertex2Comp = v_prop["conn_comp"];
  // ===============================================================================================


  // Ensures that the shortest path lengths have been computed.
  if(shortest_path_lengths.size() != nb_vertices)
  {
    survey_shortest_path_lengths();
  }

  // Ensures that the connected components have been identified.
  if(Vertex2Comp.size() != nb_vertices)
  {
    survey_connected_components();
    nb_conn_comp = g_prop["nb_conn_comp"];
  }

  // Resizes the containers.
  average_shortest_path_lengths.clear();
  average_shortest_path_lengths.resize(nb_conn_comp, 0);
  diameters.clear();
  diameters.resize(nb_conn_comp, 0);

  // Computes the average length.
  int c;
  double l;
  std::vector<int> comp_size(nb_conn_comp, 0);
  double avg_value = 0;
  for(int i(0), ii(shortest_path_lengths.size()); i<ii; ++i)
  {
    c = Vertex2Comp[i];
    for(int j(0), jj(shortest_path_lengths[i].size()-1); j<jj; ++j)
    {
      l = shortest_path_lengths[i][j];
      if(l < nb_vertices)
      {
        comp_size[c] += 1;
        average_shortest_path_lengths[c] += l;
        if(l > diameters[c])
        {
          diameters[c] = l;
        }
      }
    }
  }
  // Completes the calculation of the average values.
  for(int s(0); s<nb_conn_comp; ++s)
  {
    average_shortest_path_lengths[s] /= comp_size[s];
  }
}










// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// 9. ACCESSORS FOR GRAPH PROPERTIES
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
double undirected_graph_t::get_average_degree()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  // ===============================================================================================
  // Extracts the property.
  if(v_prop["degrees"].size() != nb_vertices)
  {
    degrees();
  }
  // Checks if the property has already been calculated. Otherwise, calculates it.
  if(g_prop.find("avg_degrees") == g_prop.end())
  {
    average_vertex_prop("degrees");
  }
  return g_prop["avg_degrees"];
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
double undirected_graph_t::get_average_local_clustering_coefficient()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  // ===============================================================================================
  // Extracts the property.
  if(v_prop["local_clust"].size() != nb_vertices)
  {
    local_clustering_coefficients();
  }
  // Checks if the property has already been calculated. Otherwise, calculates it.
  if(g_prop.find("avg_local_clust") == g_prop.end())
  {
    average_vertex_prop("local_clust");
  }
  return g_prop["avg_local_clust"];
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
std::vector<double> undirected_graph_t::get_average_shortest_path()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  int nb_conn_comp = g_prop["nb_conn_comp"];
  std::vector<double>& Vertex2Comp = v_prop["conn_comp"];
  // ===============================================================================================

  // Ensures that the shortest path lengths have been computed.
  if(shortest_path_lengths.size() != nb_vertices)
  {
    survey_shortest_path_lengths();
  }
  // Ensures that the connected components have been identified.
  if(Vertex2Comp.size() != nb_vertices)
  {
    survey_connected_components();
    nb_conn_comp = g_prop["nb_conn_comp"];
  }
  // Ensures that the average shortest path lengths have been computed.
  if(average_shortest_path_lengths.size() != g_prop["nb_conn_comp"])
  {
    topological_dimensions();
  }

  // // Gets the proper order of the components.
  std::set< std::pair<int, int> >::iterator it = ordered_connected_components.begin();
  std::set< std::pair<int, int> >::iterator end = ordered_connected_components.end();
  std::vector<int> comp_order(g_prop["nb_conn_comp"]);
  for(int c(g_prop["nb_conn_comp"] - 1); it!=end; ++it, --c)
  {
    comp_order[it->second] = c;
  }

  // Initializes the output container.
  std::vector<double> array(g_prop["nb_conn_comp"]);
  for(int c(0), cc(g_prop["nb_conn_comp"]); c<cc; ++c)
  {
    array[ comp_order[c] ] = average_shortest_path_lengths[c];
  }
  return array;
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
std::vector< std::vector<std::string> > undirected_graph_t::get_connected_components()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  std::vector<double>& Vertex2Prop = v_prop["conn_comp"];
  // ===============================================================================================
  // Checks if the property has already been calculated. Otherwise, calculates it.
  if(Vertex2Prop.size() != nb_vertices)
  {
    survey_connected_components();
  }
  // Ensures that the ID2Name has been built. Otherwise, build it.
  if(ID2Name.size() != nb_vertices)
  {
    build_ID2Name();
  }
  // // Gets the proper order of the components.
  std::set< std::pair<int, int> >::iterator it = ordered_connected_components.begin();
  std::set< std::pair<int, int> >::iterator end = ordered_connected_components.end();
  std::vector<int> comp_order(g_prop["nb_conn_comp"]);
  for(int c(g_prop["nb_conn_comp"] - 1); it!=end; ++it, --c)
  {
    comp_order[it->second] = c;
  }
  // Initializes the output container.
  std::vector< std::vector<std::string> > connected_components(g_prop["nb_conn_comp"]);
  it = ordered_connected_components.begin();
  end = ordered_connected_components.end();
  for(; it!=end; ++it)
  {
    connected_components[ comp_order[it->second] ].reserve(it->first);
  }
  // Extracts the members of the components in order of the size of the components.
  for(int v(0), c; v<nb_vertices; ++v)
  {
    c = Vertex2Prop[v];
    connected_components[ comp_order[c] ].push_back( ID2Name[v] );
  }
  return connected_components;
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
std::vector<int> undirected_graph_t::get_connected_components_size()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  int nb_conn_comp = g_prop["nb_conn_comp"];
  std::vector<double>& Vertex2Comp = v_prop["conn_comp"];
  // ===============================================================================================

  // Ensures that the connected components have been identified.
  if(Vertex2Comp.size() != nb_vertices)
  {
    survey_connected_components();
    nb_conn_comp = g_prop["nb_conn_comp"];
  }

  // Gets the proper order of the components.
  std::set< std::pair<int, int> >::iterator it = ordered_connected_components.begin();
  std::set< std::pair<int, int> >::iterator end = ordered_connected_components.end();
  std::vector<int> comp_order(g_prop["nb_conn_comp"]);
  for(int c(g_prop["nb_conn_comp"] - 1); it!=end; ++it, --c)
  {
    comp_order[it->second] = c;
  }

  // Initializes the output container.
  it = ordered_connected_components.begin();
  end = ordered_connected_components.end();
  std::vector<int> array(g_prop["nb_conn_comp"]);
  for(; it!=end; ++it)
  {
    array[ comp_order[it->second] ] = it->first;
  }
  return array;
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
std::vector<int> undirected_graph_t::get_diameter()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  int nb_conn_comp = g_prop["nb_conn_comp"];
  std::vector<double>& Vertex2Comp = v_prop["conn_comp"];
  // ===============================================================================================

  // Ensures that the shortest path lengths have been computed.
  if(shortest_path_lengths.size() != nb_vertices)
  {
    survey_shortest_path_lengths();
  }
  // Ensures that the connected components have been identified.
  if(Vertex2Comp.size() != nb_vertices)
  {
    survey_connected_components();
    nb_conn_comp = g_prop["nb_conn_comp"];
  }
  // Ensures that the average shortest path lengths have been computed.
  if(diameters.size() != g_prop["nb_conn_comp"])
  {
    topological_dimensions();
  }

  // // Gets the proper order of the components.
  std::set< std::pair<int, int> >::iterator it = ordered_connected_components.begin();
  std::set< std::pair<int, int> >::iterator end = ordered_connected_components.end();
  std::vector<int> comp_order(g_prop["nb_conn_comp"]);
  for(int c(g_prop["nb_conn_comp"] - 1); it!=end; ++it, --c)
  {
    comp_order[it->second] = c;
  }

  // Initializes the output container.
  std::vector<int> array(g_prop["nb_conn_comp"]);
  for(int c(0), cc(g_prop["nb_conn_comp"]); c<cc; ++c)
  {
    array[ comp_order[c] ] = diameters[c];
  }
  return array;
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
double undirected_graph_t::get_degree_correlation_coefficient()
{
  // Checks if the property has already been calculated. Otherwise, calculates it.
  if(g_prop.find("degree_correlation_coefficient") == g_prop.end())
  {
    degree_correlation_coefficient();
  }
  return g_prop["degree_correlation_coefficient"];
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
double undirected_graph_t::get_global_clustering_coefficient()
{
  // Checks if the property has already been calculated. Otherwise, calculates it.
  if(g_prop.find("global_clustering_coefficient") == g_prop.end())
  {
    global_clustering_coefficient();
  }
  return g_prop["global_clustering_coefficient"];
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
int undirected_graph_t::get_nb_connected_components()
{
  // Checks if the property has already been calculated. Otherwise, calculates it.
  if(g_prop.find("nb_conn_comp") == g_prop.end())
  {
    survey_connected_components();
  }
  return g_prop["nb_conn_comp"];
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
int undirected_graph_t::get_nb_edges()
{
  return g_prop["nb_edges"];
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
int undirected_graph_t::get_nb_vertices()
{
  return g_prop["nb_vertices"];
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
int undirected_graph_t::get_shortest_path(std::string v1, std::string v2)
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  // ===============================================================================================
  // Checks if the statistics has already been calculated. Otherwise, calculates it.
  if(shortest_path_lengths.size() != nb_vertices)
  {
    survey_shortest_path_lengths();
  }
  return shortest_path_length(Name2ID[v1], Name2ID[v2]);
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
std::vector< std::vector<double> > undirected_graph_t::get_shortest_path_length_distribution()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  int nb_conn_comp = g_prop["nb_conn_comp"];
  // ===============================================================================================
  // Checks if the property has already been calculated. Otherwise, calculates it.
  // Ensures that the connected components have been identified.
  if(v_prop["conn_comp"].size() != nb_vertices)
  {
    survey_connected_components();
    nb_conn_comp = g_prop["nb_conn_comp"];
  }
  if(shortest_path_length_distribution.size() != nb_conn_comp)
  {
    survey_shortest_path_length_distribution();
  }

  // Gets the proper order of the components.
  std::set< std::pair<int, int> >::iterator it = ordered_connected_components.begin();
  std::set< std::pair<int, int> >::iterator end = ordered_connected_components.end();
  std::vector<int> comp_order(nb_conn_comp);
  for(int c(nb_conn_comp - 1); it!=end; ++it, --c)
  {
    comp_order[it->second] = c;
  }
  // Initializes the output container.
  std::vector< std::vector<double> > path_length_distribution(nb_conn_comp);
  it = ordered_connected_components.begin();
  end = ordered_connected_components.end();
  for(; it!=end; ++it)
  {
    path_length_distribution[ comp_order[it->second] ] = shortest_path_length_distribution[it->second];
  }
  return path_length_distribution;
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
int undirected_graph_t::get_total_nb_triangles()
{
  // Checks if the property has already been calculated. Otherwise, calculates it.
  if(g_prop.find("nb_triangles") == g_prop.end())
  {
    number_of_triangles();
  }
  return g_prop["nb_triangles"];
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
std::vector< std::vector< std::string > > undirected_graph_t::get_triangles()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  // ===============================================================================================
  // Ensures that the ID2Name has been built. Otherwise, build it.
  if(ID2Name.size() != nb_vertices)
  {
    build_ID2Name();
  }
  // Extracts the property.
  if(triangles.size() == 0)
  {
    survey_triangles();
  }
  // Builds the map from names to property.
  std::vector< std::vector< std::string > > array(triangles.size(), std::vector<std::string>(3));
  for(int v(0), vv(triangles.size()); v<vv; ++v)
  {
    array[v][0] = ID2Name[triangles[v][0]];
    array[v][1] = ID2Name[triangles[v][1]];
    array[v][2] = ID2Name[triangles[v][2]];
  }
  return array;
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
std::vector< std::pair<std::string, std::string> > undirected_graph_t::get_edges()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  int nb_edges = g_prop["nb_edges"];
  // ===============================================================================================
  // Ensures that the ID2Name has been built. Otherwise, build it.
  if(ID2Name.size() != nb_vertices)
  {
    build_ID2Name();
  }
  // Builds the list of edges.
  std::vector< std::pair<std::string, std::string> > array(nb_edges);
  edgelist_iterator it = edgelist_begin();
  edgelist_iterator end = edgelist_end();
  for(int e(0); it!=end; ++it, ++e)
  {
    array[e] = std::make_pair(ID2Name[it->first], ID2Name[it->second]);
  }
  return array;
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
std::vector< std::string > undirected_graph_t::get_vertices()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  // ===============================================================================================
  // Ensures that the ID2Name has been built. Otherwise, build it.
  if(ID2Name.size() != nb_vertices)
  {
    build_ID2Name();
  }
  // Builds the list of vertices.
  std::vector< std::string > array(nb_vertices);
  for(int v(0); v<nb_vertices; ++v)
  {
    array[v] = ID2Name[v];
    // array[v] = ID2Name[v];
    // array[v] = ID2Name[v];
  }
  return array;
}










// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// 10. ACCESSORS FOR VERTEX STATISTICS
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
std::map<int, double> undirected_graph_t::get_average_degree_of_neighbors_spectrum(bool normalized)
{
  std::string prop = "avg_neigh_degrees";
  if(normalized)
  {
    prop = "avg_neigh_degrees_norm";
  }
  // Checks if the statistics has already been calculated. Otherwise, calculates it.
  if(v_class_average_1p.find(std::make_pair("degrees", prop)) == v_class_average_1p.end())
  {
    average_degree_of_neighbors_spectrum(normalized);
  }
  return v_class_average_1p[std::make_pair("degrees", prop)];
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
std::map<int, double> undirected_graph_t::get_clustering_spectrum()
{
  // Checks if the statistics has already been calculated. Otherwise, calculates it.
  if(v_class_average_1p.find(std::make_pair("degrees", "local_clust")) == v_class_average_1p.end())
  {
    clustering_spectrum();
  }
  return v_class_average_1p[std::make_pair("degrees", "local_clust")];
}


// // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// std::map<int, int> undirected_graph_t::get_count_vertices_in_classes(std::string prop)
// {
//   if(v_class_count_1p.find(prop) == v_class_count_1p.end())
//   {
//     count_vertices_in_classes(prop);
//   }
//   return v_class_count_1p[prop];
// }


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
std::map<int, double> undirected_graph_t::get_degree_distribution()
{
  // Checks if the statistics has already been calculated. Otherwise, calculates it.
  if(v_stat.find("degree_dist") == v_stat.end())
  {
    degree_distribution();
  }
  return v_stat["degree_dist"];
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
std::map<int, std::vector<double> > undirected_graph_t::get_onion_spectrum()
{
  // Checks if the statistics has already been calculated. Otherwise, calculates it.
  if(onion_spect.size() == 0)
  {
    onion_spectrum();
  }
  return onion_spect;
}










// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// 11. ACCESSORS FOR VERTEX PROPERTIES
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
std::map<std::string, double> undirected_graph_t::get_average_degree_of_neighbors(bool normalized)
{
  std::string prop = "avg_neigh_degrees";
  if(normalized)
  {
    prop = "avg_neigh_degrees_norm";
  }
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  std::vector<double>& Vertex2Prop = v_prop[prop];
  // ===============================================================================================
  // Ensures that the ID2Name has been built. Otherwise, build it.
  if(ID2Name.size() != nb_vertices)
  {
    build_ID2Name();
  }
  // Extracts the property.
  if(Vertex2Prop.size() != nb_vertices)
  {
    average_degree_of_neighbors(normalized);
  }
  // Builds the map from names to property.
  std::map<std::string, double> dict;
  for(int v(0); v<nb_vertices; ++v)
  {
    dict[ID2Name[v]] = Vertex2Prop[v];
  }
  return dict;
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
std::map<std::string, double> undirected_graph_t::get_closeness_centrality(bool correct_for_component_size)
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  std::string v_prop_name;
  if(correct_for_component_size)
  {
    v_prop_name = "corr_close_cent";
  }
  else
  {
    v_prop_name = "close_cent";
  }
  std::vector<double>& Vertex2Prop = v_prop[v_prop_name];
  // ===============================================================================================
  // Ensures that the ID2Name has been built. Otherwise, build it.
  if(ID2Name.size() != nb_vertices)
  {
    build_ID2Name();
  }
  // Extracts the property.
  if(Vertex2Prop.size() != nb_vertices)
  {
    closeness_centrality(correct_for_component_size);
  }
  // Builds the map from names to property.
  std::map<std::string, double> dict;
  for(int v(0); v<nb_vertices; ++v)
  {
    dict[ID2Name[v]] = Vertex2Prop[v];
  }
  return dict;
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
std::map<std::string, double> undirected_graph_t::get_degrees()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  std::vector<double>& Vertex2Prop = v_prop["degrees"];
  // ===============================================================================================
  // Ensures that the ID2Name has been built. Otherwise, build it.
  if(ID2Name.size() != nb_vertices)
  {
    build_ID2Name();
  }
  // Extracts the property.
  if(Vertex2Prop.size() != nb_vertices)
  {
    degrees();
  }
  // Builds the map from names to property.
  std::map<std::string, double> dict;
  for(int v(0); v<nb_vertices; ++v)
  {
    dict[ID2Name[v]] = Vertex2Prop[v];
  }
  return dict;
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
std::map<std::string, double> undirected_graph_t::get_k_shells()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  std::vector<double>& Vertex2Prop = v_prop["kcore"];
  // ===============================================================================================
  // Ensures that the ID2Name has been built. Otherwise, build it.
  if(ID2Name.size() != nb_vertices)
  {
    build_ID2Name();
  }
  // Extracts the property.
  if(Vertex2Prop.size() != nb_vertices)
  {
    kcore_decomposition();
  }
  // Builds the map from names to property.
  std::map<std::string, double> dict;
  for(int v(0); v<nb_vertices; ++v)
  {
    dict[ID2Name[v]] = Vertex2Prop[v];
  }
  return dict;
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
std::map<std::string, double> undirected_graph_t::get_harmonic_centrality()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  std::vector<double>& Vertex2Prop = v_prop["harmo_cent"];
  // ===============================================================================================
  // Ensures that the ID2Name has been built. Otherwise, build it.
  if(ID2Name.size() != nb_vertices)
  {
    build_ID2Name();
  }
  // Extracts the property.
  if(Vertex2Prop.size() != nb_vertices)
  {
    harmonic_centrality();
  }
  // Builds the map from names to property.
  std::map<std::string, double> dict;
  for(int v(0); v<nb_vertices; ++v)
  {
    dict[ID2Name[v]] = Vertex2Prop[v];
  }
  return dict;
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
std::map<std::string, double> undirected_graph_t::get_local_clustering_coefficients()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  std::vector<double>& Vertex2Prop = v_prop["local_clust"];
  // ===============================================================================================
  // Ensures that the ID2Name has been built. Otherwise, build it.
  if(ID2Name.size() != nb_vertices)
  {
    build_ID2Name();
  }
  // Extracts the property.
  if(Vertex2Prop.size() != nb_vertices)
  {
    local_clustering_coefficients();
  }
  // Builds the map from names to property.
  std::map<std::string, double> dict;
  for(int v(0); v<nb_vertices; ++v)
  {
    dict[ID2Name[v]] = Vertex2Prop[v];
  }
  return dict;
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
std::map<std::string, double> undirected_graph_t::get_nb_triangles()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  std::vector<double>& Vertex2Prop = v_prop["nb_triangles"];
  // ===============================================================================================
  // Ensures that the ID2Name has been built. Otherwise, build it.
  if(ID2Name.size() != nb_vertices)
  {
    build_ID2Name();
  }
  // Extracts the property.
  if(Vertex2Prop.size() != nb_vertices)
  {
    number_of_triangles();
  }
  // Builds the map from names to property.
  std::map<std::string, double> dict;
  for(int v(0); v<nb_vertices; ++v)
  {
    dict[ID2Name[v]] = Vertex2Prop[v];
  }
  return dict;
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
std::map<std::string, double> undirected_graph_t::get_od_layers()
{
  // ===============================================================================================
  // Initializes relevant objects of the class.
  int nb_vertices = g_prop["nb_vertices"];
  std::vector<double>& Vertex2Prop = v_prop["od_layer"];
  // ===============================================================================================
  // Ensures that the ID2Name has been built. Otherwise, build it.
  if(ID2Name.size() != nb_vertices)
  {
    build_ID2Name();
  }
  // Extracts the property.
  if(Vertex2Prop.size() != nb_vertices)
  {
    onion_decomposition();
  }
  // Builds the map from names to property.
  std::map<std::string, double> dict;
  for(int v(0); v<nb_vertices; ++v)
  {
    dict[ID2Name[v]] = Vertex2Prop[v];
  }
  return dict;
}










#endif
