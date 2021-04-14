// Standard Template Library
#include <chrono>
#include <cmath>
#include <numeric>
#include <random>
#include <string>
#include <vector>
// Portable Graph Library
#include "undirected_graph_t.hpp"


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
int get_root(int i, std::vector<int> &clust_id)
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
void merge_clusters(std::vector<int> &size,
                    std::vector<int> &clust_id,
                    std::vector< std::vector<int> > &adjacency_list,
                    int nb_vertices)
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
std::vector<bool> percolate(undirected_graph_t& the_graph,
                            int observation_depth_L,
                            double app_coverage,
                            double adoption_of_passports,
                            std::mt19937& engine,
                            std::uniform_real_distribution<double>& uniform_01,
                            std::discrete_distribution<int>& RandPrivacyStatus1,
                            std::discrete_distribution<int>& RandPrivacyStatus2)
{
  std::vector<double> &Vertex2PrivacyStatus1 = the_graph.v_prop["PrivacyStatus1"];
  std::vector<double> &Vertex2PrivacyStatus2 = the_graph.v_prop["PrivacyStatus2"];
  int nb_vertices = the_graph.g_prop["nb_vertices"];
  std::vector<bool> is_vertex_observed(nb_vertices, false);
  // std::vector<bool> is_vertex_observed_directly_by_the_app(nb_vertices, false);
  // std::vector<bool> is_vertex_in_second_layer(nb_vertices, false);
  std::set<int> vertices_in_current_observability_layer;
  std::set<int> vertices_in_next_observability_layer;
  // Assigns the types;
  for(int v(0); v<nb_vertices; ++v)
  {
    Vertex2PrivacyStatus1[v] = RandPrivacyStatus1(engine);
    if(Vertex2PrivacyStatus1[v] == 2 && uniform_01(engine) < adoption_of_passports) Vertex2PrivacyStatus2[v] = 2;
    else Vertex2PrivacyStatus2[v] = RandPrivacyStatus2(engine);

    if(Vertex2PrivacyStatus1[v] == 0 && uniform_01(engine) < app_coverage) //test for platform 1
    {
      is_vertex_observed[v] = true;
      // is_vertex_observed_directly_by_the_app[v] = true;
      vertices_in_current_observability_layer.insert(v);
    }
    else if(Vertex2PrivacyStatus2[v] == 0 && uniform_01(engine) < app_coverage) //test for platform 2
    {
      is_vertex_observed[v] = true;
      // is_vertex_observed_directly_by_the_app[v] = true;
      vertices_in_current_observability_layer.insert(v);
    }
  }

  for(int l(0), t11, t12, t21, t22, v2; l<observation_depth_L; ++l)
  {
    vertices_in_next_observability_layer.clear();
    for(auto v1 : vertices_in_current_observability_layer)
    {
      // v1 = *obs_v_it;
      t11 = Vertex2PrivacyStatus1[v1];
      t12 = Vertex2PrivacyStatus2[v1];
      for(int i(0), ii(the_graph.adjacency_list[v1].size()); i<ii; ++i)
      {
        v2 = the_graph.adjacency_list[v1][i];
        t21 = Vertex2PrivacyStatus1[v2];
        t22 = Vertex2PrivacyStatus2[v2];
        // if(observed_vertices[v2] == false)
        // {
          if(t21 <= (t11+1) || t22 <= (t12+1))
          {
            is_vertex_observed[v2] = true;
            vertices_in_next_observability_layer.insert(v2);
          }
        // }
      }
    }
    vertices_in_current_observability_layer.clear();
    vertices_in_current_observability_layer.insert(vertices_in_next_observability_layer.begin(),
                                                   vertices_in_next_observability_layer.end());
  }

  return is_vertex_observed;
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
std::vector< std::vector<int> > extract_adjlist(undirected_graph_t& the_graph,
                                                std::vector<bool>& is_vertex_observed,
                                                bool observable_component)
{
  // Creates the unobservable adjacency list.
  int nb_vertices = the_graph.g_prop["nb_vertices"];
  std::vector< std::vector<int> > adjlist(nb_vertices);
  int v1, v2;
  // int t1, t2;
  for(auto edge : the_graph.edgelist)
  {
    // Indentifies the first vertex.
    v1 = edge.first;
    // Indentifies the second vertex.
    v2 = edge.second;

    // Keeps the edge if none of the vertices are observed...
    if((is_vertex_observed[v1] == observable_component) && (is_vertex_observed[v2] == observable_component))
    {
      // // ...and if neither of them are of type 0. COULD BE COMMENTED
      // if((t1 > 0) && (t2 > 0))
      {
        adjlist[v1].push_back(v2);
        adjlist[v2].push_back(v1);
      }
    }
  }
  return adjlist;
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
std::vector< std::vector<int> > analyse_components(undirected_graph_t& the_graph,
                                                   std::vector< std::vector<int> >& adjlist,
                                                   std::vector<bool>& is_vertex_observed,
                                                   bool observable_component)
{
  // Find unobservable clusters
  // Starts with every vertex as an isolated cluster.
  int nb_vertices = the_graph.g_prop["nb_vertices"];
  std::vector<int> clust_id(nb_vertices);
  std::iota(clust_id.begin(), clust_id.end(), 0);
  std::vector<int> clust_size(nb_vertices, 1);
  // Merges clusters until the minimal set is obtained.
  merge_clusters(clust_size, clust_id, adjlist, nb_vertices);

  // Sizes of the connected components.
  std::vector< std::vector<int> > connected_components_composition;

  int nb_conn_comp = 0;
  int comp_id;
  int t;
  std::map<int, int> CompID;
  for(int v(0); v<nb_vertices; ++v)
  {
    if(is_vertex_observed[v] == observable_component)
    {
      comp_id = get_root(v, clust_id);
      if(CompID.find(comp_id) == CompID.end())
      {
        CompID[comp_id] = nb_conn_comp;
        connected_components_composition.push_back(std::vector<int>(3, 0));
        ++nb_conn_comp;
      }
      t = the_graph.v_prop["PrivacyStatus1"][v];
      connected_components_composition[CompID[comp_id]][t] += 1;
    }
  }
  return connected_components_composition;
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
std::vector<int>  get_components_composition(undirected_graph_t& the_graph,
                                             std::vector< std::vector<int> >& adjlist,
                                             std::vector<bool>& is_vertex_observed,
                                             bool observable_component)
{
  auto connected_components_composition = analyse_components(the_graph, adjlist, is_vertex_observed, observable_component);
  std::vector<int> comp(7, 0);
  comp[0] = connected_components_composition.size();
  int largest_component_size = 0;
  int component_size;
  for(int c(0), cc(connected_components_composition.size()); c<cc; ++c)
  {
    comp[1] += connected_components_composition[c][0];
    comp[2] += connected_components_composition[c][1];
    comp[3] += connected_components_composition[c][2];
    component_size = 0;
    component_size += connected_components_composition[c][0];
    component_size += connected_components_composition[c][1];
    component_size += connected_components_composition[c][2];
    if(component_size > largest_component_size)
    {
      largest_component_size = component_size;
      comp[4] = connected_components_composition[c][0];
      comp[5] = connected_components_composition[c][1];
      comp[6] = connected_components_composition[c][2];
    }
  }

  return comp;
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
int main(int argc, char const *argv[])
{

  // - connected component of counter-culture
  // - herd immunity? Probability for a type 1 to be observed compared to when adoption rate equals 0.
  // - what rate of adoption can prevent "cambridge analytica"-like app to successfully access

  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  // Parameters specific to the model.
  double app_coverage = std::stod(argv[3]);

  // Observation depth (number of hops) of the app (ex.: 1 gives access to the first neighbors).
  int observation_depth_L = std::stoi(argv[2]);

  // Parametrization of privacy types.
  int number_of_privacy_statuses = 3;
  double fraction_of_private_profiles = std::stod(argv[4]);
  double adoption_within_private_profiles = std::stod(argv[5]);
  double adoption_of_passports = std::stod(argv[6]);
  std::vector<double> prob_types(number_of_privacy_statuses);
    prob_types[0] = 1.0 - fraction_of_private_profiles;
    prob_types[1] = fraction_of_private_profiles * (1 - adoption_within_private_profiles);
    prob_types[2] = fraction_of_private_profiles * adoption_within_private_profiles;
  std::vector<double> prob_types2(number_of_privacy_statuses);
    double renorm = (3.0-adoption_within_private_profiles*adoption_of_passports)/3.0;
    prob_types2[0] = (1.0 - fraction_of_private_profiles)/renorm;
    prob_types2[1] = fraction_of_private_profiles * (1 - adoption_within_private_profiles) / renorm;
    prob_types2[2] = fraction_of_private_profiles * adoption_within_private_profiles * (1 - adoption_of_passports) / renorm;


  // if(adoption < 1.0/3.0)
  // {
  //   prob_types[0] = 1 - fraction_of_private_profiles;
  //   prob_types[1] = fraction_of_private_profiles - adoption;
  //   prob_types[2] = adoption;
  // }
  // else
  // {
  //   prob_types[0] = 1 - fraction_of_private_profiles - (adoption-1.0/3.0);
  //   prob_types[1] = 0.0;
  //   prob_types[2] = adoption;
  // }

  // Number of simulations.
  int nb_simulations = std::stoi(argv[7]);


  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  // Loads and sets up the graph.
  std::string edgelist_filename = argv[1];
  undirected_graph_t the_graph(edgelist_filename);
  the_graph.build_adjacency_list();
  int nb_vertices = the_graph.g_prop["nb_vertices"];

  // Creates the vertex property related to the type of vertices.
  the_graph.new_vertex_integer_property("PrivacyStatus1");
  the_graph.new_vertex_integer_property("PrivacyStatus2");
  std::vector<double> &Vertex2PrivacyStatus1 = the_graph.v_prop["PrivacyStatus1"];
  std::vector<double> &Vertex2PrivacyStatus2 = the_graph.v_prop["PrivacyStatus2"];


  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  // Random numbers.
  unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::mt19937 engine(seed);
  std::uniform_real_distribution<double> uniform_01;
  std::discrete_distribution<int> RandPrivacyStatus1(prob_types.begin(), prob_types.end());
  std::discrete_distribution<int> RandPrivacyStatus2(prob_types2.begin(), prob_types2.end());


  // =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  std::vector< std::vector<int> > obs_adjlist, nonobs_adjlist;
  std::vector<bool> is_vertex_observed;
  std::vector<int> obs_comp, nonobs_comp;
  for(int j(0); j<nb_simulations; ++j)
  {

    // Identifies the vertices that are observed (directly or indirectly).
    is_vertex_observed = percolate(the_graph, observation_depth_L, app_coverage, adoption_of_passports, engine, uniform_01, RandPrivacyStatus1, RandPrivacyStatus2);

    // Extracts the observable adjacency list.
    obs_adjlist = extract_adjlist(the_graph, is_vertex_observed, true);
    obs_comp = get_components_composition(the_graph, obs_adjlist, is_vertex_observed, true);

    // Extracts the nonobservable adjacency list.
    nonobs_adjlist = extract_adjlist(the_graph, is_vertex_observed, false);
    nonobs_comp = get_components_composition(the_graph, nonobs_adjlist, is_vertex_observed, false);




    // std::cout << edgelist_filename << " ";
    // std::cout << observation_dept h_L << " ";
    // std::cout << number_of_privacy_statuses << " ";
    std::cout << std::setw(15) << argv[8] << " ";
    std::cout << std::setw(15) << observation_depth_L << " ";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << std::setw(15) << app_coverage << " ";
    std::cout << std::setw(15) << fraction_of_private_profiles << " ";
    std::cout << std::setw(15) << adoption_of_passports <<" ";
    std::cout << std::setw(15) << adoption_within_private_profiles <<" ";
    std::cout << std::fixed << std::setprecision(0);
    std::cout << std::setw(15) << nb_vertices << " ";
    for(int n(0); n<7; ++n)
    {
      std::cout << std::setw(15) << obs_comp[n] << " ";
    }
    for(int n(0); n<7; ++n)
    {
      std::cout << std::setw(15) << nonobs_comp[n] << " ";
    }
    std::cout << std::endl;

  } //end for loop over simulations

  return 0;
}
