#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <cmath>
#include <ctime>
#include <map>
#include <random>
#include <algorithm>
#include <string>
#include <deque>


std::random_device dev;
std::mt19937 rng(dev());
std::uniform_real_distribution<double> dist(0,1); // distribution in range [0, 1)

std::random_device dev1;
std::mt19937 rng1(dev1());
std::uniform_int_distribution<int> dist1(0,1); // distribution in range [0, 1]



class Lattice
{
  // The lattice class.  This houses the link variables
  // as well as the neighbors and connectivity info.
  
  // vector containers 
  std::map<int, std::vector<std::vector<int> > > staples;
  std::vector<int> plaq;
  std::vector<int> current_staple;
  std::vector<int> spin_draw = {-1, 1};
  
  // variables int
  int num_plaqs;
  int current_spin;
  int spin_prod;
  int staple_sum;
  int current_action;
  int prop_action;
  int delta_action;
  int idx0;
  
  // variables doubles
  double probability;
  double beta;
  double r;
  double total;
    
  
 

public:
  // might want these
  int num_links;
  std::vector<int> edges;
  double acpt;
  double acpt_count;
  double total_count;
  
  Lattice(std::vector<std::vector<int> > &plaq_vec,
	  double beta_in, bool, bool);
  void metro_update(int edge_idx);
  void hb_update(int edge_idx);
  void save_config(std::string fname);
  void load_config(std::string fname);
  
  
};

Lattice::Lattice(std::vector<std::vector<int> > &plaq_vec,
		 double beta_in, bool hot=true, bool load=false)
{
  num_plaqs = plaq_vec.size();
  current_staple.resize(3);
  beta = beta_in;
  std::cout << "beta " << beta << std::endl;
  
  // build the staple map.  for each edge there is a vector
  // of staple vectors.
  for (int i = 0; i < num_plaqs; ++i)
    {
      plaq = plaq_vec[i];         
      staples[plaq[0]].push_back({plaq[1], plaq[2], plaq[3]});
      staples[plaq[1]].push_back({plaq[2], plaq[3], plaq[0]});
      staples[plaq[2]].push_back({plaq[3], plaq[0], plaq[1]});
      staples[plaq[3]].push_back({plaq[0], plaq[1], plaq[2]});
    }
  
  num_links = staples.size();
  // std::cout << num_links << std::endl;
  edges.resize(num_links);
  
  if (load)
    {
      // do nothing 
    }
  else
    {
      // the spins on the edges
      for (int i = 0; i < num_links; ++i)
        {
	  if (hot)
            {
	      idx0 = dist1(rng1);
	      edges[i] = spin_draw[idx0];
            }
	  else
            {
	      edges[i] = 1;   
            }   
        }   
    }
  
  
}


void Lattice::metro_update(int edge_idx)
{
  // The metropolis update
  current_spin = edges[edge_idx];
  
  staple_sum = 0;
  for (std::vector<std::vector<int> >::iterator it = staples[edge_idx].begin();
       it != staples[edge_idx].end(); ++it)
    {
      spin_prod = 1;
      current_staple = *it;
      for (int j = 0; j < 3; ++j)
        {
	  spin_prod *= edges[current_staple[j]];
        }
      staple_sum += spin_prod;
      // std::cout << "stape sum " << staple_sum << std::endl;
    }
  
  current_action = -1 * current_spin * staple_sum / 4;
  prop_action = -1 * current_action;
  delta_action = prop_action - current_action;
  // std::cout << delta_action << std::endl;
  
  if (delta_action < 0)
    {
      // std::cout << "taken" << std::endl;
      edges[edge_idx] = -1 * current_spin;
      acpt_count += 1;
    }
  else
    {
      r = dist(rng);
      probability = exp(-1 * beta * delta_action);
      // std::cout <<"val " << beta  << std::endl;
      // std::cout <<"prob " << probability << std::endl;
      if (r <= probability)
        {
	  edges[edge_idx] = -1 * current_spin;
	  acpt_count += 1;
        }
    }
  total_count += 1;
  acpt = acpt_count / total_count;
}

void Lattice::hb_update(int edge_idx)
{
  // the heatbath update
  current_spin = edges[edge_idx];
  
  staple_sum = 0;
  for (std::vector<std::vector<int> >::iterator it = staples[edge_idx].begin();
       it != staples[edge_idx].end(); ++it)
    {
      spin_prod = 1;
      current_staple = *it;
      for (int j = 0; j < 3; ++j)
        {
	  spin_prod *= edges[current_staple[j]];
        }
      staple_sum += spin_prod;
      // std::cout << "stape sum " << staple_sum << std::endl;
    }

  total = exp(-1 * beta * staple_sum / 4) +
    exp(beta * staple_sum / 4);
  probability = exp(beta * staple_sum / 4) / total;
  r = dist(rng);
  if (r <= probability)
    {
      edges[edge_idx] = 1;
    }
  else
    {
      edges[edge_idx] = -1;
    }
}


void Lattice::load_config(std::string fname)
{
  std::ifstream file(fname);
  std::string line;
  int k = 0;
  int v;
  while (std::getline(file, line))
    {
      std::istringstream ss(line);
      while (ss >> v)
        {
	  edges[k] = v;
        }
      k += 1;
    }
}

void Lattice::save_config(std::string fname)
{
  std::ofstream output_file(fname);
  std::ostream_iterator<int> output_iterator(output_file, "\n");
  std::copy(edges.begin(), edges.end(), output_iterator);
}
