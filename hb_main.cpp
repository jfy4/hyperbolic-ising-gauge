#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <map>
#include <random>
#include <algorithm>
#include <string>
#include <deque>
#include "lattice.h"



int main(int argc, char const *argv[])
{
  
  double beta = atof(argv[1]);
  bool load = false;
  std::string sbeta;
  std::string outname;
  std::string L = "18";
  
  // empty plaquette vector
  std::vector<std::vector<int> > plaquettes;
  std::vector<int> temp_plaq(4);
  // import plaquettes here
	// read in the flag file
  std::ifstream file("L" + L + "_facet_edges.txt");
  std::string line;
  while (std::getline(file, line))
    {
      std::istringstream ss(line);
      int v;
      int k = 0;
      while (ss >> v)
	{
	  temp_plaq[k] = v-1;
	  k += 1;
	}
      plaquettes.push_back(temp_plaq);
    }
  
  int nmeas = 100;
  int sweep_count;
  int nprint = 100;
  
  Lattice lat(plaquettes, beta, load=load);
  
  if (load)
    {
      lat.load_config("./hb_b0p7/hb_b0.700000_l5_c10500.txt");
      sweep_count = 10500;
      sweep_count += 1;
    }
  else
    {
      sweep_count = 0;
    }
  // return 0;
  
  while (true)
    {
      // lat.total_count = 0;
      // lat.acpt_count = 0;
      // lat.acpt = 0;
      for (int i = 0; i < lat.num_links; ++i)
	{
	  lat.hb_update(i);
	  // std:: cout << "acceptance = " << lat.acpt << std::endl;
	}
      if (sweep_count % nprint == 0)
	{
	  std:: cout << "sweep = " << sweep_count  << std::endl;
	}
      if (sweep_count % nmeas == 0)	
	{
	  // std::string beta_string = std::to_string(beta);
	  sbeta = std::to_string(beta).substr(0,3);
	  outname = "./hb_b" + sbeta.replace(1, 1, "p") + "/L" + L + "/hb_b" + sbeta + "_l" + L +
	    "_c" + std::to_string(sweep_count) + ".txt";
	  // std::cout << outname << std::endl;
	  lat.save_config(outname);
	}
      sweep_count += 1;
      // break;
    }
  // for (int i = 0; i < lat.num_links; ++i)
  // {
  // 	std::cout << lat.edges[i] << std::endl;
  // }
  
  return 0;
}
