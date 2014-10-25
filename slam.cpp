// SLaM: a forward population-genetic sumulation program with selection, linkage, migration, and
  // and deme-specific fitness regimes. SLaM stands for "SLiM with local adaptation and migration",
  // where "SLiM" is the software written by Philipp Messer (2013) and means "Selection on Linked
  // Mutations". SLaM was developed based on SLiM version 1.7, (C) 2013 Philipp Messer.
// version 1.0 (Month xxth, 2014)
//
// Copyright (C) 201x Simon Aeschbacher
//
// compile by:
//
// g++ -fast ./slam.cpp -lgsl -lgslcblas -o slam // iMac
// g++ -O3   ./slam.cpp -lgsl -lgslcblas -o slam // Linux
//
// You may need to specify the path to the gsl library and include path explicitly. If so, use:
// g++ -fast ./slam.cpp -L/usr/local/lib -I/usr/local/include -lgsl -lgslcblas -o slam
// Mac mini
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version (http://www.gnu.org/licenses/).
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// Modifications by Simon Aeschbacher, as compared to the original SLiM code (version 1.7):
// 05/19/2014: Mutations are now introduced in linkage equilibrium, rather than in linkage
  // disequilibrium. This concerns the method introduce_mutation() of the class 'population'.
  //
// 08/19/2014: Changed from multiplicative to additive fitness interactions among loci. This
  // concerns the method W() of the class 'subpopulation'.
  //
// xx/xx/2014: Added the possibility to define environment specific fitnesses. To this purpose, the
  // class 'environment' was introduced. Each subpopulation must be assigned a environment. For
  // each environment the user specifies how the fitnesses parameters (selection and dominance
  // coefficients) are to be modified relative to the reference environment. The reference
  // environment is implicitly given by the initial description of the mutation types, and does not
  // need to be redefined. There is also a new type of event, called "changing environment", which
  // allows assigning a environment to any subpopulation. By default, any population assumes the
  // reference environment.
  // Syntax for input spec:
  //          #ENVIRONMENTS
  //          <environment-id> <mut-type> <h> <s-modif> [<mut-type>
  //          <h> <s-modif> ...]
  //
// xx/xx/2014: Introduced a switch that allows the user to choose whether predetermined mutations
  // introduced at a given point in time should initially be in linkage disequilibrium (as it was
  // originally the case in SLiM) or at linkage equilibrium (new). The switch is added as an
  // additional entry on each line in the #PREDETERMINED MUTATIONS section, before the optional
  // parameter P <f>. The switch for initial linkage (dis)equilibrium is "e" for linkage
  // equilibrium and "d" for linkage disequilibrium.
  // Syntax:
  //          #PREDETERMINED MUTATIONS
  //          <time> <mut-type> <x> <pop-id> <nAA> <nAa> <linkage> [P <f>]
  //
// xx/xx/2014: Introduced a swith that allows the user to choose between multiplicative fitness
  // interactions between loci (as in the original SLiM) and additive fitness interactions (new).
  // This is a global setting, and has a new input section called #FITNESS INTERACTION that
  // contains one line with one entry. This entry must be "a" for additive or "m" for
  // multiplicative fitness interactions.
  // Syntax:
  //          #FITNESS INTERACTION
  //          <interaction-type>
  //
// xx/xx/2014: Changed the mutation regime at non-neutral sites, such that at most two mutations of
  // a non-neutral type are permitted at such sites. There still is not back-mutation. Neutral
  // sites can accumulate more than one neutral mutation. However, as soon as a site is hit by
  // the first mutation of a non-neutral type, no further mutations of a non-neutral type are
  // allowed. If, however, the non-neutral mutation at a given site is lost from the population,
  // the site will be considered as a neutral site again, and can accumulate at most one mutation
  // of a non-neutral type.
  //
// xx/xx/2014: Reformatted the code so that it follows more or less the GNU indent style, where it
  // did not do so before.

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <set>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include <unistd.h>


using namespace std;

const gsl_rng *rng; 


class mutation
{
public:

  int   t; // mutation type identifier
  int   x; // position
  float s; // selection coefficient
    // TODO: selection coefficients (and dominance coefficients) must become environment-specific.
    // Introduce class 'environment', and then assign a environment to each deme. environment
    // assignments can change over time. Change s to be an array of length equal to the number of
    // environment types.

  // basic constructor
  mutation(void) { ; }

  // alternative constructor
  mutation(int T, int X, float S) 
  { 
    t = T;
    x = X;
    s = S; // TODO: Will become an array of floats
  }
}; // end of class mutation


bool operator< (const mutation &M1, const mutation &M2)
{
  return M1.x < M2.x;
};


bool operator== (const mutation &M1, const mutation &M2) //
{
  return (M1.x == M2.x && M1.t == M2.t && M1.s == M2.s);
};


class event
{
  // type of events:
  //
  // t P i n [j]:  add subpopulation i of size n [drawn from j]
  // t N i n:      set size of subpopulation i to n
  // t M i j x:    set fraction x of subpopulation i that originates as migrants from j
  // t S i s:      set selfing fraction of subpopulation i to s
  // t E i e:      assign environment e to subpopulation i at time t // TODO: Implement this.
  //
  // t R i n:      output sample of n randomly drawn genomes from subpopulation i
  // t F:          output list of all mutations that have become fixed so far
  // t A [file]:   output state of entire population [into file]
  // t T m:        follow trajectory of mutation m (specified by mutation type) from generation t
  //               on

public:
  
  char   t;         // event type
  vector<string> s; // vector of strings with parameters of event
  int np;           // number of parameters

  event(char T, vector<string> S)
    {
      t = T;
      s = S;
      np = s.size();

      string options = "PNMSERFAT";
      if (options.find(t) == string::npos)
        {
          cerr << "ERROR (initialize): invalid event type \"" << t;
          for (int i = 0; i < np; i++)
            { cerr << " " << s[i]; }
          cerr << "\"" << endl;
          exit(1);
        }
    }
}; // end of class event




class mutation_type
{
  // a mutation type is specified by the DFE and the dominance coefficient. These parameters
  // implicitly define the reference environment.
  //
  // DFE options: f: fixed (s) 
  //              e: exponential (mean s)
  //              g: gamma distribution (mean s,shape)
  //
  // examples: synonymous, nonsynonymous, adaptive, etc.

public:
    
  float  h;         // dominance coefficient 
  char   d;         // DFE (f: fixed, g: gamma, e: exponential)
  vector<double> p; // DFE parameters

  mutation_type(float H, char D, vector<double> P)
  {
    h = H;
    d = D;
    p = P;

    string s = "fge";

    if (s.find(d) == string::npos)
      {
        cerr << "ERROR (initialize): invalid mutation type parameters" << endl; exit(1);
      }
    if (p.size() == 0)
      {
        cerr << "ERROR (initialize): invalid mutation type parameters" << endl; exit(1);
      }
  } // end of constructor

  float draw_s()
  {
    switch(d)
      {
      case 'f': return p[0]; 
      case 'g': return gsl_ran_gamma(rng,p[1],p[0]/p[1]); // If the mean is negative,
              // gsl_ran_gamma returns negative values.
      case 'e': return gsl_ran_exponential(rng,p[0]); // Note that if the mean is negative,
              // gsl_ran_exponential returns negative values.
      default: exit(1);
      }
  } // end of draw_s()
}; // end of class mutation_type

class environment
{
  // an environment reassigns a dominance coefficient and modifies the (mean) selection coefficient
    // with respect to the reference environment and for a specified set of mutation types

public:

  map< int,vector<double> > h; // dominance coefficients for each mutation-type that is affected in
    // this environment relative to the reference environment; the key is the mutation-type
    // identifier
  map< int,vector<double> > smodif; // modifier of (mean) selection coefficient for each mutation-
    // type that is affected in this environment relative to the reference environment; the key is
    // the mutation-type identifier

  // constructor
  environment(map< int,vector<double> >& H, map< int,vector<double> >& SMODIF)
  {

    h = H;
    smodif = SMODIF;

    if ( h.size() != smodif.size())
      { exit(1); }

  } // end of constructor

}; // end of class environment

class genomic_element
{
  // a genomic element has a genomic element type identifier (i), start (s) and end (e) position

public:

  int i, s, e;

  // constructor
  genomic_element(int I, int S, int E)
    { i = I; s = S; e = E; }

}; // end of class genomic_element


class genomic_element_type
{
  // a genomic element type is specified by a vector of the mutation type identifiers off all 
  // mutation types than can occur in such elements and a vector of their relative fractions.
  // examples: exon, intron, utr, intergenic, etc.

private:

  gsl_ran_discrete_t* LT;

public:

  vector<int>    m; // mutation types identifiers in this element
  vector<double> g; // relative fractions of each mutation type

  genomic_element_type(vector<int> M, vector<double> G)
  {
    m = M;
    g = G;  

    if (m.size() != g.size())
      { exit(1); }
    double A[m.size()];
    for (int i = 0; i < m.size(); i++)
      { A[i] = g[i]; }
    LT = gsl_ran_discrete_preproc(G.size(), A);
  }

  int draw_mutation_type()
    { return m[gsl_ran_discrete(rng, LT)]; }
}; // end of class genomic_element_type



class chromosome : public vector<genomic_element> // class chromosome inherits from
    // vector<genomic_element>
{
  // the chromosome is a vector of genomic elements (type, start, end)

private:

  // look-up tables

  gsl_ran_discrete_t* LT_M; // mutation
  gsl_ran_discrete_t* LT_R; // recombination

public:

  map<int,mutation_type>        mutation_types;
  map<int,genomic_element_type> genomic_element_types;
  vector<int>                   rec_x; // vector of end points of strata of a given recomb. rate
  vector<double>                rec_r; // vector of recombination rates for each stratum

  int    L;   // length of chromosome
  double M;   // overall mutation rate
  double R;   // overall recombination rate
  double G_f; // gene conversion fraction
  double G_l; // average stretch length

  void initialize_rng()
  {
    if (size() == 0)       { cerr << "ERROR (initialize): empty chromosome" << endl; exit(1); }
    if (rec_r.size() == 0) { cerr << "ERROR (initialize): recombination rate not specified" << endl; exit(1); }
    if (!(M>=0))           { cerr << "ERROR (initialize): invalid mutation rate" << endl; exit(1); }

    L = 0;

    for (int i=0; i<size(); i++)
      {
	if (genomic_element_types.count(operator[](i).i)==0)
	  { 
	    cerr << "ERROR (initialize): genomic element type " << operator[](i).i << " not defined" << endl; exit(1); 
	  }
      }

    for (map<int,genomic_element_type>::iterator it = genomic_element_types.begin(); it!=genomic_element_types.end(); it++)
      {
	for (int j=0; j<it->second.m.size(); j++)
	  {
	    if (mutation_types.count(it->second.m[j]) == 0)
	      {
	        cerr << "ERROR (initialize): mutation type " << it->second.m[j] << " not defined" << endl; exit(1); 
	      }
	  }
      }  

    double A[size()]; int l = 0;
    for (int i=0; i<size(); i++) 
      { 
	if (operator[](i).e > L) { L = operator[](i).e; }
	int l_i = operator[](i).e - operator[](i).s + 1.0; 
	A[i] = (double)l_i; l += l_i;
      }
    LT_M = gsl_ran_discrete_preproc(size(),A); M = M*(double)l;

    double B[rec_r.size()];
    B[0] = rec_r[0]*(double)rec_x[0]; R += B[0];
    for (int i=1; i<rec_r.size(); i++) 
      { 
	B[i] = rec_r[i]*(double)(rec_x[i]-rec_x[i-1]); R+= B[i];
	if (rec_x[i]>L) { L = rec_x[i]; }
      }
    LT_R = gsl_ran_discrete_preproc(rec_r.size(),B);
  }


  int draw_n_mut() { return gsl_ran_poisson(rng,M); } // to draw the total number of mutations
 
  mutation draw_new_mut()
  {
    int g = gsl_ran_discrete(rng,LT_M); // genomic element
    genomic_element_type ge_type = genomic_element_types.find(operator[](g).i)->second; // genomic element type

    int mut_type_id = ge_type.draw_mutation_type(); // mutation type id
    mutation_type mut_type = mutation_types.find(mut_type_id)->second; // mutation type

    int   x = operator[](g).s + gsl_rng_uniform_int(rng,operator[](g).e - operator[](g).s + 1); // position    
    float s = mut_type.draw_s(); // selection coefficient

    return mutation(mut_type_id,x,s);
  }


  vector<int> draw_breakpoints()
  {
    vector<int> r;

    // draw recombination breakpoints

    int nr = gsl_ran_poisson(rng,R);

    for (int i=0; i<nr; i++)
      {
	int x = 0;
	int j = gsl_ran_discrete(rng,LT_R);

	if (j==0) { x = gsl_rng_uniform_int(rng,rec_x[j]); }
	else     { x = rec_x[j-1] + gsl_rng_uniform_int(rng,rec_x[j]-rec_x[j-1]); }

	r.push_back(x);

	if (gsl_rng_uniform(rng)<G_f) // recombination results in gene conversion 
	  {
	    int x2 = x+gsl_ran_geometric(rng,1.0/G_l);
	    r.push_back(x2);
	  }
      }

    return r;
  }
}; // end of class chromosome


class polymorphism
{
public:

  int   i; // mutation id
  int   t; // mutation type
  float s; // selection coefficient
  int   n; // prevalence

  polymorphism(int I, int T, float S, int N)
  {
    i = I;
    t = T;
    s = S;
    n = N;
  }

  void print(int x, chromosome& chr) 
  { 
    float h = chr.mutation_types.find(t)->second.h;
    cout << i << " m" << t << " " << x+1 << " " << s << " " << h << " "<< n << endl; 
  }

  void print(ofstream& outfile, int x, chromosome& chr) 
  { 
    float h = chr.mutation_types.find(t)->second.h;
    outfile << i << " m" << t << " " << x+1 << " " << s << " " << h << " "<< n << endl; 
  }

  void print_noi(int x, chromosome& chr) 
  { 
    float h = chr.mutation_types.find(t)->second.h;
    cout << "m" << t << " " << x+1 << " " << s << " " << h << " "<< n << endl; 
  }
}; // end of class polymorphism


class substitution
{
public:

  int   t; // mutation type
  int   x; // position
  float s; // selection coefficient
  int   g; // fixation time

  substitution(mutation M, int G)
  {
    t = M.t;
    x = M.x;
    s = M.s;
    g = G;
  }

  void print(chromosome& chr) 
  { 
    float h = chr.mutation_types.find(t)->second.h;
    cout << " m" << t << " " << x+1 << " " << s << " " << h << " "<< g << endl; 
  }
}; // end of class substitution


class introduced_mutation : public mutation // inheriting from mutation
{
public:

  int i;   // subpopulation into which mutation is to be introduced
  int nAA; // number of homozygotes
  int nAa; // number of heterozygotes
  char l;  // linkage flag ('e' for linkage equilibrium, 'd' for linkage disequilibrium

  // constructor
  introduced_mutation(int T, int X, int I, int NAA, int NAa, char L)
  {
    t = T; // mutation type
    x = X; // genomic position
    i = I; // subpopulation where introduced
    nAA = NAA; // number of homozygous carriers
    nAa = NAa; // number of heterozygous carriers
    l = L; // linkage flag

  }
}; // end of class introduced_mutation


class partial_sweep
{
 public:

  int t;
  int x;
  float p;

  partial_sweep(int T, int X, float P)
  {
    t = T; x = X; p = P;
  }
}; // end of class introduced_mutation


class genome : public vector<mutation>
{
  // a genome is a vector of mutations
}; // end of class genome


// top-level function

genome fixed(genome& G1, genome& G2)
{
  // return genome G consisting only of the mutations that are present in both G1 and G2

  genome G;

  vector<mutation>::iterator g1 = G1.begin();
  vector<mutation>::iterator g2 = G2.begin();

  vector<mutation>::iterator g1_max = G1.end();
  vector<mutation>::iterator g2_max = G2.end();
  
  while (g1 != g1_max && g2 != g2_max)
      {
	// advance g1 while g1.x < g2.x

	while (g1 != g1_max && g2 != g2_max && (*g1).x < (*g2).x) { g1++; }

	// advance g2 while g1.x < g2.x

	while (g1 != g1_max && g2 != g2_max && (*g2).x < (*g1).x) { g2++; }
	   
	// identify shared mutations at positions x and add to G

	if (g2 != g2_max && g1 != g1_max && (*g2).x == (*g1).x)
	  {
	    int x = (*g1).x;

	    vector<mutation>::iterator temp;

	    while (g1 != g1_max && (*g1).x == x)
	      {
		temp = g2;
		while (temp != g2_max && (*temp).x == x)
		  {
		    if ((*temp).t==(*g1).t && (*temp).s==(*g1).s) { G.push_back(*g1); }
		    temp++;
		  }
		g1++;
	      }
	  }
      }

  return G;
} // end of method fixed()


// top-level method

genome polymorphic(genome& G1, genome& G2)
{
  // return genome G consisting only of the mutations in G1 that are not in G2

  genome G;

  vector<mutation>::iterator g1 = G1.begin();
  vector<mutation>::iterator g2 = G2.begin();

  vector<mutation>::iterator g1_max = G1.end();
  vector<mutation>::iterator g2_max = G2.end();
  
  while (g1 != g1_max && g2 != g2_max)
      {
	// advance g1 while g1.x < g2.x

	while (g1 != g1_max && g2 != g2_max && (*g1).x < (*g2).x) { G.push_back(*g1); g1++; }

	// advance g2 while g1.x < g2.x

	while (g2 != g2_max && g1 != g1_max && (*g2).x < (*g1).x) { g2++; }
	   
	// identify polymorphic mutations at positions x and add to G

	if (g2 != g2_max && g1 != g1_max && (*g2).x == (*g1).x)
	  {
	    int x = (*g1).x;

	    // go through g1 and check for those mutations that are not present in g2

	    vector<mutation>::iterator temp = g2;

	    while (g1 != g1_max && (*g1).x == x)
	      {
		bool poly = 1;

		while (temp != g2_max && (*temp).x == x)
		  {
		    if ((*g1).t==(*temp).t && (*g1).s==(*temp).s) { poly = 0; }
		    temp++;
		  }
		if (poly == 1) { G.push_back(*g1); }
		g1++;
	      }

	    while (g2 != g2_max && (*g2).x == x) { g2++; }
	  }
      }

  while (g1 != g1_max) { G.push_back(*g1); g1++; }

  return G;
} // end of method polymorphic()


class subpopulation
{
  // a subpopulation is described by the vector G of 2N genomes
  // individual i is constituted by the two genomes 2*i and 2*i+1

private:
  
  gsl_ran_discrete_t* LT; // a look-up table

public:

  int    N; // population size  
  double S; // selfing fraction

  environment E; // environment

  vector<genome> G_parent; // parent population
  vector<genome> G_child; // offspring population

  map<int,double> m; // m[i]: fraction made up of migrants from subpopulation i per generation
 
  // constructor
  subpopulation(int n)
  {
    N = n;
    S = 0.0;
    G_parent.resize(2*N); G_child.resize(2*N);
    double A[N];
    for (int i = 0; i < N; i++)
      { A[i] = 1.0; }
    LT = gsl_ran_discrete_preproc(N, A); // to draw individuals at random, which uniform weights

    // TODO: Assign the reference environment.

  } // end of constructor


  int draw_individual()
  {
    return gsl_ran_discrete(rng, LT);
  }


  void update_fitness(chromosome& chr)
  {
    // calculate fitnesses in parent population and create new lookup table
    
    gsl_ran_discrete_free(LT);
    double A[(int)(G_parent.size()/2)]; // individuals are diploid
    for (int i = 0; i < (int)(G_parent.size()/2); i++)
      {
        A[i] = W(2*i, 2*i+1, chr);
      }
    LT = gsl_ran_discrete_preproc((int)(G_parent.size()/2), A); // assigning weights to individuals
      // that correspond to their (abolute) fitness
  } // end of method update_fitness()


  double W(int i, int j, chromosome& chr)
  {
    // calculate the fitness of the individual constituted by genomes i and j in the parent
      // population, where i and j are indices to the vector population
    
    double w = 1.0;

    // a genome is a vector of mutations
    vector<mutation>::iterator pi = G_parent[i].begin(); // iterator pointing to the first mutation
      // in parent i (i.e. genome i of the parental population)
    vector<mutation>::iterator pj = G_parent[j].begin();

    vector<mutation>::iterator pi_max = G_parent[i].end(); // iterator pointing to the last
      // mutation in parent i (i.e. genome i of the parental population)
    vector<mutation>::iterator pj_max = G_parent[j].end();

    while (w > 0 && (pi != pi_max || pj != pj_max)) // while fitness is strictly positive, and
      // at least one parent has mutations not jet visited
      {

        // advance i while pi.x < pj.x (x is the physical position)
        while (pi != pi_max && (pj == pj_max || (*pi).x < (*pj).x)) // while there are unvisited
          // mutations left in genome i and (there are no unvisited mutations left in genome j or
          // the position of the current mutation in genome i is smaller than the current mutation
          // in genome j)
          // TODO: GO ON HERE. Differentiate between 'a'(dditive) and 'm'(multiplicative) fitness
            // interaction scheme. Introduce environment specific fitnesses.
          {
            // changed by SA: additive instead of mulitplicative fitness interaction
            // if ((*pi).s != 0) // if this mutation in genome i is not neutral
            //  {
            //    w = w*(1.0+chr.mutation_types.find((*pi).t)->second.h*(*pi).s);
            //  } // end of if this mutation in i is not neutral
            if ((*pi).s != 0) // if this mutation in individual i is not neutral
              {
                w = w+chr.mutation_types.find((*pi).t)->second.h*(*pi).s;
              }
          pi++;
	  }
	   
	// advance j while pj.x < pi.x

	while (pj != pj_max && (pi == pi_max || (*pj).x < (*pi).x))
	  {
          // changed by SA: additive instead of mulitplicative fitness interaction
          // if ((*pj).s != 0) { w = w*(1.0+chr.mutation_types.find((*pj).t)->second.h*(*pj).s); }
          if ((*pj).s != 0) { w = w+chr.mutation_types.find((*pj).t)->second.h*(*pj).s; }
          pj++;
	  }
	
	// check for homozygotes and heterozygotes at x

	if (pi != pi_max && pj != pj_max && (*pj).x == (*pi).x)
	  {
	    int x = (*pi).x; 
	   
	    vector<mutation>::iterator pi_start = pi;

	    // advance through pi

	    while (pi != pi_max && (*pi).x == x)
	      {
		if ((*pi).s != 0.0)
		  {
		    vector<mutation>::iterator temp_j = pj; 
		    bool homo = 0;

		    while (homo == 0 && temp_j != pj_max && (*temp_j).x == x)
		      {
			if ((*pi).t == (*temp_j).t && (*pi).s == (*temp_j).s) 
			  { 
                  // changed by SA: additive instead of mulitplicative fitness interaction
                  // w = w*(1.0+(*pi).s); homo = 1;
                  w = w+(*pi).s; homo = 1;
			  }
			temp_j++;
		      }
              // changed by SA: additive instead of mulitplicative fitness interaction
              // if (homo == 0) { w = w*(1.0+chr.mutation_types.find((*pi).t)->second.h*(*pi).s); }
              if (homo == 0) { w = w+chr.mutation_types.find((*pi).t)->second.h*(*pi).s; }
		  }
		pi++;
	      }

	    // advance through pj

	    while (pj != pj_max && (*pj).x == x)
	      {
		if ((*pj).s != 0.0)
		  {
		    vector<mutation>::iterator temp_i = pi_start; 
		    bool homo = 0;

		    while (homo == 0 && temp_i != pi_max && (*temp_i).x == x)
		      {
			if ((*pj).t == (*temp_i).t && (*pj).s == (*temp_i).s) { homo = 1; }
			temp_i++;
		      }
              // changed by SA: additive instead of mulitplicative fitness interaction
              // if (homo == 0) { w = w*(1.0+chr.mutation_types.find((*pj).t)->second.h*(*pj).s); }
              if (homo == 0) { w = w+chr.mutation_types.find((*pj).t)->second.h*(*pj).s; }
		  }
		pj++;
	      }
	  }
      }

    if (w<0) { w = 0.0; }

    return w;
  } // end of method W()

  void swap()
  {
    G_child.swap(G_parent);
    G_child.resize(2*N);
  }
}; // end of class subpopulation


class population : public map<int, subpopulation>
{
  // the population is a map of subpopulations

public: 

  vector<substitution> Substitutions;

  map<int,subpopulation>::iterator it;

  vector<string> parameters;
  // TODO: GO ON HERE. Change to include fitness interaction as argument. Change class
    // subpopulation accordingly.
  void add_subpopulation(int i, unsigned int N) 
  {
    // add new empty subpopulation i of size N (i is the key of the subpopulation)

    if (count(i) != 0) // counts number of occurrences of subpopulation with key equal to i
      {
        cerr << "ERROR (add subpopulation): subpopulation p" << i << " already exists" << endl;
        exit(1);
      }
    if (N < 1)
      {
        cerr << "ERROR (add subpopulation): subpopulation p" << i << " empty" << endl;
        exit(1);
      }

    insert(pair<int, subpopulation>(i, subpopulation(N)));
  } // end of add_subpopulation() method


  void add_subpopulation(int i, int j, unsigned int N) 
  { 
    // add new subpopulation i of size N individuals drawn from source subpopulation j

    if (count(i)!=0) { cerr << "ERROR (add subpopulation): subpopulation p"<< i << " already exists" << endl; exit(1); }
    if (count(j)==0) { cerr << "ERROR (add subpopulation): source subpopulation p"<< j << " does not exists" << endl; exit(1); }
    if (N<1)         { cerr << "ERROR (add subpopulation): subpopulation p"<< i << " empty" << endl; exit(1); }

    insert(pair<int,subpopulation>(i,subpopulation(N))); 

    for (int p=0; p<find(i)->second.N; p++)
      {
	// draw individual from subpopulation j and assign to be a parent in i  

	int m = find(j)->second.draw_individual();
	
	find(i)->second.G_parent[2*p] = find(j)->second.G_parent[2*m];
	find(i)->second.G_parent[2*p+1] = find(j)->second.G_parent[2*m+1];
      }
  }


  void set_size(int i, unsigned int N) 
  {
    // set size of subpopulation i to N

    if (count(i)==0) { cerr << "ERROR (change size): no subpopulation p"<< i << endl; exit(1); }

    if (N==0) // remove subpopulation i 
      {
	erase(i); 
	for (it = begin(); it != end(); it++) { it->second.m.erase(i); } 
      }
    else { find(i)->second.N = N; find(i)->second.G_child.resize(2*N); }
  }

  void set_selfing(int i, double s) 
  { 
    // set fraction s of i that reproduces by selfing
 
    if (count(i)==0)    { cerr << "ERROR (set selfing): no subpopulation p"<< i << endl; exit(1); }
    if (s<0.0 || s>1.0) { cerr << "ERROR (set selfing): selfing fraction has to be within [0,1]" << endl; exit(1); }

    find(i)->second.S = s; 
  }


  void set_migration(int i, int j, double m) 
  { 
    // set fraction m of i that originates as migrants from j per generation  

    if (count(i)==0)    { cerr << "ERROR (set migration): no subpopulation p"<< i << endl; exit(1); }
    if (count(j)==0)    { cerr << "ERROR (set migration): no subpopulation p"<< j << endl; exit(1); }
    if (m<0.0 || m>1.0) { cerr << "ERROR (set migration): migration fraction has to be within [0,1]" << endl; exit(1); }

    if (find(i)->second.m.count(j) !=0) { find(i)->second.m.erase(j); }

    find(i)->second.m.insert(pair<int,double>(j,m)); 
  }


  void execute_event(event& E, int g, chromosome& chr, vector<int>& FM)
  {
    char type = E.t;

    if (type == 'P') // add subpopulation
      { 
	if (E.np == 2) // empty subpopulation
	  { 
	    string sub = E.s[0]; sub.erase(0, 1);

	    int i = atoi(sub.c_str());
	    int n = (int)atof(E.s[1].c_str());
	    add_subpopulation(i,n);
	  }
	      
	if (E.np == 3) // drawn from source population
	  {
	    string sub1 = E.s[0]; sub1.erase(0,1);
	    string sub2 = E.s[2]; sub2.erase(0,1);

	    int i = atoi(sub1.c_str());
	    int j = atoi(sub2.c_str());
	    int n = (int)atof(E.s[1].c_str());
	    add_subpopulation(i,j,n);
	  } 
      }
	  
    if (type == 'N') // set subpopulation size
      { 
	string sub = E.s[0]; sub.erase(0, 1);

	int i = atoi(sub.c_str());
	int n = (int)atof(E.s[1].c_str());

	set_size(i,n);
      }

    if (type == 'S') // set selfing rate
      { 
	string sub = E.s[0]; sub.erase(0, 1);

	int i    = atoi(sub.c_str());
	double s = atof(E.s[1].c_str());

	set_selfing(i,s);
      }
	  
    if (type == 'M') // change migration rate
      {
	string sub1 = E.s[0]; sub1.erase(0,1);
	string sub2 = E.s[1]; sub2.erase(0,1);

	int    i = atoi(sub1.c_str());
	int    j = atoi(sub2.c_str());
	double m = atof(E.s[2].c_str());

	set_migration(i,j,m); 
      }

    if (type == 'A') // output state of entire population
      {
	if (E.s.size() == 0)
	  {
	    cout << "#OUT: " << g << " A" << endl;
	    print_all(chr); 
	  }	
	if (E.s.size() == 1)
	  {
	    ofstream outfile;
	    outfile.open (E.s[0].c_str());

	    for (int i=0; i<parameters.size(); i++) { outfile << parameters[i] << endl; }

	    if (outfile.is_open()) 
	      { 
		outfile << "#OUT: " << g << " A " << E.s[0].c_str() << endl;
		print_all(outfile,chr);
		outfile.close(); 
	      }
	    else { cerr << "ERROR (output): could not open "<< E.s[0].c_str() << endl; exit(1); }
	  }
      }

    if (type == 'R') // output random subpopulation sample
      {
	string sub = E.s[0]; sub.erase(0, 1);

	int    i = atoi(sub.c_str());
	int    n = atoi(E.s[1].c_str());   
	cout << "#OUT: " << g << " R p" << i << " " << n << endl;

	if (E.s.size() == 3 && E.s[2] == "MS") { print_sample_ms(i,n,chr); }
	else { print_sample(i,n,chr); }
      }

    if (type == 'F') // output list of fixed mutations
      {
	cout << "#OUT: " << g << " F " << endl;
	cout << "Mutations:" << endl;
	for (int i=0; i<Substitutions.size(); i++) { cout << i+1; Substitutions[i].print(chr); }
      }

    if (type == 'T') // track mutation-types
      {
	string sub = E.s[0]; sub.erase(0, 1);
	FM.push_back(atoi(sub.c_str()));
      }
  }


  void introduce_mutation(introduced_mutation M, chromosome& chr) 
  {
    // introduce user-defined mutation
    // TODO: Make sure that if a mutation is to be introduced at a position at which there already
        // is a mutation of a non-neutral type present, the pre-existing mutation is re-assigned
        // a neutral type (there must be only one mutation of a non-neutral type at any site
        // at any time)
    if (count(M.i)==0) { cerr << "ERROR (predetermined mutation): subpopulation "<< M.i << " does not exists" << endl; exit(1); }
    if (chr.mutation_types.count(M.t) == 0) 
      { 
	cerr << "ERROR (predetermined mutation): mutation type m"<< M.t << " has not been defined" << endl; exit(1); 
      }
    if (find(M.i)->second.G_child.size()/2 < M.nAA + M.nAa) 
      { 
	cerr << "ERROR (predetermined mutation): not enough individuals in subpopulation "<< M.i << endl; exit(1); 
      }

    mutation m;

    m.x = M.x; // position on genome
    m.t = M.t; // mutation type (to be precise, its key)
    // draw selection coefficient for novel mutation
    m.s = chr.mutation_types.find(M.t)->second.draw_s();

    subpopulation *sp = &find(M.i)->second;
    
    // test: print all children
    // print_all(chr);
      
    // shuffle genomes in subpopulation (randomisation of order)
      // this is to generate linkage equilibrium among loci and
      // has been added by Simon Aeschbacher on 05/19/2014.
    std::random_shuffle((*sp).G_child.begin(), (*sp).G_child.end());
    
    // test: print all children
    // print_all(chr);
    
    // introduce homozygotes

    for (int j=0; j<M.nAA; j++)
      {
    // recall: a genonme is a vector of mutations
    // find subpopulation, then return two of its children
	genome *g1 = &find(M.i)->second.G_child[2*j];
	genome *g2 = &find(M.i)->second.G_child[2*j+1];
    // append the new mutation to the end of the genotypes.
	(*g1).push_back(m);
	(*g2).push_back(m);
    // sort the genomes
	sort((*g1).begin(),(*g1).end());
	sort((*g2).begin(),(*g2).end());
    // unique(ForwardIterator first, ForwardIterator last): removes al but
          // the first element from every consecutive group of equivalent
          // elements in the range [first, last)
    // erase(iterator first, iterator last): removes the range
          // [first, last) of elements
    // remove *consecutive* duplicate mutations
	(*g1).erase(unique((*g1).begin(),(*g1).end()),(*g1).end());
	(*g2).erase(unique((*g2).begin(),(*g2).end()),(*g2).end());
      }

    // introduce heterozygotes

    for (int j=M.nAA; j<M.nAA+M.nAa; j++) 
      { 
	genome *g1 = &find(M.i)->second.G_child[2*j];
	(*g1).push_back(m);
	sort((*g1).begin(),(*g1).end());
    // remove duplicate mutations
	(*g1).erase(unique((*g1).begin(),(*g1).end()),(*g1).end());
      }
  }


  void track_mutations(int g, vector<int>& TM, vector<partial_sweep>& PS, chromosome& chr)
  {
    // output trajectories of followed mutations and set s=0 for partial sweeps 

    // find all polymorphism of the types that are to be tracked

    for (it = begin(); it != end(); it++) // go through all subpopulations
      {
	multimap<int,polymorphism> P;
	multimap<int,polymorphism>::iterator P_it;

	for (int i=0; i<2*it->second.N; i++) // go through all children
	  {
	    for (int k=0; k<it->second.G_child[i].size(); k++) // go through all mutations
	      {
		for (int j=0; j<TM.size(); j++)
		  {
		    if (it->second.G_child[i][k].t == TM[j]) { add_mut(P,it->second.G_child[i][k]); }
		  }
	      }
	  }

	// out put the frequencies of these mutations in each subpopulation

	for (P_it = P.begin(); P_it != P.end(); P_it++) 
	  { 
	    cout << "#OUT: " << g << " T p" << it->first << " "; P_it->second.print_noi(P_it->first,chr); 
	  }
      }

    // check partial sweeps

    multimap<int,polymorphism> P;
    multimap<int,polymorphism>::iterator P_it;

    if (PS.size()>0)
      {
	P.clear();

	int N = 0; for (it = begin(); it != end(); it++) { N += it->second.N; }

	// find all polymorphism that are supposed to undergo partial sweeps

	for (it = begin(); it != end(); it++) // go through all subpopulations
	  {
	    for (int i=0; i<2*it->second.N; i++) // go through all children
	      {
		for (int k=0; k<it->second.G_child[i].size(); k++) // go through all mutations
		  {
		    for (int j=0; j<PS.size(); j++)
		      {
			if (it->second.G_child[i][k].x == PS[j].x && it->second.G_child[i][k].t == PS[j].t) 
			  {
			    add_mut(P,it->second.G_child[i][k]); 
			  }
		      }
		  }
	      }
	  }
    
	// check whether a partial sweep has reached its target frequency

	for (P_it = P.begin(); P_it != P.end(); P_it++) 
	  { 
	    for (int j=0; j<PS.size(); j++)
	      {
		if (P_it->first == PS[j].x && P_it->second.t == PS[j].t)
		  {
		    if (((float)P_it->second.n)/(2*N) >= PS[j].p)
		      {
			// sweep has reached target frequency, set all s to zero
			
			for (it = begin(); it != end(); it++) // go through all subpopulations
			  {
			    for (int i=0; i<2*it->second.N; i++) // go through all children
			      {
				for (int k=0; k<it->second.G_child[i].size(); k++) // go through all mutations
				  {
				    if (it->second.G_child[i][k].x == PS[j].x && it->second.G_child[i][k].t == PS[j].t)
				      {
					it->second.G_child[i][k].s = 0.0;
				      }
				  }
			      }
			  }
			PS.erase(PS.begin()+j); j--;
		      }
		  }
	      }
	  }
      }
  }


  void evolve_subpopulation(int i, chromosome& chr)
  {
    int g1, g2, p1, p2, n_mut_1, n_mut_2;


    // create map of shuffled children ids

    int child_map[find(i)->second.N];          
    for (int j = 0; j < find(i)->second.N; j++) { child_map[j] = j; }
    gsl_ran_shuffle(rng,child_map,find(i)->second.N,sizeof(int));


    int c = 0; // counter over all N children (will get mapped to child_map[c])

    // migration, loop over all source populations
    
    map<int,double>::iterator it;

    for (map<int,double>::iterator it = find(i)->second.m.begin(); it != find(i)->second.m.end(); it++)
      {
	int n_migrants = (int)(it->second * find(i)->second.N + 0.5);
        
	for (int m=0; m<n_migrants; m++) 
	  {
	    if (c >= find(i)->second.N) { cerr << "ERROR (evolve subpopulation): too many migrants in subpopulation "<< i << endl; exit(1); }

	    g1 = 2*child_map[c];   // child genome 1
	    g2 = 2*child_map[c]+1; // child genome 2

	    // draw parents in source population

	    p1 = gsl_rng_uniform_int(rng,find(it->first)->second.G_parent.size()/2);
	    if (gsl_rng_uniform(rng) < find(it->first)->second.S) { p2 = p1; }
	    else { p2 = gsl_rng_uniform_int(rng,find(it->first)->second.G_parent.size()/2); }

	    // recombination, gene-conversion, mutation

	    crossover_mutation(i,g1,it->first,2*p1,2*p1+1,chr);
	    crossover_mutation(i,g2,it->first,2*p2,2*p2+1,chr);

	    c++;
	  }
      }
	    
    // remainder

    while (c<find(i)->second.N) 
      {
	g1 = 2*child_map[c];   // child genome 1
	g2 = 2*child_map[c]+1; // child genome 2

	p1 = find(i)->second.draw_individual();                 // parent 1
	if (gsl_rng_uniform(rng) < find(i)->second.S) { p2 = p1; } // parent 2
	else { p2 = find(i)->second.draw_individual(); }

	crossover_mutation(i,g1,i,2*p1,2*p1+1,chr);
	crossover_mutation(i,g2,i,2*p2,2*p2+1,chr);

	c++;
      }
  }


  void crossover_mutation(int i, int c, int j, int P1, int P2, chromosome& chr)
  {
    // child genome c in subpopulation i is assigned outcome of cross-overs at breakpoints r 
    // between parent genomes p1 and p2 from subpopulation j and new mutations added
    // 
    // example R = (r1,r2)
    // 
    // mutations (      x < r1) assigned from p1
    // mutations (r1 <= x < r2) assigned from p2
    // mutations (r2 <= x     ) assigned from p1
    //
    // p1 and p2 are swapped in half of the cases to assure random assortement

    if (gsl_rng_uniform_int(rng,2)==0) { int swap = P1; P1 = P2; P2 = swap; } // swap p1 and p2

    find(i)->second.G_child[c].clear();

    // create vector with the mutations to be added

    vector<mutation> M;
    int n_mut = chr.draw_n_mut();
    for (int k=0; k<n_mut; k++) { M.push_back(chr.draw_new_mut()); }
    sort(M.begin(),M.end());
    
    // create vector with recombination breakpoints

    vector<int> R = chr.draw_breakpoints(); 
    R.push_back(chr.L+1);
    sort(R.begin(),R.end());
    R.erase(unique(R.begin(),R.end()),R.end());
    
    vector<mutation>::iterator p1 = find(j)->second.G_parent[P1].begin();
    vector<mutation>::iterator p2 = find(j)->second.G_parent[P2].begin();

    vector<mutation>::iterator p1_max = find(j)->second.G_parent[P1].end();
    vector<mutation>::iterator p2_max = find(j)->second.G_parent[P2].end();

    vector<mutation>::iterator m     = M.begin();
    vector<mutation>::iterator m_max = M.end();

    vector<mutation>::iterator p     = p1;
    vector<mutation>::iterator p_max = p1_max;

    int r = 0; int r_max  = R.size(); int n = 0; bool present;

    while (r != r_max)
      {
	while ((p != p_max && (*p).x < R[r]) || (m != m_max && (*m).x < R[r]))
	  {
	    while (p != p_max && (*p).x < R[r] && (m == m_max || (*p).x <= (*m).x))
	      {
		present = 0;
		if (n != 0 && find(i)->second.G_child[c].back().x == (*p).x)
		  {
		    int k = n-1;
		    while (present == 0 && k >= 0)
		      {
			if (find(i)->second.G_child[c][k] == (*p)) { present = 1; }
			k--;
		      }
		  }
		if (present == 0) { find(i)->second.G_child[c].push_back(*p); n++; }
		p++;
	      }
	    while (m != m_max && (*m).x < R[r] && (p == p_max || (*m).x <= (*p).x))
	      {
	    	present = 0;
		if (n != 0 && find(i)->second.G_child[c].back().x == (*m).x)
		  {
		    int k = n-1;
		    while (present == 0 && k >= 0)
		      {
			if (find(i)->second.G_child[c][k] == (*m)) { present = 1; }
			k--;
		      }
		  }
		if (present == 0) { find(i)->second.G_child[c].push_back(*m); n++; }
		m++;
	      }
	  }

	// swap parents

	p1 = p2; p1_max = p2_max; p2 = p; p2_max = p_max; p = p1; p_max = p1_max; 

	while (p != p_max && (*p).x < R[r]) { p++; }

	r++;
      }
  }


  void swap_generations(int g, chromosome& chr)
  {
    // find and remove fixed mutations from the children in all subpopulations
    
    remove_fixed(g);

    // make children the new parents and update fitnesses

    for (it = begin(); it != end(); it++) 
      { 
	it->second.swap();
	it->second.update_fitness(chr); 
      }
  }


  void remove_fixed(int g)
  {
    // find mutations that are fixed in all child subpopulations and return vector with their ids

    genome G = begin()->second.G_child[0];

    for (it = begin(); it != end(); it++) // subpopulations
      {
	for (int i=0; i<2*it->second.N; i++) // child genomes
	  {
	    G = fixed(it->second.G_child[i],G);
	  }
      }

    if (G.size()>0)
      {
	for (it = begin(); it != end(); it++) // subpopulations
	  {
	    for (int i=0; i<2*it->second.N; i++) // child genomes
	      {
		it->second.G_child[i] = polymorphic(it->second.G_child[i],G);
	      }
	  }
	for (int i=0; i<G.size(); i++) { Substitutions.push_back(substitution(G[i],g)); } 
      }
  }


  void print_all(chromosome& chr)
  {
    // print all mutations and all genomes 

    cout << "Populations:" << endl;
    for (it = begin(); it != end(); it++) {  cout << "p" << it->first << " " << it->second.N << endl; }

    multimap<int,polymorphism> P;
    multimap<int,polymorphism>::iterator P_it;

    // add all polymorphisms

    for (it = begin(); it != end(); it++) // go through all subpopulations
      {
	for (int i=0; i<2*it->second.N; i++) // go through all children
	  {
	    for (int k=0; k<it->second.G_child[i].size(); k++) // go through all mutations
	      {
		add_mut(P,it->second.G_child[i][k]);
	      }
	  }
      }

    cout << "Mutations:"  << endl;
    
    for (P_it = P.begin(); P_it != P.end(); P_it++) { P_it->second.print(P_it->first,chr); }

    cout << "Genomes:" << endl;

    // print all genomes

    for (it = begin(); it != end(); it++) // go through all subpopulations
      {
	for (int i=0; i<2*it->second.N; i++) // go through all children
	  {
	    cout << "p" << it->first << ":" << i+1;

	    for (int k=0; k<it->second.G_child[i].size(); k++) // go through all mutations
	      {
		int id = find_mut(P,it->second.G_child[i][k]);
		cout << " " << id; 
	      }
	    cout << endl;
	  }
      }

  }


  void print_all(ofstream& outfile, chromosome& chr)
  {
    // print all mutations and all genomes 

    outfile << "Populations:" << endl;
    for (it = begin(); it != end(); it++) {  outfile << "p" << it->first << " " << it->second.N << endl; }

    multimap<int,polymorphism> P;
    multimap<int,polymorphism>::iterator P_it;

    // add all polymorphisms

    for (it = begin(); it != end(); it++) // go through all subpopulations
      {
	for (int i=0; i<2*it->second.N; i++) // go through all children
	  {
	    for (int k=0; k<it->second.G_child[i].size(); k++) // go through all mutations
	      {
		add_mut(P,it->second.G_child[i][k]);
	      }
	  }
      }

    outfile << "Mutations:"  << endl;
    
    for (P_it = P.begin(); P_it != P.end(); P_it++) { P_it->second.print(outfile,P_it->first,chr); }

    outfile << "Genomes:" << endl;

    // print all genomes

    for (it = begin(); it != end(); it++) // go through all subpopulations
      {
	for (int i=0; i<2*it->second.N; i++) // go through all children
	  {
	    outfile << "p" << it->first << ":" << i+1;

	    for (int k=0; k<it->second.G_child[i].size(); k++) // go through all mutations
	      {
		int id = find_mut(P,it->second.G_child[i][k]);
		outfile << " " << id; 
	      }
	    outfile << endl;
	  }
      }

  }


  void print_sample(int i, int n, chromosome& chr)
  {
    // print sample of n genomes from subpopulation  i

    if (count(i)==0) { cerr << "ERROR (output): subpopulation p"<< i << " does not exists" << endl; exit(1); }

    vector<int> sample; 

    multimap<int,polymorphism> P;
    multimap<int,polymorphism>::iterator P_it;
    
    for (int s=0; s<n; s++) 
      { 
	int j = gsl_rng_uniform_int(rng,find(i)->second.G_child.size());
	sample.push_back(j);

	for (int k=0; k<find(i)->second.G_child[j].size(); k++) // go through all mutations
	  {
	    add_mut(P,find(i)->second.G_child[j][k]);
	  }
      }

    cout << "Mutations:"  << endl;
    
    for (P_it = P.begin(); P_it != P.end(); P_it++) { P_it->second.print(P_it->first,chr); }

    cout << "Genomes:" << endl;

    // print all genomes

    for (int j=0; j<sample.size(); j++) // go through all children
      {
	cout << "p" << find(i)->first << ":" << sample[j]+1;

	for (int k=0; k<find(i)->second.G_child[sample[j]].size(); k++) // go through all mutations
	  {
	    int id = find_mut(P,find(i)->second.G_child[sample[j]][k]);
	    cout << " " << id; 
	  }
	cout << endl;
      }
  }


  void print_sample_ms(int i, int n, chromosome& chr)
  {
    // print sample of n genomes from subpopulation  i

    if (count(i)==0) { cerr << "ERROR (output): subpopulation p"<< i << " does not exists" << endl; exit(1); }

    vector<int> sample; 

    multimap<int,polymorphism> P;
    multimap<int,polymorphism>::iterator P_it;
    
    for (int s=0; s<n; s++) 
      { 
	int j = gsl_rng_uniform_int(rng,find(i)->second.G_child.size());
	sample.push_back(j);

	for (int k=0; k<find(i)->second.G_child[j].size(); k++) // go through all mutations
	  {
	    add_mut(P,find(i)->second.G_child[j][k]);
	  }
      }

    // print header

    cout << endl << "//" << endl << "segsites: " << P.size() << endl;

    // print all positions

    if (P.size()>0)
      {
	cout << "positions:"; 
	for (P_it = P.begin(); P_it != P.end(); P_it++) 
	  { 
	    cout << " " << fixed << setprecision(7) << (double)(P_it->first+1)/(chr.L+1); 
	  }
	cout << endl;
      }

    // print genotypes

    for (int j=0; j<sample.size(); j++) // go through all children
      {
	string genotype(P.size(),'0');

	for (int k=0; k<find(i)->second.G_child[sample[j]].size(); k++) // go through all mutations
	  {
	    int pos = 0;
	    mutation m = find(i)->second.G_child[sample[j]][k];

	    for (P_it = P.begin(); P_it != P.end(); P_it++) 
	      {
		if (P_it->first == m.x && P_it->second.t == m.t && P_it->second.s == m.s)
		  {
		    genotype.replace(pos,1,"1");
		    break;
		  }
		pos++;
	      }
	  }
	cout << genotype << endl;
      }
  }


  int find_mut(multimap<int,polymorphism>& P, mutation m)
  {
    // find m in P and return its id

    int id = 0;

    // iterate through all mutations with same position

    multimap<int,polymorphism>::iterator it;
    pair<multimap<int,polymorphism>::iterator,multimap<int,polymorphism>::iterator> range = P.equal_range(m.x);
    it = range.first;

    while (it != range.second)
      {
	if (it->second.t == m.t && it->second.s == m.s) 
	  { 
	    id = it->second.i;
	    it->second.n++;
	    it = range.second;
	  }
	else{ it++; }
      }
    
    return id;
  }


  void add_mut(multimap<int,polymorphism>& P, mutation m)
  {
    // if mutation is present in P increase prevalence, otherwise add it

    int id = 0;

    // iterate through all mutations with same position

    multimap<int,polymorphism>::iterator it;
    pair<multimap<int,polymorphism>::iterator,multimap<int,polymorphism>::iterator> range = P.equal_range(m.x);
    it = range.first;

    while (it != range.second)
      {
	if (it->second.t == m.t && it->second.s == m.s) 
	  { 
	    id = it->second.i;
	    it->second.n++;
	    it = range.second;
	  }
	else{ it++; }
      }

    // if not already present, add mutation to P

    if (id == 0)
      {
	id = P.size()+1;
	P.insert(pair<int,polymorphism>(m.x,polymorphism(id,m.t,m.s,1)));
      }
  } // end of add_mut() method

}; // end of class 'population'


void get_line(ifstream& infile, string& line)
{
  getline(infile, line);
  if (line.find("/")!= string::npos) { line.erase(line.find("/")); } // remove all after "/"; these
        // lines are interpreted as comments
  line.erase(0, line.find_first_not_of(' ')); // remove leading whitespaces
  line.erase(line.find_last_not_of(' ') + 1); // remove trailing whitespaces
} // end of method get_line()


void input_error(int type, string line)
{
  cerr << endl;

  if (type == -2) // no population defined
     {
       cerr << "ERROR (parameter file): no population to simulate:" << endl << endl;
     }
  
  else if (type == -1) // unknown parameter
     {
       cerr << "ERROR (parameter file): unknown parameter: " << line << endl << endl;
     }

  else if (type == 0) // invalid parameter file
    {
      cerr << "ERROR (parameter file): could not open: " << line << endl << endl;
    }
  
  else if (type == 1) // mutation rate
    {
      cerr << "ERROR (parameter file): invalid mutation rate: " << line << endl << endl;
      cerr << "Required syntax:" << endl << endl;
      cerr << "#MUTATION RATE" << endl;
      cerr << "<u>" << endl << endl;
      cerr << "Example:" << endl << endl;
      cerr << "#MUTATION RATE" << endl;
      cerr << "1.5e-8" << endl << endl;
    }

  else if (type == 2) // mutation type
    {
      cerr << "ERROR (parameter file): invalid mutation type: " << line << endl << endl;
      cerr << "Required syntax:" << endl << endl;
      cerr << "#MUTATION TYPES" << endl;
      cerr << "<mutation-type-id> <h> <DFE-type> [DFE parameters]" << endl;
      cerr << "..." << endl << endl;
      cerr << "Example:" << endl << endl;
      cerr << "#MUTATION TYPES" << endl;
      cerr << "m1 0.2 g -0.05 0.2" << endl;
      cerr << "m2 0.0 f 0.0" << endl;
      cerr << "m3 0.5 e 0.01" << endl << endl;
    }

  else if (type == 3) // genomic element type
    {
      cerr << "ERROR (parameter file): invalid genomic element type: " << line << endl << endl;
      cerr << "Required syntax:" << endl << endl;
      cerr << "#GENOMIC ELEMENT TYPES" << endl;
      cerr << "<element-type-id> <mut-type> <x> [<mut-type> <x>...]" << endl;
      cerr << "..." << endl << endl;
      cerr << "Example:" << endl << endl;
      cerr << "#GENOMIC ELEMENT TYPES" << endl;
      cerr << "g1 m3 0.8 m2 0.01 m1 0.19" << endl << endl;
    }

  else if (type == 4) // chromosome organization
    {
      cerr << "ERROR (parameter file): invalid chromosome organization: " << line << endl << endl;
      cerr << "Required syntax:" << endl << endl;
      cerr << "#CHROMOSOME ORGANIZATION" << endl;
      cerr << "<element-type> <start> <end>" << endl;
      cerr << "..." << endl << endl;
      cerr << "Example:" << endl << endl;
      cerr << "#CHROMOSOME ORGANIZATION" << endl;
      cerr << "g1 1000 1999" << endl << endl;
    }

  else if (type == 5) // recombination rate
    {
      cerr << "ERROR (parameter file): invalid recombination rate: " << line << endl << endl;
      cerr << "Required syntax:" << endl << endl;
      cerr << "#RECOMBINATION RATE" << endl;
      cerr << "<interval-end> <r>" << endl;
      cerr << "..." << endl << endl;
      cerr << "Example:" << endl << endl;
      cerr << "#RECOMBINATION RATE" << endl;
      cerr << "10000 1e-8" << endl;
      cerr << "20000 4.5e-8" << endl << endl;
    }
  
  else if (type == 6) // generations
    {
      cerr << "ERROR (parameter file): invalid generations: " << line << endl << endl;
      cerr << "Required syntax:" << endl << endl;
      cerr << "#GENERATIONS" << endl;
      cerr << "<t>" << endl << endl;
      cerr << "Example:" << endl << endl;
      cerr << "#GENERATIONS" << endl;
      cerr << "10000" << endl << endl;
    }

  else if (type == 7) // demography and structure
    {
      cerr << "ERROR (parameter file): invalid demography and structure: " << line << endl << endl;
      cerr << "Required syntax:" << endl << endl;
      cerr << "#DEMOGRAPHY AND STRUCTURE" << endl;
      cerr << "<time> <event-type> [event parameters]" << endl;
      cerr << "..." << endl << endl;
      cerr << "Example:" << endl << endl;
      cerr << "DEMOGRAPHY AND STRUCTURE" << endl;
      cerr << "1 P p1 1000" << endl;
      cerr << "1 S p1 0.05" << endl;
      cerr << "1000 P p2 100 p1" << endl;
      cerr << "1000 S p2 0.05" << endl;
      cerr << "2000 N p1 1e4" << endl;
      cerr << "2000 M p2 p1 0.01" << endl;
      cerr << "4000 E p1 e1" << endl;
      cerr << "4000 E p2 e2" << endl << endl;
    }

  else if (type == 8) // output
    {
      cerr << "ERROR (parameter file): invalid output: " << line << endl << endl;
      cerr << "Required syntax:" << endl << endl;
      cerr << "#OUTPUT" << endl;
      cerr << "<time> <output-type> [output parameters]" << endl;
      cerr << "..." << endl << endl;
      cerr << "Example:" << endl << endl;
      cerr << "OUTPUT" << endl;
      cerr << "2000 A outfile" << endl;
      cerr << "1000 R p1 10" << endl;
      cerr << "1000 R p1 10 MS" << endl;
      cerr << "2000 F" << endl;
      cerr << "1 T m3" << endl << endl;
    }

  else if (type == 9) // initialization
    {
      cerr << "ERROR (parameter file): invalid initialization: " << line << endl << endl;
      cerr << "Required syntax:" << endl << endl;
      cerr << "#INITIALIZATION" << endl;
      cerr << "<filename>" << endl << endl;
      cerr << "Example:" << endl << endl;
      cerr << "#INITIALIZATION" << endl;
      cerr << "outfile" << endl << endl;
    }

  else if (type == 10) // seed
    {
      cerr << "ERROR (parameter file): invalid seed: " << line << endl << endl;
      cerr << "Required syntax:" << endl << endl;
      cerr << "#SEED" << endl;
      cerr << "<seed>" << endl << endl;
      cerr << "Example:" << endl << endl;
      cerr << "#SEED" << endl;
      cerr << "141235" << endl << endl;
    }

  else if (type == 11) // predetermined mutation
    {
      cerr << "ERROR (parameter file): invalid predetermined mutations: " << line << endl << endl;
      cerr << "Required syntax:" << endl << endl;
      cerr << "#PREDETERMINED MUTATIONS" << endl;
      cerr << "<time> <mut-type> <x> <pop-id> <nAA> <nAa> <linkage> [P <f>]" << endl << endl;
      cerr << "Example:" << endl << endl;
      cerr << "#PREDETERMINED MUTATIONS" << endl;
      cerr << "5000 m7 45000 p1 0 1 e" << endl << endl;
    }

  else if (type == 12) // gene conversion
    {
      cerr << "ERROR (parameter file): invalid gene conversion: " << line << endl << endl;
      cerr << "Required syntax:" << endl << endl;
      cerr << "#GENE CONVERSION" << endl;
      cerr << "<fraction> <average-length>" << endl << endl;
      cerr << "Example:" << endl << endl;
      cerr << "#GENE CONVERSION" << endl;
      cerr << "0.5 20" << endl << endl;
    }

  else if (type == 13) // fitness interaction
    {
      cerr << "ERROR (parameter file): invalid fitness interaction: " << line << endl << endl;
      cerr << "Required syntax:" << endl << endl;
      cerr << "#FITNESS INTERACTION" << endl;
      cerr << "<interaction-type>" << endl << endl;
      cerr << "Example:" << endl << endl;
      cerr << "#FITNESS INTERACTION" << endl;
      cerr << "a" << endl << endl;
    } // end of fitness interaction

  else if (type == 14) // environment
    {
      cerr << "ERROR  (parameter file): invalid environment:" << line << endl << endl;
      cerr << "Required syntax:" << endl << endl;
      cerr << "#ENVIRONMENTS" << endl;
      cerr << "<environment-id> <mut-type> <h> <s-modif> [<mut-type> <<h>> <s-modif>...]" << endl;
      cerr << "..." << endl << endl;
      cerr << "Example:" << endl << endl;
      cerr << "#ENVIRONMENTS" << endl;
      cerr << "e1 m1 0.5 -1.0 m2 0.0 0.5 m3 1.0 1.0" << endl << endl;
    } // end of environment

  exit(1);
} // end of method input_error()


void check_input_file(char* file)
{
  int mutation_types = 0;
  int mutation_rate  = 0;
  int interaction_type = 0; // added in SLaM v. 1.0
  int genomic_element_types = 0;
  int chromosome_organization = 0;
  int recombination_rate = 0;
  int generations = 0;
  int population = 0;
  int fitness_interaction = 0;
  int environments = 0;

  ifstream infile (file);
  // check if parameter file can be opened
  if (!infile.is_open()) { input_error(0, string(file)); } // the method input_error handles
        // different types of errors, here it is case "0"

  string line; string sub;

  get_line(infile, line); // calls a modified getline function that takes care of some comment
    // lines and white spaces at the beginning and end of lines.
  // read line by line
  while (!infile.eof()) // while not hitting end of file 'infile'
    {
      if (line.find('#') != string::npos) // check for start of an input section, denoted by
        // an initial '#'; if this line starts a new parameter section
	      {
        
          // 'mutation rate' section
          if (line.find("MUTATION RATE") != string::npos)
            {
              get_line(infile, line);
              while (line.find('#') == string::npos && !infile.eof()) // while not hitting next
              // section or end of file
                {
                  if (line.length() > 0) // if line is not empty
                    {
                      if (line.find_first_not_of("1234567890.e-") != string::npos)
                        {  input_error(1, line);  }
                      else
                        {  mutation_rate++;  } // mutation_rate now 1
                    } // end of if line is not empty
                  get_line(infile, line);
                } // end of while not hitting next section or end of file
            } // end of 'mutation rate' section

          // 'fitness interaction' section
          else if (line.find("FITNESS INTERACTION") != string::npos)
            {
              get_line(infile, line);
              while (line.find('#') == string::npos && !infile.eof()) // while not hitting next
                                                                     // section or end of file
                {
                  if (line.length() > 0) // if line is not empty
                    {
                      if (line.find_first_not_of("am") != string::npos)
                        {  input_error(13, line);  }
                      else
                        {  fitness_interaction++;  } // fitness_type now 1
                    } // end of if line is not empty
                  get_line(infile, line);
                } // end of while not hitting next section or end of file
            } // end of 'mutation rate' section

          // 'mutation types' section
          else if (line.find("MUTATION TYPES") != string::npos)
            {
              get_line(infile, line);
              while (line.find('#') == string::npos && !infile.eof()) // while not hitting next
                // section or end of file
                {
                  if (line.length() > 0) // line is not empty
                    {
                      int good = 1;
                      istringstream iss(line);
                      iss >> sub;
                      if (sub.compare(0, 1, "m") != 0)
                        { good = 0; }
                      sub.erase(0, 1);
                      if (sub.find_first_not_of("1234567890") != string::npos)
                        { good = 0; } // id
                      if (iss.eof()) // note that eof() here is the end of the string representing
                          // the current 'line'
                        { good = 0; }
                      iss >> sub;
                      if (sub.find_first_not_of("1234567890.-") != string::npos)
                        { good = 0; } // h
                      if (iss.eof()) { good = 0; }
                      iss >> sub;
                      if (sub.find_first_not_of("fge") != string::npos)
                        { good = 0; } // DFE-type
                      if (sub.compare("f")==0 || sub.compare("e")==0) // one parameter
                        {
                          if (iss.eof())
                            { good = 0; }
                          iss >> sub;
                          if (sub.find_first_not_of("1234567890.-") != string::npos)
                            { good = 0; }
                          if (!iss.eof()) { good = 0; }
                        }
                      if (sub.compare("g")==0) // two parameters
                        {
                          if (iss.eof())
                            { good = 0; }
                          iss >> sub;
                          if (sub.find_first_not_of("1234567890.-") != string::npos)
                            { good = 0; }
                          if (iss.eof())
                            { good = 0; }
                          iss >> sub;
                          if (sub.find_first_not_of("1234567890.-") != string::npos)
                            { good = 0; }
                          if (!iss.eof()) { good = 0; }
                        }
                      if (good == 0)
                        { input_error(2, line); }
                      else
                        { mutation_types++; } // increase counter of mytation types
                    } // end of if line is not empty
                  get_line(infile, line);
                } // end of while not hitting next section
            } // end of mutation types section
        
          // 'environments' section

          else if (line.find("ENVIRONMENTS") != string::npos)
            {
              cout << "Entering ENVIRONMENTS" << endl;
              get_line(infile, line);
              while (line.find('#') == string::npos && !infile.eof()) // while not hitting next
                // section or end of file
                {
                  cout << "ENVIRONMENTS: while not hitting next section or end of file." << endl;
                  if (line.length() > 0) // line is not empty
                    {
                      int good = 1;
                      istringstream iss(line);
                      iss >> sub;
                      if (sub.compare(0, 1, "e") != 0)
                        { good = 0; }
                      sub.erase(0, 1); // id (of environment, numeric part)
                      if (sub.find_first_not_of("1234567890") != string::npos)
                        { good = 0; }
                      if (iss.eof())
                        { good = 0; }
                      while(!iss.eof()) // while end of file (line) not reached
                        {
                          iss >> sub;
                          if(sub.compare(0, 1, "m") != 0) // mutation-type id not starting with
                              // "m"
                            { good = 0; }
                          sub.erase(0, 1); // mutation type id (numeric part)
                          if(sub.find_first_not_of("1234567890") != string::npos)
                            { good = 0; }
                          if(iss.eof())
                            { good = 0; }
                          iss >> sub; // h (dominance)
                          if(sub.find_first_not_of("1234567890.e-") != string::npos)
                            { good = 0; }
                          if(iss.eof())
                            { good = 0; }
                          iss >> sub; // s-modif (modifier of [mean] selection coefficient)
                          if(sub.find_first_not_of("1234567890.e-") != string::npos)
                            { good = 0; }
                        } // end of while end of file (line) not reached
                      if (good == 0)
                        { input_error(14, line); }
                      else
                        { environments++; } // increase environments counter
                    } // end of if line is not empty
                  get_line(infile, line);
                } // end of while not hitting next section or end of file
            } // end of environments section
           
          // 'genomic element types' section (only testing for presence and correctness of
          else if (line.find("GENOMIC ELEMENT TYPES") != string::npos)
            {
              get_line(infile, line);
              while (line.find('#') == string::npos && !infile.eof()) // while not hitting next
                // section or end of file
                {
                  if (line.length() > 0) // line is not emtpy
                    {
                      int good = 1;
                      istringstream iss(line);
                      iss >> sub;
                      if (sub.compare(0, 1, "g") != 0)
                        { good = 0; }
                      sub.erase(0, 1); // id
                      if (sub.find_first_not_of("1234567890") != string::npos)
                        { good = 0; }
                      if (iss.eof())
                        { good = 0; }
                      while (!iss.eof()) // while end of file of line not reached
                        {
                          iss >> sub;
                          if (sub.compare(0, 1, "m") != 0)
                            { good = 0; }
                          sub.erase(0, 1); // mutation type id (numeric part)
                          if (sub.find_first_not_of("1234567890") != string::npos)
                            { good = 0; }
                          if (iss.eof())
                            { good = 0; }
                          iss >> sub; // fraction
                          if (sub.find_first_not_of("1234567890e.") != string::npos)
                            { good = 0; }
                        } // end of while not hit end of file (i.e. end of current line)
                      if (good == 0)
                        { input_error(3, line); }
                      else
                        { genomic_element_types++; } // increase counter
                    } // end of if line is not emtpy
                  get_line(infile, line);
                } // end of while not hitting next section
            } // end of 'genomic element types' section
        
          // 'chromosome organization' section
          else if (line.find("CHROMOSOME ORGANIZATION") != string::npos)
            { get_line(infile, line);
              while (line.find('#') == string::npos && !infile.eof()) // while not hitting next
                // section or end of file
                {
                  if (line.length() > 0)
                    {
                      int good = 1;
                      istringstream iss(line);
                      iss >> sub;
                      if (sub.compare(0, 1, "g") != 0)
                        { good = 0; }
                      sub.erase(0, 1);
                      if (sub.find_first_not_of("1234567890") != string::npos)
                        { good = 0; } // id
                      if (iss.eof())
                        { good = 0; }
                      iss >> sub;
                      if (sub.find_first_not_of("1234567890e") != string::npos)
                        { good = 0; } // start
                      if (iss.eof())
                        { good = 0; }
                      iss >> sub;
                      if (sub.find_first_not_of("1234567890e") != string::npos)
                        { good = 0; } // end
                      if (!iss.eof())
                        { good = 0; }
                      if (good == 0)
                        { input_error(4, line); }
                      else
                        { chromosome_organization++; } // increase counter
                    } // end of if line is not empty
                  get_line(infile, line);
                } // end of while not hitting next section or end of file
            } // end of 'chromosome organization' section
        
          // 'recombination rate' section
          else if (line.find("RECOMBINATION RATE") != string::npos)
            {
              get_line(infile, line);
              while (line.find('#') == string::npos && !infile.eof())
                {
                  if (line.length() > 0) // if line is not empty
                    {
                      int good = 1;
                      istringstream iss(line);
                      iss >> sub;
                      if (sub.find_first_not_of("1234567890e") != string::npos)
                        { good = 0; } // end
                      if (iss.eof())
                        { good = 0; }
                      iss >> sub;
                      if (sub.find_first_not_of("1234567890e.-") != string::npos)
                        { good = 0; } // rate
                      if (!iss.eof())
                        { good = 0; }
                      if (good == 0)
                        { input_error(5, line); }
                      else { recombination_rate++; } // increase counter
                    } // end of if line is not empty
                  get_line(infile, line);
                } // end of while not hitting next section or end of file
            } // end of 'recombination rate' section
        
          // 'gene conversion' section
          else if (line.find("GENE CONVERSION") != string::npos)
            {
              get_line(infile, line);
              while (line.find('#') == string::npos && !infile.eof()) // while not hitting next
                  // section or end of file
                {
                  if (line.length() > 0) // if line is not emtpy
                    {
                      int good = 1;
                      istringstream iss(line); iss >> sub;
                      if (sub.find_first_not_of("1234567890e.-") != string::npos)
                        { good = 0; } // fraction
                      if (iss.eof())
                        { good = 0; }
                      iss >> sub;
                      if (sub.find_first_not_of("1234567890e.-") != string::npos)
                        { good = 0; } // average length
                      if (!iss.eof())
                        { good = 0; }
                      if (good == 0)
                        { input_error(12, line); }
                    } // end of if line is not empty
                      get_line(infile, line);
                } // end of while not hitting next section or end of file
            } // end of 'gene conversion' section
        
          // 'generations' section
          else if (line.find("GENERATIONS") != string::npos)
            {
              get_line(infile, line);
              while (line.find('#') == string::npos && !infile.eof()) // while not hitting next
                  // section or end of file
                {
                  if (line.length() > 0) // if line is not empty
                    {
                      int good = 1;
                      istringstream iss(line);
                      iss >> sub;
                      if (sub.find_first_not_of("1234567890e") != string::npos)
                        { good = 0; } // T
                      if (!iss.eof())
                        { good = 0; }
                      iss >> sub;
                      if (good == 0)
                        { input_error(6, line); }
                      else { generations++; } // increase counter
                    } // end of line is not empty
                  get_line(infile, line);
                } // end of while not hitting next section or end of file
            } // end of 'generations' section
        
          // 'demography and structure' section
          else if (line.find("DEMOGRAPHY AND STRUCTURE") != string::npos)
            {
              get_line(infile, line);
              while (line.find('#') == string::npos && !infile.eof()) // while not hitting next
                  // section or end of file
                {
                  if (line.length() > 0) // if line is not empty
                    {
                      int good = 1;
                      istringstream iss(line);
                      iss >> sub; // t(ime)
                      if (sub.find_first_not_of("1234567890e") != string::npos)
                        { good = 0; }
                      if (iss.eof())
                        { good = 0; }
                      iss >> sub;
                      if (sub.find_first_not_of("PSMNE") != string::npos)
                        { good = 0; } // event type
                    
                      if (sub.compare("P") == 0) // adding a new population; expect two or three
                        // positive integers (the third one for the source population is optional)
                        {
                          if (iss.eof())
                            { good = 0; }
                          iss >> sub; // p1
                          if (sub.compare(0, 1, "p") != 0)
                            { good = 0; }
                          sub.erase(0, 1);
                          if (sub.find_first_not_of("1234567890") != string::npos)
                            { good = 0; }
                          if (iss.eof())
                            { good = 0; }
                          iss >> sub; // N (size of new population)
                          if (sub.find_first_not_of("1234567890e") != string::npos)
                            { good = 0; }
                          if (!iss.eof()) // p2 (source population has been specified)
                            { // if not hit end of file
                              iss >> sub;
                              if (sub.compare(0, 1, "p") != 0)
                                { good = 0; }
                              sub.erase(0, 1);
                              if (sub.find_first_not_of("1234567890") != string::npos)
                                { good = 0; }
                              if (!iss.eof())
                                { good = 0; }
                            } // end of if not hit end of file
                          population++; // increase counter
                        } // end of adding a new population
                    
                      if (sub.compare("N")==0) // changing population size; expect two positive
                        // integers
                        {
                          if (iss.eof())
                            { good = 0; }
                          iss >> sub; // p (population identifier)
                          if (sub.compare(0, 1, "p") != 0)
                            { good = 0; }
                          sub.erase(0, 1);
                          if (sub.find_first_not_of("1234567890") != string::npos)
                            { good = 0; }
                          if (iss.eof())
                            { good = 0; }
                          iss >> sub; // N (new size)
                          if (sub.find_first_not_of("1234567890e") != string::npos)
                            { good = 0; }
                          if (!iss.eof())
                            { good = 0; }
                        } // end of changing population size

                      if (sub.compare("S")==0) // changing selfing rate; expect one positive integer and a double
                        {
                          if (iss.eof())
                              { good = 0; }
                          iss >> sub; // p (population identifier)
                          if (sub.compare(0, 1, "p") != 0)
                            { good = 0; }
                          sub.erase(0, 1);
                          if (sub.find_first_not_of("1234567890") != string::npos)
                            { good = 0; }
                          if (iss.eof())
                            { good = 0; }
                          iss >> sub; // sigma (selfing rate)
                          if (sub.find_first_not_of("1234567890.-e") != string::npos)
                            { good = 0; }
                          if (!iss.eof())
                            { good = 0; }
                        } // end of changing selfing rate

                      if (sub.compare("M")==0) // changing migration rate; two positive integers
                        // and a double
                        {
                          if (iss.eof())
                            { good = 0; }
                          iss >> sub; // p (id of target population)
                          if (sub.compare(0, 1, "p") != 0)
                            { good = 0; }
                          sub.erase(0, 1);
                          if (sub.find_first_not_of("1234567890") != string::npos)
                            { good = 0; }
                          if (iss.eof())
                            { good = 0; }
                          iss >> sub; // p (id of source population)
                          if (sub.compare(0, 1, "p") != 0)
                            { good = 0; }
                          sub.erase(0, 1);
                          if (sub.find_first_not_of("1234567890") != string::npos)
                            { good = 0; }
                          if (iss.eof())
                            { good = 0; }
                          iss >> sub; // M (new migration rate)
                          if (sub.find_first_not_of("1234567890.-e") != string::npos)
                            { good = 0; }
                          if (!iss.eof())
                            { good = 0; }
                        } // end of changing migration rate

                      if (sub.compare("E") == 0) // changing environment; time, "E", population id
                        // and environment id
                        {
                          if (iss.eof())
                            { good = 0; }
                          iss >> sub; // p(opulation id)
                          if (sub.compare(0, 1, "p") != 0)
                            { good = 0; }
                          sub.erase(0, 1);
                          if (sub.find_first_not_of("1234567890") != string::npos)
                            { good = 0; }
                          if (iss.eof())
                            { good = 0; }
                          iss >> sub; // e(nvironment id)
                          if (sub.compare(0, 1, "e") != 0)
                            { good = 0; }
                          sub.erase(0, 1);
                          if (sub.find_first_not_of("1234567890") != string::npos)
                            { good = 0; }
                          if (!iss.eof())
                            { good = 0; }
                        } // end of changing environment

                      if (good == 0)
                        { input_error(7, line); }
                    } // end of if line is not emtpy
                  get_line(infile, line);
                } // end of while not hitting next section or end of file
            } // end of 'demography and structure' section
        
          // 'output section
          else if (line.find("OUTPUT") != string::npos)
            {
              cout << "reading output specs" << endl; // test
              get_line(infile, line);
              while (line.find('#') == string::npos && !infile.eof())
                { // while not hitting next section or end of file
                  if (line.length() > 0) // if line is not emtpy
                    {
                      int good = 1;
                      istringstream iss(line);
                      iss >> sub;
                      if (sub.find_first_not_of("1234567890e") != string::npos)
                        { good = 0; } // t(ime)
                      if (iss.eof())
                        { good = 0; }
                      iss >> sub;
                      if (sub.find_first_not_of("ARFT") != string::npos)
                        { good = 0; } // event type
                      if (sub.compare("A") == 0) // output state of entire population; expect no
                        // parameter or a filename
                        {
                          if (!iss.eof()) // there is a filename
                            {
                              iss >> sub;
                              if (!iss.eof())
                                { good = 0; }
                            } // end of hitting end of file
                        } // end of if this is the "A" subtype of event
                      if (sub.compare("R")==0) // output random sample; expect two parameters
                        // and an optional flag specifying whether output should be given in the ms
                        // format
                        {
                          if (iss.eof())
                            { good = 0; }
                          iss >> sub; // p (population id)
                          if (sub.compare(0, 1, "p") != 0)
                            { good = 0; }
                          sub.erase(0, 1);
                          if (sub.find_first_not_of("1234567890") != string::npos)
                            { good = 0; }
                          if (iss.eof())
                            { good = 0; }
                          iss >> sub; // size of sample
                          if (sub.find_first_not_of("1234567890") != string::npos)
                            { good = 0; }
                          if (!iss.eof()) // one more argument, the ms option, is present
                            {
                              iss >> sub; // MS
                              if (sub != "MS")
                                { good = 0; }
                            } // end of check for optional MS option
                          if (!iss.eof())
                            { good = 0; }
                        } // end of if this is the "R" subtype of event
                      if (sub.compare("F") == 0) // output list of all fixed mutations; expect no
                        // parameter
                        {
                        if (!iss.eof())
                          { good = 0; }
                        } // end of if this is the "F" subtype of event
                      if (sub.compare("T") == 0) // output track of mutations of particular type;
                        // expect one parameter (id of mutation type)
                        {
                        if (iss.eof())
                          { good = 0; }
                        iss >> sub; // mutation type
                        if (sub.compare(0, 1, "m") != 0)
                          { good = 0; }
                        sub.erase(0, 1);
                        if (sub.find_first_not_of("1234567890") != string::npos)
                          { good = 0; }
                        if (!iss.eof())
                          { good = 0; }
                        } // end of if this is the "T" subtype of event
                      if (good == 0)
                        { input_error(8, line); }
                    } // end of if line is not empty
                  get_line(infile, line);
                } // end of while not hitting next section or end of file
            } // end of 'output' section
        
          // 'initialisation' section
          else if (line.find("INITIALIZATION") != string::npos)
            {
              get_line(infile, line);
              while (line.find('#') == string::npos && !infile.eof()) // while not hitting next 
                // section or end of file
                {
                  if (line.length() > 0) // if line is not emtpy
                    {
                      int good = 1;
                      istringstream iss(line);
                      iss >> sub; // name of initialisation file
                      if (!iss.eof())
                        { good = 0; }
                      if (good == 0)
                        { input_error(9, line); }
                      population++;
                    } // enf of if line is not empty
                  get_line(infile, line);
                } // end of while not hitting next section or end of file
            } // end of 'initialisation' section
          
          // 'seed' section
          else if (line.find("SEED") != string::npos)
            { 
              get_line(infile, line);
              while (line.find('#') == string::npos && !infile.eof()) // while not hitting end of
                // section or end of file
                {
                  if (line.length() > 0) // if line is not empty
                    {
                      int good = 1;
                      istringstream iss(line);
                      iss >> sub; // seed
                      if (sub.find_first_not_of("1234567890-") != string::npos) { good = 0; }
                      if (!iss.eof()) { good = 0; }
                      if (good == 0) { input_error(10, line); }
                    }
                  get_line(infile, line);
                } // end of if line is not empty
            } // end of 'seed' section
          
          // 'predetermined mutations' section; expect seven or nine entries, depending
            // on whether or not a partial sweep is desired
          else if (line.find("PREDETERMINED MUTATIONS") != string::npos)
            { 
              get_line(infile, line);
              while (line.find('#') == string::npos && !infile.eof()) // while not hitting next
                // section or end of file
                {
                  if (line.length() > 0) // if line is not empty
                    {
                      int good = 1;
                      istringstream iss(line);
                      iss >> sub; // time
                      if (sub.find_first_not_of("1234567890e") != string::npos)
                        { good = 0; }
                      if (iss.eof())
                        { good = 0; }
                      iss >> sub; // id (of mutation type)
                      if (sub.compare(0, 1, "m") != 0)
                        { good = 0; }
                      sub.erase(0, 1);
                      if (sub.find_first_not_of("1234567890") != string::npos)
                        { good = 0; }
                      if (iss.eof()) { good = 0; }
                      iss >> sub; // x (position)
                      if (sub.find_first_not_of("1234567890e") != string::npos)
                        { good = 0; }
                      if (iss.eof())
                        { good = 0; }
                      iss >> sub; // sub (population)
                      if (sub.compare(0, 1, "p") != 0)
                        { good = 0; }
                      sub.erase(0, 1);
                      if (sub.find_first_not_of("1234567890") != string::npos)
                        { good = 0; }
                      if (iss.eof()) { good = 0; }
                      iss >> sub; // nAA (number of homozygotes)
                      if (sub.find_first_not_of("1234567890") != string::npos)
                        { good = 0; }
                      if (iss.eof())
                        { good = 0; }
                      iss >> sub; // nAa (number of heterozygotes)
                      if (sub.find_first_not_of("1234567890") != string::npos)
                        { good = 0; }
                      iss >> sub; // linkage
                      if (sub.find_first_not_of("de") != string::npos)
                        { good = 0; }
                      if (!iss.eof()) // expect two more entries (partial sweep)
                        {
                          iss >> sub;
                          if (sub.find_first_not_of("P") != string::npos)
                            { good = 0; }
                          if (iss.eof())
                            { good = 0; }
                          iss >> sub; // freq (frequency at which allele turns neutral)
                          if (sub.find_first_not_of("1234567890.-e") != string::npos)
                            { good = 0; }
                        } // end of partial sweep subsection
                      if (!iss.eof())
                        { good = 0; }
                      if (good == 0)
                      { input_error(11, line); }
                    } // end of if line is not empty
                  get_line(infile, line);
                } // end of while not hitting next section or end of file
            } // end of 'predetermined mutations' section
          else // if none of the intended settings applies
            { input_error(-1, line); } // unknown parameter (section)
        } // end of if this line starts a new parameter section; this line did not start a new
            // input section; else: read next line
      else
        {
          get_line(infile, line);
        }
    } // end of while not hitting end of file 'infile'; end of file 'infile' now reached

  if (mutation_rate != 1)
    { input_error(1, string()); }
  if (mutation_types < 1)
    { input_error(2, string()); }
  if (genomic_element_types < 1)
    { input_error(3, string()); }
  if (chromosome_organization < 1)
    { input_error(4, string()); }
  if (recombination_rate < 1)
    { input_error(5, string()); }
  if (generations < 1)
    { input_error(6, string()); }
  if (population < 1)
    { input_error(-2, string()); }
  if (fitness_interaction != 1)
    { input_error(13, string()); }
} // end of method check_input_file()


void initialize_from_file(population& P, const char* file, chromosome& chr)
{
  // initialize population from file

  map<int,mutation> M;

  string line; 
  string sub; 

  ifstream infile (file);

  if (!infile.is_open()) { cerr << "ERROR (initialize): could not open initialization file" << endl; exit(1); }

  get_line(infile, line);

  while (line.find("Populations") == string::npos && !infile.eof()) { get_line(infile, line); } 

  get_line(infile, line);

  while (line.find("Mutations") == string::npos && !infile.eof())
    { 
      istringstream iss(line); iss >> sub; sub.erase(0, 1);  
      int i = atoi(sub.c_str()); iss >> sub;  
      int n = atoi(sub.c_str());
      P.add_subpopulation(i,n);
      get_line(infile, line);      
    }

  get_line(infile, line);

  while (line.find("Genomes") == string::npos && !infile.eof()) 
    {     
      istringstream iss(line); iss >> sub; 
      int   i = atoi(sub.c_str()); iss >> sub; sub.erase(0, 1); 
      int   t = atoi(sub.c_str()); iss >> sub; 
      int   x = atoi(sub.c_str())-1; iss >> sub; 
      float s = atof(sub.c_str());

      M.insert(pair<int,mutation>(i,mutation(t,x,s)));
      get_line(infile, line); 
    }

  get_line(infile, line);

  while (!infile.eof()) // while end of file not reached
    {
      istringstream iss(line); iss >> sub; sub.erase(0, 1);
      int pos = sub.find_first_of(":"); 
      int p = atoi(sub.substr(0,pos+1).c_str()); sub.erase(0,pos+1);  
      int i = atoi(sub.c_str());

      while (iss >> sub) 
        {
          int id = atoi(sub.c_str());
          P.find(p)->second.G_parent[i-1].push_back(M.find(id)->second);
        }
      get_line(infile, line);
    } // end of while end of file not reached

  for (P.it = P.begin(); P.it != P.end(); P.it++)
    {
      P.it->second.update_fitness(chr);
    }
} // end of initialize_from_file()

/*
  Initialises the population P using parameter values to be read from file. Initialises chromosome
  chr, as well as demography and structure events E, output events O, user-defined mutations  IM
  that will be introduced, and mutations PS undergoing partial sweeps. FI denotes the fitness
  interaction
*/
void initialize(population& P, char* file, chromosome& chr, int& T, char& FI, multimap<int,event>& E, multimap<int,event>& O, multimap<int,introduced_mutation>& IM, vector<partial_sweep>& PS, vector<string>& parameters)
{
  string line; 
  string sub; 

  ifstream infile (file);

  // setting the seed based on the current time and the process id
  long pid = getpid(); // obtain process id from operating system
  time_t *tp, t;
  tp = &t;
  time(tp);
  t += pid;
  int seed = t;

  get_line(infile, line);

  while (!infile.eof()) // while not hitting end of file
    {
      if (line.find('#') != string::npos) // if hitting start of an input section
        {
          // mutation rate section
          if (line.find("MUTATION RATE") != string::npos)
            {
              get_line(infile, line);
              parameters.push_back("#MUTATION RATE");
              while (line.find('#') == string::npos && !infile.eof()) // while not hitting new
                // input section
                {
                  if (line.length() > 0) // if line is not empty
                    {
                      parameters.push_back(line);
                      istringstream iss(line);
                      iss >> sub;
                      chr.M = atof(sub.c_str()); // initialise overall mutation rate
                    } // end of if line is not empty
                  get_line(infile, line);
                } // end of while not hitting new input section
            } // end of mutation rate section

          // mutation types section
          else if (line.find("MUTATION TYPES") != string::npos)
            {
              get_line(infile, line);
              parameters.push_back("#MUTATION TYPES");
              while (line.find('#') == string::npos && !infile.eof()) // while not hitting new
                // input section
                {
                  cout << "mutation types: " << line << endl; // test
                  if (line.length() > 0) // if line is not empty
                    {
                      parameters.push_back(line);
                      // FORMAT: i h t p
                        // i: identifier; h: dominance; t [DFE-type]: f(ixed), e(xponential), or
                        // g(amma); p [DFE parameters]: one for f and e (mean s), two for g (mean
                        // s, shape)
                      int i; float h; char t;
                      vector<double> p; // one or two entries
                      istringstream iss(line);

                      iss >> sub;
                      sub.erase(0, 1); // getting rid of "m" for mutation type
                      i = atoi(sub.c_str()); // initialise mutation type id
                      if (chr.mutation_types.count(i) > 0) // assessed based on the key
                        {
                          cerr << "ERROR (initialize): mutation type "<< i << " already defined" << endl; exit(1);
                        }

                      iss >> sub;
                      h = atof(sub.c_str()); // initialise dominance coefficient

                      iss >> sub;
                      t = sub.at(0); // initialise DFE-type

                      while (iss >> sub) // there are one or two values to be read
                        {
                          p.push_back(atof(sub.c_str()));
                          // cout << "iterating over iss >> sub" << endl; // test
                        }
                      chr.mutation_types.insert(pair<int,mutation_type>(i,mutation_type(h, t, p)));
                    } // end of if line is not empty
                  get_line(infile, line);
                } // end of while not hitting new input section
            } // end of mutation types section

          // fitness interaction section
          else if (line.find("FITNESS INTERACTION") != string::npos)
            {
              // cout << "entering fitness interaction (initialize)" << endl; // test
              get_line(infile, line);
              parameters.push_back("FITNESS INTERACTION");
              while (line.find('#') == string::npos && !infile.eof()) // while not hitting next
                // input section or end of file
                {
                  if (line.length() > 0) // if line is not empty
                    {
                      parameters.push_back(line);
                      istringstream iss(line);
                      iss >> sub; // 'a' for additive; 'm' for multiplicative
                      FI = sub.at(0);
                    } // end of if line is not empty
                  get_line(infile, line);
                } // end of while not hitting next input section or end of file
              // cout << "leaving fitness interaction (initialize)" << endl; // test
            } // end of fitness interaction section

          // genomic element types section
          else if (line.find("GENOMIC ELEMENT TYPES") != string::npos)
            {
              // cout << "entering genomic element types (initialize)" << endl; // test
              get_line(infile, line);
              parameters.push_back("#GENOMIC ELEMENT TYPES");
              while (line.find('#') == string::npos && !infile.eof()) // while not hitting next
                // input section
                {
                  if (line.length() > 0) // if line is not empty
                    {
                      parameters.push_back(line);
                      // FORMAT: i m1 g1 [m2 g2 ...]
                        // (identifier, mut type class, fraction)
                      int i; vector<int> m; vector<double> g;
                      istringstream iss(line);
                      iss >> sub;
                      sub.erase(0, 1); // erasing leading "g"
                      i = atoi(sub.c_str()); // initialising genomic-element-type id
                      while (iss >> sub)
                        { // while there are multiples of two parameters to be read
                          sub.erase(0, 1); // erasing leading "m"
                          m.push_back(atoi(sub.c_str())); // initialising mutation-type id
                          iss >> sub;
                          g.push_back(atof(sub.c_str())); // initialising fraction of mutations
                            // made up of the current mutation type
                        } // end of while there are multiples of two parameters to be read
                      if (chr.genomic_element_types.count(i) > 0)
                        {
                          cerr << "ERROR (initialize): genomic element type "<< i << " already defined" << endl;
                          exit(1);
                        }
                      chr.genomic_element_types.insert(pair<int,genomic_element_type>(i,genomic_element_type(m, g))); // m is the vector of mutation types, and g
                          // the vector for fractions they make up of all mutations
                    } // end of if line is not empty
                  get_line(infile, line);
                } // end of while not hitting next input section
              // cout << "leaving genomic element types (initialize)" << endl; // test
            } // end of genomic element types section

          // chromosome organization section
          else if (line.find("CHROMOSOME ORGANIZATION") != string::npos)
            {
              get_line(infile, line);
              parameters.push_back("#CHROMOSOME ORGANIZATION");
              while (line.find('#') == string::npos && !infile.eof()) // while not hitting next
                // input section or end of file
                {
                  if (line.length() > 0) // if line is not empty
                    {
                      parameters.push_back(line);
                        // FORMAT: i s e
                          // (genomic element type identifier, start, end)
                        int i, s, e;
                        istringstream iss(line);
                        iss >> sub;
                        sub.erase(0, 1); // erasing leading "g"
                        i = atoi(sub.c_str());
                        iss >> sub;
                        s = (int) atof(sub.c_str()) - 1; // start; internally, the chromosome
                          // starts at position 0
                        iss >> sub;
                        e = (int) atof(sub.c_str()) - 1; // internal start is at 0
                        chr.push_back(genomic_element(i, s, e));
                    } // end of line is not empty
                  get_line(infile, line);
                } // end of while not hitting next input section or end of file
              // cout << "leaving chromosome organization (initialize)" << endl; // test
            } // end of chromosome organization section

          // recombination rate section
          else if (line.find("RECOMBINATION RATE") != string::npos)
            {
              get_line(infile, line);
              parameters.push_back("#RECOMBINATION RATE");
              while (line.find('#') == string::npos && !infile.eof()) // while not hitting next
                // input section or end of file
                {
                  if (line.length() > 0) // if line is not empty
                    {
                      parameters.push_back(line);
                      // FORMAT: x r
                        // (interval end, rec rate in events per bp)
                      int x; double r;
                      istringstream iss(line);
                      iss >> sub;
                      x = (int) atof(sub.c_str()) - 1; // internal start is at 0
                      iss >> sub;
                      r = atof(sub.c_str());
                      chr.rec_x.push_back(x); //
                      chr.rec_r.push_back(r);
                    } // end of if line is not empty
                  get_line(infile, line); } }

          // gene conversion section
          else if (line.find("GENE CONVERSION") != string::npos)
            {
              get_line(infile, line);
              parameters.push_back("#GENE CONVERSION");
              while (line.find('#') == string::npos && !infile.eof())
                {
                  if (line.length() > 0) // if line is not empty
                    {
                      parameters.push_back(line);
                      // FORMAT: G_f G_l
                        // (gene conversion fraction, average stretch length in bp)
                      istringstream iss(line);
                      iss >> sub;
                      chr.G_f = atof(sub.c_str());
                      iss >> sub;
                      chr.G_l = atof(sub.c_str());
                    } // end of if line is not empty
                  get_line(infile, line);
                } // end of while not hitting next input section or end of file
            } // end of gene conversion section

          // generations
          else if (line.find("GENERATIONS") != string::npos)
            {
              get_line(infile, line);
              parameters.push_back("#GENERATIONS");
              while (line.find('#') == string::npos && !infile.eof()) // while not hitting next
                // input section or end of file
                {
                  if (line.length() > 0) // if line is not empty
                    {
                      parameters.push_back(line);
                      istringstream iss(line);
                      iss >> sub;
                      T = (int)atof(sub.c_str());
                    } // end of if line is not empty
                  get_line(infile, line);
                } // end of while not hitting next input section or end of file
            } // end of generations

          // demography and structure section
          else if (line.find("DEMOGRAPHY AND STRUCTURE") != string::npos)
            {
              get_line(infile, line);
              parameters.push_back("#DEMOGRAPHY AND STRUCTURE");
              while (line.find('#') == string::npos && !infile.eof()) // while not hitting next
                // input section or end of file
                {
                  if (line.length() > 0) // if line is not emtpy
                    {
                      parameters.push_back(line);
                      // FORMAT: t c s
                        // (time, event_type, event_parameters)
                      int t; char c;
                      vector<string> s;
                      istringstream iss(line);
                      iss >> sub;
                      t = (int)atof(sub.c_str());
                      iss >> sub;
                      c = sub.at(0); // event type can be P, N, M, S, or E
                      while (iss >> sub) // while there are more entries on that line
                        {
                          s.push_back(sub.c_str()); // adding event parameters one by one
                        } // end of line reached
                      E.insert(pair<int,event>(t,event(c, s))); // add demography and structure
                        // event; note that time t is the key, and the collection is a multimap
                    } // end of if line is not empty
                  get_line(infile, line);
                } // end of while not hitting next input section or end of file
              // cout << "leaving demography and structure (initialize)" << endl; // test
            } // end of demography and structure section

          // environment section
          else if (line.find("ENVIRONMENTS") != string::npos)
            {
              cout << "entering environment (initialize)" << endl; // test
              get_line(infile, line);
              parameters.push_back("#ENVIRONMENTS");
              while (line.find('#') == string::npos && !infile.eof()) // while not hitting next
                // input section
                {
                  if (line.length() > 0) // if line is not emtpy
                    {
                      parameters.push_back(line);
                      // FORMAT: i m1 h1 s-modif1 [m2 h2 s-modif2 ...]
                        // (identifier, mut type, dominance coeff., modifier of (mean) sel. coeff
                        // , [mut type, dominance coeff., modifier of (mean) sel. coeff, ...])
                      int i; vector<int> m; vector<double> h; vector<double> smodif;
                      istringstream iss(line);

                      iss >> sub;
                      sub.erase(0, 1); // erasing leading "e"
                      i = atoi(sub.c_str()); // initialising environemt id
                      while (iss >> sub)
                        { // while there are multiples of three parameters to be read
                            sub.erase(0, 1); // erasing leading "m"
                            m.push_back(atoi(sub.c_str())); // initialising mutation-type id
                            iss >> sub;
                            h.push_back(atof(sub.c_str())); // initialising dominance coefficient
                            iss >> sub;
                            smodif.push_back(atof(sub.c_str())); // initialising modifier of
                              // selection coefficient
                        } // end of while there are multiples of three parameters to be read

                    } // end of if line is not empty
                } // end of while not hitting next input section
            } // end of environment section

          // output section
          else if (line.find("OUTPUT") != string::npos)
            {
              cout << "entering output (initialize)" << endl; // test
              get_line(infile, line);
              parameters.push_back("#OUTPUT");
              while (line.find('#') == string::npos && !infile.eof()) // while not hitting next
                // input section or end of file
                {
                  if (line.length() > 0)  // if line is not empty
                    {
                      parameters.push_back(line);
                      // FORMAT: t c s
                        // (time, event_type, event_paramaters)
                      int t; char c; vector<string> s;
                      istringstream iss(line);
                      iss >> sub;
                      t = (int)atof(sub.c_str());
                      iss >> sub;
                      c = sub.at(0); // envent type can be A, R, F
                      while (iss >> sub) // while there are more entries on that line
                        {
                          s.push_back(sub.c_str()); // adding event parameters one by one
                        } // end of line reached
                      O.insert(pair<int,event>(t,event(c, s))); // add ouptut event; note that
                        // time is the key, and the collection is a multimap
                    } // end of if line is not empty
                  get_line(infile, line);
                } // end of while not hitting next input section or end of file
              cout << "leaving output (initialize)" << endl; // test
            } // end of output section

          // predetermined mutations
          else if (line.find("PREDETERMINED MUTATIONS") != string::npos)
            {
              // cout << "entering predetermined mutations (initialize)" << endl; // test
              get_line(infile, line);
              parameters.push_back("#PREDETERMINED MUTATIONS");
              while (line.find('#') == string::npos && !infile.eof()) // while not hitting next
                // input section or end of file
                {
                  if (line.length() > 0) // if line is not empty
                    {
                      parameters.push_back(line);
                      // FORMAT: time t x i nAA nAa l [P p]
                        // (time, type, position, subpopulation, number of homozygous carriers,
                        // number of heterozygous carriers, linkage-equilibrium flag, [partial
                        // selective sweep: final frequency])
                      int time, t, x, i, nAA, nAa;
                      char l; // 'e' for linkage equilibrium, 'd' for linkage disequilibrium
                      float p;
                      istringstream iss(line);

                      iss >> sub;
                      time = (int)atof(sub.c_str());
                      iss >> sub;
                      sub.erase(0, 1); // erasing leading "m"
                      t = atoi(sub.c_str()); // mutation type
                      iss >> sub;
                      x = (int)atof(sub.c_str()) - 1; // position; internally starting at 0
                      iss >> sub;
                      sub.erase(0, 1); // erasing leading "p"
                      i = atoi(sub.c_str()); // subpopulation
                      iss >> sub;
                      nAA  = (int)atof(sub.c_str()); // number of homozygous carriers
                      iss >> sub;
                      nAa  = (int)atof(sub.c_str()); // number of heterozygous carriers
                      iss >> sub;
                      l = sub.at(0); // linkage information ("e" for LE, or "d" for LD)
                      introduced_mutation M(t, x, i, nAA, nAa, l);
                      IM.insert(pair<int,introduced_mutation>(time, M)); // time is the key;
                        // multiple events (mutations) can be associated with the same time

                      while (iss >> sub) // read final sweep frequency if present
                        {
                          if (sub.find('P') != string::npos) // if there is a 'P' in what remains
                            {
                              iss >> sub;
                              float p = atof(sub.c_str()); // assign the final frequency
                              PS.push_back(partial_sweep(t, x, p)); // time, position, final
                                // frequency
                            } // end of if there is a 'P' left
                        } // end of reading final sweep frequency
                    } // end of if line is empty
                  get_line(infile, line);
                } // end of while not hitting next input section or end of file
              // cout << "leaving predetermined mutations (initialize)" << endl; // test
            } // end of predetermined mutations section

          // seed section
          else if (line.find("SEED") != string::npos)
            {
              // cout << "entering seed (initialize)" << endl; // test
              get_line(infile, line);
              while (line.find('#') == string::npos && !infile.eof())
                {
                  if (line.length() > 0)
                    {
                      istringstream iss(line);
                      iss >> sub;
                      seed = atoi(sub.c_str());
                    }
                  get_line(infile, line);
                }
              // cout << "leaving seed (initialize)" << endl; // test
            } // end of seed section

          // initialization section
          else if (line.find("INITIALIZATION") != string::npos)
            {
              get_line(infile, line);
              parameters.push_back("#INITIALIZATION");
              while (line.find('#') == string::npos && !infile.eof()) // while not hitting next
                // section or end of file
                {
                  if (line.length() > 0) // if line is not empty
                    {
                      parameters.push_back(line);
                      istringstream iss(line);
                      iss >> sub;
                      initialize_from_file(P, sub.c_str(), chr);
                    } // end of if line is not empty
                  get_line(infile, line);
                } // end of while not hitting next section or end of file
            } // end of initialization section
        } // this line was not the start of a new input section
      else
        {
          get_line(infile, line);
          cout << "getting new line:" << endl; // test
          // TODO: GO ON HERE, testing error (endless loop)
          cout << line << endl;
        }
    // TODO: Introduce a check for the presence of lines that start with # and have not been
      // handled as intended.
      // cout << "not hitting end of file" << endl; // test
    } // end of while not hitting end of file


  // initialize chromosome

  chr.initialize_rng();

  // initialize rng

  rng = gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(rng, (long)seed);

  parameters.push_back("#SEED");
  stringstream ss; ss << seed;
  parameters.push_back(ss.str());

  // parameter output

  for (int i=0; i<P.parameters.size(); i++) { cout << parameters[i] << endl; }
} // end of method initialize()



int main(int argc,char *argv[])
{
  // initialize simulation parameters

  if (argc <= 1)
    {
      cerr << "usage: slim <parameter file>" << endl;
      exit(1);
    }

  char* input_file = argv[1];

  cout << "Starting 'check_input_file()" << endl; // test
  check_input_file(input_file); // calls method implemented above

  int T; // maximum number of generations
  chromosome chr; // chromosome
  char FI; // type of fitness interaction ('a' for additive; 'm' for multiplicative)

  population P;
  map<int, subpopulation>::iterator itP;

  P.parameters.push_back("#INPUT PARAMETER FILE");
  P.parameters.push_back(input_file);

  // demography and structure events

  multimap<int,event> E; // multiple elements can have identical keys, where the key is time;
    // events are P (adding new subpopulation), N (changing population size), M (chaning
    // migration rate), S (changing selfing rate), E (changing environment)
  multimap<int,event>::iterator itE;

  // output events (time, output)
  
  multimap<int,event> O; // multiple elements can have identical keys, where the key is time;
    // events are A (output state of entire population), R (output random sample from sub-
    // population), F (output list of all fixed mutations)
  multimap<int,event>::iterator itO;

  // user-defined mutations that will be introduced (time, mutation)

  multimap<int,introduced_mutation> IM; // multiple elements can have identical keys, where the key
    // is time
  multimap<int,introduced_mutation>::iterator itIM;

  // environments (all except the reference environment)

  map<int,environment> EV; // unique keys
  map<int,environment>::iterator itEV;

  // tracked mutation-types

  vector<int> TM; 

  // mutations undergoing partial sweeps

  vector<partial_sweep> PS;

  // GO ON HERE: understand, and extend.
  initialize(P, input_file, chr, T, FI, E, O, IM, PS, P.parameters);
 
  // evolve over t generations

  for (int g = 1; g <= T; g++)
    { 
      // execute demographic and substructure events in this generation 

      pair<multimap<int,event>::iterator,multimap<int,event>::iterator> rangeE = E.equal_range(g);
      for (itE = rangeE.first; itE != rangeE.second; itE++)
        {
          P.execute_event(itE->second, g, chr, TM);
        }
   
      // evolve all subpopulations

      for (itP = P.begin(); itP != P.end(); itP++) { P.evolve_subpopulation(itP->first,chr); }     
            
      // introduce user-defined mutations
        
      pair<multimap<int,introduced_mutation>::iterator,multimap<int,introduced_mutation>::iterator> rangeIM = IM.equal_range(g);
      for (itIM = rangeIM.first; itIM != rangeIM.second; itIM++) {
          // test
          // P.print_all(chr);
          
          P.introduce_mutation(itIM->second,chr);
          
          // test
          // P.print_all(chr);
      }

      // execute output events

      pair<multimap<int,event>::iterator,multimap<int,event>::iterator> rangeO = O.equal_range(g);
      for (itO = rangeO.first; itO != rangeO.second; itO++) { P.execute_event(itO->second,g,chr,TM); }

      // track particular mutation-types and set s=0 for partial sweeps when completed
      
      if (TM.size()>0 || PS.size()>0) { P.track_mutations(g,TM,PS,chr); }

      // swap generations

      P.swap_generations(g,chr);   
    }
}
