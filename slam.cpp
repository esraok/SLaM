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
// xx/xx/2014: Changed the mutation regime such that if a new non-neutral mutation is introduced
  // at a position at which a non-neutral mutation preexists, the preexisting mutation is assigned
  // a selection coefficient of zero (irrespective of its original mutation type) unless the
  // mutation to be introduced and the preexisting mutation are identical in state (i.e. they
  // have the same mutation type and identical selection coefficients in the reference
  // environment). This regime allows for recurrent mutation to the very same (in terms of state)
  // allele, but prevents more than two different non-neutral alleles (in terms of state) at any
  // site across the entire population. The regime grants a consistent handling of dominance, but
  // introduces a potential bias in the migration rate, because some mutations are not allowed.
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

  // default constructor
  mutation(void) { ; }

  // alternative constructor
  mutation(int T, int X, float S) 
  { 
    t = T;
    x = X;
    s = S;
  }
}; // end of class mutation


bool operator< (const mutation &M1, const mutation &M2)
// mutations are compared based on their physical position
{
  return M1.x < M2.x;
};


bool operator== (const mutation &M1, const mutation &M2)
// mutations are equal if their position, type and selection coefficient are identical
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
    } // end of constructor method

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
    h = H; // dominance coefficient; valid in the reference environment, subject to overriding in
      // user-defined environments
    d = D; // distribution of fitness effects (DFE); valid independently of environment
    p = P; // vector of parameters specifying the DFE; valid independently of environment

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

    // drawing the selection coefficient from the appropriate distribution of fitness effects (DFE)

    switch(d)
      {
        case 'f': return p[0];
        case 'g': return gsl_ran_gamma(rng, p[1], p[0]/p[1]); // If the mean is negative,
          // gsl_ran_gamma returns negative values.
        case 'e': return gsl_ran_exponential(rng, p[0]); // Note that if the mean is negative,
          // gsl_ran_exponential returns negative values.
        default: exit(1);
      } // end of switch
  } // end of draw_s()

}; // end of class mutation_type


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
  // a genomic element type is specified by a vector of the mutation type identifiers of all
  // mutation types than can occur in such elements, and a vector of their relative fractions.
  // examples: exon, intron, utr, intergenic, etc.

private:

  gsl_ran_discrete_t* LT;

public:

  vector<int>    m; // identifiers of mutation types that can occur in this element type
  vector<double> g; // relative fractions of each mutation type

  genomic_element_type(vector<int> M, vector<double> G)
  {
    m = M;
    g = G;  

    if (m.size() != g.size())
      { exit(1); }
    double A[m.size()]; // m.size() is the number of mutation types that can occur
    for (int i = 0; i < m.size(); i++)
      { A[i] = g[i]; } // weighing mutation types by their relative fractions
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

  gsl_ran_discrete_t* LT_M; // mutation (weighting genomic elements by their length)
  gsl_ran_discrete_t* LT_R; // recombination (weighting recombination strata by their length)

public:

  map<int,mutation_type>        mutation_types; // storing all possible mutation-types, the key
    // being the id of the mutation-type
  map<int,genomic_element_type> genomic_element_types; // storing the genomic element type that
    // the genomic elements of this chromosome belong to
  vector<int>                   rec_x; // vector of end points of strata of a given recomb. rate
  vector<double>                rec_r; // vector of recombination rates for each stratum

  map<int,vector<double>>       seg_nonneutr_mut; // storing the mutation-type id and selection
    // coefficient of all non-neutral mutations that are segregating in the entire population, the
    // key being their physical position

  int    L;   // length of chromosome
  double M;   // overall mutation rate, being built here as the product of the per-base pair rate
    // times the cumulative length of all genomic elements making up the chromosome; initally
    // assigned the per-base pair rate
  double R;   // overall recombination rate, being built here as the cumulative sum of the products
    // of the recombination rates and lengths of successive recombination strata in this
    // chromosome; initially assigned the per-base pair rate
  double G_f; // gene conversion fraction
  double G_l; // average stretch length

  void initialize_rng()
  {
    if (size() == 0) // if there are no genomic elements (initialisation must have failed)
      {
        cerr << "ERROR (initialize): empty chromosome" << endl;
        exit(1);
      }
    if (rec_r.size() == 0) // if there are no recombination rates (initialisation must have failed)
      {
        cerr << "ERROR (initialize): recombination rate not specified" << endl;
        exit(1);
      }
    if (!(M >= 0)) // initially, M is assigned the per-base pair recombination rate
      {
        cerr << "ERROR (initialize): invalid mutation rate" << endl;
        exit(1);
      }

    L = 0; // cumulative length of chromosome visited over the course of the algorithm

    // iterate over genomic elements

    for (int i = 0; i < size(); i++) // size(): the number of genomic elements that make up this
      // chromosome
      {
        if (genomic_element_types.count(operator[](i).i) == 0) // the first 'i' is the iterator,
          // the second 'i' is the identifier of the genomic element; if there are 0 genomic
          // element types with identifier i for the genomic element i
          {
            cerr << "ERROR (initialize): genomic element type " << operator[](i).i << " not defined" << endl;
            exit(1);
          }
      } // end of for each genomic element

    // checking if mutation types are valid

    // iterating over genomic element types
    for (map<int,genomic_element_type>::iterator it = genomic_element_types.begin(); it !=genomic_element_types.end(); it++)
      {
        // iterating over mutation types that can occur in the current genomic element type
        for (int j = 0; j < it->second.m.size(); j++)
          {
            if (mutation_types.count(it->second.m[j]) == 0)
              {
                cerr << "ERROR (initialize): mutation type " << it->second.m[j] << " not defined" << endl;
                exit(1);
              }
          }
      } // end of iterating over genomic element types

    // initialising the look-up table for recombination, and the overall recombination rate

    double A[size()]; // size() is the number of genomic elements that make up the chromosome
    int l = 0; // cumulative length of genomic elements
    // iterating over genomic elements
    for (int i = 0; i < size(); i++)
      { 
        if (operator[](i).e > L)
          {
            L = operator[](i).e; // shifting the curser to the end of the current genomic element
          }
        int l_i = operator[](i).e - operator[](i).s + 1.0; // the length of the current element
        A[i] = (double)l_i; // weight by length of current element
        l += l_i;
      } // end of iterating over genomic elements

    LT_M = gsl_ran_discrete_preproc(size(), A); // to draw genomic elements weighted by their
      // length
    M = M*(double)l; // scaling the per-base pair rate by the cumulative length of all genomic
      // elements

    // initialising the look-up table for recombination, and the overall recombination rate

    double B[rec_r.size()]; // rec_r.size() is the number of recombination strata
    B[0] = rec_r[0]*(double)rec_x[0]; // rec. rate of first stratum times length of first stratum
    R += B[0]; // upsate overall recombination rate
    // iterating over remaining recombination strata
    for (int i = 1; i < rec_r.size(); i++)
      { 
        B[i] = rec_r[i]*(double)(rec_x[i] - rec_x[i-1]); // weight of the ith recombination stratum
        R += B[i];
        if (rec_x[i] > L) // L is the cumulative length of all genomic elements; at this point, it
          // equals the length after all defined genomic elements have been visited; note that
          // the whole chromosome does not need to be covered by genomic elements, which is why
          // this update of L here is needed
        { L = rec_x[i]; }
      }
    LT_R = gsl_ran_discrete_preproc(rec_r.size(), B);
  } // end of method initialize_rng()


  int draw_n_mut()
  {
    // draw the total number of mutations
    return gsl_ran_poisson(rng, M); // M is the overall mutation rate
  } // end of method draw_n_mut()
 

  void draw_new_mut(map<int,mutation>& M_map)

  // draw a new mutation and introduce it at an appropriately drawn random position x (accounting
    // for genomic element types); a non-neutral mutation at position x is only allowed if no other
    // non-neutral mutation – except one identical in state with the one to be introduced – is
    // segregating at x in the entire population. Moreover, in any haploid genome, no more than one
    // mutation can exist at any given position. The new mutation will replace possible preexisting
    // mutations at that position. That is, if a mutation is already present in M_map at position x
    // where another one is about to be introduced, the previously existing mutation is replaced

  // M_map is a map of mutations to be modified such as to include the new mutation, if applicable
    // i.e. if the mutation is permitted; the key is the physical position

  {
    int g = gsl_ran_discrete(rng, LT_M); // draw genomic element according to weights given by
      // length of element relative to the sum of lenghts of all elements
    genomic_element_type ge_type = genomic_element_types.find(operator[](g).i)->second; // get the genomic element type of the genomic element g

    int mut_type_id = ge_type.draw_mutation_type(); // drawing mutation type id according to
      // relative weight of that type among all types allowed in genomic element type ge_type
    mutation_type mut_type = mutation_types.find(mut_type_id)->second; // mutation type
      // corresponding to mutation type id mut_type_id

    // position of new mutation, uniformly chosen over the length of the genomic element
    int x = operator[](g).s + gsl_rng_uniform_int(rng, operator[](g).e - operator[](g).s + 1);

    float s = mut_type.draw_s(); // selection coefficient

    // find potentially segregating non-neutral mutation at position x and point iterator to it
    map<int,vector<double>>::iterator sm_it = seg_nonneutr_mut.find(x);

    if (sm_it == seg_nonneutr_mut.end()) // if there is no non-neutral mutation segregating at
      // position x
      {
        // the new mutation can be added, potentially replacing a previously existing mutation at
          // the same position x
        M_map.erase(x);
        M_map.insert(pair<int,mutation>(x,mutation(mut_type_id, x, s)));
      }
    else // there is a non-neutral mutation segregating at position x
      {
        if (mut_type.p[0] == 0.0) // the mutation-type of the mutation to be introduced is neutral
          // in the reference environment (and, hence, in the entire population)
          {
            // the new neutral mutation can be added, potentially replacing a previously existing
              // mutation at the same position x
            M_map.erase(x);
            M_map.insert(pair<int,mutation>(x,mutation(mut_type_id, x, s)));
          }
        if (mut_type.t == (int) sm_it->second[0] && s == sm_it->second[1]) // the mutation to be
          // introduced is identical in state (in terms of mutation-type and selection coefficient
          // in the reference environment) with the one present at position x; this is a recurrent
          // mutation to the same allele
          {
            // the mutation is permitted, potentially replacing a previously existing mutation at
              // the sampe position x
            M_map.erase(x);
            M_map.insert(pair<int,mutation>(x,mutation(mut_type_id, x, s)));
          }
        // else: the new mutation is neither neutral nor identical in state with the segregating
          // non-neutral mutation; it cannot be introduced (which results in a downward-bias
          // of the mutation rate that is accepted here)
      } // end of there is a non-neutral mutation segregating at position x

    // OLD: return mutation(mut_type_id, x, s); // mutation type, position, selection coefficient

  } // end of method draw_new_mut()


  vector<int> draw_breakpoints()
  {
    // draw recombination breakpoints

    vector<int> r; // vector of breakpoints, given as physical positions

    int nr = gsl_ran_poisson(rng, R); // draw the random number of recombination events

    // place these events along the genome

    // iterating over recombination events
    for (int i = 0; i < nr; i++)
      {
        int x = 0;
        int j = gsl_ran_discrete(rng, LT_R);

        if (j == 0)
          {
            x = gsl_rng_uniform_int(rng, rec_x[j]);
          }
        else
          {
            x = rec_x[j - 1] + gsl_rng_uniform_int(rng, rec_x[j] - rec_x[j - 1]);
          }

        r.push_back(x);

        if (gsl_rng_uniform(rng) < G_f) // recombination results in gene conversion
          {
            int x2 = x + gsl_ran_geometric(rng, 1.0/G_l);
            r.push_back(x2);
          } // end of recombination results in gene conversion

      } // end of iterating over recombination events

    return r;

  } // end of method draw_breakpoints()

}; // end of class chromosome


class polymorphism
{
public:

  int   i; // mutation id
  int   t; // mutation type
  float s; // selection coefficient in the reference environment
  // int   n; // prevalence in the entire population (i.e. summing over all subpopulations existing
    // at the time of sampling)

  map<int,int> n; // map of prevalences in each subpopulation; the key is the subpopulation
    // identifier

  // polymorphism(int I, int T, float S, int N)
  polymorphism(int I, int T, float S, map<int,int> N)
  {
    i = I;
    t = T;
    s = S;
    n = N;
  }

  /*
  Computes total prevalence (i.e. across all subpopulations, in the entire population)
  */
  int calc_ntot()
  {
    map<int,int>::iterator nit;
    int ntot = 0; // total prevalence (i.e. sum across all subpopulations)

    // for each subpopulation
    for (nit = n.begin(); nit != n.end(); nit++)
      {
        ntot = ntot + nit->second;
      } // end of for each subpopulation

  } // end of method calc_ntot()

  void print(int x, chromosome& chr)
  { 
    float h = chr.mutation_types.find(t)->second.h; // the dominance coefficient in the reference
      // environment

    // compute total prevalence (although this could be done while printing subpopulation-specific
      // prevalences below, for better consistency with the previous output format, the total
      // prevalence is given first and hence needs to be computed first)

    int ntot = calc_ntot();

    // cout << i << " m" << t << " " << x+1 << " " << s << " " << h << " "<< n << endl;
    // print mutation id (re-assigned every generation), mutation type, physical position,
      // selection coefficient and dominance coefficient in the reference environment, and
      // total prevalence
    cout << i << " m" << t << " " << x+1 << " " << s << " " << h << " "<< n;

    // print prevalence in each subpopulation

    // for each subpopulation
    for (nit = n.begin(); nit != n.end(); nit++)
      {
        cout << " p" << nit->first << " " << nit->second;
      } // end of for each subpopulation
    cout << endl;

  } // end of method print()

  void print(ofstream& outfile, int x, chromosome& chr)
  { 
    float h = chr.mutation_types.find(t)->second.h; // the dominance coefficient in the reference
      // environment

    // compute total prevalence (although this could be done while printing subpopulation-specific
      // prevalences below, for better consistency with the previous output format, the total
      // prevalence is given first and hence needs to be computed first)

    int ntot = calc_ntot();

    // outfile << i << " m" << t << " " << x+1 << " " << s << " " << h << " "<< n << endl;
    // print mutation id (re-assigned every generation), mutation type, physical position,
      // selection coefficient and dominance coefficient in the reference environment, and
      // total prevalence
    outfile << i << " m" << t << " " << x+1 << " " << s << " " << h << " "<< n;

    // print prevalence in each subpopulation

    // for each subpopulation
    for (nit = n.begin(); nit != n.end(); nit++)
      {
        outfile << " p" << nit->first << " " << nit->second;
      } // end of for each subpopulation
    outfile << endl;

  } // end of method print()

  void print_noi(int x, chromosome& chr, int spid)

  // x is the physical position, chr a reference to the chromosome, and spid the id of the
    // subpopulation of which to print the prevalences
  { 
    float h = chr.mutation_types.find(t)->second.h; // the dominance coefficient in the reference
      // environment
      // OPTION: Change this to subpopulation-specific dominance coefficient?
      // OPTION: Change from selection coefficient in the reference environment to the one
        // modified for the respective subpopulation?

    cout << "m" << t << " " << x+1 << " " << s << " " << h << " "<< n.find(spid)->second << endl;
  } // end of method print_noi()

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
}; // end of class partial sweep


class genome : public vector<mutation>
{
  // a genome is a vector of mutations
}; // end of class genome


// top-level method

genome fixed(genome& G1, genome& G2)
{
  // return genome G consisting only of the mutations that are present in both G1 and G2

  genome G; // a genome is a vector of mutations

  vector<mutation>::iterator g1 = G1.begin();
  vector<mutation>::iterator g2 = G2.begin();

  vector<mutation>::iterator g1_max = G1.end();
  vector<mutation>::iterator g2_max = G2.end();
  
  // advance g1 and g2 while non is at their maximum
  while (g1 != g1_max && g2 != g2_max)
    {
      // advance g1 while g1.x < g2.x
      while (g1 != g1_max && g2 != g2_max && (*g1).x < (*g2).x)
        { g1++; }

      // advance g2 while g1.x < g2.x
      while (g1 != g1_max && g2 != g2_max && (*g2).x < (*g1).x)
        { g2++; }
	   
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
                  if ((*temp).t == (*g1).t && (*temp).s == (*g1).s)
                    { G.push_back(*g1); }
                  temp++;
                }
              g1++;
            }
        } // end of identify shared mutations at positions x and add to G
    } // end of advance g1 and g2

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
  
  // advance g1 and g2 while non is hitting their maximum
  while (g1 != g1_max && g2 != g2_max)
    {
      // advance g1 while g1.x < g2.x
      while (g1 != g1_max && g2 != g2_max && (*g1).x < (*g2).x)
        { G.push_back(*g1); g1++; }

      // advance g2 while g1.x < g2.x
      while (g2 != g2_max && g1 != g1_max && (*g2).x < (*g1).x)
        { g2++; }
	   
      // identify polymorphic mutations at positions x and add to G
      if (g2 != g2_max && g1 != g1_max && (*g2).x == (*g1).x)
        {
          int x = (*g1).x;

          // traverse g1 and check for those mutations that are not present in g2
          vector<mutation>::iterator temp = g2;

          while (g1 != g1_max && (*g1).x == x)
            {
              bool poly = 1;

              while (temp != g2_max && (*temp).x == x)
                {
                  if ((*g1).t==(*temp).t && (*g1).s==(*temp).s)
                    { poly = 0; }
                  temp++;
                }
              if (poly == 1)
                { G.push_back(*g1); }
              g1++;
            }

          while (g2 != g2_max && (*g2).x == x)
            { g2++; }
        } // end of identify polymorphic positions x and add to G
    } // end of advance g1 and g2

  while (g1 != g1_max)
    { G.push_back(*g1); g1++; }

  return G;

} // end of method polymorphic()


class environment
{
    // an environment reassigns a dominance coefficient and modifies the (mean) selection
    // coefficient with respect to the reference environment and for a specified set of mutation
    // types

public:

  map<int,double> h; // dominance coefficients for each mutation-type; the coefficient is by
    // default equal to the one in the reference environment for each mutation-type, but can
    // be changed for this environment by the user, for specific mutation-types; the key of h is
    // the identifier of the affected mutation-types.
  map<int,double> smodif; // modifier of selection coefficient for each mutation-type; the modifier
    // is by default equal to 1, meaning that the selection coefficient is equal to the one in the
    // reference environment. Yet, the modifier can be changed for this environment by the user,
    // for specific mutation-types; the key of smodif is the identifier of the affected mutation-
    // types.

  // default constructor
  environment(void) { ; }

  // extended constructor 1: initializes the reference environment
  environment(chromosome& chr)
    // chr is the chromosome with all mutation-types
  {

    // initialize by copying from the reference environment

    init(chr);

  } // end of extended contstructor 1

  // extended constructor 2: starts from the reference environment, then updates specified entries
  environment(chromosome& chr, vector<int> MUp, vector<double> HUp, vector<double> SMODIFUp)
    // chr is the chromosome with all the mutation-types
    // MUp is the vector of mutation-types affected by a change relative to the reference
    //    environment
    // HUp is the vector of dominance coefficients for the mutation-types affected by the change;
    //    the order of the entries must correspond to the order of the entries in MUp
    // SMODIFUp is the vector of modifiers of the selection coefficients for the mutation-types
    //    affected by the change; the order of the entries must correspond to the order of the
    //    entries in MUp
  {

    if ( HUp.size() != MUp.size() || SMODIFUp.size() != MUp.size())
      { exit(1); }

    // initialize by copying from the reference environment

    init(chr);

    // alter maps h and smodif based on updates HUp, MUp, and SMODIFUp

    // for each affected mutation type
    for (int i = 0; i < MUp.size(); i++)
      {
        h.insert(pair<int,double>(MUp[i], HUp[i]));
        smodif.insert(pair<int,double>(MUp[i], SMODIFUp[i]));
      } // end of for each affected mutation type

  } // end of extended constructor 2

private:

  // initializes the environment as (a copy of) the reference environment implicitly defined
    // by the mutation-types known for chromosome chr
  void init(chromosome& chr)
  {
    map<int,mutation_type>::iterator mti; // iterator over mutation-types in chromosome chr

    // iterate over mutation-types and retrieve dominance coefficients to initialize reference
      // environment; at the same time, set the modifier of selection to 1 for each visited
      // mutation-type (which means no modification w.r.t. the reference environment)

    for (mti = chr.mutation_types.begin(); mti != chr.mutation_types.end(); mti++)
      { // for each mutation-type in chromosome chr
        // retrieve dominance coefficient
        h.insert(pair<int,double>(mti->first, mti->second.h));
        // assign default modifier of selection
        smodif.insert(pair<int,double>(mti->first, 1.0));
      } // end of for each mutation-type in chromosome chr

  } // end of method init()

}; // end of class environment


class subpopulation
{
  // a subpopulation is described by the vector G of 2N genomes
  // individual i is constituted by the two genomes 2*i and 2*i+1

private:
  
  gsl_ran_discrete_t* LT; // a look-up table

public:

  int    N; // population size (number of diploid individuals)
  double S; // selfing fraction
  environment E; // environment

  vector<genome> G_parent; // parent population
  vector<genome> G_child; // offspring population

  map<int,double> m; // m[i]: fraction made up of migrants from subpopulation i per generation
 
  // constructor
  subpopulation(int n, const environment& e)
    // n is the population size, e a reference to the environment assigned to this subpopulation
      // (by default, the reference environment is assigned)
  {
    N = n;
    S = 0.0;
    E = e;
    G_parent.resize(2*N);
    G_child.resize(2*N);

    double A[N];

    for (int i = 0; i < N; i++)
      { A[i] = 1.0; }

    LT = gsl_ran_discrete_preproc(N, A); // to draw individuals at random, with uniform weights

  } // end of constructor

  int draw_individual()
  {
    return gsl_ran_discrete(rng, LT);
  }

  void update_fitness(chromosome& chr, char& fi)
  {
    // updating fitness given mutations on chromosome chr and fitness interaction of type fi

    // calculate fitnesses in parent population and create new lookup table
    
    gsl_ran_discrete_free(LT);
    double A[(int)(G_parent.size()/2)]; // individuals are diploid
    for (int i = 0; i < (int)(G_parent.size()/2); i++)
      {
        A[i] = W(2*i, 2*i+1, chr, fi); // recall that genomes are stored such that
          // individual i is made up of haploid genomes 2*i and and 2*i + 1
      }
    LT = gsl_ran_discrete_preproc((int)(G_parent.size()/2), A); // assigning weights to individuals
      // that correspond to their (abolute) fitness
  } // end of method update_fitness()


  double W(int i, int j, chromosome& chr, char& fi)
  {
    // calculate the fitness of the individual constituted by genomes i and j in the parent
    // population, where i and j are indices to the vector population. Consider mutations
    // on chromosome chr and assume fitness interaction type fi. Use dominance coefficients and
    // modifiers of selection coefficients with respect to the environment e of this
    // subpopulation;

  double w = 1.0;

    // a genome is a vector of mutations; recall that mutations in a genome are ordered by physical
    // position
  vector<mutation>::iterator pi = G_parent[i].begin(); // iterator pointing to the first mutation
    // in the haploid genome i of the parental population
  vector<mutation>::iterator pj = G_parent[j].begin();

  vector<mutation>::iterator pi_max = G_parent[i].end(); // iterator pointing to the end of the
    // vector with mutations in the haploid genome i of the parental population; does not point
    // to an actual element, but to the virtual past-the-end element
  vector<mutation>::iterator pj_max = G_parent[j].end();

  while (w > 0 && (pi != pi_max || pj != pj_max)) // while fitness is strictly positive, and
      // at least one parental genome has mutations not jet visited
    {

      // deal with all mutation except those residing at the same position x in the two parental
      // genomes (i.e., exclude cases where distinction must be made between homo- and hetero
      // zygotes)

      // advance i while pi.x < pj.x (x is the physical position)
      while (pi != pi_max && (pj == pj_max || (*pi).x < (*pj).x)) // while there are unvisited
          // mutations left in genome i and (there are no unvisited mutations left in genome j or
          // the position of the current mutation in genome i is smaller than that of the current
          //  mutation in genome j)
        {
          // neutrality is assessed w.r.t. the environment assigned to this subpopulation
          if (E.smodif.find((*pi).t)->second * (*pi).s != 0) // if this mutation
            // in genome i is not neutral in environment with id EI
            {
              // distinguish between two fitness regimes (additive vs. multiplicative)
              if (fi == 'a') // if fitness interaction is additive
                {
                  // DONE: Changed to environment-specific fitnesses
                  w = w + E.h.find((*pi).t)->second * E.smodif.find((*pi).t)->second * (*pi).s; // the dominance coefficient is
                    // chosen according to the environment EI assigned to this subpopulation; the
                    // selection coefficient is modified according to the environment EI assigned
                    // to this subpopulation
                  /* OLD:
                  w = w + chr.mutation_types.find((*pi).t)->second.h * (*pi).s; // the dominance
                  // coefficient is a public variable of the class mutation type; the selection
                  // coefficient is a public variable of the class mutation
                  */
                } // fitness interaction is not additive
              else if (fi == 'm') // ie fitness interaction is multiplicative
                {
                  // DONE: Changed to environment-specific fitnesses
                  w = w * (1.0 + E.h.find((*pi).t)->second * E.smodif.find((*pi).t)->second * (*pi).s);
                  /* OLD
                  w = w * (1.0 + chr.mutation_types.find((*pi).t)->second.h * (*pi).s);
                  */
                } // fitness interaction is nEIther additive nor multiplicative
              else
                { exit(1); }
            } // end of if this mutation in i is not neutral
          pi++;
        } // end of advancing i while pi.x < pj.x; pi == pi_max || (*pi).x >= (*pj).x

      // advance j while pj.x < pi.x (x is the physical position)
      while (pj != pj_max && (pi == pi_max || (*pj).x < (*pi).x)) // while there are unvisited
        // mutations left in genome j and (there are no unvisited mutations left in genome i or
        // the position of the current mutation in genome j is smaller than that of the current
        //  mutation in genome i)
        {
          // neutrality is assessed w.r.t. the environment assigned to this subpopulation
          if (E.smodif.find((*pj).t)->second * (*pj).s != 0)// if this mutation
            // in genome j is not neutral in environment EI
            {
              // distinguish between two fitness regimes (additive vs. multiplicative)
              if (fi == 'a') // if fitness interaction is additive
                {
                  // DONE: Changee to environment-specific fitnesses
                  w = w + E.h.find((*pj).t)->second * E.smodif.find((*pj).t)->second * (*pj).s;
                  /* OLD:
                  w = w + chr.mutation_types.find((*pj).t)->second.h * (*pj).s; // see above for
                    // details
                  */
                } // fitness interaction is not additive
              else if (fi == 'm')
                {
                  // DONE: Changed to environment-specific fitnesses
                  w = w * (1.0 + E.h.find((*pj).t)->second * E.smodif.find((*pj).t)->second * (*pj).s);
                  /* OLD:
                  w = w * (1.0 + chr.mutation_types.find((*pj).t)->second.h * (*pj).s);
                  */
                } // fitness interaction is neither additive nor multiplicative
              else
                { exit(1); }
            } // end of if this mutation in j is not neutral
          pj++;
        } // end of advancing j while pj.x < pi.x; pj == pj_max || (*pj).x >= (*pi).x

        // check for homozygotes and heterozygotes at position x, where homozygosity is determined
          // as identity in state in the environment assigned to this subpopulation; note that
          // this is identical to identity by descent here, because we allow for at most one non-
          // neutral mutation at any given position in the entire population

      if (pi != pi_max && pj != pj_max && (*pj).x == (*pi).x) // if there are unvisited mutations
        // left in both genomes i and j, and if the currently visited mutations are at the same
        // physical position x
        {
          int x = (*pi).x;
          int ns_i = 0; // number of mutations in genome i encountered at position x that are
            // non-neutral in this subpopulation (must be <= 1)
          int ns_j = 0; // number of mutations in genome j encountered at position x that are
            // non-neutral in this subpopulation (must be <= 1)

          vector<mutation>::iterator pi_start = pi;

          // advance through pi, i.e. through mutations in genome i that reside at the same
            // position; recall that one haploid genome can harbour more than one neutral mutation
            // at the same position; however, no more than one non-neutral mutation at position x
            // is allowed in the entire population; identity of mutations is measured in terms
            // of mutation type and selection coefficient in the reference environment
          while (pi != pi_max && (*pi).x == x) // visiting all mutations in genome i that reside
              // at a given physical position x
            {
              // neutrality is assessed w.r.t. the environment assigned to this subpopulation
              if (E.smodif.find((*pi).t)->second * (*pi).s != 0) // if this
                // mutation in genome i is not neutral in environment EI
                {
                  vector<mutation>::iterator temp_j = pj;
                  bool homo = 0;

                  while (homo == 0 && temp_j != pj_max && (*temp_j).x == x) // advance through
                    // mutations in genome j at position x, as long as no homozygote is found
                    {
                      // Note that no permanent mutation identifiers are stored for reasons of
                        // computational efficiency (memory). Hence identity is assessed based on
                        // mutation type and selection coefficient
                      // DONE: Stuck to the practice of not storing mutation identifiers, but
                        // now identity must be defined based on properties (mutation type and
                        // selection coefficient) in the *reference environment*.
                      if ((*pi).t == (*temp_j).t && E.smodif.find((*pi).t)->second * (*pi).s == E.smodif.find((*temp_j).t)->second * (*temp_j).s) // if currently
                        // visited mutations in genomes i and j are biologically the same
                        // (identical by state in environment EI)
                        {
                          // distinguish between two fitness regimes (additive vs.
                            // multiplicative)
                          if (fi == 'a')
                            {
                              // DONE: Changed to environment-specific fitnesses
                              w = w + E.smodif.find((*pi).t)->second * (*pi).s;
                              /* OLD:
                              w = w + (*pi).s;
                              */
                              homo = 1;
                            } // fitness interaction is not additive
                          else if (fi == 'm')
                            {
                              w = w * (1.0 + E.smodif.find((*pi).t)->second * (*pi).s);
                              /* OLD:
                              w = w * (1.0 + (*pi).s);
                              */
                              homo = 1;
                            } // fitness interaction is neither additive nor multiplicative
                          else
                            { exit(1); }
                        } // currently visited mutations are not the same; the parent is not
                            // homozygous
                      temp_j++;
                    } // end of advance through mutations in genome j at position x, as long as
                        // no homozygote is found
                  if (homo == 0) // if the parent is not homozygous at position x
                    // Note that the original implementation (SLiM vs. 1.x) was conceptually
                      // flawed. Dominance was applied only from the perspective of a focal
                      // mutation (*pi), irrespective of the other mutation present. However,
                      // dominance is a feature of the genotype.
                    // New: Only one mutation of a non-neutral type is allowed at any given
                      // physical position (neutrality being defined w.r.t. the reference
                      // environment). It then has to be the case that the mutation at x in gemone
                      // j must be a neutral one, given that the mutation at x in genome i is
                      // non-neutral and the individual is heterozygous.
                    {
                      // distinguish between two fitness regimes (additive vs. multiplicative)
                      if (fi == 'a')
                        {
                          // DONE: Changed to environment-specific fitnesses
                          w = w + E.h.find((*pi).t)->second * E.smodif.find((*pi).t)->second * (*pi).s;
                          /* OLD:
                          w = w + chr.mutation_types.find((*pi).t)->second.h * (*pi).s;
                          */
                        }
                      else if (fi == 'm')
                        {
                          // DONE: Changed to environment-specific fitnesses
                          w = w * (1.0 + E.h.find((*pi).t)->second * E.smodif.find((*pi).t)->second * (*pi).s);
                          /* OLD:
                          w = w * (1.0 + chr.mutation_types.find((*pi).t)->second.h * (*pi).s);
                          */
                        }
                      else
                        { exit(1); }
                    } // end of if parent is not homozygous
                  ns_i++;
                } // end of if this mutation in genome i is not neutral; hence, it is neutral
              pi++;
            } // end of visiting mutations in genome i that reside at the same physical pos. x

          // advance through pj, i.e. through mutations in genome j that reside at the same
            // position; recall that one haploid genome can harbour more than one neutral mutation
            // at the same position; however, no more than one non-neutral mutation at position x
            // is allowed in the entire population
          while (pj != pj_max && (*pj).x == x) // visiting all mutations in genome j that reside
              // at a given physical position x; note that x is the position of the focal mutation
              // currently visited in the outer loop in genome i, but we are here in the
              // situation where (*pj).x == (*pi).x holds.
            {
              // neutrality is assessed w.r.t. the environment assigned to this subpopulation
              if (E.smodif.find((*pj).t)->second * (*pj).s != 0.0) // if this mutation in genome j is not neutral in environment EI
                {
                  vector<mutation>::iterator temp_i = pi_start; // recall that pi_start was
                    // assigned pi above, i.e. pi_start = pi, but in the meantime, pi has been
                    // increased
                  bool homo = 0;

                  while (homo == 0 && temp_i != pi_max && (*temp_i).x == x) // advance through
                      // mutations in genome i at position x, as long as no homozygote is found
                    {
                      // recall that homozygosity is defined as identity by state in the
                        // environment assigned to this subpopulation
                      if ((*pj).t == (*temp_i).t && E.smodif.find((*pj).t)->second * (*pj).s == E.smodif.find((*temp_i).t)->second * (*temp_i).s) // if currently visited mutations in genomes i and j are
                        // biologically the same (identical by state) in environment EI
                        {
                          // this genotype was encountered before and its contribution to fitness
                          // has been incorporated
                          homo = 1;
                        }
                      temp_i++;
                    } // end of advance through mutations in genome i at position x, as long as
                        // no homozygote is found
                  if (homo == 0) // if the parent is not homozygous at position x
                    // Note that this implementation (the original one in SLiM) is conceptually
                      // flawed. Dominance is applied only from the perspective of a focal
                      // mutation (*pj), irrespective of the other mutation present. However,
                      // dominance is a feature of the genotype.
                    // Only one mutation of a non-neutral type is allowed at any given physical
                      // position (neutrality being defined w.r.t. the reference environment). It
                      // then has to be the case that the mutation at position x in gemone
                      // j must be a neutral one, given that the mutation at x in genome i is
                      // non-neutral and the individual is heterozygous.
                    {
                      // distinguish between two fitness regimes (additive vs. multiplicative)
                      if (fi == 'a')
                        {
                          // DONE: Changed to environment-specific fitnesses
                          w = w + E.h.find((*pj).t)->second * E.smodif.find((*pj).t)->second * (*pj).s;
                          /* OLD:
                          w = w + chr.mutation_types.find((*pj).t)->second.h * (*pj).s;
                          */
                        }
                      else if (fi == 'm')
                        {
                          // DONE: Changed to environment-specific fitnesses
                          w = w * (1.0 + E.h.find((*pj).t)->second * E.smodif.find((*pj).t)->second * (*pj).s);
                          /* OLD:
                          w = w * (1.0 + chr.mutation_types.find((*pj).t)->second.h * (*pj).s);
                          */
                        }
                      else
                        { exit(1); }
                    } // end of if parent is not homozygous
                  ns_j++;
                } // end of if this mutation in genome j is not neutral; hence, it is neutral
              pj++;
            } // end of visiting mutations in genome j that reside at the same physical pos. x

          if (ns_i > 1 || ns_j > 1) // encountered more than one non-neutral mutation at pos. x
            // in at least one parental haplotype
            {
              if (ns_i > 1)
                {
                  cerr << "ERROR (compute fitness): encountered " << ns_i << " non-neutral mutations at position " << x << " in parental genome " << i << "." << endl;
                }
              if (ns_j > 1)
                {
                  cerr << "ERROR (compute fitness): encountered " << ns_j << " non-neutral mutations at position " << x << " in parental genome " << j << "." << endl;
                }
              exit(1);
            }
        } // end of if there are unvisited mutations left in both genomes i and j, and if the
            // currently visited mutations are at the same physical position x; there are either
            // no unvisited mutations left in at least one genome, or physical positions are not
            // the same
    } // end of while fitness is strictly positive and at least one parental genome has mutations
        // not yet visited; hence, fitness is either not strictly positive or at least one of the
        // parental genomes has no more mutations left to be visited.

  if (w < 0)
    { w = 0.0; }

  return w;

  } // end of method W()


  double oldW(int i, int j, chromosome& chr, char& fi)
  {
    // calculate the fitness of the individual constituted by genomes i and j in the parent
      // population, where i and j are indices to the vector population. Consider mutations
      // on chromosome chr and assume fitness interaction type fi

    double w = 1.0;

    // a genome is a vector of mutations; recall that mutations in a genome are ordered by physical
      // position
    vector<mutation>::iterator pi = G_parent[i].begin(); // iterator pointing to the first mutation
      // in the haploid genome i of the parental population
    vector<mutation>::iterator pj = G_parent[j].begin();

    vector<mutation>::iterator pi_max = G_parent[i].end(); // iterator pointing to the end of the
      // vector with mutations in the haploid genome i of the parental population; does not point
      // to an actual element, but to the virtual past-the-end element
    vector<mutation>::iterator pj_max = G_parent[j].end();

    while (w > 0 && (pi != pi_max || pj != pj_max)) // while fitness is strictly positive, and
      // at least one parental genome has mutations not jet visited
      {

        // deal with all mutation except those residing at the same position x in the two parental
          // genomes (i.e., exclude cases where distinction must be made between homo- and hetero
          // zygotes)

        // advance i while pi.x < pj.x (x is the physical position)
        while (pi != pi_max && (pj == pj_max || (*pi).x < (*pj).x)) // while there are unvisited
          // mutations left in genome i and (there are no unvisited mutations left in genome j or
          // the position of the current mutation in genome i is smaller than that of the current
          //  mutation in genome j)
          {
            if ((*pi).s != 0) // if this mutation in genome i is not neutral
              {
                // distinguish between two fitness regimes (additive vs. multiplicative)
                if (fi == 'a') // if fitness interaction is additive
                  {
                    // TODO: Change to environment-specific fitnesses
                    w = w + chr.mutation_types.find((*pi).t)->second.h * (*pi).s; // the dominance
                      // coefficient is a public variable of the class mutation type; the selection
                      // coefficient is a public variable of the class mutation
                  } // fitness interaction is not additive
                else if (fi == 'm') // ie fitness interaction is multiplicative
                  {
                    // TODO: Change to environment-specific fitnesses
                    w = w * (1.0 + chr.mutation_types.find((*pi).t)->second.h * (*pi).s);
                  } // fitness interaction is neither additive nor multiplicative
                else
                  { exit(1); }
              } // end of if this mutation in i is not neutral
            pi++;
          } // end of advancing i while pi.x < pj.x; pi == pi_max || (*pi).x >= (*pj).x
	   
        // advance j while pj.x < pi.x (x is the physical position)
        while (pj != pj_max && (pi == pi_max || (*pj).x < (*pi).x)) // while there are unvisited
          // mutations left in genome j and (there are no unvisited mutations left in genome i or
          // the position of the current mutation in genome j is smaller than that of the current
          //  mutation in genome i)
          {
            if ((*pj).s != 0) // if this mutation in genome j is not neutral
              {
                // distinguish between two fitness regimes (additive vs. multiplicative)
                if (fi == 'a') // if fitness interaction is additive
                  {
                    // TODO: Change to environment-specific fitnesses
                    w = w + chr.mutation_types.find((*pj).t)->second.h * (*pj).s; // see above for
                      // details
                  } // fitness interaction is not additive
                else if (fi == 'm')
                  {
                    // TODO: Change to environment-specific fitnesses
                    w = w * (1.0 + chr.mutation_types.find((*pj).t)->second.h * (*pj).s);
                  } // fitness interaction is neither additive nor multiplicative
                else
                  { exit(1); }
              } // end of if this mutation in j is not neutral
            pj++;
          } // end of advancing j while pj.x < pi.x; pj == pj_max || (*pj).x >= (*pi).x
	
        // check for homozygotes and heterozygotes at position x

        if (pi != pi_max && pj != pj_max && (*pj).x == (*pi).x) // if there are unvisited mutations
          // left in both genomes i and j, and if the currently visited mutations are at the same
            // physical position x
          {
            int x = (*pi).x;
	   
            vector<mutation>::iterator pi_start = pi;

            // advance through pi, i.e. through mutations in genome i that reside at the same
              // position; recall that one (haploid!) 'genome' can harbour more than one mutation
              // at the same position, as long as the mutations do not have the same type AND
              // selection coefficient.
            // TODO: Change this. Allow for only two segregating mutations of a non-neutral type at
              // a given physical position in the whole population. Probably best to handle this
              // case when mutations are created, but then do a check here, and exit if the rule is
              // broken.

            while (pi != pi_max && (*pi).x == x) // visiting all mutations in genome i that reside
              // at a given physical position x
              {
                if ((*pi).s != 0.0) // if this mutation in genome i is not neutral
                  {
                    vector<mutation>::iterator temp_j = pj;
                    bool homo = 0;

                    while (homo == 0 && temp_j != pj_max && (*temp_j).x == x) // advance through
                      // mutations in genome j at position x, as long as no homozygote is found
                      {
                        // Note that no permanent mutation identifiers are stored for reasons of
                          // computational efficiency (memory). hence identity is assessed based on
                          // mutation type and selection coefficient
                        // TODO: Stick to the practice of not storing mutation identifiers, but
                          // then identity must be defined based on properties (mutation type and
                          // selection coefficient) in the reference environment.
                        if ((*pi).t == (*temp_j).t && (*pi).s == (*temp_j).s) // if currently
                          // visited mutations in genomes i and j are biologically the same
                            // (identical by state)
                          {
                            // distinguish between two fitness regimes (additive vs.
                              // multiplicative)
                            if (fi == 'a')
                              {
                                // TODO: Change to environment-specific fitnesses
                                w = w + (*pi).s;
                                homo = 1;
                              } // fitness interaction is not additive
                            else if (fi == 'm')
                              {
                                w = w * (1.0 + (*pi).s);
                                homo = 1;
                              } // fitness interaction is neither additive nor multiplicative
                            else
                              { exit(1); }
                          } // currently visited mutations are not the same; the parent is not
                              // homozygous
                        temp_j++;
                      } // end of advance through mutations in genome j at position x, as long as
                          // no homozygote is found
                    if (homo == 0) // if the parent is not homozygous at position x
                      // Note that this implementation (the original one in SLiM) is conceptually
                        // flawed. Dominance is applied here only from the perspective of a focal
                        // mutation (*pi), irrespective of the other mutation present. However,
                        // dominance is a feature of the genotype.
                      // TODO: Fix this by allowing only one mutation of a non-neutral type at any
                        // physical position (neutrality being defined w.r.t. the reference
                        // environment). It then has to be the case that the mutation at x in
                        // gemone j must be a neutral one, given that the mutation at x in genome i
                        // is non-neutral and the individual is heterozygous.
                      {
                        // distinguish between two fitness regimes (additive vs. multiplicative)
                        if (fi == 'a')
                          {
                            // TODO: Change to environment-specific fitnesses
                            w = w + chr.mutation_types.find((*pi).t)->second.h * (*pi).s;
                          }
                        else if (fi == 'm')
                          {
                            // TODO: Change to environment-specific fitnesses
                            w = w * (1.0 + chr.mutation_types.find((*pi).t)->second.h * (*pi).s);
                          }
                        else
                          { exit(1); }
                      } // enf of if parent is not homozygous
                  } // end of if this mutation in genome i is not neutral; hence, it is neutral
                pi++;
              } // end of visiting mutations in genome i that reside at the same physical pos. x

            // advance through pj, i.e. through mutations in genome j that reside at the same
              // position; recall that one (haploid!) 'genome' can harbour more than one mutation
              // at the same position, as long as the mutations do not have the same type AND
              // selection coefficient.
           while (pj != pj_max && (*pj).x == x) // visiting all mutations in genome j that reside
              // at a given physical position x; note that x is the position of the focal mutation
                // currently visited in the outer loop in genome i, but we are here in the
                // situation where (*pj).x == (*pi).x holds.
              {
                if ((*pj).s != 0.0) // if this mutation in genome j is not neutral
                  {
                    vector<mutation>::iterator temp_i = pi_start; // recall that pi_start was
                      // assigned pi above, i.e. pi_start = pi, but in the meantime, pi has been
                      // increased
                    bool homo = 0;

                    while (homo == 0 && temp_i != pi_max && (*temp_i).x == x) // advance through
                      // mutations in genome i at position x, as long as no homozygote is found
                      {
                        // homozygosity is defined as identity by state in the environment assigned
                          // to this subpopulation
                        if ((*pj).t == (*temp_i).t && (*pj).s == (*temp_i).s) // if currently
                          // visited mutations in genomes i and j are biologically the same
                            // (identical by state)
                          {
                            // we encountered this genotype before and have incorporated its
                              // contribution to fitness
                            homo = 1;
                          }
                        temp_i++;
                      } // end of advance through mutations in genome i at position x, as long as
                          // no homozygote is found
                    if (homo == 0) // if the parent is not homozygous at position x
                      // Note that this implementation (the original one in SLiM) is conceptually
                        // flawed. Dominance is applied only from the perspective of a focal
                        // mutation (*pj), irrespective of the other mutation present. However,
                        // dominance is a feature of the genotype.
                      {
                        // distinguish between two fitness regimes (additive vs. multiplicative)
                        if (fi == 'a')
                          {
                            w = w + chr.mutation_types.find((*pj).t)->second.h * (*pj).s;
                          }
                        else if (fi == 'm')
                          {
                            w = w * (1.0 + chr.mutation_types.find((*pj).t)->second.h * (*pj).s);
                          }
                        else
                          { exit(1); }
                      } // end of if parent is not homozygous
                  } // end of if this mutation in genome j is not neutral; hence, it is neutral
                pj++;
              } // end of visiting mutations in genome j that reside at the same physical pos. x
          } // end of if there are unvisited mutations left in both genomes i and j, and if the
              // currently visited mutations are at the same physical position x; there are either
              // no unvisited mutations left in at least one genome, or physical positions are not
              // the same
      } // end of while fitness is strictly positive and at least one parental genome has mutations
          // not yet visited; hence, fitness is either not strictly positive or at least one of the
          // parental genomes has no more mutations left to be visited.

    if (w < 0)
      { w = 0.0; }

    return w;

  } // end of method oldW()

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

  void add_subpopulation(int i, unsigned int N, const environment& env)
  {
    // add new empty subpopulation i of size N (i is the key of the subpopulation) and assign
      // environment env to it (by default, the reference environment should be assigned)

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

    insert(pair<int, subpopulation>(i, subpopulation(N, env)));

  } // end of add_subpopulation() method

  void add_subpopulation(int i, int j, unsigned int N, const environment& env)
  { 
    // add new subpopulation i of size N individuals drawn from source subpopulation j and assign
      // environment env to it (by default, the reference environment should be assigned)

    if (count(i) != 0)
      {
        cerr << "ERROR (add subpopulation): subpopulation p"<< i << " already exists" << endl;
        exit(1);
      }
    if (count(j) == 0)
      {
        cerr << "ERROR (add subpopulation): source subpopulation p"<< j << " does not exists" << endl;
        exit(1);
      }
    if (N < 1)
      {
        cerr << "ERROR (add subpopulation): subpopulation p"<< i << " empty" << endl;
        exit(1);
      }

    insert(pair<int,subpopulation>(i,subpopulation(N, env)));

    for (int p = 0; p < find(i)->second.N; p++)
      {
        // draw individual from subpopulation j and assign to be a parent in i

        int m = find(j)->second.draw_individual(); // returns an id of an individual in
          // subpopulation j at random
	
        find(i)->second.G_parent[2*p] = find(j)->second.G_parent[2*m];
        find(i)->second.G_parent[2*p+1] = find(j)->second.G_parent[2*m+1];
      }
  } // end of add_subpopulation() [alt] method


  void set_size(int i, unsigned int N) 
  {
    // set size of subpopulation i to N

    if (count(i) == 0)
      {
        cerr << "ERROR (change size): no subpopulation p"<< i << endl;
        exit(1);
      }

    if (N == 0) // remove subpopulation i
      {
        erase(i);

        // visit all subpopulations and erase the migration rate from the removed subpopulation i
          // to the current subpopulation
        for (it = begin(); it != end(); it++)
          { it->second.m.erase(i); }
      }
    else
      {
        find(i)->second.N = N;
        find(i)->second.G_child.resize(2*N);
      }
  } // end of method set_size()

  void set_selfing(int i, double s) 
  { 
    // set fraction s of i that reproduces by selfing
 
    if (count(i) == 0)
      {
        cerr << "ERROR (set selfing): no subpopulation p"<< i << endl;
        exit(1);
      }
    if (s < 0.0 || s > 1.0)
      {
        cerr << "ERROR (set selfing): selfing fraction has to be within [0,1]" << endl;
        exit(1);
      }
    find(i)->second.S = s;
  } // end of method set_selfing()


  void set_migration(int i, int j, double m) 
  { 
    // set fraction m of i that originates as migrants from j per generation  

    if (count(i) == 0)
      {
        cerr << "ERROR (set migration): no subpopulation p"<< i << endl;
        exit(1);
      }
    if (count(j) == 0)
      {
        cerr << "ERROR (set migration): no subpopulation p"<< j << endl;
        exit(1);
      }
    if (m < 0.0 || m > 1.0)
      {
        cerr << "ERROR (set migration): migration fraction has to be within [0,1]" << endl;
        exit(1);
      }

    // erase previously existing entries for migration fraction from j to i
    if (find(i)->second.m.count(j) != 0) // m is a map
      {
        find(i)->second.m.erase(j);
      }

      find(i)->second.m.insert(pair<int,double>(j, m));
  } // end of method set_migration()


  void assign_environment(int i, int e, const map<int, environment>& envs)

    // assign environment e to subpopulation i
  {

    if (count(i) == 0)
      {
        cerr << "ERROR (assign environment): no subpopulation p"<< i << endl;
        exit(1);
      }

    if (envs.count(e) == 0)
      {
        cerr << "ERROR (assign environment): no environment e"<< e << endl;
        exit(1);
      }

    find(i)->second.E = envs.find(e)->second;

  } // end of method assign_environment()

  /*
   Executes event E in generation g for chromosome chr. The vector FM stores the id's of mutation
   types to be tracked. The reference environment renv and the map of user-specified environments
   env are taken as further arguments.
   */
  void execute_event(event& E, int g, chromosome& chr, vector<int>& FM, const environment& renv, map<int, environment>& envs)
  {
    char type = E.t;

    if (type == 'P') // add subpopulation
      { 
        if (E.np == 2) // empty subpopulation (i.e. not drawing from another existing
          // subpopulation); np is the number of event-specific parameters
          {
            string sub = E.s[0]; // recall: s of class event is a vector of strings
            sub.erase(0, 1); // erasing leading "p" for subpopulation

            int i = atoi(sub.c_str()); // subpopulation identifier
            int n = (int)atof(E.s[1].c_str()); // initial size of new subpopulation
            add_subpopulation(i, n, renv);
          } // end of if empty subpopulation
	      
        if (E.np == 3) // drawn from source population
          {
            string sub1 = E.s[0]; sub1.erase(0,1);
            string sub2 = E.s[2]; sub2.erase(0,1);

            int i = atoi(sub1.c_str());
            int j = atoi(sub2.c_str());
            int n = (int)atof(E.s[1].c_str());
            add_subpopulation(i, j, n, renv); // note that the reference environment is assigned
              // to the new subpopulation even though founder individuals are drawn from
              // subpopulation j, which may have been assigned an environment other than the
              // reference environment
          } // end of if draw from source population
      } // end of if add subpopultion
	  
    if (type == 'N') // set subpopulation size
      { 
        string sub = E.s[0];
        sub.erase(0, 1); // erasing leading "p" for subpopulation

        int i = atoi(sub.c_str());
        int n = (int)atof(E.s[1].c_str());

        set_size(i, n);
      } // end of set subpopulation size

    if (type == 'S') // set selfing rate
      { 
        string sub = E.s[0];
        sub.erase(0, 1); // erasing leading "p" for subpopulation

        int i    = atoi(sub.c_str());
        double s = atof(E.s[1].c_str());

        set_selfing(i,s);
      } // end of set selfing rate
	  
    if (type == 'M') // change migration rate
      {
        string sub1 = E.s[0]; sub1.erase(0, 1);
        string sub2 = E.s[1]; sub2.erase(0, 1);

        int    i = atoi(sub1.c_str());
        int    j = atoi(sub2.c_str());
        double m = atof(E.s[2].c_str()); // the proportion of subpopulation i made up by
          // individuals from subpopulation j every generation (m_ij)

        set_migration(i, j, m);
      } // end of set migration rate

    if (type == 'A') // output state of entire population
      {
        if (E.s.size() == 0) // if no filename given (i.e. if printing to screen)
          {
            cout << "#OUT: " << g << " A" << endl;
            print_all(chr);
          } // end of if printing to screen
        if (E.s.size() == 1) // if filename given (i.e. if writing to file)
          {
            ofstream outfile;
            outfile.open (E.s[0].c_str()); // open outstream to specified output file

            // start by printing all parameter settings
            for (int i = 0; i<parameters.size(); i++)
              { outfile << parameters[i] << endl; }

            // if outstream is open, write polymorphsims and close outstream
            if (outfile.is_open())
              {
                outfile << "#OUT: " << g << " A " << E.s[0].c_str() << endl;
                print_all(outfile, chr);
                outfile.close();
              }
            else // outstream not open
              {
                cerr << "ERROR (output): could not open "<< E.s[0].c_str() << endl;
                exit(1);
              }
          } // end of if filename given
      } // end of output state of entire population

    if (type == 'R') // output random subpopulation sample
      {
        string sub = E.s[0];
        sub.erase(0, 1); // remove leading "p" of subpopulation id

        int i = atoi(sub.c_str()); // subpopulation id
        int n = atoi(E.s[1].c_str()); // sample size (number of haploid genomes)

        cout << "#OUT: " << g << " R p" << i << " " << n << endl;

        if (E.s.size() == 3 && E.s[2] == "MS") // if ms format requested
          {
            print_sample_ms(i, n, chr);
          }
        else
          {
            print_sample(i, n, chr);
          }
      } // end of output random subpopulation sample

    if (type == 'F') // output list of fixed mutations
      {
        cout << "#OUT: " << g << " F " << endl;
        cout << "Mutations:" << endl;
        for (int i = 0; i < Substitutions.size(); i++) // for each substitution
          {
            cout << i+1;
            // print mutation type, position, selection coefficient and dominance coefficient (both
              // in the reference environment), and generation at which fixed
            Substitutions[i].print(chr);
          } // end of for each substitution
      } // end of output list of fixed mutations

    if (type == 'T') // track mutation-types
      {
        string sub = E.s[0];
        sub.erase(0, 1);
        FM.push_back(atoi(sub.c_str()));
      } // end of track mutation-types

    if (type == 'E') // assign environment to subpopulation
      {
        // event has two parameters, the subpopulation (e.g. p1) and an environment (e.g. e1) to
          // be assigned to that subpopulation
        string sub1 = E.s[0]; sub1.erase(0, 1); // numeric id of subpopulation
        string sub2 = E.s[1]; sub2.erase(0, 1); // numeric id of environment

        int p = atoi(sub1.c_str());
        int e = atoi(sub2.c_str());

        assign_environment(p, e, envs);

      } // end of if assign environment

  } // end of execute_event() method


  void introduce_mutation(introduced_mutation M, chromosome& chr) 
  {
    // introduce user-defined mutation, and handle potential conflicts arising from the fact
      // that no more than one non-neutral mutation is allowed to segregate at any given position
      // in the entire population, and that only one mutation may be present at any given position
      // in any haploid genome

    // M    the mutation to be introduced (class introduced_mutation inherits from mutation)
    // chr  the chromosome

    // some tests

    if (count(M.i) == 0)
      {
        cerr << "ERROR (predetermined mutation): subpopulation "<< M.i << " does not exists" << endl;
        exit(1); }
    if (chr.mutation_types.count(M.t) == 0) 
      {
        cerr << "ERROR (predetermined mutation): mutation type m"<< M.t << " has not been defined" << endl;
        exit(1);
      }
    if (find(M.i)->second.G_child.size()/2 < M.nAA + M.nAa) 
      {
        cerr << "ERROR (predetermined mutation): not enough individuals in subpopulation "<< M.i << endl; exit(1);
      }

    mutation m;

    m.x = M.x; // position on genome
    m.t = M.t; // mutation type (to be precise, its key in the map mutation_types)

    // draw selection coefficient for novel mutation
    m.s = chr.mutation_types.find(M.t)->second.draw_s();

    // OLD:
    // subpopulation* sp = &find(M.i)->second; // sp is a pointer to the subpopulation in which the
      // new mutation should be introduced
    
    // test: print all children
    // print_all(chr);

    int n = find(M.i)->second.G_child.size()/2; // the number of diploid childred in subpopulation
      // i
    int ix_all[n], ix_hom[M.nAA], ix_het[M.nAa]; // indices of all children, and of those to become
      // homozygous and heterozygous carriers, respectively

    // assign current index to each child
    for (int i = 0; i < n; i++)
      {
        ix_all[i] = i;
      }

    if (M.l == 'e') // if successive introduced mutations should be in linkage equilibrium
      {

        // shuffle the individuals' indices
        std::random_suffle(std::begin(ix_all), std:end(ix_all));

        // OLD:
        // shuffle genomes in subpopulation (*sp) to generate linkage equilibrium
        // std::random_shuffle((*sp).G_child.begin(), (*sp).G_child.end());
        // std::random_shuffle(find(M.i)->second.G_child.begin(), find(M.i)->second.G_child.end());

      }
    else if (M.l != 'd')
      {
        cerr << "ERROR (predetermined mutation): linkage flag is "<< M.l << "instead of either e for linkage equilibrium or d for linkage disequilibrium among introduced mutations" << endl;
        exit(1);
      }

    // assign first M.nAA indices to homozygotes, and successive M.nAa to heterozygotes
    for (int i = 0; i < M.nAA; i++)
      {
        ix_hom[i] = ix_all[i];
      }
    for (int i = 0; i < M.nAa; i++)
      {
        ix_het[i] = ix_all[M.nAA + i];
      }

    // test: print all children
    // print_all(chr);

    // test against chr.seg_nonneutr_mut and clear entire population from preexisting non-
      // neutral mutation(s) at position m.x, unless the existing non-neutral mutation and the
      // one to be introduced are identical in state (in terms of mutation-type and selection
      // coefficient in the reference environment); if the latter is the case, the desired
      // genotypes are introduced; the implementation below ensures that a given haploid genome has
      // only one mutation at any given position x; previously existing mutations at x are
      // overridden;

    // find potential segregating non-neutral mutation at position x and point to it
    map<int,vector<double>>::iterator snnm_it = chr.seg_nonneutr_mut.find(m.x);

    // if a non-neutral mutation is segregating at position m.x
    if (snnm_it != chr.seg_nonneutr_mut.end())
      {
        // if the preexisting segregating non-neutral mutation is not identical in state with
          // the proposed mutation M (w.r.t. mutation-type and selection coefficient)
        if (snnm_it->second[0] != m.t || snnm_it->second[1] != m.s)
          {
            // clear the entire population from the preexisting non-neutral mutation at
              // position m.x
            remove_seg_mut(m.x);
          } // no preexisting segregating non-neutral mutation that is not also identical in
              // state with the mutation to be introduced is segregating in the entire population
      } // no preexisting non-neutral mutation other than one identical in state with m is
          // segregating in the entire population at position m.x

    vector<mutation>::iterator m_it;
    int found;

    // introduce homozygotes

    for (int j = 0; j < M.nAA; j++) // for each homozygote mutant individual to be created
      {
        // recall: a genonme is a vector of mutations

        // find subpopulation, then return two of its child genomes
        genome* g1 = &find(M.i)->second.G_child[2*ix_hom[j]]; // recall: a diploid individual k is
          // made up of haploid genomes 2*k and 2*k + 1; returns a pointer to child genome 2*k in
          // subpopulation M.i; here, k = ix_hom[j]
        genome* g2 = &find(M.i)->second.G_child[2*ix_hom[j] + 1]; // returns a pointer to child
          // genome 2*k + 1 in subpopulation M.i; here, k = ix_hom[j]

        // introduce the desired homozygote, but ensure that previously existing mutations
          // (neutral ones or those identical in state with the one to be introduced) at position
          // m.x are removed from both haploid genomes

        // find mutations present at position m.x in genomes g1 and g2 and remove them

        // advance over mutations in genome g1
        m_it = (*g1).begin();
        found = 0;

        while (m_it != (*g1).end() && found == 0) // in the new mutation scheme, there may be at
          // most one mutation at any given position in a haploid genome; once one mutation at
          // position m.x has been found, all have been found
          {
            if ((*m_it).x == m.x)
              {
                // remove currently visited mutation and point iterator to next following mutation
                  // (in principle, the pointing-to-the-next mutation is not necessary)
                m_it = (*g1).erase(m_it);
                found = 1;
              }
            else // currently visited mutation is not located at position m.x
              {
                // increase iterator
                m_it++;
              }
          } // end of advance over mutations in genome g1; all mutations in g1 visited or found the
              // one at position m.x

        // advance over mutations in genome g2
        m_it = (*g2).begin();
        found = 0;

        while (m_it != (*g2).end() && found == 0) // in the new mutation scheme, there may be at
          // most one mutation at any given position in a haploid genome; once one mutation at
          // position m.x has been found, all have been found
          {
            if ((*m_it).x == m.x)
              {
                // remove currently visited mutation and point iterator to next following mutation
                m_it = (*g2).erase(m_it);
                found = 1;
              }
            else // currently visited mutation is not located at position m.x
              {
                // increase iterator
                m_it++;
              }
          } // end of advance over mutations in genome g2; all mutations in g2 visited or found the
              // one at position m.x

        // add new mutation m to both genomes
        (*g1).push_back(m);
        (*g2).push_back(m);

        // NEW: sorting genomes and removing duplicate mutations is no longer needed under the
          // new mutation scheme, as it allows for at most one mutation at any given position
          // in any haploid genome
        /*
        // OLD:
        // sort the genomes
        sort((*g1).begin(),(*g1).end()); // sorting according to physical position
        sort((*g2).begin(),(*g2).end());

        // unique(ForwardIterator first, ForwardIterator last): removes al but
          // the first element from every consecutive group of equivalent
          // elements in the range [first, last)
        // erase(iterator first, iterator last): removes the range
          // [first, last) of elements

        // remove *consecutive* duplicate mutations, where consecutive referst to physical position
        (*g1).erase(unique((*g1).begin(),(*g1).end()),(*g1).end()); // operator ==: same position,
          // same mutation type, and same selection coefficient (defined for comparison of
          // instances of class mutation
        (*g2).erase(unique((*g2).begin(),(*g2).end()),(*g2).end());
        */

      } // end of for each homozygote mutant individual to be created

    // introduce heterozygotes

    for (int j = 0; j < M.nAa; j++) // for each heterozygote mutant individual to be created
      { 
        // find subpopulation, then return one of its child genomes (randomise over the two
          // genomes in an individual)
        int l = gsl_rng_uniform_int(rng, 2);
        genome *g1 = &find(M.i)->second.G_child[2*ix_het[j] + l]; // recall: a diploid individual k
          // is made up of haploid genomes 2*k and 2*k + 1; returns a pointer to child genome 2*k
          // in subpopulation M.i; here, k = ix_het[j]

        // introduce the desired heterozygote, but ensure that previously existing mutations
          // (neutral ones or those identical in state with the one to be introduced) at position
          // m.x are removed from the haploid genome referred to by g1

        // find mutations present at position m.x in genome g1 and remove them

        // advance over mutations in genome g1
        m_it = (*g1).begin();
        found = 0;

        while (m_it != (*g1).end() && found == 0) // in the new mutation scheme, there may be at
          // most one mutation at any given position in a haploid genome; once one mutation at
          // position m.x. has been found, all have been found
          {
            if ((*m_it).x == m.x)
              {
                // remove currently visited mutation and point iterator to next following mutation
                m_it = (*g1).erase(m_it);
                found = 1;
              }
            else
              {
                // increase iterator
                m_it++;
              }
          } // end of advance over mutations in g1; all mutations in g1 vistited or found the one
              // at position m.x

        // add new mutation m to genome g1
        (*g1).push_back(m);

        // NEW: sorting genomes and removing duplicate mutations is no longer needed under the
          // new mutation scheme, as it allows for at most one mutation at any given position
          // in any haploid genome
        /*
        // OLD:
        // sort the genomes
        sort((*g1).begin(),(*g1).end()); // sorting by physical position

        // remove duplicate mutations
        (*g1).erase(unique((*g1).begin(),(*g1).end()),(*g1).end());
         */

      } // end of for each heterozygote mutant individual to be created

      // add information about currently introduced mutation to chr.seg_nonneut_mut
      vector<double> v;
      v.push_back((double) m.t);
      v.push_back(m.s);
      chr.seg_nonneutr_mut.insert(pair<int,vector<double>>(m.x,v));

  } // end of method introduce_mutation()


  void track_mutations(int g, vector<int>& TM, vector<partial_sweep>& PS, chromosome& chr)
  {
    // output trajectories of followed mutations and set s = 0 for partial sweeps when completed

    // g is the generation, TM refers to the vector of ids of tracked mutation-types, PS to the
      // vector of partial sweeps, and chr to the chromosome

    // find all polymorphism of the types that are to be tracked

    for (it = begin(); it != end(); it++) // iterate over all subpopulations
      {
        multimap<int,polymorphism> P;
        multimap<int,polymorphism>::iterator P_it;

        for (int i = 0; i < 2*it->second.N; i++) // iterate over all children
          {
            for (int k = 0; k < it->second.G_child[i].size(); k++) // iterate over all mutations
              {
                for (int j = 0; j < TM.size(); j++) // iterate over tracked mutation-types
                  {
                    if (it->second.G_child[i][k].t == TM[j]) // if current mutation-type is equal
                      // to mutation-type of currently visited mutation
                      {
                        add_mut(P, it->second.G_child[i][k], it->first);
                      }
                  } // end of iterate over tracked mutation-types
              } // end of iterate over all mutations
          } // end of iterate over all children

        // output the frequencies of these mutations in each subpopulation separately

        for (P_it = P.begin(); P_it != P.end(); P_it++) // iterate over polymorphisms
          {
            cout << "#OUT: " << g << " T p" << it->first << " ";
            P_it->second.print_noi(P_it->first, chr, it->first);
          } // end of iterate over polymorphisms
      } // end of iterate over all subpopulations

    // check partial sweeps

    multimap<int,polymorphism> P; // for the entire population (above, P was specific to a given
      // subpopulation); P will be filled by population::add_mut() such that the key corresponds to
      // the physical position; note that more than one mutation are allowed at a given position in
      // general, but no more than one non-neutral mutation other than those identical in state
    multimap<int,polymorphism>::iterator P_it;

    if (PS.size() > 0) // if there are partial sweeps
      {
        P.clear();

        int N = 0;
        // compute the size of the entire population (summing diploid individuals across all
          // subpopulations)
        for (it = begin(); it != end(); it++)
          { N += it->second.N; }

        // find all polymorphism that are supposed to undergo partial sweeps

        for (it = begin(); it != end(); it++) // iterate over all subpopulations
          {
            for (int i = 0; i < 2*it->second.N; i++) // iterate over all children
              {
                for (int k = 0; k < it->second.G_child[i].size(); k++) // iterate over all
                  // mutations; recall that they are ordered by physical position in a genome
                  {
                    for (int j = 0; j < PS.size(); j++) // iterate over all partial sweeps
                      {
                        // if the currently visited mutation and the current partially sweeping
                          // are identical in mutation-type and position; to prevent confusion
                          // with homoplasy via recurrent mutation, the user shoud specify
                          // exclusive mutation-types for partial sweeps and predetermined
                          // mutations
                        if (it->second.G_child[i][k].x == PS[j].x && it->second.G_child[i][k].t == PS[j].t)
                          {
                            add_mut(P, it->second.G_child[i][k], it.first);
                          }
                      } // end of iterate over all partial sweeps
                  } // end of iterate over all mutations
              } // end of iterate over all children
          } // end of interate over all subpopulations
    
        // check whether a partial sweep has reached its target frequency

        for (P_it = P.begin(); P_it != P.end(); P_it++) // iterate over all polymorphisms
          {
            for (int j = 0; j < PS.size(); j++) // iterate over all partial sweeps
              {
                // if currently visited polymorphism and current partial sweep are identical in
                  // terms of physical position and mutation-type
                if (P_it->first == PS[j].x && P_it->second.t == PS[j].t)
                  {
                    // if the global frequency (i.e. in the entire population) of the (partially)
                      // sweeping mutation has reached the threshold frequency
                    if (((float)P_it->second.calc_ntot())/(2*N) >= PS[j].p)
                      {

                        // sweep has reached target frequency, set all s to zero
			
                        for (it = begin(); it != end(); it++) // iterate over all subpopulations
                          {
                            for (int i = 0; i < 2*it->second.N; i++) // iterate over all children
                              {
                                for (int k = 0; k < it->second.G_child[i].size(); k++) // iterate
                                  // over all mutations
                                  {
                                    // if currently visited mutation and current partial sweep
                                      // are identical in terms of position and mutation-type,
                                      // set the respective selection coefficient to 0.
                                    if (it->second.G_child[i][k].x == PS[j].x && it->second.G_child[i][k].t == PS[j].t)
                                      {
                                        it->second.G_child[i][k].s = 0.0;
                                      }
                                  } // end of for each mutation
                              } // end of for each child
                          } // end of for each subpopulation
                        // erase currently visited (and completed) (partial) sweep and decrease
                          // iterator
                        PS.erase(PS.begin() + j);
                        j--;
                      } // end of if current (partial) sweep has reached its target frequency
                  } // end of if current polymorphism and current partial sweep are identical
              } // end of for all partial sweeps
          } // end of for all polymorphisms
      } // end of if there are partial sweeps
  } // end of method track_mutations()


  void evolve_subpopulation(int i, chromosome& chr)

  // evolving subpopulation indexed by i

  {
    int g1, g2; // indices of the two haploid genomes belonging to a child
    int p1, p2; // indices of the two diploid parents
    int n_mut_1, n_mut_2;

    // create map of shuffled children ids

    int child_map[find(i)->second.N]; // declare array child_map of size equal to size of
      // subpopulation i
    // fill child_map with children ids and shuffle them
    for (int j = 0; j < find(i)->second.N; j++)
      {
        child_map[j] = j;
      }
    gsl_ran_shuffle(rng, child_map, find(i)->second.N, sizeof(int));


    int c = 0; // counter over all N children (will get mapped to child_map[c])

    // migration; loop over all subpopulation from which subpopulation i receives migrants
    
    map<int,double>::iterator it;

    // for each entry in the map of migration rates pertaining to subpopulation i
    for (/* map<int,double>::iterator*/ it = find(i)->second.m.begin(); it != find(i)->second.m.end(); it++)
      {

        // deterministically determine the number of immigrants from subpopulation corresponding
          // to it to subpopulation i
        // Note that in a new version, the number of immigrants will be a random number
        int n_migrants = (int)(it->second * find(i)->second.N + 0.5);

        // for each immigrant to be drawn
        for (int m = 0; m < n_migrants; m++)
          {
            if (c >= find(i)->second.N)
              {
                cerr << "ERROR (evolve subpopulation): too many migrants in subpopulation "<< i << endl;
                exit(1);
              }

            g1 = 2*child_map[c];   // child genome 1
            g2 = 2*child_map[c]+1; // child genome 2

            // draw parents in source population

            p1 = gsl_rng_uniform_int(rng, find(it->first)->second.G_parent.size()/2); // draws the
              // index of a parent in the currently relevant subpopulation at random
            if (gsl_rng_uniform(rng) < find(it->first)->second.S) // if this reproduction event is
              // one of selfing
              {
                p2 = p1;
              }
            else // no selfing; draw index of second parent
              {
                p2 = gsl_rng_uniform_int(rng,find(it->first)->second.G_parent.size()/2);
              }

            // recombination, gene-conversion, mutation

            crossover_mutation(i, g1, it->first, 2*p1, 2*p1+1, chr);
            crossover_mutation(i, g2, it->first, 2*p2, 2*p2+1, chr);

            c++;
          } // end of for each immigrant to be drawn
      } // end of for each entry in the map of migration rates pertaining to subpopulation i
	    
    // remainder, i.e. children not resulting from immigration, but from matings within
      // subpopulation i

    // for each remaining child
    while (c < find(i)->second.N)
      {
        g1 = 2*child_map[c];   // child genome 1
        g2 = 2*child_map[c]+1; // child genome 2

        p1 = find(i)->second.draw_individual(); // parent 1
        if (gsl_rng_uniform(rng) < find(i)->second.S) // if selfing
          {
            p2 = p1;
          }
        else // no selfing
          {
            p2 = find(i)->second.draw_individual();
          }

        crossover_mutation(i, g1, i, 2*p1, 2*p1+1, chr);
        crossover_mutation(i, g2, i, 2*p2, 2*p2+1, chr);

        c++;
      } // end of for each remaining child
  } // end of method evolve_subpopulation()


  void crossover_mutation(int i, int c, int j, int P1, int P2, chromosome& chr)
  {
    // child genome c in subpopulation i is assigned outcome of cross-overs at breakpoints r 
    // between parent genomes P1 and P2 from subpopulation j and new mutations are added
    // 
    // example R = (r1, r2)
    // 
    // mutations (      x < r1) assigned from P1
    // mutations (r1 <= x < r2) assigned from P2
    // mutations (r2 <= x     ) assigned from P1
    //
    // p1 and p2 are swapped in half of the cases to assure random assortement

    if (gsl_rng_uniform_int(rng, 2) == 0)
      { int swap = P1; P1 = P2; P2 = swap; } // swap P1 and P2

    // clear child genome c in subpopulation i
    find(i)->second.G_child[c].clear();

    // create vector with the mutations to be added

    vector<mutation> M; // storing de-novo mutations arising in child c
    map<int,mutation> M_map; // helper construct to store de-novo mutations arising in child c,
      // making it easier to guarantee that there is at most one mutation at any physical positin;
      // the physical position is the key;
    map<int,mutation>::iterator M_map_it;

    int n_mut = chr.draw_n_mut(); // draw random number of mutations

    // introduce n_mut new mutations according to genomic elements and relative importance of
      // respective mutation-types
    for (int k = 0; k < n_mut; k++)
      {
        // OLD: M.push_back(chr.draw_new_mut());
        // draw a new mutation and modify M such as to include the new mutation; note that the
          // mutation will not be introduced if it were to choose a position at which there already
          // is a segregating non-neutral mutation, unless that segregating non-neutral mutation is
          // identical in state with the new mutation (mutation-type and identical selection
          // coefficient w.r.t. the reference environment); tis is ensured by the method
          // chromosome::draw_new_mut()
        chr.draw_new_mut(M_map); // the method draw_new_mut() does not allow for more than
          // one mutation of any type at a given physical position; older mutations are overridden
          // by newer ones
      }

    // copy values from M_map to M, recalling that the values (mutations) in M_map are sorted by
      // the key corresponding to the physical position

    M_map_it = M_map.begin();

    // for each element in M_map
    while (M_map_it != M_map.end())
      {
        M.push_back(M_map_it->second());
        M_map_it++;
      } // end of for each element in M_map
    // mutations in M are now sorted according to physical position

    // OLD:
    // sort(M.begin(), M.end()); // sort by physical position

    // TODO: GO ON HERE, accounting for the fact that M contains at most one mutation at any given
      // physical position x

    // create vector with recombination breakpoints

    vector<int> R = chr.draw_breakpoints(); // breakpoints are physical positions
    R.push_back(chr.L + 1);
    sort(R.begin(), R.end());
    R.erase(unique(R.begin(), R.end()), R.end()); // erase consecutive duplicate breakpoints

    // initialize various iterators

    vector<mutation>::iterator p1 = find(j)->second.G_parent[P1].begin(); // pointing to the first
      // mutation in parent genome P1
    vector<mutation>::iterator p2 = find(j)->second.G_parent[P2].begin(); // ditto for parent
      // genome P2

    vector<mutation>::iterator p1_max = find(j)->second.G_parent[P1].end(); // pointing to the last
      // mutation in parent genome P1
    vector<mutation>::iterator p2_max = find(j)->second.G_parent[P2].end(); // ditto for parent
      // genome P2

    vector<mutation>::iterator m     = M.begin(); // pointing to the first mutation added above
    vector<mutation>::iterator m_max = M.end(); // pointing to the last mutation added above

    vector<mutation>::iterator p     = p1; // pointing to the first mutation in parent genome P1;
      // later updated to pointing to currently visited mutation in P1
    vector<mutation>::iterator p_max = p1_max; // pointing to the last mutation in parent genome P1

    int r = 0; // counter of visited recombination breakpoints
    int r_max = R.size(); // total number of recombination breakpoints
    int n = 0; // number of mutations (those transmitted from a parent or de-novo ones) added
      // to child c
    bool present;

    while (r != r_max) // while there are further recombination breakpoints
      {
        while ((p != p_max && (*p).x < R[r]) || (m != m_max && (*m).x < R[r])) // while there are
          // previously existing or de-novo mutations prior to the current recombination breakpoint
          {
            // advance p, and check if the previously existing (parental) mutation referred to by p
              // can be added to child genome c; it can be added only if there is not yet a
              // mutation at the respective position in c, as such a mutation would necessarily
              // have to be a de-novo mutation arising in this generation, which is given priority
              // over previously existing parental mutations
            while (p != p_max && (*p).x < R[r] && (m == m_max || (*p).x <= (*m).x)) // while there
              // are previously existing mutations prior to the current recombination breakpoint
              // and (all de-novo mutations have been visited or the position of the currently
              // visited previously existing mutation is not larger than the one of the current de-
              // novo mutation)

              {
                present = 0;
                if (n != 0 && find(i)->second.G_child[c].back().x == (*p).x) // if at least one
                  // previously existing or de-novo mutation has been added and the position of
                  // the last mutation in child genome c is identical with the position of the
                  // currently visited parental mutation
                  {

                    // NEW:

                    // there is at most one mutation at any given position x in M, P1, and P2;
                      // hence, if .back() of child c is at position x, it is the only one
                    present = 1;

                    // OLD:
                    /*
                    // search the child genome c for mutations at the same position as (*p).x;
                      // according to the new mutation regime, there may be at most one such
                      // mutation

                    int k = n - 1;
                    while (present == 0 && k >= 0) // while no mutation has been found in child
                      // genome c at position (*p).x and there are further mutations in child
                      // genome c that have been added in previous iterations of this loop
                      {
                        if (find(i)->second.G_child[c][k].x == (*p).x) // if the currently visited
                          // parental mutation (p) is at the same position as mutation k in child c

                          // mutation k in child c has priority over mutation p from the parent, as
                            // it is more recent; p must not be added

                          {
                            present = 1;
                          }
                        // decrease iterator, as we want to check the preceding mutation, added
                          // to child genome c in a previous iteration
                        k--;
                      } // encountered a mutation at position (*p).x in child genome c or there
                          // are no more mutations added to child genome c to be checked

                  } // there is no previously existing or de-novo mutation in child genome c or the
                      // position of the last mutation in child genome c is not identical with the
                      // position of the currently visited parental mutation or
                  */

                  } // if a de-novo mutation in child genome c was present at position x, it has
                    // been identified

                if (present == 0) // if there is no mutation present in child c at the same
                  // position as the currently visited parental mutation resides, the parental
                  // mutation can be added to c
                  {
                    find(i)->second.G_child[c].push_back(*p);
                    n++;
                  }
                p++; // shift pointer to next parental mutation
              } // there are no more previously existing (parental) mutations prior to the current
                  // recombination breakpoint or all de-novo mutations have been visited or the
                  // position of the currently visited parental mutation is not larger than the one
                  // of the current de-novo mutation

            // advance m and add the de-novo mutation referred to by m

            while (m != m_max && (*m).x < R[r] && (p == p_max || (*m).x <= (*p).x)) // while there
              // are de-novo mutations prior to the current recombination breakpoint and (all
              // previously existing mutations have been visited or the position of the currently
              // visited de-novo mutation is not larger than the one of the current parental
              // mutation
              {
                present = 0;
                if (n != 0 && find(i)->second.G_child[c].back().x == (*m).x) // if at least one
                  // previously existing or de-novo mutation has been added and the position of
                  // the last mutation in child genome c is identical with the position of the
                  // currently visited de-novo mutation
                  {

                    int k = n - 1;

                    // NEW:

                    // there is at most one mutation at any given position x in M, P1, and P2;
                      // hence, if .back() of child c is at position x, it is the only one
                    present = 1;

                    // the previously existing mutation at .back().x in c must be replaced
                      // by the currently visited de-novo mutation

                    find(i)->second.G_child[c][k] = (*m);

                    // no need to increase n, as we only replaced an existing mutation

                    // OLD:
                    /*
                    // search the child genome c for mutations at the sampe osition as (*m).x
                    while (present == 0 && k >= 0) // while no mutation has been found in child
                      // genome c at position (*m).x and there are further mutations in child
                      // genome c that have been added in previous iterations of the loop
                      {
                        // OLD:
                        // if (find(i)->second.G_child[c][k] == (*m))
                        if (find(i)->second.G_child[c][k].x == (*m).x) // if the currently visited
                          // de-novo mutation (m) is at the same position as mutation k in child c

                          // mutation k in child c must be replaced by the de-novo mutation
                            // m, as de-novo mutations have priority over previously existing
                            // (parental) mutations
                          {
                            present = 1;
                            find(i)->second.G_child[c][k] = (*m);
                            // as we just replace a mutation, we do not need to increase n
                          }
                        // decrease iterator, as we want to check the preceding mutation, added
                          // to child genome c in a previous step
                        k--;
                      } // encountered a mutation at position (*m).x in child genome c or there
                          // are no more mutations added to child genome c to be checked
                  } // there is no previously existing or de-novo mutation in child genome c or the
                      // position of the last mutation in child genome c is not identical with the
                      // position of the currently visited parental mutation
                  */

                  } // if a previously existing mutation in child genome c was present at position
                    // x, it has been identified and replaced by the currently visited de-novo
                    // muation

                if (present == 0) // if there is no mutation present in child c at the same
                  // position as the currently visited de-novo mutation resides, the de-novo
                  // mutation can be added to c
                  {
                    find(i)->second.G_child[c].push_back(*m);
                    n++;
                  }
                m++; // shift pointer to next de-novo mutation
              } // there are no more de-novo mutations prior to the current recombination
                  // breakpoint or all parental mutations have been visited or the position of the
                  // currently visited de-novo mutation is not larger than the one of the current
                  // parental mutation
          } // end of while there are previously existing or new mutations prior to the current
              // recombination breakpoint

        // swap parents

        p1 = p2; p1_max = p2_max;
        p2 = p; p2_max = p_max;
        p = p1; p_max = p1_max;

        // shift iterator p (which now refers to mutations in the other parent) to current
          // breakpoint
        while (p != p_max && (*p).x < R[r])
          { p++; }

        r++;

      } // end of while there are further recombination breakpoints

  } // end of method crossover_mutation()


  void swap_generations(int g, chromosome& chr, char& fi, map<int,environment>& env)
  {

    // perform change of generation at time g and update fitnesses given mutations on chromosome
      // chr assuming fitness interaction fi (additive or multiplicative)

    // find and remove fixed mutations from the children in all subpopulations, and update the
      // register of positions at which non-neutral mutations are segregating w.r.t. the entire
      // population
    
    remove_fixed(g, chr);

    // make children the new parents and update fitnesses

    for (it = begin(); it != end(); it++) // iterating over subpopulations
      { 
        it->second.swap();
        it->second.update_fitness(chr, fi);
      } // end of iterating over subpopulations
  } // end of method swap_generations()


  void remove_fixed(int g, chromosome& chr)
  {
    // find mutations that are fixed in all subpopulations and update the substitutions vector

    // g is the generation; chr is the chromosome

    genome G = begin()->second.G_child[0]; // copy the first genome in the first subpopulation
      // as a reference to which all other genomes across the entire population are compared; this
      // reference will be updated genome by genome, and accumulate fixed mutations

    for (it = begin(); it != end(); it++) // iterate over subpopulations
      {
        for (int i = 0; i < 2*it->second.N; i++) // iterate over child genomes
          {
            // returning a genome consisting only of mutations that are present in both genomes
              // passed as arguments; one is the reference G (which is cumulatively updated), and
              // the other genome is the child currently visited
            G = fixed(it->second.G_child[i], G);
          } // end of for each child genome
      } // end of for each subpopulation

    // G is now a vector of those mutations that are fixed across the entire population

    // if the cumulative genome made up of mutations fixed across the entire population is not
      // empty
    if (G.size() > 0)
      {
        for (it = begin(); it != end(); it++) // iterate over subpopulations
          {
            for (int i = 0; i < 2*it->second.N; i++) // iterate over child genomes
              {
                // returning a genome consisting only of mutations present in the genome given
                  // as the first argument, but not in the genome given as the second argument;
                  // this clears the currently visited child genome of mutations fixed across
                  // the entire population
                it->second.G_child[i] = polymorphic(it->second.G_child[i], G);
              } // end of for each child genome
          } // end of for each subpopulation

        // add substitution to book keeping (a substitution is a fixed mutation); remove respective
          // entry from register that holds position and type of non-neutral mutations segregating
          // in the entire population
        for (int i = 0; i < G.size(); i++) // for each fixed mutation
          {
            // recall that genome G is a vector of mutations

            Substitutions.push_back(substitution(G[i], g));

            // remove from register of segregating non-neutral mutations
            int check = chr.seg_nonneutr_mut.erase(G[i].x);
            if (check != 1 && G[i].s != 0.0) // if no segregating non-neutral mutation has been
              // removed and G[i] is a non-neutral mutation
              {
                cerr << "ERROR (remove_fixed): non-neutral substitution not found in seg_nonneutr_mut" << endl;
                exit(1);
              }
          } // end of for each fixed mutation

      } // end of if genome G is not empty
  } // end of method remove_fixed()


  void print_all(chromosome& chr)
  {
    // print all mutations and all genomes 

    cout << "Populations:" << endl;
    // iterating over subpopulations
    for (it = begin(); it != end(); it++)
      {
        cout << "p" << it->first << " " << it->second.N << endl; // print size of subpopulation
      }

    multimap<int,polymorphism> P;
    multimap<int,polymorphism>::iterator P_it;

    // add all polymorphisms
    for (it = begin(); it != end(); it++) // iterate over all subpopulations
      {
        for (int i = 0; i < 2*it->second.N; i++) // iterate over all children
          {
            for (int k = 0; k < it->second.G_child[i].size(); k++) // iterate over all mutations;
              // note that mutations in a genome (here, G_child[i]) are ordered by the key
              // identical to their physical position
              {
                add_mut(P, it->second.G_child[i][k], it->first); // add mutation k to P or
                  // increase its prevalence if already contained in P; prevalence are made
                  // specifically to the subpopulation
              } // end of iterate over all mutations
          } // end of iterate over all children
      } // end of iterate over all subpopulations

    cout << "Mutations:"  << endl;

    // print mutation id, mutation-type, position, selection coefficient and dominance in the
      // reference environment, total prevalence in the entire population, and subpopulation-
      // specific prevalences
    for (P_it = P.begin(); P_it != P.end(); P_it++)
      {
        P_it->second.print(P_it->first, chr); // the key is the physical position
      }

    cout << "Genomes:" << endl;

    // print all genomes

    for (it = begin(); it != end(); it++) // iterate over all subpopulations
      {
        for (int i = 0; i < 2*it->second.N; i++) // iterate over all children
          {
            cout << "p" << it->first << ":" << i+1;
            // print rest of line
            for (int k = 0; k < it->second.G_child[i].size(); k++) // iterate over all mutations
              {
                // returns id of mutation found in P; if no mutation is found in P, that mutation
                  // is added to P as a polymorphism and the id returned
                int id = find_mut(P, it->second.G_child[i][k]);
                cout << " " << id;
              }
            cout << endl;
          } // end of iterate over children
      } // end of iterate over subpopulations

  } // end of method print_all() 1


  void print_all(ofstream& outfile, chromosome& chr)
  {
    // print all mutations and all genomes 

    outfile << "Populations:" << endl;
    // iterating over subpopulations
    for (it = begin(); it != end(); it++)
      {
        outfile << "p" << it->first << " " << it->second.N << endl;
      }

    multimap<int,polymorphism> P;
    multimap<int,polymorphism>::iterator P_it;

    // add all polymorphisms
    for (it = begin(); it != end(); it++) // iterate over all subpopulations
      {
        for (int i = 0; i < 2*it->second.N; i++) // iterate over all children
          {
            for (int k = 0; k < it->second.G_child[i].size(); k++) // iterate over all mutations
              {
                add_mut(P, it->second.G_child[i][k], it.first);
              } // end of iterate over all mutations
          } // end of iterate over all children
      } // end of iterate over all subpopulations

    outfile << "Mutations:"  << endl;

    // write mutation id, mutation-type, position, selection coefficient and dominance in the
      // reference environment, total prevalence in the entire population, and subpopulation-
      // specific prevalences
    for (P_it = P.begin(); P_it != P.end(); P_it++)
      {
        P_it->second.print(outfile, P_it->first, chr);
      }

    outfile << "Genomes:" << endl;

    // write all genomes

    for (it = begin(); it != end(); it++) // iterate over all subpopulations
      {
        for (int i = 0; i < 2*it->second.N; i++) // iterate over all children
          {
            outfile << "p" << it->first << ":" << i+1;
            // write rest of line
            for (int k = 0; k < it->second.G_child[i].size(); k++) // iterate over all mutations
              {
                // returns id of mutation found in P; if no mutation is found in P, that mutation
                  // is added to P as a polymorphism and the id returned
                int id = find_mut(P, it->second.G_child[i][k]);
                outfile << " " << id;
              }
            outfile << endl;
          } // end of iterate over all children
      } // end of iterate over all subpopulations

  } // end of method print_all() 2


  void print_sample(int i, int n, chromosome& chr)
  {
    // print sample of n genomes from subpopulation  i

    if (count(i) == 0)
      {
        cerr << "ERROR (output): subpopulation p"<< i << " does not exists" << endl; exit(1);
      }

    vector<int> sample; // collecting indices of sampled genomes in subpopulation i

    multimap<int,polymorphism> P;
    multimap<int,polymorphism>::iterator P_it;
    
    // iterate over random genomes to be drawn
    for (int s = 0; s < n; s++)
      {
        // draw a random index for a haploid genome in subpopulation i
        int j = gsl_rng_uniform_int(rng, find(i)->second.G_child.size());
        sample.push_back(j);

        for (int k = 0; k < find(i)->second.G_child[j].size(); k++) // iterate over all mutations
          // present in genome j
          {
            add_mut(P, find(i)->second.G_child[j][k], find(i)->first);
          } // end of for each mutation
      } // end of for each random genome to be drawn

    cout << "Mutations:"  << endl;

    // print mutation id, mutation-type, position, selection coefficient and dominance in the
      // reference environment, total prevalence in the entire population, and subpopulation-
      // specific prevalences
    for (P_it = P.begin(); P_it != P.end(); P_it++)
      {
        P_it->second.print(P_it->first, chr);
      }

    cout << "Genomes:" << endl;

    // print all genomes

    for (int j = 0; j < sample.size(); j++) // iterate over all children (i.e. haploid genomes
      // drawn)
      {
        cout << "p" << find(i)->first << ":" << sample[j] + 1; // write subpopulation id and
          // genome index

        for (int k=0; k<find(i)->second.G_child[sample[j]].size(); k++) // iterate over all
          // mutations
          {
            int id = find_mut(P, find(i)->second.G_child[sample[j]][k]);
            cout << " " << id;
          } // end of iterate over all mutations
        cout << endl;
      } // end of iterate over all children
  } // end of method print_sample()


  void print_sample_ms(int i, int n, chromosome& chr)
  {
    // print sample of n genomes from subpopulation  i

    if (count(i) == 0)
    {
      cerr << "ERROR (output): subpopulation p"<< i << " does not exists" << endl;
      exit(1);
    }

    vector<int> sample; // collecting indices of sampled genomes in subpopulation i

    multimap<int,polymorphism> P;
    multimap<int,polymorphism>::iterator P_it;

    // iterate over random genomes to be drawn
    for (int s = 0; s < n; s++)
      { 
        // draw a random index for a haploid genome in subpopulation i
        int j = gsl_rng_uniform_int(rng, find(i)->second.G_child.size());
        sample.push_back(j);

        for (int k = 0; k < find(i)->second.G_child[j].size(); k++) // iterate over all mutations
          {
            add_mut(P, find(i)->second.G_child[j][k], find(i)->first);
          } // end of for each mutation
      } // end of for each random genome to be drawn

    // print header

    cout << endl << "//" << endl << "segsites: " << P.size() << endl;

    // print all positions

    if (P.size() > 0) // if there are polymorphisms
      {
        cout << "positions:";

        for (P_it = P.begin(); P_it != P.end(); P_it++) // iterate over polymorphisms
          {
            cout << " " << fixed << setprecision(7) << (double)(P_it->first+1)/(chr.L + 1);
          } // end of iterate over polymorphisms
        cout << endl;
      } // end of if there are polymorphisms

    // print genotypes

    for (int j = 0; j < sample.size(); j++) // iterate over all children
      {
        string genotype(P.size(), '0');

        for (int k = 0; k < find(i)->second.G_child[sample[j]].size(); k++) // iterate over all
          // mutations
          {
            int pos = 0;
            mutation m = find(i)->second.G_child[sample[j]][k];

            for (P_it = P.begin(); P_it != P.end(); P_it++)
              {
                if (P_it->first == m.x && P_it->second.t == m.t && P_it->second.s == m.s)
                  {
                    genotype.replace(pos, 1, "1");
                    break;
                  }
                pos++;
              } // end of for each polymorphism
          } // end of for each mutation
        cout << genotype << endl;
      } // end of for each child
  } // end of method print_sample_ms()

  int find_mut(multimap<int,polymorphism>& P, mutation m)
  {
    // find m in P and return its id

    int id = 0;

    // iterate through all mutations with same position

    multimap<int,polymorphism>::iterator it;
    pair<multimap<int,polymorphism>::iterator,multimap<int,polymorphism>::iterator> range = P.equal_range(m.x);
    it = range.first;
    // advance through polymorphisms
    while (it != range.second)
      {
        // if mutation m and currently visited polymorphism are identical in terms of mutation-type
          // and selection coefficient
        if (it->second.t == m.t && it->second.s == m.s)
          {
            id = it->second.i;
            // it->second.n++; // A mistake in the original version of SLiM
            it = range.second;
          }
        else
          { it++; }
      } // end of advance through polymorphisms
    
    return id;
  } // end of method find_mut()


  void add_mut(multimap<int,polymorphism>& P, mutation m, int spid)
  {
    // if mutation is present in P (i.e. in the entire population) increase prevalence in the
      // appropriate subpopulation pop, otherwise add it to P for the subpopulation with key spid

    int id = 0; // the value assigned to no existing mutation; ids start at 1

    // iterate over all mutations with same position

    multimap<int,polymorphism>::iterator it;
    pair<multimap<int,polymorphism>::iterator,multimap<int,polymorphism>::iterator> range = P.equal_range(m.x); // establish the range of polymorphisms in P with the same physical
        // position x;
    it = range.first; // initialise iterator to first mutation at position x

    // iterate over mutations (polymorphisms) at the same physical position
    while (it != range.second)
      {
        // if mutation m is identical in state (i.e. type and selection coefficient in the
          // reference environment) with the currently visited polymorphism
        if (it->second.t == m.t && it->second.s == m.s)
          {
            id = it->second.i; // obtain id of mutation already present at x
            // OLD
            // it->second.n++; // increase counter for subpopulation with key spid
            // NEW
            // if key spid exists, increase corresponding counter, otherwise insert that counter
            map<int,int>::iterator count_it;
            count_it = it->second.n.find(spid);

            if (count_it == second.n.end()) // counter does not yet exist
              {
                // insert counter
                it->second.n.insert(pair<int,int>(spid,1));
              }
            else // counter exists
              {
                // update counter for subpopulation with id spid
                count_it->second++;
              }

            it = range.second; // set iterator equal to end of last mutation at x
          }
        else // mutation m is not identical in state with currently visited polymorphism
          { it++; } // advance to next polymorphism
      } // end of iterate over polymorphism at the same physical position

    // if not already present in any subpopulation, add polymorphism to P

    if (id == 0)
      {
        id = P.size() + 1;
        P.insert(pair<int,polymorphism>(m.x, polymorphism(id, m.t, m.s, pair<int,int>(spid,1)))); // the last argument to the constructor of polymorphism is the map of prevalences;
          // here, we initalise the polymorphism with just one count in subpopulation spid
      }
  } // end of add_mut() method

  void remove_seg_mut(int x)
  {

    // removes mutations segregating at position x from the entire population

    vector<mutation>::iterator g;

    // iterate over all subpopulations
    for (it = begin(); it != end(); it++)
      {
        // iterate over all child genomes in current subpopulation
        for (int j = 0; j < 2*it->second.N; j++)
          {
            g = it->second.G_child[j].begin();
            int found = 0;

            // advance over all mutations in child genome j
            while (g != it->second.G_child[j].end() && found == 0) // in the new mutation scheme,
              // there may be at most one mutation at any given position in a haploid genome; once
              // one mutation at position m.x has been found, all have been found
              {
                // if currently visited mutation is located at position x
                if ((*g).x == m.x)
                  {
                    // remove currently visited mutation and point iterator to next
                      // following mutation
                    g = it->second.G_child[j].erase(g);
                    found = 1;
                  }
                else // currently visited mutation is not located at position m.x
                  {
                    // increase iterator
                    g++;
                  }
              } // end of while there are mutations to be checked
          } // end of for each child genome
      } // end of for each subpopulation
  } // end of method remove_seg_mut()

}; // end of class 'population'



void get_line(ifstream& infile, string& line)
{
  getline(infile, line);
  if (line.find("/")!= string::npos)
    { line.erase(line.find("/")); } // remove all after "/"; these
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
                          iss >> sub; // s-modif (modifier of selection coefficient)
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


void initialize_from_file(population& P, const char* file, chromosome& chr, char& fi, const environment& renv)
{
  // initialize population from file

  map<int,mutation> M;

  string line; 
  string sub; 

  ifstream infile (file);

  if (!infile.is_open()) { cerr << "ERROR (initialize): could not open initialization file" << endl; exit(1); }

  get_line(infile, line);

  // advancing to populations section
  while (line.find("Populations") == string::npos && !infile.eof())
    {
      get_line(infile, line);
    } // populations section reached

  get_line(infile, line);

  // advancing to mutations, reading subpopulations
  while (line.find("Mutations") == string::npos && !infile.eof())
    { 
      istringstream iss(line);
      iss >> sub;
      sub.erase(0, 1); // erasing leading "p"
      int i = atoi(sub.c_str()); // subpopulation identifier
      iss >> sub;
      int n = atoi(sub.c_str()); // subpopulation size
      P.add_subpopulation(i, n, renv); // assigning reference environment renv by default
      get_line(infile, line);      
    } // subpopulations read, mutations section reached

  get_line(infile, line);

  // advancing to genomes, reading mutations
  while (line.find("Genomes") == string::npos && !infile.eof()) 
    {     
      istringstream iss(line);
      iss >> sub;
      int   i = atoi(sub.c_str()); // id
      iss >> sub;
      sub.erase(0, 1); // erasing leading "m"
      int   t = atoi(sub.c_str()); // mutation type
      iss >> sub;
      int   x = atoi(sub.c_str()) - 1; // physical position (internally starting at 0)
      iss >> sub;
      float s = atof(sub.c_str()); // selection coefficient

      M.insert(pair<int,mutation>(i,mutation(t, x, s)));
      get_line(infile, line); 
    } // mutations read, genomes section reached

  get_line(infile, line);

  // advancing to end of file, reading genomes
  while (!infile.eof()) // while end of file not reached
    {
      istringstream iss(line);
      iss >> sub;
      sub.erase(0, 1); // erasing leading "p"
      int pos = sub.find_first_of(":"); // helper giving position of ":"
      int p = atoi(sub.substr(0, pos).c_str()); // numeric subpopulation identifier
      sub.erase(0, pos+1);
      int i = atoi(sub.c_str()); //

      while (iss >> sub) 
        {
          int id = atoi(sub.c_str());
          P.find(p)->second.G_parent[i-1].push_back(M.find(id)->second);
        }
      get_line(infile, line);
    } // end of while end of file not reached

  for (P.it = P.begin(); P.it != P.end(); P.it++)
    {
      P.it->second.update_fitness(chr, fi);
    }
} // end of initialize_from_file()

/*
  Initialises the population P using parameter values to be read from file. Initialises chromosome
  chr, as well as demography and structure events E, output events O, user-defined mutations  IM
  that will be introduced, and mutations PS undergoing partial sweeps. Moreover, fi denotes the
  fitness interaction, env the vector of environments, and renv the reference environment
*/
void initialize(population& P, char* file, chromosome& chr, int& T, char& fi, multimap<int,event>& E, multimap<int,event>& O, multimap<int,introduced_mutation>& IM, vector<partial_sweep>& PS, map<int,environment>& env, environment& renv, vector<string>& parameters)
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
                          cerr << "ERROR (initialize): mutation type "<< i << " already defined" << endl;
                          exit(1);
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
                      fi = sub.at(0); // assign fitness interaction, with fi being a reference to a
                        // global variable
                    } // end of if line is not empty
                  get_line(infile, line);
                } // end of while not hitting next input section or end of file
              // TODO:
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
                      E.insert(pair<int,event>(t, event(c, s))); // add demography and structure
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
                        // (identifier, mut type, dominance coeff., modifier of selection coeff.,
                        // [mut type, dominance coeff., modifier of selection coeff., ...])
                      int i; vector<int> m; vector<double> h; vector<double> smodif;
                      istringstream iss(line);

                      iss >> sub;
                      sub.erase(0, 1); // erasing leading "e"
                      i = atoi(sub.c_str()); // initialising environemt id
                      if (env.count(i) > 0) // assessed based on the key; if environment already
                        // present
                        {
                          cerr << "ERROR (initialize): environment " << i << " already defined" << endl;
                          exit(1);
                        } // environment not yet present
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
                      env.insert(pair<int,environment>(i, environment(chr, m, h, smodif)));
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
                      c = sub.at(0); // event type can be A, R, F
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
                      initialize_from_file(P, sub.c_str(), chr, fi, renv);
                    } // end of if line is not empty
                  get_line(infile, line);
                } // end of while not hitting next section or end of file
            } // end of initialization section
        } // this line was not the start of a new input section
      else
        {
          get_line(infile, line);
          cout << "getting new line:" << endl; // test
          // TODO: GO ON HERE (1), testing error (endless loop)
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
  stringstream ss;
  ss << seed;
  parameters.push_back(ss.str());

  // parameter output

  for (int i = 0; i < P.parameters.size(); i++)
    {
      cout << parameters[i] << endl;
    }
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

  char FI; // fitness interaction ('a' for additive; 'm' for multiplicative)

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

  // user-defined environments (all except the reference environment)

  map<int,environment> EV; // unique keys
  map<int,environment>::iterator itEV;

  // reference environment
  environment REV;

  // tracked mutation-types

  vector<int> TM; 

  // mutations undergoing partial sweeps

  vector<partial_sweep> PS;

  initialize(P, input_file, chr, T, FI, E, O, IM, PS, EV, REV, P.parameters);
 
  // evolve over t generations

  for (int g = 1; g <= T; g++)
    { 
      // execute demographic and substructure events in this generation 

      pair<multimap<int,event>::iterator,multimap<int,event>::iterator> rangeE = E.equal_range(g);

      for (itE = rangeE.first; itE != rangeE.second; itE++)
        {
          P.execute_event(itE->second, g, chr, TM, REV, EV);
        }
   
      // evolve all subpopulations

      for (itP = P.begin(); itP != P.end(); itP++) // for each subpopulation
        {
          // GO ON HERE (2): understand, and extend.
          P.evolve_subpopulation(itP->first, chr);
        } // end of for each subpopulation

      // introduce user-defined mutations
        
      pair<multimap<int,introduced_mutation>::iterator,multimap<int,introduced_mutation>::iterator> rangeIM = IM.equal_range(g);

      // for all mutations to be introduced in the current generation
      for (itIM = rangeIM.first; itIM != rangeIM.second; itIM++)
        {
          // test
          // P.print_all(chr);
          
          P.introduce_mutation(itIM->second, chr);
          
          // test
          // P.print_all(chr);
        } // end of for all mutations to be introduce in the current generation

      // execute output events

      // range of output events happening in generation g
      pair<multimap<int,event>::iterator,multimap<int,event>::iterator> rangeO = O.equal_range(g);

      for (itO = rangeO.first; itO != rangeO.second; itO++)
        {
          P.execute_event(itO->second, g, chr, TM, REV, EV);
        }

      // track particular mutation-types and set s = 0 for partial sweeps when completed
      
      if (TM.size() > 0 || PS.size() > 0)
      {
        P.track_mutations(g, TM, PS, chr);
      }

      // swap generations

      P.swap_generations(g, chr, FI, EV);
    }
} // end of method main()
