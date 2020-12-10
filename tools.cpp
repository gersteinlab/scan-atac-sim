#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <float.h>
#include <algorithm>

#include <bits/stdc++.h> 
using namespace std;

int countlines(char *filename)
{
  // count the number of lines in the file called filename                                    
  FILE *fp = fopen(filename,"r");
  int ch=0;
  int lines=0;

  if (fp == NULL){
    printf("file pointer NULL\n");
    return 0;
  }

  lines++;
  while ((ch = fgetc(fp)) != EOF)
  {
    if (ch == '\n')
    lines++;
  }
  fclose(fp);
  return lines;
}


static unsigned int g_seed;

// Used to seed the generator.           
inline void fast_srand(int seed) {
    g_seed = seed;
}

// Compute a pseudorandom integer.
// Output value in range [0, 32767]
inline int fast_rand(void) {
    g_seed = (214013*g_seed+2531011);
    return (g_seed>>16)&0x7FFF;
}

long double rexp(long double scale)
{
    long double u;
    u = fast_rand() / (32767 + 1.0);
    return -log(1-u) / scale;
}

vector<int> sample_int_expj(int n, int size, vector<long double> prob, int trial_seed) {

  fast_srand(trial_seed);

  // Corner case
  if (size == 0) {
    vector<int> empty;
    return empty;
  }

  // Step 1: The first m items of V are inserted into R
  // Step 2: For each item v_i ∈ R: Calculate a key k_i = u_i^(1/w),
  // where u_i = random(0, 1)
  // (Modification: Calculate and store -log k_i = e_i / w where e_i = exp(1),
  //  reservoir is a priority queue that pops the *maximum* elements)
  priority_queue<pair<long double, int> > R;

  for (vector<long double>::iterator iprob = prob.begin(); 
    iprob != prob.begin() + size; ++iprob) {
    long double k_i = rexp(1.0) / (*iprob);
    // long double k_i = rexp(*iprob);
    // printf("%f\n", k_i);
    // printf("%d\n", iprob - prob.begin() + 1);
    R.push(make_pair(k_i, iprob - prob.begin() + 1));
  }

  // Step 4: Repeat Steps 5–10 until the population is exhausted
  {
    // Step 3: The threshold T_w is the minimum key of R
    // (Modification: This is now the logarithm)
    // Step 10: The new threshold T w is the new minimum key of R
    const pair<long double, int>& T_w = R.top();


    // Incrementing iprob is part of Step 7
    for (vector<long double>::iterator iprob = prob.begin() + size; iprob != prob.end(); ++iprob) {

      // Step 5: Let r = random(0, 1) and X_w = log(r) / log(T_w)
      // (Modification: Use e = -exp(1) instead of log(r))
      long double X_w = rexp(1.0) / T_w.first;

      // Step 6: From the current item v_c skip items until item v_i, such that:
      long double w = 0.0;

      // Step 7: w_c + w_{c+1} + ··· + w_{i−1} < X_w <= w_c + w_{c+1} + ··· + w_{i−1} + w_i
      for (; iprob != prob.end(); ++iprob) {
        w += *iprob;
        if (X_w <= w)
          break;
      }

      // Step 7: No such item, terminate
      if (iprob == prob.end())
        break;

      // Step 9: Let t_w = T_w^{w_i}, r_2 = random(t_w, 1) and v_i’s key: k_i = (r_2)^{1/w_i}
      // (Mod: Let t_w = log(T_w) * {w_i}, e_2 = log(random(e^{t_w}, 1)) and v_i’s key: k_i = -e_2 / w_i)
      long double t_w = -T_w.first * *iprob;
      long double exp_tw = exp(t_w);
      long double e_2 = log((((long double)fast_rand()/(long double) 32767) * (1.0-exp_tw)) + exp_tw);
      long double k_i = -e_2 / *iprob;

      // Step 8: The item in R with the minimum key is replaced by item v_i
      R.pop();
      R.push(make_pair(k_i, iprob - prob.begin() + 1));

    }
  }

  vector<int> ret;
  ret.resize(size);

  for (vector<int>::iterator iret = ret.end(); iret != ret.begin(); ) {
    --iret;

    if (R.empty()) {
       cout << "Reservoir empty before all elements have been filled";
       exit(0);
    }

    *iret = R.top().second;
    R.pop();
  }

  if (!R.empty()) {
    cout << "Reservoir empty before all elements have been filled";
    exit(0);
  }

  return ret;
}


char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}

void cmdline_error_message()
{
  printf("./weighted_sampling -f <foreground.bed> -b <background.bed>");
  printf(" -of <output foreground> -ob <output background> -n <frag_number> (-nv <fragment number variance>) -c <cell_number>");
  printf(" -s <signal to noise ratio> (-u/-g) (-min <min threshold>) (-max <max threshold>) (-h)");
  exit(1);
}
