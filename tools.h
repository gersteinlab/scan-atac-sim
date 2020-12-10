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

int countlines(char *filename);
inline void fast_srand(int seed);
inline int fast_rand(void);
long double rexp(long double scale);
vector<int> sample_int_expj(int n, int size, vector<long double> prob, int trial_seed);

char* getCmdOption(char ** begin, char ** end, const std::string & option);
bool cmdOptionExists(char** begin, char** end, const std::string& option);
void cmdline_error_message();
