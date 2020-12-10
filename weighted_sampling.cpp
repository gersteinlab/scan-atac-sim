#include "tools.h"

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

int main(int argc, char *argv[])
{

  //detect for mandatory command line options
  if (!cmdOptionExists(argv, argv+argc, "-f") || 
    !cmdOptionExists(argv, argv+argc, "-b") || 
    !cmdOptionExists(argv, argv+argc, "-of") || 
    !cmdOptionExists(argv, argv+argc, "-ob") || 
    !cmdOptionExists(argv, argv+argc, "-n") || 
    !cmdOptionExists(argv, argv+argc, "-c") || 
    !cmdOptionExists(argv, argv+argc, "-s") || 
    cmdOptionExists(argv, argv+argc, "-h")) {
    cmdline_error_message();
  }

  //detect for gaussian fragment distribution sampling
  bool uniform_sampling = true;
  unsigned long min_threshold = 0;
  unsigned long max_threshold = 0;
  float frag_num_var = 0.0;
  if(!cmdOptionExists(argv, argv+argc, "-u"))
  {
    uniform_sampling = false;
    min_threshold = atoi(getCmdOption(argv, argv + argc, "-min"));
    max_threshold = atoi(getCmdOption(argv, argv + argc, "-max"));
    frag_num_var = atof(getCmdOption(argv, argv + argc, "-nv"));
  }

  //temporary storage variables
  char chrom[6];
  string chrom_s;
  unsigned long start;
  unsigned long end;
  char name[20];
  string name_s;
  unsigned long number;

  //load command line arguments
  char *foreground_fn = getCmdOption(argv, argv + argc, "-f");
  char *background_fn = getCmdOption(argv, argv + argc, "-b");
  char *output_foreground_name = getCmdOption(argv, argv + argc, "-of");
  char *output_background_name = getCmdOption(argv, argv + argc, "-ob");
  unsigned long fragment_number = atoi(getCmdOption(argv, argv + argc, "-n"));
  unsigned long cell_number = atoi(getCmdOption(argv, argv + argc, "-c"));
  float signal2noise = atof(getCmdOption(argv, argv + argc, "-s"));

  // read in foreground peaks
  unsigned long peak_file_len = countlines(foreground_fn) - 1;
  fprintf(stdout, "%zu\n", peak_file_len);

  vector<long double> WF;
  WF.resize(peak_file_len);
  vector<tuple<string, unsigned long, unsigned long, string>> Foreground;
  Foreground.resize(peak_file_len);

  FILE* foreground = fopen(foreground_fn, "r+");
  if (foreground == NULL) {
    fprintf(stderr, "fail to open file\n");
    exit(1);
  }

  unsigned long filterd_file_len = 0;
  for (unsigned long i = 0; i < peak_file_len; i++) {

      int got = fscanf(foreground, "%s %zu %zu %s %zu", 
        &chrom, &start, &end, &name, &number);
      if (number != 0) {
        WF[filterd_file_len] = number;

        chrom_s = chrom;
        name_s = name;

        Foreground[filterd_file_len] = tuple<string, unsigned long, unsigned long, string> {chrom_s, start, end, name_s};
        filterd_file_len++;
      } 
  }
  fprintf(stdout, "%zu\n", filterd_file_len);
  WF.resize(filterd_file_len);
  Foreground.resize(filterd_file_len);

  fclose(foreground);
  // printf("finish reading in foreground of %zu lines\n", filterd_file_len);

  // calculated the average region weight
  // FILE *foreground_weights=fopen("foreground_weights.txt", "w");
  // if (foreground_weights == NULL) {
  //   fprintf(stderr, "could not open output file\n");
  // }
  unsigned long sum = 0;
  for (vector<long double>::iterator i = WF.begin(); i != WF.end(); ++i) {
      sum += *i;
      // fprintf(foreground_weights, "%Lf\n", *i);
  }
  // fprintf(stdout, "%zu\n", sum);
  for (vector<long double>::iterator i = WF.begin(); i != WF.end(); ++i) {
      *i = *i / sum;
      //fprintf(foreground_weights, "%.15Lf\n", *i);
  }
  //fclose(foreground_weights);

  // read in foreground peaks
  peak_file_len = countlines(background_fn) - 1;
  fprintf(stdout, "%zu\n", peak_file_len);

  vector<long double> WB;
  WB.resize(peak_file_len);
  vector<tuple<string, unsigned long, unsigned long, string>> Background;
  Background.resize(peak_file_len);

  FILE* background = fopen(background_fn, "r+");
  if (background == NULL) {
    fprintf(stderr, "fail to open file\n");
    exit(1);
  }

  filterd_file_len = 0;
  for (unsigned long i = 0; i < peak_file_len; i++) {

      int got = fscanf(background, "%s %zu %zu %s %zu", 
        &chrom, &start, &end, &name, &number);
      if (number != 0) {
        WB[filterd_file_len] = number;

        chrom_s = chrom;
        name_s = name;

        Background[filterd_file_len] = tuple<string, unsigned long, unsigned long, string> {chrom_s, start, end, name_s};
        filterd_file_len++;
      } 
  }

  fprintf(stdout, "%zu\n", filterd_file_len);
  WB.resize(filterd_file_len);
  Background.resize(filterd_file_len);

  fclose(background);
  // printf("finish reading in background of %zu lines\n", filterd_file_len);

  // FILE *background_weights=fopen("background_weights.txt", "w");
  // if (background_weights == NULL) {
  //   fprintf(stderr, "could not open output file\n");
  // }
  // calculated the average region weight
  sum = 0;
  for (vector<long double>::iterator i = WB.begin(); i != WB.end(); ++i) {
      sum += *i;
      // fprintf(background_weights, "%Lf\n", *i);
  }
  // fprintf(stdout, "%zu\n", sum);
  for (vector<long double>::iterator i = WB.begin(); i != WB.end(); ++i) {
      *i = *i / sum;
      // fprintf(background_weights, "%.15Lf\n", *i);
  }
  //fclose(background_weights);

  // printf("read in both foreground and background peaks\n");

  // open output file
  FILE *output_foreground=fopen(output_foreground_name, "w");
  if (output_foreground == NULL) {
    fprintf(stderr, "could not open output file\n");
  }

  FILE *output_background=fopen(output_background_name, "w");
  if (output_background == NULL) {
    fprintf(stderr, "could not open output file\n");
  }

  //instantiate fragment number of each cell if not uniform
  unsigned long cell_frag_num = fragment_number;

  //instantiate parallel workers
  #pragma omp parallel shared(WF, WB, fragment_number, frag_num_var,\
  cell_number, signal2noise, max_threshold, min_threshold, \
  output_foreground, output_background, Foreground, Background, uniform_sampling) \
  private(chrom_s, start, end, name_s, cell_frag_num, g_seed) \
  default(none)
  {
    #pragma omp for //parallelize the simulation of each cell
    for (unsigned long i = 0; i < cell_number; i++) {

      fast_srand(i * (omp_get_thread_num() + 1));

      // if (i % 1000 == 0) printf("cell %d\n", i);

      cell_frag_num = fragment_number; 
      //sample the number of fragments from a log normal distribution
      if (!uniform_sampling) {
        //create log normal distribution for gaussian fragment number distribution
        default_random_engine generator(fast_rand());
        normal_distribution<long double> distribution(log(fragment_number), frag_num_var); // mean and variance

        cell_frag_num = exp(distribution(generator));
        // printf("%zu\n", cell_frag_num);
        if (cell_frag_num > max_threshold) cell_frag_num = max_threshold;
        if (cell_frag_num < min_threshold) cell_frag_num = min_threshold;
      }

      //sampling foreground fragments twice
      unsigned long peak_round_1 = round((float) cell_frag_num * signal2noise / 2.0);
      unsigned long peak_round_2 = round((float) cell_frag_num * signal2noise) - peak_round_1;

      // printf("%zu %zu %zu\n", cell_frag_num, peak_round_1, peak_round_2);

      vector<int> sampled = sample_int_expj(WF.size(), peak_round_1, WF, fast_rand());
      for (unsigned long j = 0; j < peak_round_1; ++j) {
        tie(chrom_s, start, end, name_s) = Foreground[sampled[j] - 1];
        fprintf(output_foreground, "%s\t%zu\t%zu\t%s\tcell_%d\n", 
          chrom_s.c_str(), start, end, name_s.c_str(), i);
      }

      sampled.clear();
      sampled = sample_int_expj(WF.size(), peak_round_2, WF, fast_rand());
      for (unsigned long j = 0; j < peak_round_2; ++j) {
        tie(chrom_s, start, end, name_s) = Foreground[sampled[j] - 1];
        fprintf(output_foreground, "%s\t%zu\t%zu\t%s\tcell_%d\n", 
          chrom_s.c_str(), start, end, name_s.c_str(), i);
      }

      //if needed, sampling background fragments twice
      unsigned long bg_round_1 = 0;
      unsigned long bg_round_2 = 0;
      if(signal2noise<1.0){

        bg_round_1 = round(cell_frag_num * (1.0 - signal2noise) / 2);
        bg_round_2 = cell_frag_num - peak_round_1 - peak_round_2 - bg_round_1;

        sampled.clear();
        sampled = sample_int_expj(WB.size(), bg_round_1, WB, fast_rand());
        for (unsigned long j = 0; j < sampled.size(); ++j) {
          tie(chrom_s, start, end, name_s) = Background[sampled[j] - 1];
          fprintf(output_background, "%s\t%zu\t%zu\t%s\tcell_%d\n", 
            chrom_s.c_str(), start, end, name_s.c_str(), i);
        }

        sampled.clear();
        sampled = sample_int_expj(WB.size(), bg_round_2, WB, fast_rand());
        for (unsigned long j = 0; j < sampled.size(); ++j) {
          tie(chrom_s, start, end, name_s) = Background[sampled[j] - 1];
          fprintf(output_background, "%s\t%zu\t%zu\t%s\tcell_%d\n", 
            chrom_s.c_str(), start, end, name_s.c_str(), i);
        }
      }
    }
  }


  fclose(output_foreground);
  fclose(output_background);
  exit(0);
}
