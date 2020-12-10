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
  //this script assumes: 
  //the input is in standard BED file format

  if (argc != 4) {
    printf("please input ./uniform_sampling <peak_read_intersect.bed> <sampled_peaks.bed> <output>\n");
    exit(1);
  }

  char chrom[6];
  string chrom_s;
  unsigned long start;
  unsigned long end;
  char name[30];
  string name_s;

  char pre_name[30] = "p_0";

  char rchrom[6];
  string rchrom_s;
  unsigned long rstart;
  unsigned long rend;

  char id[30];
  string id_s;

  char cell_name[20];
  string cell_name_s;

  unsigned long score;
  char strand[1];
  string strand_s;

  // declare a map called Peak_Reads
  map <
    tuple<string, unsigned long, unsigned long>,
    vector<tuple<string, unsigned long, unsigned long, string, unsigned long, string>>
  > Peak_Reads;

  //read in reads files
  FILE* reads_file = fopen(argv[1], "r+");
  if (reads_file == NULL) {
    fprintf(stderr, "fail to open file\n");
    exit(1);
  }
  // printf("file opened\n");

  int got = 10;
  do {
    got = fscanf(reads_file, "%s %zu %zu %s %s %zu %zu %s %zu %s", 
      &chrom, &start, &end, &name,
      &rchrom, &rstart, &rend, &id,
      &score, &strand);


    chrom_s = chrom;
    name_s = name;
    rchrom_s = rchrom;
    id_s = id;
    strand_s = strand;

    if (strcmp(pre_name, name) != 0) {
      vector<tuple<string, unsigned long, unsigned long, string, unsigned long, string>> vect;
      vect.push_back(make_tuple(rchrom_s, rstart, rend, id_s, score, strand_s));
      Peak_Reads[make_tuple(chrom_s, start, end)] = vect;
      strcpy(pre_name, name);

    } else {
      Peak_Reads[make_tuple(chrom_s, start, end)].push_back(make_tuple(rchrom_s, rstart, rend, id_s, score, strand_s));
    }
  }
  while (got == 10); //iterate until EOF
  fclose(reads_file);

  // printf("peak-reads dictionary constructed\n");

  FILE* sampled_file = fopen(argv[2], "r+");
  if (sampled_file == NULL) {
    fprintf(stderr, "fail to open file\n");
    exit(1);
  }

  got = 5;
  vector <tuple<string, unsigned long, unsigned long, string, string>> sampled_peaks;
  do{
    //not thread safe
    got = fscanf(sampled_file, "%s\t%zu\t%zu\t%s\t%s", 
      &chrom, &start, &end, &name, &cell_name);
    chrom_s = chrom;
    name_s = name;
    cell_name_s = cell_name;

    sampled_peaks.push_back(make_tuple(chrom_s, start, end, name_s, cell_name_s));
  } while (got == 5);

  FILE *fp=fopen(argv[3],"w");
  if (fp == NULL) {
    fprintf(stderr, "fail to open file\n");
    exit(1);
  }

  unsigned long rand_choice = 0;
  unsigned long peak_read_size;
  unsigned long peak_file_len = sampled_peaks.size()-1;
  // printf("all %zu sampled peaks IO'ed\n", peak_file_len-1);
  
  #pragma omp parallel for \
  shared(peak_file_len, sampled_peaks, Peak_Reads, fp) \
  private(chrom_s, start, end, name_s, cell_name_s, rand_choice, rchrom_s, rstart, rend, id, id_s, score, strand_s, peak_read_size) \
  default(none)
  for (unsigned long i = 0; i < peak_file_len; i++) {

    tie(chrom_s, start, end, name_s, cell_name_s) = sampled_peaks[i];
    peak_read_size = Peak_Reads[make_tuple(chrom_s, start, end)].size();

    if (peak_read_size != 0) {
      rand_choice = fast_rand() % peak_read_size;
    } else {
      printf("sampled peaks with no reads!\n");
      printf("%s %zu %zu %zu\n", chrom_s.c_str(), start, end, peak_read_size);
      // cout << chrom_s << " " << start << " " << end << " " << peak_read_size << "\n";
      continue;
    }
    tie(rchrom_s, rstart, rend, id_s, score, strand_s) = Peak_Reads[make_tuple(chrom_s, start, end)][rand_choice];

    fprintf(fp, "%s\t%zu\t%zu\t%s\t%s\t%zu\t%zu\t%s\t%zu\t%s\t%s\n", 
      chrom_s.c_str(), start, end, name_s.c_str(),
      rchrom_s.c_str(), rstart, rend, id_s.c_str(), 
      score, strand_s.c_str(), cell_name_s.c_str());
    
  }
  // while (got == 5); //iterate until EOF
  fclose(sampled_file);
  fclose(fp);

  return 0;

}
