#include <iostream>
#include "common.h"
#include <dirent.h>

/*
 *  Recompute the full clvs and likelihood for the given tree and sequences
 *  iterations times, and print the elapsed time in ms
 */
void partitioned_full_traversal(int argc, char *params[])
{
  if (argc != 6) {
    std::cerr << "Error : syntax is" << std::endl;
    std::cerr << "newick partitions_dir use_repeats update_repeats iterations arch" << std::endl;
    return ;
  }
  unsigned int i = 0;
  const char *newick = params[i++];
  const char *lbdir = params[i++];
  unsigned int use_repeats = atoi(params[i++]);
  unsigned int update_repeats = atoi(params[i++]);
  unsigned int iterations = atoi(params[i++]);
  const char *arch = params[i++];


  unsigned int attribute = PLLHelper::compute_attribute(use_repeats, arch);
  if (INVALID_ATTRIBUTE == attribute) {
    return;
  }

  std::vector<PLLHelper *> helpers;


  char buf[500];
  DIR *mydir;
  struct dirent *myfile;
  mydir = opendir(lbdir);
  while((myfile = readdir(mydir)) != NULL)
  {
    const char * seqdir = myfile->d_name;
    if (seqdir[0] == 'c') {
      sprintf(buf, "%s/%s", lbdir, myfile->d_name);
      char buf2[500];
      DIR *mydir2;
      struct dirent *myfile2;
      mydir2 = opendir(buf);
      while((myfile2 = readdir(mydir2)) != NULL)
      {
        const char * seq = myfile2->d_name;
        if (seq[strlen(seq) - 1] == 'y') {
          sprintf(buf2, "%s/%s", buf, myfile2->d_name);
          std::cout << buf2 << std::endl;
          PLLHelper *helper = new PLLHelper(newick, buf2, attribute);
          helper->update_all_partials();
          helpers.push_back(helper);
        }
      }
      closedir(mydir);
    }
  }
  Timer t;
  for (i = 0; i < iterations; ++i) {
    for (unsigned int j = 0; j < helpers.size(); ++j) {
      helpers[j]->update_all_partials(update_repeats);
      helpers[j]->get_likelihood();
    }
  }
  std::cout << t.get_time() << "ms" << std::endl; 
  for (i = 0; i < helpers.size(); ++i) {
    delete helpers[i];
  }
}


