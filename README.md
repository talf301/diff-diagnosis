diff-diagnosis
==============

Improving phenotips' differential diagnosis tool

These scripts all take a directory of patient phenotype descriptions, and write out suggestions to a directory, running some sort of differential diagnosis algorithm on those patients.

diagnose_phenomizer_boring_pval.py: finds differential diagnoses for every patient in directory, uses ancestor overlap and information content
  $1 -> directory containing patients. Each patient is a comma-separated list of HPO terms.
  $2 -> file containing information content for each HPO term. Two tab-separated columns: HPO term on left, IC on right.
  $3 -> out directory. Files are 20 lines, the top disease suggestions. They are tab-separated: p-value, score, disease code.
  $4 -> directory containing p-values. Contains one subfolder for each disease. Each disease subfolder contains 10 files, for each query size from 1 to 10. Each file has one float per line, smallest number at the top.
  $5 -> normalization option. 'True' runs like simgic with p-values, 'half' will normalize by query information content, 'False' does not normalize.
  
diagnose_phenomizer_boring.py: same as above, without p-values, and no normalization.
  $1, $2, $3 -> in dir, ic file, out dir, same as above.

diagnose_phenomizer_simgic_pval.py: finds differential diagnoses for every patient in directory, uses simgic score - information content of ancestor overlap, intersection over union
  $1, $2, $3, $4 -> in dir, ic file, out dir, pval dir, same as above.
  
diagnose_phenomizer_simgic.py: same as above, without p-values.
  $1, $2, $3 -> in dir, ic file, out dir, same as above.

diagnose_phenomizer_mica.py: finds differential diagnoses for every patient in directory, uses most informative common ancestor algorithm (a la Phenomizer)
  $1 -> in directory, as above
  $2 -> out directory, as above

The following scripts are used to make the p-value distributions.

make_distribution.py: creates simgic distribution for a given disease.
  $1 -> disease code.
  $2 -> out dir. Writes 10 files, one for each query size from 1-10, with 100000 sorted random observations for each.
  $3 -> ic file, as above.

make_distribution_boring.py: creates unnormalized simgic distribution for given disease.
  $1, $2, $3 -> same as above.

dispatch_pval_simgic.sh, dispatch_pval_boring.sh: dispatches jobs across cluster to create above distributions in parallel.
  $1 -> determines which ic file to use.

The following scripts are used to get results from Phenomizer or Phenotips.

phenomizer_get.sh: uses wget queries to get diagnosis results from Phenomizer at compbio.charite.de.
phenotips_get.sh: uses wget queries to get diagnosis results from Phenotips at playground.phenotips.org.
phenotips_get_diff.sh: uses wget queries to get diff-diagnosis phenotype suggestions from Phenotips at playground.phenotips.org.

parse_phenotips.py, parse_phenotips_diff.py: parse the returns from their respective shell scripts. Take $1, $2 -> in dir, out dir.

The following scripts score a results directory from a disease suggestion script.

score_phenomizer.py: Creates score.log in results directory, containing info of how many top 1, top 5, and top 20 hits this directory contains.
  $1 -> in dir. Filenames are formatted 'orphanetnumber_indexnumber_hpo.txt'.
  
score_results.sh: Prints out distribution of hits from 1 to 20, and the number of misses.
  $1 -> in dir. There is a diagnosis file which maps patient names to correct diagnoses.

These files are not run to create any results, but they contain functions used by other scripts.

get_hp_ic.py
hpo_lib.py: has many HPO manipulation functions.
hpo.py: loads HPO, has many HPO manipulation functions.
omim.py: does parsing of phenotype annotation file.
orpha.py: does mapping between Orphanet and OMIM.
