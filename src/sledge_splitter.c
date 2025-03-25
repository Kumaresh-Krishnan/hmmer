/* phmmer: search a protein sequence against a protein database
 */
#include <p7_config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_exponential.h"
#include "esl_getopts.h"
#include "esl_gumbel.h"
#include "esl_vectorops.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_random.h"
#include "esl_scorematrix.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"

#ifdef HMMER_MPI
#include "mpi.h"
#include "esl_mpi.h"
#endif

#ifdef HMMER_THREADS
#include <unistd.h>
#include "esl_threads.h"
#include "esl_workqueue.h"
#endif

#include "hmmer.h"
#include "sledge_dev.h"


static ESL_OPTIONS options[] = {
  /* name             type              default   env  range   toggles   reqs   incomp                             help                                       docgroup*/
  { "-h",             eslARG_NONE,        FALSE, NULL, NULL,      NULL,  NULL,  NULL,              "show brief help on version and usage",                         1 },
/* Control of output */
  { "-o",             eslARG_OUTFILE,      NULL, NULL, NULL,      NULL,  NULL,  NULL,              "direct output to file <f>, not stdout",                        2 },
  { "-A",             eslARG_OUTFILE,      NULL, NULL, NULL,      NULL,  NULL,  NULL,              "save multiple alignment of hits to file <f>",                  2 },
  { "--tblout",       eslARG_OUTFILE,      NULL, NULL, NULL,      NULL,  NULL,  NULL,              "save parseable table of per-sequence hits to file <f>",        2 },
  { "--domtblout",    eslARG_OUTFILE,      NULL, NULL, NULL,      NULL,  NULL,  NULL,              "save parseable table of per-domain hits to file <f>",          2 },
  { "--pfamtblout",   eslARG_OUTFILE,      NULL, NULL, NULL,      NULL,  NULL,  NULL,              "save table of hits and domains to file, in Pfam format <f>",   2 },
  { "--acc",          eslARG_NONE,        FALSE, NULL, NULL,      NULL,  NULL,  NULL,              "prefer accessions over names in output",                       2 },
  { "--noali",        eslARG_NONE,        FALSE, NULL, NULL,      NULL,  NULL,  NULL,              "don't output alignments, so output is smaller",                2 },
  { "--notextw",      eslARG_NONE,         NULL, NULL, NULL,      NULL,  NULL, "--textw",          "unlimit ASCII text output line width",                         2 },
  { "--textw",        eslARG_INT,         "120", NULL, "n>=120",  NULL,  NULL, "--notextw",        "set max width of ASCII text output lines",                     2 },
/* Control of scoring system */
  { "--popen",        eslARG_REAL,       "0.02", NULL, "0<=x<0.5",NULL,  NULL,  NULL,              "gap open probability",                                         3 },
  { "--pextend",      eslARG_REAL,        "0.4", NULL, "0<=x<1",  NULL,  NULL,  NULL,              "gap extend probability",                                       3 },
  { "--mx",           eslARG_STRING, "BLOSUM62", NULL, NULL,      NULL,  NULL,  "--mxfile",        "substitution score matrix choice (of some built-in matrices)", 3 },
  { "--mxfile",       eslARG_INFILE,       NULL, NULL, NULL,      NULL,  NULL,  "--mx",            "read substitution score matrix from file <f>",                 3 },
/* Control of reporting thresholds */
  { "-E",             eslARG_REAL,       "10.0", NULL, "x>0",     NULL,  NULL,  REPOPTS,           "report sequences <= this E-value threshold in output",         4 },
  { "-T",             eslARG_REAL,        FALSE, NULL,  NULL,     NULL,  NULL,  REPOPTS,           "report sequences >= this score threshold in output",           4 },
  { "--domE",         eslARG_REAL,       "10.0", NULL, "x>0",     NULL,  NULL,  DOMREPOPTS,        "report domains <= this E-value threshold in output",           4 },
  { "--domT",         eslARG_REAL,        FALSE, NULL,  NULL,     NULL,  NULL,  DOMREPOPTS,        "report domains >= this score cutoff in output",                4 },
/* Control of inclusion thresholds */
  { "--incE",         eslARG_REAL,       "0.01", NULL, "x>0",     NULL,  NULL,  INCOPTS,           "consider sequences <= this E-value threshold as significant",  5 },
  { "--incT",         eslARG_REAL,        FALSE, NULL,  NULL,     NULL,  NULL,  INCOPTS,           "consider sequences >= this score threshold as significant",    5 },
  { "--incdomE",      eslARG_REAL,       "0.01", NULL, "x>0",     NULL,  NULL,  INCDOMOPTS,        "consider domains <= this E-value threshold as significant",    5 },
  { "--incdomT",      eslARG_REAL,        FALSE, NULL,  NULL,     NULL,  NULL,  INCDOMOPTS,        "consider domains >= this score threshold as significant",      5 },
/* Model-specific thresholding for both reporting and inclusion (unused in phmmer)*/
  { "--cut_ga",       eslARG_NONE,        FALSE, NULL, NULL,      NULL,  NULL,  THRESHOPTS,        "use profile's GA gathering cutoffs to set all thresholding",  99 },
  { "--cut_nc",       eslARG_NONE,        FALSE, NULL, NULL,      NULL,  NULL,  THRESHOPTS,        "use profile's NC noise cutoffs to set all thresholding",      99 },
  { "--cut_tc",       eslARG_NONE,        FALSE, NULL, NULL,      NULL,  NULL,  THRESHOPTS,        "use profile's TC trusted cutoffs to set all thresholding",    99 },
/* Control of filter pipeline */
  { "--max",          eslARG_NONE,        FALSE, NULL, NULL,      NULL,  NULL, "--F1,--F2,--F3",   "Turn all heuristic filters off (less speed, more power)",      7 },
  { "--F1",           eslARG_REAL,       "0.02", NULL, NULL,      NULL,  NULL, "--max",            "Stage 1 (MSV) threshold: promote hits w/ P <= F1",             7 },
  { "--F2",           eslARG_REAL,       "1e-3", NULL, NULL,      NULL,  NULL, "--max",            "Stage 2 (Vit) threshold: promote hits w/ P <= F2",             7 },
  { "--F3",           eslARG_REAL,       "1e-5", NULL, NULL,      NULL,  NULL, "--max",            "Stage 3 (Fwd) threshold: promote hits w/ P <= F3",             7 },
  { "--nobias",       eslARG_NONE,        NULL,  NULL, NULL,      NULL,  NULL, "--max",            "turn off composition bias filter",                             7 },
/* Control of E-value calibration */
  { "--EmL",          eslARG_INT,         "200", NULL,"n>0",      NULL,  NULL,  NULL,              "length of sequences for MSV Gumbel mu fit",                   11 },   
  { "--EmN",          eslARG_INT,         "200", NULL,"n>0",      NULL,  NULL,  NULL,              "number of sequences for MSV Gumbel mu fit",                   11 },   
  { "--EvL",          eslARG_INT,         "200", NULL,"n>0",      NULL,  NULL,  NULL,              "length of sequences for Viterbi Gumbel mu fit",               11 },   
  { "--EvN",          eslARG_INT,         "200", NULL,"n>0",      NULL,  NULL,  NULL,              "number of sequences for Viterbi Gumbel mu fit",               11 },   
  { "--EfL",          eslARG_INT,         "100", NULL,"n>0",      NULL,  NULL,  NULL,              "length of sequences for Forward exp tail tau fit",            11 },   
  { "--EfN",          eslARG_INT,         "200", NULL,"n>0",      NULL,  NULL,  NULL,              "number of sequences for Forward exp tail tau fit",            11 },   
  { "--Eft",          eslARG_REAL,       "0.04", NULL,"0<x<1",    NULL,  NULL,  NULL,              "tail mass for Forward exponential tail tau fit",              11 },   
/* other options */
  { "--nonull2",      eslARG_NONE,        NULL,  NULL, NULL,      NULL,  NULL,  NULL,              "turn off biased composition score corrections",               12 },
  { "-Z",             eslARG_REAL,       FALSE,  NULL, "x>0",     NULL,  NULL,  NULL,              "set # of comparisons done, for E-value calculation",          12 },
  { "--domZ",         eslARG_REAL,       FALSE,  NULL, "x>0",     NULL,  NULL,  NULL,              "set # of significant seqs, for domain E-value calculation",   12 },
  { "--seed",         eslARG_INT,         "42",  NULL, "n>=0",    NULL,  NULL,  NULL,              "set RNG seed to <n> (if 0: one-time arbitrary seed)",         12 },
  { "--qformat",      eslARG_STRING,      NULL,  NULL, NULL,      NULL,  NULL,  NULL,              "assert query <seqfile> is in format <s>: no autodetection",   12 },
  { "--tformat",      eslARG_STRING,      NULL,  NULL, NULL,      NULL,  NULL,  NULL,              "assert target <seqdb> is in format <s>>: no autodetection",   12 },
/* Sledgehmmer specific options */  
  { "--qblock",       eslARG_INT,        "10",   NULL, NULL,      NULL,  NULL,  NULL,              "max queries handled at a time",                               12 },
  { "--dbblock",      eslARG_INT,        "1000", NULL, NULL,      NULL,  NULL,  NULL,              "max targets handled at a time",                               12 },
  { "--filter",       eslARG_NONE,       FALSE,  NULL, NULL,      NULL,  NULL,  NULL,              "filter mode toggle",                                          12 },
  { "--halt",         eslARG_INT,        "-1",   NULL, NULL,      NULL,  NULL,  NULL,              "stop after n sequences - for debugging",                         12 },
  { "--train_frac",   eslARG_REAL,        "0.9",   NULL, NULL,      NULL,  NULL,  NULL,              "train proportion after test limit reached",                    12 },
  { "--test_limit",   eslARG_INT,        "75000",   NULL, NULL,      NULL,  NULL,  NULL,              "required minimum number of test sequences",                  12 },
  { "--init_chunk",   eslARG_INT,        "20",   NULL, NULL,      NULL,  NULL,  NULL,              "query chunk size till test limit",                              12 },
  { "--train_only",   eslARG_NONE,        FALSE,   NULL, NULL,      NULL,  NULL,  NULL,              "growing only train set",                                     12 },
  { "--test_only",    eslARG_NONE,        FALSE,   NULL, NULL,      NULL,  NULL,  NULL,              "growing only test set",                                      12 },
  { "--suppress",     eslARG_NONE,        FALSE,   NULL, NULL,      NULL,  NULL,  NULL,              "turn off progress bar",                                      12 },
  { "--task_id",      eslARG_INT,        "0",   NULL, NULL,      NULL,  NULL,  NULL,              "slurm array task ID or shard number",                          12 },
  { "--plow",         eslARG_REAL,       "0.00", NULL, NULL,      NULL,  NULL,  NULL,              "pid lower limit for accepting a sequence",                    12 },
  { "--phigh",        eslARG_REAL,       "1.00", NULL, NULL,      NULL,  NULL,  NULL,              "pid upper limit for accepting a sequence",                    12 },
  { "--qsize",        eslARG_INT,        "1",    NULL, NULL,      NULL,  NULL,  NULL,              "num queries per thread",                                      12 },
  { "--train_path",   eslARG_STRING,        "../results/train_seqs",    NULL, NULL,      NULL,  NULL,  NULL,              "filename for train seqs with path",     12 },
  { "--test_path",    eslARG_STRING,        "../results/test_seqs",    NULL, NULL,      NULL,  NULL,  NULL,              "filename for test seqs with path",       12 },
  { "--discard_path", eslARG_STRING,        "../results/discard_seqs",    NULL, NULL,      NULL,  NULL,  NULL,              "filename for discard seqs with path", 12 },
  { "--load_tr",    eslARG_STRING,        "-",    NULL, NULL,      NULL,  NULL,  NULL,              "Resume from this train db", 12 },
  { "--load_te",    eslARG_STRING,        "-",    NULL, NULL,      NULL,  NULL,  NULL,              "Resume from this test db",  12 },
#ifdef HMMER_THREADS
  { "--cpu",        eslARG_INT,        "1","HMMER_NCPU", "n>=0",NULL,  NULL,  CPUOPTS,           "number of CPU cores to use",                                  12 },
#endif
#ifdef HMMER_MPI
  { "--stall",      eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--mpi", NULL,                   "arrest after start: for debugging MPI under gdb",             12 },  
  { "--mpi",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  MPIOPTS,                "run as an MPI parallel program",                              12 },
#endif
 {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <seqdb>";
static char banner[] = "split a database into train and test set";

/* process_commandline()
 * Take argc, argv, and options; parse the command line;
 * display help/usage info.
 */
static int
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_qfile, char **ret_dbfile)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  int          status;

  if (esl_opt_ProcessEnvironment(go)         != eslOK)  { if (printf("Failed to process environment: %s\n", go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { if (printf("Failed to parse command line: %s\n",  go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_VerifyConfig(go)               != eslOK)  { if (printf("Failed to parse command line: %s\n",  go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      p7_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      if (puts("\nBasic options:")                                           < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 120=textwidth*/

      if (puts("\nOptions directing output:")                                < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 

      if (puts("\nOptions controlling scoring system:")                      < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80); 

      if (puts("\nOptions controlling reporting thresholds:")                < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 

      if (puts("\nOptions controlling inclusion (significance) thresholds:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80); 

      if (puts("\nOptions controlling acceleration heuristics:")             < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80); 

      if (puts("\nOptions controlling E value calibration:")                 < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 11, 2, 80); 

      if (puts("\nOther expert options:")                                    < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 12, 2, 80); 
      exit(0);
    }

  if (esl_opt_ArgNumber(go)                 != 1)    { if (puts("Incorrect number of command line arguments.")      < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((*ret_dbfile = esl_opt_GetArg(go, 1)) == NULL) { if (puts("Failed to get <seqdb> argument on command line")   < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  /* Validate any attempted use of stdin streams */
  if (strcmp(*ret_dbfile, "-") == 0) 
    { if (puts("Either <seqfile> or <seqdb> may be '-' (to read from stdin), but not both.") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");  goto FAILURE; }

  *ret_go = go;
  return eslOK;
  
 FAILURE:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, argv[0], usage);
  if (puts("\nwhere basic options are:")                                       < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 120); /* 1= group; 2 = indentation; 120=textwidth*/
  if (printf("\nTo see more help on available options, do %s -h\n\n", argv[0]) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  esl_getopts_Destroy(go);

  exit(1);

 ERROR:
  if (go) esl_getopts_Destroy(go);
  exit(status);
}

int
main(int argc, char **argv)
{ 
  int status = eslOK;

  ESL_GETOPTS     *go               = NULL;	              /* command line processing    */
  ESL_SQFILE      *dbfp             = NULL;               /* open dbfile                */
  ESL_SQ_BLOCK    *train_db         = NULL;               /* train database             */
  ESL_SQ_BLOCK    *test_db          = NULL;               /* test database              */
  ESL_SQ_BLOCK    *qdb              = NULL;               /* query database             */
  ESL_SQ_BLOCK    *sqdb             = NULL;               /* sequence database to split */
  int              dbformat         = eslSQFILE_UNKNOWN;  /* format of dbfile           */
  ESL_ALPHABET    *abc              = NULL;               /* sequence alphabet          */
  char            *results          = NULL;               /* results of seq inclusion   */
  struct cfg_s     cfg;                                   /* configuration data         */
  SLEDGE_INFO      si;                                    /* Sledgehmmer info           */
  int              num_records;                           /* Total sequences in db      */
  ESL_RANDOMNESS  *rnd;                                   /* Random number generator    */
  int              seq_ctr;                               /* Count sequences processed  */
  int              tmp_ctr;                               /* Temp counter               */
  int              max_qsize;                             /* Queries processed per core */
  int              max_cores;                             /* Number of cores to use     */
  bool             assign_train;                          /* Flag for train assignment  */
  bool             low;                                   /* Flag for low throughput    */
  bool             write_train;                           /* Flag to save train seqs    */
  bool             write_test;                            /* Flag to save test seqs     */
  bool             suppress;                              /* Flag to suppress outputs   */
  int              Q_SIZE;                                /* Query chunk size           */
  double           current_train_frac;                    /* Current train proportion   */
  char            *train_file;                            /* Train sequence file name   */
  char            *test_file;                             /* Test sequence file name    */
  char            *discard_file;                          /* Discard sequence file name */
  FILE            *train_fp;                              /* Train sequence file ptr    */
  FILE            *test_fp;                               /* Test sequence file ptr     */
  FILE            *discard_fp;                            /* Discard sequence file ptr  */
  int              len;                                   /* Calculate file name length */
  clock_t          start_time;                            /* Keep track of start time   */
  clock_t          elapsed_time;                          /* Elapsed time               */

  /* Set processor specific flags */
  impl_Init();

  /* Initialize what we can in the config structure (without knowing the alphabet yet) */
  cfg.qfile        = NULL;
  cfg.dbfile       = NULL;
  cfg.do_mpi       = FALSE;	         /* this gets reset below, if we init MPI */
  cfg.nproc        = 0;		           /* this gets reset below, if we init MPI */
  cfg.my_rank      = 0;		           /* this gets reset below, if we init MPI */
  cfg.firstseq_key = NULL;
  cfg.n_targetseq  = -1;

  /* Initializations */

  p7_FLogsumInit();		/* we're going to use table-driven Logsum() approximations at times */
  process_commandline(argc, argv, &go, &cfg.qfile, &cfg.dbfile);

  /* Initialize parameters for Sledgehmmer*/
  si.pid_low      = esl_opt_GetReal    (go, "--plow");
  si.pid_high     = esl_opt_GetReal    (go, "--phigh");
  si.filter       = esl_opt_GetBoolean (go, "--filter");
  si.qsize        = esl_opt_GetInteger (go, "--qsize");
  si.cores        = esl_opt_GetInteger (go, "--cpu");
  si.incE         = esl_opt_GetReal    (go, "--incE");
  si.db_size      = esl_opt_GetReal    (go, "-Z");
  si.train_only   = esl_opt_GetBoolean (go, "--train_only");
  si.test_only    = esl_opt_GetBoolean (go, "--test_only");
  si.init_chunk   = esl_opt_GetInteger (go, "--init_chunk");
  si.train_frac   = esl_opt_GetReal    (go, "--train_frac");
  si.test_limit   = esl_opt_GetInteger (go, "--test_limit");
  si.task_id      = esl_opt_GetInteger (go, "--task_id");
  si.train_path   = esl_opt_GetString  (go, "--train_path");
  si.test_path    = esl_opt_GetString  (go, "--test_path");
  si.discard_path = esl_opt_GetString  (go, "--discard_path");
  si.out_path     = NULL; /* This is used only for the filtering code */
  si.out_fp       = NULL; /* This is used only for the filtering code */
  suppress        = esl_opt_GetBoolean (go, "--suppress"); /* Store suppress flag */

  max_cores = si.cores; /* Store max possible cores for high throughput */
  max_qsize = si.qsize; /* Store max possible qsize for high throughput */

  /* If caller declared input formats, decode them */
  if (esl_opt_IsOn(go, "--tformat"))
  {
    dbformat = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--tformat"));
    if (dbformat == eslSQFILE_UNKNOWN) p7_Fail("%s is not a recognized sequence database file format\n", esl_opt_GetString(go, "--tformat"));
  }
  
  /* Initialize alphabet */
  abc = esl_alphabet_Create(eslAMINO);

  /* Open the sequence database file to split */
  status =  esl_sqfile_OpenDigital(abc, cfg.dbfile, dbformat, p7_SEQDBENV, &dbfp);
  if      (status == eslENOTFOUND) p7_Fail("Failed to open target sequence database %s for reading\n",      cfg.dbfile);
  else if (status == eslEFORMAT)   p7_Fail("Target sequence database file %s is empty or misformatted\n",   cfg.dbfile);
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        p7_Fail("Unexpected error %d opening target sequence database file %s\n", status, cfg.dbfile);
  
  sqdb = esl_sq_CreateDigitalBlock(esl_opt_GetInteger(go, "--dbblock"), abc);
  if (sqdb == NULL) p7_Fail("Failed to allocate sequence database block");

  status = read_cust(dbfp, sqdb);
  if (status != eslOK) p7_Fail("Failed to read target sequences from database");
  esl_sqfile_Close(dbfp);

  /* Calculate total number of sequences to process */
  num_records = sqdb->count;
  if (esl_opt_GetInteger(go, "--halt") > 0) num_records = esl_opt_GetInteger(go, "--halt");

  /* Allocate memory for the train SQ_BLOCK */
  train_db = esl_sq_CreateDigitalBlock((int) (si.train_frac * sqdb->count), abc);
  if (train_db == NULL) p7_Fail("Failed to allocate train database block");

  /* Check if we are loading an existing train db initially */
  if (strcmp(esl_opt_GetString(go, "--load_tr"), "-") != 0)
  {
    char *fname = NULL;
    fname = esl_opt_GetString(go, "--load_tr");
    ESL_SQFILE *fp = NULL;
    status =  esl_sqfile_OpenDigital(abc, fname, dbformat, p7_SEQDBENV, &fp);
    if      (status == eslENOTFOUND) p7_Fail("Failed to open train sequence database %s for reading\n",      fname);
    else if (status == eslEFORMAT)   p7_Fail("Train sequence database file %s is empty or misformatted\n",   fname);
    else if (status == eslEINVAL)    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
    else if (status != eslOK)        p7_Fail("Unexpected error %d opening train sequence database file %s\n", status, fname);

    read_cust(fp, train_db);
    if (status != eslOK) p7_Fail("Failed to read train sequences from database");
      esl_sqfile_Close(fp);
    free(fp);
  }

  /* Allocate memory for the test SQ_BLOCK */
  test_db = esl_sq_CreateDigitalBlock((int) ((1.0-si.train_frac) * sqdb->count), abc);
  if (test_db == NULL) p7_Fail("Failed to allocate test database block");

  /* Check if we are loading an existing test db initially */
  if (strcmp(esl_opt_GetString(go, "--load_te"), "-") != 0)
  {
    char *fname = NULL;
    fname = esl_opt_GetString(go, "--load_te");
    ESL_SQFILE *fp = NULL;
    status =  esl_sqfile_OpenDigital(abc, fname, dbformat, p7_SEQDBENV, &fp);
    if      (status == eslENOTFOUND) p7_Fail("Failed to open test sequence database %s for reading\n",      fname);
    else if (status == eslEFORMAT)   p7_Fail("Test sequence database file %s is empty or misformatted\n",   fname);
    else if (status == eslEINVAL)    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
    else if (status != eslOK)        p7_Fail("Unexpected error %d opening test sequence database file %s\n", status, fname);

    read_cust(fp, test_db);
    if (status != eslOK) p7_Fail("Failed to read test sequences from database");
    esl_sqfile_Close(fp);
  }

  /* Allocate memory for the query SQ_BLOCK */
  qdb = esl_sq_CreateDigitalBlock(si.cores * si.qsize, abc); /* Allocate for max possible, not just Q_SIZE */
  if (qdb == NULL) p7_Fail("Failed to allocate query database block");

  if (si.test_limit == -1) si.test_limit = test_db->listSize; /* If test limit not defined then decide based on train fraction */
  
  /* Set initial cores and query chunk sizes */
  Q_SIZE = esl_opt_GetInteger(go, "--init_chunk");  /* Initial query set size        */
  low = TRUE;                                       /* Low throughput mode initially */
  current_train_frac = 0.5;                         /* 50-50 split to start with     */

  if ((si.train_only) || (si.test_only))
  {
    Q_SIZE = si.qsize * si.cores;
    low = FALSE;
    write_train = (si.test_only) ? FALSE : TRUE;
    write_test = (si.train_only) ? FALSE : TRUE;
  }

  else if (si.cores < Q_SIZE)
  {
    if (Q_SIZE % si.cores == 0)
    {
      si.qsize = Q_SIZE / si.cores;
      si.cores;
    }

    else
    {
      printf("Since cores < init chunk, they must be exactly divisible\n");
      exit(-20);
    }
  }

  else
  {
    si.qsize = 1;       /* In low throughput mode, qsize 1 and cores as needed */
    si.cores = Q_SIZE;  /* Cores as much as needed for Q_SIZE                  */
  }

  /* Allocate max possible memory for results + 1 for \0 */
  ESL_ALLOC(results, (qdb->listSize+1)*sizeof(char));

  /* Create seeded random generator for the rest of the code to use */
  rnd = esl_randomness_Create(42);
  
  /* Open the discard file for permanent writing */

  len = snprintf(NULL, 0, "%s_%d.fasta", si.discard_path, si.task_id) + 1;
  ESL_ALLOC(discard_file, len*sizeof(char));
  
  snprintf(discard_file, len, "%s_%d.fasta", si.discard_path, si.task_id);
  if ((discard_fp = fopen(discard_file, "w")) == NULL) p7_Fail("Failed to open discard file %s for writing\n", discard_file);

  /* Loop through all sequences in database */
  start_time = get_coarse_time();
  
  for (seq_ctr = 0; seq_ctr < num_records; seq_ctr++)
  {
    if (seq_ctr < 0) continue; /* Useful for resuming the code from a certain point */
    
    else
    {
      add_seq(qdb, sqdb->list + seq_ctr, abc);

      if (((qdb->count == Q_SIZE) || (seq_ctr == num_records-1)) && ((si.train_only) || (si.test_only)))
      {
        if (si.train_only)
          status = assign_master(go, qdb, test_db, &si, results);
        else
          status = assign_master(go, qdb, train_db, &si, results);

        for (tmp_ctr=0; tmp_ctr<qdb->count; tmp_ctr++)
        {
          if (results[tmp_ctr] == '1')
          {
            if (si.train_only) add_seq(train_db, qdb->list+tmp_ctr, abc);
            else add_seq(test_db, qdb->list+tmp_ctr, abc);
          }
          else
            esl_sqio_Write(discard_fp, qdb->list+tmp_ctr, eslSQFILE_FASTA, FALSE);  
        }
      }

      else if ((qdb->count== Q_SIZE) || (seq_ctr == num_records-1))
      {
        if (esl_random(rnd) < current_train_frac) assign_train = TRUE;
        else assign_train = FALSE;

        if (((assign_train) && (test_db->count == 0)) || ((!assign_train) && (train_db->count == 0)))
        {
          memset(results, '1', Q_SIZE*sizeof(char));
          results[Q_SIZE] = '\0';
        }

        else if (assign_train)
          status = assign_master(go, qdb, test_db, &si, results);

        else
          status = assign_master(go, qdb, train_db, &si, results);

        ESL_SQ_BLOCK *cdb = esl_sq_CreateDigitalBlock(Q_SIZE, abc); /* Allocate similar to qdb for checking other side */
        if (cdb == NULL) p7_Fail("Failed to allocate candidate database block");
        
        for (tmp_ctr=0; tmp_ctr<qdb->count; tmp_ctr++)
        {
          if (results[tmp_ctr] == '1')
          {
            if (assign_train) add_seq(train_db, qdb->list+tmp_ctr, abc);
            else add_seq(test_db, qdb->list+tmp_ctr, abc);
          }
          else
            add_seq(cdb, qdb->list+tmp_ctr, abc);
        }

        /* Process candidates for other side */

        if (cdb->count > 0)
        {
          int old_cores = si.cores;                 /* Must revert to this after checking other side */
          int old_qsize = si.qsize;                 /* Must revert to this after checking other side */

          if (si.cores > cdb->count)
          {
            si.cores = cdb->count;                  /* Only as many cores as candidate seqs          */
            si.qsize = 1;                           /* One query per core for candidate seqs         */
          }
          else
          {
            si.qsize = ceil(cdb->count / si.cores); /* Too many candidates, assign qsize accordingly */
            si.cores = ceil(cdb->count / si.qsize); /* Now adjust number of cores based on qsize     */
          }

          if (assign_train)
            status = assign_master(go, cdb, train_db, &si, results);
          else
            status = assign_master(go, cdb, test_db, &si, results);
          
          for (tmp_ctr=0; tmp_ctr<cdb->count; tmp_ctr++)
          {
            if (results[tmp_ctr] == '1')
            {
              if (assign_train) add_seq(test_db, cdb->list+tmp_ctr, abc);
              else add_seq(train_db, cdb->list+tmp_ctr, abc);
            }
            else
              esl_sqio_Write(discard_fp, cdb->list+tmp_ctr, eslSQFILE_FASTA, FALSE); 
          }
          si.cores = old_cores;                   /* Restore cores to prior value                   */
          si.qsize = old_qsize;                   /* Restore qsize to prior value                   */
        }

        esl_sq_DestroyBlock(cdb);                 /* Free the candidate db created here */
      }
      
      if (qdb->count == Q_SIZE)
      {
        qdb->count = 0;

        if (num_records - seq_ctr - 1 < Q_SIZE)
        {
          Q_SIZE = num_records - seq_ctr - 1;

          if (si.cores >= Q_SIZE)
          {
            si.qsize = 1;                         /* Last chunk, use only required qsize           */
            si.cores = Q_SIZE;                    /* Last chunk, use only required cores           */
          }

          else
          {
            si.qsize = ceil(Q_SIZE / si.cores);   /* Core limited, so choose qsize wisely          */
            si.cores = ceil(Q_SIZE / si.qsize);   /* New cores based on adjusted qsize             */
          }
        }
      }

      if ((low) && (test_db->count >= si.test_limit))
      {
        si.cores = max_cores;                     /* High throughput mode with available cores     */
        si.qsize = max_qsize;                     /* High throughput with max qsize                */
        Q_SIZE = si.cores * si.qsize;             /* Recompute queries processed at a time         */
        current_train_frac = si.train_frac;       /* Switch to desired train proportion            */
        si.train_only = 1;                        /* Assign only to train set since limit is met   */
        low = FALSE;                              /* Turn off low throughput flag                  */
      }
    }

    if ((seq_ctr % 1000 == 0) && !suppress) /* How often do we want to update the progress bar? It has to be db dependent*/
    {
      elapsed_time = get_coarse_time() - start_time;
      print_progress(seq_ctr+1, num_records, elapsed_time);

      // if (seq_ctr %1000 == 0) fflush(stdout);
      fflush(stderr);
    }
  }

  /* Print progress for one final time */
  if (!suppress)
  {
    elapsed_time = get_coarse_time() - start_time;;
    print_progress(seq_ctr, num_records, elapsed_time);
  }

  /* Writing train sequences to a file */
  if (write_train)
  {
    len = snprintf(NULL, 0, "%s_%d.fasta", si.train_path, si.task_id) + 1;
    ESL_ALLOC(train_file, len*sizeof(char));
    snprintf(train_file, len, "%s_%d.fasta", si.train_path, si.task_id);
    if ((train_fp = fopen(train_file, "w")) == NULL) p7_Fail("Failed to open train file %s for writing\n", train_file);
    
    for(tmp_ctr=0; tmp_ctr < train_db->count; tmp_ctr++)
      esl_sqio_Write(train_fp, train_db->list+tmp_ctr, eslSQFILE_FASTA, FALSE);
  }

  /* Writing test sequences to a file */
  if (write_test)
  {
    len = snprintf(NULL, 0, "%s_%d.fasta", si.test_path, si.task_id) + 1;
    ESL_ALLOC(test_file, len*sizeof(char));
    snprintf(test_file, len, "%s_%d.fasta", si.test_path, si.task_id);
    if ((test_fp = fopen(test_file, "w")) == NULL) p7_Fail("Failed to open test file %s for writing\n", test_file);

    for(tmp_ctr=0; tmp_ctr < test_db->count; tmp_ctr++)
      esl_sqio_Write(test_fp, test_db->list+tmp_ctr, eslSQFILE_FASTA, FALSE);
  }

  /* Close all files */

  if (write_train) fclose(train_fp);
  if (write_test) fclose(test_fp);
  fclose(discard_fp);

  /* Mop the floor, empty the trash, free the souls and leave */

  if (write_train) free(train_file);
  if (write_test) free(test_file);
  free(discard_file);
  free(results);

  esl_sq_DestroyBlock(train_db);
  esl_sq_DestroyBlock(test_db);
  esl_sq_DestroyBlock(qdb);
  esl_sq_DestroyBlock(sqdb);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(rnd);

  ERROR:
    exit(status);

  exit(status);
}
