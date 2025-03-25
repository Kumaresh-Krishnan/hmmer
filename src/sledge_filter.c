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
  { "-h",             eslARG_NONE,        FALSE, NULL, NULL,        NULL,  NULL,  NULL,              "show brief help on version and usage",                         1 },
/* Control of output */
  { "-A",             eslARG_OUTFILE,      NULL, NULL, NULL,        NULL,  NULL,  NULL,              "save multiple alignment of hits to file <f>",                  2 },
  { "--tblout",       eslARG_OUTFILE,      NULL, NULL, NULL,        NULL,  NULL,  NULL,              "save parseable table of per-sequence hits to file <f>",        2 },
  { "--domtblout",    eslARG_OUTFILE,      NULL, NULL, NULL,        NULL,  NULL,  NULL,              "save parseable table of per-domain hits to file <f>",          2 },
  { "--pfamtblout",   eslARG_OUTFILE,      NULL, NULL, NULL,        NULL,  NULL,  NULL,              "save table of hits and domains to file, in Pfam format <f>",   2 },
  { "--acc",          eslARG_NONE,        FALSE, NULL, NULL,        NULL,  NULL,  NULL,              "prefer accessions over names in output",                       2 },
  { "--noali",        eslARG_NONE,        FALSE, NULL, NULL,        NULL,  NULL,  NULL,              "don't output alignments, so output is smaller",                2 },
  { "--notextw",      eslARG_NONE,         NULL, NULL, NULL,        NULL,  NULL, "--textw",          "unlimit ASCII text output line width",                         2 },
  { "--textw",        eslARG_INT,         "120", NULL, "n>=120",    NULL,  NULL, "--notextw",        "set max width of ASCII text output lines",                     2 },
/* Control of scoring system */
  { "--popen",        eslARG_REAL,       "0.02", NULL, "0<=x<0.5",  NULL,  NULL,  NULL,              "gap open probability",                                         3 },
  { "--pextend",      eslARG_REAL,        "0.4", NULL, "0<=x<1",    NULL,  NULL,  NULL,              "gap extend probability",                                       3 },
  { "--mx",           eslARG_STRING, "BLOSUM62", NULL, NULL,        NULL,  NULL,  "--mxfile",        "substitution score matrix choice (of some built-in matrices)", 3 },
  { "--mxfile",       eslARG_INFILE,       NULL, NULL, NULL,        NULL,  NULL,  "--mx",            "read substitution score matrix from file <f>",                 3 },
/* Control of reporting thresholds */
  { "-E",             eslARG_REAL,       "10.0", NULL, "x>0",       NULL,  NULL,  REPOPTS,           "report sequences <= this E-value threshold in output",         4 },
  { "-T",             eslARG_REAL,        FALSE, NULL,  NULL,       NULL,  NULL,  REPOPTS,           "report sequences >= this score threshold in output",           4 },
  { "--domE",         eslARG_REAL,       "10.0", NULL, "x>0",       NULL,  NULL,  DOMREPOPTS,        "report domains <= this E-value threshold in output",           4 },
  { "--domT",         eslARG_REAL,        FALSE, NULL,  NULL,       NULL,  NULL,  DOMREPOPTS,        "report domains >= this score cutoff in output",                4 },
/* Control of inclusion thresholds */
  { "--incE",         eslARG_REAL,       "0.01", NULL, "x>0",       NULL,  NULL,  INCOPTS,           "consider sequences <= this E-value threshold as significant",  5 },
  { "--incT",         eslARG_REAL,        FALSE, NULL,  NULL,       NULL,  NULL,  INCOPTS,           "consider sequences >= this score threshold as significant",    5 },
  { "--incdomE",      eslARG_REAL,       "0.01", NULL, "x>0",       NULL,  NULL,  INCDOMOPTS,        "consider domains <= this E-value threshold as significant",    5 },
  { "--incdomT",      eslARG_REAL,        FALSE, NULL,  NULL,       NULL,  NULL,  INCDOMOPTS,        "consider domains >= this score threshold as significant",      5 },
/* Model-specific thresholding for both reporting and inclusion (unused in phmmer)*/
  { "--cut_ga",       eslARG_NONE,        FALSE, NULL, NULL,        NULL,  NULL,  THRESHOPTS,        "use profile's GA gathering cutoffs to set all thresholding",  99 },
  { "--cut_nc",       eslARG_NONE,        FALSE, NULL, NULL,        NULL,  NULL,  THRESHOPTS,        "use profile's NC noise cutoffs to set all thresholding",      99 },
  { "--cut_tc",       eslARG_NONE,        FALSE, NULL, NULL,        NULL,  NULL,  THRESHOPTS,        "use profile's TC trusted cutoffs to set all thresholding",    99 },
/* Control of filter pipeline */
  { "--max",          eslARG_NONE,        FALSE, NULL, NULL,        NULL,  NULL, "--F1,--F2,--F3",   "Turn all heuristic filters off (less speed, more power)",      7 },
  { "--F1",           eslARG_REAL,       "0.02", NULL, NULL,        NULL,  NULL, "--max",            "Stage 1 (MSV) threshold: promote hits w/ P <= F1",             7 },
  { "--F2",           eslARG_REAL,       "1e-3", NULL, NULL,        NULL,  NULL, "--max",            "Stage 2 (Vit) threshold: promote hits w/ P <= F2",             7 },
  { "--F3",           eslARG_REAL,       "1e-5", NULL, NULL,        NULL,  NULL, "--max",            "Stage 3 (Fwd) threshold: promote hits w/ P <= F3",             7 },
  { "--nobias",       eslARG_NONE,        NULL,  NULL, NULL,        NULL,  NULL, "--max",            "turn off composition bias filter",                             7 },
/* Control of E-value calibration */
  { "--EmL",          eslARG_INT,         "200", NULL,"n>0",        NULL,  NULL,  NULL,              "length of sequences for MSV Gumbel mu fit",                   11 },   
  { "--EmN",          eslARG_INT,         "200", NULL,"n>0",        NULL,  NULL,  NULL,              "number of sequences for MSV Gumbel mu fit",                   11 },   
  { "--EvL",          eslARG_INT,         "200", NULL,"n>0",        NULL,  NULL,  NULL,              "length of sequences for Viterbi Gumbel mu fit",               11 },   
  { "--EvN",          eslARG_INT,         "200", NULL,"n>0",        NULL,  NULL,  NULL,              "number of sequences for Viterbi Gumbel mu fit",               11 },   
  { "--EfL",          eslARG_INT,         "100", NULL,"n>0",        NULL,  NULL,  NULL,              "length of sequences for Forward exp tail tau fit",            11 },   
  { "--EfN",          eslARG_INT,         "200", NULL,"n>0",        NULL,  NULL,  NULL,              "number of sequences for Forward exp tail tau fit",            11 },   
  { "--Eft",          eslARG_REAL,       "0.04", NULL,"0<x<1",      NULL,  NULL,  NULL,              "tail mass for Forward exponential tail tau fit",              11 },   
/* other options */
  { "--nonull2",      eslARG_NONE,        NULL,  NULL, NULL,        NULL,  NULL,  NULL,              "turn off biased composition score corrections",               12 },
  { "-Z",             eslARG_REAL,       FALSE,  NULL, "x>0",       NULL,  NULL,  NULL,              "set # of comparisons done, for E-value calculation",          12 },
  { "--domZ",         eslARG_REAL,       FALSE,  NULL, "x>0",       NULL,  NULL,  NULL,              "set # of significant seqs, for domain E-value calculation",   12 },
  { "--seed",         eslARG_INT,         "42",  NULL, "n>=0",      NULL,  NULL,  NULL,              "set RNG seed to <n> (if 0: one-time arbitrary seed)",         12 },
  { "--qformat",      eslARG_STRING,      NULL,  NULL, NULL,        NULL,  NULL,  NULL,              "assert query <seqfile> is in format <s>: no autodetection",   12 },
  { "--tformat",      eslARG_STRING,      NULL,  NULL, NULL,        NULL,  NULL,  NULL,              "assert target <seqdb> is in format <s>>: no autodetection",   12 },
/* Sledgehmmer specific options */  
  { "--qblock",       eslARG_INT,        "10",    NULL, NULL,       NULL,  NULL,  NULL,              "max queries handled at a time",                               12 },
  { "--dbblock",      eslARG_INT,        "1000",  NULL, NULL,       NULL,  NULL,  NULL,              "max targets handled at a time",                               12 },
  { "--filter",       eslARG_NONE,       FALSE,   NULL, NULL,       NULL,  NULL,  NULL,              "filter mode toggle",                                          12 },
  { "--halt",         eslARG_INT,        "-1",    NULL, NULL,       NULL,  NULL,  NULL,              "stop after n sequences - for debugging",                      12 },
  { "--train_frac",   eslARG_REAL,        "0.9",  NULL, NULL,       NULL,  NULL,  NULL,              "train proportion after test limit reached",                   12 },
  { "--test_limit",   eslARG_INT,        "75000", NULL, NULL,       NULL,  NULL,  NULL,              "required minimum number of test sequences",                   12 },
  { "--init_chunk",   eslARG_INT,        "20",    NULL, NULL,       NULL,  NULL,  NULL,              "query chunk size till test limit",                            12 },
  { "--train_only",   eslARG_NONE,        FALSE,  NULL, NULL,       NULL,  NULL,  NULL,              "growing only train set",                                      12 },
  { "--test_only",    eslARG_NONE,        FALSE,  NULL, NULL,       NULL,  NULL,  NULL,              "growing only test set",                                       12 },
  { "--suppress",     eslARG_NONE,        FALSE,   NULL, NULL,      NULL,  NULL,  NULL,              "turn off progress bar",                                       12 },
  { "--task_id",      eslARG_INT,        "0",     NULL, NULL,       NULL,  NULL,  NULL,              "slurm array task ID or shard number",                         12 },
  { "--shard_id",     eslARG_STRING,     "",      NULL, NULL,       NULL,  NULL,  NULL,              "second level shard number as a string",                       12 },
  { "--plow",         eslARG_REAL,       "0.00",  NULL, NULL,       NULL,  NULL,  NULL,              "pid lower limit for accepting a sequence",                    12 },
  { "--phigh",        eslARG_REAL,       "1.00",  NULL, NULL,       NULL,  NULL,  NULL,              "pid upper limit for accepting a sequence",                    12 },
  { "--qsize",        eslARG_INT,        "1",     NULL, NULL,       NULL,  NULL,  NULL,              "num queries per thread",                                      12 },
  { "--train_path",   eslARG_STRING,     "../results/train_seqs",   NULL,  NULL,  NULL, NULL, NULL,  "filename for train seqs with path",                           12 },
  { "--test_path",    eslARG_STRING,     "../results/test_seqs",    NULL,  NULL,  NULL, NULL, NULL,  "filename for test seqs with path",                            12 },
  { "--discard_path", eslARG_STRING,     "../results/discard_seqs", NULL,  NULL,  NULL, NULL, NULL,  "filename for discard seqs with path",                         12 },
  { "-o",             eslARG_STRING,     "-",     NULL, NULL,       NULL,  NULL,  NULL,              "output file name",                                             2 },
#ifdef HMMER_THREADS
  { "--cpu",          eslARG_INT,        "1","HMMER_NCPU", "n>=0",  NULL,  NULL,  CPUOPTS,           "number of CPU cores to use",                                  12 },
#endif
#ifdef HMMER_MPI
  { "--stall",        eslARG_NONE,       FALSE, NULL, NULL,         NULL,  "--mpi", NULL,            "arrest after start: for debugging MPI under gdb",             12 },  
  { "--mpi",          eslARG_NONE,       FALSE, NULL, NULL,         NULL,  NULL,  MPIOPTS,           "run as an MPI parallel program",                              12 },
#endif

  /* Restrict search to subset of database - hidden because these flags are
   *   (a) currently for internal use
   *   (b) probably going to change
   * Doesn't work with MPI
   */
  { "--restrictdb_stkey", eslARG_STRING, "0",   NULL,  NULL,        NULL,  NULL,  NULL,            "Search starts at the sequence with name <s> (not with MPI)",    99 },
  { "--restrictdb_n",eslARG_INT,         "-1",  NULL,  NULL,        NULL,  NULL,  NULL,            "Search <j> target sequences (starting at --restrictdb_stkey)",  99 },
  { "--ssifile",    eslARG_STRING,       NULL,  NULL,  NULL,        NULL,  NULL,  NULL,            "restrictdb_x values require ssi file. Override default to <s>", 99 },

 {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <qdb> <seqdb>";
static char banner[] = "filter hits from qdb to seqdb";

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

  if (esl_opt_ArgNumber(go)                 != 2)    { if (puts("Incorrect number of command line arguments.")      < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((*ret_qfile  = esl_opt_GetArg(go, 1)) == NULL) { if (puts("Failed to get <qdb> argument on command line")     < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((*ret_dbfile = esl_opt_GetArg(go, 2)) == NULL) { if (puts("Failed to get <seqdb> argument on command line")   < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  /* Validate any attempted use of stdin streams */
  if (strcmp(*ret_qfile, "-") == 0 && strcmp(*ret_dbfile, "-") == 0) 
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
  ESL_SQFILE      *qfp              = NULL;               /* open qfile                 */
  ESL_SQ_BLOCK    *qdb              = NULL;               /* query database             */
  ESL_SQ_BLOCK    *sqdb             = NULL;               /* target database            */
  ESL_SQ_BLOCK    *cdb              = NULL;               /* chunk database             */
  int              qformat          = eslSQFILE_UNKNOWN;  /* format of qfile            */
  int              dbformat         = eslSQFILE_UNKNOWN;  /* format of dbfile           */
  ESL_ALPHABET    *abc              = NULL;               /* sequence alphabet          */
  char            *results          = NULL;               /* results of seq inclusion   */
  struct cfg_s     cfg;                                   /* configuration data         */
  SLEDGE_INFO      si;                                    /* Sledgehmmer info           */
  int              num_records;                           /* Total sequences in qdb     */
  int              seq_ctr;                               /* Count sequences processed  */
  char            *out_file;                              /* Discard sequence file name */
  int              len;                                   /* Calculate file name length */
  clock_t          start_time;                            /* Keep track of start time   */
  clock_t          elapsed_time;                          /* Elapsed time               */
  bool             suppress;                              /* Suppress progress update   */

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
  si.qsize        = esl_opt_GetInteger (go, "--qsize");
  si.cores        = esl_opt_GetInteger (go, "--cpu");
  si.incE         = esl_opt_GetReal    (go, "--incE");
  si.db_size      = esl_opt_GetReal    (go, "-Z");
  si.out_path     = esl_opt_GetString  (go, "-o");
  si.task_id      = esl_opt_GetInteger (go, "--task_id");
  suppress        = esl_opt_GetBoolean (go, "--suppress");
  si.filter       = TRUE;
  si.out_fp       = NULL;
  si.offset       = 0;

  /* Set unused Sledgehmmer parameters to NULL */
  bool   train_only   = FALSE;
  bool   test_only    = FALSE;
  float  train_frac   = 0.0;
  int    init_chunk   = 0;
  int    test_limit   = 0;
  char  *train_path   = NULL;
  char  *test_path    = NULL;
  char  *discard_path = NULL;

  /* If caller declared input formats, decode them */
  if (esl_opt_IsOn(go, "--qformat"))
  {  
    qformat = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--qformat"));
    if (qformat == eslSQFILE_UNKNOWN) p7_Fail("%s is not a recognized input sequence file format\n", esl_opt_GetString(go, "--qformat"));
  }
  if (esl_opt_IsOn(go, "--tformat"))
  {
    dbformat = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--tformat"));
    if (dbformat == eslSQFILE_UNKNOWN) p7_Fail("%s is not a recognized sequence database file format\n", esl_opt_GetString(go, "--tformat"));
  }
  
  /* Initialize alphabet */
  abc = esl_alphabet_Create(eslAMINO);

  /* Open and read the target database file */
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

  /* Open and read the query database file */
  status =  esl_sqfile_OpenDigital(abc, cfg.qfile, qformat, p7_SEQDBENV, &qfp);
  if      (status == eslENOTFOUND) p7_Fail("Failed to open query sequence database %s for reading\n",       cfg.qfile);
  else if (status == eslEFORMAT)   p7_Fail("Query sequence database file %s is empty or misformatted\n",    cfg.qfile);
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        p7_Fail("Unexpected error %d opening target sequence database file %s\n", status, cfg.qfile);
  
  qdb = esl_sq_CreateDigitalBlock(esl_opt_GetInteger(go, "--qblock"), abc);
  if (qdb == NULL) p7_Fail("Failed to allocate sequence database block");

  status = read_cust(qfp, qdb);
  if (status != eslOK) p7_Fail("Failed to read target sequences from database");
  esl_sqfile_Close(qfp);

  /* Allocate a database to process each chunk of queries */
  cdb = esl_sq_CreateDigitalBlock(si.cores * si.qsize, abc); /* Allocate for max possible */
  if (cdb == NULL) p7_Fail("Failed to allocate chunk database block");

  /* Calculate total number of queries to process */
  num_records = qdb->count;

  /* Allocate large buffer for filter results */
  ESL_ALLOC(results, BUFFER_MAX * sizeof(char));
  
  /* Open the output file (or stdout) for permanent writing */

  if (strcmp(si.out_path, "-") == 0)
    si.out_fp = stdout;
  else
  {
    len = snprintf(NULL, 0, "%s_%d_%s.txt", si.out_path, si.task_id, esl_opt_GetString(go, "--shard_id")) + 1;
    ESL_ALLOC(out_file, len*sizeof(char));
    snprintf(out_file, len, "%s_%d_%s.txt", si.out_path, si.task_id, esl_opt_GetString(go, "--shard_id"));
    if ((si.out_fp = fopen(out_file, "w")) == NULL) p7_Fail("Failed to open output file %s for writing\n", out_file);
  }

  /* Header in output file/stream */
  fprintf(si.out_fp, "%-15s %-15s %6s %7s %12s %8s\n",
    "Query", "Target", "PID", "Length", "E-value", "Decision");

  fprintf(si.out_fp, "--------------- --------------- ------ ------- ------------ --------\n"); /* Just pretty printing */

  /* Loop through all sequences in database */
  start_time = get_coarse_time();

  for (seq_ctr = 0; seq_ctr < num_records; seq_ctr++)
  {
    add_seq(cdb, qdb->list + seq_ctr, abc);

    if (((cdb->count == si.cores * si.qsize) || (seq_ctr == num_records-1)))
    {
      if (seq_ctr >= num_records)
        si.cores = ceil(cdb->count / si.qsize); /* Not optimized, but this is called only for last chunk */

      status = assign_master(go, cdb, sqdb, &si, results);

      cdb->count = 0;
    }

    if ((!suppress) && (seq_ctr % 1000 == 0)) /* How often do we want to update the progress bar? It has to be db dependent*/
    {
      elapsed_time = get_coarse_time() - start_time;
      print_progress(seq_ctr+1, num_records, elapsed_time);

      fflush(stderr);
    }
  }

  /* Flush the buffer one last time if buffer isn't empty */
  if (si.offset > 0)
    fwrite(results, 1, si.offset, si.out_fp);

  /* Print progress for one final time */
  elapsed_time = get_coarse_time() - start_time;;
  if (!suppress) print_progress(seq_ctr, num_records, elapsed_time);

  /* Close all files if used */

  if (si.out_fp != stdout)
  {
    fclose(si.out_fp);
    free(out_file);
  }

  /* Mop the floor, empty the trash, free the souls and leave */
  
  free(results);
  
  esl_sq_DestroyBlock(qdb);
  esl_sq_DestroyBlock(sqdb);
  esl_sq_DestroyBlock(cdb);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);

  ERROR:
    exit(status);

  exit(status);
}
