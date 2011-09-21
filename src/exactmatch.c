#include "hmmer.h"
#include <sys/times.h>

#include "easel.h"
#include <string.h>

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles   reqs   incomp              help                                                      docgroup*/
  { "-h",           eslARG_NONE,      FALSE, NULL, NULL,    NULL,  NULL,  NULL,    "show brief help on version and usage",                      1 },

  { "--out",      eslARG_STRING,     "none", NULL, NULL,    NULL,  NULL,  NULL,    "save list of hits to file <s>  ('-' writes to stdout)",     2 },
  { "--count_only", eslARG_NONE,      FALSE, NULL, NULL,    NULL,  NULL,  NULL,    "compute just counts, not locations",                        2 },

  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};


static char usage[]  = "[options] <fmfile> <qfile>";
static char banner[] = "Find all instances of each <qfile> sequence in the database represented by <fmfile>";


static void
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_fmfile, char **ret_qfile)
{
  ESL_GETOPTS *go = NULL;

  if ((go = esl_getopts_Create(options))     == NULL)     p7_Die("Internal failure creating options object");
  if (esl_opt_ProcessEnvironment(go)         != eslOK)  { printf("Failed to process environment: %s\n", go->errbuf); goto ERROR; }
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { printf("Failed to parse command line: %s\n", go->errbuf); goto ERROR; }
  if (esl_opt_VerifyConfig(go)               != eslOK)  { printf("Failed to parse command line: %s\n", go->errbuf); goto ERROR; }

  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) {
      p7_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      puts("\nBasic options:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 120=textwidth*/

      puts("\nSpecial options:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); /* 2= group; 2 = indentation; 120=textwidth*/

      exit(0);
  }

  if (esl_opt_ArgNumber(go)                  != 2)     { puts("Incorrect number of command line arguments.");      goto ERROR; }
  if ((*ret_fmfile = esl_opt_GetArg(go, 1)) == NULL)  { puts("Failed to get <fmfile> argument on command line"); goto ERROR; }
  if ((*ret_qfile  = esl_opt_GetArg(go, 2)) == NULL)  { puts("Failed to get <qfile> argument on command line");   goto ERROR; }

  /* Validate any attempted use of stdin streams */
  if (esl_strcmp(*ret_fmfile, "-") == 0 && esl_strcmp(*ret_qfile, "-") == 0) {
    puts("Either <fmfile> or <qfile> may be '-' (to read from stdin), but not both.");
    goto ERROR;
  }

  *ret_go = go;
  return;

ERROR:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, argv[0], usage);
  puts("\nwhere basic options are:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/
  printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
  exit(1);
}

//see: http://c-faq.com/stdio/commaprint.html
char *
commaprint(unsigned long n)
{
    static int comma = ',';
    static char retbuf[30];
    char *p = &retbuf[sizeof(retbuf)-1];
    int i = 0;
    *p = '\0';

    do {
      if(i%3 == 0 && i != 0)
              *--p = comma;
      *--p = '0' + n % 10;
      n /= 10;
      i++;
    } while(n != 0);

    return p;
}

static int
output_header(FM_METADATA *meta, FILE *ofp, const ESL_GETOPTS *go, char *fmfile, char *qfile)
{
  p7_banner(ofp, go->argv[0], banner);

  fprintf(ofp, "# input binary-formatted HMMER database:   %s\n", fmfile);
  fprintf(ofp, "# input file of query sequences:           %s\n", qfile);

  if (esl_opt_IsUsed(go, "--out")) {
    fprintf(ofp, "# output file containing list of hits:     ");
    char *outfile = esl_opt_GetString(go, "--out");
    if (esl_strcmp(outfile, "-"))
      fprintf(ofp, "stdout\n");
    else
      fprintf(ofp, "%s\n", outfile);
  }

  if (esl_opt_IsUsed(go, "--count_only"))
    fprintf(ofp, "# output only counts, not hit locations\n");

  char *alph;
  if (meta->alph_type == fm_DNA)
    alph = "dna";
  else if (meta->alph_type == fm_DNA_full)
    alph = "dna_full";
  else if (meta->alph_type == fm_AMINO)
      alph = "amino";
  fprintf(ofp, "# alphabet     :                           %s\n", alph);

  fprintf(ofp, "# bin_length   :                           %d\n", meta->freq_cnt_b);
  fprintf(ofp, "# suffix array sample rate:                %d\n", meta->freq_SA);
  fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
  return eslOK;
}



/* Function:  getFMHits()
 * Synopsis:  For a given interval, identify the position in original text for each element
 *            of interval
 * Purpose:   Implement Algorithm 3.7 (p17) of Firth paper (A Comparison of BWT Approaches
 *            to String Pattern Matching). Most of the meat is in the method of counting
 *            characters - bwt_getOccCount, which depends on compilation choices.
 */
//#ifndef FMDEBUG
//inline
//#endif
int
getFMHits( FM_DATA *fm, FM_CFG *cfg, FM_INTERVAL *interval, int block_id, int hit_offset, int hit_length, FM_HIT *hits_ptr, int fm_direction) {

  int i, j, len = 0;

  for (i = interval->lower;  i<= interval->upper; i++) {
    j = i;
    len = 0;

    while ( j != fm->term_loc && (j & cfg->maskSA)) { //go until we hit a position in the full SA that was sampled during FM index construction
      uint8_t c = getChar( cfg->meta->alph_type, j, fm->BWT);
      j = fm_getOccCount (fm, cfg, j-1, c);
      j += abs(fm->C[c]);
      len++;
    }


    hits_ptr[hit_offset + i - interval->lower].block     = block_id;
    hits_ptr[hit_offset + i - interval->lower].direction = fm_direction;
    hits_ptr[hit_offset + i - interval->lower].length    = hit_length;

    hits_ptr[hit_offset + i - interval->lower].start     = len + (j==fm->term_loc ? 0 : fm->SA[ j >> cfg->shiftSA ]) ; // len is how many backward steps we had to take to find a sampled SA position
    if (fm_direction == fm_backward)
      hits_ptr[hit_offset + i - interval->lower].start  +=  hit_length - 1 ;

  }

  return eslOK;

}

/* hit_sorter(): qsort's pawn, below */
static int
hit_sorter(const void *a, const void *b)
{
  FM_HIT *h1 = (FM_HIT*)a;
  FM_HIT *h2 = (FM_HIT*)b;

  if      (h1->sortkey > h2->sortkey) return  1;
  else if (h1->sortkey < h2->sortkey) return -1;
  else {
    if      (h1->direction > h2->direction) return 1;
    else if (h1->direction < h2->direction) return -1;
    else {
      if  (h1->start > h2->start) return  1;
      else                        return -1;
    }
  }
}

/* Function:  main()
 * Synopsis:  Run set of queries against an FM
 * Incept:    TJW, Fri Dec 24 21:30:51 MST 2010 [Tucson]
 * Purpose:   Read in a FM and a file of query sequences.
 *            For each query, find matching FM interval, then collect positions in
 *            the original text T for the corresponding occurences. These positions
 *            are 0-based (so first character is position 0).
 */
int
main(int argc,  char *argv[]) {

  void* tmp; // used for RALLOC calls
  clock_t t1, t2;
  struct tms ts1, ts2;
  char *fname_fm      = NULL;
  char *fname_queries = NULL;
  char *inv_alph      = NULL;
  char *alph          = NULL;
  FM_HIT *hits        = NULL;
  char *line          = NULL;
  int status        = eslOK;
  int hit_cnt       = 0;
  int hit_indiv_cnt = 0;
  int miss_cnt      = 0;
  int hit_num       = 0;
  int hit_num2       = 0;
  int hits_size     = 0;
  int i;
  int count_only    = 0;

  FM_DATA *fmsf;
  FM_DATA *fmsb;
  FM_INTERVAL interval;
  FILE* fp_fm = NULL;
  FILE* fp = NULL;
  FILE* out = NULL;
  char *outname = NULL;

  ESL_GETOPTS     *go  = NULL;    /* command line processing                 */
  FM_CFG *cfg;
  FM_METADATA *meta;

  //start timer
  t1 = times(&ts1);

  process_commandline(argc, argv, &go, &fname_fm, &fname_queries);


  if (esl_opt_IsOn(go, "--out")) {
    outname = esl_opt_GetString(go, "--out");
    if ( esl_strcmp ("-", outname) == 0 ) {
      out = stdout;
      outname = "stdout";
    } else {
      out = fopen(optarg,"w");
    }
  }

  if (esl_opt_IsOn(go, "--count_only"))
    count_only = 1;




  if((fp_fm = fopen(fname_fm, "rb")) == NULL)
    esl_fatal("Cannot open file `%s': ", fname_fm);


  ESL_ALLOC(cfg, sizeof(FM_CFG));
  ESL_ALLOC(cfg->meta, sizeof(FM_METADATA));
  meta = cfg->meta;
  cfg->occCallCnt = 0;

  readFMmeta( fp_fm, cfg->meta);

  //read in FM-index blocks
  ESL_ALLOC(fmsf, cfg->meta->block_count * sizeof(FM_DATA) );
  for (i=0; i<meta->block_count; i++) {
    readFM( fp_fm, fmsf+i, cfg->meta, 1 );

    if (!meta->fwd_only) {
      ESL_ALLOC(fmsb, meta->block_count * sizeof(FM_DATA) );
      readFM( fp_fm, fmsb+i, cfg->meta, 0 );
      fmsb[i].SA = fmsf[i].SA;
      fmsb[i].T = fmsf[i].T;
    }
  }
  fclose(fp_fm);

  output_header(cfg->meta, stdout, go, fname_fm, fname_queries);


  /* initialize a few global variables, then call initGlobals
   * to do architecture-specific initialization
   */
  cfg->maskSA       =  meta->freq_SA - 1;
  cfg->shiftSA      =  meta->SA_shift;
  fm_initConfig(cfg);


  fm_createAlphabet(meta->alph_type, &alph, &inv_alph, &(meta->alph_size), NULL); // don't override charBits


  fp = fopen(fname_queries,"r");
  if (fp == NULL)
    esl_fatal("Unable to open file %s\n", fname_queries);

  ESL_ALLOC(line, FM_MAX_LINE * sizeof(char));

  hits_size = 200;
  ESL_ALLOC(hits, hits_size * sizeof(FM_HIT));


  while(fgets(line, FM_MAX_LINE, fp) ) {
    int qlen=0;
    while (line[qlen] != '\0' && line[qlen] != '\n')  qlen++;
    if (line[qlen] == '\n')  line[qlen] = '\0';

    hit_num = 0;


    for (i=0; i<meta->block_count; i++) {

      getSARangeForward(fmsb+i, cfg, line, inv_alph, &interval);// yes, use the backward fm to produce a forward search on the forward fm
      if (interval.lower>0 && interval.lower <= interval.upper) {
        int new_hit_num =  interval.upper - interval.lower + 1;
        hit_num += new_hit_num;
        if (!count_only) {
          if (hit_num > hits_size) {
            hits_size = 2*hit_num;
            ESL_RALLOC(hits, tmp, hits_size * sizeof(FM_HIT));
          }
          //even though I used fmsb above, use fmsf here, since we'll now do a backward trace
          //in the FM-index to find the next sampled SA position
          getFMHits(fmsf+i, cfg, &interval, i, hit_num-new_hit_num, qlen, hits, fm_forward);
        }
      }


      /* find reverse hits, using backward search on the forward FM*/
      if (!meta->fwd_only) {
        getSARangeReverse(fmsf+i, cfg, line, inv_alph, &interval);
        if (interval.lower>0 && interval.lower <= interval.upper) {
          int new_hit_num =  interval.upper - interval.lower + 1;
          hit_num += new_hit_num;
          if (!count_only) {
            if (hit_num > hits_size) {
              hits_size = 2*hit_num;
              ESL_RALLOC(hits, tmp, hits_size * sizeof(FM_HIT));
            }
            getFMHits(fmsf+i, cfg, &interval, i, hit_num-new_hit_num, qlen, hits, fm_backward);
          }
        }
      }
    }



    if (hit_num > 0) {
      if (count_only) {
        hit_cnt++;
        hit_indiv_cnt += hit_num;
      } else {
        hit_num2 = 0;

        //for each hit, identify the sequence id and position within that sequence
        for (i = 0; i< hit_num; i++) {
          int block = hits[i].block;
          int seq_offset = computeSequenceOffset( fmsf, meta, block, hits[i].start);
          int pos =  ( hits[i].start - meta->seq_data[ seq_offset ].offset) + meta->seq_data[ seq_offset ].start - 1;

          //verify that the hit doesn't extend beyond the bounds of the target sequence
          if (hits[i].direction == fm_forward) {
            if (pos + hits[i].length > meta->seq_data[ seq_offset ].length ) {
              hits[i].block  = hits[i].sortkey  = hits[i].start  = -1;  // goes into the next sequence, so it should be ignored
              continue;
            }
          } else { //backward
            if (pos - hits[i].length + 1 < 0 ) {
              hits[i].block  = hits[i].sortkey  = hits[i].start  = -1; // goes into the previous sequence, so it should be ignored
              continue;
            }
          }
          hit_num2++; // legitimate hit

          //reuse hit variables.  Now "block" has the index into the matching sequence (in meta), and "start" has the pos within that sequence
          hits[i].block   = seq_offset;
          hits[i].start   = pos;
          hits[i].sortkey = meta->seq_data[ seq_offset ].id;
        }
        if (hit_num2 > 0)
          hit_cnt++;


        //now sort according the the sequence_id corresponding to that seq_offset
        qsort(hits, hit_num, sizeof(FM_HIT), hit_sorter);

        //skim past the skipped entries
        i = 0;
        while ( i < hit_num ) {
          if (hits[i].block != -1 )
            break;  //
          i++;
        }

        if (i < hit_num) {
          if (out != NULL) {
            fprintf (out, "%s\n",line);
            //fprintf (out, "\t%10s (%8d %s)\n",meta->seq_data[ hits[i].block ].name, hits[i].start, (hits[i].direction==fm_forward?"+":"-"));
            fprintf (out, "    %8d %s %10s\n", hits[i].start, (hits[i].direction==fm_forward?"1":"0"), meta->seq_data[ hits[i].block ].name);
          }
          hit_indiv_cnt++;
          i++; // skip the first one, since I'll be comparing each to the previous
          for (  ; i< hit_num; i++) {
            if ( //meta->seq_data[ hits[i].block ].id != meta->seq_data[ hits[i-1].block ].id ||
                 hits[i].sortkey   != hits[i-1].sortkey ||  //sortkey is seq_data[].id
                 hits[i].direction != hits[i-1].direction ||
                 hits[i].start     != hits[i-1].start )
            {
              if (out != NULL)
                //fprintf (out, "\t%10s (%8d %s)\n",meta->seq_data[ hits[i].block ].name, hits[i].start, (hits[i].direction==fm_forward?"+":"-"));
                fprintf (out, "    %8d %s %10s\n", hits[i].start, (hits[i].direction==fm_forward?"1":"0"), meta->seq_data[ hits[i].block ].name);
              hit_indiv_cnt++;
            }
          }
          if (out != NULL)
            fprintf (out, "\n");
        }
      }
    } else {
      miss_cnt++;
    }


  }

  for (i=0; i<meta->block_count; i++) {
    freeFM( fmsb+i, 1 );
    if (!meta->fwd_only)
      freeFM( fmsf+i, 0 );
  }

  for (i=0; i<meta->seq_count; i++)
    free (meta->seq_data[i].name);

  free (meta->seq_data);
  free (meta);
  free (hits);
  free (line);

  fclose(fp);
  fm_destroyConfig(cfg);
  free(cfg);

  // compute and print the elapsed time in millisec
  t2 = times(&ts2);
  {
    double clk_ticks = sysconf(_SC_CLK_TCK);
    double elapsedTime = (t2-t1)/clk_ticks;
    double throughput = cfg->occCallCnt/elapsedTime;

    fprintf (stderr, "hit: %-10d  (%d)\n", hit_cnt, hit_indiv_cnt);
    fprintf (stderr, "miss:%-10d\n", miss_cnt);
    fprintf (stderr, "run time:  %.2f seconds\n", elapsedTime);
    fprintf (stderr, "occ calls: %12s\n", commaprint(cfg->occCallCnt));
    fprintf (stderr, "occ/sec:   %12s\n", commaprint(throughput));
  }

  exit(eslOK);


ERROR:
  printf ("failure allocating memory for hits\n");
  exit(status);


}





/*****************************************************************
 * @LICENSE@
 *****************************************************************/
