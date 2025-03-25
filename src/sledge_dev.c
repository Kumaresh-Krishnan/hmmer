/* phmmer: search a protein sequence against a protein database
 */
#include <p7_config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_exponential.h"
#include "esl_getopts.h"
#include "esl_gumbel.h"
#include "esl_vectorops.h"
#include "esl_msa.h"
#include "esl_msafile.h"
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

int
read_cust(ESL_SQFILE *sqfp, ESL_SQ_BLOCK *sqBlock)
{
  int     status = eslOK;

  sqBlock->count = 0;

  while((status = esl_sqio_Read(sqfp, sqBlock->list + sqBlock->count)) == eslOK)
    ++sqBlock->count;
  
  /* EOF will be returned only in the case were no sequences were read */
  if (status == eslEOF && sqBlock->count > 0) status = eslOK;
  
  sqBlock->complete = TRUE;

  return status;
}

int
append_target(RESULT_INFO *result, const char *target_id, float pid, float eval, int length)
{    
  /* Reallocate memory and double the capacity if current max is reached */
  if (result->n_targets >= result->capacity)
  {  
    int new_capacity = result->capacity * 2;
    char **new_targets = realloc(result->targets, new_capacity * sizeof(char*));
    float *new_pids    = realloc(result->pids, new_capacity * sizeof(float));
    float *new_evals   = realloc(result->evals, new_capacity * sizeof(float));
    int *new_lengths   = realloc(result->lengths, new_capacity * sizeof(int));
    
    if(new_targets == NULL || new_pids == NULL || new_evals == NULL) {
      return eslENORESULT;
    }
    
    /* Place old values in the newly allocated object */
    result->targets  = new_targets;
    result->pids     = new_pids;
    result->evals    = new_evals;
    result->capacity = new_capacity;
    result->lengths  = new_lengths;
  }
  
  /* Append new data */
  result->targets[result->n_targets] = strdup(target_id);
  
  if(result->targets[result->n_targets] == NULL) {
      return eslENORESULT;
  }
  
  result->pids[result->n_targets]  = pid;
  result->evals[result->n_targets] = eval;
  result->lengths[result->n_targets] = length;
  result->n_targets++;
  result->decision = REJECT; /* Only storing the rejected sequences */

  return eslOK;
}

int
destroy_info(WORKER_INFO *info, int num_cores)
{
  /* Free all the objects created inside info across all cores and queries */
  for (int i = 0; i < num_cores; ++i)
  { 
    for (int j=0; j<info[i].num_queries; j++)
    {
      for (int k=0; k<info[i].result[j].n_targets; k++)
        free((info[i].result[j]).targets[k]);

      free(info[i].result[j].targets);
      free(info[i].result[j].pids);
      free(info[i].result[j].evals);
      free(info[i].result[j].lengths);
      free(info[i].result[j].source);
    }

    free(info[i].result);   /* Free the result array */
  }

  free(info);

  return eslOK;
}

int
add_seq(ESL_SQ_BLOCK *db, ESL_SQ *seq, const ESL_ALPHABET *abc)
{

  if (db->count+1 > db->listSize)
    esl_sq_BlockGrowTo(db, db->listSize*2, 1, abc);

  esl_sq_Copy(seq, db->list + db->count);
  db->count++;

  return eslOK;
}

void
add_to_buffer(char *buffer, char *line, size_t line_length, SLEDGE_INFO *si)
{  
  if (si->offset + line_length >= BUFFER_MAX)
  {
    fwrite(buffer, 1, si->offset, si->out_fp);
    si->offset = 0;
  }

  memcpy(buffer + si->offset, line, line_length);
  si->offset += line_length;
}

int
store_results(ESL_THREADS *threadObj, SLEDGE_INFO *si, char *results)
{
  WORKER_INFO *info       = NULL;
  int          result_ctr = 0;
  char         line[128];
  
  for(int i=0; i<si->cores; i++)
  {
    info = (WORKER_INFO*) (threadObj->data[i]);
    
    for(int j=0; j<info->num_queries; j++)
    {
      if (si->filter) /* Want information about the failed hits */
      {
        for(int hit=0; hit<(info->result[j]).n_targets; hit++)
        {
            snprintf(line, sizeof(line), "%-15s %-15s %6.2f %7d %12.2E %8c\n",
              (info->result[j]).source,
              (info->result[j]).targets[hit],
              (info->result[j]).pids[hit],
              (info->result[j]).lengths[hit],
              exp((info->result[j]).evals[hit]),
              (info->result[j]).decision);

              add_to_buffer(results, line, strlen(line), si);
        }    
      }

      else
      {
        results[result_ctr] = (info->result[j]).decision; /* Only store ACCEPT/REJECT based on PID/E-val */
        result_ctr++;
      }
    }
  }

  if (!(si->filter)) /* Terminate ACCEPT/REJECT string */
    results[result_ctr] = '\0';

  return eslOK;
}

double get_coarse_time()
{
  struct timespec ts;
  
  clock_gettime(CLOCK_MONOTONIC_COARSE, &ts);
  
  return ts.tv_sec + (ts.tv_nsec / 1.0e9);
}

void
format_time(double seconds, char *buffer, size_t buffer_size)
{
  int hours = (int)(seconds / 3600);
  int minutes = (int)((seconds - (hours * 3600)) / 60);
  int secs = (int)(seconds - (hours * 3600) - (minutes * 60));
  snprintf(buffer, buffer_size, "%02d:%02d:%02d", hours, minutes, secs);
}

void
print_progress(int current, int total, double elapsed)
{
  int   bar_width = 50;                        /* Progress bar width            */
  float progress  = (float) current / total;   /* Compute the progress          */
  int   pos       = bar_width * progress;      /* Find position                 */
  char  eta_str[10];                           /* Buffer for formatted ETA time */
  char  total_str[10];                         /* Buffer for elapsed time       */
  
  format_time((elapsed / progress - elapsed),
              eta_str,
              sizeof(eta_str));                /* Convert ETA to HH:MM:SS       */
  format_time(elapsed,
              total_str,
              sizeof(total_str));              /* Convert ETA to HH:MM:SS       */

  fprintf(stderr, "\rProgress: [");            /* Move cursor to beginning      */
  for (int i = 0; i < bar_width; i++) {
      if (i < pos)
        fprintf(stderr, "\xe2\x96\x88");       /* Progress block                */
      else fprintf(stderr, " ");
  }

  fprintf(stderr, "] %.2f%% | %.2f it/s | Elapsed: %s | ETA: %s ", (progress * 100), (elapsed == 0) ? 0.0 : (current/elapsed), total_str, eta_str);

  if (current == total)
    fprintf(stderr, "\n");
}

int
assign_master(ESL_GETOPTS *go, ESL_SQ_BLOCK *qdb, ESL_SQ_BLOCK *sqdb, SLEDGE_INFO *si, char *results)
{
  ESL_STOPWATCH   *timer     = NULL;               /* for timing                                       */
  ESL_THREADS     *threadObj = NULL;               /* stores information about all threads             */ 
  WORKER_INFO     *info      = NULL;               /* Info array to hold data for all threads          */
  int              status    = eslOK;              /* for status                                       */
  int              i;                              /* Counter for loops                                */

  /* Initialize thread object and relevant data */

  ESL_ALLOC(info, (ptrdiff_t) sizeof(*info) * si->cores);

  threadObj = esl_threads_Create(&pipeline_thread);
  int q_ctr = 0;
  
  /* Assign data for each thread along with query chunk information */
  for(i=0; i<si->cores; i++)
  {    
    info[i].si    = si;       /* Sledgehmmer info           */
    info[i].sqdb  = sqdb;     /* Target database            */
    info[i].qdb   = qdb;      /* Query database             */
    info[i].go    = go;       /* Commandline options        */
    info[i].start = q_ctr;    /* Query subset start in file */
    
    /* Assign only overflow queries for last thread */
    if (si->qsize < qdb->count - q_ctr)
    {
      q_ctr = q_ctr + si->qsize;
      info[i].num_queries = si->qsize; 
    }
    else info[i].num_queries = qdb->count - q_ctr;
    
    ESL_ALLOC(info[i].result, (info[i].num_queries)*sizeof(RESULT_INFO));
    
    esl_threads_AddThread(threadObj, &info[i]);
  }

  esl_threads_WaitForStart(threadObj);
  esl_threads_WaitForFinish(threadObj);

  status = store_results(threadObj, si, results);
  
  /* Free objects */
  destroy_info(info, si->cores);
  esl_threads_Destroy(threadObj);

  return status;

  ERROR:
    return status;
}

void
pipeline_thread(void *arg)
{  
  int i;
  int q_ctr;
  int status;
  int workeridx;
  int seed;
  ESL_ALPHABET *abc;
  WORKER_INFO *info;
  ESL_THREADS *obj;
  P7_BUILDER *bld = NULL;

  impl_Init();
  abc = esl_alphabet_Create(eslAMINO);
  obj = (ESL_THREADS *) arg;
  esl_threads_Started(obj, &workeridx);

  info = (WORKER_INFO *) esl_threads_GetData(obj, workeridx);
  ESL_GETOPTS *go = info->go;
  P7_BG *bg = p7_bg_Create(abc);

  bld = p7_builder_Create(NULL, abc);
  if ((seed = esl_opt_GetInteger(go, "--seed")) > 0)
    {				
      /* a little wasteful - we're blowing a couple of usec by reinitializing */
      esl_randomness_Init(bld->r, seed);
      bld->do_reseeding = TRUE;
    }
  bld->EmL = esl_opt_GetInteger(go, "--EmL");
  bld->EmN = esl_opt_GetInteger(go, "--EmN");
  bld->EvL = esl_opt_GetInteger(go, "--EvL");
  bld->EvN = esl_opt_GetInteger(go, "--EvN");
  bld->EfL = esl_opt_GetInteger(go, "--EfL");
  bld->EfN = esl_opt_GetInteger(go, "--EfN");
  bld->Eft = esl_opt_GetReal   (go, "--Eft");
    
  status = p7_builder_LoadScoreSystem(bld, esl_opt_GetString(go, "--mx"), esl_opt_GetReal(go, "--popen"), esl_opt_GetReal(go, "--pextend"), bg);
  if (status != eslOK) p7_Fail("Failed to set single query seq score system:\n%s\n", bld->errbuf);

  for(q_ctr=0; q_ctr<info->num_queries; q_ctr++)
  {
    P7_OPROFILE *om = NULL;
    P7_TOPHITS *th = p7_tophits_Create();
    P7_PIPELINE *pli = NULL;

    ESL_SQ *qsq = info->qdb->list + info->start + q_ctr;

    /* Initialize all members of the RESULT_INFO struct */
    (info->result[q_ctr]).n_targets = 0;
    (info->result[q_ctr]).max_pid   = -1.0;
    (info->result[q_ctr]).decision  = ACCEPT;
    (info->result[q_ctr]).capacity  = INITIAL_HITS_CAPACITY;
    (info->result[q_ctr]).source    = strdup(qsq->name);
    (info->result[q_ctr]).targets   = malloc(INITIAL_HITS_CAPACITY * sizeof(char *));
    (info->result[q_ctr]).pids      = malloc(INITIAL_HITS_CAPACITY * sizeof(float));
    (info->result[q_ctr]).evals     = malloc(INITIAL_HITS_CAPACITY * sizeof(float));
    (info->result[q_ctr]).lengths   = malloc(INITIAL_HITS_CAPACITY * sizeof(int));
    
    /* Build the query model */
    p7_SingleBuilder(bld, qsq, bg, NULL, NULL, NULL, &om); /* bypass HMM - only need model */
    pli = p7_pipeline_Create(info->go, om->M, 100, FALSE, p7_SEARCH_SEQS);
    p7_pli_NewModel(pli, om, bg);

    for (i = 0; i < info->sqdb->count; i++)
    {
      ESL_SQ *tsq = info->sqdb->list + i;
      p7_pli_NewSeq(pli, tsq);
      p7_bg_SetLength(bg, tsq->n);
      p7_oprofile_ReconfigLength(om, tsq->n);

      status = p7_Pipeline_cust(pli, om, bg, tsq, qsq, th, &(info->result[q_ctr]), info->si);


      if ((info->result[q_ctr].decision == REJECT) && (!((info->si)->filter))) break;
    }
    
    p7_tophits_Destroy(th);
    p7_pipeline_Destroy(pli);
    p7_oprofile_Destroy(om);    
  }
  
  esl_threads_Finished(obj, workeridx);

  p7_bg_Destroy(bg);
  p7_builder_Destroy(bld);

  return;  
}

/* Original function modified to keep track of decision and qsq instead of ntsq since we want access to query */
int
p7_Pipeline_cust(P7_PIPELINE *pli, P7_OPROFILE *om, P7_BG *bg, const ESL_SQ *sq, const ESL_SQ *qsq, P7_TOPHITS *hitlist, RESULT_INFO *result, SLEDGE_INFO *si)
{
  P7_HIT          *hit     = NULL;        /* ptr to the current hit output data             */
  float            usc, vfsc, fwdsc;      /* filter scores                                  */
  float            filtersc;              /* HMM null filter score                          */
  float            nullsc;                /* null model score                               */
  float            seqbias;               /* squence bias                                   */
  float            seq_score;             /* the corrected per-seq bit score                */
  float            sum_score;             /* the corrected reconstruction score for the seq */
  float            pre_score, pre2_score; /* uncorrected bit scores for seq                 */
  double           P;                     /* P-value of a hit                               */
  double           lnP;                   /* log P-value of a hit                           */
  int              Ld;                    /* # of residues in envelopes                     */
  int              d;
  int              status;

  if (sq->n == 0) return eslOK;    /* silently skip length 0 seqs; they'd cause us all sorts of weird problems */
  if (sq->n > 100000) ESL_EXCEPTION(eslETYPE, "Target sequence length > 100K, over comparison pipeline limit.\n(Did you mean to use nhmmer/nhmmscan?)");

  p7_omx_GrowTo(pli->oxf, om->M, 0, sq->n);    /* expand the one-row omx if needed */

  /* Base null model score (we could calculate this in NewSeq(), for a scan pipeline) */
  p7_bg_NullOne  (bg, sq->dsq, sq->n, &nullsc);

  /* First level filter: the MSV filter, multihit with <om> */
  p7_MSVFilter(sq->dsq, sq->n, om, pli->oxf, &usc);
  seq_score = (usc - nullsc) / eslCONST_LOG2;
  P = esl_gumbel_surv(seq_score,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
  if (P > pli->F1) {
    if (result->n_targets == 0)
      result->decision = ACCEPT;

    return eslOK;
  }
  pli->n_past_msv++;

  /* biased composition HMM filtering */
  if (pli->do_biasfilter)
  {
    p7_bg_FilterScore(bg, sq->dsq, sq->n, &filtersc);
    seq_score = (usc - filtersc) / eslCONST_LOG2;
    P = esl_gumbel_surv(seq_score,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
    if (P > pli->F1) {
      if (result->n_targets == 0)
        result->decision = ACCEPT;

      return eslOK;
    }
  }
  else filtersc = nullsc;
  pli->n_past_bias++;

  /* In scan mode, if it passes the MSV filter, read the rest of the profile */
  if (pli->mode == p7_SCAN_MODELS)
  {
    if (pli->hfp) p7_oprofile_ReadRest(pli->hfp, om);
    p7_oprofile_ReconfigRestLength(om, sq->n);
    if ((status = p7_pli_NewModelThresholds(pli, om)) != eslOK) return status; /* pli->errbuf has err msg set */
  }

  /* Second level filter: ViterbiFilter(), multihit with <om> */
  if (P > pli->F2)
  {
    p7_ViterbiFilter(sq->dsq, sq->n, om, pli->oxf, &vfsc);  
    seq_score = (vfsc-filtersc) / eslCONST_LOG2;
    P  = esl_gumbel_surv(seq_score,  om->evparam[p7_VMU],  om->evparam[p7_VLAMBDA]);
    if (P > pli->F2) {
      if (result->n_targets == 0)
        result->decision = ACCEPT;

      return eslOK;
    }
  }
  pli->n_past_vit++;

  /* Parse it with Forward and obtain its real Forward score. */
  p7_ForwardParser(sq->dsq, sq->n, om, pli->oxf, &fwdsc);
  seq_score = (fwdsc-filtersc) / eslCONST_LOG2;
  P = esl_exp_surv(seq_score,  om->evparam[p7_FTAU],  om->evparam[p7_FLAMBDA]);
  if (P > pli->F3) {
    if (result->n_targets == 0)
      result->decision = ACCEPT;

    return eslOK;
  }
  pli->n_past_fwd++;
  
  /* ok, it's for real. Now a Backwards parser pass, and hand it to domain definition workflow */
  p7_omx_GrowTo(pli->oxb, om->M, 0, sq->n);
  p7_BackwardParser(sq->dsq, sq->n, om, pli->oxf, pli->oxb, NULL);
  status = p7_domaindef_cust(sq, qsq, om, pli->oxf, pli->oxb, pli->fwd, pli->bck,
                              pli->ddef, bg, FALSE, NULL, NULL, NULL,
                              result, &nullsc, si);

  return status;
}

int
p7_domaindef_cust(const ESL_SQ *sq, const ESL_SQ *qsq, P7_OPROFILE *om,
				   P7_OMX *oxf, P7_OMX *oxb, P7_OMX *fwd, P7_OMX *bck, 
				   P7_DOMAINDEF *ddef, P7_BG *bg, int long_target,
				   P7_BG *bg_tmp, float *scores_arr, float *fwd_emissions_arr,
           RESULT_INFO *result, float *nullsc, SLEDGE_INFO *si)
{
  int i, j;
  int triggered;
  int d;
  int i2,j2;
  int last_j2;
  int nc;
  int saveL     = om->L;	                          /* Save the length config of <om>; will restore upon return                               */
  int save_mode = om->mode;	                        /* Likewise for the mode.                                                                 */
  int status;
  
  if ((status = p7_domaindef_GrowTo(ddef, sq->n))      != eslOK) return eslENORESULT;  /* ddef's btot,etot,mocc now ready for seq of length n */
  if ((status = p7_DomainDecoding(om, oxf, oxb, ddef)) != eslOK) return eslENORESULT;  /* ddef->{btot,etot,mocc} now made.                    */

  esl_vec_FSet(ddef->n2sc, sq->n+1, 0.0);          /* ddef->n2sc null2 scores are initialized                                                 */
  ddef->nexpected = ddef->btot[sq->n];             /* posterior expectation for # of domains (same as etot[sq->n])                            */

  p7_oprofile_ReconfigUnihit(om, saveL);	         /* process each domain in unihit mode, regardless of om->mode                              */
  i = -1;
  triggered = FALSE;

  for (j = 1; j <= sq->n; j++)
  {

    if (! triggered)
    {			/* xref J2/101 for what the logic below is: */
      if       (ddef->mocc[j] - (ddef->btot[j] - ddef->btot[j-1]) <  ddef->rt2) i = j;
      else if  (i == -1)                                                        i = j;
      if       (ddef->mocc[j]                                     >= ddef->rt1) triggered = TRUE;
    }
    else if (ddef->mocc[j] - (ddef->etot[j] - ddef->etot[j-1])  <  ddef->rt2)
    {
        /* We have a region i..j to evaluate. */
        p7_omx_GrowTo(fwd, om->M, j-i+1, j-i+1);
        p7_omx_GrowTo(bck, om->M, j-i+1, j-i+1);
        ddef->nregions++;
        if (is_multidomain_region(ddef, i, j))
        {
            /* This region appears to contain more than one domain, so we have to
             * resolve it by cluster analysis of posterior trace samples, to define
             * one or more domain envelopes.
             */
            ddef->nclustered++;

            /* Resolve the region into domains by stochastic trace
             * clustering; assign position-specific null2 model by
             * stochastic trace clustering; there is redundancy
             * here; we will consolidate later if null2 strategy
             * works
             */
            p7_oprofile_ReconfigMultihit(om, saveL);
            p7_Forward(sq->dsq+i-1, j-i+1, om, fwd, NULL);

            region_trace_ensemble(ddef, om, sq->dsq, i, j, fwd, bck, &nc);
            p7_oprofile_ReconfigUnihit(om, saveL);
            /* ddef->n2sc is now set on i..j by the traceback-dependent method */

            last_j2 = 0;
            for (d = 0; d < nc; d++)
            {
                  p7_spensemble_GetClusterCoords(ddef->sp, d, &i2, &j2, NULL, NULL, NULL);
                  if (i2 <= last_j2) ddef->noverlaps++;

                  /* Note that k..m coords on model are available, but
                     * we're currently ignoring them.  This leads to a
                     * rare clustering bug that we eventually need to fix
                     * properly [xref J3/32]: two different regions in one
                     * profile HMM might have hit same seq domain, and
                     * when we now go to calculate an OA trace, nothing
                     * constrains us to find the two different alignments
                     * to the HMM; in fact, because OA is optimal, we'll
                     * find one and the *same* alignment, leading to an
                     * apparent duplicate alignment in the output.
                     *
                     * Registered as #h74, Dec 2009, after EBI finds and
                     * reports it.  #h74 is worked around in p7_tophits.c
                     * by hiding all but one envelope with an identical
                     * alignment, in the rare event that this
                     * happens. [xref J5/130].
                  */
                  ddef->nenvelopes++;

                  /*the !long_target argument will cause the function to recompute null2
                   * scores if this is part of a long_target (nhmmer) pipeline */
                  if ((status = rescore_domain_cust(ddef, om, sq, qsq, fwd, bck, i2, j2, TRUE, bg, long_target, bg_tmp, scores_arr, fwd_emissions_arr, result, nullsc, si)) == eslOK)
                    last_j2 = j2;
                  if ((result->decision == REJECT) && (!(si->filter))) break; // breaking out of first loop over domains if early stopping

            }

            p7_spensemble_Reuse(ddef->sp);
            p7_trace_Reuse(ddef->tr);
            if ((result->decision == REJECT) && (!(si->filter))) break; // must break out of second loop if early stopping!!
        }
        else
        {
            /* The region looks simple, single domain; convert the region to an envelope. */
            ddef->nenvelopes++;
            status = rescore_domain_cust(ddef, om, sq, qsq, fwd, bck, i, j, FALSE, bg, long_target, bg_tmp, scores_arr, fwd_emissions_arr, result, nullsc, si);
            if ((result->decision == REJECT) && (!(si->filter))) break; // break from single domain case if rejected
        }
        i = -1;
        triggered = FALSE;
    }
  }

  /* Restore model to uni/multihit mode, and to its original length model */
  if (p7_IsMulti(save_mode)) p7_oprofile_ReconfigMultihit(om, saveL); 
  else                       p7_oprofile_ReconfigUnihit  (om, saveL); 

  return status;
}

static int
is_multidomain_region(P7_DOMAINDEF *ddef, int i, int j)
{
  int   z;
  float max;
  float expected_n;

  max = -1.0;
  for (z = i; z <= j; z++)
    {
      expected_n = ESL_MIN( (ddef->etot[z] - ddef->etot[i-1]), (ddef->btot[j] - ddef->btot[z-1]));
      max        = ESL_MAX(max, expected_n);
    }

  return ( (max >= ddef->rt3) ? TRUE : FALSE);
}

static int
region_trace_ensemble(P7_DOMAINDEF *ddef, const P7_OPROFILE *om, const ESL_DSQ *dsq, int ireg, int jreg, 
		      const P7_OMX *fwd, P7_OMX *wrk, int *ret_nc)
{
  int    Lr  = jreg-ireg+1;
  int    t, d, d2;
  int    nov, n;
  int    nc;
  int    pos;
  float  null2[p7_MAXCODE];

  esl_vec_FSet(ddef->n2sc+ireg, Lr, 0.0); /* zero the null2 scores in region */

  /* By default, we make results reproducible by forcing a reset of
   * the RNG to its originally seeded state.
   */
  if (ddef->do_reseeding) 
    esl_randomness_Init(ddef->r, esl_randomness_GetSeed(ddef->r));

  /* Collect an ensemble of sampled traces; calculate null2 odds ratios from these */
  for (t = 0; t < ddef->nsamples; t++)
    {
      p7_StochasticTrace(ddef->r, dsq+ireg-1, Lr, om, fwd, ddef->tr);
      p7_trace_Index(ddef->tr);

      pos = 1;
      for (d = 0; d < ddef->tr->ndom; d++)
	{
	  p7_spensemble_Add(ddef->sp, t, ddef->tr->sqfrom[d]+ireg-1, ddef->tr->sqto[d]+ireg-1, ddef->tr->hmmfrom[d], ddef->tr->hmmto[d]);

	  p7_Null2_ByTrace(om, ddef->tr, ddef->tr->tfrom[d], ddef->tr->tto[d], wrk, null2);
	  
	  /* residues outside domains get bumped +1: because f'(x) = f(x), so f'(x)/f(x) = 1 in these segments */
	  for (; pos <= ddef->tr->sqfrom[d]; pos++) ddef->n2sc[ireg+pos-1] += 1.0;

	  /* Residues inside domains get bumped by their null2 ratio */
	  for (; pos <= ddef->tr->sqto[d];   pos++) ddef->n2sc[ireg+pos-1] += null2[dsq[ireg+pos-1]];
	}
      /* the remaining residues in the region outside any domains get +1 */
      for (; pos <= Lr; pos++)  ddef->n2sc[ireg+pos-1] += 1.0;

      p7_trace_Reuse(ddef->tr);        
    }

  /* Convert the accumulated n2sc[] ratios in this region to log odds null2 scores on each residue. */
  for (pos = ireg; pos <= jreg; pos++)
    ddef->n2sc[pos] = logf(ddef->n2sc[pos] / (float) ddef->nsamples);

  /* Cluster the ensemble of traces to break region into envelopes. */
  p7_spensemble_Cluster(ddef->sp, ddef->min_overlap, ddef->of_smaller, ddef->max_diagdiff, ddef->min_posterior, ddef->min_endpointp, &nc);

  /* A little hacky now. Remove "dominated" domains relative to seq coords. */
  for (d = 0; d < nc; d++) 
    ddef->sp->assignment[d] = 0; /* overload <assignment> to flag that a domain is dominated */

  /* who dominates who? (by post prob) */
  for (d = 0; d < nc; d++)
    {
      for (d2 = d+1; d2 < nc; d2++)
	{
	  nov = ESL_MIN(ddef->sp->sigc[d].j, ddef->sp->sigc[d2].j) - ESL_MAX(ddef->sp->sigc[d].i, ddef->sp->sigc[d2].i) + 1;
	  if (nov == 0) break;
	  n   = ESL_MIN(ddef->sp->sigc[d].j - ddef->sp->sigc[d].i + 1,  ddef->sp->sigc[d2].j - ddef->sp->sigc[d2].i + 1);
	  if ((float) nov / (float) n >= 0.8) /* overlap */
	    {
	      if (ddef->sp->sigc[d].prob > ddef->sp->sigc[d2].prob) ddef->sp->assignment[d2] = 1;
	      else                                                  ddef->sp->assignment[d]  = 1;
	    }
	}
    }
      
  /* shrink the sigc list, removing dominated domains */
  d = 0;
  for (d2 = 0; d2 < nc; d2++)
    {
      if (ddef->sp->assignment[d2]) continue; /* skip domain d2, it's dominated. */
      if (d != d2) memcpy(ddef->sp->sigc + d, ddef->sp->sigc + d2, sizeof(struct p7_spcoord_s));
      d++;
    }
  ddef->sp->nc = d;
  *ret_nc = d;
  return eslOK;
}

static int
rescore_domain_cust(P7_DOMAINDEF *ddef, P7_OPROFILE *om, const ESL_SQ *sq, const ESL_SQ *qsq,
			P7_OMX *ox1, P7_OMX *ox2, int i, int j, int null2_is_done, P7_BG *bg, int long_target,
			P7_BG *bg_tmp, float *scores_arr, float *fwd_emissions_arr, 
      RESULT_INFO *result, float *nullsc, SLEDGE_INFO *si)
{
  P7_DOMAIN     *dom           = NULL;
  int            Ld            = j-i+1;
  float          domcorrection = 0.0;
  float          envsc, oasc;
  int            z;
  int            pos;
  float          null2[p7_MAXCODE];
  int            status;
  int            max_env_extra = 20;
  int            orig_L;
  float          eval;
  float          pid;
  char           decision;

  
  p7_Forward (sq->dsq + i-1, Ld, om, ox1, &envsc);
  p7_Backward(sq->dsq + i-1, Ld, om, ox1, ox2, NULL);

  status = p7_Decoding(om, ox1, ox2, ox2);      /* <ox2> is now overwritten with post probabilities     */
  if (status == eslERANGE) { /* rare: numeric overflow; domain is assumed to be repetitive garbage [J3/119-121] */
    goto ERROR;
  }

  /* Find an optimal accuracy alignment */
  p7_OptimalAccuracy(om, ox2, ox1, &oasc);      /* <ox1> is now overwritten with OA scores              */
  p7_OATrace        (om, ox2, ox1, ddef->tr);   /* <tr>'s seq coords are offset by i-1, rel to orig dsq */

  /* hack the trace's sq coords to be correct w.r.t. original dsq */
  for (z = 0; z < ddef->tr->N; z++)
    if (ddef->tr->i[z] > 0) ddef->tr->i[z] += i-1;

  /* Find the percent identity and accept/reject the sequence based on threshold */
  find_pid(ddef->tr, sq, qsq, &pid, &decision, si);
  
  /////////// IS THIS BOXED SECTION NEEDED? REMOVE, TEST AND UPDATE ///////////
  /* get ptr to next empty domain structure in domaindef's results */
  if (ddef->ndom == ddef->nalloc) {
    ESL_REALLOC(ddef->dcl, sizeof(P7_DOMAIN) * (ddef->nalloc*2));
    ddef->nalloc *= 2;
  }
  /////////////////////////////////////////////////////////////////////////////

  if (decision == REJECT)
  {
    dom = &(ddef->dcl[ddef->ndom]);
    dom->ad             = p7_alidisplay_Create(ddef->tr, 0, om, sq, qsq);
    dom->scores_per_pos = NULL;

    if (!null2_is_done) {
      p7_Null2_ByExpectation(om, ox2, null2);
      for (pos = i; pos <= j; pos++)
        ddef->n2sc[pos]  = logf(null2[sq->dsq[pos]]);
    }
    for (pos = i; pos <= j; pos++)
      domcorrection += ddef->n2sc[pos]; /* domcorrection is in units of NATS */

    dom->domcorrection = domcorrection; /* in units of NATS */
    dom->ienv          = i;
    dom->jenv          = j;
    dom->envsc         = envsc;         /* in units of NATS */
    dom->oasc          = oasc;          /* in units of expected # of correctly aligned residues */
    dom->dombias       = 0.0;           /* gets set later, using bg->omega and dombias */
    dom->bitscore      = 0.0;           /* gets set later by caller, using envsc, null score, and dombias */
    dom->lnP           = 0.0;           /* gets set later by caller, using bitscore */
    dom->is_reported   = FALSE;         /* gets set later by caller */
    dom->is_included   = FALSE;         /* gets set later by caller */

    ddef->ndom++;
    
    /* E-value calculation */
    Ld = dom->jenv - dom->ienv + 1;
    dom->bitscore = dom->envsc + (sq->n-Ld) * log((float) sq->n / (float) (sq->n+3));                 /* NATS, for the moment...   */
    dom->dombias  = p7_FLogsum(0.0, log(bg->omega) + dom->domcorrection);                             /* NATS, and will stay so    */
    dom->bitscore = (dom->bitscore - (*nullsc + dom->dombias)) / eslCONST_LOG2;                       /* now BITS, as it should be */
    dom->lnP      = esl_exp_logsurv (dom->bitscore,  om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);  /* p-value in log space      */
    eval = dom->lnP + log(si->db_size);                                                               /* E-value based on db size  */

    if (eval < log(si->incE))                                                                         /* E-value threshold check   */
      append_target(result, sq->name, pid, eval, ddef->tr->N);
  }
  
  p7_trace_Reuse(ddef->tr);

  return eslOK;

 ERROR:
  p7_trace_Reuse(ddef->tr);
  return status;
}

static int
find_pid(P7_TRACE *tr, const ESL_SQ *sq, const ESL_SQ *qsq, float *pid, char *decision, SLEDGE_INFO *si)
{
  int j;
  double match = 0.0;

  double q_length = 0;
  double t_length = 0;

  /* Walk through the trace, check for match state and whether canonical residue */
  for(j=0; j<tr->N; j++)
  {
    if ((tr->st[j] < p7T_M) || (tr->st[j] > p7T_I)) continue;
    else 
    {
      if ((tr->st[j] == p7T_M) && ((qsq->dsq[tr->k[j]]) == (sq->dsq[tr->i[j]])) && 
            (qsq->dsq[tr->k[j]] < 20) && (sq->dsq[tr->i[j]] < 20)) { match++; q_length++;  t_length++; }
      else if (tr->st[j] == p7T_D) t_length++;
      else q_length++;

    }
  }
  
  /* Store the pid and decision for this domain comparison */
  *pid = match / q_length;
  *decision = ((*pid >= si->pid_low) && (*pid < si->pid_high)) ? ACCEPT : REJECT;

  return eslOK;
}