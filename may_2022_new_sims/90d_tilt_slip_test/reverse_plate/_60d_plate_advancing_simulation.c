#line 1 "60d_plate_advancing_simulation-cpp.c"
#line 1 "<built-in>"
#line 1 "<command-line>"
#line 1 "/usr/include/stdc-predef.h"
#line 1 "<command-line>"
#line 1 "60d_plate_advancing_simulation-cpp.c"
#if _XOPEN_SOURCE < 700
#undef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif
#if _GNU_SOURCE
#include <stdint.h>
#include <string.h>
#include <fenv.h>
#endif



#line 1 "/home/fpl/softwares/basilisk/src/common.h"
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#ifndef assert
# include <assert.h>
#endif
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#if _OPENMP
# include <omp.h>
# define OMP(x) _Pragma(#x)
#elif 1

# define OMP(x)

# include <mpi.h>
static int mpi_rank, mpi_npe;
# define tid() mpi_rank
# define pid() mpi_rank
# define npe() mpi_npe

#else

# define OMP(x)

#endif

#if _CADNA
# include <cadna.h>
#endif

#if __cplusplus
# define delete delete_qcc
# define right right_qcc
# define left left_qcc
# define norm norm_qcc
# define new new_qcc
#endif

#define pi 3.14159265358979
#undef HUGE
#define HUGE ((double)1e30)
#define nodata HUGE
#define _NVARMAX 65536
#define is_constant(v) ((v).i >= _NVARMAX)
#define constant(v) (is_constant(v) ? _constant[(v).i - _NVARMAX] : nodata)

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))
#define sq(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))
#define sign(x) ((x) > 0 ? 1 : -1)
#define noise() (1. - 2.*rand()/(double)RAND_MAX)
#define clamp(x,a,b) ((x) < (a) ? (a) : (x) > (b) ? (b) : (x))
#define swap(type,a,b) do { type __tmp = a; a = b; b = __tmp; } while(0)
#define unmap(x,y)

#define trash(x)


#define systderr stderr
#define systdout stdout
#define sysfprintf fprintf

#if 1
FILE * qstderr (void);
FILE * qstdout (void);
FILE * ferr = NULL, * fout = NULL;
#define not_mpi_compatible()\
do {\
  if (npe() > 1) {\
    fprintf (ferr, "%s() is not compatible with MPI (yet)\n", __func__);\
    exit (1);\
  }\
} while(0)\

#line 84

# define system(command) (pid() == 0 ? system(command) : 0)
#else
# define qstderr() stderr
# define qstdout() stdout
# define ferr stderr
# define fout stdout
# define not_mpi_compatible()
#endif



static inline void qassert (const char * file, int line, const char * cond) {
  fprintf (ferr, "%s:%d: Assertion `%s' failed.\n", file, line, cond);
  abort();
}
#line 108 "/home/fpl/softwares/basilisk/src/common.h"
#define sysmalloc malloc
#define syscalloc calloc
#define sysrealloc realloc
#define sysfree free
#define systrdup strdup

#if MTRACE

struct {
  FILE * fp;
  size_t total, max;
  size_t overhead, maxoverhead;
  size_t nr;
  size_t startrss, maxrss;
  char * fname;
} pmtrace;

typedef struct {
  char * func, * file;
  size_t max, total;
  int line, id;
} pmfunc;

typedef struct {
  size_t id, size;
} pmdata;

static pmfunc * pmfuncs = NULL;
static int pmfuncn = 0;

static int pmfunc_index (const char * func, const char * file, int line)
{
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn; i++, p++)
    if (p->line == line && !strcmp(func, p->func) && !strcmp(file, p->file))
      return p->id;
  pmfuncn++;
  pmfuncs = (pmfunc *) sysrealloc (pmfuncs, pmfuncn*sizeof(pmfunc));
  p = &pmfuncs[pmfuncn - 1];
  memset (p, 0, sizeof(pmfunc));
  p->func = systrdup(func);
  p->file = systrdup(file);
  p->line = line;
  p->id = pmfuncn;
  if (pmtrace.fp)
    fprintf (pmtrace.fp, "@ %d %s %s %d\n", pmfuncn, func, file, line);
  return pmfuncn;
}

static void pmfunc_trace (pmfunc * f, char c)
{
  if (pmtrace.fp)
    fprintf (pmtrace.fp, "%c %d %ld %ld %ld",
      c, f->id, pmtrace.nr, pmtrace.total, f->total);
#if _GNU_SOURCE
  if (pmtrace.nr % 1 == 0) {
    struct rusage usage;
    getrusage (RUSAGE_SELF, &usage);
    if (pmtrace.fp)
      fprintf (pmtrace.fp, " %ld", usage.ru_maxrss*1024);
    if (!pmtrace.nr)
      pmtrace.startrss = usage.ru_maxrss;
    if (usage.ru_maxrss > pmtrace.maxrss)
      pmtrace.maxrss = usage.ru_maxrss;
  }
#endif
  if (pmtrace.fp)
    fputc ('\n', pmtrace.fp);
  pmtrace.nr++;
}

static void * pmfunc_alloc (pmdata * d, size_t size,
       const char * func, const char * file, int line,
       char c)
{
  if (!(d != NULL)) qassert ("/home/fpl/softwares/basilisk/src/common.h", 183, "d != NULL");
  OMP (omp critical)
  {
    d->id = pmfunc_index(func, file, line);
    d->size = size;
    pmfunc * f = &pmfuncs[d->id - 1];
    f->total += size;
    if (f->total > f->max)
      f->max = f->total;
    pmtrace.total += size;
    pmtrace.overhead += sizeof(pmdata);
    if (pmtrace.total > pmtrace.max) {
      pmtrace.max = pmtrace.total;
      pmtrace.maxoverhead = pmtrace.overhead;
    }
    pmfunc_trace (f, c);
  }
  return ((char *)d) + sizeof(pmdata);
}

static void * pmfunc_free (void * ptr, char c)
{
  if (!ptr)
    return ptr;
  pmdata * d = (pmdata *) (((char *)ptr) - sizeof(pmdata));
  if (d->id < 1 || d->id > pmfuncn) {
    fputs ("*** MTRACE: ERROR!: corrupted free()", ferr);
    if (d->size == 0)
      fputs (", possible double free()", ferr);
    else
      fputs (", not traced?", ferr);
    fputs (", aborting...\n", ferr);
    abort();
    return ptr;
  }
  else
  OMP (omp critical)
  {
    pmfunc * f = &pmfuncs[d->id - 1];
    if (f->total < d->size) {
      fprintf (ferr, "*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\n",
        f->total, d->size);
      abort();
    }
    else
      f->total -= d->size;
    if (pmtrace.total < d->size) {
      fprintf (ferr, "*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\n",
        pmtrace.total, d->size);
      abort();
    }
    else {
      pmtrace.total -= d->size;
      pmtrace.overhead -= sizeof(pmdata);
    }
    d->id = 0;
    d->size = 0;
    pmfunc_trace (f, c);
  }
  return d;
}

static void * pmalloc (size_t size,
         const char * func, const char * file, int line)
{
  return pmfunc_alloc ((pmdata *) sysmalloc (sizeof(pmdata) + size),
         size, func, file, line, '+');
}

static void * pcalloc (size_t nmemb, size_t size,
         const char * func, const char * file, int line)
{
  void * p = pmalloc (nmemb*size, func, file, line);
  return memset (p, 0, nmemb*size);
}

static void * prealloc (void * ptr, size_t size,
   const char * func, const char * file, int line)
{
  return pmfunc_alloc ((pmdata *) sysrealloc (pmfunc_free(ptr, '<'),
           sizeof(pmdata) + size),
         size, func, file, line, '>');
}

static void pfree (void * ptr,
     const char * func, const char * file, int line)
{
  sysfree (pmfunc_free (ptr, '-'));
}

static char * pstrdup (const char * s,
         const char * func, const char * file, int line)
{
  char * d = (char *) pmalloc (strlen(s) + 1, func, file, line);
  return strcpy (d, s);
}

#if MTRACE < 3
static int pmaxsort (const void * a, const void * b) {
  const pmfunc * p1 = a, * p2 = b;
  return p1->max < p2->max;
}
#endif

static int ptotalsort (const void * a, const void * b) {
  const pmfunc * p1 = (const pmfunc *) a, * p2 = (const pmfunc *) b;
  return p1->total < p2->total;
}

static void pmfuncs_free ()
{
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn; i++, p++) {
    sysfree (p->func);
    sysfree (p->file);
  }
  sysfree (pmfuncs);
}

void pmuntrace (void)
{
#if MTRACE < 3
  fprintf (ferr,
    "*** MTRACE: max resident  set size: %10ld bytes\n"
    "*** MTRACE: max traced memory size: %10ld bytes"
    " (tracing overhead %.1g%%)\n"
    "%10s    %20s   %s\n",
    pmtrace.maxrss*1024,
    pmtrace.max, pmtrace.maxoverhead*100./pmtrace.max,
    "max bytes", "function", "file");
  qsort (pmfuncs, pmfuncn, sizeof(pmfunc), pmaxsort);
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn && p->max > 0; i++, p++)
    fprintf (ferr, "%10ld    %20s   %s:%d\n",
      p->max, p->func, p->file, p->line);

  if (pmtrace.fp) {
    char * fname = pmtrace.fname, * s;
    while ((s = strchr(fname,'/')))
      fname = s + 1;

    fputs ("load(\"`echo $BASILISK`/mtrace.plot\")\n", pmtrace.fp);
    fprintf (pmtrace.fp,
      "plot '%s' u 3:($6-%g) w l t 'ru_maxrss - %.3g',"
      "total(\"%s\") w l t 'total'",
      fname,
      pmtrace.startrss*1024.,
      pmtrace.startrss*1024.,
      fname);
    pmfunc * p = pmfuncs;
    for (int i = 0; i < pmfuncn && p->max > 0.01*pmtrace.max; i++, p++)
      fprintf (pmtrace.fp,
        ",func(\"%s\",%d) w l t '%s'",
        fname, p->id, p->func);
    fputc ('\n', pmtrace.fp);
    fprintf (ferr,
      "*** MTRACE: To get a graph use: tail -n 2 %s | gnuplot -persist\n",
      fname);
    fclose (pmtrace.fp);
    pmtrace.fp = NULL;
    sysfree (pmtrace.fname);
  }
#endif

  if (pmtrace.total > 0) {
    qsort (pmfuncs, pmfuncn, sizeof(pmfunc), ptotalsort);
    pmfunc * p = pmfuncs;
    for (int i = 0; i < pmfuncn && p->total > 0; i++, p++)
      fprintf (ferr, "%s:%d: error: %ld bytes leaked here\n",
        p->file, p->line, p->total);
    pmfuncs_free();
    exit(1);
  }
  else {
#if MTRACE < 3
    fputs ("*** MTRACE: No memory leaks\n", ferr);
#endif
    pmfuncs_free();
  }
}

#else
# define pmalloc(s,func,file,line) malloc(s)
# define pcalloc(n,s,func,file,line) calloc(n,s)
# define prealloc(p,s,func,file,line) realloc(p,s)
# define pfree(p,func,file,line) free(p)
# define pstrdup(s,func,file,line) strdup(s)
#endif







typedef struct {
  void * p;
  long max, len;
} Array;

Array * array_new()
{
  Array * a = ((Array *) pmalloc ((1)*sizeof(Array),__func__,__FILE__,__LINE__));
  a->p = NULL;
  a->max = a->len = 0;
  return a;
}

void array_free (Array * a)
{
  pfree (a->p,__func__,__FILE__,__LINE__);
  pfree (a,__func__,__FILE__,__LINE__);
}

void array_append (Array * a, void * elem, size_t size)
{
  if (a->len + size >= a->max) {
    a->max += max (size, 4096);
    a->p = prealloc (a->p, a->max,__func__,__FILE__,__LINE__);
  }
  memcpy (((char *)a->p) + a->len, elem, size);
  a->len += size;
}

void * array_shrink (Array * a)
{
  void * p = prealloc (a->p, a->len,__func__,__FILE__,__LINE__);
  pfree (a,__func__,__FILE__,__LINE__);
  return p;
}



#if TRACE == 1
#include <extrae_user_events.h>

typedef struct {
  Array index, stack;
  extrae_type_t type;
} Trace;

Trace trace_func = {
  {NULL, 0, 0}, {NULL, 0, 0},
  60000010,
};

Trace trace_mpi_func = {
  {NULL, 0, 0}, {NULL, 0, 0},
  60000011,
};

static int lookup_func (Array * a, const char * func)
{
  for (int i = 0; i < a->len/sizeof(char *); i++) {
    char * s = ((char **)a->p)[i];
    if (!strcmp (func, s))
      return i + 1;
  }
  char * s = pstrdup (func,__func__,__FILE__,__LINE__);
  array_append (a, &s, sizeof(char *));
  return a->len;
}

static void trace_push (Trace * t, const char * func)
{
  int value = lookup_func (&t->index, func);
  Extrae_eventandcounters (t->type, value);
  array_append (&t->stack, &value, sizeof(int));
}

static void trace_pop (Trace * t, const char * func)
{
  if (!(t->stack.len > 0)) qassert ("/home/fpl/softwares/basilisk/src/common.h", 455, "t->stack.len > 0");
  t->stack.len -= sizeof(int);
  int value = t->stack.len > 0 ?
    ((int *)t->stack.p)[t->stack.len/sizeof(int) - 1] : 0;
  Extrae_eventandcounters (t->type, value);
}

static void trace_define (Trace * t, char * description)
{
  if (t->index.len > 0) {
    extrae_value_t values[t->index.len/sizeof(char *) + 1];
    char * names[t->index.len/sizeof(char *) + 1],
      ** func = (char **) t->index.p;
    names[0] = "OTHER";
    values[0] = 0;
    unsigned len = 1;
    for (int i = 0; i < t->index.len/sizeof(char *); i++, func++) {
      names[len] = *func;
      values[len++] = i + 1;
    }
    Extrae_define_event_type (&t->type, description, &len, values, names);
  }
}

static void trace_free (Trace * t)
{
  char ** func = (char **) t->index.p;
  for (int i = 0; i < t->index.len/sizeof(char *); i++, func++)
    pfree (*func,__func__,__FILE__,__LINE__);
  pfree (t->index.p,__func__,__FILE__,__LINE__);
  pfree (t->stack.p,__func__,__FILE__,__LINE__);
}

static void trace_off ()
{
  trace_define (&trace_func, "Basilisk functions");
  trace_define (&trace_mpi_func, "Basilisk functions (MPI-related)");
  trace_free (&trace_func);
  trace_free (&trace_mpi_func);
}






# define trace(func, file, line) trace_push (&trace_func, func)
# define end_trace(func, file, line) trace_pop (&trace_func, func)

#elif TRACE

typedef struct {
  char * func, * file;
  int line, calls;
  double total, self;
#if 1
  double min, max;
#endif
} TraceIndex;

struct {
  Array stack, index;
  double t0;
} Trace = {
  {NULL, 0, 0}, {NULL, 0, 0},
  -1
};

static void trace_add (const char * func, const char * file, int line,
         double total, double self)
{
  TraceIndex * t = (TraceIndex *) Trace.index.p;
  int i, len = Trace.index.len/sizeof(TraceIndex);
  for (i = 0; i < len; i++, t++)
    if (t->line == line && !strcmp (func, t->func) && !strcmp (file, t->file))
      break;
  if (i == len) {
    TraceIndex t = {pstrdup(func,__func__,__FILE__,__LINE__), pstrdup(file,__func__,__FILE__,__LINE__), line, 1, total, self};
    array_append (&Trace.index, &t, sizeof(TraceIndex));
  }
  else
    t->calls++, t->total += total, t->self += self;
}

static void trace (const char * func, const char * file, int line)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  if (Trace.t0 < 0)
    Trace.t0 = tv.tv_sec + tv.tv_usec/1e6;
  double t[2] = {(tv.tv_sec - Trace.t0) + tv.tv_usec/1e6, 0.};
  array_append (&Trace.stack, t, 2*sizeof(double));




}

static void end_trace (const char * func, const char * file, int line)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  double te = (tv.tv_sec - Trace.t0) + tv.tv_usec/1e6;
  double * t = (double *) Trace.stack.p;
  if (!(Trace.stack.len >= 2*sizeof(double))) qassert ("/home/fpl/softwares/basilisk/src/common.h", 559, "Trace.stack.len >= 2*sizeof(double)");
  t += Trace.stack.len/sizeof(double) - 2;
  Trace.stack.len -= 2*sizeof(double);
  double dt = te - t[0];




  trace_add (func, file, line, dt, dt - t[1]);
  if (Trace.stack.len >= 2*sizeof(double)) {
    t -= 2;
    t[1] += dt;
  }
}

static int compar_self (const void * p1, const void * p2)
{
  const TraceIndex * t1 = p1, * t2 = p2;
  return t1->self < t2->self;
}

#if 1
static int compar_func (const void * p1, const void * p2)
{
  const TraceIndex * t1 = p1, * t2 = p2;
  if (t1->line != t2->line)
    return t1->line < t2->line;
  return strcmp (t1->file, t2->file);
}
#endif

void trace_print (FILE * fp, double threshold)
{
  int i, len = Trace.index.len/sizeof(TraceIndex);
  double total = 0.;
  TraceIndex * t;
  Array * index = array_new();
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    array_append (index, t, sizeof(TraceIndex)), total += t->self;
#if 1
  qsort (index->p, len, sizeof(TraceIndex), compar_func);
  double tot[len], self[len], min[len], max[len];
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    tot[i] = t->total, self[i] = t->self;
  MPI_Reduce (self, min, len, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce (self, max, len, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce (pid() ? self : MPI_IN_PLACE,
       self, len, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (pid() ? tot : MPI_IN_PLACE,
       tot, len, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  total = 0.;
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    t->total = tot[i]/npe(), t->self = self[i]/npe(),
      t->max = max[i], t->min = min[i], total += t->self;
#endif
  qsort (index->p, len, sizeof(TraceIndex), compar_self);
  fprintf (fp, "   calls    total     self   %% total   function\n");
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    if (t->self*100./total > threshold) {
      fprintf (fp, "%8d   %6.2f   %6.2f     %4.1f%%",
        t->calls, t->total, t->self, t->self*100./total);
#if 1
      fprintf (fp, " (%4.1f%% - %4.1f%%)", t->min*100./total, t->max*100./total);
#endif
      fprintf (fp, "   %s():%s:%d\n", t->func, t->file, t->line);
    }
  fflush (fp);
  array_free (index);
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    t->calls = t->total = t->self = 0.;
}

static void trace_off ()
{
  trace_print (fout, 0.);

  int i, len = Trace.index.len/sizeof(TraceIndex);
  TraceIndex * t;
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    pfree (t->func,__func__,__FILE__,__LINE__), pfree (t->file,__func__,__FILE__,__LINE__);

  pfree (Trace.index.p,__func__,__FILE__,__LINE__);
  Trace.index.p = NULL;
  Trace.index.len = Trace.index.max = 0;

  pfree (Trace.stack.p,__func__,__FILE__,__LINE__);
  Trace.stack.p = NULL;
  Trace.stack.len = Trace.stack.max = 0;
}

#else
# define trace(...)
# define end_trace(...)
#endif



#if _OPENMP

#define tid() omp_get_thread_num()
#define pid() 0
#define npe() omp_get_num_threads()
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_array(v,type,op,elem)

#elif 1

static bool in_prof = false;
static double prof_start, _prof;
#define prof_start(name)\
  if (!(!in_prof)) qassert ("/home/fpl/softwares/basilisk/src/common.h", 669, "!in_prof"); in_prof = true;\
  prof_start = MPI_Wtime();\

#line 671

#define prof_stop()\
  if (!(in_prof)) qassert ("/home/fpl/softwares/basilisk/src/common.h", 673, "in_prof"); in_prof = false;\
  _prof = MPI_Wtime();\
  mpi_time += _prof - prof_start;\

#line 676


#if FAKE_MPI
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_array(v,type,op,elem)
#else

int mpi_all_reduce0 (void *sendbuf, void *recvbuf, int count,
       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{ trace ("mpi_all_reduce0", "/home/fpl/softwares/basilisk/src/common.h", 685);
  { int _ret =  MPI_Allreduce (sendbuf, recvbuf, count, datatype, op, comm); end_trace("mpi_all_reduce0", "/home/fpl/softwares/basilisk/src/common.h", 686);  return _ret; }
 end_trace("mpi_all_reduce0", "/home/fpl/softwares/basilisk/src/common.h", 687); }
#define mpi_all_reduce(v,type,op) {\
  prof_start ("mpi_all_reduce");\
  union { int a; float b; double c;} global;\
  mpi_all_reduce0 (&(v), &global, 1, type, op, MPI_COMM_WORLD);\
  memcpy (&(v), &global, sizeof (v));\
  prof_stop();\
}\

#line 695

#define mpi_all_reduce_array(v,type,op,elem) {\
  prof_start ("mpi_all_reduce");\
  type global[elem], tmp[elem];\
  for (int i = 0; i < elem; i++)\
    tmp[i] = (v)[i];\
  MPI_Datatype datatype;\
  if (!strcmp(#type, "double")) datatype = MPI_DOUBLE;\
  else if (!strcmp(#type, "int")) datatype = MPI_INT;\
  else if (!strcmp(#type, "long")) datatype = MPI_LONG;\
  else if (!strcmp(#type, "bool")) datatype = MPI_C_BOOL;\
  else {\
    fprintf (ferr, "unknown reduction type '%s'\n", #type);\
    fflush (ferr);\
    abort();\
  }\
  mpi_all_reduce0 (tmp, global, elem, datatype, op, MPI_COMM_WORLD);\
  for (int i = 0; i < elem; i++)\
    (v)[i] = global[i];\
  prof_stop();\
}\

#line 716


#endif

#define QFILE FILE

FILE * qstderr (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "log-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systderr;
  }
  return fp;
}

FILE * qstdout (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "out-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systdout;
  }
  return fp;
}

static void finalize (void)
{
  MPI_Finalize();
}

void mpi_init ()
{
  int initialized;
  MPI_Initialized (&initialized);
  if (!initialized) {
    MPI_Init (NULL, NULL);
    MPI_Comm_set_errhandler (MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
    atexit (finalize);
  }
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_npe);
  srand (mpi_rank + 1);
  if (ferr == NULL) {
    if (mpi_rank > 0) {
      ferr = fopen ("/dev/null", "w");
      fout = fopen ("/dev/null", "w");
    }
    else {
      ferr = systderr;
      fout = systdout;
    }
    char * etrace = getenv ("MALLOC_TRACE"), name[80];
    if (etrace && mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      setenv ("MALLOC_TRACE", name, 1);
    }
#if MTRACE == 1
    etrace = getenv ("MTRACE");
    if (!etrace)
      etrace = "mtrace";
    if (mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      pmtrace.fp = fopen (name, "w");
      pmtrace.fname = systrdup(name);
    }
    else {
      pmtrace.fp = fopen (etrace, "w");
      pmtrace.fname = systrdup(etrace);
    }
#endif
  }
}

#else

#define tid() 0
#define pid() 0
#define npe() 1
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_array(v,type,op,elem)

#endif

#define OMP_PARALLEL() OMP(omp parallel)

#define NOT_UNUSED(x) (void)(x)

#define VARIABLES ;
#define _index(a,m) (a.i)
#define val(a,k,l,m) data(k,l,m)[_index(a,m)]

double _val_higher_dimension = 0.;
#define _val_higher_dimension(x,a,b,c) _val_higher_dimension
#line 828 "/home/fpl/softwares/basilisk/src/common.h"
#if (_GNU_SOURCE || __APPLE__) && !_OPENMP && !_CADNA
double undefined;
# if __APPLE__
# include <stdint.h>
# include "fp_osx.h"
# endif
# define enable_fpe(flags) feenableexcept (flags)
# define disable_fpe(flags) fedisableexcept (flags)
static void set_fpe (void) {
  int64_t lnan = 0x7ff0000000000001;
  if (!(sizeof (int64_t) == sizeof (double))) qassert ("/home/fpl/softwares/basilisk/src/common.h", 838, "sizeof (int64_t) == sizeof (double)");
  memcpy (&undefined, &lnan, sizeof (double));
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}
#else
# define undefined ((double) DBL_MAX)
# define enable_fpe(flags)
# define disable_fpe(flags)
static void set_fpe (void) {}
#endif


typedef struct {
  long n;
  long tn;
  int depth;
  int maxdepth;
} Grid;
Grid * grid = NULL;

double X0 = 0., Y0 = 0., Z0 = 0.;

double L0 = 1.;


int N = 64;




typedef struct { int i; } scalar;

typedef struct {
  scalar x;

  scalar y;




} vector;

typedef struct {
  scalar * x;

  scalar * y;




} vectorl;

typedef struct {
  vector x;

  vector y;




} tensor;

struct { int x, y, z; } Period = {false, false, false};

typedef struct {
  double x, y, z;
} coord;

OMP(omp declare reduction (+ : coord :
      omp_out.x += omp_in.x,
      omp_out.y += omp_in.y,
      omp_out.z += omp_in.z))
#line 922 "/home/fpl/softwares/basilisk/src/common.h"
void normalize (coord * n)
{
  double norm = 0.;
  {
#line 925

    norm += sq(n->x);
#line 925

    norm += sq(n->y);}
  norm = sqrt(norm);
  {
#line 928

    n->x /= norm;
#line 928

    n->y /= norm;}
}

struct _origin { double x, y, z; };

void origin (struct _origin p) {
  X0 = p.x; Y0 = p.y; Z0 = p.z;
}

void size (double L) {
  L0 = L;
}

double zero (double s0, double s1, double s2) { return 0.; }






  enum { right, left, top, bottom };



int nboundary = 2*2;



#define dirichlet(expr) (2.*(expr) - val(_s,0,0,0))
#define dirichlet_homogeneous() (- val(_s,0,0,0))
#define neumann(expr) (Delta*(expr) + val(_s,0,0,0))
#define neumann_homogeneous() (val(_s,0,0,0))

double * _constant = NULL;
extern size_t datasize;
typedef struct _Point Point;



#define strongif(x) if(x)
#define IF(x) if((x)||1)

typedef struct { int i[2]; } IJK;

void _stencil_fprintf (const char * file, int line,
         FILE * stream, const char *format, ...) {}
void _stencil_printf (const char * file, int line,
        const char *format, ...) {}
void _stencil_fputc (const char * file, int line,
       int c, FILE * stream) {}
void _stencil_fputs (const char * file, int line,
       const char * s, FILE * stream) {}
#define _stencil_qassert(...) ((void)0)

int _stencil_access (scalar s, IJK i, const char * file, int line);

#if 2 == 1
#define _stencil_val(file, line, s, i1, i2, i3)\
  (_attribute[is_constant(s) || s.i < 0 ? -1 : s.i].write\
   [_stencil_access(s, (IJK){{point.i + i1}}, file, line)])\

#line 989

#elif 2 == 2
#define _stencil_val(file, line, s, i1, i2, i3)\
  (_attribute[is_constant(s) || s.i < 0 ? -1 : s.i].write\
   [_stencil_access(s, (IJK){{point.i + i1, point.j + i2}}, file, line)])\

#line 994

#else
#define _stencil_val(file, line, s, i1, i2, i3)\
  (_attribute[is_constant(s) || s.i < 0 ? -1 : s.i].write\
   [_stencil_access(s, (IJK){{point.i + i1, point.j + i2, point.k + i3}},\
      file, line)])\

#line 1000

#endif




#define _stencil_fine(file, line, s, i1, i2, i3)\
  _stencil_val(file, line, s, i1, i2, i3)\

#line 1008

#define _stencil_coarse(file, line, s, i1, i2, i3)\
  _stencil_val(file, line, s, i1, i2, i3)\

#line 1011


#line 1 "/home/fpl/softwares/basilisk/src/grid/boundaries.h"


typedef struct _Boundary Boundary;

struct _Boundary {
  void (* destroy) (Boundary * b);
  void (* level) (const Boundary * b, scalar * list, int l);

  void (* restriction) (const Boundary * b, scalar * list, int l);
};

static Boundary ** boundaries = NULL;

void add_boundary (Boundary * b) {
  int len = 0;
  if (boundaries) {
    Boundary ** i = boundaries;
    while (*i++) len++;
  }
  boundaries = (Boundary * *) prealloc (boundaries, (len + 2)*sizeof(Boundary *),__func__,__FILE__,__LINE__);
  boundaries[len] = b;
  boundaries[len+1] = NULL;
}

void free_boundaries () {
  if (!boundaries)
    return;
  Boundary ** i = boundaries, * b;
  while ((b = *i++))
    if (b->destroy)
      b->destroy (b);
    else
      pfree (b,__func__,__FILE__,__LINE__);
  pfree (boundaries,__func__,__FILE__,__LINE__);
  boundaries = NULL;
}
#line 47 "/home/fpl/softwares/basilisk/src/grid/boundaries.h"
typedef struct {
  Boundary parent;
  int d;
} BoxBoundary;
#line 1014 "/home/fpl/softwares/basilisk/src/common.h"



typedef struct {

#line 1019 "/home/fpl/softwares/basilisk/src/common.h"

  double (** boundary) (Point, Point, scalar, void *);
  double (** boundary_homogeneous) (Point, Point, scalar, void *);
  double (* gradient) (double, double, double);
  void (* delete) (scalar);
  char * name;
  struct {
    int x;

    int y;




  } d;
  vector v;
  int face;
  bool nodump, freed;
  int block;
  scalar * depends;

#line 18 "/home/fpl/softwares/basilisk/src/grid/stencils.h"

  double * write;
  int * read;
  int dirty;




#line 17 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

  void (* prolongation) (Point, scalar);
  void (* restriction) (Point, scalar);

#line 8 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"

  void (* refine) (Point, scalar);

#line 94 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"

  void (* coarsen) (Point, scalar);

#line 81 "/home/fpl/softwares/basilisk/src/fractions.h"

  vector n;

#line 454 "/home/fpl/softwares/basilisk/src/heights.h"

  vector height;

#line 20 "/home/fpl/softwares/basilisk/src/iforce.h"

  scalar phi;

#line 21 "/home/fpl/softwares/basilisk/src/tension.h"

  double sigma;

#line 27 "/home/fpl/softwares/basilisk/src/vof.h"

  scalar * tracers, c;
  bool inverse;

} _Attributes;
_Attributes * _attribute;
#define _call_vof_concentration_gradient_x 1
#define _call_vof_concentration_gradient_y 1
#define _call_height_position 1
#define _call_interfacial 1
#define _call_centroids_curvature_fit 1
#define _call_height_curvature_fit 1
#define _call_height_curvature 1
#define _call_facet_normal 1
#define _call_mycs 1
#define _call_interface_normal 1
#if _call_vof_concentration_refine
#  undef _call_vof_concentration_gradient_x
#  define _call_vof_concentration_gradient_x 1
#  undef _call_vof_concentration_gradient_y
#  define _call_vof_concentration_gradient_y 1
#endif
#if _call_interface_normal
#  undef _call_mycs
#  define _call_mycs 1
#  undef _call_height_normal
#  define _call_height_normal 1
#endif
#if _call_height_position
#  undef _call_pos_x
#  define _call_pos_x 1
#  undef _call_pos_y
#  define _call_pos_y 1
#endif
#if _call_centroids_curvature_fit
#  undef _call_mycs
#  define _call_mycs 1
#endif
#if _call_height_curvature_fit
#  undef _call_mycs
#  define _call_mycs 1
#endif
#if _call_facet_normal
#  undef _call_interface_normal
#  define _call_interface_normal 1
#endif
#if _call_fraction_refine
#  undef _call_mycs
#  define _call_mycs 1
#endif
#if _call_locals_pids
#  undef _call_is_local_prolongation
#  define _call_is_local_prolongation 1
#endif
#if _call_rcv_pid_append
#  undef _call_rcv_append
#  define _call_rcv_append 1
#endif
#if _call_treey
#  undef _call_treey
#  define _call_treey 1
#endif
#if _call_treex
#  undef _call_treex
#  define _call_treex 1
#endif
#if _call_refine_face_solenoidal
#  undef _call_refine_face
#  define _call_refine_face 1
#endif
#if _call_coarsen_cell_recursive
#  undef _call_coarsen_cell
#  define _call_coarsen_cell 1
#  undef _call_coarsen_cell_recursive
#  define _call_coarsen_cell_recursive 1
#endif
#if _call_coarsen_cell
#  undef _call_decrement_neighbors
#  define _call_decrement_neighbors 1
#endif
#if _call_refine_cell
#  undef _call_increment_neighbors
#  define _call_increment_neighbors 1
#  undef _call_refine_cell
#  define _call_refine_cell 1
#endif
#if _call_multigrid_debug
#  undef _call_cartesian_debug
#  define _call_cartesian_debug 1
#endif
#if _call_refine_biquadratic
#  undef _call_biquadratic
#  define _call_biquadratic 1
#endif
#if _call_refine_bilinear
#  undef _call_bilinear
#  define _call_bilinear 1
#endif
#if _call_restriction_face
#  undef _call_face_average
#  define _call_face_average 1
#endif
#if _call_decrement_neighbors
#  undef _call_free_children
#  define _call_free_children 1
#endif
#if _call_increment_neighbors
#  undef _call_alloc_children
#  define _call_alloc_children 1
#endif
#line 1017 "/home/fpl/softwares/basilisk/src/common.h"





























int list_len (scalar * list)
{
  if (!list) return 0;
  int ns = 0;
  strongif (list) for (scalar s = *list, *_i0 = list; ((scalar *)&s)->i >= 0; s = *++_i0) ns++;
  return ns;
}

scalar * list_append (scalar * list, scalar s)
{
  int len = list_len (list);
  list = (scalar *) prealloc (list, (len + 2)*sizeof(scalar),__func__,__FILE__,__LINE__);
  list[len] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_prepend (scalar * list, scalar s)
{
  int len = list_len (list);
  list = (scalar *) prealloc (list, (len + 2)*sizeof(scalar),__func__,__FILE__,__LINE__);
  for (int i = len; i >= 1; i--)
    list[i] = list[i-1];
  list[0] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_add (scalar * list, scalar s)
{
  strongif (list) for (scalar t = *list, *_i1 = list; ((scalar *)&t)->i >= 0; t = *++_i1)
    if (t.i == s.i)
      return list;
  return list_append (list, s);
}

int list_lookup (scalar * l, scalar s)
{
  if (l != NULL)
    strongif (l) for (scalar s1 = *l, *_i2 = l; ((scalar *)&s1)->i >= 0; s1 = *++_i2)
      if (s1.i == s.i)
 return true;
  return false;
}

scalar * list_copy (scalar * l)
{
  scalar * list = NULL;
  if (l != NULL)
    strongif (l) for (scalar s = *l, *_i3 = l; ((scalar *)&s)->i >= 0; s = *++_i3)
      list = list_append (list, s);
  return list;
}

scalar * list_concat (scalar * l1, scalar * l2)
{
  scalar * l3 = list_copy (l1);
  strongif (l2) for (scalar s = *l2, *_i4 = l2; ((scalar *)&s)->i >= 0; s = *++_i4)
    l3 = list_append (l3, s);
  return l3;
}

void list_print (scalar * l, FILE * fp)
{
  int i = 0;
  strongif (l) for (scalar s = *l, *_i5 = l; ((scalar *)&s)->i >= 0; s = *++_i5)
    fprintf (fp, "%s%s", i++ == 0 ? "{" : ",", _attribute[s.i].name);
  fputs (i > 0 ? "}\n" : "{}\n", fp);
}

int vectors_len (vector * list)
{
  if (!list) return 0;
  int nv = 0;
  strongif (list) for (vector v = *list, *_i6 = list; ((scalar *)&v)->i >= 0; v = *++_i6) nv++;
  return nv;
}

vector * vectors_append (vector * list, vector v)
{
  int len = vectors_len (list);
  list = (vector *) prealloc (list, (len + 2)*sizeof(vector),__func__,__FILE__,__LINE__);
  list[len] = v;
  list[len + 1] = (vector){{-1}};
  return list;
}

vector * vectors_add (vector * list, vector v)
{
  strongif (list) for (vector w = *list, *_i7 = list; ((scalar *)&w)->i >= 0; w = *++_i7) {
    bool id = true;
    {
#line 1137

      if (w.x.i != v.x.i)
 id = false;
#line 1137

      if (w.y.i != v.y.i)
 id = false;}
    if (id)
      return list;
  }
  return vectors_append (list, v);
}

vector * vectors_copy (vector * l)
{
  vector * list = NULL;
  if (l != NULL)
    strongif (l) for (vector v = *l, *_i8 = l; ((scalar *)&v)->i >= 0; v = *++_i8)
      list = vectors_append (list, v);
  return list;
}

vector * vectors_from_scalars (scalar * s)
{
  vector * list = NULL;
  while (s->i >= 0) {
    vector v;
    {
#line 1160
 {
      if (!(s->i >= 0)) qassert ("/home/fpl/softwares/basilisk/src/common.h", 1161, "s->i >= 0");
      v.x = *s++;
    }
#line 1160
 {
      if (!(s->i >= 0)) qassert ("/home/fpl/softwares/basilisk/src/common.h", 1161, "s->i >= 0");
      v.y = *s++;
    }}
    list = vectors_append (list, v);
  }
  return list;
}

int tensors_len (tensor * list)
{
  if (!list) return 0;
  int nt = 0;
  strongif (list) for (tensor t = *list, *_i9 = list; ((scalar *)&t)->i >= 0; t = *++_i9) nt++;
  return nt;
}

tensor * tensors_append (tensor * list, tensor t)
{
  int len = tensors_len (list);
  list = (tensor *) prealloc (list, (len + 2)*sizeof(tensor),__func__,__FILE__,__LINE__);
  list[len] = t;
  list[len + 1] = (tensor){{{-1}}};
  return list;
}

tensor * tensors_from_vectors (vector * v)
{
  tensor * list = NULL;
  while (v->x.i >= 0) {
    tensor t;
    {
#line 1191
 {
      if (!(v->x.i >= 0)) qassert ("/home/fpl/softwares/basilisk/src/common.h", 1192, "v->x.i >= 0");
      t.x = *v++;
    }
#line 1191
 {
      if (!(v->y.i >= 0)) qassert ("/home/fpl/softwares/basilisk/src/common.h", 1192, "v->x.i >= 0");
      t.y = *v++;
    }}
    list = tensors_append (list, t);
  }
  return list;
}

scalar * all = NULL;
scalar * baseblock = NULL;



scalar (* init_scalar) (scalar, const char *);
scalar (* init_vertex_scalar) (scalar, const char *);
vector (* init_vector) (vector, const char *);
tensor (* init_tensor) (tensor, const char *);
vector (* init_face_vector) (vector, const char *);





typedef struct _Event Event;
typedef int (* Expr) (int *, double *, Event *);

struct _Event {
  int last, nexpr;
  int (* action) (const int, const double, Event *);
  Expr expr[3];
  int * arrayi;
  double * arrayt;
  char * file;
  int line;
  char * name;
  double t;
  int i, a;
  void * data;
  Event * next;
};

static Event * Events = NULL;

int iter = 0, inext = 0;
double t = 0, tnext = 0;
void init_events (void);
void event_register (Event event);
void _init_solver (void);



#if 1
static double mpi_time = 0.;
#endif

typedef struct {
  clock_t c;
  struct timeval tv;
  double tm;
} timer;

timer timer_start (void)
{
  timer t;
  t.c = clock();
  gettimeofday (&t.tv, NULL);
#if 1
  t.tm = mpi_time;
#endif
  return t;
}

double timer_elapsed (timer t)
{
  struct timeval tvend;
  gettimeofday (&tvend, NULL);
  return ((tvend.tv_sec - t.tv.tv_sec) +
   (tvend.tv_usec - t.tv.tv_usec)/1e6);
}



vector zerof= {{_NVARMAX + 0},{_NVARMAX + 1}};
vector unityf= {{_NVARMAX + 2},{_NVARMAX + 3}};
scalar unity= {_NVARMAX + 4};
scalar zeroc= {_NVARMAX + 5};



 vector fm = {{_NVARMAX + 2},{_NVARMAX + 3}};
 scalar cm = {(_NVARMAX + 4)};
#line 1296 "/home/fpl/softwares/basilisk/src/common.h"
static FILE ** qpopen_pipes = NULL;

FILE * qpopen (const char * command, const char * type)
{
  if (pid() > 0)
    return fopen ("/dev/null", type);
  FILE * fp = popen (command, type);
  if (fp) {
    FILE ** i = qpopen_pipes;
    int n = 0;
    while (i && *i) { n++; i++; }
    qpopen_pipes = (FILE * *) prealloc (qpopen_pipes, (n + 2)*sizeof(FILE *),__func__,__FILE__,__LINE__);
    qpopen_pipes[n] = fp;
    qpopen_pipes[n+1] = NULL;
  }
  return fp;
}

int qpclose (FILE * fp)
{
  if (pid() > 0)
    return fclose (fp);
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i == fp)
      *i = (FILE *) 1;
    i++;
  }
  return pclose (fp);
}

static void qpclose_all ()
{
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i != (FILE *) 1)
      pclose (*i);
    i++;
  }
  pfree (qpopen_pipes,__func__,__FILE__,__LINE__);
  qpopen_pipes = NULL;
}






FILE * lfopen (const char * name, const char * mode)
{
  char fname[80];
  sprintf (fname, "%s-%d", name, pid());
  return fopen (fname, mode);
}



void * matrix_new (int n, int p, size_t size)
{
  void ** m = ((void * *) pmalloc ((n)*sizeof(void *),__func__,__FILE__,__LINE__));
  char * a = ((char *) pmalloc ((n*p*size)*sizeof(char),__func__,__FILE__,__LINE__));
  for (int i = 0; i < n; i++)
    m[i] = a + i*p*size;
  return m;
}

double matrix_inverse (double ** m, int n, double pivmin)
{
  int indxc[n], indxr[n], ipiv[n];
  int i, icol = 0, irow = 0, j, k, l, ll;
  double big, dum, pivinv, minpiv = HUGE;

  for (j = 0; j < n; j++)
    ipiv[j] = -1;

  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++)
      if (ipiv[j] != 0)
 for (k = 0; k < n; k++) {
   if (ipiv[k] == -1) {
     if (fabs (m[j][k]) >= big) {
       big = fabs (m[j][k]);
       irow = j;
       icol = k;
     }
   }
 }
    ipiv[icol]++;
    if (irow != icol)
      for (l = 0; l < n; l++)
 swap (double, m[irow][l], m[icol][l]);
    indxr[i] = irow;
    indxc[i] = icol;
    if (fabs (m[icol][icol]) <= pivmin)
      return 0.;
    if (fabs (m[icol][icol]) < minpiv)
      minpiv = fabs (m[icol][icol]);
    pivinv = 1.0/m[icol][icol];
    m[icol][icol] = 1.0;
    for (l = 0; l < n; l++) m[icol][l] *= pivinv;
    for (ll = 0; ll < n; ll++)
      if (ll != icol) {
 dum = m[ll][icol];
 m[ll][icol] = 0.0;
 for (l = 0; l < n; l++)
   m[ll][l] -= m[icol][l]*dum;
      }
  }
  for (l = n - 1; l >= 0; l--) {
    if (indxr[l] != indxc[l])
      for (k = 0; k < n; k++)
 swap (double, m[k][indxr[l]], m[k][indxc[l]]);
  }
  return minpiv;
}

void matrix_free (void * m)
{
  pfree (((void **) m)[0],__func__,__FILE__,__LINE__);
  pfree (m,__func__,__FILE__,__LINE__);
}



void init_solver ()
{
#if _CADNA
  cadna_init (-1);
#endif
#if 1
  mpi_init();
#elif MTRACE == 1
  char * etrace = getenv ("MTRACE");
  pmtrace.fp = fopen (etrace ? etrace : "mtrace", "w");
  pmtrace.fname = systrdup (etrace ? etrace : "mtrace");
#endif
}

void allocate_globals (int nvar)
{
  _attribute = (_Attributes *) pcalloc (nvar + 1, sizeof (_Attributes),__func__,__FILE__,__LINE__);
  _attribute[0].write = (double *) pcalloc (1, sizeof(double),__func__,__FILE__,__LINE__);
  _attribute++;
  all = (scalar *) pmalloc (sizeof (scalar)*(nvar + 1),__func__,__FILE__,__LINE__);
  baseblock = (scalar *) pmalloc (sizeof (scalar)*(nvar + 1),__func__,__FILE__,__LINE__);
  for (int i = 0; i < nvar; i++)
    baseblock[i].i = all[i].i = i;
  baseblock[nvar].i = all[nvar].i = -1;
}

typedef void (* free_solver_func) (void);

static Array * free_solver_funcs = NULL;

void free_solver_func_add (free_solver_func func)
{
  if (!free_solver_funcs)
    free_solver_funcs = array_new();
  array_append (free_solver_funcs, &func, sizeof(free_solver_func));
}



static char * display_defaults = NULL;

struct _display {
  const char * commands;
  bool overwrite;
};

static void free_display_defaults () {
  pfree (display_defaults,__func__,__FILE__,__LINE__);
}

void display (struct _display p)
{
  if (display_defaults == NULL)
    free_solver_func_add (free_display_defaults);
  if (p.overwrite) {
    pfree (display_defaults,__func__,__FILE__,__LINE__);
    display_defaults = pmalloc (strlen(p.commands) + 2,__func__,__FILE__,__LINE__);
    strcpy (display_defaults, "@");
    strcat (display_defaults, p.commands);
  }
  else {
    if (!display_defaults)
      display_defaults = pstrdup ("@",__func__,__FILE__,__LINE__);
    display_defaults =
      prealloc (display_defaults,
        strlen(display_defaults) + strlen(p.commands) + 1,__func__,__FILE__,__LINE__);
    strcat (display_defaults, p.commands);
  }
}
#line 14 "60d_plate_advancing_simulation-cpp.c"
#line 1 "grid/quadtree.h"
#line 1 "/home/fpl/softwares/basilisk/src/grid/quadtree.h"


#line 1 "grid/tree.h"
#line 1 "/home/fpl/softwares/basilisk/src/grid/tree.h"
#line 1 "grid/mempool.h"
#line 1 "/home/fpl/softwares/basilisk/src/grid/mempool.h"





typedef struct _Pool Pool;

struct _Pool {
  Pool * next;
};

typedef struct {
  char * first, * lastb;
  size_t size;
  size_t poolsize;
  Pool * pool, * last;
} Mempool;

typedef struct {
  char * next;
} FreeBlock;

Mempool * mempool_new (size_t poolsize, size_t size)
{

  if (!(poolsize % 8 == 0)) qassert ("/home/fpl/softwares/basilisk/src/grid/mempool.h", 26, "poolsize % 8 == 0");
  if (!(size >= sizeof(FreeBlock))) qassert ("/home/fpl/softwares/basilisk/src/grid/mempool.h", 27, "size >= sizeof(FreeBlock)");


  poolsize = min(1 << 20, poolsize + sizeof(Pool));
  Mempool * m = ((Mempool *) pcalloc (1, sizeof(Mempool),__func__,__FILE__,__LINE__));
  m->poolsize = poolsize;
  m->size = size;
  return m;
}

void mempool_destroy (Mempool * m)
{
  Pool * p = m->pool;
  while (p) {
    Pool * next = p->next;
    pfree (p,__func__,__FILE__,__LINE__);
    p = next;
  }
  pfree (m,__func__,__FILE__,__LINE__);
}

void * mempool_alloc (Mempool * m)
{
  if (!m->first) {

    Pool * p = (Pool *) pmalloc (m->poolsize,__func__,__FILE__,__LINE__);
    p->next = NULL;
    if (m->last)
      m->last->next = p;
    else
      m->pool = p;
    m->last = p;
    m->first = m->lastb = ((char *)m->last) + sizeof(Pool);
    FreeBlock * b = (FreeBlock *) m->first;
    b->next = NULL;
  }
  void * ret = m->first;
  FreeBlock * b = (FreeBlock *) ret;
  char * next = b->next;
  if (!next) {
    m->lastb += m->size;
    next = m->lastb;
    if (next + m->size > ((char *) m->last) + m->poolsize)
      next = NULL;
    else {
      FreeBlock * b = (FreeBlock *) next;
      b->next = NULL;
    }
  }
  m->first = next;
#if TRASH
  double * v = (double *) ret;
  for (int i = 0; i < m->size/sizeof(double); i++)
    v[i] = undefined;
#endif
  return ret;
}

void * mempool_alloc0 (Mempool * m)
{
  void * ret = mempool_alloc (m);
  memset (ret, 0, m->size);
  return ret;
}

void mempool_free (Mempool * m, void * p)
{
#if TRASH
  double * v = (double *) p;
  for (int i = 0; i < m->size/sizeof(double); i++)
    v[i] = undefined;
#endif
  FreeBlock * b = (FreeBlock *) p;
  b->next = m->first;
  m->first = (char *) p;
}
#line 2 "/home/fpl/softwares/basilisk/src/grid/tree.h"




#line 1 "grid/memindex/range.h"
#line 1 "/home/fpl/softwares/basilisk/src/grid/memindex/range.h"
#line 15 "/home/fpl/softwares/basilisk/src/grid/memindex/range.h"
typedef struct {
  void ** p;
  int size;
} Memalloc;

typedef struct {
  int start, end;
} Memrange;
#line 34 "/home/fpl/softwares/basilisk/src/grid/memindex/range.h"
void memrange_alloc (Memrange * r, Memalloc * mem, int i)
{
  if (r->start == r->end) {
    r->start = i;
    r->end = i + 1;
    for (Memalloc * m = mem; m->p; m++) {
      *m->p = pcalloc (1, m->size,__func__,__FILE__,__LINE__);
      *m->p = (char *)(*m->p) - i*m->size;
    }
  }
  else if (i >= r->end) {
    for (Memalloc * m = mem; m->p; m++) {
      *m->p = prealloc ((char *)(*m->p) + r->start*m->size,
         m->size*(i + 1 - r->start),__func__,__FILE__,__LINE__);
      *m->p = (char *)(*m->p) - r->start*m->size;
      memset ((char *)(*m->p) + r->end*m->size, 0, (i - r->end + 1)*m->size);
    }
    r->end = i + 1;
  }
  else if (i < r->start) {
    for (Memalloc * m = mem; m->p; m++) {
      *m->p = prealloc ((char *)(*m->p) + r->start*m->size, m->size*(r->end - i),__func__,__FILE__,__LINE__);
      memmove ((char *)(*m->p) + (r->start - i)*m->size, *m->p,
        m->size*(r->end - r->start));
      memset ((char *)(*m->p), 0, (r->start - i)*m->size);
      *m->p = (char *)(*m->p) - i*m->size;
    }
    r->start = i;
  }
}
#line 73 "/home/fpl/softwares/basilisk/src/grid/memindex/range.h"
bool memrange_free (Memrange * r, Memalloc * mem, int i)
{
  if (i == r->start) {
    if (i == r->end - 1) {
      for (Memalloc * m = mem; m->p; m++) {
 pfree ((char *)(*m->p) + r->start*m->size,__func__,__FILE__,__LINE__);
 *m->p = NULL;
      }
      r->start = r->end = 0;
      return true;
    }
    else {
      for (i = i + 1; i < r->end &&
      !*(void **)((char *)(*mem->p) + i*mem->size); i++);
      for (Memalloc * m = mem; m->p; m++) {
 memmove ((char *)(*m->p) + r->start*m->size,
   (char *)(*m->p) + i*m->size, m->size*(r->end - i));
 *m->p = prealloc ((char *)(*m->p) + r->start*m->size,
    m->size*(r->end - i),__func__,__FILE__,__LINE__);
 *m->p = (char *)(*m->p) - i*m->size;
      }
      r->start = i;
    }
  }
  else if (i == r->end - 1) {
    for (i = i - 1; i >= r->start &&
    !*(void **)((char *)(*mem->p) + i*mem->size); i--);
    r->end = i + 1;
    for (Memalloc * m = mem; m->p; m++) {
      *m->p = prealloc ((char *)(*m->p) + r->start*m->size,
         m->size*(r->end - r->start),__func__,__FILE__,__LINE__);
      *m->p = (char *)(*m->p) - r->start*m->size;
    }
  }
  else {
    if (!(i > r->start && i < r->end)) qassert ("/home/fpl/softwares/basilisk/src/grid/memindex/range.h", 108, "i > r->start && i < r->end");
    for (Memalloc * m = mem; m->p; m++)
      memset ((char *)(*m->p) + i*m->size, 0, m->size);
  }
  return false;
}







struct _Memindex {
  Memrange r1;

  Memrange * r2;







  char *** b;



};
#line 171 "/home/fpl/softwares/basilisk/src/grid/memindex/range.h"
struct _Memindex * mem_new (int len)
{
  struct _Memindex * m = pcalloc (1, sizeof (struct _Memindex),__func__,__FILE__,__LINE__);
  return m;
}





void mem_destroy (struct _Memindex * m, int len)
{

  for (int i = m->r1.start; i < m->r1.end; i++)
    if (m->b[i]) {






      pfree (m->b[i] + m->r2[i].start,__func__,__FILE__,__LINE__);
    }
  if (m->b) {
    pfree (m->r2 + m->r1.start,__func__,__FILE__,__LINE__);



  }

  if (m->b)
    pfree (m->b + m->r1.start,__func__,__FILE__,__LINE__);
  pfree (m,__func__,__FILE__,__LINE__);
}
#line 218 "/home/fpl/softwares/basilisk/src/grid/memindex/range.h"
void mem_assign (struct _Memindex * m, int i, int j, int len, void * b)
{
  Memalloc mem[] = {{(void **)&m->b, sizeof(char **)},
      {(void **)&m->r2, sizeof(Memrange)},
      {NULL}};
  memrange_alloc (&m->r1, mem, i);
  Memalloc mem1[] = {{(void **)&m->b[i], sizeof(char *)},
       {NULL}};
  memrange_alloc (&m->r2[i], mem1, j);
  ((m)->b[i][j]) = b;
}
#line 259 "/home/fpl/softwares/basilisk/src/grid/memindex/range.h"
void mem_free (struct _Memindex * m, int i, int j, int len)
{
  Memalloc mem[] = {{(void **)&m->b[i], sizeof(char *)},
      {NULL}};
  if (memrange_free (&m->r2[i], mem, j)) {
    Memalloc mem[] = {{(void **)&m->b, sizeof(char **)},
        {(void **)&m->r2, sizeof(Memrange)},
        {NULL}};
    memrange_free (&m->r1, mem, i);
  }
}
#line 305 "/home/fpl/softwares/basilisk/src/grid/memindex/range.h"
#define foreach_mem(_m, _len, _i) {\
  Point point = {0};\
  for (point.i = max(Period.x*2, (_m)->r1.start);\
       point.i < min(_len - Period.x*2, (_m)->r1.end);\
       point.i += _i)\
    if ((_m)->b[point.i])\
      for (point.j = max(Period.y*2, (_m)->r2[point.i].start);\
    point.j < min(_len - Period.y*2, (_m)->r2[point.i].end);\
    point.j += _i)\
 if ((_m)->b[point.i][point.j]) {\

#line 315

#define end_foreach_mem() }}
#line 7 "/home/fpl/softwares/basilisk/src/grid/tree.h"
#line 24 "/home/fpl/softwares/basilisk/src/grid/tree.h"
typedef struct {
  unsigned short flags;

  unsigned short neighbors;
  int pid;
} Cell;

enum {
  active = 1 << 0,
  leaf = 1 << 1,
  border = 1 << 2,
  vertex = 1 << 3,
  user = 4,

  face_x = 1 << 0

  , face_y = 1 << 1




};

#define is_active(cell) ((cell).flags & active)
#define is_leaf(cell) ((cell).flags & leaf)
#define is_coarse() ((cell).neighbors > 0)
#define is_border(cell) ((cell).flags & border)
#define is_local(cell) ((cell).pid == pid())
#define is_vertex(cell) ((cell).flags & vertex)



typedef struct {
  int i;

  int j;




} IndexLevel;

typedef struct {
  IndexLevel * p;
  int n, nm;
} CacheLevel;

typedef struct {
  int i;

  int j;




  int level, flags;
} Index;

typedef struct {
  Index * p;
  int n, nm;
} Cache;



typedef struct {
  struct _Memindex * m;
  Mempool * pool;
  long nc;
  int len;
} Layer;

static size_t _size (size_t depth)
{
  return (1 << depth) + 2*2;
}

static size_t poolsize (size_t depth, size_t size)
{




  return sq(_size(depth))*size;



}

static Layer * new_layer (int depth)
{
  Layer * l = ((Layer *) pmalloc ((1)*sizeof(Layer),__func__,__FILE__,__LINE__));
  l->len = _size (depth);
  if (depth == 0)
    l->pool = NULL;
  else {
    size_t size = sizeof(Cell) + datasize;


    l->pool = mempool_new (poolsize (depth, size), (1 << 2)*size);
  }
  l->m = mem_new (l->len);
  l->nc = 0;
  return l;
}

static void destroy_layer (Layer * l)
{
  if (l->pool)
    mempool_destroy (l->pool);
  mem_destroy (l->m, l->len);
  pfree (l,__func__,__FILE__,__LINE__);
}



typedef struct {
  Grid g;
  Layer ** L;

  Cache leaves;
  Cache faces;
  Cache vertices;
  Cache refined;
  CacheLevel * active;
  CacheLevel * prolongation;
  CacheLevel * boundary;

  CacheLevel * restriction;

  bool dirty;
} Tree;



struct _Point {

  int i;

  int j;




  int level;
#ifdef foreach_block
  int l;
#define _BLOCK_INDEX , point.l
#else
#define _BLOCK_INDEX
#endif
};
static Point last_point;



static void cache_level_append (CacheLevel * c, Point p)
{
  if (c->n >= c->nm) {
    c->nm += 128;
    c->p = (IndexLevel *) prealloc (c->p, (c->nm)*sizeof(IndexLevel),__func__,__FILE__,__LINE__);
  }
  c->p[c->n].i = p.i;

  c->p[c->n].j = p.j;




  c->n++;
}

static void cache_level_shrink (CacheLevel * c)
{
  if (c->nm > (c->n/128 + 1)*128) {
    c->nm = (c->n/128 + 1)*128;
    if (!(c->nm > c->n)) qassert ("/home/fpl/softwares/basilisk/src/grid/tree.h", 200, "c->nm > c->n");
    c->p = (IndexLevel *) prealloc (c->p, sizeof (Index)*c->nm,__func__,__FILE__,__LINE__);
  }
}

static void cache_append (Cache * c, Point p, unsigned short flags)
{
  if (c->n >= c->nm) {
    c->nm += 128;
    c->p = (Index *) prealloc (c->p, (c->nm)*sizeof(Index),__func__,__FILE__,__LINE__);
  }
  c->p[c->n].i = p.i;

  c->p[c->n].j = p.j;




  c->p[c->n].level = p.level;
  c->p[c->n].flags = flags;
  c->n++;
}

void cache_shrink (Cache * c)
{
  cache_level_shrink ((CacheLevel *)c);
}
#line 243 "/home/fpl/softwares/basilisk/src/grid/tree.h"
#define allocated(k,l,n) (((point.i+k) >= (((Tree *)grid)->L[point.level]->m)->r1.start && (point.i+k) < (((Tree *)grid)->L[point.level]->m->r1.end) && (((Tree *)grid)->L[point.level]->m)->b[point.i+k] && (point.j+l) >= (((Tree *)grid)->L[point.level]->m)->r2[point.i+k].start && (point.j+l) < (((Tree *)grid)->L[point.level]->m)->r2[point.i+k].end && (((Tree *)grid)->L[point.level]->m)->b[point.i+k][point.j+l])\
                               )\

#line 245

#define NEIGHBOR(k,l,n) (((((Tree *)grid)->L[point.level]->m)->b[point.i+k][point.j+l])\
                            )\

#line 248

#define PARENT(k,l,n) (((((Tree *)grid)->L[point.level-1]->m)->b[(point.i+2)/2+k][(point.j+2)/2+l])\
                                                    )\

#line 251

#define allocated_child(k,l,n) (level < depth() &&\
         ((2*point.i-2 +k) >= (((Tree *)grid)->L[point.level+1]->m)->r1.start && (2*point.i-2 +k) < (((Tree *)grid)->L[point.level+1]->m->r1.end) && (((Tree *)grid)->L[point.level+1]->m)->b[2*point.i-2 +k] && (2*point.j-2 +l) >= (((Tree *)grid)->L[point.level+1]->m)->r2[2*point.i-2 +k].start && (2*point.j-2 +l) < (((Tree *)grid)->L[point.level+1]->m)->r2[2*point.i-2 +k].end && (((Tree *)grid)->L[point.level+1]->m)->b[2*point.i-2 +k][2*point.j-2 +l])\
\
                             )\

#line 256

#define CHILD(k,l,n) (((((Tree *)grid)->L[point.level+1]->m)->b[2*point.i-2 +k][2*point.j-2 +l])\
                                                )\

#line 259

#line 284 "/home/fpl/softwares/basilisk/src/grid/tree.h"
#define CELL(m) (*((Cell *)(m)))


#define depth() (grid->depth)
#define aparent(k,l,n) CELL(PARENT(k,l,n))
#define child(k,l,n) CELL(CHILD(k,l,n))


#define cell CELL(NEIGHBOR(0,0,0))
#define neighbor(k,l,n) CELL(NEIGHBOR(k,l,n))
#define neighborp(l,m,n) (Point) {\
    point.i + l,\
\
    point.j + m,\
\
\
\
\
    point.level\
    _BLOCK_INDEX\
}\

#line 305



#define data(k,l,n) ((double *) (NEIGHBOR(k,l,n) + sizeof(Cell)))
#define fine(a,k,p,n) ((double *) (CHILD(k,p,n) + sizeof(Cell)))[_index(a,n)]
#define coarse(a,k,p,n) ((double *) (PARENT(k,p,n) + sizeof(Cell)))[_index(a,n)]

#define POINT_VARIABLES\
  VARIABLES\
  int level = point.level; NOT_UNUSED(level);\
\
\
\
  struct { int x, y; } child = {\
    2*((point.i+2)%2)-1, 2*((point.j+2)%2)-1\
  };\
\
\
\
\
\
  NOT_UNUSED(child);\
  Point parent = point; NOT_UNUSED(parent);\
  parent.level--;\
  parent.i = (point.i + 2)/2;\
\
  parent.j = (point.j + 2)/2;\
\
\
\
\

#line 336


#line 1 "grid/foreach_cell.h"
#line 1 "/home/fpl/softwares/basilisk/src/grid/foreach_cell.h"
#line 66 "/home/fpl/softwares/basilisk/src/grid/foreach_cell.h"
#define foreach_cell_root(root)\
  {\
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);\
    Point point = {0};\
\
\
\
    struct { int l, i, j, stage; } stack[20];\
\
\
\
\
    int _s = -1;\
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].stage = 0; };\
    while (_s >= 0) {\
      int stage;\
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; stage = stack[_s].stage; _s--; };\
      if (!allocated (0,0,0))\
 continue;\
      switch (stage) {\
      case 0: {\
 POINT_VARIABLES;\
\

#line 89

#define end_foreach_cell_root()\
        if (point.level < grid->depth) {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 1; };\
          { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; };\
        }\
        break;\
      }\
\
\
\
      case 1: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 2; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; }; break;\
      case 2: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 3; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; }; break;\
      case 3: { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; }; break;\
\
      }\
    }\
  }\

#line 123


#define foreach_cell() {\
\
\
\
  Point root = {2,2,0};\
\
\
\
  foreach_cell_root (root)\

#line 134

#define end_foreach_cell() end_foreach_cell_root() }

#define foreach_cell_all() {\
  Point root = {0};\
  for (root.i = 2*Period.x; root.i <= 2*(2 - Period.x); root.i++)\
\
    for (root.j = 2*Period.y; root.j <= 2*(2 - Period.y); root.j++)\
\
\
\
\
 foreach_cell_root (root)\

#line 147

#define end_foreach_cell_all() end_foreach_cell_root() }

#define foreach_cell_post_root(condition, root)\
  {\
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);\
    Point point = {0};\
\
\
\
    struct { int l, i, j, stage; } stack[20];\
\
\
\
\
    int _s = -1;\
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].stage = 0; };\
    while (_s >= 0) {\
      int stage;\
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; stage = stack[_s].stage; _s--; };\
      if (!allocated (0,0,0))\
 continue;\
      switch (stage) {\
      case 0: {\
        POINT_VARIABLES;\
 if (point.level == grid->depth) {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 8; };\
 }\
 else {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 1; };\
   if (condition)\
     { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; };\
 }\
 break;\
      }\
\
\
\
\
\
\
\
      case 1:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 2; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; };\
 break;\
      case 2:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 3; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; };\
 break;\
      case 3:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 4; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; };\
 break;\
\
      default: {\
        POINT_VARIABLES;\
\

#line 244

#define end_foreach_cell_post_root()\
      }\
      }\
    }\
  }\

#line 250


#define foreach_cell_post(condition)\
  {\
\
\
\
    Point root = {2,2,0};\
\
\
\
    foreach_cell_post_root(condition, root)\

#line 262

#define end_foreach_cell_post() end_foreach_cell_post_root() }

#define foreach_cell_post_all(condition) {\
  Point root = {0};\
  for (root.i = 0; root.i <= 2*2; root.i++)\
\
    for (root.j = 0; root.j <= 2*2; root.j++)\
\
\
\
\
 foreach_cell_post_root (condition, root)\

#line 275

#define end_foreach_cell_post_all() end_foreach_cell_post_root() }

#define foreach_leaf() foreach_cell()\
  if (is_leaf (cell)) {\
    if (is_active(cell) && is_local(cell)) {\

#line 281

#define end_foreach_leaf() } continue; } end_foreach_cell()
#line 339 "/home/fpl/softwares/basilisk/src/grid/tree.h"
#line 356 "/home/fpl/softwares/basilisk/src/grid/tree.h"
#define foreach_child() {\
  int _i = 2*point.i - 2, _j = 2*point.j - 2;\
  point.level++;\
  for (int _k = 0; _k < 2; _k++) {\
    point.i = _i + _k;\
    for (int _l = 0; _l < 2; _l++) {\
      point.j = _j + _l;\
      POINT_VARIABLES;\

#line 364

#define end_foreach_child()\
    }\
  }\
  point.i = (_i + 2)/2; point.j = (_j + 2)/2;\
  point.level--;\
}\

#line 371

#define foreach_child_break() _k = _l = 2
#line 402 "/home/fpl/softwares/basilisk/src/grid/tree.h"
#define is_refined_check() ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) &&\
    point.i > 0 && point.i < (1 << level) + 2*2 - 1\
\
    && point.j > 0 && point.j < (1 << level) + 2*2 - 1\
\
\
\
\
    )\

#line 411


#define foreach_cache(_cache) {\
  OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.i = 2;\
\
  point.j = 2;\
\
\
\
\
  int _k; unsigned short _flags; NOT_UNUSED(_flags);\
  OMP(omp for schedule(static))\
  for (_k = 0; _k < _cache.n; _k++) {\
    point.i = _cache.p[_k].i;\
\
    point.j = _cache.p[_k].j;\
\
\
\
\
    point.level = _cache.p[_k].level;\
    _flags = _cache.p[_k].flags;\
    POINT_VARIABLES;\

#line 437

#define end_foreach_cache() } } }

#define foreach_cache_level(_cache,_l) {\
  OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.i = 2;\
\
  point.j = 2;\
\
\
\
\
  point.level = _l;\
  int _k;\
  OMP(omp for schedule(static))\
  for (_k = 0; _k < _cache.n; _k++) {\
    point.i = _cache.p[_k].i;\
\
    point.j = _cache.p[_k].j;\
\
\
\
\
    POINT_VARIABLES;\

#line 463

#define end_foreach_cache_level() } } }

#define foreach_boundary_level(_l) {\
  if (_l <= depth()) {\
    { if (((Tree *)grid)->dirty) update_cache_f(); };\
    CacheLevel _boundary = ((Tree *)grid)->boundary[_l];\
    foreach_cache_level (_boundary,_l)\

#line 471

#define end_foreach_boundary_level() end_foreach_cache_level(); }}



#define foreach_boundary(_b) {\
  for (int _l = depth(); _l >= 0; _l--)\
    foreach_boundary_level(_l) {\
      if ((- cell.pid - 1) == _b)\
 for (int _d = 0; _d < 2; _d++) {\
   for (int _i = -1; _i <= 1; _i += 2) {\
     if (_d == 0) ig = _i; else if (_d == 1) jg = _i; else kg = _i;\
     if (allocated(-ig,-jg,-kg) &&\
  is_leaf (neighbor(-ig,-jg,-kg)) &&\
  !(neighbor(-ig,-jg,-kg).pid < 0) &&\
  is_local(neighbor(-ig,-jg,-kg))) {\
       point.i -= ig; x -= ig*Delta/2.;\
\
       point.j -= jg; y -= jg*Delta/2.;\
\
\
\
\

#line 494

#define end_foreach_boundary()\
       point.i += ig; x += ig*Delta/2.;\
\
       point.j += jg; y += jg*Delta/2.;\
\
\
\
\
            }\
   }\
   ig = jg = kg = 0;\
 }\
    } end_foreach_boundary_level(); }\

#line 508


#define foreach_halo(_name,_l) {\
  if (_l <= depth()) {\
    { if (((Tree *)grid)->dirty) update_cache_f(); };\
    CacheLevel _cache = ((Tree *)grid)->_name[_l];\
    foreach_cache_level (_cache, _l)\

#line 515

#define end_foreach_halo() end_foreach_cache_level(); }}

#line 1 "grid/neighbors.h"
#line 1 "/home/fpl/softwares/basilisk/src/grid/neighbors.h"
#line 17 "/home/fpl/softwares/basilisk/src/grid/neighbors.h"
#define foreach_neighbor(_s) {\
  int _nn = _s + 0 ? _s + 0 : 2;\
  int _i = point.i, _j = point.j;\
  for (int _k = - _nn; _k <= _nn; _k++) {\
    point.i = _i + _k;\
    for (int _l = - _nn; _l <= _nn; _l++) {\
      point.j = _j + _l;\
      POINT_VARIABLES;\

#line 25

#define end_foreach_neighbor()\
    }\
  }\
  point.i = _i; point.j = _j;\
}\

#line 31

#define foreach_neighbor_break() _k = _l = _nn + 1
#line 519 "/home/fpl/softwares/basilisk/src/grid/tree.h"

static inline bool has_local_children (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 521 "/home/fpl/softwares/basilisk/src/grid/tree.h"

   { foreach_child()
    if (is_local(cell))
      return true; end_foreach_child(); }
  return false;

#if _call_has_local_children
}
#define _IN_STENCIL 1

#line 520
static bool _has_local_children (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 521 "/home/fpl/softwares/basilisk/src/grid/tree.h"

   { foreach_child()
    IF (is_local(cell))
      return true; end_foreach_child(); }
  return false;

#undef _IN_STENCIL

#endif

#line 526
}

static inline void cache_append_face (Point point, unsigned short flags)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 529 "/home/fpl/softwares/basilisk/src/grid/tree.h"

  Tree * q = ((Tree *)grid);
  cache_append (&q->faces, point, flags);

  if (!is_vertex(cell)) {
    cache_append (&q->vertices, point, 0);
    cell.flags |= vertex;
  }
  {
#line 537

    if ((flags & face_y) && !is_vertex(neighbor(1,0,0))) {
      cache_append (&q->vertices, neighborp(1,0,0), 0);
      neighbor(1,0,0).flags |= vertex;
    }
#line 537

    if ((flags & face_x) && !is_vertex(neighbor(0,1,0))) {
      cache_append (&q->vertices, neighborp(0,1,0), 0);
      neighbor(0,1,0).flags |= vertex;
    }}
#line 552 "/home/fpl/softwares/basilisk/src/grid/tree.h"

#if _call_cache_append_face
}
#define _IN_STENCIL 1

#line 528
static void _cache_append_face (Point point, unsigned short flags)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 529 "/home/fpl/softwares/basilisk/src/grid/tree.h"

  Tree * q = ((Tree *)grid);
  cache_append (&q->faces, point, flags);

  IF (!is_vertex(cell)) {
    cache_append (&q->vertices, point, 0);
    cell.flags |= vertex;
  }
  {
#line 537

    IF ((flags & face_y) && !is_vertex(neighbor(1,0,0))) {
      cache_append (&q->vertices, neighborp(1,0,0), 0);
      neighbor(1,0,0).flags |= vertex;
    }
#line 537

    IF ((flags & face_x) && !is_vertex(neighbor(0,1,0))) {
      cache_append (&q->vertices, neighborp(0,1,0), 0);
      neighbor(0,1,0).flags |= vertex;
    }}
#line 552 "/home/fpl/softwares/basilisk/src/grid/tree.h"

#undef _IN_STENCIL

#endif

#line 552
}



static void update_cache_f (void)
{
  Tree * q = ((Tree *)grid);

   { foreach_cache (q->vertices){

#line 560 "/home/fpl/softwares/basilisk/src/grid/tree.h"

    if (level <= depth() && allocated(0,0,0))
      cell.flags &= ~vertex; } end_foreach_cache(); }


  q->leaves.n = q->faces.n = q->vertices.n = 0;
  for (int l = 0; l <= depth(); l++)
    q->active[l].n = q->prolongation[l].n =
      q->boundary[l].n = q->restriction[l].n = 0;

  const unsigned short fboundary = 1 << user;
   { foreach_cell(){

#line 571 "/home/fpl/softwares/basilisk/src/grid/tree.h"
 {



    if (is_local(cell) && is_active(cell)) {


      cache_level_append (&q->active[level], point);
    }
#line 595 "/home/fpl/softwares/basilisk/src/grid/tree.h"
    if (!(cell.pid < 0)) {

       { foreach_neighbor (2)
 if (allocated(0,0,0) && (cell.pid < 0) && !(cell.flags & fboundary)) {
   cache_level_append (&q->boundary[level], point);
   cell.flags |= fboundary;
 } end_foreach_neighbor(); }
    }

    else if (level > 0 && is_local(aparent(0,0,0)))
      cache_level_append (&q->restriction[level], point);

    if (is_leaf (cell)) {
      if (is_local(cell)) {
 cache_append (&q->leaves, point, 0);

 unsigned short flags = 0;
 {
#line 612

   if ((neighbor(-1,0,0).pid < 0) || (!is_leaf(neighbor(-1,0,0)) && !neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) ||
       is_leaf(neighbor(-1,0,0)))
     flags |= face_x;
#line 612

   if ((neighbor(0,-1,0).pid < 0) || (!is_leaf(neighbor(0,-1,0)) && !neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) ||
       is_leaf(neighbor(0,-1,0)))
     flags |= face_y;}
 if (flags)
   cache_append (&q->faces, point, flags);
 {
#line 618

   if ((neighbor(1,0,0).pid < 0) || (!is_leaf(neighbor(1,0,0)) && !neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0) ||
       (!is_local(neighbor(1,0,0)) && is_leaf(neighbor(1,0,0))))
     cache_append (&q->faces, neighborp(1,0,0), face_x);
#line 618

   if ((neighbor(0,1,0).pid < 0) || (!is_leaf(neighbor(0,1,0)) && !neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0) ||
       (!is_local(neighbor(0,1,0)) && is_leaf(neighbor(0,1,0))))
     cache_append (&q->faces, neighborp(0,1,0), face_y);}

 for (int i = 0; i <= 1; i++)

   for (int j = 0; j <= 1; j++)




       if (!is_vertex(neighbor(i,j,k))) {
  cache_append (&q->vertices, neighborp(i,j,k), 0);
  neighbor(i,j,k).flags |= vertex;
       }

        if (cell.neighbors > 0)
   cache_level_append (&q->prolongation[level], point);
      }
      else if (!(cell.pid < 0) || is_local(aparent(0,0,0))) {

 unsigned short flags = 0;
 {
#line 641

   if (allocated(-1,0,0) &&
       is_local(neighbor(-1,0,0)) && (!is_leaf(neighbor(-1,0,0)) && !neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0))
     flags |= face_x;
#line 641

   if (allocated(0,-1,0) &&
       is_local(neighbor(0,-1,0)) && (!is_leaf(neighbor(0,-1,0)) && !neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0))
     flags |= face_y;}
 if (flags)
   cache_append_face (point, flags);
 {
#line 647

   if (allocated(1,0,0) && is_local(neighbor(1,0,0)) &&
       (!is_leaf(neighbor(1,0,0)) && !neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0))
     cache_append_face (neighborp(1,0,0), face_x);
#line 647

   if (allocated(0,1,0) && is_local(neighbor(0,1,0)) &&
       (!is_leaf(neighbor(0,1,0)) && !neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0))
     cache_append_face (neighborp(0,1,0), face_y);}
      }

      continue;

    }
  } } end_foreach_cell(); }


  cache_shrink (&q->leaves);
  cache_shrink (&q->faces);
  cache_shrink (&q->vertices);
  for (int l = 0; l <= depth(); l++) {
    cache_level_shrink (&q->active[l]);
    cache_level_shrink (&q->prolongation[l]);
    cache_level_shrink (&q->boundary[l]);
    cache_level_shrink (&q->restriction[l]);
}

  q->dirty = false;


  for (int l = depth(); l >= 0; l--)
     { foreach_boundary_level (l){

#line 673 "/home/fpl/softwares/basilisk/src/grid/tree.h"

      cell.flags &= ~fboundary; } end_foreach_boundary_level(); }



  grid->n = q->leaves.n;

#if !1
  grid->tn = grid->n;
  grid->maxdepth = grid->depth;
#endif
}

#define foreach() { if (((Tree *)grid)->dirty) update_cache_f(); }; foreach_cache(((Tree *)grid)->leaves)
#define end_foreach() end_foreach_cache()

#define foreach_face_generic()\
  { if (((Tree *)grid)->dirty) update_cache_f(); };\
  foreach_cache(((Tree *)grid)->faces) 
#line 690

#define end_foreach_face_generic() end_foreach_cache()

#define is_face_x() (_flags & face_x)

#define is_face_y() (_flags & face_y)





#define foreach_vertex()\
  { if (((Tree *)grid)->dirty) update_cache_f(); };\
  foreach_cache(((Tree *)grid)->vertices) {\
    x -= Delta/2.;\
\
    y -= Delta/2.;\
\
\
\
\

#line 712

#define end_foreach_vertex() } end_foreach_cache()
#line 724 "/home/fpl/softwares/basilisk/src/grid/tree.h"
#define foreach_level(l) {\
  if (l <= depth()) {\
    { if (((Tree *)grid)->dirty) update_cache_f(); };\
    CacheLevel _active = ((Tree *)grid)->active[l];\
    foreach_cache_level (_active,l)\

#line 729

#define end_foreach_level() end_foreach_cache_level(); }}

#define foreach_coarse_level(l) foreach_level(l) if (!is_leaf(cell)) {
#define end_foreach_coarse_level() } end_foreach_level()

#define foreach_level_or_leaf(l) {\
  for (int _l1 = l; _l1 >= 0; _l1--)\
    foreach_level(_l1)\
      if (_l1 == l || is_leaf (cell)) {\

#line 739

#define end_foreach_level_or_leaf() } end_foreach_level(); }

#if TRASH
# undef trash
# define trash(list) reset(list, undefined)
#endif

void reset (void * alist, double val)
{
  scalar * list = (scalar *) alist;
  Tree * q = ((Tree *)grid);

  for (int l = 0; l <= depth(); l++) {
    Layer * L = q->L[l];
     { foreach_mem (L->m, L->len, 1){

#line 754 "/home/fpl/softwares/basilisk/src/grid/tree.h"
 {
      point.level = l;
      strongif (list) for (scalar s = *list, *_i10 = list; ((scalar *)&s)->i >= 0; s = *++_i10) {
 if (!is_constant(s))
   for (int b = 0; b < _attribute[s.i].block; b++)
     data(0,0,0)[s.i + b] = val;
      }
    } } end_foreach_mem(); }
  }
}

#define cache_level_resize(name, a)\
{\
  for (int i = 0; i <= depth() - a; i++)\
    pfree (q->name[i].p,__func__,__FILE__,__LINE__);\
  pfree (q->name,__func__,__FILE__,__LINE__);\
  q->name = ((CacheLevel *) pcalloc (depth() + 1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));\
}\

#line 772


static void update_depth (int inc)
{
  Tree * q = ((Tree *)grid);
  grid->depth += inc;
  q->L = &(q->L[-1]);
  q->L = (Layer * *) prealloc (q->L, (grid->depth + 2)*sizeof(Layer *),__func__,__FILE__,__LINE__);
  q->L = &(q->L[1]);
  if (inc > 0)
    q->L[grid->depth] = new_layer (grid->depth);
  cache_level_resize (active, inc);
  cache_level_resize (prolongation, inc);
  cache_level_resize (boundary, inc);
  cache_level_resize (restriction, inc);
}
#line 814 "/home/fpl/softwares/basilisk/src/grid/tree.h"
typedef void (* PeriodicFunction) (struct _Memindex *, int, int, int, void *);

static void periodic_function (struct _Memindex * m, int i, int j, int len, void * b,
          PeriodicFunction f)
{
  f(m, i, j, len, b);
  if (Period.x) {
    int nl = len - 2*2;
    for (int l = - 1; l <= 1; l += 2)
      for (int n = i + l*nl; n >= 0 && n < len; n += l*nl)
 f(m, n, j, len, b);
    if (Period.y)
      for (int l = - 1; l <= 1; l += 2)
 for (int n = j + l*nl; n >= 0 && n < len; n += l*nl) {
   f(m, i, n, len, b);
   for (int o = - 1; o <= 1; o += 2)
     for (int p = i + o*nl; p >= 0 && p < len; p += o*nl)
       f(m, p, n, len, b);
 }
  }
  else if (Period.y) {
    int nl = len - 2*2;
    for (int l = - 1; l <= 1; l += 2)
      for (int n = j + l*nl; n >= 0 && n < len; n += l*nl)
 f(m, i, n, len, b);
  }
}

static void assign_periodic (struct _Memindex * m, int i, int j, int len, void * b)
{
  periodic_function (m, i, j, len, b, mem_assign);
}

static void free_periodic (struct _Memindex * m, int i, int j, int len)
{
  periodic_function (m, i, j, len, NULL, (PeriodicFunction) mem_free);
}
#line 929 "/home/fpl/softwares/basilisk/src/grid/tree.h"
static void alloc_children (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 930 "/home/fpl/softwares/basilisk/src/grid/tree.h"

  if (point.level == grid->depth)
    update_depth (+1);
  else if (allocated_child(0,0,0))
    return;


  Layer * L = ((Tree *)grid)->L[point.level + 1];
  L->nc++;
  size_t len = sizeof(Cell) + datasize;
  char * b = (char *) mempool_alloc0 (L->pool);
  int i = 2*point.i - 2;
  for (int k = 0; k < 2; k++, i++) {




    int j = 2*point.j - 2;
    for (int l = 0; l < 2; l++, j++) {
      assign_periodic (L->m, i, j, L->len, b);
      b += len;
    }
#line 962 "/home/fpl/softwares/basilisk/src/grid/tree.h"
  }

  int pid = cell.pid;
   { foreach_child() {
    cell.pid = pid;
#if TRASH
    strongif (all) for (scalar s = *all, *_i11 = all; ((scalar *)&s)->i >= 0; s = *++_i11)
      val(s,0,0,0) = undefined;
#endif
  } end_foreach_child(); }

#if _call_alloc_children
}
#define _IN_STENCIL 1

#line 929
static void _alloc_children (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 930 "/home/fpl/softwares/basilisk/src/grid/tree.h"

  IF (point.level == grid->depth)
    update_depth (+1);
  IF (allocated_child(0,0,0))
    return;


  Layer * L = ((Tree *)grid)->L[point.level + 1];
  L->nc++;
  size_t len = sizeof(Cell) + datasize;
  char * b = (char *) mempool_alloc0 (L->pool);
  int i = 2*point.i - 2;
  for (int k = 0; k < 2; k++, i++) {




    int j = 2*point.j - 2;
    for (int l = 0; l < 2; l++, j++) {
      assign_periodic (L->m, i, j, L->len, b);
      b += len;
    }
#line 962 "/home/fpl/softwares/basilisk/src/grid/tree.h"
  }

  int pid = cell.pid;
   { foreach_child() {
    cell.pid = pid;
#if TRASH
    strongif (all) for (scalar s = *all, *_i11 = all; ((scalar *)&s)->i >= 0; s = *++_i11)
      _stencil_val(__FILE__,__LINE__,s,0,0,0) = undefined;
#endif
  } end_foreach_child(); }

#undef _IN_STENCIL

#endif

#line 972
}
#line 991 "/home/fpl/softwares/basilisk/src/grid/tree.h"
static void free_children (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 992 "/home/fpl/softwares/basilisk/src/grid/tree.h"


  Layer * L = ((Tree *)grid)->L[point.level + 1];
  int i = 2*point.i - 2, j = 2*point.j - 2;
  if (!(((L->m)->b[i][j]))) qassert ("/home/fpl/softwares/basilisk/src/grid/tree.h", 996, "mem_data (L->m,i,j)");
  mempool_free (L->pool, ((L->m)->b[i][j]));
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++)
      free_periodic (L->m, i + k, j + l, L->len);
  if (--L->nc == 0) {
    destroy_layer (L);
    if (!(point.level + 1 == grid->depth)) qassert ("/home/fpl/softwares/basilisk/src/grid/tree.h", 1003, "point.level + 1 == grid->depth");
    update_depth (-1);
  }

#if _call_free_children
}
#define _IN_STENCIL 1

#line 991
static void _free_children (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 992 "/home/fpl/softwares/basilisk/src/grid/tree.h"


  Layer * L = ((Tree *)grid)->L[point.level + 1];
  int i = 2*point.i - 2, j = 2*point.j - 2;
  IF (!(((L->m)->b[i][j]))) _stencil_qassert (__FILE__,__LINE__,"/home/fpl/softwares/basilisk/src/grid/tree.h", 996, "mem_data (L->m,i,j)");
  mempool_free (L->pool, ((L->m)->b[i][j]));
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++)
      free_periodic (L->m, i + k, j + l, L->len);
  IF (--L->nc == 0) {
    destroy_layer (L);
    IF (!(point.level + 1 == grid->depth)) _stencil_qassert (__FILE__,__LINE__,"/home/fpl/softwares/basilisk/src/grid/tree.h", 1003, "point.level + 1 == grid->depth");
    update_depth (-1);
  }

#undef _IN_STENCIL

#endif

#line 1006
}
#line 1032 "/home/fpl/softwares/basilisk/src/grid/tree.h"
void increment_neighbors (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1033 "/home/fpl/softwares/basilisk/src/grid/tree.h"

  ((Tree *)grid)->dirty = true;
  if (cell.neighbors++ == 0)
    alloc_children (point);
   { foreach_neighbor (2/2)
    if (cell.neighbors++ == 0)
      alloc_children (point); end_foreach_neighbor(); }
  cell.neighbors--;

#if _call_increment_neighbors
}
#define _IN_STENCIL 1
#define alloc_children _alloc_children

#line 1032
static void _increment_neighbors (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1033 "/home/fpl/softwares/basilisk/src/grid/tree.h"

  ((Tree *)grid)->dirty = true;
  IF (cell.neighbors++ == 0)
    alloc_children (point);
   { foreach_neighbor (2/2)
    IF (cell.neighbors++ == 0)
      alloc_children (point); end_foreach_neighbor(); }
  cell.neighbors--;

#undef alloc_children
#undef _IN_STENCIL

#endif

#line 1041
}

void decrement_neighbors (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1044 "/home/fpl/softwares/basilisk/src/grid/tree.h"

  ((Tree *)grid)->dirty = true;
   { foreach_neighbor (2/2)
    if (allocated(0,0,0)) {
      cell.neighbors--;
      if (cell.neighbors == 0)
 free_children (point);
    } end_foreach_neighbor(); }
  if (cell.neighbors) {
    int pid = cell.pid;
     { foreach_child() {
      cell.flags = 0;
      cell.pid = pid;
    } end_foreach_child(); }
  }

#if _call_decrement_neighbors
}
#define _IN_STENCIL 1
#define free_children _free_children

#line 1043
static void _decrement_neighbors (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1044 "/home/fpl/softwares/basilisk/src/grid/tree.h"

  ((Tree *)grid)->dirty = true;
   { foreach_neighbor (2/2)
    IF (allocated(0,0,0)) {
      cell.neighbors--;
      IF (cell.neighbors == 0)
 free_children (point);
    } end_foreach_neighbor(); }
  IF (cell.neighbors) {
    int pid = cell.pid;
     { foreach_child() {
      cell.flags = 0;
      cell.pid = pid;
    } end_foreach_child(); }
  }

#undef free_children
#undef _IN_STENCIL

#endif

#line 1059
}

void realloc_scalar (int size)
{

  Tree * q = ((Tree *)grid);
  size_t oldlen = sizeof(Cell) + datasize;
  size_t newlen = oldlen + size;
  datasize += size;

  Layer * L = q->L[0];
   { foreach_mem (L->m, L->len, 1){

#line 1070 "/home/fpl/softwares/basilisk/src/grid/tree.h"
 {




    char * p = (char *) prealloc (((L->m)->b[point.i][point.j]),
     newlen*sizeof(char),__func__,__FILE__,__LINE__);
    assign_periodic (L->m, point.i, point.j, L->len, p);





  } } end_foreach_mem(); }

  for (int l = 1; l <= depth(); l++) {
    Layer * L = q->L[l];
    Mempool * oldpool = L->pool;
    L->pool = mempool_new (poolsize (l, newlen), (1 << 2)*newlen);
     { foreach_mem (L->m, L->len, 2){

#line 1089 "/home/fpl/softwares/basilisk/src/grid/tree.h"
 {
      char * new = (char *) mempool_alloc (L->pool);







      for (int k = 0; k < 2; k++)
 for (int o = 0; o < 2; o++) {
   memcpy (new, ((L->m)->b[point.i + k][point.j + o]), oldlen);
   assign_periodic (L->m, point.i + k, point.j + o, L->len, new);
   new += newlen;
 }
#line 1115 "/home/fpl/softwares/basilisk/src/grid/tree.h"
    } } end_foreach_mem(); }
    mempool_destroy (oldpool);
  }
}



#define VN v.x
#define VT v.y
#define VR v.z




#if 1
# define disable_fpe_for_mpi() disable_fpe (FE_DIVBYZERO|FE_INVALID)
# define enable_fpe_for_mpi() enable_fpe (FE_DIVBYZERO|FE_INVALID)
#else
# define disable_fpe_for_mpi()
# define enable_fpe_for_mpi()
#endif

static inline void no_restriction (Point point, scalar s);
#if _call_no_restriction
static void _no_restriction (Point point, scalar s);
#endif

#line 1137

static bool normal_neighbor (Point point, scalar * scalars, vector * vectors)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1140 "/home/fpl/softwares/basilisk/src/grid/tree.h"

  for (int k = 1; k <= 2; k++)
    {
#line 1142

      for (int i = -k; i <= k; i += 2*k)
 if ((allocated(i,0,0) && !(neighbor(i,0,0).pid < 0))) {
   Point neighbor = neighborp(i,0,0);
   int id = (- cell.pid - 1);
   strongif (scalars) for (scalar s = *scalars, *_i12 = scalars; ((scalar *)&s)->i >= 0; s = *++_i12)
    
       val(s,0,0,0) = _attribute[s.i].boundary[id](neighbor, point, s, NULL);
   strongif (vectors) for (vector v = *vectors, *_i13 = vectors; ((scalar *)&v)->i >= 0; v = *++_i13)
     {
       scalar vn = VN;
       val(v.x,0,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.x, NULL);

       scalar vt = VT;
       val(v.y,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.y, NULL);





     }
   return true;
 }
#line 1142

      for (int i = -k; i <= k; i += 2*k)
 if ((allocated(0,i,0) && !(neighbor(0,i,0).pid < 0))) {
   Point neighbor = neighborp(0,i,0);
   int id = (- cell.pid - 1);
   strongif (scalars) for (scalar s = *scalars, *_i12 = scalars; ((scalar *)&s)->i >= 0; s = *++_i12)
    
       val(s,0,0,0) = _attribute[s.i].boundary[id](neighbor, point, s, NULL);
   strongif (vectors) for (vector v = *vectors, *_i13 = vectors; ((scalar *)&v)->i >= 0; v = *++_i13)
     {
       scalar vn = VN;
       val(v.y,0,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.y, NULL);

       scalar vt = VT;
       val(v.x,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.x, NULL);





     }
   return true;
 }}
  return false;

#if _call_normal_neighbor
}
#define _IN_STENCIL 1

#line 1139
static bool _normal_neighbor (Point point, scalar * scalars, vector * vectors)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1140 "/home/fpl/softwares/basilisk/src/grid/tree.h"

  for (int k = 1; k <= 2; k++)
    {
#line 1142

      for (int i = -k; i <= k; i += 2*k)
 IF ((allocated(i,0,0) && !(neighbor(i,0,0).pid < 0))) {
   Point neighbor = neighborp(i,0,0);
   int id = (- cell.pid - 1);
   strongif (scalars) for (scalar s = *scalars, *_i12 = scalars; ((scalar *)&s)->i >= 0; s = *++_i12)
    
       _stencil_val(__FILE__,__LINE__,s,0,0,0) = _attribute[s.i].boundary[id](neighbor, point, s, NULL);
   strongif (vectors) for (vector v = *vectors, *_i13 = vectors; ((scalar *)&v)->i >= 0; v = *++_i13)
     {
       scalar vn = VN;
       _stencil_val(__FILE__,__LINE__,v.x,0,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.x, NULL);

       scalar vt = VT;
       _stencil_val(__FILE__,__LINE__,v.y,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.y, NULL);





     }
   return true;
 }
#line 1142

      for (int i = -k; i <= k; i += 2*k)
 IF ((allocated(0,i,0) && !(neighbor(0,i,0).pid < 0))) {
   Point neighbor = neighborp(0,i,0);
   int id = (- cell.pid - 1);
   strongif (scalars) for (scalar s = *scalars, *_i12 = scalars; ((scalar *)&s)->i >= 0; s = *++_i12)
    
       _stencil_val(__FILE__,__LINE__,s,0,0,0) = _attribute[s.i].boundary[id](neighbor, point, s, NULL);
   strongif (vectors) for (vector v = *vectors, *_i13 = vectors; ((scalar *)&v)->i >= 0; v = *++_i13)
     {
       scalar vn = VN;
       _stencil_val(__FILE__,__LINE__,v.y,0,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.y, NULL);

       scalar vt = VT;
       _stencil_val(__FILE__,__LINE__,v.x,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.x, NULL);





     }
   return true;
 }}
  return false;

#undef _IN_STENCIL

#endif

#line 1166
}

static bool diagonal_neighbor_2D (Point point,
      scalar * scalars, vector * vectors)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1170 "/home/fpl/softwares/basilisk/src/grid/tree.h"


  for (int k = 1; k <= 2; k++)



      for (int i = -k; i <= k; i += 2*k)
 for (int j = -k; j <= k; j += 2*k)
   if (allocated(i,j,0) && (allocated(i,j,0) && !(neighbor(i,j,0).pid < 0)) &&
       allocated(i,0,0) && (neighbor(i,0,0).pid < 0) &&
       allocated(0,j,0) && (neighbor(0,j,0).pid < 0)) {
     Point n = neighborp(i,j,0),
       n1 = neighborp(i,0,0), n2 = neighborp(0,j,0);
     int id1 = (- neighbor(i,0,0).pid - 1), id2 = (- neighbor(0,j,0).pid - 1);
     strongif (scalars) for (scalar s = *scalars, *_i14 = scalars; ((scalar *)&s)->i >= 0; s = *++_i14)
      
  val(s,0,0,0) = (_attribute[s.i].boundary[id1](n,n1,s,NULL) +
         _attribute[s.i].boundary[id2](n,n2,s,NULL) -
         val(s,i,j,0));
     strongif (vectors) for (vector v = *vectors, *_i15 = vectors; ((scalar *)&v)->i >= 0; v = *++_i15)
       {
  scalar vt = VT, vn = VN;
  val(v.x,0,0,0) = (_attribute[vt.i].boundary[id1](n,n1,v.x,NULL) +
    _attribute[vn.i].boundary[id2](n,n2,v.x,NULL) -
    val(v.x,i,j,0));
  val(v.y,0,0,0) = (_attribute[vn.i].boundary[id1](n,n1,v.y,NULL) +
    _attribute[vt.i].boundary[id2](n,n2,v.y,NULL) -
    val(v.y,i,j,0));






       }
     return true;
   }

  return false;

#if _call_diagonal_neighbor_2D
}
#define _IN_STENCIL 1

#line 1168
static bool _diagonal_neighbor_2D (Point point,
      scalar * scalars, vector * vectors)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1170 "/home/fpl/softwares/basilisk/src/grid/tree.h"


  for (int k = 1; k <= 2; k++)



      for (int i = -k; i <= k; i += 2*k)
 for (int j = -k; j <= k; j += 2*k)
   IF (allocated(i,j,0) && (allocated(i,j,0) && !(neighbor(i,j,0).pid < 0)) &&
       allocated(i,0,0) && (neighbor(i,0,0).pid < 0) &&
       allocated(0,j,0) && (neighbor(0,j,0).pid < 0)) {
     Point n = neighborp(i,j,0),
       n1 = neighborp(i,0,0), n2 = neighborp(0,j,0);
     int id1 = (- neighbor(i,0,0).pid - 1), id2 = (- neighbor(0,j,0).pid - 1);
     strongif (scalars) for (scalar s = *scalars, *_i14 = scalars; ((scalar *)&s)->i >= 0; s = *++_i14)
      
  _stencil_val(__FILE__,__LINE__,s,0,0,0) = (_attribute[s.i].boundary[id1](n,n1,s,NULL) +
         _attribute[s.i].boundary[id2](n,n2,s,NULL) -
         _stencil_val(__FILE__,__LINE__,s,i,j,0));
     strongif (vectors) for (vector v = *vectors, *_i15 = vectors; ((scalar *)&v)->i >= 0; v = *++_i15)
       {
  scalar vt = VT, vn = VN;
  _stencil_val(__FILE__,__LINE__,v.x,0,0,0) = (_attribute[vt.i].boundary[id1](n,n1,v.x,NULL) +
    _attribute[vn.i].boundary[id2](n,n2,v.x,NULL) -
    _stencil_val(__FILE__,__LINE__,v.x,i,j,0));
  _stencil_val(__FILE__,__LINE__,v.y,0,0,0) = (_attribute[vn.i].boundary[id1](n,n1,v.y,NULL) +
    _attribute[vt.i].boundary[id2](n,n2,v.y,NULL) -
    _stencil_val(__FILE__,__LINE__,v.y,i,j,0));






       }
     return true;
   }

  return false;

#undef _IN_STENCIL

#endif

#line 1209
}

static bool diagonal_neighbor_3D (Point point,
      scalar * scalars, vector * vectors)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1213 "/home/fpl/softwares/basilisk/src/grid/tree.h"

#line 1257 "/home/fpl/softwares/basilisk/src/grid/tree.h"
  return false;

#if _call_diagonal_neighbor_3D
}
#define _IN_STENCIL 1

#line 1211
static bool _diagonal_neighbor_3D (Point point,
      scalar * scalars, vector * vectors)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1213 "/home/fpl/softwares/basilisk/src/grid/tree.h"

#line 1257 "/home/fpl/softwares/basilisk/src/grid/tree.h"
  return false;

#undef _IN_STENCIL

#endif

#line 1258
}



#line 1261

static Point tangential_neighbor_x (Point point, bool * zn)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1263 "/home/fpl/softwares/basilisk/src/grid/tree.h"

  for (int k = 1; k <= 2; k++)
    for (int j = -k; j <= k; j += 2*k) {
      if ((allocated(0,j,0) && !(neighbor(0,j,0).pid < 0)) || (allocated(-1,j,0) && !(neighbor(-1,j,0).pid < 0))) {
 *zn = false;
 return neighborp(0,j,0);
      }







    }
  return (Point){.level = -1};

#if _call_tangential_neighbor_x
}
#define _IN_STENCIL 1

#line 1262
static Point _tangential_neighbor_x (Point point, bool * zn)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1263 "/home/fpl/softwares/basilisk/src/grid/tree.h"

  for (int k = 1; k <= 2; k++)
    for (int j = -k; j <= k; j += 2*k) {
      IF ((allocated(0,j,0) && !(neighbor(0,j,0).pid < 0)) || (allocated(-1,j,0) && !(neighbor(-1,j,0).pid < 0))) {
 *zn = false;
 return neighborp(0,j,0);
      }







    }
  return (Point){.level = -1};

#undef _IN_STENCIL

#endif

#line 1279
}
#line 1261

static Point tangential_neighbor_y (Point point, bool * zn)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1263 "/home/fpl/softwares/basilisk/src/grid/tree.h"

  for (int k = 1; k <= 2; k++)
    for (int j = -k; j <= k; j += 2*k) {
      if ((allocated(j,0,0) && !(neighbor(j,0,0).pid < 0)) || (allocated(j,-1,0) && !(neighbor(j,-1,0).pid < 0))) {
 *zn = false;
 return neighborp(j,0,0);
      }







    }
  return (Point){.level = -1};

#if _call_tangential_neighbor_y
}
#define _IN_STENCIL 1

#line 1262
static Point _tangential_neighbor_y (Point point, bool * zn)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1263 "/home/fpl/softwares/basilisk/src/grid/tree.h"

  for (int k = 1; k <= 2; k++)
    for (int j = -k; j <= k; j += 2*k) {
      IF ((allocated(j,0,0) && !(neighbor(j,0,0).pid < 0)) || (allocated(j,-1,0) && !(neighbor(j,-1,0).pid < 0))) {
 *zn = false;
 return neighborp(j,0,0);
      }







    }
  return (Point){.level = -1};

#undef _IN_STENCIL

#endif

#line 1279
}


static inline bool is_boundary_point (Point point) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1282 "/home/fpl/softwares/basilisk/src/grid/tree.h"

  return (cell.pid < 0);

#if _call_is_boundary_point
}
#define _IN_STENCIL 1

#line 1282
static bool _is_boundary_point (Point point) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1282 "/home/fpl/softwares/basilisk/src/grid/tree.h"

  return (cell.pid < 0);

#undef _IN_STENCIL

#endif

#line 1284
}

static void box_boundary_level (const Boundary * b, scalar * list, int l)
{
  disable_fpe_for_mpi();
  scalar * scalars = NULL;
  vector * vectors = NULL, * faces = NULL;
  strongif (list) for (scalar s = *list, *_i16 = list; ((scalar *)&s)->i >= 0; s = *++_i16)
    if (!is_constant(s) && _attribute[s.i].refine != no_restriction) {
      if (_attribute[s.i].v.x.i == s.i) {
 if (_attribute[s.i].face)
   faces = vectors_add (faces, _attribute[s.i].v);
 else
   vectors = vectors_add (vectors, _attribute[s.i].v);
      }
      else if (_attribute[s.i].v.x.i < 0 && _attribute[s.i].boundary[0])
 scalars = list_add (scalars, s);
    }

   { foreach_boundary_level (l){

#line 1303 "/home/fpl/softwares/basilisk/src/grid/tree.h"
 {
    if (!normal_neighbor (point, scalars, vectors) &&
 !diagonal_neighbor_2D (point, scalars, vectors) &&
 !diagonal_neighbor_3D (point, scalars, vectors)) {

      strongif (scalars) for (scalar s = *scalars, *_i17 = scalars; ((scalar *)&s)->i >= 0; s = *++_i17)

   val(s,0,0,0) = undefined;
      strongif (vectors) for (vector v = *vectors, *_i18 = vectors; ((scalar *)&v)->i >= 0; v = *++_i18)

   {
#line 1313

     val(v.x,0,0,0) = undefined;
#line 1313

     val(v.y,0,0,0) = undefined;}
    }
    if (faces) {
      int id = (- cell.pid - 1);
      {
#line 1318

 for (int i = -1; i <= 1; i += 2) {

   if ((allocated(i,0,0) && !(neighbor(i,0,0).pid < 0))) {
     Point neighbor = neighborp(i,0,0);
     strongif (faces) for (vector v = *faces, *_i19 = faces; ((scalar *)&v)->i >= 0; v = *++_i19) {
       scalar vn = VN;
       if (_attribute[vn.i].boundary[id])
 
    val(v.x,(i + 1)/2,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.x, NULL);
     }
   }

   else if (i == -1) {

     bool zn;
     Point neighbor = tangential_neighbor_x (point, &zn);
     if (neighbor.level >= 0) {
       int id = is_boundary_point (neighbor) ?
  (- neighbor(-1,0,0).pid - 1) : (- cell.pid - 1);
       strongif (faces) for (vector v = *faces, *_i20 = faces; ((scalar *)&v)->i >= 0; v = *++_i20) {

  scalar vt = VT;



 
    val(v.x,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.x, NULL);
       }
     }
     else

       strongif (faces) for (vector v = *faces, *_i21 = faces; ((scalar *)&v)->i >= 0; v = *++_i21)
 
    val(v.x,0,0,0) = 0.;
   }

 }
#line 1318

 for (int i = -1; i <= 1; i += 2) {

   if ((allocated(0,i,0) && !(neighbor(0,i,0).pid < 0))) {
     Point neighbor = neighborp(0,i,0);
     strongif (faces) for (vector v = *faces, *_i19 = faces; ((scalar *)&v)->i >= 0; v = *++_i19) {
       scalar vn = VN;
       if (_attribute[vn.i].boundary[id])
 
    val(v.y,0,(i + 1)/2,0) = _attribute[vn.i].boundary[id](neighbor, point, v.y, NULL);
     }
   }

   else if (i == -1) {

     bool zn;
     Point neighbor = tangential_neighbor_y (point, &zn);
     if (neighbor.level >= 0) {
       int id = is_boundary_point (neighbor) ?
  (- neighbor(0,-1,0).pid - 1) : (- cell.pid - 1);
       strongif (faces) for (vector v = *faces, *_i20 = faces; ((scalar *)&v)->i >= 0; v = *++_i20) {

  scalar vt = VT;



 
    val(v.y,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.y, NULL);
       }
     }
     else

       strongif (faces) for (vector v = *faces, *_i21 = faces; ((scalar *)&v)->i >= 0; v = *++_i21)
 
    val(v.y,0,0,0) = 0.;
   }

 }}
    }
  } } end_foreach_boundary_level(); }

  pfree (scalars,__func__,__FILE__,__LINE__);
  pfree (vectors,__func__,__FILE__,__LINE__);
  pfree (faces,__func__,__FILE__,__LINE__);
  enable_fpe_for_mpi();
}



#undef VN
#undef VT
#define VN _attribute[s.i].v.x
#define VT _attribute[s.i].v.y

static double masked_average (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1373 "/home/fpl/softwares/basilisk/src/grid/tree.h"

  double sum = 0., n = 0.;
   { foreach_child()
    if (!(cell.pid < 0) && val(s,0,0,0) != nodata)
      sum += val(s,0,0,0), n++; end_foreach_child(); }
  return n ? sum/n : nodata;

#if _call_masked_average
}
#define _IN_STENCIL 1

#line 1372
static double _masked_average (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1373 "/home/fpl/softwares/basilisk/src/grid/tree.h"

  double sum = 0., n = 0.;
   { foreach_child()
    IF (!(cell.pid < 0) && _stencil_val(__FILE__,__LINE__,s,0,0,0) != nodata)
      sum += _stencil_val(__FILE__,__LINE__,s,0,0,0), n++; end_foreach_child(); }
  return n ? sum/n : nodata;

#undef _IN_STENCIL

#endif

#line 1379
}


#line 1381

static double masked_average_x (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1383 "/home/fpl/softwares/basilisk/src/grid/tree.h"

  double sum = 0., n = 0.;
   { foreach_child()
    if (child.x < 0 && (!(cell.pid < 0) || !(neighbor(1,0,0).pid < 0)) &&
 val(s,1,0,0) != nodata)
      sum += val(s,1,0,0), n++; end_foreach_child(); }
  return n ? sum/n : nodata;

#if _call_masked_average_x
}
#define _IN_STENCIL 1

#line 1382
static double _masked_average_x (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1383 "/home/fpl/softwares/basilisk/src/grid/tree.h"

  double sum = 0., n = 0.;
   { foreach_child()
    IF (child.x < 0 && (!(cell.pid < 0) || !(neighbor(1,0,0).pid < 0)) &&
 _stencil_val(__FILE__,__LINE__,s,1,0,0) != nodata)
      sum += _stencil_val(__FILE__,__LINE__,s,1,0,0), n++; end_foreach_child(); }
  return n ? sum/n : nodata;

#undef _IN_STENCIL

#endif

#line 1390
}
#line 1381

static double masked_average_y (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1383 "/home/fpl/softwares/basilisk/src/grid/tree.h"

  double sum = 0., n = 0.;
   { foreach_child()
    if (child.y < 0 && (!(cell.pid < 0) || !(neighbor(0,1,0).pid < 0)) &&
 val(s,0,1,0) != nodata)
      sum += val(s,0,1,0), n++; end_foreach_child(); }
  return n ? sum/n : nodata;

#if _call_masked_average_y
}
#define _IN_STENCIL 1

#line 1382
static double _masked_average_y (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1383 "/home/fpl/softwares/basilisk/src/grid/tree.h"

  double sum = 0., n = 0.;
   { foreach_child()
    IF (child.y < 0 && (!(cell.pid < 0) || !(neighbor(0,1,0).pid < 0)) &&
 _stencil_val(__FILE__,__LINE__,s,1,0,0) != nodata)
      sum += _stencil_val(__FILE__,__LINE__,s,1,0,0), n++; end_foreach_child(); }
  return n ? sum/n : nodata;

#undef _IN_STENCIL

#endif

#line 1390
}

static void masked_boundary_restriction (const Boundary * b,
      scalar * list, int l)
{
  scalar * scalars = NULL;
  vector * faces = NULL;
  strongif (list) for (scalar s = *list, *_i22 = list; ((scalar *)&s)->i >= 0; s = *++_i22)
    if (!is_constant(s) && _attribute[s.i].refine != no_restriction) {
      if (_attribute[s.i].v.x.i == s.i && _attribute[s.i].face)
 faces = vectors_add (faces, _attribute[s.i].v);
      else
 scalars = list_add (scalars, s);
    }

   { foreach_halo (restriction, l){

#line 1405 "/home/fpl/softwares/basilisk/src/grid/tree.h"
 {
    strongif (scalars) for (scalar s = *scalars, *_i23 = scalars; ((scalar *)&s)->i >= 0; s = *++_i23)
      val(s,0,0,0) = masked_average (parent, s);
    strongif (faces) for (vector v = *faces, *_i24 = faces; ((scalar *)&v)->i >= 0; v = *++_i24)
      {
#line 1409
 {
 double average = masked_average_x (parent, v.x);
 if ((neighbor(-1,0,0).pid < 0))
   val(v.x,0,0,0) = average;
 if ((neighbor(1,0,0).pid < 0))
   val(v.x,1,0,0) = average;
      }
#line 1409
 {
 double average = masked_average_y (parent, v.y);
 if ((neighbor(0,-1,0).pid < 0))
   val(v.y,0,0,0) = average;
 if ((neighbor(0,1,0).pid < 0))
   val(v.y,0,1,0) = average;
      }}
  } } end_foreach_halo(); }

  pfree (scalars,__func__,__FILE__,__LINE__);
  pfree (faces,__func__,__FILE__,__LINE__);
}
#line 1445 "/home/fpl/softwares/basilisk/src/grid/tree.h"
static void free_cache (CacheLevel * c)
{
  for (int l = 0; l <= depth(); l++)
    pfree (c[l].p,__func__,__FILE__,__LINE__);
  pfree (c,__func__,__FILE__,__LINE__);
}

void free_grid (void)
{
  if (!grid)
    return;
  free_boundaries();
  Tree * q = ((Tree *)grid);
  pfree (q->leaves.p,__func__,__FILE__,__LINE__);
  pfree (q->faces.p,__func__,__FILE__,__LINE__);
  pfree (q->vertices.p,__func__,__FILE__,__LINE__);
  pfree (q->refined.p,__func__,__FILE__,__LINE__);


  Layer * L = q->L[0];
   { foreach_mem (L->m, L->len, 1){

#line 1465 "/home/fpl/softwares/basilisk/src/grid/tree.h"
 {



    pfree (((L->m)->b[point.i][point.j]),__func__,__FILE__,__LINE__);



  } } end_foreach_mem(); }
  for (int l = 0; l <= depth(); l++)
    destroy_layer (q->L[l]);
  q->L = &(q->L[-1]);
  pfree (q->L,__func__,__FILE__,__LINE__);
  free_cache (q->active);
  free_cache (q->prolongation);
  free_cache (q->boundary);
  free_cache (q->restriction);
  pfree (q,__func__,__FILE__,__LINE__);
  grid = NULL;
}

static void refine_level (int depth);


void init_grid (int n)
{ trace ("init_grid", "/home/fpl/softwares/basilisk/src/grid/tree.h", 1490);

  if (!(sizeof(Cell) % 8 == 0)) qassert ("/home/fpl/softwares/basilisk/src/grid/tree.h", 1492, "sizeof(Cell) % 8 == 0");

  free_grid();
  int depth = 0;
  while (n > 1) {
    if (n % 2) {
      fprintf (ferr, "tree: N must be a power-of-two\n");
      exit (1);
    }
    n /= 2;
    depth++;
  }
  Tree * q = ((Tree *) pcalloc (1, sizeof(Tree),__func__,__FILE__,__LINE__));
  grid = (Grid *) q;
  grid->depth = 0;


  q->L = ((Layer * *) pmalloc ((2)*sizeof(Layer *),__func__,__FILE__,__LINE__));

  q->L[0] = NULL; q->L = &(q->L[1]);

  Layer * L = new_layer (0);
  q->L[0] = L;
#line 1528 "/home/fpl/softwares/basilisk/src/grid/tree.h"
  for (int i = Period.x*2; i < L->len - Period.x*2; i++)
    for (int j = Period.y*2; j < L->len - Period.y*2; j++)
      assign_periodic (L->m, i, j, L->len,
         (char *) pcalloc (1, sizeof(Cell) + datasize,__func__,__FILE__,__LINE__));
  CELL(((L->m)->b[2][2])).flags |= leaf;
  if (pid() == 0)
    CELL(((L->m)->b[2][2])).flags |= active;
  for (int k = - 2*(1 - Period.x); k <= 2*(1 - Period.x); k++)
    for (int l = -2*(1 - Period.y); l <= 2*(1 - Period.y); l++)
      CELL(((L->m)->b[2 +k][2 +l])).pid =
 (k < 0 ? -1 - left :
  k > 0 ? -1 - right :
  l > 0 ? -1 - top :
  l < 0 ? -1 - bottom :
  0);
  CELL(((L->m)->b[2][2])).pid = 0;
#line 1566 "/home/fpl/softwares/basilisk/src/grid/tree.h"
  q->active = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));
  q->prolongation = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));
  q->boundary = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));
  q->restriction = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));
  q->dirty = true;
  N = 1 << depth;
#if 1
  void mpi_boundary_new();
  mpi_boundary_new();
#endif

  Boundary * b = ((Boundary *) pcalloc (1, sizeof(Boundary),__func__,__FILE__,__LINE__));
  b->level = box_boundary_level;
  b->restriction = masked_boundary_restriction;
  add_boundary (b);
  refine_level (depth);
  reset (all, 0.);
  { if (((Tree *)grid)->dirty) update_cache_f(); };
 end_trace("init_grid", "/home/fpl/softwares/basilisk/src/grid/tree.h", 1584); }


void check_two_one (void)
{
   { foreach_leaf(){

#line 1589 "/home/fpl/softwares/basilisk/src/grid/tree.h"

    if (level > 0)
      for (int k = -1; k <= 1; k++)
 for (int l = -1; l <= 1; l++) {

   int i = (point.i + 2)/2 + k;
   int j = (point.j + 2)/2 + l;
   double x = ((i - 2 + 0.5)*(1./(1 << point.level))*2. - 0.5);
   double y = ((j - 2 + 0.5)*(1./(1 << point.level))*2. - 0.5);
   if (x > -0.5 && x < 0.5 && y > -0.5 && y < 0.5 &&
       !(aparent(k,l,0).flags & active)) {
     FILE * fp = fopen("check_two_one_loc", "w");
     fprintf (fp,
       "# %d %d\n"
       "%g %g\n%g %g\n",
       k, l,
       (((point.i - 2) + 0.5)*(1./(1 << point.level)) - 0.5),
       (((point.j - 2) + 0.5)*(1./(1 << point.level)) - 0.5),
       x, y);
     fclose (fp);





     if (!(false)) qassert ("/home/fpl/softwares/basilisk/src/grid/tree.h", 1614, "false");
   }
 } } end_foreach_leaf(); }
}


struct _locate { double x, y, z; };

Point locate (struct _locate p)
{
  for (int l = depth(); l >= 0; l--) {
    Point point = {0};  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1625 "/home/fpl/softwares/basilisk/src/grid/tree.h"

    point.level = l;
    int n = 1 << point.level;
    point.i = (p.x - X0)/L0*n + 2;

    point.j = (p.y - Y0)/L0*n + 2;




    if (point.i >= 0 && point.i < n + 2*2

 && point.j >= 0 && point.j < n + 2*2




 ) {
      if (allocated(0,0,0) && is_local(cell) && is_leaf(cell))
 return point;
    }
    else
      break;
  }
  Point point = {0};  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 1649 "/home/fpl/softwares/basilisk/src/grid/tree.h"

  point.level = -1;
  return point;
}



bool tree_is_full ()
{
  { if (((Tree *)grid)->dirty) update_cache_f(); };
  return (grid->tn == 1L << grid->maxdepth*2);
}

#line 1 "grid/tree-common.h"
#line 1 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"



#line 1 "grid/multigrid-common.h"
#line 1 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"


#line 1 "grid/cartesian-common.h"
#line 1 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"
#line 1 "grid/events.h"
#line 1 "/home/fpl/softwares/basilisk/src/grid/events.h"







static void event_error (Event * ev, const char * s)
{
  fprintf (ferr, "%s:%d: error: %s\n", ev->file, ev->line, s);
  exit (1);
}

static void init_event (Event * ev)
{
  if (ev->arrayi || ev->arrayt) {
    ev->i = ev->t = -1;
    if (ev->arrayi)
      ev->i = ev->arrayi[0];
    else
      ev->t = ev->arrayt[0];
    ev->a = 1;
    ev->expr[1] = NULL;
  }
  else {
    if (ev->nexpr > 0) {
      Expr init = NULL, cond = NULL, inc = NULL;
      for (int j = 0; j < ev->nexpr; j++) {
 int i = -123456; double t = i;
 (* ev->expr[j]) (&i, &t, ev);
 if (i == -123456 && t == -123456) {

   if (cond)
     event_error (ev, "events can only use a single condition");
   cond = ev->expr[j];
 }
 else {

   int i1 = i; double t1 = t;
   (* ev->expr[j]) (&i1, &t1, ev);
   if (i1 == i && t1 == t) {


     if (init)
       event_error (ev, "events can only use a single initialisation");
     init = ev->expr[j];
   }
   else {

     if (inc)
       event_error (ev, "events can only use a single increment");
     inc = ev->expr[j];
   }
 }
      }
      ev->expr[0] = init;
      ev->expr[1] = cond;
      ev->expr[2] = inc;
      ev->nexpr = 0;
    }
    ev->i = ev->t = -1;
    if (ev->expr[0]) {
      (* ev->expr[0]) (&ev->i, &ev->t, ev);
      if (ev->i == 1234567890 || ev->t == 1234567890) {
 ev->i = 1234567890; ev->t = -1;
      }
    }
    else if (ev->expr[2]) {
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (ev->i != -1)
 ev->i = 0;
      if (ev->t != -1)
 ev->t = 0;
    }
  }
}

enum { event_done, event_alive, event_stop };

static int event_finished (Event * ev)
{
  ev->t = ev->i = -1;
  return event_done;
}

void event_register (Event event) {
  if (!(Events)) qassert ("/home/fpl/softwares/basilisk/src/grid/events.h", 87, "Events");
  if (!(!event.last)) qassert ("/home/fpl/softwares/basilisk/src/grid/events.h", 88, "!event.last");
  int n = 0, parent = -1;
  for (Event * ev = Events; !ev->last; ev++) {
    if (!strcmp (event.name, ev->name)) {
      if (!(parent < 0)) qassert ("/home/fpl/softwares/basilisk/src/grid/events.h", 92, "parent < 0");
      parent = n;
    }
    n++;
  }
  if (parent < 0) {
    Events = (Event *) prealloc (Events, (n + 2)*sizeof(Event),__func__,__FILE__,__LINE__);
    Events[n] = event;
    Events[n].next = NULL;
    Events[n + 1].last = true;
    init_event (&Events[n]);
  }
  else {
    Event * ev = ((Event *) pcalloc (1, sizeof(Event),__func__,__FILE__,__LINE__));
    *ev = Events[parent];
    Events[parent] = event;
    Events[parent].next = ev;
    init_event (&Events[parent]);
  }
}

static int event_cond (Event * ev, int i, double t)
{
  if (!ev->expr[1])
    return true;
  return (* ev->expr[1]) (&i, &t, ev);
}
#line 131 "/home/fpl/softwares/basilisk/src/grid/events.h"
static int event_do (Event * ev, bool action)
{
  if ((iter > ev->i && t > ev->t) || !event_cond (ev, iter, t))
    return event_finished (ev);
  if (iter == ev->i || fabs (t - ev->t) <= 1e-9) {
    if (action) {
      bool finished = false;
      for (Event * e = ev; e; e = e->next) {



 if ((* e->action) (iter, t, e))
   finished = true;
      }
      if (finished) {
 event_finished (ev);
 return event_stop;
      }
    }
    if (ev->arrayi) {
      ev->i = ev->arrayi[ev->a++];
      if (ev->i < 0)
 return event_finished (ev);
    }
    if (ev->arrayt) {
      ev->t = ev->arrayt[ev->a++];
      if (ev->t < 0)
 return event_finished (ev);
    }
    else if (ev->expr[2]) {
      int i0 = ev->i;
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (i0 == -1 && ev->i != i0)
 ev->i += iter + 1;
      if (!event_cond (ev, iter + 1, ev->t))
 return event_finished (ev);
    }
    else if (ev->expr[0] && !ev->expr[1])
      return event_finished (ev);
  }
  return event_alive;
}

static void end_event_do (bool action)
{




  for (Event * ev = Events; !ev->last; ev++)
    if (ev->i == 1234567890 && action)
      for (Event * e = ev; e; e = e->next) {



 e->action (iter, t, e);
      }
}

int events (bool action)
{





  if (iter == 0)
    for (Event * ev = Events; !ev->last; ev++)
      init_event (ev);

  int cond = 0, cond1 = 0;
  inext = 1234567890; tnext = HUGE;
  for (Event * ev = Events; !ev->last && !cond; ev++)
    if (ev->i != 1234567890 &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond = 1;
  for (Event * ev = Events; !ev->last; ev++) {
    int status = event_do (ev, action);
    if (status == event_stop) {
      end_event_do (action);
      return 0;
    }
    if (status == event_alive && ev->i != 1234567890 &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond1 = 1;
    if (ev->t > t && ev->t < tnext)
      tnext = ev->t;
    if (ev->i > iter && ev->i < inext)
      inext = ev->i;
  }
  if ((!cond || cond1) && (tnext != HUGE || inext != 1234567890)) {
    inext = iter + 1;
    return 1;
  }
  end_event_do (action);
  return 0;
}

void event (const char * name)
{
  for (Event * ev = Events; !ev->last; ev++)
    if (!strcmp (ev->name, name))
      for (Event * e = ev; e; e = e->next) {



 (* e->action) (0, 0, e);
      }
}

double dtnext (double dt)
{
  if (tnext != HUGE && tnext > t) {
    unsigned int n = (tnext - t)/dt;
    if (!(n < INT_MAX)) qassert ("/home/fpl/softwares/basilisk/src/grid/events.h", 245, "n < INT_MAX");
    if (n == 0)
      dt = tnext - t;
    else {
      double dt1 = (tnext - t)/n;
      if (dt1 > dt + 1e-9)
 dt = (tnext - t)/(n + 1);
      else if (dt1 < dt)
 dt = dt1;
      tnext = t + dt;
    }
  }
  else
    tnext = t + dt;
  return dt;
}
#line 2 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"

void (* debug) (Point);

#define _val_constant(a,k,l,m) ((const double) _constant[a.i -_NVARMAX])
#define diagonalize(a)
#define val_diagonal(a,k,l,m) ((k) == 0 && (l) == 0 && (m) == 0)

#undef VARIABLES
#define VARIABLES\
  double Delta = L0*(1./(1 << point.level));\
  \
    double Delta_x = Delta;\
    double Delta_y = Delta;\
\
  double x = (ig/2. + (point.i - 2) + 0.5)*Delta + X0; NOT_UNUSED(x);\
\
  double y = (jg/2. + (point.j - 2) + 0.5)*Delta + Y0;\
\
\
\
 NOT_UNUSED(y);\
\
\
\
  double z = 0.;\
\
  NOT_UNUSED(z);\
\
  NOT_UNUSED(Delta);\
  \
    NOT_UNUSED(Delta_x);\
    NOT_UNUSED(Delta_y);\
\
  ;\

#line 34


#line 1 "grid/fpe.h"
#line 1 "/home/fpl/softwares/basilisk/src/grid/fpe.h"


#include <signal.h>
#include <unistd.h>

static int gdb ()
{
  if (last_point.level >= 0) {
    debug (last_point);
    fputc ('\n', ferr);
    fflush (ferr);
  }
  char command[80];
  sprintf (command, "exec xterm -e 'gdb -p %d' & xterm -e 'gnuplot plot -'",
    getpid());
  return system (command);
}

static void caught_abort (int sig)
{
  fprintf (ferr, "Caught signal %d (Aborted)\n", sig);
  gdb();
}

static void caught_fpe (int sig)
{
  fprintf (ferr, "Caught signal %d (Floating Point Exception)\n", sig);
  gdb();
  exit (1);
}

static void caught_segfault (int sig)
{
  fprintf (ferr, "Caught signal %d (Segmentation Fault)\n", sig);
  gdb();
  exit (2);
}

void catch_fpe (void)
{
  struct sigaction act;
  act.sa_handler = caught_fpe;
  sigemptyset (&act.sa_mask);
  act.sa_flags = 0;
  last_point.level = -1;
  sigaction (8, &act, NULL);
  act.sa_handler = caught_segfault;
  sigaction (11, &act, NULL);
  act.sa_handler = caught_abort;
  act.sa_flags = SA_RESETHAND;
  sigaction (6, &act, NULL);
}
#line 37 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"

#define end_foreach_face()

#line 1 "grid/stencils.h"
#line 1 "/home/fpl/softwares/basilisk/src/grid/stencils.h"
#line 18 "/home/fpl/softwares/basilisk/src/grid/stencils.h"








#line 43 "/home/fpl/softwares/basilisk/src/grid/stencils.h"
typedef struct {
  const char * fname;
  int line;
  int first;
  const char * each;
  int face;
  bool vertex;
} ForeachData;
#line 59 "/home/fpl/softwares/basilisk/src/grid/stencils.h"
static inline int stencil_index (scalar s, IJK i)
{
  int len = 1, index = 0;
  for (int d = 0; d < 2; d++) {
    if (i.i[d] < 0 || i.i[d] >= 5)
      return -1;
    index += len*i.i[d], len *= 5;
  }
  return index;
}




static inline IJK stencil_ijk (int index)
{
  IJK i;
  int len = (5*5);
  for (int d = 2 - 1; d >= 0; d--) {
    len /= 5;
    i.i[d] = index/len;
    index -= len*i.i[d];
    i.i[d] -= 5/2;
  }
  return i;
}




static void write_stencil_index (IJK i, int shift)
{
  sysfprintf (qstderr(), "[%d", i.i[0] - shift);
  for (int d = 1; d < 2; d++)
    sysfprintf (qstderr(), ",%d", i.i[d] - shift);
  sysfprintf (qstderr(), "]");
}






int _stencil_access (scalar s, IJK i, const char * file, int line)
{
  if (is_constant(s) || s.i < 0)
    return 0;
  int index = stencil_index (s, i);
  if (index < 0) {
    sysfprintf (qstderr(), "%s:%d: error: stencil overflow: %s",
  file, line, _attribute[s.i].name);
    write_stencil_index (i, 5/2);
    sysfprintf (qstderr(), "\n");
    fflush (qstderr());
    abort();
  }
  _attribute[s.i].read[index]++;
  return index;
}
#line 127 "/home/fpl/softwares/basilisk/src/grid/stencils.h"
#define foreach_stencil() {\
  strongif (baseblock) for (scalar _s = *baseblock, *_i25 = baseblock; ((scalar *)&_s)->i >= 0; _s = *++_i25) {\
    for (int _i = 0; _i < (5*5); _i++) {\
      _attribute[_s.i].read[_i] = 0;\
      _attribute[_s.i].write[_i] = 1.7759437274 + _s.i + _i;\
    }\
  }\
  _foreach_data.face = 0;\
  int ig = 0, jg = 0, kg = 0;\
  NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
\
\
\
\
  for (int _i = 0; _i < 2; _i++)\
    ((int *)&point)[_i] = 5/2;\
  if (sizeof(Point) >= (2 + 2)*sizeof(int))\
    ((int *)&point)[2 + 1] = 1;\
  POINT_VARIABLES\

#line 147






#define end_foreach_stencil()\
  strongif (baseblock) for (scalar _s = *baseblock, *_i26 = baseblock; ((scalar *)&_s)->i >= 0; _s = *++_i26) {\
    for (int _i = 0; _i < (5*5); _i++)\
      _attribute[_s.i].write[_i] = (_attribute[_s.i].write[_i] != 1.7759437274 + _s.i + _i);\
  }\
  end_stencil (&_foreach_data);\
}\

#line 160


#define foreach_face_stencil foreach_stencil
#define end_foreach_face_stencil() end_foreach_stencil()

#define foreach_vertex_stencil foreach_stencil
#define end_foreach_vertex_stencil() end_foreach_stencil()

#define is_stencil_face_x() ((_foreach_data.face |= (1 << 0)))
#define is_stencil_face_y() ((_foreach_data.face |= (1 << 1)))
#define is_stencil_face_z() ((_foreach_data.face |= (1 << 2)))




void reduction_warning (const char * fname, int line, const char * var)
{
  fprintf (ferr,
  "%s:%d: warning: variable '%s' is modified by this foreach loop:\n"
  "%s:%d: warning: use a loop-local variable, a reduction operation\n"
  "%s:%d: warning: or a serial loop to get rid of this warning\n",
    fname, line, var, fname, line, fname, line);
}







static inline bool scalar_is_dirty (scalar s)
{
  if (_attribute[s.i].dirty)
    return true;
  scalar * depends = _attribute[s.i].depends;
  strongif (depends) for (scalar d = *depends, *_i27 = depends; ((scalar *)&d)->i >= 0; d = *++_i27)
    if (_attribute[d.i].dirty)
      return true;
  return false;
}




static inline bool scalar_depends_from (scalar a, scalar b)
{
  scalar * depends = _attribute[a.i].depends;
  strongif (depends) for (scalar s = *depends, *_i28 = depends; ((scalar *)&s)->i >= 0; s = *++_i28)
    if (s.i == b.i)
      return true;
  return false;
}







void boundary_internal (scalar * list, const char * fname, int line);
void (* boundary_face) (vectorl);







void end_stencil (ForeachData * loop)
{
  scalar * listc = NULL, * dirty = NULL;
  vectorl listf = {NULL};
  bool flux = false;
  loop->vertex = !strcmp (loop->each, "foreach_vertex");




  strongif (baseblock) for (scalar s = *baseblock, *_i29 = baseblock; ((scalar *)&s)->i >= 0; s = *++_i29) {
    bool write = false, read = false;
    int max = 0;




    for (int n = 0; n < (5*5); n++)
      if (_attribute[s.i].write[n] || _attribute[s.i].read[n]) {
 IJK i = stencil_ijk (n);





 if (_attribute[s.i].write[n]) {
   for (int d = 0; d < 2; d++)
     if (i.i[d] != 0) {
       fprintf (ferr,
         "%s:%d: error: illegal write within this loop: %s",
         loop->fname, loop->line, _attribute[s.i].name);
       write_stencil_index (i, 0);
       fprintf (ferr, "\n");
       fflush (ferr);
       abort();
     }
   write = true;
 }





 if (_attribute[s.i].read[n]) {
   read = true;
   int d = 0;
   {
#line 274
 {
     if ((!_attribute[s.i].face || _attribute[s.i].v.x.i != s.i) && abs(i.i[d]) > max)
       max = abs(i.i[d]);
     d++;
   }
#line 274
 {
     if ((!_attribute[s.i].face || _attribute[s.i].v.y.i != s.i) && abs(i.i[d]) > max)
       max = abs(i.i[d]);
     d++;
   }}
 }
      }







    {





      if (read && scalar_is_dirty (s)) {





 if (_attribute[s.i].face) {
   if (max > 0)
     listc = list_append (listc, s);
   else if (!write) {
     scalar sn = _attribute[s.i].v.x.i >= 0 ? _attribute[s.i].v.x : s;
     {
#line 305

       if (_attribute[s.i].v.x.i == s.i) {




  if (_attribute[sn.i].boundary[left] || _attribute[sn.i].boundary[right])
    listc = list_append (listc, s);
  else if (_attribute[s.i].dirty != 2) {
    listf.x = list_append (listf.x, s);
    flux = true;
  }
       }
#line 305

       if (_attribute[s.i].v.y.i == s.i) {




  if (_attribute[sn.i].boundary[bottom] || _attribute[sn.i].boundary[top])
    listc = list_append (listc, s);
  else if (_attribute[s.i].dirty != 2) {
    listf.y = list_append (listf.y, s);
    flux = true;
  }
       }}
   }
 }





 else if (max > 0)
   listc = list_append (listc, s);
      }





      if (write) {
 if (2 > 1 && !loop->vertex && loop->first) {
   bool vertex = true;
   {
#line 336

     if (_attribute[s.i].d.x != -1)
       vertex = false;
#line 336

     if (_attribute[s.i].d.y != -1)
       vertex = false;}
   if (vertex)
     fprintf (ferr,
       "%s:%d: warning: vertex scalar '%s' should be assigned with"
       " a foreach_vertex() loop\n",
       loop->fname, loop->line, _attribute[s.i].name);
 }
 if (_attribute[s.i].face) {
   if (loop->face == 0 && loop->first)
     fprintf (ferr,
       "%s:%d: warning: face vector '%s' should be assigned with"
       " a foreach_face() loop\n",
       loop->fname, loop->line, _attribute[s.i].name);
 }
 else if (loop->face) {
   if (_attribute[s.i].v.x.i < 0) {
     int d = 1, i = 0;
     {
#line 355
 {
       if (loop->face == d) {
  _attribute[s.i].face = 2, _attribute[s.i].v.x.i = s.i;
  _attribute[s.i].boundary[left] = _attribute[s.i].boundary[right] = NULL;





       }
       d *= 2, i++;
     }
#line 355
 {
       if (loop->face == d) {
  _attribute[s.i].face = 2, _attribute[s.i].v.y.i = s.i;
  _attribute[s.i].boundary[bottom] = _attribute[s.i].boundary[top] = NULL;





       }
       d *= 2, i++;
     }}
     if (!_attribute[s.i].face && loop->first)
       fprintf (ferr,
         "%s:%d: warning: scalar '%s' should be assigned with "
         "a foreach_face(x|y|z) loop\n",
         loop->fname, loop->line, _attribute[s.i].name);
   }
   else {
     char * name = NULL;
     if (_attribute[s.i].name) {
       name = pstrdup (_attribute[s.i].name,__func__,__FILE__,__LINE__);
       char * s = name + strlen(name) - 1;
       while (s != name && *s != '.') s--;
       if (s != name) *s = '\0';
     }
     init_face_vector (_attribute[s.i].v, name);




     pfree (name,__func__,__FILE__,__LINE__);
   }
 }
 else if (loop->vertex) {
   bool vertex = true;
   {
#line 391

     if (_attribute[s.i].d.x != -1)
       vertex = false;
#line 391

     if (_attribute[s.i].d.y != -1)
       vertex = false;}
   if (!vertex) {
     char * name = NULL;
     if (_attribute[s.i].name) name = pstrdup (_attribute[s.i].name,__func__,__FILE__,__LINE__);
     init_vertex_scalar (s, name);
     {
#line 398

       _attribute[s.i].v.x.i = -1;
#line 398

       _attribute[s.i].v.y.i = -1;}




     pfree (name,__func__,__FILE__,__LINE__);
   }
 }





 dirty = list_append (dirty, s);
 strongif (baseblock) for (scalar d = *baseblock, *_i30 = baseblock; ((scalar *)&d)->i >= 0; d = *++_i30)
   if (scalar_depends_from (d, s))
     dirty = list_append (dirty, d);
      }
    }
  }




  if (flux) {
#line 436 "/home/fpl/softwares/basilisk/src/grid/stencils.h"
    boundary_face (listf);
    {
#line 437

      pfree (listf.x,__func__,__FILE__,__LINE__);
#line 437

      pfree (listf.y,__func__,__FILE__,__LINE__);}
  }




  if (listc) {






    boundary_internal (listc, loop->fname, loop->line);
    pfree (listc,__func__,__FILE__,__LINE__);
  }





  if (dirty) {






    strongif (dirty) for (scalar s = *dirty, *_i31 = dirty; ((scalar *)&s)->i >= 0; s = *++_i31)
      _attribute[s.i].dirty = true;
    pfree (dirty,__func__,__FILE__,__LINE__);
  }
}
#line 41 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"



static void init_block_scalar (scalar sb, const char * name, const char * ext,
          int n, int block)
{
  char bname[strlen(name) + strlen(ext) + 10];
  if (n == 0) {
    sprintf (bname, "%s%s", name, ext);
    _attribute[sb.i].block = block;
    init_scalar (sb, bname);
    baseblock = list_append (baseblock, sb);
  }
  else {
    sprintf (bname, "%s%d%s", name, n, ext);
    _attribute[sb.i].block = - n;
    init_scalar (sb, bname);
  }
  all = list_append (all, sb);
}

scalar new_block_scalar (const char * name, const char * ext, int block)
{
  int nvar = datasize/sizeof(double);

  scalar s = {0};
  while (s.i < nvar) {
    int n = 0;
    scalar sb = s;
    while (sb.i < nvar && n < block && _attribute[sb.i].freed)
      n++, sb.i++;
    if (n >= block) {
      for (sb.i = s.i, n = 0; n < block; n++, sb.i++)
 init_block_scalar (sb, name, ext, n, block);
      trash (((scalar []){s, {-1}}));
      return s;
    }
    s.i = sb.i + 1;
  }


  s = (scalar){nvar};
  if (!(nvar + block <= _NVARMAX)) qassert ("/home/fpl/softwares/basilisk/src/grid/cartesian-common.h", 83, "nvar + block <= _NVARMAX");





  if (_attribute == NULL) {
    _attribute = (_Attributes *) pcalloc (nvar + block + 1, sizeof (_Attributes),__func__,__FILE__,__LINE__);
    _attribute[0].write = (double *) pcalloc (1, sizeof(double),__func__,__FILE__,__LINE__);
  }
  else
    _attribute = (_Attributes *)
      prealloc (_attribute - 1, (nvar + block + 1)*sizeof (_Attributes),__func__,__FILE__,__LINE__);
  _attribute++;
  memset (&_attribute[nvar], 0, block*sizeof (_Attributes));
  for (int n = 0; n < block; n++, nvar++) {
    scalar sb = (scalar){nvar};
    init_block_scalar (sb, name, ext, n, block);
  }

  realloc_scalar (block*sizeof(double));
  trash (((scalar []){s, {-1}}));
  return s;
}

scalar new_scalar (const char * name)
{
  return new_block_scalar (name, "", 1);
}

scalar new_vertex_scalar (const char * name)
{
  return init_vertex_scalar (new_scalar (name), name);
}

static vector alloc_block_vector (const char * name, int block)
{
  vector v;
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  {
#line 122

    v.x = new_block_scalar (name, ext.x, block);
#line 122

    v.y = new_block_scalar (name, ext.y, block);}
  return v;
}

vector new_vector (const char * name)
{
  vector v = alloc_block_vector (name, 1);
  init_vector (v, NULL);
  return v;
}

vector new_face_vector (const char * name)
{
  vector v = alloc_block_vector (name, 1);
  init_face_vector (v, NULL);
  return v;
}

vector new_block_vector (const char * name, int block)
{
  vector v = alloc_block_vector (name, block);
  for (int i = 0; i < block; i++) {
    vector vb;
    {
#line 146

      vb.x.i = v.x.i + i;
#line 146

      vb.y.i = v.y.i + i;}
    init_vector (vb, NULL);
    {
#line 149

      _attribute[vb.x.i].block = - i;
#line 149

      _attribute[vb.y.i].block = - i;}
  }
  {
#line 152

    _attribute[v.x.i].block = block;
#line 152

    _attribute[v.y.i].block = block;}
  return v;
}

vector new_block_face_vector (const char * name, int block)
{
  vector v = alloc_block_vector (name, block);
  for (int i = 0; i < block; i++) {
    vector vb;
    {
#line 162

      vb.x.i = v.x.i + i;
#line 162

      vb.y.i = v.y.i + i;}
    init_face_vector (vb, NULL);
    {
#line 165

      _attribute[vb.x.i].block = - i;
#line 165

      _attribute[vb.y.i].block = - i;}
  }
  {
#line 168

    _attribute[v.x.i].block = block;
#line 168

    _attribute[v.y.i].block = block;}
  return v;
}

tensor new_tensor (const char * name)
{
  char cname[strlen(name) + 3];
  struct { char * x, * y, * z; } ext = {"%s.x", "%s.y", "%s.z"};
  tensor t;
  {
#line 178
 {
    sprintf (cname, ext.x, name);
    t.x = new_vector (cname);
  }
#line 178
 {
    sprintf (cname, ext.y, name);
    t.y = new_vector (cname);
  }}
  init_tensor (t, NULL);
  return t;
}

tensor new_symmetric_tensor (const char * name)
{
  char cname[strlen(name) + 5];
  struct { char * x, * y, * z; } ext = {"%s.x.x", "%s.y.y", "%s.z.z"};
  tensor t;
  {
#line 191
 {
    sprintf (cname, ext.x, name);
    t.x.x = new_scalar(cname);
  }
#line 191
 {
    sprintf (cname, ext.y, name);
    t.y.y = new_scalar(cname);
  }}

    sprintf (cname, "%s.x.y", name);
    t.x.y = new_scalar(cname);
    t.y.x = t.x.y;
#line 211 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"
  init_tensor (t, NULL);
  return t;
}

static int nconst = 0;

void init_const_scalar (scalar s, const char * name, double val)
{
  if (s.i - _NVARMAX >= nconst) {
    nconst = s.i - _NVARMAX + 1;
    _constant = (double *) prealloc (_constant, (nconst)*sizeof(double),__func__,__FILE__,__LINE__);
  }
  _constant[s.i - _NVARMAX] = val;
}

scalar new_const_scalar (const char * name, int i, double val)
{
  scalar s = (scalar){i + _NVARMAX};
  init_const_scalar (s, name, val);
  return s;
}

void init_const_vector (vector v, const char * name, double * val)
{
  {
#line 235

    init_const_scalar (v.x, name, *val++);
#line 235

    init_const_scalar (v.y, name, *val++);}
}

vector new_const_vector (const char * name, int i, double * val)
{
  vector v;
  {
#line 242

    v.x.i = _NVARMAX + i++;
#line 242

    v.y.i = _NVARMAX + i++;}
  init_const_vector (v, name, val);
  return v;
}

void scalar_clone (scalar a, scalar b)
{
  char * name = _attribute[a.i].name;
  double * write = _attribute[a.i].write;
  int * read = _attribute[a.i].read;
  double (** boundary) (Point, Point, scalar, void *) = _attribute[a.i].boundary;
  double (** boundary_homogeneous) (Point, Point, scalar, void *) =
    _attribute[a.i].boundary_homogeneous;
  if (!(_attribute[b.i].block > 0 && _attribute[a.i].block == _attribute[b.i].block)) qassert ("/home/fpl/softwares/basilisk/src/grid/cartesian-common.h", 256, "b.block > 0 && a.block == b.block");
  pfree (_attribute[a.i].depends,__func__,__FILE__,__LINE__);
  _attribute[a.i] = _attribute[b.i];
  _attribute[a.i].name = name;
  _attribute[a.i].write = write;
  _attribute[a.i].read = read;
  _attribute[a.i].boundary = boundary;
  _attribute[a.i].boundary_homogeneous = boundary_homogeneous;
  for (int i = 0; i < nboundary; i++) {
    _attribute[a.i].boundary[i] = _attribute[b.i].boundary[i];
    _attribute[a.i].boundary_homogeneous[i] = _attribute[b.i].boundary_homogeneous[i];
  }
  _attribute[a.i].depends = list_copy (_attribute[b.i].depends);
}

scalar * list_clone (scalar * l)
{
  scalar * list = NULL;
  int nvar = datasize/sizeof(double), map[nvar];
  for (int i = 0; i < nvar; i++)
    map[i] = -1;
  strongif (l) for (scalar s = *l, *_i32 = l; ((scalar *)&s)->i >= 0; s = *++_i32) {
    scalar c = _attribute[s.i].block > 1 ? new_block_scalar("c", "", _attribute[s.i].block) :
      new_scalar("c");
    scalar_clone (c, s);
    map[s.i] = c.i;
    list = list_append (list, c);
  }
  strongif (list) for (scalar s = *list, *_i33 = list; ((scalar *)&s)->i >= 0; s = *++_i33)
    {
#line 285

      if (_attribute[s.i].v.x.i >= 0 && map[_attribute[s.i].v.x.i] >= 0)
 _attribute[s.i].v.x.i = map[_attribute[s.i].v.x.i];
#line 285

      if (_attribute[s.i].v.y.i >= 0 && map[_attribute[s.i].v.y.i] >= 0)
 _attribute[s.i].v.y.i = map[_attribute[s.i].v.y.i];}
  return list;
}

void delete (scalar * list)
{
  if (all == NULL)
    return;

  strongif (list) for (scalar f = *list, *_i34 = list; ((scalar *)&f)->i >= 0; f = *++_i34) {
    if (_attribute[f.i].block > 0) {
      pfree (_attribute[f.i].write,__func__,__FILE__,__LINE__); _attribute[f.i].write = NULL;
      pfree (_attribute[f.i].read,__func__,__FILE__,__LINE__); _attribute[f.i].read = NULL;
    }
    for (int i = 0; i < _attribute[f.i].block; i++) {
      scalar fb = {f.i + i};
      if (_attribute[f.i].delete)
 _attribute[f.i].delete (fb);
      pfree (_attribute[fb.i].name,__func__,__FILE__,__LINE__); _attribute[fb.i].name = NULL;
      _attribute[fb.i].read = NULL;
      _attribute[fb.i].write = NULL;
      pfree (_attribute[fb.i].boundary,__func__,__FILE__,__LINE__); _attribute[fb.i].boundary = NULL;
      pfree (_attribute[fb.i].boundary_homogeneous,__func__,__FILE__,__LINE__); _attribute[fb.i].boundary_homogeneous = NULL;
      pfree (_attribute[fb.i].depends,__func__,__FILE__,__LINE__); _attribute[fb.i].depends = NULL;
      _attribute[fb.i].freed = true;
    }
  }

  if (list == all) {
    all[0].i = -1;
    baseblock[0].i = -1;
    return;
  }

  trash (list);
  strongif (list) for (scalar f = *list, *_i35 = list; ((scalar *)&f)->i >= 0; f = *++_i35) {
    if (_attribute[f.i].block > 0) {
      scalar * s;
      for (s = all; s->i >= 0 && s->i != f.i; s++);
      if (s->i == f.i) {
 for (; s[_attribute[f.i].block].i >= 0; s++)
   s[0] = s[_attribute[f.i].block];
 s->i = -1;
      }
      for (s = baseblock; s->i >= 0 && s->i != f.i; s++);
      if (s->i == f.i) {
 for (; s[1].i >= 0; s++)
   s[0] = s[1];
 s->i = -1;
      }
    }
  }
}

void free_solver ()
{
  if (!(_val_higher_dimension == 0.)) qassert ("/home/fpl/softwares/basilisk/src/grid/cartesian-common.h", 343, "_val_higher_dimension == 0.");

  if (free_solver_funcs) {
    free_solver_func * a = (free_solver_func *) free_solver_funcs->p;
    for (int i = 0; i < free_solver_funcs->len/sizeof(free_solver_func); i++)
      a[i] ();
    array_free (free_solver_funcs);
  }

  delete (all);
  pfree (all,__func__,__FILE__,__LINE__); all = NULL;
  pfree (baseblock,__func__,__FILE__,__LINE__); baseblock = NULL;
  for (Event * ev = Events; !ev->last; ev++) {
    Event * e = ev->next;
    while (e) {
      Event * next = e->next;
      pfree (e,__func__,__FILE__,__LINE__);
      e = next;
    }
  }

  pfree (Events,__func__,__FILE__,__LINE__); Events = NULL;
  pfree ((_attribute - 1)->write,__func__,__FILE__,__LINE__);
  pfree (_attribute - 1,__func__,__FILE__,__LINE__); _attribute = NULL;
  pfree (_constant,__func__,__FILE__,__LINE__); _constant = NULL;
  free_grid();
  qpclose_all();
#if TRACE
  trace_off();
#endif
#if MTRACE
  pmuntrace();
#endif
#if _CADNA
  cadna_end();
#endif
}



void (* boundary_level) (scalar *, int l);
void (* boundary_face) (vectorl);




void boundary_flux (vector * list) __attribute__ ((deprecated));

void boundary_flux (vector * list)
{
  vectorl list1 = {NULL};
  strongif (list) for (vector v = *list, *_i36 = list; ((scalar *)&v)->i >= 0; v = *++_i36)
    {
#line 395

      list1.x = list_append (list1.x, v.x);
#line 395

      list1.y = list_append (list1.y, v.y);}
  boundary_face (list1);
  {
#line 398

    pfree (list1.x,__func__,__FILE__,__LINE__);
#line 398

    pfree (list1.y,__func__,__FILE__,__LINE__);}
}

static scalar * list_add_depends (scalar * list, scalar s)
{
  strongif (list) for (scalar t = *list, *_i37 = list; ((scalar *)&t)->i >= 0; t = *++_i37)
    if (t.i == s.i)
      return list;
  scalar * list1 = list;
  strongif (_attribute[s.i].depends) for (scalar d = *_attribute[s.i].depends, *_i38 = _attribute[s.i].depends; ((scalar *)&d)->i >= 0; d = *++_i38)
    if (_attribute[d.i].dirty)
      list1 = list_add_depends (list1, d);
  return list_append (list1, s);
}


void boundary_internal (scalar * list, const char * fname, int line)
{ trace ("boundary_internal", "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h", 416);
  if (list == NULL)
    { ; end_trace("boundary_internal", "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h", 418);  return; }
  scalar * listc = NULL;
  vectorl listf = {NULL};
  bool flux = false;
  strongif (list) for (scalar s = *list, *_i39 = list; ((scalar *)&s)->i >= 0; s = *++_i39)
    if (!is_constant(s) && _attribute[s.i].block > 0) {
      if (scalar_is_dirty (s)) {
 if (_attribute[s.i].face && _attribute[s.i].dirty != 2)
   {
#line 426

     if (_attribute[s.i].v.x.i == s.i)
       listf.x = list_add (listf.x, s), flux = true;
#line 426

     if (_attribute[s.i].v.y.i == s.i)
       listf.y = list_add (listf.y, s), flux = true;}
 if (!is_constant(cm) && _attribute[cm.i].dirty)
   listc = list_add_depends (listc, cm);
 if (_attribute[s.i].face != 2)
   listc = list_add_depends (listc, s);
      }




    }
  if (flux) {
    boundary_face (listf);
    {
#line 441

      pfree (listf.x,__func__,__FILE__,__LINE__);
#line 441

      pfree (listf.y,__func__,__FILE__,__LINE__);}
  }
  if (listc) {
    boundary_level (listc, -1);
    strongif (listc) for (scalar s = *listc, *_i40 = listc; ((scalar *)&s)->i >= 0; s = *++_i40)
      _attribute[s.i].dirty = false;
    pfree (listc,__func__,__FILE__,__LINE__);
  }
 end_trace("boundary_internal", "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h", 450); }

void cartesian_boundary_level (scalar * list, int l)
{
  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, l); };
}

void cartesian_boundary_face (vectorl list)
{
  {
#line 459

    strongif (list.x) for (scalar s = *list.x, *_i41 = list.x; ((scalar *)&s)->i >= 0; s = *++_i41)
      _attribute[s.i].dirty = 2;
#line 459

    strongif (list.y) for (scalar s = *list.y, *_i41 = list.y; ((scalar *)&s)->i >= 0; s = *++_i41)
      _attribute[s.i].dirty = 2;}
}

static double symmetry (Point point, Point neighbor, scalar s, void * data)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 465 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"

  return val(s,0,0,0);

#if _call_symmetry
}
#define _IN_STENCIL 1

#line 464
static double _symmetry (Point point, Point neighbor, scalar s, void * data)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 465 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"

  return _stencil_val(__FILE__,__LINE__,s,0,0,0);

#undef _IN_STENCIL

#endif

#line 467
}

static double antisymmetry (Point point, Point neighbor, scalar s, void * data)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 470 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"

  return -val(s,0,0,0);

#if _call_antisymmetry
}
#define _IN_STENCIL 1

#line 469
static double _antisymmetry (Point point, Point neighbor, scalar s, void * data)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 470 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"

  return -_stencil_val(__FILE__,__LINE__,s,0,0,0);

#undef _IN_STENCIL

#endif

#line 472
}

double (* default_scalar_bc[]) (Point, Point, scalar, void *) = {
  symmetry, symmetry, symmetry, symmetry, symmetry, symmetry
};

static double centered (double s0, double s1, double s2) {
  return (s2 - s0)/2.;
}

scalar cartesian_init_scalar (scalar s, const char * name)
{

  char * pname;
  if (name) {
    pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
    pname = pstrdup (name,__func__,__FILE__,__LINE__);
  }
  else
    pname = _attribute[s.i].name;
  int block = _attribute[s.i].block;
  double * write = _attribute[s.i].write;
  int * read = _attribute[s.i].read;
  double (** boundary) (Point, Point, scalar, void *) = _attribute[s.i].boundary;
  double (** boundary_homogeneous) (Point, Point, scalar, void *) =
    _attribute[s.i].boundary_homogeneous;

  _attribute[s.i] = (const _Attributes){0};
  _attribute[s.i].name = pname;
  if (block < 0) {
    scalar base = {s.i + block};
    _attribute[s.i].block = block;
    _attribute[s.i].write = _attribute[base.i].write;
    _attribute[s.i].read = _attribute[base.i].read;
  }
  else {
    _attribute[s.i].block = block > 0 ? block : 1;
    _attribute[s.i].write = write ? write : pmalloc((5*5)*sizeof(double),__func__,__FILE__,__LINE__);
    _attribute[s.i].read = read ? read : pmalloc((5*5)*sizeof(int),__func__,__FILE__,__LINE__);
  }

  _attribute[s.i].boundary = boundary ? boundary :
    (double (**)(Point, Point, scalar, void *))
    pmalloc (nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  _attribute[s.i].boundary_homogeneous = boundary_homogeneous ? boundary_homogeneous :
    (double (**)(Point, Point, scalar, void *))
    pmalloc (nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  for (int b = 0; b < nboundary; b++)
    _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] =
      b < 2*2 ? default_scalar_bc[b] : symmetry;
  _attribute[s.i].gradient = centered;
  {
#line 523
 {
    _attribute[s.i].d.x = 0;
    _attribute[s.i].v.x.i = -1;
  }
#line 523
 {
    _attribute[s.i].d.y = 0;
    _attribute[s.i].v.y.i = -1;
  }}
  _attribute[s.i].face = false;
  return s;
}

scalar cartesian_init_vertex_scalar (scalar s, const char * name)
{
  {
#line 533

    _attribute[s.i].d.x = -1;
#line 533

    _attribute[s.i].d.y = -1;}
  for (int d = 0; d < nboundary; d++)
    _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = NULL;
  return s;
}

double (* default_vector_bc[]) (Point, Point, scalar, void *) = {
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry
};

vector cartesian_init_vector (vector v, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  {
#line 549
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.x);
      init_scalar (v.x, cname);
    }
    else
      init_scalar (v.x, NULL);
    _attribute[v.x.i].v = v;
  }
#line 549
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.y);
      init_scalar (v.y, cname);
    }
    else
      init_scalar (v.y, NULL);
    _attribute[v.y.i].v = v;
  }}

  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] =
      d < 2*2 ? default_vector_bc[d] : antisymmetry;
  return v;
}

vector cartesian_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_vector (v, name);
  {
#line 569
 {
    _attribute[v.x.i].d.x = -1;
    _attribute[v.x.i].face = true;
  }
#line 569
 {
    _attribute[v.y.i].d.y = -1;
    _attribute[v.y.i].face = true;
  }}
  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
  return v;
}

tensor cartesian_init_tensor (tensor t, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  {
#line 581
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.x);
      init_vector (t.x, cname);
    }
    else
      init_vector (t.x, NULL);
  }
#line 581
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.y);
      init_vector (t.y, cname);
    }
    else
      init_vector (t.y, NULL);
  }}






    for (int b = 0; b < nboundary; b++) {
      _attribute[t.x.x.i].boundary[b] = _attribute[t.y.x.i].boundary[b] =
 _attribute[t.x.x.i].boundary_homogeneous[b] = _attribute[t.y.y.i].boundary_homogeneous[b] =
 b < 2*2 ? default_scalar_bc[b] : symmetry;
      _attribute[t.x.y.i].boundary[b] = _attribute[t.y.y.i].boundary[b] =
 _attribute[t.x.y.i].boundary_homogeneous[b] = _attribute[t.y.x.i].boundary_homogeneous[b] =
 b < 2*2 ? default_vector_bc[b] : antisymmetry;
    }



  return t;
}

struct OutputCells {
  FILE * fp;
  coord c;
  double size;
};

void output_cells (struct OutputCells p)
{
  if (!p.fp) p.fp = fout;
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h", .line = 619,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 619 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"
 {
    bool inside = true;
    coord o = {x,y,z};
    {
#line 622

      IF (inside && p.size > 0. &&
   (o.x > p.c.x + p.size || o.x < p.c.x - p.size))
 inside = false;
#line 622

      IF (inside && p.size > 0. &&
   (o.y > p.c.y + p.size || o.y < p.c.y - p.size))
 inside = false;}
    IF (inside) {
      Delta /= 2.;



      _stencil_fprintf (__FILE__,__LINE__,p.fp, "%g %g\n%g %g\n%g %g\n%g %g\n%g %g\n\n",
        x - Delta, y - Delta,
        x - Delta, y + Delta,
        x + Delta, y + Delta,
        x + Delta, y - Delta,
        x - Delta, y - Delta);
#line 651 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"
    }
  } } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 652
foreach(){

#line 619 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"
 {
    bool inside = true;
    coord o = {x,y,z};
    {
#line 622

      if (inside && p.size > 0. &&
   (o.x > p.c.x + p.size || o.x < p.c.x - p.size))
 inside = false;
#line 622

      if (inside && p.size > 0. &&
   (o.y > p.c.y + p.size || o.y < p.c.y - p.size))
 inside = false;}
    if (inside) {
      Delta /= 2.;



      fprintf (p.fp, "%g %g\n%g %g\n%g %g\n%g %g\n%g %g\n\n",
        x - Delta, y - Delta,
        x - Delta, y + Delta,
        x + Delta, y + Delta,
        x + Delta, y - Delta,
        x - Delta, y - Delta);
#line 651 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"
    }
  } } end_foreach(); }
  fflush (p.fp);
}


static void output_cells_internal (FILE * fp)
{
  output_cells ((struct OutputCells){fp});
}


static char * replace_ (const char * vname)
{
  char * name = pstrdup (vname,__func__,__FILE__,__LINE__), * c = name;
  while (*c != '\0') {
    if (*c == '.')
      *c = '_';
    c++;
  }
  return name;
}

static void debug_plot (FILE * fp, const char * name, const char * cells,
   const char * stencil)
{
  char * vname = replace_ (name);
  fprintf (fp,
    "  load 'debug.plot'\n"
    "  v=%s\n"




    "  plot '%s' w l lc 0, "
    "'%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 1 title columnhead(3+3*v)",





    vname, cells, stencil);
  pfree (vname,__func__,__FILE__,__LINE__);
}

void cartesian_debug (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 697 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"

  char name[80] = "cells";
  if (pid() > 0)
    sprintf (name, "cells-%d", pid());
  FILE * fp = fopen (name, "w");
  output_cells ((struct OutputCells){fp, (coord){x,y,z}, 4.*Delta});
  fclose (fp);

  char stencil[80] = "stencil";
  if (pid() > 0)
    sprintf (stencil, "stencil-%d", pid());
  fp = fopen (stencil, "w");
  strongif (all) for (scalar v = *all, *_i42 = all; ((scalar *)&v)->i >= 0; v = *++_i42)



    fprintf (fp, "x y %s ", _attribute[v.i].name);



  fputc ('\n', fp);
#line 730 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"
    for (int k = -2; k <= 2; k++)
      for (int l = -2; l <= 2; l++) {
 strongif (all) for (scalar v = *all, *_i43 = all; ((scalar *)&v)->i >= 0; v = *++_i43) {
   fprintf (fp, "%g %g ",
     x + k*Delta + _attribute[v.i].d.x*Delta/2.,
     y + l*Delta + _attribute[v.i].d.y*Delta/2.);
   if (allocated(k,l,0))
     fprintf (fp, "%g ", val(v,k,l,0));
   else
     fputs ("n/a ", fp);
 }
 fputc ('\n', fp);
      }
#line 760 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"
  fclose (fp);

  fp = fopen ("debug.plot", "w");
  fprintf (fp,
    "set term x11\n"
    "set size ratio -1\n"
    "set key outside\n");
  strongif (all) for (scalar s = *all, *_i44 = all; ((scalar *)&s)->i >= 0; s = *++_i44) {
    char * name = replace_ (_attribute[s.i].name);
    fprintf (fp, "%s = %d\n", name, s.i);
    pfree (name,__func__,__FILE__,__LINE__);
  }
  fclose (fp);

  fprintf (ferr, "Last point stencils can be displayed using (in gnuplot)\n");
  debug_plot (ferr, _attribute[0].name, name, stencil);
  fflush (ferr);

  fp = fopen ("plot", "w");
  debug_plot (fp, _attribute[0].name, name, stencil);
  fclose (fp);

#if _call_cartesian_debug
}
#define _IN_STENCIL 1

#line 696
static void _cartesian_debug (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 697 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"

  char name[80] = "cells";
  IF (pid() > 0)
    sprintf (name, "cells-%d", pid());
  FILE * fp = fopen (name, "w");
  output_cells ((struct OutputCells){fp, (coord){x,y,z}, 4.*Delta});
  fclose (fp);

  char stencil[80] = "stencil";
  IF (pid() > 0)
    sprintf (stencil, "stencil-%d", pid());
  fp = fopen (stencil, "w");
  strongif (all) for (scalar v = *all, *_i42 = all; ((scalar *)&v)->i >= 0; v = *++_i42)



    _stencil_fprintf (__FILE__,__LINE__,fp, "x y %s ", _attribute[v.i].name);



  _stencil_fputc (__FILE__,__LINE__,'\n', fp);
#line 730 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"
    for (int k = -2; k <= 2; k++)
      for (int l = -2; l <= 2; l++) {
 strongif (all) for (scalar v = *all, *_i43 = all; ((scalar *)&v)->i >= 0; v = *++_i43) {
   _stencil_fprintf (__FILE__,__LINE__,fp, "%g %g ",
     x + k*Delta + _attribute[v.i].d.x*Delta/2.,
     y + l*Delta + _attribute[v.i].d.y*Delta/2.);
   IF (allocated(k,l,0))
     _stencil_fprintf (__FILE__,__LINE__,fp, "%g ", _stencil_val(__FILE__,__LINE__,v,k,l,0));
   
     _stencil_fputs (__FILE__,__LINE__,"n/a ", fp);
 }
 _stencil_fputc (__FILE__,__LINE__,'\n', fp);
      }
#line 760 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"
  fclose (fp);

  fp = fopen ("debug.plot", "w");
  _stencil_fprintf (__FILE__,__LINE__,fp,
    "set term x11\n"
    "set size ratio -1\n"
    "set key outside\n");
  strongif (all) for (scalar s = *all, *_i44 = all; ((scalar *)&s)->i >= 0; s = *++_i44) {
    char * name = replace_ (_attribute[s.i].name);
    _stencil_fprintf (__FILE__,__LINE__,fp, "%s = %d\n", name, s.i);
    pfree (name,__func__,__FILE__,__LINE__);
  }
  fclose (fp);

  _stencil_fprintf (__FILE__,__LINE__,ferr, "Last point stencils can be displayed using (in gnuplot)\n");
  debug_plot (ferr, _attribute[0].name, name, stencil);
  fflush (ferr);

  fp = fopen ("plot", "w");
  debug_plot (fp, _attribute[0].name, name, stencil);
  fclose (fp);

#undef _IN_STENCIL

#endif

#line 781
}

void cartesian_methods ()
{
  init_scalar = cartesian_init_scalar;
  init_vertex_scalar = cartesian_init_vertex_scalar;
  init_vector = cartesian_init_vector;
  init_tensor = cartesian_init_tensor;
  init_face_vector = cartesian_init_face_vector;
  boundary_level = cartesian_boundary_level;
  boundary_face = cartesian_boundary_face;
  debug = cartesian_debug;
}

struct _interpolate {
  scalar v;
  double x, y, z;
};

static double interpolate_linear (Point point, struct _interpolate p)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 801 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"

  scalar v = p.v;







  x = (p.x - x)/Delta - _attribute[v.i].d.x/2.;
  y = (p.y - y)/Delta - _attribute[v.i].d.y/2.;
  int i = sign(x), j = sign(y);
  x = fabs(x); y = fabs(y);

  return ((val(v,0,0,0)*(1. - x) + val(v,i,0,0)*x)*(1. - y) +
   (val(v,0,j,0)*(1. - x) + val(v,i,j,0)*x)*y);
#line 829 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"

#if _call_interpolate_linear
}
#define _IN_STENCIL 1

#line 800
static double _interpolate_linear (Point point, struct _interpolate p)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 801 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"

  scalar v = p.v;







  x = (p.x - x)/Delta - _attribute[v.i].d.x/2.;
  y = (p.y - y)/Delta - _attribute[v.i].d.y/2.;
  int i = sign(x), j = sign(y);
  x = fabs(x); y = fabs(y);

  return ((_stencil_val(__FILE__,__LINE__,v,0,0,0)*(1. - x) + _stencil_val(__FILE__,__LINE__,v,i,0,0)*x)*(1. - y) +
   (_stencil_val(__FILE__,__LINE__,v,0,j,0)*(1. - x) + _stencil_val(__FILE__,__LINE__,v,i,j,0)*x)*y);
#line 829 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"

#undef _IN_STENCIL

#endif

#line 829
}


double interpolate (struct _interpolate p)
{ trace ("interpolate", "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h", 833);
  scalar v = p.v;
  boundary_internal ((scalar *)(((scalar []){v,{-1}})), "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h", 835);
  Point point = locate ((struct _locate){p.x, p.y, p.z});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 836 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"

  if (point.level < 0)
    { double _ret =  nodata; end_trace("interpolate", "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h", 838);  return _ret; }
  { double _ret =  interpolate_linear (point, p); end_trace("interpolate", "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h", 839);  return _ret; }
 end_trace("interpolate", "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h", 840); }


void interpolate_array (scalar * list, coord * a, int n, double * v, bool linear)
{ trace ("interpolate_array", "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h", 844);
  boundary_internal ((scalar *)(list), "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h", 845);
  int j = 0;
  for (int i = 0; i < n; i++) {
    Point point = locate ((struct _locate){a[i].x, a[i].y, a[i].z});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 848 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"

    if (point.level >= 0) {
      strongif (list) for (scalar s = *list, *_i45 = list; ((scalar *)&s)->i >= 0; s = *++_i45)
 v[j++] = !linear ? val(s,0,0,0) :
   interpolate_linear (point,
         (struct _interpolate){s, a[i].x, a[i].y, a[i].z});
    }
    else
      strongif (list) for (scalar s = *list, *_i46 = list; ((scalar *)&s)->i >= 0; s = *++_i46)
 v[j++] = nodata;
  }
#if 1
  if (pid() == 0)
    MPI_Reduce (MPI_IN_PLACE, v, n*list_len(list), MPI_DOUBLE,
  MPI_MIN, 0, MPI_COMM_WORLD);
  else
    MPI_Reduce (v, v, n*list_len(list), MPI_DOUBLE,
  MPI_MIN, 0, MPI_COMM_WORLD);
#endif
 end_trace("interpolate_array", "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h", 867); }



typedef int bid;

bid new_bid ()
{
  int b = nboundary++;
  strongif (all) for (scalar s = *all, *_i47 = all; ((scalar *)&s)->i >= 0; s = *++_i47) {
    _attribute[s.i].boundary = (double (**)(Point, Point, scalar, void *))
      prealloc (_attribute[s.i].boundary, nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
    _attribute[s.i].boundary_homogeneous = (double (**)(Point, Point, scalar, void *))
      prealloc (_attribute[s.i].boundary_homogeneous, nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  }
  strongif (all) for (scalar s = *all, *_i48 = all; ((scalar *)&s)->i >= 0; s = *++_i48) {
    if (_attribute[s.i].v.x.i < 0)
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] = symmetry;
    else if (_attribute[s.i].v.x.i == s.i) {
      vector v = _attribute[s.i].v;
      {
#line 887

 _attribute[v.y.i].boundary[b] = _attribute[v.y.i].boundary_homogeneous[b] = symmetry;
#line 887

 _attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] = symmetry;}
      _attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] =
 _attribute[v.x.i].face ? NULL : antisymmetry;
    }
  }
  return b;
}



static double periodic_bc (Point point, Point neighbor, scalar s, void * data)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 899 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"

  return nodata;

#if _call_periodic_bc
}
#define _IN_STENCIL 1

#line 898
static double _periodic_bc (Point point, Point neighbor, scalar s, void * data)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 899 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"

  return nodata;

#undef _IN_STENCIL

#endif

#line 901
}

static void periodic_boundary (int d)
{

  strongif (all) for (scalar s = *all, *_i49 = all; ((scalar *)&s)->i >= 0; s = *++_i49)
    _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = periodic_bc;

  strongif (all) for (scalar s = *all, *_i50 = all; ((scalar *)&s)->i >= 0; s = *++_i50)
    if (_attribute[s.i].face) {
      vector v = _attribute[s.i].v;
      _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
    }

  default_scalar_bc[d] = periodic_bc;
  default_vector_bc[d] = periodic_bc;
}

void periodic (int dir)
{



    if (!(dir <= bottom)) qassert ("/home/fpl/softwares/basilisk/src/grid/cartesian-common.h", 924, "dir <= bottom");




  int c = dir/2;
  periodic_boundary (2*c);
  periodic_boundary (2*c + 1);
  (&Period.x)[c] = true;
}


double getvalue (Point point, scalar s, int i, int j, int k)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 937 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"

  return val(s,i,j,k);

#if _call_getvalue
}
#define _IN_STENCIL 1

#line 936
static double _getvalue (Point point, scalar s, int i, int j, int k)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 937 "/home/fpl/softwares/basilisk/src/grid/cartesian-common.h"

  return _stencil_val(__FILE__,__LINE__,s,i,j,k);

#undef _IN_STENCIL

#endif

#line 939
}
#line 4 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

#ifndef foreach_level_or_leaf
# define foreach_level_or_leaf foreach_level
# define end_foreach_level_or_leaf end_foreach_level
#endif

#ifndef foreach_coarse_level
# define foreach_coarse_level foreach_level
# define end_foreach_coarse_level end_foreach_level
#endif










void (* restriction) (scalar *);

static inline void restriction_average (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 27 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

  double sum = 0.;
   { foreach_child()
    sum += val(s,0,0,0); end_foreach_child(); }
  val(s,0,0,0) = sum/(1 << 2);

#if _call_restriction_average
}
#define _IN_STENCIL 1

#line 26
static void _restriction_average (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 27 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

  double sum = 0.;
   { foreach_child()
    sum += _stencil_val(__FILE__,__LINE__,s,0,0,0); end_foreach_child(); }
  _stencil_val(__FILE__,__LINE__,s,0,0,0) = sum/(1 << 2);

#undef _IN_STENCIL

#endif

#line 32
}

static inline void restriction_volume_average (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 35 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 35

  double sum = 0.;
   { foreach_child()
    sum += val_cm(cm,0,0,0)*val(s,0,0,0); end_foreach_child(); }
  val(s,0,0,0) = sum/(1 << 2)/(val_cm(cm,0,0,0) + 1e-30);
 }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 35

  double sum = 0.;
   { foreach_child()
    sum += val_cm(cm,0,0,0)*val(s,0,0,0); end_foreach_child(); }
  val(s,0,0,0) = sum/(1 << 2)/(val_cm(cm,0,0,0) + 1e-30);
 }
#if _call_restriction_volume_average
}
#define _IN_STENCIL 1

#line 34
static void _restriction_volume_average (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 35 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 35

  double sum = 0.;
   { foreach_child()
    sum += val_cm(cm,0,0,0)*_stencil_val(__FILE__,__LINE__,s,0,0,0); end_foreach_child(); }
  _stencil_val(__FILE__,__LINE__,s,0,0,0) = sum/(1 << 2)/(val_cm(cm,0,0,0) + 1e-30);
 }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 35

  double sum = 0.;
   { foreach_child()
    sum += val_cm(cm,0,0,0)*_stencil_val(__FILE__,__LINE__,s,0,0,0); end_foreach_child(); }
  _stencil_val(__FILE__,__LINE__,s,0,0,0) = sum/(1 << 2)/(val_cm(cm,0,0,0) + 1e-30);
 }
#undef _IN_STENCIL

#endif

#line 40
}

static inline void face_average (Point point, vector v)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 43 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

  {
#line 44
 {




      val(v.x,0,0,0) = (fine(v.x,0,0,0) + fine(v.x,0,1,0))/2.;
      val(v.x,1,0,0) = (fine(v.x,2,0,0) + fine(v.x,2,1,0))/2.;






  }
#line 44
 {




      val(v.y,0,0,0) = (fine(v.y,0,0,0) + fine(v.y,1,0,0))/2.;
      val(v.y,0,1,0) = (fine(v.y,0,2,0) + fine(v.y,1,2,0))/2.;






  }}

#if _call_face_average
}
#define _IN_STENCIL 1

#line 42
static void _face_average (Point point, vector v)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 43 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

  {
#line 44
 {




      _stencil_val(__FILE__,__LINE__,v.x,0,0,0) = (_stencil_fine(__FILE__,__LINE__,v.x,0,0,0) + _stencil_fine(__FILE__,__LINE__,v.x,0,1,0))/2.;
      _stencil_val(__FILE__,__LINE__,v.x,1,0,0) = (_stencil_fine(__FILE__,__LINE__,v.x,2,0,0) + _stencil_fine(__FILE__,__LINE__,v.x,2,1,0))/2.;






  }
#line 44
 {




      _stencil_val(__FILE__,__LINE__,v.y,0,0,0) = (_stencil_fine(__FILE__,__LINE__,v.y,0,0,0) + _stencil_fine(__FILE__,__LINE__,v.y,1,0,0))/2.;
      _stencil_val(__FILE__,__LINE__,v.y,0,1,0) = (_stencil_fine(__FILE__,__LINE__,v.y,0,2,0) + _stencil_fine(__FILE__,__LINE__,v.y,1,2,0))/2.;






  }}

#undef _IN_STENCIL

#endif

#line 58
}

static inline void restriction_face (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 61 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

  face_average (point, _attribute[s.i].v);

#if _call_restriction_face
}
#define _IN_STENCIL 1
#define face_average _face_average

#line 60
static void _restriction_face (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 61 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

  face_average (point, _attribute[s.i].v);

#undef face_average
#undef _IN_STENCIL

#endif

#line 63
}

static inline void restriction_vertex (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 66 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

  for (int i = 0; i <= 1; i++) {
    val(s,i,0,0) = fine(s,2*i,0,0);

    val(s,i,1,0) = fine(s,2*i,2,0);





  }

#if _call_restriction_vertex
}
#define _IN_STENCIL 1

#line 65
static void _restriction_vertex (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 66 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

  for (int i = 0; i <= 1; i++) {
    _stencil_val(__FILE__,__LINE__,s,i,0,0) = _stencil_fine(__FILE__,__LINE__,s,2*i,0,0);

    _stencil_val(__FILE__,__LINE__,s,i,1,0) = _stencil_fine(__FILE__,__LINE__,s,2*i,2,0);





  }

#undef _IN_STENCIL

#endif

#line 77
}

static inline void no_restriction (Point point, scalar s) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 78 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

#if _call_no_restriction
}
#define _IN_STENCIL 1

#line 79
static void _no_restriction (Point point, scalar s) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 78 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

#undef _IN_STENCIL

#endif

#line 79
}

static inline void no_data (Point point, scalar s) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 81 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    val(s,0,0,0) = nodata; end_foreach_child(); }

#if _call_no_data
}
#define _IN_STENCIL 1

#line 81
static void _no_data (Point point, scalar s) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 81 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    _stencil_val(__FILE__,__LINE__,s,0,0,0) = nodata; end_foreach_child(); }

#undef _IN_STENCIL

#endif

#line 84
}

void wavelet (scalar s, scalar w)
{
  restriction (((scalar []){s,{-1}}));
  for (int l = depth() - 1; l >= 0; l--) {
     { foreach_coarse_level (l){

#line 90 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"
 {
       { foreach_child()
        val(w,0,0,0) = val(s,0,0,0); end_foreach_child(); }
      _attribute[s.i].prolongation (point, s);
       { foreach_child() {
        double sp = val(s,0,0,0);
        val(s,0,0,0) = val(w,0,0,0);

        val(w,0,0,0) -= sp;
      } end_foreach_child(); }
    } } end_foreach_coarse_level(); }
    boundary_level (((scalar []){w,{-1}}), l + 1);
  }

   { foreach_level(0){

#line 104 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

    val(w,0,0,0) = val(s,0,0,0); } end_foreach_level(); }
  boundary_level (((scalar []){w,{-1}}), 0);
}

void inverse_wavelet (scalar s, scalar w)
{
   { foreach_level(0){

#line 111 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

    val(s,0,0,0) = val(w,0,0,0); } end_foreach_level(); }
  boundary_level (((scalar []){s,{-1}}), 0);
  for (int l = 0; l <= depth() - 1; l++) {
     { foreach_coarse_level (l){

#line 115 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"
 {
      _attribute[s.i].prolongation (point, s);
       { foreach_child()
        val(s,0,0,0) += val(w,0,0,0); end_foreach_child(); }
    } } end_foreach_coarse_level(); }
    boundary_level (((scalar []){s,{-1}}), l + 1);
  }
}

static inline double bilinear (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 125 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"




    return (9.*coarse(s,0,0,0) +
     3.*(coarse(s,child.x,0,0) + coarse(s,0,child.y,0)) +
     coarse(s,child.x,child.y,0))/16.;
#line 140 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

#if _call_bilinear
}
#define _IN_STENCIL 1

#line 124
static double _bilinear (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 125 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"




    return (9.*_stencil_coarse(__FILE__,__LINE__,s,0,0,0) +
     3.*(_stencil_coarse(__FILE__,__LINE__,s,child.x,0,0) + _stencil_coarse(__FILE__,__LINE__,s,0,child.y,0)) +
     _stencil_coarse(__FILE__,__LINE__,s,child.x,child.y,0))/16.;
#line 140 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

#undef _IN_STENCIL

#endif

#line 140
}

static inline void refine_bilinear (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 143 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    val(s,0,0,0) = bilinear (point, s); end_foreach_child(); }

#if _call_refine_bilinear
}
#define _IN_STENCIL 1
#define bilinear _bilinear

#line 142
static void _refine_bilinear (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 143 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    _stencil_val(__FILE__,__LINE__,s,0,0,0) = bilinear (point, s); end_foreach_child(); }

#undef bilinear
#undef _IN_STENCIL

#endif

#line 146
}

static inline double quadratic (double a, double b, double c)
{
  return (30.*a + 5.*b - 3.*c)/32.;
}

static inline double biquadratic (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 154 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"




  return
    quadratic (quadratic (coarse(s,0,0,0),
     coarse(s,child.x,0,0),
     coarse(s,-child.x,0,0)),
        quadratic (coarse(s,0,child.y,0),
     coarse(s,child.x,child.y,0),
     coarse(s,-child.x,child.y,0)),
        quadratic (coarse(s,0,-child.y,0),
     coarse(s,child.x,-child.y,0),
     coarse(s,-child.x,-child.y,0)));





#if _call_biquadratic
}
#define _IN_STENCIL 1

#line 153
static double _biquadratic (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 154 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"




  return
    quadratic (quadratic (_stencil_coarse(__FILE__,__LINE__,s,0,0,0),
     _stencil_coarse(__FILE__,__LINE__,s,child.x,0,0),
     _stencil_coarse(__FILE__,__LINE__,s,-child.x,0,0)),
        quadratic (_stencil_coarse(__FILE__,__LINE__,s,0,child.y,0),
     _stencil_coarse(__FILE__,__LINE__,s,child.x,child.y,0),
     _stencil_coarse(__FILE__,__LINE__,s,-child.x,child.y,0)),
        quadratic (_stencil_coarse(__FILE__,__LINE__,s,0,-child.y,0),
     _stencil_coarse(__FILE__,__LINE__,s,child.x,-child.y,0),
     _stencil_coarse(__FILE__,__LINE__,s,-child.x,-child.y,0)));





#undef _IN_STENCIL

#endif

#line 172
}

static inline double biquadratic_vertex (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 175 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"




  return (36.*val(s,0,0,0) + 18.*(val(s,-1,0,0) + val(s,0,-1,0)) - 6.*(val(s,1,0,0) + val(s,0,1,0)) +
   9.*val(s,-1,-1,0) - 3.*(val(s,1,-1,0) + val(s,-1,1,0)) + val(s,1,1,0))/64.;





#if _call_biquadratic_vertex
}
#define _IN_STENCIL 1

#line 174
static double _biquadratic_vertex (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 175 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"




  return (36.*_stencil_val(__FILE__,__LINE__,s,0,0,0) + 18.*(_stencil_val(__FILE__,__LINE__,s,-1,0,0) + _stencil_val(__FILE__,__LINE__,s,0,-1,0)) - 6.*(_stencil_val(__FILE__,__LINE__,s,1,0,0) + _stencil_val(__FILE__,__LINE__,s,0,1,0)) +
   9.*_stencil_val(__FILE__,__LINE__,s,-1,-1,0) - 3.*(_stencil_val(__FILE__,__LINE__,s,1,-1,0) + _stencil_val(__FILE__,__LINE__,s,-1,1,0)) + _stencil_val(__FILE__,__LINE__,s,1,1,0))/64.;





#undef _IN_STENCIL

#endif

#line 185
}

static inline void refine_biquadratic (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 188 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    val(s,0,0,0) = biquadratic (point, s); end_foreach_child(); }

#if _call_refine_biquadratic
}
#define _IN_STENCIL 1
#define biquadratic _biquadratic

#line 187
static void _refine_biquadratic (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 188 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    _stencil_val(__FILE__,__LINE__,s,0,0,0) = biquadratic (point, s); end_foreach_child(); }

#undef biquadratic
#undef _IN_STENCIL

#endif

#line 191
}

static inline void refine_linear (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 194 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 194

  coord g;
  if (_attribute[s.i].gradient)
    {
#line 197

      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));
#line 197

      g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));}
  else
    {
#line 200

      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;
#line 200

      g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*val_cm(cm,0,0,0), sum = val_cm(cm,0,0,0)*(1 << 2);
   { foreach_child() {
    val(s,0,0,0) = sc;
    {
#line 206

      val(s,0,0,0) += child.x*g.x*val_cm(cm,-child.x,0,0)/cmc;
#line 206

      val(s,0,0,0) += child.y*g.y*val_cm(cm,0,-child.y,0)/cmc;}
    sum -= val_cm(cm,0,0,0);
  } end_foreach_child(); }
  if (!(fabs(sum) < 1e-10)) qassert ("/home/fpl/softwares/basilisk/src/grid/multigrid-common.h", 210, "fabs(sum) < 1e-10");
 }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 194

  coord g;
  if (_attribute[s.i].gradient)
    {
#line 197

      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));
#line 197

      g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));}
  else
    {
#line 200

      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;
#line 200

      g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*val_cm(cm,0,0,0), sum = val_cm(cm,0,0,0)*(1 << 2);
   { foreach_child() {
    val(s,0,0,0) = sc;
    {
#line 206

      val(s,0,0,0) += child.x*g.x*val_cm(cm,-child.x,0,0)/cmc;
#line 206

      val(s,0,0,0) += child.y*g.y*val_cm(cm,0,-child.y,0)/cmc;}
    sum -= val_cm(cm,0,0,0);
  } end_foreach_child(); }
  if (!(fabs(sum) < 1e-10)) qassert ("/home/fpl/softwares/basilisk/src/grid/multigrid-common.h", 210, "fabs(sum) < 1e-10");
 }
#if _call_refine_linear
}
#define _IN_STENCIL 1

#line 193
static void _refine_linear (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 194 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 194

  coord g;
  IF (_attribute[s.i].gradient)
    {
#line 197

      g.x = _attribute[s.i].gradient (_stencil_val(__FILE__,__LINE__,s,-1,0,0), _stencil_val(__FILE__,__LINE__,s,0,0,0), _stencil_val(__FILE__,__LINE__,s,1,0,0));
#line 197

      g.y = _attribute[s.i].gradient (_stencil_val(__FILE__,__LINE__,s,0,-1,0), _stencil_val(__FILE__,__LINE__,s,0,0,0), _stencil_val(__FILE__,__LINE__,s,0,1,0));}
  
    {
#line 200

      g.x = (_stencil_val(__FILE__,__LINE__,s,1,0,0) - _stencil_val(__FILE__,__LINE__,s,-1,0,0))/2.;
#line 200

      g.y = (_stencil_val(__FILE__,__LINE__,s,0,1,0) - _stencil_val(__FILE__,__LINE__,s,0,-1,0))/2.;}

  double sc = _stencil_val(__FILE__,__LINE__,s,0,0,0), cmc = 4.*val_cm(cm,0,0,0), sum = val_cm(cm,0,0,0)*(1 << 2);
   { foreach_child() {
    _stencil_val(__FILE__,__LINE__,s,0,0,0) = sc;
    {
#line 206

      _stencil_val(__FILE__,__LINE__,s,0,0,0) += child.x*g.x*val_cm(cm,-child.x,0,0)/cmc;
#line 206

      _stencil_val(__FILE__,__LINE__,s,0,0,0) += child.y*g.y*val_cm(cm,0,-child.y,0)/cmc;}
    sum -= val_cm(cm,0,0,0);
  } end_foreach_child(); }
  IF (!(fabs(sum) < 1e-10)) _stencil_qassert (__FILE__,__LINE__,"/home/fpl/softwares/basilisk/src/grid/multigrid-common.h", 210, "fabs(sum) < 1e-10");
 }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 194

  coord g;
  IF (_attribute[s.i].gradient)
    {
#line 197

      g.x = _attribute[s.i].gradient (_stencil_val(__FILE__,__LINE__,s,-1,0,0), _stencil_val(__FILE__,__LINE__,s,0,0,0), _stencil_val(__FILE__,__LINE__,s,1,0,0));
#line 197

      g.y = _attribute[s.i].gradient (_stencil_val(__FILE__,__LINE__,s,0,-1,0), _stencil_val(__FILE__,__LINE__,s,0,0,0), _stencil_val(__FILE__,__LINE__,s,0,1,0));}
  
    {
#line 200

      g.x = (_stencil_val(__FILE__,__LINE__,s,1,0,0) - _stencil_val(__FILE__,__LINE__,s,-1,0,0))/2.;
#line 200

      g.y = (_stencil_val(__FILE__,__LINE__,s,0,1,0) - _stencil_val(__FILE__,__LINE__,s,0,-1,0))/2.;}

  double sc = _stencil_val(__FILE__,__LINE__,s,0,0,0), cmc = 4.*val_cm(cm,0,0,0), sum = val_cm(cm,0,0,0)*(1 << 2);
   { foreach_child() {
    _stencil_val(__FILE__,__LINE__,s,0,0,0) = sc;
    {
#line 206

      _stencil_val(__FILE__,__LINE__,s,0,0,0) += child.x*g.x*val_cm(cm,-child.x,0,0)/cmc;
#line 206

      _stencil_val(__FILE__,__LINE__,s,0,0,0) += child.y*g.y*val_cm(cm,0,-child.y,0)/cmc;}
    sum -= val_cm(cm,0,0,0);
  } end_foreach_child(); }
  IF (!(fabs(sum) < 1e-10)) _stencil_qassert (__FILE__,__LINE__,"/home/fpl/softwares/basilisk/src/grid/multigrid-common.h", 210, "fabs(sum) < 1e-10");
 }
#undef _IN_STENCIL

#endif

#line 211
}

static inline void refine_reset (Point point, scalar v)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 214 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    val(v,0,0,0) = 0.; end_foreach_child(); }

#if _call_refine_reset
}
#define _IN_STENCIL 1

#line 213
static void _refine_reset (Point point, scalar v)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 214 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    _stencil_val(__FILE__,__LINE__,v,0,0,0) = 0.; end_foreach_child(); }

#undef _IN_STENCIL

#endif

#line 217
}

static inline void refine_injection (Point point, scalar v)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 220 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

  double val = val(v,0,0,0);
   { foreach_child()
    val(v,0,0,0) = val; end_foreach_child(); }

#if _call_refine_injection
}
#define _IN_STENCIL 1

#line 219
static void _refine_injection (Point point, scalar v)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 220 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

  double val = _stencil_val(__FILE__,__LINE__,v,0,0,0);
   { foreach_child()
    _stencil_val(__FILE__,__LINE__,v,0,0,0) = val; end_foreach_child(); }

#undef _IN_STENCIL

#endif

#line 224
}

static scalar multigrid_init_scalar (scalar s, const char * name)
{
  s = cartesian_init_scalar (s, name);
  _attribute[s.i].prolongation = refine_bilinear;
  _attribute[s.i].restriction = restriction_average;
  return s;
}

static scalar multigrid_init_vertex_scalar (scalar s, const char * name)
{
  s = cartesian_init_vertex_scalar (s, name);
  _attribute[s.i].restriction = restriction_vertex;
  return s;
}

static vector multigrid_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  {
#line 244

    _attribute[v.y.i].restriction = no_restriction;
#line 244

    _attribute[v.x.i].restriction = no_restriction;}
  _attribute[v.x.i].restriction = restriction_face;
  return v;
}

void multigrid_debug (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 251 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

  cartesian_debug (point);

  FILE * plot = fopen ("plot", "a");
  if (point.level > 0) {
    char name[80] = "coarse";
    if (pid() > 0)
      sprintf (name, "coarse-%d", pid());
    FILE * fp = fopen (name, "w");
#line 271 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"
      double xc = x - child.x*Delta/2., yc = y - child.y*Delta/2.;
      for (int k = 0; k <= 1; k++)
 for (int l = 0; l <= 1; l++) {
   strongif (all) for (scalar v = *all, *_i51 = all; ((scalar *)&v)->i >= 0; v = *++_i51)
     fprintf (fp, "%g %g %g ",
       xc + k*child.x*Delta*2. + _attribute[v.i].d.x*Delta,
       yc + l*child.y*Delta*2. + _attribute[v.i].d.y*Delta,
       coarse(v,k*child.x,l*child.y,0));
   fputc ('\n', fp);
 }
      fprintf (ferr, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 3 t ''", name);
      fprintf (plot, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 3 t ''", name);
#line 302 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"
    fclose (fp);
  }

  if (is_coarse()) {
    char name[80] = "fine";
    if (pid() > 0)
      sprintf (name, "fine-%d", pid());
    FILE * fp = fopen (name, "w");
#line 324 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"
      double xf = x - Delta/4., yf = y - Delta/4.;
      for (int k = -2; k <= 3; k++)
 for (int l = -2; l <= 3; l++) {
   strongif (all) for (scalar v = *all, *_i52 = all; ((scalar *)&v)->i >= 0; v = *++_i52) {
     fprintf (fp, "%g %g ",
       xf + k*Delta/2. + _attribute[v.i].d.x*Delta/4.,
       yf + l*Delta/2. + _attribute[v.i].d.y*Delta/4.);
     if (allocated_child(k,l,0))
       fprintf (fp, "%g ", fine(v,k,l,0));
     else
       fputs ("n/a ", fp);
   }
   fputc ('\n', fp);
 }
      fprintf (ferr, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 2 t ''", name);
      fprintf (plot, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 2 t ''", name);
#line 362 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"
    fclose (fp);
  }
  fflush (ferr);
  fclose (plot);

#if _call_multigrid_debug
}
#define _IN_STENCIL 1
#define cartesian_debug _cartesian_debug

#line 250
static void _multigrid_debug (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 251 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

  cartesian_debug (point);

  FILE * plot = fopen ("plot", "a");
  IF (point.level > 0) {
    char name[80] = "coarse";
    IF (pid() > 0)
      sprintf (name, "coarse-%d", pid());
    FILE * fp = fopen (name, "w");
#line 271 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"
      double xc = x - child.x*Delta/2., yc = y - child.y*Delta/2.;
      for (int k = 0; k <= 1; k++)
 for (int l = 0; l <= 1; l++) {
   strongif (all) for (scalar v = *all, *_i51 = all; ((scalar *)&v)->i >= 0; v = *++_i51)
     _stencil_fprintf (__FILE__,__LINE__,fp, "%g %g %g ",
       xc + k*child.x*Delta*2. + _attribute[v.i].d.x*Delta,
       yc + l*child.y*Delta*2. + _attribute[v.i].d.y*Delta,
       _stencil_coarse(__FILE__,__LINE__,v,k*child.x,l*child.y,0));
   _stencil_fputc (__FILE__,__LINE__,'\n', fp);
 }
      _stencil_fprintf (__FILE__,__LINE__,ferr, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 3 t ''", name);
      _stencil_fprintf (__FILE__,__LINE__,plot, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 3 t ''", name);
#line 302 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"
    fclose (fp);
  }

  IF (is_coarse()) {
    char name[80] = "fine";
    IF (pid() > 0)
      sprintf (name, "fine-%d", pid());
    FILE * fp = fopen (name, "w");
#line 324 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"
      double xf = x - Delta/4., yf = y - Delta/4.;
      for (int k = -2; k <= 3; k++)
 for (int l = -2; l <= 3; l++) {
   strongif (all) for (scalar v = *all, *_i52 = all; ((scalar *)&v)->i >= 0; v = *++_i52) {
     _stencil_fprintf (__FILE__,__LINE__,fp, "%g %g ",
       xf + k*Delta/2. + _attribute[v.i].d.x*Delta/4.,
       yf + l*Delta/2. + _attribute[v.i].d.y*Delta/4.);
     IF (allocated_child(k,l,0))
       _stencil_fprintf (__FILE__,__LINE__,fp, "%g ", _stencil_fine(__FILE__,__LINE__,v,k,l,0));
     
       _stencil_fputs (__FILE__,__LINE__,"n/a ", fp);
   }
   _stencil_fputc (__FILE__,__LINE__,'\n', fp);
 }
      _stencil_fprintf (__FILE__,__LINE__,ferr, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 2 t ''", name);
      _stencil_fprintf (__FILE__,__LINE__,plot, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 2 t ''", name);
#line 362 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"
    fclose (fp);
  }
  fflush (ferr);
  fclose (plot);

#undef cartesian_debug
#undef _IN_STENCIL

#endif

#line 366
}

static void multigrid_restriction (scalar * list)
{
  scalar * listdef = NULL, * listc = NULL, * list2 = NULL;
  strongif (list) for (scalar s = *list, *_i53 = list; ((scalar *)&s)->i >= 0; s = *++_i53)
    if (!is_constant (s) && _attribute[s.i].block > 0) {
      if (_attribute[s.i].restriction == restriction_average) {
 listdef = list_add (listdef, s);
 list2 = list_add (list2, s);
      }
      else if (_attribute[s.i].restriction != no_restriction) {
 listc = list_add (listc, s);
 if (_attribute[s.i].face)
   {
#line 380

     list2 = list_add (list2, _attribute[s.i].v.x);
#line 380

     list2 = list_add (list2, _attribute[s.i].v.y);}
 else
   list2 = list_add (list2, s);
      }
    }

  if (listdef || listc) {
    for (int l = depth() - 1; l >= 0; l--) {
       { foreach_coarse_level(l){

#line 389 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"
 {
 strongif (listdef) for (scalar s = *listdef, *_i54 = listdef; ((scalar *)&s)->i >= 0; s = *++_i54)
  
     restriction_average (point, s);
 strongif (listc) for (scalar s = *listc, *_i55 = listc; ((scalar *)&s)->i >= 0; s = *++_i55) {
  
     _attribute[s.i].restriction (point, s);
 }
      } } end_foreach_coarse_level(); }
      { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list2, l); };
    }
    pfree (listdef,__func__,__FILE__,__LINE__);
    pfree (listc,__func__,__FILE__,__LINE__);
    pfree (list2,__func__,__FILE__,__LINE__);
  }
}

void multigrid_methods ()
{
  cartesian_methods();
  debug = multigrid_debug;
  init_scalar = multigrid_init_scalar;
  init_vertex_scalar = multigrid_init_vertex_scalar;
  init_face_vector = multigrid_init_face_vector;
  restriction = multigrid_restriction;
}







void subtree_size (scalar size, bool leaves)
{




   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h", .line = 428,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 428 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

    _stencil_val(__FILE__,__LINE__,size,0,0,0) = 1; } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 429
foreach(){

#line 428 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"

    val(size,0,0,0) = 1; } end_foreach(); }





  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar []){size,{-1}}), depth()); };
  for (int l = depth() - 1; l >= 0; l--) {
     { foreach_coarse_level(l){

#line 437 "/home/fpl/softwares/basilisk/src/grid/multigrid-common.h"
 {
      double sum = !leaves;
       { foreach_child()
 sum += val(size,0,0,0); end_foreach_child(); }
      val(size,0,0,0) = sum;
    } } end_foreach_coarse_level(); }
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar []){size,{-1}}), l); };
  }
}
#line 5 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"






#line 21 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"
int refine_cell (Point point, scalar * list, int flag, Cache * refined)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 22 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"

  int nr = 0;


  if (level > 0)
    for (int k = 0; k != 2*child.x; k += child.x)

      for (int l = 0; l != 2*child.y; l += child.y)




   if (aparent(k,l,m).pid >= 0 && is_leaf(aparent(k,l,m))) {
     Point p = point;


     p.level = point.level - 1;
     p.i = (point.i + 2)/2 + k;
     do { if (p.i < 2) p.i += 1 << p.level; else if (p.i >= 2 + (1 << p.level)) p.i -= 1 << p.level; } while(0);

       p.j = (point.j + 2)/2 + l;
       do { if (p.j < 2) p.j += 1 << p.level; else if (p.j >= 2 + (1 << p.level)) p.j -= 1 << p.level; } while(0);





     nr += refine_cell (p, list, flag, refined);
     aparent(k,l,m).flags |= flag;
   }



  increment_neighbors (point);

  int cflag = is_active(cell) ? (active|leaf) : leaf;
   { foreach_child()
    cell.flags |= cflag; end_foreach_child(); }


  strongif (list) for (scalar s = *list, *_i56 = list; ((scalar *)&s)->i >= 0; s = *++_i56)
    if (is_local(cell) || _attribute[s.i].face)
      _attribute[s.i].refine (point, s);


  cell.flags &= ~leaf;

#if 1
  if (is_border(cell)) {
     { foreach_child() {
      bool bord = false;
       { foreach_neighbor() {
 if (!is_local(cell) || (level > 0 && !is_local(aparent(0,0,0))))
   bord = true, foreach_neighbor_break();
 if (is_refined_check())
    { foreach_child()
     if (!is_local(cell))
       bord = true, foreach_child_break(); end_foreach_child(); }
 if (bord)
   foreach_neighbor_break();
      } end_foreach_neighbor(); }
      if (bord)
 cell.flags |= border;
    } end_foreach_child(); }
    if (refined)
      cache_append (refined, point, cell.flags);
    nr++;
  }
#endif
  return nr;

#if _call_refine_cell
}
#define _IN_STENCIL 1
#define increment_neighbors _increment_neighbors
#define refine_cell _refine_cell

#line 21
static int _refine_cell (Point point, scalar * list, int flag, Cache * refined)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 22 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"

  int nr = 0;


  IF (level > 0)
    for (int k = 0; k != 2*child.x; k += child.x)

      for (int l = 0; l != 2*child.y; l += child.y)




   IF (aparent(k,l,m).pid >= 0 && is_leaf(aparent(k,l,m))) {
     Point p = point;


     p.level = point.level - 1;
     p.i = (point.i + 2)/2 + k;
     do { IF (p.i < 2) p.i += 1 << p.level; IF (p.i >= 2 + (1 << p.level)) p.i -= 1 << p.level; } while(0);

       p.j = (point.j + 2)/2 + l;
       do { IF (p.j < 2) p.j += 1 << p.level; IF (p.j >= 2 + (1 << p.level)) p.j -= 1 << p.level; } while(0);





     nr += refine_cell (p, list, flag, refined);
     aparent(k,l,m).flags |= flag;
   }



  increment_neighbors (point);

  int cflag = is_active(cell) ? (active|leaf) : leaf;
   { foreach_child()
    cell.flags |= cflag; end_foreach_child(); }


  strongif (list) for (scalar s = *list, *_i56 = list; ((scalar *)&s)->i >= 0; s = *++_i56)
    IF (is_local(cell) || _attribute[s.i].face)
      _attribute[s.i].refine (point, s);


  cell.flags &= ~leaf;

#if 1
  IF (is_border(cell)) {
     { foreach_child() {
      bool bord = false;
       { foreach_neighbor() {
 IF (!is_local(cell) || (level > 0 && !is_local(aparent(0,0,0))))
   bord = true, foreach_neighbor_break();
 IF (is_refined_check())
    { foreach_child()
     IF (!is_local(cell))
       bord = true, foreach_child_break(); end_foreach_child(); }
 IF (bord)
   foreach_neighbor_break();
      } end_foreach_neighbor(); }
      IF (bord)
 cell.flags |= border;
    } end_foreach_child(); }
    IF (refined)
      cache_append (refined, point, cell.flags);
    nr++;
  }
#endif
  return nr;

#undef increment_neighbors
#undef refine_cell
#undef _IN_STENCIL

#endif

#line 92
}





bool coarsen_cell (Point point, scalar * list)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 99 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"




  int pid = cell.pid;
   { foreach_child()
    if (cell.neighbors || (cell.pid < 0 && cell.pid != pid))
      return false; end_foreach_child(); }



  strongif (list) for (scalar s = *list, *_i57 = list; ((scalar *)&s)->i >= 0; s = *++_i57) {
    _attribute[s.i].restriction (point, s);
    if (_attribute[s.i].coarsen)
      _attribute[s.i].coarsen (point, s);
  }


  cell.flags |= leaf;


  decrement_neighbors (point);

#if 1
  if (!is_local(cell)) {
    cell.flags &= ~(active|border);
     { foreach_neighbor(1)
      if (cell.neighbors)
  { foreach_child()
   if (allocated(0,0,0) && is_local(cell) && !is_border(cell))
     cell.flags |= border; end_foreach_child(); } end_foreach_neighbor(); }
  }
#endif

  return true;

#if _call_coarsen_cell
}
#define _IN_STENCIL 1
#define decrement_neighbors _decrement_neighbors

#line 98
static bool _coarsen_cell (Point point, scalar * list)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 99 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"




  int pid = cell.pid;
   { foreach_child()
    IF (cell.neighbors || (cell.pid < 0 && cell.pid != pid))
      return false; end_foreach_child(); }



  strongif (list) for (scalar s = *list, *_i57 = list; ((scalar *)&s)->i >= 0; s = *++_i57) {
    _attribute[s.i].restriction (point, s);
    IF (_attribute[s.i].coarsen)
      _attribute[s.i].coarsen (point, s);
  }


  cell.flags |= leaf;


  decrement_neighbors (point);

#if 1
  IF (!is_local(cell)) {
    cell.flags &= ~(active|border);
     { foreach_neighbor(1)
      IF (cell.neighbors)
  { foreach_child()
   IF (allocated(0,0,0) && is_local(cell) && !is_border(cell))
     cell.flags |= border; end_foreach_child(); } end_foreach_neighbor(); }
  }
#endif

  return true;

#undef decrement_neighbors
#undef _IN_STENCIL

#endif

#line 134
}

void coarsen_cell_recursive (Point point, scalar * list)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 137 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"



   { foreach_child()
    if (cell.neighbors)
       { foreach_neighbor(1)
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
   coarsen_cell_recursive (point, list); end_foreach_neighbor(); } end_foreach_child(); }

  if (!(coarsen_cell (point, list))) qassert ("/home/fpl/softwares/basilisk/src/grid/tree-common.h", 146, "coarsen_cell (point, list)");

#if _call_coarsen_cell_recursive
}
#define _IN_STENCIL 1
#define coarsen_cell _coarsen_cell
#define coarsen_cell_recursive _coarsen_cell_recursive

#line 136
static void _coarsen_cell_recursive (Point point, scalar * list)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 137 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"



   { foreach_child()
    IF (cell.neighbors)
       { foreach_neighbor(1)
 IF ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
   coarsen_cell_recursive (point, list); end_foreach_neighbor(); } end_foreach_child(); }

  IF (!(coarsen_cell (point, list))) _stencil_qassert (__FILE__,__LINE__,"/home/fpl/softwares/basilisk/src/grid/tree-common.h", 146, "coarsen_cell (point, list)");

#undef coarsen_cell
#undef coarsen_cell_recursive
#undef _IN_STENCIL

#endif

#line 147
}

void mpi_boundary_refine (scalar *);
void mpi_boundary_coarsen (int, int);
void mpi_boundary_update (scalar *);

typedef struct {
  int nc, nf;
} astats;

struct Adapt {
  scalar * slist;
  double * max;
  int maxlevel;
  int minlevel;
  scalar * list;
};


astats adapt_wavelet (struct Adapt p)
{ trace ("adapt_wavelet", "/home/fpl/softwares/basilisk/src/grid/tree-common.h", 167);
  scalar * list = p.list;

  if (is_constant(cm)) {
    if (list == NULL || list == all)
      list = list_copy (all);
    boundary_internal ((scalar *)(list), "/home/fpl/softwares/basilisk/src/grid/tree-common.h", 173);
    restriction (p.slist);
  }
  else {
    if (list == NULL || list == all) {
      list = list_copy (((scalar []){cm,fm.x,fm.y,{-1}}));
      strongif (all) for (scalar s = *all, *_i58 = all; ((scalar *)&s)->i >= 0; s = *++_i58)
 list = list_add (list, s);
    }
    boundary_internal ((scalar *)(list), "/home/fpl/softwares/basilisk/src/grid/tree-common.h", 182);
    scalar * listr = list_concat (p.slist, ((scalar []){cm,{-1}}));
    restriction (listr);
    pfree (listr,__func__,__FILE__,__LINE__);
  }

  astats st = {0, 0};
  scalar * listc = NULL;
  strongif (list) for (scalar s = *list, *_i59 = list; ((scalar *)&s)->i >= 0; s = *++_i59)
    if (!is_constant(s) && _attribute[s.i].restriction != no_restriction)
      listc = list_add (listc, s);


  if (p.minlevel < 1)
    p.minlevel = 1;
  ((Tree *)grid)->refined.n = 0;
  static const int refined = 1 << user, too_fine = 1 << (user + 1);
   { foreach_cell(){

#line 199 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"
 {
    if (is_active(cell)) {
      static const int too_coarse = 1 << (user + 2);
      if (is_leaf (cell)) {
 if (cell.flags & too_coarse) {
   cell.flags &= ~too_coarse;
   refine_cell (point, listc, refined, &((Tree *)grid)->refined);
   st.nf++;
 }
 continue;
      }
      else {
 if (cell.flags & refined) {

   cell.flags &= ~too_coarse;
   continue;
 }

 bool local = is_local(cell);
 if (!local)
    { foreach_child()
     if (is_local(cell))
       local = true, foreach_child_break(); end_foreach_child(); }
 if (local) {
   int i = 0;
   static const int just_fine = 1 << (user + 3);
   strongif (p.slist) for (scalar s = *p.slist, *_i60 = p.slist; ((scalar *)&s)->i >= 0; s = *++_i60) {
     double max = p.max[i++], sc[1 << 2];
     int c = 0;
      { foreach_child()
       sc[c++] = val(s,0,0,0); end_foreach_child(); }
     _attribute[s.i].prolongation (point, s);
     c = 0;
      { foreach_child() {
       double e = fabs(sc[c] - val(s,0,0,0));
       if (e > max && level < p.maxlevel) {
  cell.flags &= ~too_fine;
  cell.flags |= too_coarse;
       }
       else if ((e <= max/1.5 || level > p.maxlevel) &&
         !(cell.flags & (too_coarse|just_fine))) {
  if (level >= p.minlevel)
    cell.flags |= too_fine;
       }
       else if (!(cell.flags & too_coarse)) {
  cell.flags &= ~too_fine;
  cell.flags |= just_fine;
       }
       val(s,0,0,0) = sc[c++];
     } end_foreach_child(); }
   }
    { foreach_child() {
     cell.flags &= ~just_fine;
     if (!is_leaf(cell)) {
       cell.flags &= ~too_coarse;
       if (level >= p.maxlevel)
  cell.flags |= too_fine;
     }
     else if (!is_active(cell))
       cell.flags &= ~too_coarse;
   } end_foreach_child(); }
 }
      }
    }
    else
      continue;
  } } end_foreach_cell(); }
  mpi_boundary_refine (listc);



  for (int l = depth(); l >= 0; l--) {
     { foreach_cell(){

#line 271 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"

      if (!(cell.pid < 0)) {
 if (level == l) {
   if (!is_leaf(cell)) {
     if (cell.flags & refined)

       cell.flags &= ~(refined|too_fine);
     else if (cell.flags & too_fine) {
       if (is_local(cell) && coarsen_cell (point, listc))
  st.nc++;
       cell.flags &= ~too_fine;
     }
   }
   if (cell.flags & too_fine)
     cell.flags &= ~too_fine;
   else if (level > 0 && (aparent(0,0,0).flags & too_fine))
     aparent(0,0,0).flags &= ~too_fine;
   continue;
 }
 else if (is_leaf(cell))
   continue;
      } } end_foreach_cell(); }
    mpi_boundary_coarsen (l, too_fine);
  }
  pfree (listc,__func__,__FILE__,__LINE__);

  mpi_all_reduce (st.nf, MPI_INT, MPI_SUM);
  mpi_all_reduce (st.nc, MPI_INT, MPI_SUM);
  if (st.nc || st.nf)
    mpi_boundary_update (list);

  if (list != p.list)
    pfree (list,__func__,__FILE__,__LINE__);

  { astats _ret =  st; end_trace("adapt_wavelet", "/home/fpl/softwares/basilisk/src/grid/tree-common.h", 305);  return _ret; }
 end_trace("adapt_wavelet", "/home/fpl/softwares/basilisk/src/grid/tree-common.h", 306); }
#line 328 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"
static void refine_level (int depth)
{
  int refined;
  do {
    refined = 0;
    ((Tree *)grid)->refined.n = 0;
     { foreach_leaf(){

#line 334 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"

      if (level < depth) {
 refine_cell (point, NULL, 0, &((Tree *)grid)->refined);
 refined++;
 continue;
      } } end_foreach_leaf(); }
    mpi_all_reduce (refined, MPI_INT, MPI_SUM);
    if (refined) {
      mpi_boundary_refine (NULL);
      mpi_boundary_update (NULL);
    }
  } while (refined);
}
#line 373 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"

static void halo_face (vectorl vl)
{ trace ("halo_face", "/home/fpl/softwares/basilisk/src/grid/tree-common.h", 375);
  {
#line 376

    strongif (vl.x) for (scalar s = *vl.x, *_i61 = vl.x; ((scalar *)&s)->i >= 0; s = *++_i61)
      _attribute[s.i].dirty = 2;
#line 376

    strongif (vl.y) for (scalar s = *vl.y, *_i61 = vl.y; ((scalar *)&s)->i >= 0; s = *++_i61)
      _attribute[s.i].dirty = 2;}

  for (int l = depth() - 1; l >= 0; l--)
     { foreach_halo (prolongation, l){

#line 381 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"

      {
#line 382

        if (vl.x) {
#line 392 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"
   if ((!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0))
     strongif (vl.x) for (scalar s = *vl.x, *_i62 = vl.x; ((scalar *)&s)->i >= 0; s = *++_i62)
       val(s,0,0,0) = (fine(s,0,0,0) + fine(s,0,1,0))/2.;
   if ((!is_leaf (neighbor(1,0,0)) && neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0))
     strongif (vl.x) for (scalar s = *vl.x, *_i63 = vl.x; ((scalar *)&s)->i >= 0; s = *++_i63)
       val(s,1,0,0) = (fine(s,2,0,0) + fine(s,2,1,0))/2.;
#line 408 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"
 }
#line 382

        if (vl.y) {
#line 392 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"
   if ((!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0))
     strongif (vl.y) for (scalar s = *vl.y, *_i62 = vl.y; ((scalar *)&s)->i >= 0; s = *++_i62)
       val(s,0,0,0) = (fine(s,0,0,0) + fine(s,1,0,0))/2.;
   if ((!is_leaf (neighbor(0,1,0)) && neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0))
     strongif (vl.y) for (scalar s = *vl.y, *_i63 = vl.y; ((scalar *)&s)->i >= 0; s = *++_i63)
       val(s,0,1,0) = (fine(s,0,2,0) + fine(s,1,2,0))/2.;
#line 408 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"
 }} } end_foreach_halo(); }
 end_trace("halo_face", "/home/fpl/softwares/basilisk/src/grid/tree-common.h", 409); }



static scalar tree_init_scalar (scalar s, const char * name)
{
  s = multigrid_init_scalar (s, name);
  _attribute[s.i].refine = _attribute[s.i].prolongation;
  return s;
}

static void prolongation_vertex (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 421 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"


  fine(s,1,1,0) = (val(s,0,0,0) + val(s,1,0,0) + val(s,0,1,0) + val(s,1,1,0))/4.;





  for (int i = 0; i <= 1; i++) {
    for (int j = 0; j <= 1; j++)





      if (allocated_child(2*i,2*j,0))
 fine(s,2*i,2*j,0) = val(s,i,j,0);


    {
#line 440

      if (neighbor(i,0,0).neighbors) {

 fine(s,2*i,1,0) = (val(s,i,0,0) + val(s,i,1,0))/2.;
#line 453 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"
      }
#line 440

      if (neighbor(0,i,0).neighbors) {

 fine(s,1,2*i,0) = (val(s,0,i,0) + val(s,1,i,0))/2.;
#line 453 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"
      }}
  }

#if _call_prolongation_vertex
}
#define _IN_STENCIL 1

#line 420
static void _prolongation_vertex (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 421 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"


  _stencil_fine(__FILE__,__LINE__,s,1,1,0) = (_stencil_val(__FILE__,__LINE__,s,0,0,0) + _stencil_val(__FILE__,__LINE__,s,1,0,0) + _stencil_val(__FILE__,__LINE__,s,0,1,0) + _stencil_val(__FILE__,__LINE__,s,1,1,0))/4.;





  for (int i = 0; i <= 1; i++) {
    for (int j = 0; j <= 1; j++)





      IF (allocated_child(2*i,2*j,0))
 _stencil_fine(__FILE__,__LINE__,s,2*i,2*j,0) = _stencil_val(__FILE__,__LINE__,s,i,j,0);


    {
#line 440

      IF (neighbor(i,0,0).neighbors) {

 _stencil_fine(__FILE__,__LINE__,s,2*i,1,0) = (_stencil_val(__FILE__,__LINE__,s,i,0,0) + _stencil_val(__FILE__,__LINE__,s,i,1,0))/2.;
#line 453 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"
      }
#line 440

      IF (neighbor(0,i,0).neighbors) {

 _stencil_fine(__FILE__,__LINE__,s,1,2*i,0) = (_stencil_val(__FILE__,__LINE__,s,0,i,0) + _stencil_val(__FILE__,__LINE__,s,1,i,0))/2.;
#line 453 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"
      }}
  }

#undef _IN_STENCIL

#endif

#line 455
}

static scalar tree_init_vertex_scalar (scalar s, const char * name)
{
  s = multigrid_init_vertex_scalar (s, name);
  _attribute[s.i].refine = _attribute[s.i].prolongation = prolongation_vertex;
  return s;
}


#line 464

static void refine_face_x (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 466 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"

  vector v = _attribute[s.i].v;

  if (!(!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) &&
      (is_local(cell) || is_local(neighbor(-1,0,0)))) {
    double g1 = (val(v.x,0,+1,0) - val(v.x,0,-1,0))/8.;
    for (int j = 0; j <= 1; j++)
      fine(v.x,0,j,0) = val(v.x,0,0,0) + (2*j - 1)*g1;
  }
  if (!(!is_leaf (neighbor(1,0,0)) && neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0) && neighbor(1,0,0).neighbors &&
      (is_local(cell) || is_local(neighbor(1,0,0)))) {
    double g1 = (val(v.x,1,+1,0) - val(v.x,1,-1,0))/8.;
    for (int j = 0; j <= 1; j++)
      fine(v.x,2,j,0) = val(v.x,1,0,0) + (2*j - 1)*g1;
  }
  if (is_local(cell)) {
    double g1 = (val(v.x,0,+1,0) - val(v.x,0,-1,0) + val(v.x,1,+1,0) - val(v.x,1,-1,0))/16.;
    for (int j = 0; j <= 1; j++)
      fine(v.x,1,j,0) = (val(v.x,0,0,0) + val(v.x,1,0,0))/2. + (2*j - 1)*g1;
  }
#line 511 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"

#if _call_refine_face_x
}
#define _IN_STENCIL 1

#line 465
static void _refine_face_x (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 466 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"

  vector v = _attribute[s.i].v;

  IF (!(!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) &&
      (is_local(cell) || is_local(neighbor(-1,0,0)))) {
    double g1 = (_stencil_val(__FILE__,__LINE__,v.x,0,+1,0) - _stencil_val(__FILE__,__LINE__,v.x,0,-1,0))/8.;
    for (int j = 0; j <= 1; j++)
      _stencil_fine(__FILE__,__LINE__,v.x,0,j,0) = _stencil_val(__FILE__,__LINE__,v.x,0,0,0) + (2*j - 1)*g1;
  }
  IF (!(!is_leaf (neighbor(1,0,0)) && neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0) && neighbor(1,0,0).neighbors &&
      (is_local(cell) || is_local(neighbor(1,0,0)))) {
    double g1 = (_stencil_val(__FILE__,__LINE__,v.x,1,+1,0) - _stencil_val(__FILE__,__LINE__,v.x,1,-1,0))/8.;
    for (int j = 0; j <= 1; j++)
      _stencil_fine(__FILE__,__LINE__,v.x,2,j,0) = _stencil_val(__FILE__,__LINE__,v.x,1,0,0) + (2*j - 1)*g1;
  }
  IF (is_local(cell)) {
    double g1 = (_stencil_val(__FILE__,__LINE__,v.x,0,+1,0) - _stencil_val(__FILE__,__LINE__,v.x,0,-1,0) + _stencil_val(__FILE__,__LINE__,v.x,1,+1,0) - _stencil_val(__FILE__,__LINE__,v.x,1,-1,0))/16.;
    for (int j = 0; j <= 1; j++)
      _stencil_fine(__FILE__,__LINE__,v.x,1,j,0) = (_stencil_val(__FILE__,__LINE__,v.x,0,0,0) + _stencil_val(__FILE__,__LINE__,v.x,1,0,0))/2. + (2*j - 1)*g1;
  }
#line 511 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"

#undef _IN_STENCIL

#endif

#line 511
}
#line 464

static void refine_face_y (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 466 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"

  vector v = _attribute[s.i].v;

  if (!(!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) &&
      (is_local(cell) || is_local(neighbor(0,-1,0)))) {
    double g1 = (val(v.y,+1,0,0) - val(v.y,-1,0,0))/8.;
    for (int j = 0; j <= 1; j++)
      fine(v.y,j,0,0) = val(v.y,0,0,0) + (2*j - 1)*g1;
  }
  if (!(!is_leaf (neighbor(0,1,0)) && neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0) && neighbor(0,1,0).neighbors &&
      (is_local(cell) || is_local(neighbor(0,1,0)))) {
    double g1 = (val(v.y,+1,1,0) - val(v.y,-1,1,0))/8.;
    for (int j = 0; j <= 1; j++)
      fine(v.y,j,2,0) = val(v.y,0,1,0) + (2*j - 1)*g1;
  }
  if (is_local(cell)) {
    double g1 = (val(v.y,+1,0,0) - val(v.y,-1,0,0) + val(v.y,+1,1,0) - val(v.y,-1,1,0))/16.;
    for (int j = 0; j <= 1; j++)
      fine(v.y,j,1,0) = (val(v.y,0,0,0) + val(v.y,0,1,0))/2. + (2*j - 1)*g1;
  }
#line 511 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"

#if _call_refine_face_y
}
#define _IN_STENCIL 1

#line 465
static void _refine_face_y (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 466 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"

  vector v = _attribute[s.i].v;

  IF (!(!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) &&
      (is_local(cell) || is_local(neighbor(0,-1,0)))) {
    double g1 = (_stencil_val(__FILE__,__LINE__,v.y,0,+1,0) - _stencil_val(__FILE__,__LINE__,v.y,0,-1,0))/8.;
    for (int j = 0; j <= 1; j++)
      _stencil_fine(__FILE__,__LINE__,v.y,0,j,0) = _stencil_val(__FILE__,__LINE__,v.y,0,0,0) + (2*j - 1)*g1;
  }
  IF (!(!is_leaf (neighbor(0,1,0)) && neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0) && neighbor(0,1,0).neighbors &&
      (is_local(cell) || is_local(neighbor(0,1,0)))) {
    double g1 = (_stencil_val(__FILE__,__LINE__,v.y,1,+1,0) - _stencil_val(__FILE__,__LINE__,v.y,1,-1,0))/8.;
    for (int j = 0; j <= 1; j++)
      _stencil_fine(__FILE__,__LINE__,v.y,2,j,0) = _stencil_val(__FILE__,__LINE__,v.y,1,0,0) + (2*j - 1)*g1;
  }
  IF (is_local(cell)) {
    double g1 = (_stencil_val(__FILE__,__LINE__,v.y,0,+1,0) - _stencil_val(__FILE__,__LINE__,v.y,0,-1,0) + _stencil_val(__FILE__,__LINE__,v.y,1,+1,0) - _stencil_val(__FILE__,__LINE__,v.y,1,-1,0))/16.;
    for (int j = 0; j <= 1; j++)
      _stencil_fine(__FILE__,__LINE__,v.y,1,j,0) = (_stencil_val(__FILE__,__LINE__,v.y,0,0,0) + _stencil_val(__FILE__,__LINE__,v.y,1,0,0))/2. + (2*j - 1)*g1;
  }
#line 511 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"

#undef _IN_STENCIL

#endif

#line 511
}

void refine_face (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 514 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"

  vector v = _attribute[s.i].v;
  {
#line 516

    _attribute[v.x.i].prolongation (point, v.x);
#line 516

    _attribute[v.y.i].prolongation (point, v.y);}

#if _call_refine_face
}
#define _IN_STENCIL 1

#line 513
static void _refine_face (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 514 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"

  vector v = _attribute[s.i].v;
  {
#line 516

    _attribute[v.x.i].prolongation (point, v.x);
#line 516

    _attribute[v.y.i].prolongation (point, v.y);}

#undef _IN_STENCIL

#endif

#line 518
}

void refine_face_solenoidal (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 521 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"

  refine_face (point, s);

  if (is_local(cell)) {

    vector v = _attribute[s.i].v;
    double d[1 << 2], p[1 << 2];
    int i = 0;
     { foreach_child() {
      d[i] = 0.;
      {
#line 531

 d[i] += val(v.x,1,0,0) - val(v.x,0,0,0);
#line 531

 d[i] += val(v.y,0,1,0) - val(v.y,0,0,0);}
      i++;
    } end_foreach_child(); }

    p[0] = 0.;
    p[1] = (3.*d[3] + d[0])/4. + d[2]/2.;
    p[2] = (d[3] + 3.*d[0])/4. + d[2]/2.;
    p[3] = (d[3] + d[0])/2. + d[2];
    fine(v.x,1,1,0) += p[1] - p[0];
    fine(v.x,1,0,0) += p[3] - p[2];
    fine(v.y,0,1,0) += p[0] - p[2];
    fine(v.y,1,1,0) += p[1] - p[3];
#line 571 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"
  }


#if _call_refine_face_solenoidal
}
#define _IN_STENCIL 1
#define refine_face _refine_face

#line 520
static void _refine_face_solenoidal (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 521 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"

  refine_face (point, s);

  IF (is_local(cell)) {

    vector v = _attribute[s.i].v;
    double d[1 << 2], p[1 << 2];
    int i = 0;
     { foreach_child() {
      d[i] = 0.;
      {
#line 531

 d[i] += _stencil_val(__FILE__,__LINE__,v.x,1,0,0) - _stencil_val(__FILE__,__LINE__,v.x,0,0,0);
#line 531

 d[i] += _stencil_val(__FILE__,__LINE__,v.y,0,1,0) - _stencil_val(__FILE__,__LINE__,v.y,0,0,0);}
      i++;
    } end_foreach_child(); }

    p[0] = 0.;
    p[1] = (3.*d[3] + d[0])/4. + d[2]/2.;
    p[2] = (d[3] + 3.*d[0])/4. + d[2]/2.;
    p[3] = (d[3] + d[0])/2. + d[2];
    _stencil_fine(__FILE__,__LINE__,v.x,1,1,0) += p[1] - p[0];
    _stencil_fine(__FILE__,__LINE__,v.x,1,0,0) += p[3] - p[2];
    _stencil_fine(__FILE__,__LINE__,v.y,0,1,0) += p[0] - p[2];
    _stencil_fine(__FILE__,__LINE__,v.y,1,1,0) += p[1] - p[3];
#line 571 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"
  }


#undef refine_face
#undef _IN_STENCIL

#endif

#line 573
}

vector tree_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  {
#line 578

    _attribute[v.x.i].restriction = _attribute[v.x.i].refine = no_restriction;
#line 578

    _attribute[v.y.i].restriction = _attribute[v.y.i].refine = no_restriction;}
  _attribute[v.x.i].restriction = restriction_face;
  _attribute[v.x.i].refine = refine_face;
  {
#line 582

    _attribute[v.x.i].prolongation = refine_face_x;
#line 582

    _attribute[v.y.i].prolongation = refine_face_y;}
  return v;
}


static void tree_boundary_level (scalar * list, int l)
{ trace ("tree_boundary_level", "/home/fpl/softwares/basilisk/src/grid/tree-common.h", 589);
  int depth = l < 0 ? depth() : l;

  if (tree_is_full()) {
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, depth); };
    { ; end_trace("tree_boundary_level", "/home/fpl/softwares/basilisk/src/grid/tree-common.h", 594);  return; }
  }

  scalar * listdef = NULL, * listc = NULL, * list2 = NULL, * vlist = NULL;
  strongif (list) for (scalar s = *list, *_i64 = list; ((scalar *)&s)->i >= 0; s = *++_i64)
    if (!is_constant (s)) {
      if (_attribute[s.i].restriction == restriction_average) {
 listdef = list_add (listdef, s);
 list2 = list_add (list2, s);
      }
      else if (_attribute[s.i].restriction != no_restriction) {
 listc = list_add (listc, s);
 if (_attribute[s.i].face)
   {
#line 607

     list2 = list_add (list2, _attribute[s.i].v.x);
#line 607

     list2 = list_add (list2, _attribute[s.i].v.y);}
 else {
   list2 = list_add (list2, s);
   if (_attribute[s.i].restriction == restriction_vertex)
     vlist = list_add (vlist, s);
 }
      }
    }

  if (vlist)






     { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/grid/tree-common.h", .line = 624,
    .each = "foreach_vertex", .first = _first_call
  };
foreach_vertex_stencil(){

#line 624 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"
 {
      IF ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) || (!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) ||
   (!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) || (!is_leaf (neighbor(-1,-1,0)) && neighbor(-1,-1,0).neighbors && neighbor(-1,-1,0).pid >= 0)) {

 strongif (vlist) for (scalar s = *vlist, *_i65 = vlist; ((scalar *)&s)->i >= 0; s = *++_i65)
   _stencil_val(__FILE__,__LINE__,s,0,0,0) = is_vertex (child(0,0,0)) ? _stencil_fine(__FILE__,__LINE__,s,0,0,0) : nodata;
      }
      
 {
#line 632

   IF (child.y == 1 &&
       ((!is_leaf(cell) && !cell.neighbors && cell.pid >= 0) || (!is_leaf(neighbor(-1,0,0)) && !neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0))) {

     strongif (vlist) for (scalar s = *vlist, *_i66 = vlist; ((scalar *)&s)->i >= 0; s = *++_i66)
       _stencil_val(__FILE__,__LINE__,s,0,0,0) = is_vertex(neighbor(0,-1,0)) && is_vertex(neighbor(0,1,0)) ?
  (_stencil_val(__FILE__,__LINE__,s,0,-1,0) + _stencil_val(__FILE__,__LINE__,s,0,1,0))/2. : nodata;
   }
#line 632

   IF (child.x == 1 &&
       ((!is_leaf(cell) && !cell.neighbors && cell.pid >= 0) || (!is_leaf(neighbor(0,-1,0)) && !neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0))) {

     strongif (vlist) for (scalar s = *vlist, *_i66 = vlist; ((scalar *)&s)->i >= 0; s = *++_i66)
       _stencil_val(__FILE__,__LINE__,s,0,0,0) = is_vertex(neighbor(-1,0,0)) && is_vertex(neighbor(1,0,0)) ?
  (_stencil_val(__FILE__,__LINE__,s,-1,0,0) + _stencil_val(__FILE__,__LINE__,s,1,0,0))/2. : nodata;
   }}
    } } end_foreach_vertex_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 640
foreach_vertex(){

#line 624 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"
 {
      if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) || (!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) ||
   (!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) || (!is_leaf (neighbor(-1,-1,0)) && neighbor(-1,-1,0).neighbors && neighbor(-1,-1,0).pid >= 0)) {

 strongif (vlist) for (scalar s = *vlist, *_i65 = vlist; ((scalar *)&s)->i >= 0; s = *++_i65)
   val(s,0,0,0) = is_vertex (child(0,0,0)) ? fine(s,0,0,0) : nodata;
      }
      else
 {
#line 632

   if (child.y == 1 &&
       ((!is_leaf(cell) && !cell.neighbors && cell.pid >= 0) || (!is_leaf(neighbor(-1,0,0)) && !neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0))) {

     strongif (vlist) for (scalar s = *vlist, *_i66 = vlist; ((scalar *)&s)->i >= 0; s = *++_i66)
       val(s,0,0,0) = is_vertex(neighbor(0,-1,0)) && is_vertex(neighbor(0,1,0)) ?
  (val(s,0,-1,0) + val(s,0,1,0))/2. : nodata;
   }
#line 632

   if (child.x == 1 &&
       ((!is_leaf(cell) && !cell.neighbors && cell.pid >= 0) || (!is_leaf(neighbor(0,-1,0)) && !neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0))) {

     strongif (vlist) for (scalar s = *vlist, *_i66 = vlist; ((scalar *)&s)->i >= 0; s = *++_i66)
       val(s,0,0,0) = is_vertex(neighbor(-1,0,0)) && is_vertex(neighbor(1,0,0)) ?
  (val(s,-1,0,0) + val(s,1,0,0))/2. : nodata;
   }}
    } } end_foreach_vertex(); }
#line 673 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"
  pfree (vlist,__func__,__FILE__,__LINE__);

  if (listdef || listc) {
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, list2, depth); };
    for (int l = depth - 1; l >= 0; l--) {
       { foreach_coarse_level(l){

#line 678 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"
 {
 strongif (listdef) for (scalar s = *listdef, *_i67 = listdef; ((scalar *)&s)->i >= 0; s = *++_i67)
   restriction_average (point, s);
 strongif (listc) for (scalar s = *listc, *_i68 = listc; ((scalar *)&s)->i >= 0; s = *++_i68)
   _attribute[s.i].restriction (point, s);
      } } end_foreach_coarse_level(); }
      { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, list2, l); };
    }
    pfree (listdef,__func__,__FILE__,__LINE__);
    pfree (listc,__func__,__FILE__,__LINE__);
    pfree (list2,__func__,__FILE__,__LINE__);
  }

  scalar * listr = NULL;
  vector * listf = NULL;
  strongif (list) for (scalar s = *list, *_i69 = list; ((scalar *)&s)->i >= 0; s = *++_i69)
    if (!is_constant (s) && _attribute[s.i].refine != no_restriction) {
      if (_attribute[s.i].face)
 listf = vectors_add (listf, _attribute[s.i].v);
      else
 listr = list_add (listr, s);
    }

  if (listr || listf) {
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, 0); };
    for (int i = 0; i < depth; i++) {
       { foreach_halo (prolongation, i){

#line 704 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"
 {
 strongif (listr) for (scalar s = *listr, *_i70 = listr; ((scalar *)&s)->i >= 0; s = *++_i70)
          _attribute[s.i].prolongation (point, s);
 strongif (listf) for (vector v = *listf, *_i71 = listf; ((scalar *)&v)->i >= 0; v = *++_i71)
   {
#line 708

     _attribute[v.x.i].prolongation (point, v.x);
#line 708

     _attribute[v.y.i].prolongation (point, v.y);}
      } } end_foreach_halo(); }
      { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, i + 1); };
    }
    pfree (listr,__func__,__FILE__,__LINE__);
    pfree (listf,__func__,__FILE__,__LINE__);
  }
 end_trace("tree_boundary_level", "/home/fpl/softwares/basilisk/src/grid/tree-common.h", 716); }

double treex (Point point) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 718 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"

  if (level == 0)
    return 0;

  double i = 2*child.x - child.y;
  if (i <= 1 && i >= -1) i = -i;




  return treex(parent) + i/(1 << 2*(level - 1));

#if _call_treex
}
#define _IN_STENCIL 1
#define treex _treex

#line 718
static double _treex (Point point) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 718 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"

  IF (level == 0)
    return 0;

  double i = 2*child.x - child.y;
  IF (i <= 1 && i >= -1) i = -i;




  return treex(parent) + i/(1 << 2*(level - 1));

#undef treex
#undef _IN_STENCIL

#endif

#line 729
}

double treey (Point point) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 731 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"

  if (level == 0)
    return 0;
  return treey(parent) + 4./(1 << 2*(level - 1));

#if _call_treey
}
#define _IN_STENCIL 1
#define treey _treey

#line 731
static double _treey (Point point) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 731 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"

  IF (level == 0)
    return 0;
  return treey(parent) + 4./(1 << 2*(level - 1));

#undef treey
#undef _IN_STENCIL

#endif

#line 735
}

void output_tree (FILE * fp)
{
   { foreach_cell(){

#line 739 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"

    if (cell.neighbors)
       { foreach_child()
 if (is_local(cell))
   fprintf (fp, "%g %g\n%g %g\n\n",
     treex(parent), treey(parent), treex(point), treey(point)); end_foreach_child(); }; } end_foreach_cell(); }
}


void tree_check ()
{ trace ("tree_check", "/home/fpl/softwares/basilisk/src/grid/tree-common.h", 749);


  long nleaves = 0, nactive = 0;
   { foreach_cell_all(){

#line 753 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"
 {
    if (is_leaf(cell)) {
      if (!(cell.pid >= 0)) qassert ("/home/fpl/softwares/basilisk/src/grid/tree-common.h", 755, "cell.pid >= 0");
      nleaves++;
    }
    if (is_local(cell))
      if (!(is_active(cell) || (!is_leaf(cell) && !cell.neighbors && cell.pid >= 0))) qassert ("/home/fpl/softwares/basilisk/src/grid/tree-common.h", 759, "is_active(cell) || is_prolongation(cell)");
    if (is_active(cell))
      nactive++;

    int neighbors = 0;
     { foreach_neighbor(1)
      if (allocated(0,0,0) && (!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
 neighbors++; end_foreach_neighbor(); }
    if (!(cell.neighbors == neighbors)) qassert ("/home/fpl/softwares/basilisk/src/grid/tree-common.h", 767, "cell.neighbors == neighbors");


    if (!cell.neighbors)
      if (!(!allocated_child(0,0,0))) qassert ("/home/fpl/softwares/basilisk/src/grid/tree-common.h", 771, "!allocated_child(0)");
  } } end_foreach_cell_all(); }


  long reachable = 0;
   { foreach_cell(){

#line 776 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"
 {
    if (is_active(cell))
      reachable++;
    else
      continue;
  } } end_foreach_cell(); }
  if (!(nactive == reachable)) qassert ("/home/fpl/softwares/basilisk/src/grid/tree-common.h", 782, "nactive == reachable");


  reachable = 0;
   { foreach_cell(){

#line 786 "/home/fpl/softwares/basilisk/src/grid/tree-common.h"

    if (is_leaf(cell)) {
      reachable++;
      continue;
    } } end_foreach_cell(); }
  if (!(nleaves == reachable)) qassert ("/home/fpl/softwares/basilisk/src/grid/tree-common.h", 791, "nleaves == reachable");
 end_trace("tree_check", "/home/fpl/softwares/basilisk/src/grid/tree-common.h", 792); }


static void tree_restriction (scalar * list) { trace ("tree_restriction", "/home/fpl/softwares/basilisk/src/grid/tree-common.h", 795);
  boundary_internal ((scalar *)(list), "/home/fpl/softwares/basilisk/src/grid/tree-common.h", 796);
  if (tree_is_full())
    multigrid_restriction (list);
 end_trace("tree_restriction", "/home/fpl/softwares/basilisk/src/grid/tree-common.h", 799); }

void tree_methods ()
{
  multigrid_methods();
  init_scalar = tree_init_scalar;
  init_vertex_scalar = tree_init_vertex_scalar;
  init_face_vector = tree_init_face_vector;
  boundary_level = tree_boundary_level;
  boundary_face = halo_face;
  restriction = tree_restriction;
}
#line 1663 "/home/fpl/softwares/basilisk/src/grid/tree.h"


void tree_periodic (int dir)
{
  int depth = grid ? depth() : -1;
  if (grid)
    free_grid();
  periodic (dir);
  if (depth >= 0)
    init_grid (1 << depth);
}


#if 1
#line 1 "grid/tree-mpi.h"
#line 1 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

int debug_iteration = -1;

void debug_mpi (FILE * fp1);

typedef struct {
  CacheLevel * halo;
  void * buf;
  MPI_Request r;
  int depth;
  int pid;
  int maxdepth;
} Rcv;

typedef struct {
  Rcv * rcv;
  char * name;
  int npid;
} RcvPid;

typedef struct {
  RcvPid * rcv, * snd;
} SndRcv;

typedef struct {
  Boundary parent;

  SndRcv mpi_level, mpi_level_root, restriction;
  Array * send, * receive;
} MpiBoundary;

static void cache_level_init (CacheLevel * c)
{
  c->p = NULL;
  c->n = c->nm = 0;
}

static void rcv_append (Point point, Rcv * rcv)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 39 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

  if (level > rcv->depth) {
    rcv->halo = (CacheLevel *) prealloc (rcv->halo, (level + 1)*sizeof(CacheLevel),__func__,__FILE__,__LINE__);
    for (int j = rcv->depth + 1; j <= level; j++)
      cache_level_init (&rcv->halo[j]);
    rcv->depth = level;
  }
  cache_level_append (&rcv->halo[level], point);
  if (level > rcv->maxdepth)
    rcv->maxdepth = level;

#if _call_rcv_append
}
#define _IN_STENCIL 1

#line 38
static void _rcv_append (Point point, Rcv * rcv)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 39 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

  IF (level > rcv->depth) {
    rcv->halo = (CacheLevel *) prealloc (rcv->halo, (level + 1)*sizeof(CacheLevel),__func__,__FILE__,__LINE__);
    for (int j = rcv->depth + 1; j <= level; j++)
      cache_level_init (&rcv->halo[j]);
    rcv->depth = level;
  }
  cache_level_append (&rcv->halo[level], point);
  IF (level > rcv->maxdepth)
    rcv->maxdepth = level;

#undef _IN_STENCIL

#endif

#line 49
}

void rcv_print (Rcv * rcv, FILE * fp, const char * prefix)
{
  for (int l = 0; l <= rcv->depth; l++)
    if (rcv->halo[l].n > 0)
       { foreach_cache_level(rcv->halo[l], l){

#line 55 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

 fprintf (fp, "%s%g %g %g %d %d\n", prefix, x, y, z, rcv->pid, level); } end_foreach_cache_level(); }
}

static void rcv_free_buf (Rcv * rcv)
{
  if (rcv->buf) {
    prof_start ("rcv_pid_receive");
    MPI_Wait (&rcv->r, MPI_STATUS_IGNORE);
    pfree (rcv->buf,__func__,__FILE__,__LINE__);
    rcv->buf = NULL;
    prof_stop();
  }
}

static void rcv_destroy (Rcv * rcv)
{
  rcv_free_buf (rcv);
  for (int i = 0; i <= rcv->depth; i++)
    if (rcv->halo[i].n > 0)
      pfree (rcv->halo[i].p,__func__,__FILE__,__LINE__);
  pfree (rcv->halo,__func__,__FILE__,__LINE__);
}

static RcvPid * rcv_pid_new (const char * name)
{
  RcvPid * r = ((RcvPid *) pcalloc (1, sizeof(RcvPid),__func__,__FILE__,__LINE__));
  r->name = pstrdup (name,__func__,__FILE__,__LINE__);
  return r;
}

static Rcv * rcv_pid_pointer (RcvPid * p, int pid)
{
  if (!(pid >= 0 && pid < npe())) qassert ("/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 88, "pid >= 0 && pid < npe()");

  int i;
  for (i = 0; i < p->npid; i++)
    if (pid == p->rcv[i].pid)
      break;

  if (i == p->npid) {
    p->rcv = (Rcv *) prealloc (p->rcv, (++p->npid)*sizeof(Rcv),__func__,__FILE__,__LINE__);
    Rcv * rcv = &p->rcv[p->npid-1];
    rcv->pid = pid;
    rcv->depth = rcv->maxdepth = 0;
    rcv->halo = ((CacheLevel *) pmalloc ((1)*sizeof(CacheLevel),__func__,__FILE__,__LINE__));
    rcv->buf = NULL;
    cache_level_init (&rcv->halo[0]);
  }
  return &p->rcv[i];
}

static void rcv_pid_append (RcvPid * p, int pid, Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 108 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

  rcv_append (point, rcv_pid_pointer (p, pid));

#if _call_rcv_pid_append
}
#define _IN_STENCIL 1
#define rcv_append _rcv_append

#line 107
static void _rcv_pid_append (RcvPid * p, int pid, Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 108 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

  rcv_append (point, rcv_pid_pointer (p, pid));

#undef rcv_append
#undef _IN_STENCIL

#endif

#line 110
}

static void rcv_pid_append_pids (RcvPid * p, Array * pids)
{

  for (int i = 0; i < p->npid; i++) {
    int pid = p->rcv[i].pid, j, * a;
    for (j = 0, a = pids->p; j < pids->len/sizeof(int); j++,a++)
      if (*a == pid)
 break;
    if (j == pids->len/sizeof(int))
      array_append (pids, &pid, sizeof(int));
  }
}

void rcv_pid_write (RcvPid * p, const char * name)
{
  for (int i = 0; i < p->npid; i++) {
    Rcv * rcv = &p->rcv[i];
    char fname[80];
    sprintf (fname, "%s-%d-%d", name, pid(), rcv->pid);
    FILE * fp = fopen (fname, "w");
    rcv_print (rcv, fp, "");
    fclose (fp);
  }
}

static void rcv_pid_print (RcvPid * p, FILE * fp, const char * prefix)
{
  for (int i = 0; i < p->npid; i++)
    rcv_print (&p->rcv[i], fp, prefix);
}

static void rcv_pid_destroy (RcvPid * p)
{
  for (int i = 0; i < p->npid; i++)
    rcv_destroy (&p->rcv[i]);
  pfree (p->rcv,__func__,__FILE__,__LINE__);
  pfree (p->name,__func__,__FILE__,__LINE__);
  pfree (p,__func__,__FILE__,__LINE__);
}

static Boundary * mpi_boundary = NULL;






void debug_mpi (FILE * fp1);

static void apply_bc (Rcv * rcv, scalar * list, scalar * listv,
        vector * listf, int l, MPI_Status s)
{
  double * b = rcv->buf;
   { foreach_cache_level(rcv->halo[l], l){

#line 165 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"
 {
    strongif (list) for (scalar s = *list, *_i72 = list; ((scalar *)&s)->i >= 0; s = *++_i72) {
      memcpy (&val(s,0,0,0), b, sizeof(double)*_attribute[s.i].block);
      b += _attribute[s.i].block;
    }
    strongif (listf) for (vector v = *listf, *_i73 = listf; ((scalar *)&v)->i >= 0; v = *++_i73)
      {
#line 171
 {
 memcpy (&val(v.x,0,0,0), b, sizeof(double)*_attribute[v.x.i].block);
 b += _attribute[v.x.i].block;
 if (*b != nodata && allocated(1,0,0))
   memcpy (&val(v.x,1,0,0), b, sizeof(double)*_attribute[v.x.i].block);
 b += _attribute[v.x.i].block;
      }
#line 171
 {
 memcpy (&val(v.y,0,0,0), b, sizeof(double)*_attribute[v.y.i].block);
 b += _attribute[v.y.i].block;
 if (*b != nodata && allocated(0,1,0))
   memcpy (&val(v.y,0,1,0), b, sizeof(double)*_attribute[v.y.i].block);
 b += _attribute[v.y.i].block;
      }}
    strongif (listv) for (scalar s = *listv, *_i74 = listv; ((scalar *)&s)->i >= 0; s = *++_i74) {
      for (int i = 0; i <= 1; i++)
 for (int j = 0; j <= 1; j++)







          {
     if (*b != nodata && allocated(i,j,0))
       memcpy (&val(s,i,j,0), b, sizeof(double)*_attribute[s.i].block);
     b += _attribute[s.i].block;
          }

    }
  } } end_foreach_cache_level(); }
  size_t size = b - (double *) rcv->buf;
  pfree (rcv->buf,__func__,__FILE__,__LINE__);
  rcv->buf = NULL;

  int rlen;
  MPI_Get_count (&s, MPI_DOUBLE, &rlen);
  if (rlen != size) {
    fprintf (ferr,
      "rlen (%d) != size (%ld), %d receiving from %d at level %d\n"
      "Calling debug_mpi(NULL)...\n"
      "Aborting...\n",
      rlen, size, pid(), rcv->pid, l);
    fflush (ferr);
    debug_mpi (NULL);
    MPI_Abort (MPI_COMM_WORLD, -2);
  }
}
#line 234 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"
static void mpi_recv_check (void * buf, int count, MPI_Datatype datatype,
       int source, int tag,
       MPI_Comm comm, MPI_Status * status,
       const char * name)
{
#line 269 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"
  int errorcode = MPI_Recv (buf, count, datatype, source, tag, comm, status);
  if (errorcode != MPI_SUCCESS) {
    char string[MPI_MAX_ERROR_STRING];
    int resultlen;
    MPI_Error_string (errorcode, string, &resultlen);
    fprintf (ferr,
      "ERROR MPI_Recv \"%s\" (count = %d, source = %d, tag = %d):\n%s\n"
      "Calling debug_mpi(NULL)...\n"
      "Aborting...\n",
      name, count, source, tag, string);
    fflush (ferr);
    debug_mpi (NULL);
    MPI_Abort (MPI_COMM_WORLD, -1);
  }





}


static int mpi_waitany (int count, MPI_Request array_of_requests[], int *indx,
   MPI_Status *status)
{ trace ("mpi_waitany", "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 293);
  { int _ret =  MPI_Waitany (count, array_of_requests, indx, status); end_trace("mpi_waitany", "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 294);  return _ret; }
 end_trace("mpi_waitany", "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 295); }

static int list_lenb (scalar * list) {
  int len = 0;
  strongif (list) for (scalar s = *list, *_i75 = list; ((scalar *)&s)->i >= 0; s = *++_i75)
    len += _attribute[s.i].block;
  return len;
}

static int vectors_lenb (vector * list) {
  int len = 0;
  strongif (list) for (vector v = *list, *_i76 = list; ((scalar *)&v)->i >= 0; v = *++_i76)
    len += _attribute[v.x.i].block;
  return len;
}

static void rcv_pid_receive (RcvPid * m, scalar * list, scalar * listv,
        vector * listf, int l)
{
  if (m->npid == 0)
    return;

  prof_start ("rcv_pid_receive");

  int len = list_lenb (list) + 2*2*vectors_lenb (listf) +
    (1 << 2)*list_lenb (listv);

  MPI_Request r[m->npid];
  Rcv * rrcv[m->npid];
  int nr = 0;
  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0) {
      if (!(!rcv->buf)) qassert ("/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 328, "!rcv->buf");
      rcv->buf = pmalloc (sizeof (double)*rcv->halo[l].n*len,__func__,__FILE__,__LINE__);






      MPI_Irecv (rcv->buf, rcv->halo[l].n*len, MPI_DOUBLE, rcv->pid,
   (l), MPI_COMM_WORLD, &r[nr]);
      rrcv[nr++] = rcv;






    }
  }


  if (nr > 0) {
    int i;
    MPI_Status s;
    mpi_waitany (nr, r, &i, &s);
    while (i != MPI_UNDEFINED) {
      Rcv * rcv = rrcv[i];
      if (!(l <= rcv->depth && rcv->halo[l].n > 0)) qassert ("/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 355, "l <= rcv->depth && rcv->halo[l].n > 0");
      if (!(rcv->buf)) qassert ("/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 356, "rcv->buf");
      apply_bc (rcv, list, listv, listf, l, s);
      mpi_waitany (nr, r, &i, &s);
    }
  }

  prof_stop();
}


static void rcv_pid_wait (RcvPid * m)
{ trace ("rcv_pid_wait", "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 367);

  for (int i = 0; i < m->npid; i++)
    rcv_free_buf (&m->rcv[i]);
 end_trace("rcv_pid_wait", "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 371); }

static void rcv_pid_send (RcvPid * m, scalar * list, scalar * listv,
     vector * listf, int l)
{
  if (m->npid == 0)
    return;

  prof_start ("rcv_pid_send");

  int len = list_lenb (list) + 2*2*vectors_lenb (listf) +
    (1 << 2)*list_lenb (listv);


  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0) {
      if (!(!rcv->buf)) qassert ("/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 388, "!rcv->buf");
      rcv->buf = pmalloc (sizeof (double)*rcv->halo[l].n*len,__func__,__FILE__,__LINE__);
      double * b = rcv->buf;
       { foreach_cache_level(rcv->halo[l], l){

#line 391 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"
 {
 strongif (list) for (scalar s = *list, *_i77 = list; ((scalar *)&s)->i >= 0; s = *++_i77) {
   memcpy (b, &val(s,0,0,0), sizeof(double)*_attribute[s.i].block);
   b += _attribute[s.i].block;
 }
 strongif (listf) for (vector v = *listf, *_i78 = listf; ((scalar *)&v)->i >= 0; v = *++_i78)
   {
#line 397
 {
     memcpy (b, &val(v.x,0,0,0), sizeof(double)*_attribute[v.x.i].block);
     b += _attribute[v.x.i].block;
     if (allocated(1,0,0))
       memcpy (b, &val(v.x,1,0,0), sizeof(double)*_attribute[v.x.i].block);
     else
       *b = nodata;
     b += _attribute[v.x.i].block;
   }
#line 397
 {
     memcpy (b, &val(v.y,0,0,0), sizeof(double)*_attribute[v.y.i].block);
     b += _attribute[v.y.i].block;
     if (allocated(0,1,0))
       memcpy (b, &val(v.y,0,1,0), sizeof(double)*_attribute[v.y.i].block);
     else
       *b = nodata;
     b += _attribute[v.y.i].block;
   }}
 strongif (listv) for (scalar s = *listv, *_i79 = listv; ((scalar *)&s)->i >= 0; s = *++_i79) {
   for (int i = 0; i <= 1; i++)
     for (int j = 0; j <= 1; j++)
#line 418 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"
       {
  if (allocated(i,j,0))
    memcpy (b, &val(s,i,j,0), sizeof(double)*_attribute[s.i].block);
  else
    *b = nodata;
  b += _attribute[s.i].block;
       }

 }
      } } end_foreach_cache_level(); }





      MPI_Isend (rcv->buf, (b - (double *) rcv->buf),
   MPI_DOUBLE, rcv->pid, (l), MPI_COMM_WORLD,
   &rcv->r);
    }
  }

  prof_stop();
}

static void rcv_pid_sync (SndRcv * m, scalar * list, int l)
{
  scalar * listr = NULL, * listv = NULL;
  vector * listf = NULL;
  strongif (list) for (scalar s = *list, *_i80 = list; ((scalar *)&s)->i >= 0; s = *++_i80)
    if (!is_constant(s) && _attribute[s.i].block > 0) {
      if (_attribute[s.i].face)
 listf = vectors_add (listf, _attribute[s.i].v);
      else if (_attribute[s.i].restriction == restriction_vertex)
 listv = list_add (listv, s);
      else
 listr = list_add (listr, s);
    }
  rcv_pid_send (m->snd, listr, listv, listf, l);
  rcv_pid_receive (m->rcv, listr, listv, listf, l);
  rcv_pid_wait (m->snd);
  pfree (listr,__func__,__FILE__,__LINE__);
  pfree (listf,__func__,__FILE__,__LINE__);
  pfree (listv,__func__,__FILE__,__LINE__);
}

static void snd_rcv_destroy (SndRcv * m)
{
  rcv_pid_destroy (m->rcv);
  rcv_pid_destroy (m->snd);
}

static void snd_rcv_init (SndRcv * m, const char * name)
{
  char s[strlen(name) + 5];
  strcpy (s, name);
  strcat (s, ".rcv");
  m->rcv = rcv_pid_new (s);
  strcpy (s, name);
  strcat (s, ".snd");
  m->snd = rcv_pid_new (s);
}

static void mpi_boundary_destroy (Boundary * b)
{
  MpiBoundary * m = (MpiBoundary *) b;
  snd_rcv_destroy (&m->mpi_level);
  snd_rcv_destroy (&m->mpi_level_root);
  snd_rcv_destroy (&m->restriction);
  array_free (m->send);
  array_free (m->receive);
  pfree (m,__func__,__FILE__,__LINE__);
}


static void mpi_boundary_level (const Boundary * b, scalar * list, int l)
{ trace ("mpi_boundary_level", "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 493);
  MpiBoundary * m = (MpiBoundary *) b;
  rcv_pid_sync (&m->mpi_level, list, l);
  rcv_pid_sync (&m->mpi_level_root, list, l);
 end_trace("mpi_boundary_level", "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 497); }


static void mpi_boundary_restriction (const Boundary * b, scalar * list, int l)
{ trace ("mpi_boundary_restriction", "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 501);
  MpiBoundary * m = (MpiBoundary *) b;
  rcv_pid_sync (&m->restriction, list, l);
 end_trace("mpi_boundary_restriction", "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 504); }

void mpi_boundary_new ()
{
  mpi_boundary = (Boundary *) ((MpiBoundary *) pcalloc (1, sizeof(MpiBoundary),__func__,__FILE__,__LINE__));
  mpi_boundary->destroy = mpi_boundary_destroy;
  mpi_boundary->level = mpi_boundary_level;
  mpi_boundary->restriction = mpi_boundary_restriction;
  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;
  snd_rcv_init (&mpi->mpi_level, "mpi_level");
  snd_rcv_init (&mpi->mpi_level_root, "mpi_level_root");
  snd_rcv_init (&mpi->restriction, "restriction");
  mpi->send = array_new();
  mpi->receive = array_new();
  add_boundary (mpi_boundary);
}

static FILE * fopen_prefix (FILE * fp, const char * name, char * prefix)
{
  if (fp) {
    sprintf (prefix, "%s-%d ", name, pid());
    return fp;
  }
  else {
    strcpy (prefix, "");
    char fname[80];
    if (debug_iteration >= 0)
      sprintf (fname, "%s-%d-%d", name, debug_iteration, pid());
    else
      sprintf (fname, "%s-%d", name, pid());
    return fopen (fname, "w");
  }
}

void debug_mpi (FILE * fp1)
{
  void output_cells_internal (FILE * fp);

  char prefix[80];
  FILE * fp;


  if (fp1 == NULL) {
    char name[80];
    sprintf (name, "halo-%d", pid()); remove (name);
    sprintf (name, "cells-%d", pid()); remove (name);
    sprintf (name, "faces-%d", pid()); remove (name);
    sprintf (name, "vertices-%d", pid()); remove (name);
    sprintf (name, "neighbors-%d", pid()); remove (name);
    sprintf (name, "mpi-level-rcv-%d", pid()); remove (name);
    sprintf (name, "mpi-level-snd-%d", pid()); remove (name);
    sprintf (name, "mpi-level-root-rcv-%d", pid()); remove (name);
    sprintf (name, "mpi-level-root-snd-%d", pid()); remove (name);
    sprintf (name, "mpi-restriction-rcv-%d", pid()); remove (name);
    sprintf (name, "mpi-restriction-snd-%d", pid()); remove (name);
    sprintf (name, "mpi-border-%d", pid()); remove (name);
    sprintf (name, "exterior-%d", pid()); remove (name);
    sprintf (name, "depth-%d", pid()); remove (name);
    sprintf (name, "refined-%d", pid()); remove (name);
  }


  fp = fopen_prefix (fp1, "halo", prefix);
  for (int l = 0; l < depth(); l++)
     { foreach_halo (prolongation, l){

#line 568 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

       { foreach_child()
        fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level); end_foreach_child(); }; } end_foreach_halo(); }
  if (!fp1)
    fclose (fp);

  if (!fp1) {
    fp = fopen_prefix (fp1, "cells", prefix);
    output_cells_internal (fp);
    fclose (fp);
  }

  fp = fopen_prefix (fp1, "faces", prefix);
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  char _prefix[80];
 memcpy (_prefix, prefix, (80)*sizeof(char));
{ char * prefix = _prefix; NOT_UNUSED(prefix);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", .line = 581,
    .each = "foreach_face", .first = _first_call
  };
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 581
{

#line 581 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

    _stencil_fprintf (__FILE__,__LINE__,fp, "%s%g %g %g %d\n", prefix, x, y, z, level); }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 581
{

#line 581 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

    _stencil_fprintf (__FILE__,__LINE__,fp, "%s%g %g %g %d\n", prefix, x, y, z, level); }  }}  end_foreach_face_stencil()
#line 582
 if (_first_call) {
 for (int i = 0; i < (80); i++)
   if (_prefix[i] != prefix[i]) {
     reduction_warning ("/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 581, "prefix");
     break; }
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 582
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 581
{

#line 581 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level); }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 581
{

#line 581 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level); }  }}  end_foreach_face_generic()
#line 582
 end_foreach_face(); }
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "vertices", prefix);
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  char _prefix[80];
 memcpy (_prefix, prefix, (80)*sizeof(char));
{ char * prefix = _prefix; NOT_UNUSED(prefix);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", .line = 587,
    .each = "foreach_vertex", .first = _first_call
  };
foreach_vertex_stencil(){

#line 587 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

    _stencil_fprintf (__FILE__,__LINE__,fp, "%s%g %g %g %d\n", prefix, x, y, z, level); } end_foreach_vertex_stencil(); if (_first_call) {
 for (int i = 0; i < (80); i++)
   if (_prefix[i] != prefix[i]) {
     reduction_warning ("/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 587, "prefix");
     break; }
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 588
foreach_vertex(){

#line 587 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level); } end_foreach_vertex(); }
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "neighbors", prefix);
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  char _prefix[80];
 memcpy (_prefix, prefix, (80)*sizeof(char));
{ char * prefix = _prefix; NOT_UNUSED(prefix);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", .line = 593,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 593 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"
 {
    int n = 0;
     { foreach_neighbor(1)
      IF ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
 n++; end_foreach_neighbor(); }
    _stencil_fprintf (__FILE__,__LINE__,fp, "%s%g %g %g %d\n", prefix, x, y, z, cell.neighbors);
    IF (!(cell.neighbors == n)) _stencil_qassert (__FILE__,__LINE__,"/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 599, "cell.neighbors == n");
  } } end_foreach_stencil(); if (_first_call) {
 for (int i = 0; i < (80); i++)
   if (_prefix[i] != prefix[i]) {
     reduction_warning ("/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 593, "prefix");
     break; }
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 600
foreach(){

#line 593 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"
 {
    int n = 0;
     { foreach_neighbor(1)
      if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
 n++; end_foreach_neighbor(); }
    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, cell.neighbors);
    if (!(cell.neighbors == n)) qassert ("/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 599, "cell.neighbors == n");
  } } end_foreach(); }
  if (!fp1)
    fclose (fp);

  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;

  fp = fopen_prefix (fp1, "mpi-level-rcv", prefix);
  rcv_pid_print (mpi->mpi_level.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-level-root-rcv", prefix);
  rcv_pid_print (mpi->mpi_level_root.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-restriction-rcv", prefix);
  rcv_pid_print (mpi->restriction.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-level-snd", prefix);
  rcv_pid_print (mpi->mpi_level.snd, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-level-root-snd", prefix);
  rcv_pid_print (mpi->mpi_level_root.snd, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-restriction-snd", prefix);
  rcv_pid_print (mpi->restriction.snd, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-border", prefix);
   { foreach_cell(){

#line 637 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"
 {
    if (is_border(cell))
      fprintf (fp, "%s%g %g %g %d %d %d\n",
        prefix, x, y, z, level, cell.neighbors, cell.pid);
    else
      continue;
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "exterior", prefix);
   { foreach_cell(){

#line 650 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"
 {
    if (!is_local(cell))
      fprintf (fp, "%s%g %g %g %d %d %d %d\n",
        prefix, x, y, z, level, cell.neighbors,
        cell.pid, cell.flags & leaf);






  } } end_foreach_cell(); }
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "depth", prefix);
  fprintf (fp, "depth: %d %d\n", pid(), depth());
  fprintf (fp, "======= mpi_level.snd ======\n");
  RcvPid * snd = mpi->mpi_level.snd;
  for (int i = 0; i < snd->npid; i++)
    fprintf (fp, "%d %d %d\n", pid(), snd->rcv[i].pid, snd->rcv[i].maxdepth);
  fprintf (fp, "======= mpi_level.rcv ======\n");
  snd = mpi->mpi_level.rcv;
  for (int i = 0; i < snd->npid; i++)
    fprintf (fp, "%d %d %d\n", pid(), snd->rcv[i].pid, snd->rcv[i].maxdepth);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "refined", prefix);
   { foreach_cache (((Tree *)grid)->refined){

#line 679 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level); } end_foreach_cache(); }
  if (!fp1)
    fclose (fp);
}

static void snd_rcv_free (SndRcv * p)
{
  char name[strlen(p->rcv->name) + 1];
  strcpy (name, p->rcv->name);
  rcv_pid_destroy (p->rcv);
  p->rcv = rcv_pid_new (name);
  strcpy (name, p->snd->name);
  rcv_pid_destroy (p->snd);
  p->snd = rcv_pid_new (name);
}

static bool is_root (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 697 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

  if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
     { foreach_child()
      if (is_local(cell))
 return true; end_foreach_child(); }
  return false;

#if _call_is_root
}
#define _IN_STENCIL 1

#line 696
static bool _is_root (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 697 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

  IF ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
     { foreach_child()
      IF (is_local(cell))
 return true; end_foreach_child(); }
  return false;

#undef _IN_STENCIL

#endif

#line 703
}


static bool is_local_prolongation (Point point, Point p)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 707 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"


  struct { int x, y; } dp = {p.i - point.i, p.j - point.j};



  {
#line 713
 {
    if (dp.x == 0 && ((!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) || (!is_leaf (neighbor(1,0,0)) && neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0)))
      return true;
    if ((!is_leaf (neighbor(dp.x,0,0)) && neighbor(dp.x,0,0).neighbors && neighbor(dp.x,0,0).pid >= 0))
      return true;
  }
#line 713
 {
    if (dp.y == 0 && ((!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) || (!is_leaf (neighbor(0,1,0)) && neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0)))
      return true;
    if ((!is_leaf (neighbor(0,dp.y,0)) && neighbor(0,dp.y,0).neighbors && neighbor(0,dp.y,0).pid >= 0))
      return true;
  }}
  return false;

#if _call_is_local_prolongation
}
#define _IN_STENCIL 1

#line 706
static bool _is_local_prolongation (Point point, Point p)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 707 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"


  struct { int x, y; } dp = {p.i - point.i, p.j - point.j};



  {
#line 713
 {
    IF (dp.x == 0 && ((!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) || (!is_leaf (neighbor(1,0,0)) && neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0)))
      return true;
    IF ((!is_leaf (neighbor(dp.x,0,0)) && neighbor(dp.x,0,0).neighbors && neighbor(dp.x,0,0).pid >= 0))
      return true;
  }
#line 713
 {
    IF (dp.y == 0 && ((!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) || (!is_leaf (neighbor(0,1,0)) && neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0)))
      return true;
    IF ((!is_leaf (neighbor(0,dp.y,0)) && neighbor(0,dp.y,0).neighbors && neighbor(0,dp.y,0).pid >= 0))
      return true;
  }}
  return false;

#undef _IN_STENCIL

#endif

#line 720
}



static void append_pid (Array * pids, int pid)
{
  for (int i = 0, * p = (int *) pids->p; i < pids->len/sizeof(int); i++, p++)
    if (*p == pid)
      return;
  array_append (pids, &pid, sizeof(int));
}

static int locals_pids (Point point, Array * pids)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 733 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

  if (is_leaf(cell)) {
    if (is_local(cell)) {
      Point p = point;
       { foreach_neighbor(1) {
 if ((cell.pid >= 0 && cell.pid != pid()) &&
     ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) || is_local_prolongation (point, p)))
   append_pid (pids, cell.pid);
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
    { foreach_child()
     if ((cell.pid >= 0 && cell.pid != pid()))
       append_pid (pids, cell.pid); end_foreach_child(); }
      } end_foreach_neighbor(); }
    }
  }
  else
     { foreach_neighbor(1) {
      if ((cell.pid >= 0 && cell.pid != pid()))
 append_pid (pids, cell.pid);
      if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
  { foreach_child()
   if ((cell.pid >= 0 && cell.pid != pid()))
     append_pid (pids, cell.pid); end_foreach_child(); }
    } end_foreach_neighbor(); }
  return pids->len/sizeof(int);

#if _call_locals_pids
}
#define _IN_STENCIL 1
#define is_local_prolongation _is_local_prolongation

#line 732
static int _locals_pids (Point point, Array * pids)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 733 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

  IF (is_leaf(cell)) {
    IF (is_local(cell)) {
      Point p = point;
       { foreach_neighbor(1) {
 IF ((cell.pid >= 0 && cell.pid != pid()) &&
     ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) || is_local_prolongation (point, p)))
   append_pid (pids, cell.pid);
 IF ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
    { foreach_child()
     IF ((cell.pid >= 0 && cell.pid != pid()))
       append_pid (pids, cell.pid); end_foreach_child(); }
      } end_foreach_neighbor(); }
    }
  }
  
     { foreach_neighbor(1) {
      IF ((cell.pid >= 0 && cell.pid != pid()))
 append_pid (pids, cell.pid);
      IF ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
  { foreach_child()
   IF ((cell.pid >= 0 && cell.pid != pid()))
     append_pid (pids, cell.pid); end_foreach_child(); }
    } end_foreach_neighbor(); }
  return pids->len/sizeof(int);

#undef is_local_prolongation
#undef _IN_STENCIL

#endif

#line 758
}

static int root_pids (Point point, Array * pids)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 761 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

   { foreach_child()
    if ((cell.pid >= 0 && cell.pid != pid()))
      append_pid (pids, cell.pid); end_foreach_child(); }
  return pids->len/sizeof(int);

#if _call_root_pids
}
#define _IN_STENCIL 1

#line 760
static int _root_pids (Point point, Array * pids)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 761 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

   { foreach_child()
    IF ((cell.pid >= 0 && cell.pid != pid()))
      append_pid (pids, cell.pid); end_foreach_child(); }
  return pids->len/sizeof(int);

#undef _IN_STENCIL

#endif

#line 766
}







static void rcv_pid_row (RcvPid * m, int l, int * row)
{
  for (int i = 0; i < npe(); i++)
    row[i] = 0;
  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0)
      row[rcv->pid] = rcv->halo[l].n;
  }
}

void check_snd_rcv_matrix (SndRcv * sndrcv, const char * name)
{
  int maxlevel = depth();
  mpi_all_reduce (maxlevel, MPI_INT, MPI_MAX);
  int * row = ((int *) pmalloc ((npe())*sizeof(int),__func__,__FILE__,__LINE__));
  for (int l = 0; l <= maxlevel; l++) {
    int status = 0;
    if (pid() == 0) {


      int ** send = matrix_new (npe(), npe(), sizeof(int));
      int ** receive = matrix_new (npe(), npe(), sizeof(int));
      rcv_pid_row (sndrcv->snd, l, row);
      MPI_Gather (row, npe(), MPI_INT, &send[0][0], npe(), MPI_INT, 0,
    MPI_COMM_WORLD);
      rcv_pid_row (sndrcv->rcv, l, row);
      MPI_Gather (row, npe(), MPI_INT, &receive[0][0], npe(), MPI_INT, 0,
    MPI_COMM_WORLD);

      int * astatus = ((int *) pmalloc ((npe())*sizeof(int),__func__,__FILE__,__LINE__));
      for (int i = 0; i < npe(); i++)
 astatus[i] = 0;
      for (int i = 0; i < npe(); i++)
 for (int j = 0; j < npe(); j++)
   if (send[i][j] != receive[j][i]) {
     fprintf (ferr, "%s: %d sends    %d to   %d at level %d\n",
       name, i, send[i][j], j, l);
     fprintf (ferr, "%s: %d receives %d from %d at level %d\n",
       name, j, receive[j][i], i, l);
     fflush (ferr);
     for (int k = i - 2; k <= i + 2; k++)
       if (k >= 0 && k < npe())
  astatus[k] = 1;
     for (int k = j - 2; k <= j + 2; k++)
       if (k >= 0 && k < npe())
  astatus[k] = 1;
   }
      MPI_Scatter (astatus, 1, MPI_INT, &status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      pfree (astatus,__func__,__FILE__,__LINE__);

      matrix_free (send);
      matrix_free (receive);
    }
    else {
      rcv_pid_row (sndrcv->snd, l, row);
      MPI_Gather (row, npe(), MPI_INT, NULL, npe(), MPI_INT, 0, MPI_COMM_WORLD);
      rcv_pid_row (sndrcv->rcv, l, row);
      MPI_Gather (row, npe(), MPI_INT, NULL, npe(), MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Scatter (NULL, 1, MPI_INT, &status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    if (status) {
      fprintf (ferr,
        "check_snd_rcv_matrix \"%s\" failed\n"
        "Calling debug_mpi(NULL)...\n"
        "Aborting...\n",
        name);
      fflush (ferr);
      debug_mpi (NULL);
      MPI_Abort (MPI_COMM_WORLD, -3);
    }
  }
  pfree (row,__func__,__FILE__,__LINE__);
}

static bool has_local_child (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 850 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

   { foreach_child()
    if (is_local(cell))
      return true; end_foreach_child(); }
  return false;

#if _call_has_local_child
}
#define _IN_STENCIL 1

#line 849
static bool _has_local_child (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 850 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

   { foreach_child()
    IF (is_local(cell))
      return true; end_foreach_child(); }
  return false;

#undef _IN_STENCIL

#endif

#line 855
}


void mpi_boundary_update_buffers ()
{ trace ("mpi_boundary_update_buffers", "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 859);
  if (npe() == 1)
    { ; end_trace("mpi_boundary_update_buffers", "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 861);  return; }

  prof_start ("mpi_boundary_update_buffers");

  MpiBoundary * m = (MpiBoundary *) mpi_boundary;
  SndRcv * mpi_level = &m->mpi_level;
  SndRcv * mpi_level_root = &m->mpi_level_root;
  SndRcv * restriction = &m->restriction;

  snd_rcv_free (mpi_level);
  snd_rcv_free (mpi_level_root);
  snd_rcv_free (restriction);

  static const unsigned short used = 1 << user;
   { foreach_cell(){

#line 875 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"
 {
    if (is_active(cell) && !is_border(cell))



      continue;

    if (cell.neighbors) {

      Array pids = {NULL, 0, 0};
      int n = locals_pids (point, &pids);
      if (n) {
  { foreach_child()
   if (is_local(cell))
     for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
       rcv_pid_append (mpi_level->snd, *p, point); end_foreach_child(); }
 pfree (pids.p,__func__,__FILE__,__LINE__);
      }

      bool locals = false;
      if (is_leaf(cell)) {
 if ((cell.pid >= 0 && cell.pid != pid())) {
   Point p = point;
    { foreach_neighbor(1)
     if ((is_local(cell) &&
   ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) || is_local_prolongation (point, p))) ||
  is_root(point))
       locals = true, foreach_neighbor_break(); end_foreach_neighbor(); }
 }
      }
      else
  { foreach_neighbor(1)
   if (is_local(cell) || is_root(point))
     locals = true, foreach_neighbor_break(); end_foreach_neighbor(); }
      if (locals)
  { foreach_child()
   if ((cell.pid >= 0 && cell.pid != pid()))
            rcv_pid_append (mpi_level->rcv, cell.pid, point),
       cell.flags |= used; end_foreach_child(); }


      if (!is_leaf(cell)) {

 if (is_local(cell)) {
   Array pids = {NULL, 0, 0};

   int n = root_pids (point, &pids);
   if (n) {
      { foreach_neighbor()
       for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
  if (cell.pid >= 0 && cell.pid != *p)
    rcv_pid_append (mpi_level_root->snd, *p, point); end_foreach_neighbor(); }

     for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
       rcv_pid_append (restriction->snd, *p, point);
     pfree (pids.p,__func__,__FILE__,__LINE__);
   }
 }

 else if ((cell.pid >= 0 && cell.pid != pid())) {
   bool root = false;
    { foreach_child()
     if (is_local(cell))
       root = true, foreach_child_break(); end_foreach_child(); }
   if (root) {
     int pid = cell.pid;
      { foreach_neighbor()
       if ((cell.pid >= 0 && cell.pid != pid()))
  rcv_pid_append (mpi_level_root->rcv, pid, point),
    cell.flags |= used; end_foreach_neighbor(); }

     rcv_pid_append (restriction->rcv, pid, point);
   }
 }
      }
    }


    if (level > 0) {
      if (is_local(cell)) {

 Array pids = {NULL, 0, 0};
 if ((aparent(0,0,0).pid >= 0 && aparent(0,0,0).pid != pid()))
   append_pid (&pids, aparent(0,0,0).pid);
 int n = root_pids (parent, &pids);
 if (n) {
   for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
     rcv_pid_append (restriction->snd, *p, point);
   pfree (pids.p,__func__,__FILE__,__LINE__);
 }
      }
      else if ((cell.pid >= 0 && cell.pid != pid())) {

 if (is_local(aparent(0,0,0)) || has_local_child (parent))
   rcv_pid_append (restriction->rcv, cell.pid, point);
      }
    }
  } } end_foreach_cell(); }





  static const unsigned short keep = 1 << (user + 1);
  for (int l = depth(); l >= 0; l--)
     { foreach_cell(){

#line 980 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

      if (level == l) {
 if (level > 0 && (cell.pid < 0 || is_local(cell) || (cell.flags & used)))
   aparent(0,0,0).flags |= keep;
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) && !(cell.flags & keep))
   coarsen_cell (point, NULL);
 cell.flags &= ~(used|keep);
 continue;
      } } end_foreach_cell(); }


  m->send->len = m->receive->len = 0;
  rcv_pid_append_pids (mpi_level->snd, m->send);
  rcv_pid_append_pids (mpi_level_root->snd, m->send);
  rcv_pid_append_pids (mpi_level->rcv, m->receive);
  rcv_pid_append_pids (mpi_level_root->rcv, m->receive);

  prof_stop();
#line 1012 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"
 end_trace("mpi_boundary_update_buffers", "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 1012); }


void mpi_boundary_refine (scalar * list)
{ trace ("mpi_boundary_refine", "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 1016);
  prof_start ("mpi_boundary_refine");

  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;


  Array * snd = mpi->send;
  MPI_Request r[2*snd->len/sizeof(int)];
  int nr = 0;
  for (int i = 0, * dest = snd->p; i < snd->len/sizeof(int); i++,dest++) {
    int len = ((Tree *)grid)->refined.n;
    MPI_Isend (&((Tree *)grid)->refined.n, 1, MPI_INT, *dest,
        (128), MPI_COMM_WORLD, &r[nr++]);
    if (len > 0)
      MPI_Isend (((Tree *)grid)->refined.p, sizeof(Index)/sizeof(int)*len,
   MPI_INT, *dest, (128), MPI_COMM_WORLD, &r[nr++]);
  }



  Array * rcv = mpi->receive;
  Cache rerefined = {NULL, 0, 0};
  for (int i = 0, * source = rcv->p; i < rcv->len/sizeof(int); i++,source++) {
    int len;
    mpi_recv_check (&len, 1, MPI_INT, *source, (128),
      MPI_COMM_WORLD, MPI_STATUS_IGNORE,
      "mpi_boundary_refine (len)");
    if (len > 0) {
      Index p[len];
      mpi_recv_check (p, sizeof(Index)/sizeof(int)*len,
        MPI_INT, *source, (128),
        MPI_COMM_WORLD, MPI_STATUS_IGNORE,
        "mpi_boundary_refine (p)");
      Cache refined = {p, len, len};
       { foreach_cache (refined){

#line 1050 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

 if (level <= depth() && allocated(0,0,0)) {
   if (is_leaf(cell)) {
     bool neighbors = false;
      { foreach_neighbor()
       if (allocated(0,0,0) && (is_active(cell) || is_local(aparent(0,0,0))))
  neighbors = true, foreach_neighbor_break(); end_foreach_neighbor(); }

     if (neighbors)
       refine_cell (point, list, 0, &rerefined);
   }
 } } end_foreach_cache(); }
    }
  }


  if (nr)
    MPI_Waitall (nr, r, MPI_STATUSES_IGNORE);


  pfree (((Tree *)grid)->refined.p,__func__,__FILE__,__LINE__);
  ((Tree *)grid)->refined = rerefined;

  prof_stop();



  mpi_all_reduce (rerefined.n, MPI_INT, MPI_SUM);
  if (rerefined.n)
    mpi_boundary_refine (list);
  strongif (list) for (scalar s = *list, *_i81 = list; ((scalar *)&s)->i >= 0; s = *++_i81)
    _attribute[s.i].dirty = true;
 end_trace("mpi_boundary_refine", "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 1082); }

static void check_depth ()
{
#line 1117 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"
}

typedef struct {
  int refined, leaf;
} Remote;




void mpi_boundary_coarsen (int l, int too_fine)
{ trace ("mpi_boundary_coarsen", "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 1127);
  if (npe() == 1)
    { ; end_trace("mpi_boundary_coarsen", "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 1129);  return; }

  check_depth();

  if (!(sizeof(Remote) == sizeof(double))) qassert ("/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 1133, "sizeof(Remote) == sizeof(double)");

  scalar remote= new_scalar("remote");
   { foreach_cell(){

#line 1136 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"
 {
    if (level == l) {
      if (is_local(cell)) {
 ((Remote *)&val(remote,0,0,0))->refined = (!is_leaf (cell) && cell.neighbors && cell.pid >= 0);
 ((Remote *)&val(remote,0,0,0))->leaf = is_leaf(cell);
      }
      else {
 ((Remote *)&val(remote,0,0,0))->refined = true;
 ((Remote *)&val(remote,0,0,0))->leaf = false;
      }
      continue;
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }
  mpi_boundary_level (mpi_boundary, ((scalar []){remote,{-1}}), l);

   { foreach_cell(){

#line 1153 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"
 {
    if (level == l) {
      if (!is_local(cell)) {
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) && !((Remote *)&val(remote,0,0,0))->refined)
   coarsen_cell_recursive (point, NULL);
 else if (is_leaf(cell) && cell.neighbors && ((Remote *)&val(remote,0,0,0))->leaf) {
   int pid = cell.pid;
    { foreach_child()
     cell.pid = pid; end_foreach_child(); }
 }
      }
      continue;
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

  check_depth();

  if (l > 0) {
     { foreach_cell(){

#line 1173 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"
 {
      if (level == l) {
 val(remote,0,0,0) = is_local(cell) ? cell.neighbors : 0;
 continue;
      }
      if (is_leaf(cell))
 continue;
    } } end_foreach_cell(); }
    mpi_boundary_level (mpi_boundary, ((scalar []){remote,{-1}}), l);
     { foreach_cell(){

#line 1182 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"
 {
      if (level == l)
 if (!is_local(cell) && is_local(aparent(0,0,0)) && val(remote,0,0,0)) {
   aparent(0,0,0).flags &= ~too_fine;
   continue;
 }
      if (is_leaf(cell))
 continue;
    } } end_foreach_cell(); }
  }
 delete (((scalar []){remote,{-1}}));  end_trace("mpi_boundary_coarsen", "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 1192); }

static void flag_border_cells ()
{
   { foreach_cell(){

#line 1196 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"
 {
    if (is_active(cell)) {
      short flags = cell.flags & ~border;
       { foreach_neighbor() {
 if (!is_local(cell) || (level > 0 && !is_local(aparent(0,0,0))))
   flags |= border, foreach_neighbor_break();

 if (is_refined_check())
    { foreach_child()
     if (!is_local(cell))
       flags |= border, foreach_child_break(); end_foreach_child(); }
 if (flags & border)
   foreach_neighbor_break();
      } end_foreach_neighbor(); }
      cell.flags = flags;
    }
    else {
      cell.flags &= ~border;

    }
    if (is_leaf(cell)) {
      if (cell.neighbors) {
  { foreach_child()
   cell.flags &= ~border; end_foreach_child(); }
 if (is_border(cell)) {
   bool remote = false;
    { foreach_neighbor (2/2)
     if (!is_local(cell))
       remote = true, foreach_neighbor_break(); end_foreach_neighbor(); }
   if (remote)
      { foreach_child()
       cell.flags |= border; end_foreach_child(); }
 }
      }
      continue;
    }
  } } end_foreach_cell(); }
}

static int balanced_pid (long index, long nt, int nproc)
{
  long ne = max(1, nt/nproc), nr = nt % nproc;
  int pid = index < nr*(ne + 1) ?
    index/(ne + 1) :
    nr + (index - nr*(ne + 1))/ne;
  return min(nproc - 1, pid);
}



void mpi_partitioning ()
{ trace ("mpi_partitioning", "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 1247);
  prof_start ("mpi_partitioning");

  long nt = 0;
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  long _nt = nt;
{ long nt = _nt; NOT_UNUSED(nt);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", .line = 1251,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 1251 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

    nt++; } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 1252

#if _OPENMP
  #undef OMP
  #define OMP(x)
#endif
#line 1251
foreach (){

#line 1251 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

    nt++; } end_foreach();
#if _OPENMP
  #undef OMP
  #define OMP(x) _Pragma(#x)
#endif
#line 1252
 }


  long i = 0;
  ((Tree *)grid)->dirty = true;
   { foreach_cell_post (is_active (cell)){

#line 1257 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

    if (is_active (cell)) {
      if (is_leaf (cell)) {
 cell.pid = balanced_pid (i++, nt, npe());
 if (cell.neighbors > 0) {
   int pid = cell.pid;
    { foreach_child()
     cell.pid = pid; end_foreach_child(); }
 }
 if (!is_local(cell))
   cell.flags &= ~active;
      }
      else {
 cell.pid = child(0,0,0).pid;
 bool inactive = true;
  { foreach_child()
   if (is_active(cell))
     inactive = false, foreach_child_break(); end_foreach_child(); }
 if (inactive)
   cell.flags &= ~active;
      }
    } } end_foreach_cell_post(); }

  flag_border_cells();

  prof_stop();

  mpi_boundary_update_buffers();
 end_trace("mpi_partitioning", "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 1285); }

void restore_mpi (FILE * fp, scalar * list1)
{
  long index = 0, nt = 0, start = ftell (fp);
  scalar size= new_scalar("size"), * list = list_concat (((scalar []){size,{-1}}), list1);;
  long offset = sizeof(double)*list_len(list);


  static const unsigned short set = 1 << user;
  scalar * listm = is_constant(cm) ? NULL : (scalar *)((vector []){{fm.x,fm.y},{{-1},{-1}}});
   { foreach_cell(){

#line 1296 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

    if (balanced_pid (index, nt, npe()) <= pid()) {
      unsigned flags;
      if (fread (&flags, sizeof(unsigned), 1, fp) != 1) {
 fprintf (ferr, "restore(): error: expecting 'flags'\n");
 exit (1);
      }
      strongif (list) for (scalar s = *list, *_i82 = list; ((scalar *)&s)->i >= 0; s = *++_i82) {
 double val;
 if (fread (&val, sizeof(double), 1, fp) != 1) {
   fprintf (ferr, "restore(): error: expecting scalar\n");
   exit (1);
 }
 if (s.i != INT_MAX)
   val(s,0,0,0) = val;
      }
      if (level == 0)
 nt = val(size,0,0,0);
      cell.pid = balanced_pid (index, nt, npe());
      cell.flags |= set;
      if (!(flags & leaf) && is_leaf(cell)) {
 if (balanced_pid (index + val(size,0,0,0) - 1, nt, npe()) < pid()) {
   fseek (fp, (sizeof(unsigned) + offset)*(val(size,0,0,0) - 1), SEEK_CUR);
   index += val(size,0,0,0);
   continue;
 }
 refine_cell (point, listm, 0, NULL);
      }
      index++;
      if (is_leaf(cell))
 continue;
    } } end_foreach_cell(); }


  fseek (fp, start, SEEK_SET);
  index = 0;
   { foreach_cell(){

#line 1332 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"
 {
    unsigned flags;
    if (fread (&flags, sizeof(unsigned), 1, fp) != 1) {
      fprintf (ferr, "restore(): error: expecting 'flags'\n");
      exit (1);
    }
    if (cell.flags & set)
      fseek (fp, offset, SEEK_CUR);
    else {
      strongif (list) for (scalar s = *list, *_i83 = list; ((scalar *)&s)->i >= 0; s = *++_i83) {
 double val;
 if (fread (&val, sizeof(double), 1, fp) != 1) {
   fprintf (ferr, "restore(): error: expecting a scalar\n");
   exit (1);
 }
 if (s.i != INT_MAX)
   val(s,0,0,0) = val;
      }
      cell.pid = balanced_pid (index, nt, npe());
      if (is_leaf(cell) && cell.neighbors) {
 int pid = cell.pid;
  { foreach_child()
   cell.pid = pid; end_foreach_child(); }
      }
    }
    if (!(flags & leaf) && is_leaf(cell)) {
      bool locals = false;
       { foreach_neighbor(1)
 if ((cell.flags & set) && (is_local(cell) || is_root(point)))
   locals = true, foreach_neighbor_break(); end_foreach_neighbor(); }
      if (locals)
 refine_cell (point, listm, 0, NULL);
      else {
 fseek (fp, (sizeof(unsigned) + offset)*(val(size,0,0,0) - 1), SEEK_CUR);
 index += val(size,0,0,0);
 continue;
      }
    }
    index++;
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }


   { foreach_cell_post (is_active (cell)){

#line 1376 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"
 {
    cell.flags &= ~set;
    if (is_active (cell)) {
      if (is_leaf (cell)) {
 if (cell.neighbors > 0) {
   int pid = cell.pid;
    { foreach_child()
     cell.pid = pid; end_foreach_child(); }
 }
 if (!is_local(cell))
   cell.flags &= ~active;
      }
      else if (!is_local(cell)) {
 bool inactive = true;
  { foreach_child()
   if (is_active(cell))
     inactive = false, foreach_child_break(); end_foreach_child(); }
 if (inactive)
   cell.flags &= ~active;
      }
    }
  } } end_foreach_cell_post(); }

  flag_border_cells();

  mpi_boundary_update (list);
  pfree (list,__func__,__FILE__,__LINE__);
 delete (((scalar []){size,{-1}})); }
#line 1425 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

double z_indexing (scalar index, bool leaves)
{ trace ("z_indexing", "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 1427);



  scalar size= new_scalar("size");
  subtree_size (size, leaves);






  double maxi = -1.;
  if (pid() == 0)
     { foreach_level(0){

#line 1441 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

      maxi = val(size,0,0,0) - 1.; } end_foreach_level(); }




   { foreach_level(0){

#line 1447 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"

    val(index,0,0,0) = 0; } end_foreach_level(); }
  for (int l = 0; l < depth(); l++) {
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar []){index,{-1}}), l); };
     { foreach_cell(){

#line 1451 "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h"
 {
      if (level == l) {
 if (is_leaf(cell)) {
   if (is_local(cell) && cell.neighbors) {
     int i = val(index,0,0,0);
      { foreach_child()
       val(index,0,0,0) = i; end_foreach_child(); }
   }
 }
 else {
   bool loc = is_local(cell);
   if (!loc)
      { foreach_child()
       if (is_local(cell))
  loc = true, foreach_child_break(); end_foreach_child(); }
   if (loc) {
     int i = val(index,0,0,0) + !leaves;
      { foreach_child() {
       val(index,0,0,0) = i;
       i += val(size,0,0,0);
     } end_foreach_child(); }
   }
 }
 continue;
      }
      if (is_leaf(cell))
 continue;
    } } end_foreach_cell(); }
  }
  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar []){index,{-1}}), depth()); };

  { double _ret =  maxi; delete (((scalar []){size,{-1}}));  end_trace("z_indexing", "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 1482);  return _ret; }
 delete (((scalar []){size,{-1}}));  end_trace("z_indexing", "/home/fpl/softwares/basilisk/src/grid/tree-mpi.h", 1483); }
#line 1678 "/home/fpl/softwares/basilisk/src/grid/tree.h"
#line 1 "grid/balance.h"
#line 1 "/home/fpl/softwares/basilisk/src/grid/balance.h"


typedef struct {
  short leaf, prolongation;
  int pid;
} NewPid;



#if TRASH
# define is_newpid() (!isnan(val(newpid,0,0,0)) && ((NewPid *)&val(newpid,0,0,0))->pid > 0)
#else
# define is_newpid() (((NewPid *)&val(newpid,0,0,0))->pid > 0)
#endif

Array * linear_tree (size_t size, scalar newpid)
{
  const unsigned short sent = 1 << user, next = 1 << (user + 1);
  Array * a = array_new();

   { foreach_cell_post_all (true){

#line 21 "/home/fpl/softwares/basilisk/src/grid/balance.h"

    if (level > 0 && (cell.flags & (sent|next)))
      aparent(0,0,0).flags |= next; } end_foreach_cell_post_all(); }

  bool empty = true;
   { foreach_cell_all(){

#line 26 "/home/fpl/softwares/basilisk/src/grid/balance.h"
 {
    if (cell.flags & sent) {
      array_append (a, &cell, size);
      cell.flags &= ~sent;
      empty = false;
    }
    else {
      if (cell.pid >= 0 && ((NewPid *)&val(newpid,0,0,0))->leaf)
 if (!(is_leaf(cell))) qassert ("/home/fpl/softwares/basilisk/src/grid/balance.h", 34, "is_leaf(cell)");
      if (is_refined_check()) {


 bool prolo = false;
  { foreach_child()
   if (((NewPid *)&val(newpid,0,0,0))->prolongation)
     prolo = true; end_foreach_child(); }
 if (prolo) {

   cell.flags |= leaf;
   array_append (a, &cell, sizeof(Cell));
   cell.flags &= ~leaf;
 }
 else
   array_append (a, &cell, sizeof(Cell));
      }
      else
 array_append (a, &cell, sizeof(Cell));
    }
    if (cell.flags & next)
      cell.flags &= ~next;
    else
      continue;
  } } end_foreach_cell_all(); }

  if (empty)
    a->len = 0;
  return a;
}

#define foreach_tree(t, size, list)\
{\
  const unsigned short _sent = 1 << user, _next = 1 << (user + 1);\
  scalar * _list = list;\
  char * _i = (char *) (t)->p;\
  foreach_cell_all() {\
    Cell * c = (Cell *) _i;\
    if (c->flags & _sent) {\
      _i += size;\

#line 74


#define end_foreach_tree()\
    }\
    else\
      _i += sizeof(Cell);\
    if (c->flags & _next) {\
      if (!(c->neighbors)) qassert ("/home/fpl/softwares/basilisk/src/grid/balance.h", 81, "c->neighbors");\
      if (!(c->flags & leaf) && is_leaf(cell) &&\
   (!is_newpid() || !((NewPid *)&val(newpid,0,0,0))->leaf))\
\
 refine_cell (point, _list, 0, NULL);\
      else if (!cell.neighbors)\
\
 alloc_children (point);\
    }\
    else\
      continue;\
  } end_foreach_cell_all();\
}\

#line 94


Array * neighborhood (scalar newpid, int nextpid, FILE * fp)
{
  const unsigned short sent = 1 << user;
   { foreach_cell(){

#line 99 "/home/fpl/softwares/basilisk/src/grid/balance.h"
 {

    bool root = false;
    if ((!is_local(cell) || ((NewPid *)&val(newpid,0,0,0))->pid - 1 != nextpid) && (!is_leaf (cell) && cell.neighbors && cell.pid >= 0)) {
       { foreach_child()
 if (is_local(cell) && ((NewPid *)&val(newpid,0,0,0))->pid - 1 == nextpid)
   root = true, foreach_child_break(); end_foreach_child(); }
      if (root && cell.pid != nextpid) {
  { foreach_neighbor()
   if (cell.pid != nextpid && is_newpid()) {
     if (fp)
       fprintf (fp, "%g %g %g %d %d root\n",
         x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
     cell.flags |= sent;
   } end_foreach_neighbor(); }
      }
    }

    if ((is_local(cell) && ((NewPid *)&val(newpid,0,0,0))->pid - 1 == nextpid) || root) {
       { foreach_neighbor(1)
 if (cell.neighbors && cell.pid != nextpid)
    { foreach_child()
     if (cell.pid != nextpid && is_newpid()) {
       if (fp)
  fprintf (fp, "%g %g %g %d %d nextpid\n",
    x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
       cell.flags |= sent;
     } end_foreach_child(); } end_foreach_neighbor(); }
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

  return linear_tree (sizeof(Cell) + datasize, newpid);
}

static void send_tree (Array * a, int to, MPI_Request * r)
{
  MPI_Isend (&a->len, 1, MPI_LONG, to, (256), MPI_COMM_WORLD, &r[0]);
  if (a->len > 0) {
    MPI_Isend (a->p, a->len, MPI_BYTE, to, (256), MPI_COMM_WORLD, &r[1]);
    ((Tree *)grid)->dirty = true;
  }
}

static void receive_tree (int from, scalar newpid, FILE * fp)
{
  Array a;
  mpi_recv_check (&a.len, 1, MPI_LONG, from, (256),
    MPI_COMM_WORLD, MPI_STATUS_IGNORE, "receive_tree (len)");
  if (a.len > 0) {
    a.p = pmalloc (a.len,__func__,__FILE__,__LINE__);
    if (fp)
      fprintf (fp, "receiving %ld from %d\n", a.len, from);
    mpi_recv_check (a.p, a.len, MPI_BYTE, from, (256),
      MPI_COMM_WORLD, MPI_STATUS_IGNORE, "receive_tree (p)");

     { foreach_tree (&a, sizeof(Cell) + datasize, NULL){

#line 156 "/home/fpl/softwares/basilisk/src/grid/balance.h"
 {
      memcpy (((char *)&cell) + sizeof(Cell), ((char *)c) + sizeof(Cell),
       datasize);
      if (!(((NewPid *)&val(newpid,0,0,0))->pid > 0)) qassert ("/home/fpl/softwares/basilisk/src/grid/balance.h", 159, "NEWPID()->pid > 0");
      if (fp)
 fprintf (fp, "%g %g %g %d %d %d %d %d %d recv\n",
   x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid,
   c->flags & leaf,
   cell.flags & leaf, from, ((NewPid *)&val(newpid,0,0,0))->leaf);
    } } end_foreach_tree(); }
    pfree (a.p,__func__,__FILE__,__LINE__);
    ((Tree *)grid)->dirty = true;
  }
}

static void wait_tree (Array * a, MPI_Request * r)
{
  MPI_Wait (&r[0], MPI_STATUS_IGNORE);
  if (a->len > 0)
    MPI_Wait (&r[1], MPI_STATUS_IGNORE);
}

static void check_flags ()
{







}

struct {
  int min;
  bool leaves;

  int npe;
} mpi = {
  1,
  true
};


bool balance ()
{ trace ("balance", "/home/fpl/softwares/basilisk/src/grid/balance.h", 201);
  if (npe() == 1)
    { bool _ret =  false; end_trace("balance", "/home/fpl/softwares/basilisk/src/grid/balance.h", 203);  return _ret; }

  if (!(sizeof(NewPid) == sizeof(double))) qassert ("/home/fpl/softwares/basilisk/src/grid/balance.h", 205, "sizeof(NewPid) == sizeof(double)");

  check_flags();

  long nl = 0, nt = 0;
   { foreach_cell(){

#line 210 "/home/fpl/softwares/basilisk/src/grid/balance.h"
 {
    if (is_local(cell)) {
      nt++;
      if (is_leaf(cell))
 nl++;
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

  grid->n = grid->tn = nl;
  grid->maxdepth = depth();
  long nmin = nl, nmax = nl;

  mpi_all_reduce (nmax, MPI_LONG, MPI_MAX);
  mpi_all_reduce (nmin, MPI_LONG, MPI_MIN);
  mpi_all_reduce (grid->tn, MPI_LONG, MPI_SUM);
  mpi_all_reduce (grid->maxdepth, MPI_INT, MPI_MAX);
  if (mpi.leaves)
    nt = grid->tn;
  else
    mpi_all_reduce (nt, MPI_LONG, MPI_SUM);

  long ne = max(1, nt/npe());

  if (ne < mpi.min) {
    mpi.npe = max(1, nt/mpi.min);
    ne = max(1, nt/mpi.npe);
  }
  else
    mpi.npe = npe();

  if (nmax - nmin <= 1)
    { bool _ret =  false; end_trace("balance", "/home/fpl/softwares/basilisk/src/grid/balance.h", 243);  return _ret; }

  scalar newpid= new_scalar("newpid");
  double zn = z_indexing (newpid, mpi.leaves);
  if (pid() == 0)
    if (!(zn + 1 == nt)) qassert ("/home/fpl/softwares/basilisk/src/grid/balance.h", 248, "zn + 1 == nt");

  FILE * fp = NULL;
#line 260 "/home/fpl/softwares/basilisk/src/grid/balance.h"
  bool next = false, prev = false;
   { foreach_cell_all(){

#line 261 "/home/fpl/softwares/basilisk/src/grid/balance.h"
 {
    if (is_local(cell)) {
      int pid = balanced_pid (val(newpid,0,0,0), nt, mpi.npe);
      pid = clamp (pid, cell.pid - 1, cell.pid + 1);
      if (pid == pid() + 1)
 next = true;
      else if (pid == pid() - 1)
 prev = true;
      ((NewPid *)&val(newpid,0,0,0))->pid = pid + 1;
      ((NewPid *)&val(newpid,0,0,0))->leaf = is_leaf(cell);
      ((NewPid *)&val(newpid,0,0,0))->prolongation = (!is_leaf(cell) && !cell.neighbors && cell.pid >= 0);
      if (fp)
 fprintf (fp, "%g %g %d %d newpid\n", x, y, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
    }
    else
      val(newpid,0,0,0) = 0;
  } } end_foreach_cell_all(); }
  for (int l = 0; l <= depth(); l++)
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, ((scalar []){newpid,{-1}}), l); };
#line 304 "/home/fpl/softwares/basilisk/src/grid/balance.h"
  Array * anext = next ? neighborhood (newpid, pid() + 1, fp) : array_new();
  Array * aprev = prev ? neighborhood (newpid, pid() - 1, fp) : array_new();

  if (fp)
    fflush (fp);

  check_flags();


  MPI_Request rprev[2], rnext[2];
  if (pid() > 0)
    send_tree (aprev, pid() - 1, rprev);
  if (pid() < npe() - 1)
    send_tree (anext, pid() + 1, rnext);


  if (pid() < npe() - 1)
    receive_tree (pid() + 1, newpid, fp);
  if (pid() > 0)
    receive_tree (pid() - 1, newpid, fp);


  if (pid() > 0)
    wait_tree (aprev, rprev);
  array_free (aprev);
  if (pid() < npe() - 1)
    wait_tree (anext, rnext);
  array_free (anext);

  if (fp)
    fflush (fp);


  int pid_changed = false;
   { foreach_cell_all(){

#line 338 "/home/fpl/softwares/basilisk/src/grid/balance.h"
 {
    if (cell.pid >= 0) {
      if (is_newpid()) {
 if (fp)
   fprintf (fp, "%g %g %g %d %d %d %d %d new\n",
     x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid,
     is_leaf(cell), cell.neighbors, ((NewPid *)&val(newpid,0,0,0))->leaf);
 if (cell.pid != ((NewPid *)&val(newpid,0,0,0))->pid - 1) {
   cell.pid = ((NewPid *)&val(newpid,0,0,0))->pid - 1;
   cell.flags &= ~(active|border);
   if (is_local(cell))
     cell.flags |= active;
   pid_changed = true;
 }
 if (((NewPid *)&val(newpid,0,0,0))->leaf && !is_leaf(cell) && cell.neighbors)
   coarsen_cell_recursive (point, NULL);
      }
      else if (level > 0 && ((NewPid *)&coarse(newpid,0,0,0))->leaf)
 cell.pid = aparent(0,0,0).pid;
    }

    if (!cell.neighbors && allocated_child(0,0,0)) {
      if (fp)
 fprintf (fp, "%g %g %g %d %d freechildren\n",
   x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
      free_children (point);
    }
  } } end_foreach_cell_all(); }

  if (((Tree *)grid)->dirty || pid_changed) {


     { foreach_cell_post (!is_leaf (cell)){

#line 370 "/home/fpl/softwares/basilisk/src/grid/balance.h"

      if (!is_leaf(cell) && !is_local(cell)) {
 unsigned short flags = cell.flags & ~active;
  { foreach_child()
   if (is_active(cell))
     flags |= active, foreach_child_break(); end_foreach_child(); }
 cell.flags = flags;
      } } end_foreach_cell_post(); }

    flag_border_cells();
    pid_changed = true;
  }

  if (fp)
    fclose (fp);

  mpi_all_reduce (pid_changed, MPI_INT, MPI_MAX);
  if (pid_changed)
    mpi_boundary_update_buffers();

  { bool _ret =  pid_changed; delete (((scalar []){newpid,{-1}}));  end_trace("balance", "/home/fpl/softwares/basilisk/src/grid/balance.h", 390);  return _ret; }
 delete (((scalar []){newpid,{-1}}));  end_trace("balance", "/home/fpl/softwares/basilisk/src/grid/balance.h", 391); }

void mpi_boundary_update (scalar * list)
{
  mpi_boundary_update_buffers();
  strongif (list) for (scalar s = *list, *_i84 = list; ((scalar *)&s)->i >= 0; s = *++_i84)
    _attribute[s.i].dirty = true;
  grid->tn = 0;
  boundary_internal ((scalar *)(list), "/home/fpl/softwares/basilisk/src/grid/balance.h", 399);
  while (balance());
}
#line 1679 "/home/fpl/softwares/basilisk/src/grid/tree.h"
#else
void mpi_boundary_refine (scalar * list){}
void mpi_boundary_coarsen (int a, int b){}
void mpi_boundary_update (scalar * list) {
  strongif (list) for (scalar s = *list, *_i85 = list; ((scalar *)&s)->i >= 0; s = *++_i85)
    _attribute[s.i].dirty = true;
  boundary_internal ((scalar *)(list), "/home/fpl/softwares/basilisk/src/grid/tree.h", 1685);
}
#endif
#line 4 "/home/fpl/softwares/basilisk/src/grid/quadtree.h"

void quadtree_methods () {
  tree_methods();
}
#line 15 "60d_plate_advancing_simulation-cpp.c"
static double _boundary0 (Point point, Point neighbor, scalar _s, void * data);
static double _boundary0_homogeneous (Point point, Point neighbor, scalar _s, void * data);
static double _boundary1 (Point point, Point neighbor, scalar _s, void * data);
static double _boundary1_homogeneous (Point point, Point neighbor, scalar _s, void * data);
static double _boundary2 (Point point, Point neighbor, scalar _s, void * data);
static double _boundary2_homogeneous (Point point, Point neighbor, scalar _s, void * data);
static double _boundary3 (Point point, Point neighbor, scalar _s, void * data);
static double _boundary3_homogeneous (Point point, Point neighbor, scalar _s, void * data);
static double _boundary4 (Point point, Point neighbor, scalar _s, void * data);
static double _boundary4_homogeneous (Point point, Point neighbor, scalar _s, void * data);
static double _boundary5 (Point point, Point neighbor, scalar _s, void * data);
static double _boundary5_homogeneous (Point point, Point neighbor, scalar _s, void * data);
static double _boundary6 (Point point, Point neighbor, scalar _s, void * data);
static double _boundary6_homogeneous (Point point, Point neighbor, scalar _s, void * data);
static double _boundary7 (Point point, Point neighbor, scalar _s, void * data);
static double _boundary7_homogeneous (Point point, Point neighbor, scalar _s, void * data);
static double _boundary8 (Point point, Point neighbor, scalar _s, void * data);
static double _boundary8_homogeneous (Point point, Point neighbor, scalar _s, void * data);
static double _boundary9 (Point point, Point neighbor, scalar _s, void * data);
static double _boundary9_homogeneous (Point point, Point neighbor, scalar _s, void * data);
static double _boundary10 (Point point, Point neighbor, scalar _s, void * data);
static double _boundary10_homogeneous (Point point, Point neighbor, scalar _s, void * data);
static double _boundary11 (Point point, Point neighbor, scalar _s, void * data);
static double _boundary11_homogeneous (Point point, Point neighbor, scalar _s, void * data);
static double _boundary12 (Point point, Point neighbor, scalar _s, void * data);
static double _boundary12_homogeneous (Point point, Point neighbor, scalar _s, void * data);
static double _boundary13 (Point point, Point neighbor, scalar _s, void * data);
static double _boundary13_homogeneous (Point point, Point neighbor, scalar _s, void * data);
static double _boundary14 (Point point, Point neighbor, scalar _s, void * data);
static double _boundary14_homogeneous (Point point, Point neighbor, scalar _s, void * data);
static double _boundary15 (Point point, Point neighbor, scalar _s, void * data);
static double _boundary15_homogeneous (Point point, Point neighbor, scalar _s, void * data);
static double _boundary16 (Point point, Point neighbor, scalar _s, void * data);
static double _boundary16_homogeneous (Point point, Point neighbor, scalar _s, void * data);
static double _boundary17 (Point point, Point neighbor, scalar _s, void * data);
static double _boundary17_homogeneous (Point point, Point neighbor, scalar _s, void * data);
#line 1 "60d_plate_advancing_simulation.c"
#line 18 "60d_plate_advancing_simulation.c"
#line 1 "navier-stokes/centered.h"
#line 1 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
#line 27 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
#line 1 "./run.h"
#line 1 "/home/fpl/softwares/basilisk/src/run.h"
#line 9 "/home/fpl/softwares/basilisk/src/run.h"
double dt = 1.;

#line 1 "./utils.h"
#line 1 "/home/fpl/softwares/basilisk/src/utils.h"







double DT = 1e10, CFL = 0.5;




struct {

  long nc;

  long tnc;

  double t;

  double speed;

  timer gt;
} perf;





void update_perf () {
  perf.nc += grid->n;
  perf.tnc += grid->tn;
  perf.t = timer_elapsed (perf.gt);
  perf.speed = perf.tnc/perf.t;
}






typedef struct {
  double cpu;
  double real;
  double speed;
  double min;
  double avg;
  double max;
  size_t tnc;
  long mem;
} timing;






timing timer_timing (timer t, int i, size_t tnc, double * mpi)
{
  timing s;
#if 1
  s.avg = mpi_time - t.tm;
#endif
  clock_t end = clock();
  s.cpu = ((double) (end - t.c))/CLOCKS_PER_SEC;
  s.real = timer_elapsed (t);
  if (tnc == 0) {
    double n = 0;
     { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _n = n;
{ double n = _n; NOT_UNUSED(n);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/utils.h", .line = 69,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 69 "/home/fpl/softwares/basilisk/src/utils.h"
 n++; } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 69

#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:n)) {

#line 69
foreach(){

#line 69 "/home/fpl/softwares/basilisk/src/utils.h"
 n++; } end_foreach();mpi_all_reduce_array (&n, double, MPI_SUM, 1);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 69
 }
    s.tnc = n;
    tnc = n*i;
  }
  else
    s.tnc = tnc;
#if _GNU_SOURCE
  struct rusage usage;
  getrusage (RUSAGE_SELF, &usage);
  s.mem = usage.ru_maxrss;
#else
  s.mem = 0;
#endif
#if 1
  if (mpi)
    MPI_Allgather (&s.avg, 1, MPI_DOUBLE, mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
  s.max = s.min = s.avg;
  mpi_all_reduce (s.max, MPI_DOUBLE, MPI_MAX);
  mpi_all_reduce (s.min, MPI_DOUBLE, MPI_MIN);
  mpi_all_reduce (s.avg, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.real, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.mem, MPI_LONG, MPI_SUM);
  s.real /= npe();
  s.avg /= npe();
  s.mem /= npe();
#else
  s.min = s.max = s.avg = 0.;
#endif
  s.speed = s.real > 0. ? tnc/s.real : -1.;
  return s;
}




void timer_print (timer t, int i, size_t tnc)
{
  timing s = timer_timing (t, i, tnc, NULL);
  fprintf (fout,
    "\n# " "Quadtree"
    ", %d steps, %g CPU, %.4g real, %.3g points.step/s, %d var\n",
    i, s.cpu, s.real, s.speed, (int) (datasize/sizeof(double)));
#if 1
  fprintf (fout,
    "# %d procs, MPI: min %.2g (%.2g%%) "
    "avg %.2g (%.2g%%) max %.2g (%.2g%%)\n",
    npe(),
    s.min, 100.*s.min/s.real,
    s.avg, 100.*s.avg/s.real,
    s.max, 100.*s.max/s.real);
#endif
}







typedef struct {
  double avg, rms, max, volume;
} norm;

norm normf (scalar f)
{
  double avg = 0., rms = 0., max = 0., volume = 0.;
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _max = max;
 double _avg = avg;
 double _rms = rms;
 double _volume = volume;
{ double max = _max; NOT_UNUSED(max);
 double avg = _avg; NOT_UNUSED(avg);
 double rms = _rms; NOT_UNUSED(rms);
 double volume = _volume; NOT_UNUSED(volume);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/utils.h", .line = 135,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 135
foreach_stencil(){

#line 136 "/home/fpl/softwares/basilisk/src/utils.h"

    IF (_stencil_val(__FILE__,__LINE__,f,0,0,0) != nodata && (sq(Delta)*val_cm(cm,0,0,0)) > 0.) {
      double v = fabs(_stencil_val(__FILE__,__LINE__,f,0,0,0));
      IF (v > max) max = v;
      volume += (sq(Delta)*val_cm(cm,0,0,0));
      avg += (sq(Delta)*val_cm(cm,0,0,0))*v;
      rms += (sq(Delta)*val_cm(cm,0,0,0))*sq(v);
    } } end_foreach_stencil(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 135
foreach_stencil(){

#line 136 "/home/fpl/softwares/basilisk/src/utils.h"

    IF (_stencil_val(__FILE__,__LINE__,f,0,0,0) != nodata && (sq(Delta)*val_cm(cm,0,0,0)) > 0.) {
      double v = fabs(_stencil_val(__FILE__,__LINE__,f,0,0,0));
      IF (v > max) max = v;
      volume += (sq(Delta)*val_cm(cm,0,0,0));
      avg += (sq(Delta)*val_cm(cm,0,0,0))*v;
      rms += (sq(Delta)*val_cm(cm,0,0,0))*sq(v);
    } } end_foreach_stencil(); }  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 143

#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:max)  reduction(+:avg) 
   reduction(+:rms)  reduction(+:volume)) {

#line 135

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 135
foreach(){

#line 136 "/home/fpl/softwares/basilisk/src/utils.h"

    if (val(f,0,0,0) != nodata && (sq(Delta)*val_cm(cm,0,0,0)) > 0.) {
      double v = fabs(val(f,0,0,0));
      if (v > max) max = v;
      volume += (sq(Delta)*val_cm(cm,0,0,0));
      avg += (sq(Delta)*val_cm(cm,0,0,0))*v;
      rms += (sq(Delta)*val_cm(cm,0,0,0))*sq(v);
    } } end_foreach(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 135
foreach(){

#line 136 "/home/fpl/softwares/basilisk/src/utils.h"

    if (val(f,0,0,0) != nodata && (sq(Delta)*val_cm(cm,0,0,0)) > 0.) {
      double v = fabs(val(f,0,0,0));
      if (v > max) max = v;
      volume += (sq(Delta)*val_cm(cm,0,0,0));
      avg += (sq(Delta)*val_cm(cm,0,0,0))*v;
      rms += (sq(Delta)*val_cm(cm,0,0,0))*sq(v);
    } } end_foreach(); }mpi_all_reduce_array (&max, double, MPI_MAX, 1);
mpi_all_reduce_array (&avg, double, MPI_SUM, 1);
mpi_all_reduce_array (&rms, double, MPI_SUM, 1);
mpi_all_reduce_array (&volume, double, MPI_SUM, 1);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 143
 }
  norm n;
  n.avg = volume ? avg/volume : 0.;
  n.rms = volume ? sqrt(rms/volume) : 0.;
  n.max = max;
  n.volume = volume;
  return n;
}





typedef struct {
  double min, max, sum, stddev, volume;
} stats;

stats statsf (scalar f)
{
  double min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _sum = sum;
 double _sum2 = sum2;
 double _volume = volume;
 double _max = max;
 double _min = min;
{ double sum = _sum; NOT_UNUSED(sum);
 double sum2 = _sum2; NOT_UNUSED(sum2);
 double volume = _volume; NOT_UNUSED(volume);
 double max = _max; NOT_UNUSED(max);
 double min = _min; NOT_UNUSED(min);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/utils.h", .line = 163,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 163
foreach_stencil(){

#line 164 "/home/fpl/softwares/basilisk/src/utils.h"

    IF ((sq(Delta)*val_cm(cm,0,0,0)) > 0. && _stencil_val(__FILE__,__LINE__,f,0,0,0) != nodata) {
      volume += (sq(Delta)*val_cm(cm,0,0,0));
      sum += (sq(Delta)*val_cm(cm,0,0,0))*_stencil_val(__FILE__,__LINE__,f,0,0,0);
      sum2 += (sq(Delta)*val_cm(cm,0,0,0))*sq(_stencil_val(__FILE__,__LINE__,f,0,0,0));
      IF (_stencil_val(__FILE__,__LINE__,f,0,0,0) > max) max = _stencil_val(__FILE__,__LINE__,f,0,0,0);
      IF (_stencil_val(__FILE__,__LINE__,f,0,0,0) < min) min = _stencil_val(__FILE__,__LINE__,f,0,0,0);
    } } end_foreach_stencil(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 163
foreach_stencil(){

#line 164 "/home/fpl/softwares/basilisk/src/utils.h"

    IF ((sq(Delta)*val_cm(cm,0,0,0)) > 0. && _stencil_val(__FILE__,__LINE__,f,0,0,0) != nodata) {
      volume += (sq(Delta)*val_cm(cm,0,0,0));
      sum += (sq(Delta)*val_cm(cm,0,0,0))*_stencil_val(__FILE__,__LINE__,f,0,0,0);
      sum2 += (sq(Delta)*val_cm(cm,0,0,0))*sq(_stencil_val(__FILE__,__LINE__,f,0,0,0));
      IF (_stencil_val(__FILE__,__LINE__,f,0,0,0) > max) max = _stencil_val(__FILE__,__LINE__,f,0,0,0);
      IF (_stencil_val(__FILE__,__LINE__,f,0,0,0) < min) min = _stencil_val(__FILE__,__LINE__,f,0,0,0);
    } } end_foreach_stencil(); }  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 171

#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:sum)  reduction(+:sum2)  reduction(+:volume) 
   reduction(max:max)  reduction(min:min)) {

#line 163

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 163
foreach(){

#line 164 "/home/fpl/softwares/basilisk/src/utils.h"

    if ((sq(Delta)*val_cm(cm,0,0,0)) > 0. && val(f,0,0,0) != nodata) {
      volume += (sq(Delta)*val_cm(cm,0,0,0));
      sum += (sq(Delta)*val_cm(cm,0,0,0))*val(f,0,0,0);
      sum2 += (sq(Delta)*val_cm(cm,0,0,0))*sq(val(f,0,0,0));
      if (val(f,0,0,0) > max) max = val(f,0,0,0);
      if (val(f,0,0,0) < min) min = val(f,0,0,0);
    } } end_foreach(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 163
foreach(){

#line 164 "/home/fpl/softwares/basilisk/src/utils.h"

    if ((sq(Delta)*val_cm(cm,0,0,0)) > 0. && val(f,0,0,0) != nodata) {
      volume += (sq(Delta)*val_cm(cm,0,0,0));
      sum += (sq(Delta)*val_cm(cm,0,0,0))*val(f,0,0,0);
      sum2 += (sq(Delta)*val_cm(cm,0,0,0))*sq(val(f,0,0,0));
      if (val(f,0,0,0) > max) max = val(f,0,0,0);
      if (val(f,0,0,0) < min) min = val(f,0,0,0);
    } } end_foreach(); }mpi_all_reduce_array (&sum, double, MPI_SUM, 1);
mpi_all_reduce_array (&sum2, double, MPI_SUM, 1);
mpi_all_reduce_array (&volume, double, MPI_SUM, 1);
mpi_all_reduce_array (&max, double, MPI_MAX, 1);
mpi_all_reduce_array (&min, double, MPI_MIN, 1);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 171
 }
  stats s;
  s.min = min, s.max = max, s.sum = sum, s.volume = volume;
  if (volume > 0.)
    sum2 -= sum*sum/volume;
  s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;
  return s;
}
#line 187 "/home/fpl/softwares/basilisk/src/utils.h"
static double generic_limiter (double r, double beta)
{
  double v1 = min (r, beta), v2 = min (beta*r, 1.);
  v1 = max (0., v1);
  return max (v1, v2);
}

double minmod (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 1.)*(s1 - s0);
}

double superbee (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 2.)*(s1 - s0);
}

double sweby (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 1.5)*(s1 - s0);
}
#line 213 "/home/fpl/softwares/basilisk/src/utils.h"
double theta = 1.3;

double minmod2 (double s0, double s1, double s2)
{
  if (s0 < s1 && s1 < s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 < d1) d1 = d2;
    return min(d1, d3);
  }
  if (s0 > s1 && s1 > s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 > d1) d1 = d2;
    return max(d1, d3);
  }
  return 0.;
}
#line 237 "/home/fpl/softwares/basilisk/src/utils.h"
void gradients (scalar * f, vector * g)
{
  if (!(list_len(f) == vectors_len(g))) qassert ("/home/fpl/softwares/basilisk/src/utils.h", 239, "list_len(f) == vectors_len(g)");
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/utils.h", .line = 240,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 240 "/home/fpl/softwares/basilisk/src/utils.h"
 {
    scalar s; vector v;
    scalar * _i0 = f; vector * _i1 = g; strongif (f) for (s = *f, v = *g; ((scalar *)&s)->i >= 0; s = *++_i0, v = *++_i1)
      {
#line 243
 {





   _stencil_val(__FILE__,__LINE__,v.x,0,0,0) = _attribute[s.i].gradient (_stencil_val(__FILE__,__LINE__,s,-1,0,0), _stencil_val(__FILE__,__LINE__,s,0,0,0), _stencil_val(__FILE__,__LINE__,s,1,0,0))/Delta;
      }
#line 243
 {





   _stencil_val(__FILE__,__LINE__,v.y,0,0,0) = _attribute[s.i].gradient (_stencil_val(__FILE__,__LINE__,s,0,-1,0), _stencil_val(__FILE__,__LINE__,s,0,0,0), _stencil_val(__FILE__,__LINE__,s,0,1,0))/Delta;
      }}
  } } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 251
foreach(){

#line 240 "/home/fpl/softwares/basilisk/src/utils.h"
 {
    scalar s; vector v;
    scalar * _i0 = f; vector * _i1 = g; strongif (f) for (s = *f, v = *g; ((scalar *)&s)->i >= 0; s = *++_i0, v = *++_i1)
      {
#line 243
 {





   val(v.x,0,0,0) = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0))/Delta;
      }
#line 243
 {





   val(v.y,0,0,0) = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0))/Delta;
      }}
  } } end_foreach(); }
}
#line 269 "/home/fpl/softwares/basilisk/src/utils.h"
void vorticity (const vector u, scalar omega)
{
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/utils.h", .line = 271,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 271
foreach_stencil(){

#line 271 "/home/fpl/softwares/basilisk/src/utils.h"

    _stencil_val(__FILE__,__LINE__,omega,0,0,0) = ((val_fm_x(fm.x,1,0,0) - val_fm_x(fm.x,0,0,0))*_stencil_val(__FILE__,__LINE__,u.y,0,0,0) +
        val_fm_x(fm.x,1,0,0)*_stencil_val(__FILE__,__LINE__,u.y,1,0,0) - val_fm_x(fm.x,0,0,0)*_stencil_val(__FILE__,__LINE__,u.y,-1,0,0) -
        (val_fm_y(fm.y,0,1,0) - val_fm_y(fm.y,0,0,0))*_stencil_val(__FILE__,__LINE__,u.x,0,0,0) +
        val_fm_y(fm.y,0,0,0)*_stencil_val(__FILE__,__LINE__,u.x,0,-1,0) - val_fm_y(fm.y,0,1,0)*_stencil_val(__FILE__,__LINE__,u.x,0,1,0))/(2.*val_cm(cm,0,0,0)*Delta + 0.); } end_foreach_stencil(); }
strongif (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 271
foreach_stencil(){

#line 271 "/home/fpl/softwares/basilisk/src/utils.h"

    _stencil_val(__FILE__,__LINE__,omega,0,0,0) = ((val_fm_x(fm.x,1,0,0) - val_fm_x(fm.x,0,0,0))*_stencil_val(__FILE__,__LINE__,u.y,0,0,0) +
        val_fm_x(fm.x,1,0,0)*_stencil_val(__FILE__,__LINE__,u.y,1,0,0) - val_fm_x(fm.x,0,0,0)*_stencil_val(__FILE__,__LINE__,u.y,-1,0,0) -
        (val_fm_y(fm.y,0,1,0) - val_fm_y(fm.y,0,0,0))*_stencil_val(__FILE__,__LINE__,u.x,0,0,0) +
        val_fm_y(fm.y,0,0,0)*_stencil_val(__FILE__,__LINE__,u.x,0,-1,0) - val_fm_y(fm.y,0,1,0)*_stencil_val(__FILE__,__LINE__,u.x,0,1,0))/(2.*val_cm(cm,0,0,0)*Delta + 0.); } end_foreach_stencil(); }
strongif (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 271
foreach_stencil(){

#line 271 "/home/fpl/softwares/basilisk/src/utils.h"

    _stencil_val(__FILE__,__LINE__,omega,0,0,0) = ((val_fm_x(fm.x,1,0,0) - val_fm_x(fm.x,0,0,0))*_stencil_val(__FILE__,__LINE__,u.y,0,0,0) +
        val_fm_x(fm.x,1,0,0)*_stencil_val(__FILE__,__LINE__,u.y,1,0,0) - val_fm_x(fm.x,0,0,0)*_stencil_val(__FILE__,__LINE__,u.y,-1,0,0) -
        (val_fm_y(fm.y,0,1,0) - val_fm_y(fm.y,0,0,0))*_stencil_val(__FILE__,__LINE__,u.x,0,0,0) +
        val_fm_y(fm.y,0,0,0)*_stencil_val(__FILE__,__LINE__,u.x,0,-1,0) - val_fm_y(fm.y,0,1,0)*_stencil_val(__FILE__,__LINE__,u.x,0,1,0))/(2.*val_cm(cm,0,0,0)*Delta + 0.); } end_foreach_stencil(); }
strongif (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 271
foreach_stencil(){

#line 271 "/home/fpl/softwares/basilisk/src/utils.h"

    _stencil_val(__FILE__,__LINE__,omega,0,0,0) = ((val_fm_x(fm.x,1,0,0) - val_fm_x(fm.x,0,0,0))*_stencil_val(__FILE__,__LINE__,u.y,0,0,0) +
        val_fm_x(fm.x,1,0,0)*_stencil_val(__FILE__,__LINE__,u.y,1,0,0) - val_fm_x(fm.x,0,0,0)*_stencil_val(__FILE__,__LINE__,u.y,-1,0,0) -
        (val_fm_y(fm.y,0,1,0) - val_fm_y(fm.y,0,0,0))*_stencil_val(__FILE__,__LINE__,u.x,0,0,0) +
        val_fm_y(fm.y,0,0,0)*_stencil_val(__FILE__,__LINE__,u.x,0,-1,0) - val_fm_y(fm.y,0,1,0)*_stencil_val(__FILE__,__LINE__,u.x,0,1,0))/(2.*val_cm(cm,0,0,0)*Delta + 0.); } end_foreach_stencil(); }  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 275

strongif (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 271
foreach(){

#line 271 "/home/fpl/softwares/basilisk/src/utils.h"

    val(omega,0,0,0) = ((val_fm_x(fm.x,1,0,0) - val_fm_x(fm.x,0,0,0))*val(u.y,0,0,0) +
        val_fm_x(fm.x,1,0,0)*val(u.y,1,0,0) - val_fm_x(fm.x,0,0,0)*val(u.y,-1,0,0) -
        (val_fm_y(fm.y,0,1,0) - val_fm_y(fm.y,0,0,0))*val(u.x,0,0,0) +
        val_fm_y(fm.y,0,0,0)*val(u.x,0,-1,0) - val_fm_y(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*val_cm(cm,0,0,0)*Delta + 0.); } end_foreach(); }
strongif (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 271
foreach(){

#line 271 "/home/fpl/softwares/basilisk/src/utils.h"

    val(omega,0,0,0) = ((val_fm_x(fm.x,1,0,0) - val_fm_x(fm.x,0,0,0))*val(u.y,0,0,0) +
        val_fm_x(fm.x,1,0,0)*val(u.y,1,0,0) - val_fm_x(fm.x,0,0,0)*val(u.y,-1,0,0) -
        (val_fm_y(fm.y,0,1,0) - val_fm_y(fm.y,0,0,0))*val(u.x,0,0,0) +
        val_fm_y(fm.y,0,0,0)*val(u.x,0,-1,0) - val_fm_y(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*val_cm(cm,0,0,0)*Delta + 0.); } end_foreach(); }
strongif (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 271
foreach(){

#line 271 "/home/fpl/softwares/basilisk/src/utils.h"

    val(omega,0,0,0) = ((val_fm_x(fm.x,1,0,0) - val_fm_x(fm.x,0,0,0))*val(u.y,0,0,0) +
        val_fm_x(fm.x,1,0,0)*val(u.y,1,0,0) - val_fm_x(fm.x,0,0,0)*val(u.y,-1,0,0) -
        (val_fm_y(fm.y,0,1,0) - val_fm_y(fm.y,0,0,0))*val(u.x,0,0,0) +
        val_fm_y(fm.y,0,0,0)*val(u.x,0,-1,0) - val_fm_y(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*val_cm(cm,0,0,0)*Delta + 0.); } end_foreach(); }
strongif (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 271
foreach(){

#line 271 "/home/fpl/softwares/basilisk/src/utils.h"

    val(omega,0,0,0) = ((val_fm_x(fm.x,1,0,0) - val_fm_x(fm.x,0,0,0))*val(u.y,0,0,0) +
        val_fm_x(fm.x,1,0,0)*val(u.y,1,0,0) - val_fm_x(fm.x,0,0,0)*val(u.y,-1,0,0) -
        (val_fm_y(fm.y,0,1,0) - val_fm_y(fm.y,0,0,0))*val(u.x,0,0,0) +
        val_fm_y(fm.y,0,0,0)*val(u.x,0,-1,0) - val_fm_y(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*val_cm(cm,0,0,0)*Delta + 0.); } end_foreach(); } }
}





double change (scalar s, scalar sn)
{
  double max = 0.;
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _max = max;
{ double max = _max; NOT_UNUSED(max);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/utils.h", .line = 285,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 285
foreach_stencil(){

#line 285 "/home/fpl/softwares/basilisk/src/utils.h"
 {
    IF ((sq(Delta)*val_cm(cm,0,0,0)) > 0.) {
      double ds = fabs (_stencil_val(__FILE__,__LINE__,s,0,0,0) - _stencil_val(__FILE__,__LINE__,sn,0,0,0));
      IF (ds > max)
 max = ds;
    }
    _stencil_val(__FILE__,__LINE__,sn,0,0,0) = _stencil_val(__FILE__,__LINE__,s,0,0,0);
  } } end_foreach_stencil(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 285
foreach_stencil(){

#line 285 "/home/fpl/softwares/basilisk/src/utils.h"
 {
    IF ((sq(Delta)*val_cm(cm,0,0,0)) > 0.) {
      double ds = fabs (_stencil_val(__FILE__,__LINE__,s,0,0,0) - _stencil_val(__FILE__,__LINE__,sn,0,0,0));
      IF (ds > max)
 max = ds;
    }
    _stencil_val(__FILE__,__LINE__,sn,0,0,0) = _stencil_val(__FILE__,__LINE__,s,0,0,0);
  } } end_foreach_stencil(); }  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 292

#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:max)) {

#line 285

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 285
foreach(){

#line 285 "/home/fpl/softwares/basilisk/src/utils.h"
 {
    if ((sq(Delta)*val_cm(cm,0,0,0)) > 0.) {
      double ds = fabs (val(s,0,0,0) - val(sn,0,0,0));
      if (ds > max)
 max = ds;
    }
    val(sn,0,0,0) = val(s,0,0,0);
  } } end_foreach(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 285
foreach(){

#line 285 "/home/fpl/softwares/basilisk/src/utils.h"
 {
    if ((sq(Delta)*val_cm(cm,0,0,0)) > 0.) {
      double ds = fabs (val(s,0,0,0) - val(sn,0,0,0));
      if (ds > max)
 max = ds;
    }
    val(sn,0,0,0) = val(s,0,0,0);
  } } end_foreach(); }mpi_all_reduce_array (&max, double, MPI_MAX, 1);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 292
 }
  return max;
}





scalar lookup_field (const char * name)
{
  if (name)
    strongif (all) for (scalar s = *all, *_i86 = all; ((scalar *)&s)->i >= 0; s = *++_i86)
      if (!strcmp (_attribute[s.i].name, name))
 return s;
  return (scalar){-1};
}

vector lookup_vector (const char * name)
{
  if (name) {
    char component[strlen(name) + 3];
    strcpy (component, name);
    strcat (component, ".x");
    strongif (all) for (scalar s = *all, *_i87 = all; ((scalar *)&s)->i >= 0; s = *++_i87)
      if (!strcmp (_attribute[s.i].name, component))
 return _attribute[s.i].v;
  }
  return (vector){{-1}};
}







#define foreach_segment(_S,_p) {\
  coord t = {(_S)[1].x - (_S)[0].x, (_S)[1].y - (_S)[0].y};\
  double norm = sqrt(sq(t.x) + sq(t.y));\
  if (!(norm > 0.)) qassert ("/home/fpl/softwares/basilisk/src/utils.h", 331, "norm > 0.");\
  t.x = t.x/norm + 1e-6, t.y = t.y/norm - 1.5e-6;\
  double alpha = ((_S)[0].x*((_S)[1].y - (_S)[0].y) -\
    (_S)[0].y*((_S)[1].x - (_S)[0].x))/norm;\
  foreach()\
    if (fabs(t.y*x - t.x*y - alpha) < 0.708*Delta) {\
      coord _o = {x,y}, _p[2];\
      int _n = 0;\
      {\
 if (t.x)\
   for (int _i = -1; _i <= 1 && _n < 2; _i += 2) {\
     _p[_n].x = _o.x + _i*Delta/2.;\
     double a = (_p[_n].x - (_S)[0].x)/t.x;\
     _p[_n].y = (_S)[0].y + a*t.y;\
     if (fabs(_p[_n].y - _o.y) <= Delta/2.) {\
       a = clamp (a, 0., norm);\
       _p[_n].x = (_S)[0].x + a*t.x, _p[_n].y = (_S)[0].y + a*t.y;\
       if (fabs(_p[_n].x - _o.x) <= Delta/2. &&\
    fabs(_p[_n].y - _o.y) <= Delta/2.)\
  _n++;\
     }\
   }\
 if (t.y)\
   for (int _i = -1; _i <= 1 && _n < 2; _i += 2) {\
     _p[_n].y = _o.y + _i*Delta/2.;\
     double a = (_p[_n].y - (_S)[0].y)/t.y;\
     _p[_n].x = (_S)[0].x + a*t.x;\
     if (fabs(_p[_n].x - _o.x) <= Delta/2.) {\
       a = clamp (a, 0., norm);\
       _p[_n].y = (_S)[0].y + a*t.y, _p[_n].x = (_S)[0].x + a*t.x;\
       if (fabs(_p[_n].y - _o.y) <= Delta/2. &&\
    fabs(_p[_n].x - _o.x) <= Delta/2.)\
  _n++;\
     }\
   }}\
      if (_n == 2) {\

#line 354

#define end_foreach_segment() } } end_foreach(); }




void fields_stats ()
{
  fprintf (ferr, "# t = %g, fields = {", t);
  strongif (all) for (scalar s = *all, *_i88 = all; ((scalar *)&s)->i >= 0; s = *++_i88)
    fprintf (ferr, " %s", _attribute[s.i].name);
  fputs (" }\n", ferr);
  fprintf (ferr, "# %12s: %12s %12s %12s %12s\n",
    "name", "min", "avg", "stddev", "max");
  strongif (all) for (scalar s = *all, *_i89 = all; ((scalar *)&s)->i >= 0; s = *++_i89) {
    stats ss = statsf (s);
    fprintf (ferr, "# %12s: %12g %12g %12g %12g\n",
      _attribute[s.i].name, ss.min, ss.sum/ss.volume, ss.stddev, ss.max);
  }
}

#line 1 "./output.h"
#line 1 "/home/fpl/softwares/basilisk/src/output.h"
#line 37 "/home/fpl/softwares/basilisk/src/output.h"
struct OutputField {
  scalar * list;
  FILE * fp;
  int n;
  bool linear;
  double box[2][2];
};


void output_field (struct OutputField p)
{ trace ("output_field", "/home/fpl/softwares/basilisk/src/output.h", 47);
  if (!p.list) p.list = all;
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = fout;
  p.n++;
  if (p.box[0][0] == 0. && p.box[0][1] == 0. &&
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0; p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
  }

  boundary_internal ((scalar *)(p.list), "/home/fpl/softwares/basilisk/src/output.h", 58);
  int len = list_len(p.list);
  double Delta = 0.999999*(p.box[1][0] - p.box[0][0])/(p.n - 1);
  int ny = (p.box[1][1] - p.box[0][1])/Delta + 1;
  double ** field = (double **) matrix_new (p.n, ny, len*sizeof(double));
  for (int i = 0; i < p.n; i++) {
    double x = Delta*i + p.box[0][0];
    for (int j = 0; j < ny; j++) {
      double y = Delta*j + p.box[0][1];
      if (p.linear) {
 int k = 0;
 strongif (p.list) for (scalar s = *p.list, *_i90 = p.list; ((scalar *)&s)->i >= 0; s = *++_i90)
   field[i][len*j + k++] = interpolate ((struct _interpolate){s, x, y});
      }
      else {
 Point point = locate ((struct _locate){x, y});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 73 "/home/fpl/softwares/basilisk/src/output.h"

 int k = 0;
 strongif (p.list) for (scalar s = *p.list, *_i91 = p.list; ((scalar *)&s)->i >= 0; s = *++_i91)
   field[i][len*j + k++] = point.level >= 0 ? val(s,0,0,0) : nodata;
      }
    }
  }

  if (pid() == 0) {
#if 1
    MPI_Reduce (MPI_IN_PLACE, field[0], len*p.n*ny, MPI_DOUBLE, MPI_MIN, 0,
  MPI_COMM_WORLD);
#endif
    fprintf (p.fp, "# 1:x 2:y");
    int i = 3;
    strongif (p.list) for (scalar s = *p.list, *_i92 = p.list; ((scalar *)&s)->i >= 0; s = *++_i92)
      fprintf (p.fp, " %d:%s", i++, _attribute[s.i].name);
    fputc('\n', p.fp);
    for (int i = 0; i < p.n; i++) {
      double x = Delta*i + p.box[0][0];
      for (int j = 0; j < ny; j++) {
 double y = Delta*j + p.box[0][1];

 fprintf (p.fp, "%g %g", x, y);
 int k = 0;
 strongif (p.list) for (scalar s = *p.list, *_i93 = p.list; ((scalar *)&s)->i >= 0; s = *++_i93)
   fprintf (p.fp, " %g", field[i][len*j + k++]);
 fputc ('\n', p.fp);
      }
      fputc ('\n', p.fp);
    }
    fflush (p.fp);
  }
#if 1
  else
    MPI_Reduce (field[0], NULL, len*p.n*ny, MPI_DOUBLE, MPI_MIN, 0,
  MPI_COMM_WORLD);
#endif

  matrix_free (field);
 end_trace("output_field", "/home/fpl/softwares/basilisk/src/output.h", 113); }
#line 141 "/home/fpl/softwares/basilisk/src/output.h"
struct OutputMatrix {
  scalar f;
  FILE * fp;
  int n;
  bool linear;
};


void output_matrix (struct OutputMatrix p)
{ trace ("output_matrix", "/home/fpl/softwares/basilisk/src/output.h", 150);
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = fout;
  if (p.linear) {
    scalar f = p.f;
    boundary_internal ((scalar *)(((scalar []){f,{-1}})), "/home/fpl/softwares/basilisk/src/output.h", 155);
  }
  float fn = p.n;
  float Delta = (float) L0/fn;
  fwrite (&fn, sizeof(float), 1, p.fp);
  for (int j = 0; j < p.n; j++) {
    float yp = (float) (Delta*j + X0 + Delta/2.);
    fwrite (&yp, sizeof(float), 1, p.fp);
  }
  for (int i = 0; i < p.n; i++) {
    float xp = (float) (Delta*i + X0 + Delta/2.);
    fwrite (&xp, sizeof(float), 1, p.fp);
    for (int j = 0; j < p.n; j++) {
      float yp = (float)(Delta*j + Y0 + Delta/2.), v;
      if (p.linear)
 v = interpolate ((struct _interpolate){p.f, xp, yp});
      else {
 Point point = locate ((struct _locate){xp, yp});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 172 "/home/fpl/softwares/basilisk/src/output.h"

 if (!(point.level >= 0)) qassert ("/home/fpl/softwares/basilisk/src/output.h", 173, "point.level >= 0");
 v = val(p.f,0,0,0);
      }
      fwrite (&v, sizeof(float), 1, p.fp);
    }
  }
  fflush (p.fp);
 end_trace("output_matrix", "/home/fpl/softwares/basilisk/src/output.h", 180); }
#line 189 "/home/fpl/softwares/basilisk/src/output.h"
typedef void (* colormap) (double cmap[127][3]);

void jet (double cmap[127][3])
{
  for (int i = 0; i < 127; i++) {
    cmap[i][0] =
      i <= 46 ? 0. :
      i >= 111 ? -0.03125*(i - 111) + 1. :
      i >= 78 ? 1. :
      0.03125*(i - 46);
    cmap[i][1] =
      i <= 14 || i >= 111 ? 0. :
      i >= 79 ? -0.03125*(i - 111) :
      i <= 46 ? 0.03125*(i - 14) :
      1.;
    cmap[i][2] =
      i >= 79 ? 0. :
      i >= 47 ? -0.03125*(i - 79) :
      i <= 14 ? 0.03125*(i - 14) + 1.:
      1.;
  }
}

void cool_warm (double cmap[127][3])
{






  static double basemap[33][3] = {
    {0.2298057, 0.298717966, 0.753683153},
    {0.26623388, 0.353094838, 0.801466763},
    {0.30386891, 0.406535296, 0.84495867},
    {0.342804478, 0.458757618, 0.883725899},
    {0.38301334, 0.50941904, 0.917387822},
    {0.424369608, 0.558148092, 0.945619588},
    {0.46666708, 0.604562568, 0.968154911},
    {0.509635204, 0.648280772, 0.98478814},
    {0.552953156, 0.688929332, 0.995375608},
    {0.596262162, 0.726149107, 0.999836203},
    {0.639176211, 0.759599947, 0.998151185},
    {0.681291281, 0.788964712, 0.990363227},
    {0.722193294, 0.813952739, 0.976574709},
    {0.761464949, 0.834302879, 0.956945269},
    {0.798691636, 0.849786142, 0.931688648},
    {0.833466556, 0.860207984, 0.901068838},
    {0.865395197, 0.86541021, 0.865395561},
    {0.897787179, 0.848937047, 0.820880546},
    {0.924127593, 0.827384882, 0.774508472},
    {0.944468518, 0.800927443, 0.726736146},
    {0.958852946, 0.769767752, 0.678007945},
    {0.96732803, 0.734132809, 0.628751763},
    {0.969954137, 0.694266682, 0.579375448},
    {0.966811177, 0.650421156, 0.530263762},
    {0.958003065, 0.602842431, 0.481775914},
    {0.943660866, 0.551750968, 0.434243684},
    {0.923944917, 0.49730856, 0.387970225},
    {0.89904617, 0.439559467, 0.343229596},
    {0.869186849, 0.378313092, 0.300267182},
    {0.834620542, 0.312874446, 0.259301199},
    {0.795631745, 0.24128379, 0.220525627},
    {0.752534934, 0.157246067, 0.184115123},
    {0.705673158, 0.01555616, 0.150232812}
  };

  for (int i = 0; i < 127; i++) {
    double x = i*(32 - 1e-10)/(127 - 1);
    int j = x; x -= j;
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (1. - x)*basemap[j][k] + x*basemap[j+1][k];
  }
}

void gray (double cmap[127][3])
{
  for (int i = 0; i < 127; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = i/(127 - 1.);
}

void randomap (double cmap[127][3])
{
  srand(0);
  for (int i = 0; i < 127; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (noise() + 1.)/2.;
}

void blue_white_red (double cmap[127][3])
{
  for (int i = 0; i < (127 + 1)/2; i++) {
    cmap[i][0] = i/((127 - 1)/2.);
    cmap[i][1] = i/((127 - 1)/2.);
    cmap[i][2] = 1.;
  }
  for (int i = 0; i < (127 - 1)/2; i++) {
    cmap[i + (127 + 1)/2][0] = 1.;
    cmap[i + (127 + 1)/2][1] = cmap[(127 - 3)/2 - i][1];
    cmap[i + (127 + 1)/2][2] = cmap[(127 - 3)/2 - i][1];
  }
}





typedef struct {
  unsigned char r, g, b;
} color;

color colormap_color (double cmap[127][3],
        double val, double min, double max)
{
  color c;
  if (val == nodata) {
    c.r = c.g = c.b = 0;
    return c;
  }
  int i;
  double coef;
  if (max != min)
    val = (val - min)/(max - min);
  else
    val = 0.;
  if (val <= 0.) i = 0, coef = 0.;
  else if (val >= 1.) i = 127 - 2, coef = 1.;
  else {
    i = val*(127 - 1);
    coef = val*(127 - 1) - i;
  }
  if (!(i < 127 - 1)) qassert ("/home/fpl/softwares/basilisk/src/output.h", 321, "i < NCMAP - 1");
  unsigned char * c1 = (unsigned char *) &c;
  for (int j = 0; j < 3; j++)
    c1[j] = 255*(cmap[i][j]*(1. - coef) + cmap[i + 1][j]*coef);
  return c;
}
#line 340 "/home/fpl/softwares/basilisk/src/output.h"
static const char * extension (const char * file, const char * ext) {
  int len = strlen(file);
  return len > 4 && !strcmp (file + len - 4, ext) ? file + len - 4 : NULL;
}

static const char * is_animation (const char * file) {
  const char * ext;
  if ((ext = extension (file, ".mp4")) ||
      (ext = extension (file, ".ogv")) ||
      (ext = extension (file, ".gif")))
    return ext;
  return NULL;
}

static struct {
  FILE ** fp;
  char ** names;
  int n;
} open_image_data = {NULL, NULL, 0};

static void open_image_cleanup ()
{
  for (int i = 0; i < open_image_data.n; i++) {
    qpclose (open_image_data.fp[i]);
    pfree (open_image_data.names[i],__func__,__FILE__,__LINE__);
  }
  pfree (open_image_data.fp,__func__,__FILE__,__LINE__);
  pfree (open_image_data.names,__func__,__FILE__,__LINE__);
  open_image_data.fp = NULL;
  open_image_data.names = NULL;
  open_image_data.n = 0;
}

static FILE * open_image_lookup (const char * file)
{
  for (int i = 0; i < open_image_data.n; i++)
    if (!strcmp (file, open_image_data.names[i]))
      return open_image_data.fp[i];
  return NULL;
}

static bool which (const char * command)
{
  char * s = getenv ("PATH");
  if (!s)
    return false;
  char path[strlen(s) + 1];
  strcpy (path, s);
  s = strtok (path, ":");
  while (s) {
    char f[strlen(s) + strlen(command) + 2];
    strcpy (f, s);
    strcat (f, "/");
    strcat (f, command);
    FILE * fp = fopen (f, "r");
    if (fp) {
      fclose (fp);
      return true;
    }
    s = strtok (NULL, ":");
  }
  return false;
}

static FILE * ppm_fallback (const char * file, const char * mode)
{
  char filename[strlen(file) + 5];
  strcpy (filename, file);
  strcat (filename, ".ppm");
  FILE * fp = fopen (filename, mode);
  if (!fp) {
    perror (file);

    MPI_Abort (MPI_COMM_WORLD, 1);

    exit (1);
  }
  return fp;
}

FILE * open_image (const char * file, const char * options)
{
  if (!(pid() == 0)) qassert ("/home/fpl/softwares/basilisk/src/output.h", 422, "pid() == 0");
  const char * ext;
  if ((ext = is_animation (file))) {
    FILE * fp = open_image_lookup (file);
    if (fp)
      return fp;

    int len = strlen ("ppm2???    ") + strlen (file) +
      (options ? strlen (options) : 0);
    char command[len];
    strcpy (command, "ppm2"); strcat (command, ext + 1);

    static int has_ffmpeg = -1;
    if (has_ffmpeg < 0) {
      if (which (command) && (which ("ffmpeg") || which ("avconv")))
 has_ffmpeg = true;
      else {
 fprintf (ferr,
   "open_image(): cannot find '%s' or 'ffmpeg'/'avconv'\n"
   "  falling back to raw PPM outputs\n", command);
 has_ffmpeg = false;
      }
    }
    if (!has_ffmpeg)
      return ppm_fallback (file, "a");

    static bool added = false;
    if (!added) {
      free_solver_func_add (open_image_cleanup);
      added = true;
    }
    open_image_data.n++;
    open_image_data.names = (char * *) prealloc (open_image_data.names, (open_image_data.n)*sizeof(char *),__func__,__FILE__,__LINE__);
    open_image_data.names[open_image_data.n - 1] = pstrdup (file,__func__,__FILE__,__LINE__);

    if (options) {
      strcat (command, " ");
      strcat (command, options);
    }
    strcat (command, !strcmp (ext, ".mp4") ? " " : " > ");
    strcat (command, file);
    open_image_data.fp = (FILE * *) prealloc (open_image_data.fp, (open_image_data.n)*sizeof(FILE *),__func__,__FILE__,__LINE__);
    return open_image_data.fp[open_image_data.n - 1] = qpopen (command, "w");
  }
  else {
    static int has_convert = -1;
    if (has_convert < 0) {
      if (which ("convert"))
 has_convert = true;
      else {
 fprintf (ferr,
   "open_image(): cannot find 'convert'\n"
   "  falling back to raw PPM outputs\n");
 has_convert = false;
      }
    }
    if (!has_convert)
      return ppm_fallback (file, "w");

    int len = strlen ("convert ppm:-   ") + strlen (file) +
      (options ? strlen (options) : 0);
    char command[len];
    strcpy (command, "convert ppm:- ");
    if (options) {
      strcat (command, options);
      strcat (command, " ");
    }
    strcat (command, file);
    return qpopen (command, "w");
  }
}

void close_image (const char * file, FILE * fp)
{
  if (!(pid() == 0)) qassert ("/home/fpl/softwares/basilisk/src/output.h", 496, "pid() == 0");
  if (is_animation (file)) {
    if (!open_image_lookup (file))
      fclose (fp);
  }
  else if (which ("convert"))
    qpclose (fp);
  else
    fclose (fp);
}
#line 571 "/home/fpl/softwares/basilisk/src/output.h"
struct OutputPPM {
  scalar f;
  FILE * fp;
  int n;
  char * file;
  double min, max, spread, z;
  bool linear;
  double box[2][2];
  scalar mask;
  colormap map;
  char * opt;
};


void output_ppm (struct OutputPPM p)
{ trace ("output_ppm", "/home/fpl/softwares/basilisk/src/output.h", 586);

  if (p.n == 0) p.n = N;
  if (p.min == 0 && p.max == 0) {
    stats s = statsf (p.f);
    if (p.spread < 0.)
      p.min = s.min, p.max = s.max;
    else {
      double avg = s.sum/s.volume, spread = (p.spread ? p.spread : 5.)*s.stddev;
      p.min = avg - spread; p.max = avg + spread;
    }
  }
  if (p.box[0][0] == 0. && p.box[0][1] == 0. &&
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0; p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
  }
  if (!p.map)
    p.map = jet;
  if (p.linear) {
    scalar f = p.f, mask = p.mask;
    if (mask.i)
      boundary_internal ((scalar *)(((scalar []){f,mask,{-1}})), "/home/fpl/softwares/basilisk/src/output.h", 608);
    else
      boundary_internal ((scalar *)(((scalar []){f,{-1}})), "/home/fpl/softwares/basilisk/src/output.h", 610);
  }

  double fn = p.n;
  double Delta = (p.box[1][0] - p.box[0][0])/fn;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;
  if (ny % 2) ny++;

  color ** ppm = (color **) matrix_new (ny, p.n, sizeof(color));
  double cmap[127][3];
  p.map (cmap);
  OMP_PARALLEL() {
    OMP(omp for schedule(static))
      for (int j = 0; j < ny; j++) {
 double yp = Delta*j + p.box[0][1] + Delta/2.;
 for (int i = 0; i < p.n; i++) {
   double xp = Delta*i + p.box[0][0] + Delta/2., v;
   if (p.mask.i) {
     if (p.linear) {
       double m = interpolate ((struct _interpolate){p.mask, xp, yp, p.z});
       if (m < 0.)
  v = nodata;
       else
  v = interpolate ((struct _interpolate){p.f, xp, yp, p.z});
     }
     else {
       Point point = locate ((struct _locate){xp, yp, p.z});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 636 "/home/fpl/softwares/basilisk/src/output.h"

       if (point.level < 0 || val(p.mask,0,0,0) < 0.)
  v = nodata;
       else
  v = val(p.f,0,0,0);
     }
   }
   else if (p.linear)
     v = interpolate ((struct _interpolate){p.f, xp, yp, p.z});
   else {
     Point point = locate ((struct _locate){xp, yp, p.z});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 646 "/home/fpl/softwares/basilisk/src/output.h"

     v = point.level >= 0 ? val(p.f,0,0,0) : nodata;
   }
   ppm[ny - 1 - j][i] = colormap_color (cmap, v, p.min, p.max);
 }
      }
  }

  if (pid() == 0) {
#if 1
    MPI_Reduce (MPI_IN_PLACE, ppm[0], 3*ny*p.n, MPI_UNSIGNED_CHAR, MPI_MAX, 0,
  MPI_COMM_WORLD);
#endif
    if (!p.fp) p.fp = fout;
    if (p.file)
      p.fp = open_image (p.file, p.opt);

    fprintf (p.fp, "P6\n%u %u 255\n", p.n, ny);
    fwrite (((void **) ppm)[0], sizeof(color), ny*p.n, p.fp);

    if (p.file)
      close_image (p.file, p.fp);
    else
      fflush (p.fp);
  }
#if 1
  else
    MPI_Reduce (ppm[0], NULL, 3*ny*p.n, MPI_UNSIGNED_CHAR, MPI_MAX, 0,
  MPI_COMM_WORLD);
#endif

  matrix_free (ppm);
 end_trace("output_ppm", "/home/fpl/softwares/basilisk/src/output.h", 678); }
#line 710 "/home/fpl/softwares/basilisk/src/output.h"
struct OutputGRD {
  scalar f;
  FILE * fp;
  double Delta;
  bool linear;
  double box[2][2];
  scalar mask;
};


void output_grd (struct OutputGRD p)
{ trace ("output_grd", "/home/fpl/softwares/basilisk/src/output.h", 721);

  if (!p.fp) p.fp = fout;
  if (p.box[0][0] == 0. && p.box[0][1] == 0. &&
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0; p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
    if (p.Delta == 0) p.Delta = L0/N;
  }
  if (p.linear) {
    scalar f = p.f, mask = p.mask;
    if (mask.i)
      boundary_internal ((scalar *)(((scalar []){f,mask,{-1}})), "/home/fpl/softwares/basilisk/src/output.h", 733);
    else
      boundary_internal ((scalar *)(((scalar []){f,{-1}})), "/home/fpl/softwares/basilisk/src/output.h", 735);
  }

  double Delta = p.Delta;
  int nx = (p.box[1][0] - p.box[0][0])/Delta;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;


  fprintf (p.fp, "ncols          %d\n", nx);
  fprintf (p.fp, "nrows          %d\n", ny);
  fprintf (p.fp, "xllcorner      %g\n", p.box[0][0]);
  fprintf (p.fp, "yllcorner      %g\n", p.box[0][1]);
  fprintf (p.fp, "cellsize       %g\n", Delta);
  fprintf (p.fp, "nodata_value   -9999\n");


  for (int j = ny-1; j >= 0; j--) {
    double yp = Delta*j + p.box[0][1] + Delta/2.;
    for (int i = 0; i < nx; i++) {
      double xp = Delta*i + p.box[0][0] + Delta/2., v;
      if (p.mask.i) {
 if (p.linear) {
   double m = interpolate ((struct _interpolate){p.mask, xp, yp});
   if (m < 0.)
     v = nodata;
   else
     v = interpolate ((struct _interpolate){p.f, xp, yp});
 }
 else {
   Point point = locate ((struct _locate){xp, yp});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 764 "/home/fpl/softwares/basilisk/src/output.h"

   if (point.level < 0 || val(p.mask,0,0,0) < 0.)
     v = nodata;
   else
     v = val(p.f,0,0,0);
 }
      }
      else if (p.linear)
 v = interpolate ((struct _interpolate){p.f, xp, yp});
      else {
 Point point = locate ((struct _locate){xp, yp});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 774 "/home/fpl/softwares/basilisk/src/output.h"

 v = point.level >= 0 ? val(p.f,0,0,0) : nodata;
      }
      if (v == nodata)
 fprintf (p.fp, "-9999 ");
      else
 fprintf (p.fp, "%f ", v);
    }
    fprintf (p.fp, "\n");
  }

  fflush (p.fp);
 end_trace("output_grd", "/home/fpl/softwares/basilisk/src/output.h", 786); }
#line 813 "/home/fpl/softwares/basilisk/src/output.h"
struct OutputGfs {
  FILE * fp;
  scalar * list;
  double t;
  char * file;
  bool translate;
};

static char * replace (const char * input, int target, int with,
         bool translate)
{
  if (translate) {
    if (!strcmp (input, "u.x"))
      return pstrdup ("U",__func__,__FILE__,__LINE__);
    if (!strcmp (input, "u.y"))
      return pstrdup ("V",__func__,__FILE__,__LINE__);
    if (!strcmp (input, "u.z"))
      return pstrdup ("W",__func__,__FILE__,__LINE__);
  }
  char * name = pstrdup (input,__func__,__FILE__,__LINE__), * i = name;
  while (*i != '\0') {
    if (*i == target)
      *i = with;
    i++;
  }
  return name;
}


void output_gfs (struct OutputGfs p)
{ trace ("output_gfs", "/home/fpl/softwares/basilisk/src/output.h", 843);
  char * fname = p.file;

#if 1



  FILE * fp = p.fp;
  if (p.file == NULL) {
    long pid = getpid();
    MPI_Bcast (&pid, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    fname = ((char *) pmalloc ((80)*sizeof(char),__func__,__FILE__,__LINE__));
    snprintf (fname, 80, ".output-%ld", pid);
    p.fp = NULL;
  }
#endif

  bool opened = false;
  if (p.fp == NULL) {
    if (fname == NULL)
      p.fp = fout;
    else if (!(p.fp = fopen (fname, "w"))) {
      perror (fname);
      exit (1);
    }
    else
      opened = true;
  }

  scalar * list = p.list ? p.list : list_copy (all);

  restriction (list);
  fprintf (p.fp,
    "1 0 GfsSimulation GfsBox GfsGEdge { binary = 1"
    " x = %g y = %g ",
    0.5 + X0/L0, 0.5 + Y0/L0);




  if (list != NULL && list[0].i != -1) {
    scalar s = list[0];
    char * name = replace (_attribute[s.i].name, '.', '_', p.translate);
    fprintf (p.fp, "variables = %s", name);
    pfree (name,__func__,__FILE__,__LINE__);
    for (int i = 1; i < list_len(list); i++) {
      scalar s = list[i];
      if (_attribute[s.i].name) {
 char * name = replace (_attribute[s.i].name, '.', '_', p.translate);
 fprintf (p.fp, ",%s", name);
 pfree (name,__func__,__FILE__,__LINE__);
      }
    }
    fprintf (p.fp, " ");
  }
  fprintf (p.fp, "} {\n");
  fprintf (p.fp, "  Time { t = %g }\n", t);
  if (L0 != 1.)
    fprintf (p.fp, "  PhysicalParams { L = %g }\n", L0);
  fprintf (p.fp, "  VariableTracerVOF f\n");
  fprintf (p.fp, "}\nGfsBox { x = 0 y = 0 z = 0 } {\n");

#if 1
  long header;
  if ((header = ftell (p.fp)) < 0) {
    perror ("output_gfs(): error in header");
    exit (1);
  }
  int cell_size = sizeof(unsigned) + sizeof(double);
  strongif (list) for (scalar s = *list, *_i94 = list; ((scalar *)&s)->i >= 0; s = *++_i94)
    if (_attribute[s.i].name)
      cell_size += sizeof(double);
  scalar index = new_scalar("index");
  size_t total_size = header + (z_indexing (index, false) + 1)*cell_size;
#endif



   { foreach_cell(){

#line 921 "/home/fpl/softwares/basilisk/src/output.h"
 {
#if 1
    if (is_local(cell))
#endif
    {
#if 1
      if (fseek (p.fp, header + val(index,0,0,0)*cell_size, SEEK_SET) < 0) {
 perror ("output_gfs(): error while seeking");
 exit (1);
      }
#endif
      unsigned flags =
 level == 0 ? 0 :



      child.x == -1 && child.y == -1 ? 0 :
 child.x == -1 && child.y == 1 ? 1 :
 child.x == 1 && child.y == -1 ? 2 :
 3;
#line 951 "/home/fpl/softwares/basilisk/src/output.h"
      if (is_leaf(cell))
 flags |= (1 << 4);
      fwrite (&flags, sizeof (unsigned), 1, p.fp);
      double a = -1;
      fwrite (&a, sizeof (double), 1, p.fp);
      strongif (list) for (scalar s = *list, *_i95 = list; ((scalar *)&s)->i >= 0; s = *++_i95)
 if (_attribute[s.i].name) {
   if (_attribute[s.i].v.x.i >= 0) {




     if (_attribute[s.i].v.x.i == s.i) {
       s = _attribute[s.i].v.y;
       a = is_local(cell) && val(s,0,0,0) != nodata ? val(s,0,0,0) : (double) DBL_MAX;
     }
     else if (_attribute[s.i].v.y.i == s.i) {
       s = _attribute[s.i].v.x;
       a = is_local(cell) && val(s,0,0,0) != nodata ? - val(s,0,0,0) : (double) DBL_MAX;
     }





   }
   else
     a = is_local(cell) && val(s,0,0,0) != nodata ? val(s,0,0,0) : (double) DBL_MAX;
   fwrite (&a, sizeof (double), 1, p.fp);
 }
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

#if 1
  delete (((scalar []){index,{-1}}));
  if (!pid() && fseek (p.fp, total_size, SEEK_SET) < 0) {
    perror ("output_gfs(): error while finishing");
    exit (1);
  }
  if (!pid())
#endif
    fputs ("}\n", p.fp);
  fflush (p.fp);

  if (!p.list)
    pfree (list,__func__,__FILE__,__LINE__);
  if (opened)
    fclose (p.fp);

#if 1
  if (p.file == NULL) {
    MPI_Barrier (MPI_COMM_WORLD);
    if (pid() == 0) {
      if (fp == NULL)
 fp = fout;
      p.fp = fopen (fname, "r");
      size_t l;
      unsigned char buffer[8192];
      while ((l = fread (buffer, 1, 8192, p.fp)) > 0)
 fwrite (buffer, 1, l, fp);
      fflush (fp);
      remove (fname);
    }
    pfree (fname,__func__,__FILE__,__LINE__);
  }
#endif
 end_trace("output_gfs", "/home/fpl/softwares/basilisk/src/output.h", 1019); }
#line 1043 "/home/fpl/softwares/basilisk/src/output.h"
struct Dump {
  char * file;
  scalar * list;
  FILE * fp;
  bool unbuffered;
};

struct DumpHeader {
  double t;
  long len;
  int i, depth, npe, version;
  coord n;
};

static const int dump_version =

  170901;

static scalar * dump_list (scalar * lista)
{
  scalar * list = is_constant(cm) ? NULL : list_concat (((scalar []){cm,{-1}}), NULL);
  strongif (lista) for (scalar s = *lista, *_i96 = lista; ((scalar *)&s)->i >= 0; s = *++_i96)
    if (!_attribute[s.i].face && !_attribute[s.i].nodump && s.i != cm.i)
      list = list_add (list, s);
  return list;
}

static void dump_header (FILE * fp, struct DumpHeader * header, scalar * list)
{
  if (fwrite (header, sizeof(struct DumpHeader), 1, fp) < 1) {
    perror ("dump(): error while writing header");
    exit (1);
  }
  strongif (list) for (scalar s = *list, *_i97 = list; ((scalar *)&s)->i >= 0; s = *++_i97) {
    unsigned len = strlen(_attribute[s.i].name);
    if (fwrite (&len, sizeof(unsigned), 1, fp) < 1) {
      perror ("dump(): error while writing len");
      exit (1);
    }
    if (fwrite (_attribute[s.i].name, sizeof(char), len, fp) < len) {
      perror ("dump(): error while writing s.name");
      exit (1);
    }
  }
  double o[4] = {X0,Y0,Z0,L0};
  if (fwrite (o, sizeof(double), 4, fp) < 4) {
    perror ("dump(): error while writing coordinates");
    exit (1);
  }
}

#if !1

void dump (struct Dump p)
{ trace ("dump", "/home/fpl/softwares/basilisk/src/output.h", 1097);
  FILE * fp = p.fp;
  char def[] = "dump", * file = p.file ? p.file : p.fp ? NULL : def;

  char * name = NULL;
  if (file) {
    name = (char *) pmalloc (strlen(file) + 2,__func__,__FILE__,__LINE__);
    strcpy (name, file);
    if (!p.unbuffered)
      strcat (name, "~");
    if ((fp = fopen (name, "w")) == NULL) {
      perror (name);
      exit (1);
    }
  }
  if (!(fp)) qassert ("/home/fpl/softwares/basilisk/src/output.h", 1112, "fp");

  scalar * dlist = dump_list (p.list ? p.list : all);
  scalar size= new_scalar("size");
  scalar * list = list_concat (((scalar []){size,{-1}}), dlist); pfree (dlist,__func__,__FILE__,__LINE__);
  struct DumpHeader header = { t, list_len(list), iter, depth(), npe(),
          dump_version };
  dump_header (fp, &header, list);

  subtree_size (size, false);

   { foreach_cell(){

#line 1123 "/home/fpl/softwares/basilisk/src/output.h"
 {
    unsigned flags = is_leaf(cell) ? leaf : 0;
    if (fwrite (&flags, sizeof(unsigned), 1, fp) < 1) {
      perror ("dump(): error while writing flags");
      exit (1);
    }
    strongif (list) for (scalar s = *list, *_i98 = list; ((scalar *)&s)->i >= 0; s = *++_i98)
      if (fwrite (&val(s,0,0,0), sizeof(double), 1, fp) < 1) {
 perror ("dump(): error while writing scalars");
 exit (1);
      }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

  pfree (list,__func__,__FILE__,__LINE__);
  if (file) {
    fclose (fp);
    if (!p.unbuffered)
      rename (name, file);
    pfree (name,__func__,__FILE__,__LINE__);
  }
 delete (((scalar []){size,{-1}}));  end_trace("dump", "/home/fpl/softwares/basilisk/src/output.h", 1145); }
#else

void dump (struct Dump p)
{ trace ("dump", "/home/fpl/softwares/basilisk/src/output.h", 1149);
  FILE * fp = p.fp;
  char def[] = "dump", * file = p.file ? p.file : p.fp ? NULL : def;

  if (fp != NULL || file == NULL) {
    fprintf (ferr, "dump(): must specify a file name when using MPI\n");
    exit(1);
  }

  char name[strlen(file) + 2];
  strcpy (name, file);
  if (!p.unbuffered)
    strcat (name, "~");
  FILE * fh = fopen (name, "w");
  if (fh == NULL) {
    perror (name);
    exit (1);
  }

  scalar * dlist = dump_list (p.list ? p.list : all);
  scalar size= new_scalar("size");
  scalar * list = list_concat (((scalar []){size,{-1}}), dlist); pfree (dlist,__func__,__FILE__,__LINE__);
  struct DumpHeader header = { t, list_len(list), iter, depth(), npe(),
          dump_version };







  if (pid() == 0)
    dump_header (fh, &header, list);

  scalar index = {-1};

  index = new_scalar("index");
  z_indexing (index, false);
  int cell_size = sizeof(unsigned) + header.len*sizeof(double);
  int sizeofheader = sizeof(header) + 4*sizeof(double);
  strongif (list) for (scalar s = *list, *_i99 = list; ((scalar *)&s)->i >= 0; s = *++_i99)
    sizeofheader += sizeof(unsigned) + sizeof(char)*strlen(_attribute[s.i].name);
  long pos = pid() ? 0 : sizeofheader;

  subtree_size (size, false);

   { foreach_cell(){

#line 1195 "/home/fpl/softwares/basilisk/src/output.h"
 {

    if (is_local(cell)) {
      long offset = sizeofheader + val(index,0,0,0)*cell_size;
      if (pos != offset) {
 fseek (fh, offset, SEEK_SET);
 pos = offset;
      }
      unsigned flags = is_leaf(cell) ? leaf : 0;
      fwrite (&flags, 1, sizeof(unsigned), fh);
      strongif (list) for (scalar s = *list, *_i100 = list; ((scalar *)&s)->i >= 0; s = *++_i100)
 fwrite (&val(s,0,0,0), 1, sizeof(double), fh);
      pos += cell_size;
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

  delete (((scalar []){index,{-1}}));

  pfree (list,__func__,__FILE__,__LINE__);
  fclose (fh);
  if (!p.unbuffered && pid() == 0)
    rename (name, file);
 delete (((scalar []){size,{-1}}));  end_trace("dump", "/home/fpl/softwares/basilisk/src/output.h", 1219); }
#endif


bool restore (struct Dump p)
{ trace ("restore", "/home/fpl/softwares/basilisk/src/output.h", 1224);
  FILE * fp = p.fp;
  char * file = p.file;
  if (file && (fp = fopen (file, "r")) == NULL)
    { bool _ret =  false; end_trace("restore", "/home/fpl/softwares/basilisk/src/output.h", 1228);  return _ret; }
  if (!(fp)) qassert ("/home/fpl/softwares/basilisk/src/output.h", 1229, "fp");

  struct DumpHeader header;
  if (fread (&header, sizeof(header), 1, fp) < 1) {
    fprintf (ferr, "restore(): error: expecting header\n");
    exit (1);
  }


  init_grid (1);
   { foreach_cell(){

#line 1239 "/home/fpl/softwares/basilisk/src/output.h"
 {
    cell.pid = pid();
    cell.flags |= active;
  } } end_foreach_cell(); }
  ((Tree *)grid)->dirty = true;
#line 1264 "/home/fpl/softwares/basilisk/src/output.h"
  bool restore_all = (p.list == all);
  scalar * list = dump_list (p.list ? p.list : all);
  if (header.version == 161020) {
    if (header.len - 1 != list_len (list)) {
      fprintf (ferr,
        "restore(): error: the list lengths don't match: "
        "%ld (file) != %d (code)\n",
        header.len - 1, list_len (list));
      exit (1);
    }
  }
  else {
    if (header.version != dump_version) {
      fprintf (ferr,
        "restore(): error: file version mismatch: "
        "%d (file) != %d (code)\n",
        header.version, dump_version);
      exit (1);
    }

    scalar * input = NULL;
    for (int i = 0; i < header.len; i++) {
      unsigned len;
      if (fread (&len, sizeof(unsigned), 1, fp) < 1) {
 fprintf (ferr, "restore(): error: expecting len\n");
 exit (1);
      }
      char name[len + 1];
      if (fread (name, sizeof(char), len, fp) < 1) {
 fprintf (ferr, "restore(): error: expecting s.name\n");
 exit (1);
      }
      name[len] = '\0';

      if (i > 0) {
 bool found = false;
 strongif (list) for (scalar s = *list, *_i101 = list; ((scalar *)&s)->i >= 0; s = *++_i101)
   if (!strcmp (_attribute[s.i].name, name)) {
     input = list_append (input, s);
     found = true; break;
   }
 if (!found) {
   if (restore_all) {
     scalar s = new_scalar("s");
     pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
     _attribute[s.i].name = pstrdup (name,__func__,__FILE__,__LINE__);
     input = list_append (input, s);
   }
   else
     input = list_append (input, (scalar){INT_MAX});
 }
      }
    }
    pfree (list,__func__,__FILE__,__LINE__);
    list = input;

    double o[4];
    if (fread (o, sizeof(double), 4, fp) < 4) {
      fprintf (ferr, "restore(): error: expecting coordinates\n");
      exit (1);
    }
    origin ((struct _origin){o[0], o[1], o[2]});
    size (o[3]);
  }
#line 1339 "/home/fpl/softwares/basilisk/src/output.h"
  scalar * listm = is_constant(cm) ? NULL : (scalar *)((vector []){{fm.x,fm.y},{{-1},{-1}}});

  restore_mpi (fp, list);
#line 1369 "/home/fpl/softwares/basilisk/src/output.h"
  scalar * other = NULL;
  strongif (all) for (scalar s = *all, *_i102 = all; ((scalar *)&s)->i >= 0; s = *++_i102)
    if (!list_lookup (list, s) && !list_lookup (listm, s))
      other = list_append (other, s);
  reset (other, 0.);
  pfree (other,__func__,__FILE__,__LINE__);

  pfree (list,__func__,__FILE__,__LINE__);
  if (file)
    fclose (fp);


  while (iter < header.i && events (false))
    iter = inext;
  events (false);
  while (t < header.t && events (false))
    t = tnext;
  t = header.t;
  events (false);

  { bool _ret =  true; end_trace("restore", "/home/fpl/softwares/basilisk/src/output.h", 1389);  return _ret; }
 end_trace("restore", "/home/fpl/softwares/basilisk/src/output.h", 1390); }
#line 376 "/home/fpl/softwares/basilisk/src/utils.h"
#line 12 "/home/fpl/softwares/basilisk/src/run.h"


void run (void)
{ trace ("run", "/home/fpl/softwares/basilisk/src/run.h", 15);
  iter = 0, t = 0., dt = 1.;
  init_grid (N);

  perf.nc = perf.tnc = 0;
  perf.gt = timer_start();
  while (events (true)) {





    update_perf();
    iter = inext, t = tnext;
  }




  timer_print (perf.gt, iter, perf.tnc);

  free_grid();
 end_trace("run", "/home/fpl/softwares/basilisk/src/run.h", 37); }




static int defaults_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults (const int i, const double t, Event * _ev) { trace ("defaults", "/home/fpl/softwares/basilisk/src/run.h", 42);  {
  display ((struct _display){"box();"});
 end_trace("defaults", "/home/fpl/softwares/basilisk/src/run.h", 44); } return 0; } 





static int cleanup_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (t = 1234567890);   *ip = i; *tp = t;   return ret; } static int cleanup (const int i, const double t, Event * _ev) { trace ("cleanup", "/home/fpl/softwares/basilisk/src/run.h", 50);  {
  display ((struct _display){"", true});
 end_trace("cleanup", "/home/fpl/softwares/basilisk/src/run.h", 52); } return 0; } 
#line 28 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
#line 1 "./timestep.h"
#line 1 "/home/fpl/softwares/basilisk/src/timestep.h"

double timestep (const vector u, double dtmax)
{
  static double previous = 0.;
  dtmax /= CFL;
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dtmax = dtmax;
{ double dtmax = _dtmax; NOT_UNUSED(dtmax);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/timestep.h", .line = 6,
    .each = "foreach_face", .first = _first_call
  };

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 6
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 6
{

#line 6 "/home/fpl/softwares/basilisk/src/timestep.h"

    IF (_stencil_val(__FILE__,__LINE__,u.x,0,0,0) != 0.) {
      double dt = Delta/fabs(_stencil_val(__FILE__,__LINE__,u.x,0,0,0));




      dt *= val_cm(cm,0,0,0);

      IF (dt < dtmax) dtmax = dt;
    } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 6
{

#line 6 "/home/fpl/softwares/basilisk/src/timestep.h"

    IF (_stencil_val(__FILE__,__LINE__,u.y,0,0,0) != 0.) {
      double dt = Delta/fabs(_stencil_val(__FILE__,__LINE__,u.y,0,0,0));




      dt *= val_cm(cm,0,0,0);

      IF (dt < dtmax) dtmax = dt;
    } }  }}  end_foreach_face_stencil()
#line 16
 }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 6
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 6
{

#line 6 "/home/fpl/softwares/basilisk/src/timestep.h"

    IF (_stencil_val(__FILE__,__LINE__,u.x,0,0,0) != 0.) {
      double dt = Delta/fabs(_stencil_val(__FILE__,__LINE__,u.x,0,0,0));




      dt *= val_cm(cm,0,0,0);

      IF (dt < dtmax) dtmax = dt;
    } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 6
{

#line 6 "/home/fpl/softwares/basilisk/src/timestep.h"

    IF (_stencil_val(__FILE__,__LINE__,u.y,0,0,0) != 0.) {
      double dt = Delta/fabs(_stencil_val(__FILE__,__LINE__,u.y,0,0,0));




      dt *= val_cm(cm,0,0,0);

      IF (dt < dtmax) dtmax = dt;
    } }  }}  end_foreach_face_stencil()
#line 16
 }  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 16

#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(min:dtmax)) {

#line 6

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 6
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 6
{

#line 6 "/home/fpl/softwares/basilisk/src/timestep.h"

    if (val(u.x,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.x,0,0,0));




      dt *= val_cm(cm,0,0,0);

      if (dt < dtmax) dtmax = dt;
    } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 6
{

#line 6 "/home/fpl/softwares/basilisk/src/timestep.h"

    if (val(u.y,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.y,0,0,0));




      dt *= val_cm(cm,0,0,0);

      if (dt < dtmax) dtmax = dt;
    } }  }}  end_foreach_face_generic()
#line 16
 end_foreach_face(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 6
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 6
{

#line 6 "/home/fpl/softwares/basilisk/src/timestep.h"

    if (val(u.x,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.x,0,0,0));




      dt *= val_cm(cm,0,0,0);

      if (dt < dtmax) dtmax = dt;
    } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 6
{

#line 6 "/home/fpl/softwares/basilisk/src/timestep.h"

    if (val(u.y,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.y,0,0,0));




      dt *= val_cm(cm,0,0,0);

      if (dt < dtmax) dtmax = dt;
    } }  }}  end_foreach_face_generic()
#line 16
 end_foreach_face(); }mpi_all_reduce_array (&dtmax, double, MPI_MIN, 1);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 16
 }
  dtmax *= CFL;
  if (dtmax > previous)
    dtmax = (previous + 0.1*dtmax)/1.1;
  previous = dtmax;
  return dtmax;
}
#line 29 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
#line 1 "./bcg.h"
#line 1 "/home/fpl/softwares/basilisk/src/bcg.h"
#line 11 "/home/fpl/softwares/basilisk/src/bcg.h"
void tracer_fluxes (scalar f,
      vector uf,
      vector flux,
      double dt,
       scalar src)
{





  vector g= new_vector("g");
  gradients (((scalar []){f,{-1}}), ((vector []){{g.x,g.y},{{-1},{-1}}}));




   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dt = dt;
{ double dt = _dt; NOT_UNUSED(dt);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/bcg.h", .line = 28,
    .each = "foreach_face", .first = _first_call
  };

strongif (!is_constant(fm.x) && !is_constant(src)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_src
#define val_src(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_src
#define fine_src(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_src
#define coarse_src(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 28
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 28
{

#line 28 "/home/fpl/softwares/basilisk/src/bcg.h"
 {







    double un = dt*_stencil_val(__FILE__,__LINE__,uf.x,0,0,0)/(val_fm_x(fm.x,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = _stencil_val(__FILE__,__LINE__,f,i,0,0) + (val_src(src,0,0,0) + val_src(src,-1,0,0))*dt/4. + s*(1. - s*un)*_stencil_val(__FILE__,__LINE__,g.x,i,0,0)*Delta/2.;





    IF (val_fm_y(fm.y,i,0,0) && val_fm_y(fm.y,i,1,0)) {
      double vn = (_stencil_val(__FILE__,__LINE__,uf.y,i,0,0) + _stencil_val(__FILE__,__LINE__,uf.y,i,1,0))/(val_fm_y(fm.y,i,0,0) + val_fm_y(fm.y,i,1,0));
      double fyy = vn < 0. ? _stencil_val(__FILE__,__LINE__,f,i,1,0) - _stencil_val(__FILE__,__LINE__,f,i,0,0) : _stencil_val(__FILE__,__LINE__,f,i,0,0) - _stencil_val(__FILE__,__LINE__,f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/fpl/softwares/basilisk/src/bcg.h"
    _stencil_val(__FILE__,__LINE__,flux.x,0,0,0) = f2*_stencil_val(__FILE__,__LINE__,uf.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 28
{

#line 28 "/home/fpl/softwares/basilisk/src/bcg.h"
 {







    double un = dt*_stencil_val(__FILE__,__LINE__,uf.y,0,0,0)/(val_fm_y(fm.y,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = _stencil_val(__FILE__,__LINE__,f,0,i,0) + (val_src(src,0,0,0) + val_src(src,0,-1,0))*dt/4. + s*(1. - s*un)*_stencil_val(__FILE__,__LINE__,g.y,0,i,0)*Delta/2.;





    IF (val_fm_x(fm.x,0,i,0) && val_fm_x(fm.x,1,i,0)) {
      double vn = (_stencil_val(__FILE__,__LINE__,uf.x,0,i,0) + _stencil_val(__FILE__,__LINE__,uf.x,1,i,0))/(val_fm_x(fm.x,0,i,0) + val_fm_x(fm.x,1,i,0));
      double fyy = vn < 0. ? _stencil_val(__FILE__,__LINE__,f,1,i,0) - _stencil_val(__FILE__,__LINE__,f,0,i,0) : _stencil_val(__FILE__,__LINE__,f,0,i,0) - _stencil_val(__FILE__,__LINE__,f,-1,i,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/fpl/softwares/basilisk/src/bcg.h"
    _stencil_val(__FILE__,__LINE__,flux.y,0,0,0) = f2*_stencil_val(__FILE__,__LINE__,uf.y,0,0,0);
  } }  }}  end_foreach_face_stencil()
#line 59
 }
strongif (is_constant(fm.x) && !is_constant(src)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_src
#define val_src(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_src
#define fine_src(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_src
#define coarse_src(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 28
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 28
{

#line 28 "/home/fpl/softwares/basilisk/src/bcg.h"
 {







    double un = dt*_stencil_val(__FILE__,__LINE__,uf.x,0,0,0)/(val_fm_x(fm.x,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = _stencil_val(__FILE__,__LINE__,f,i,0,0) + (val_src(src,0,0,0) + val_src(src,-1,0,0))*dt/4. + s*(1. - s*un)*_stencil_val(__FILE__,__LINE__,g.x,i,0,0)*Delta/2.;





    IF (val_fm_y(fm.y,i,0,0) && val_fm_y(fm.y,i,1,0)) {
      double vn = (_stencil_val(__FILE__,__LINE__,uf.y,i,0,0) + _stencil_val(__FILE__,__LINE__,uf.y,i,1,0))/(val_fm_y(fm.y,i,0,0) + val_fm_y(fm.y,i,1,0));
      double fyy = vn < 0. ? _stencil_val(__FILE__,__LINE__,f,i,1,0) - _stencil_val(__FILE__,__LINE__,f,i,0,0) : _stencil_val(__FILE__,__LINE__,f,i,0,0) - _stencil_val(__FILE__,__LINE__,f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/fpl/softwares/basilisk/src/bcg.h"
    _stencil_val(__FILE__,__LINE__,flux.x,0,0,0) = f2*_stencil_val(__FILE__,__LINE__,uf.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 28
{

#line 28 "/home/fpl/softwares/basilisk/src/bcg.h"
 {







    double un = dt*_stencil_val(__FILE__,__LINE__,uf.y,0,0,0)/(val_fm_y(fm.y,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = _stencil_val(__FILE__,__LINE__,f,0,i,0) + (val_src(src,0,0,0) + val_src(src,0,-1,0))*dt/4. + s*(1. - s*un)*_stencil_val(__FILE__,__LINE__,g.y,0,i,0)*Delta/2.;





    IF (val_fm_x(fm.x,0,i,0) && val_fm_x(fm.x,1,i,0)) {
      double vn = (_stencil_val(__FILE__,__LINE__,uf.x,0,i,0) + _stencil_val(__FILE__,__LINE__,uf.x,1,i,0))/(val_fm_x(fm.x,0,i,0) + val_fm_x(fm.x,1,i,0));
      double fyy = vn < 0. ? _stencil_val(__FILE__,__LINE__,f,1,i,0) - _stencil_val(__FILE__,__LINE__,f,0,i,0) : _stencil_val(__FILE__,__LINE__,f,0,i,0) - _stencil_val(__FILE__,__LINE__,f,-1,i,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/fpl/softwares/basilisk/src/bcg.h"
    _stencil_val(__FILE__,__LINE__,flux.y,0,0,0) = f2*_stencil_val(__FILE__,__LINE__,uf.y,0,0,0);
  } }  }}  end_foreach_face_stencil()
#line 59
 }
strongif (!is_constant(fm.x) && is_constant(src)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
const double _const_src = _constant[src.i -_NVARMAX];
NOT_UNUSED(_const_src);
#undef val_src
#define val_src(a,i,j,k) _const_src
#undef fine_src
#define fine_src(a,i,j,k) _const_src
#undef coarse_src
#define coarse_src(a,i,j,k) _const_src
#line 28
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 28
{

#line 28 "/home/fpl/softwares/basilisk/src/bcg.h"
 {







    double un = dt*_stencil_val(__FILE__,__LINE__,uf.x,0,0,0)/(val_fm_x(fm.x,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = _stencil_val(__FILE__,__LINE__,f,i,0,0) + (val_src(src,0,0,0) + val_src(src,-1,0,0))*dt/4. + s*(1. - s*un)*_stencil_val(__FILE__,__LINE__,g.x,i,0,0)*Delta/2.;





    IF (val_fm_y(fm.y,i,0,0) && val_fm_y(fm.y,i,1,0)) {
      double vn = (_stencil_val(__FILE__,__LINE__,uf.y,i,0,0) + _stencil_val(__FILE__,__LINE__,uf.y,i,1,0))/(val_fm_y(fm.y,i,0,0) + val_fm_y(fm.y,i,1,0));
      double fyy = vn < 0. ? _stencil_val(__FILE__,__LINE__,f,i,1,0) - _stencil_val(__FILE__,__LINE__,f,i,0,0) : _stencil_val(__FILE__,__LINE__,f,i,0,0) - _stencil_val(__FILE__,__LINE__,f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/fpl/softwares/basilisk/src/bcg.h"
    _stencil_val(__FILE__,__LINE__,flux.x,0,0,0) = f2*_stencil_val(__FILE__,__LINE__,uf.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 28
{

#line 28 "/home/fpl/softwares/basilisk/src/bcg.h"
 {







    double un = dt*_stencil_val(__FILE__,__LINE__,uf.y,0,0,0)/(val_fm_y(fm.y,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = _stencil_val(__FILE__,__LINE__,f,0,i,0) + (val_src(src,0,0,0) + val_src(src,0,-1,0))*dt/4. + s*(1. - s*un)*_stencil_val(__FILE__,__LINE__,g.y,0,i,0)*Delta/2.;





    IF (val_fm_x(fm.x,0,i,0) && val_fm_x(fm.x,1,i,0)) {
      double vn = (_stencil_val(__FILE__,__LINE__,uf.x,0,i,0) + _stencil_val(__FILE__,__LINE__,uf.x,1,i,0))/(val_fm_x(fm.x,0,i,0) + val_fm_x(fm.x,1,i,0));
      double fyy = vn < 0. ? _stencil_val(__FILE__,__LINE__,f,1,i,0) - _stencil_val(__FILE__,__LINE__,f,0,i,0) : _stencil_val(__FILE__,__LINE__,f,0,i,0) - _stencil_val(__FILE__,__LINE__,f,-1,i,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/fpl/softwares/basilisk/src/bcg.h"
    _stencil_val(__FILE__,__LINE__,flux.y,0,0,0) = f2*_stencil_val(__FILE__,__LINE__,uf.y,0,0,0);
  } }  }}  end_foreach_face_stencil()
#line 59
 }
strongif (is_constant(fm.x) && is_constant(src)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_src = _constant[src.i -_NVARMAX];
NOT_UNUSED(_const_src);
#undef val_src
#define val_src(a,i,j,k) _const_src
#undef fine_src
#define fine_src(a,i,j,k) _const_src
#undef coarse_src
#define coarse_src(a,i,j,k) _const_src
#line 28
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 28
{

#line 28 "/home/fpl/softwares/basilisk/src/bcg.h"
 {







    double un = dt*_stencil_val(__FILE__,__LINE__,uf.x,0,0,0)/(val_fm_x(fm.x,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = _stencil_val(__FILE__,__LINE__,f,i,0,0) + (val_src(src,0,0,0) + val_src(src,-1,0,0))*dt/4. + s*(1. - s*un)*_stencil_val(__FILE__,__LINE__,g.x,i,0,0)*Delta/2.;





    IF (val_fm_y(fm.y,i,0,0) && val_fm_y(fm.y,i,1,0)) {
      double vn = (_stencil_val(__FILE__,__LINE__,uf.y,i,0,0) + _stencil_val(__FILE__,__LINE__,uf.y,i,1,0))/(val_fm_y(fm.y,i,0,0) + val_fm_y(fm.y,i,1,0));
      double fyy = vn < 0. ? _stencil_val(__FILE__,__LINE__,f,i,1,0) - _stencil_val(__FILE__,__LINE__,f,i,0,0) : _stencil_val(__FILE__,__LINE__,f,i,0,0) - _stencil_val(__FILE__,__LINE__,f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/fpl/softwares/basilisk/src/bcg.h"
    _stencil_val(__FILE__,__LINE__,flux.x,0,0,0) = f2*_stencil_val(__FILE__,__LINE__,uf.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 28
{

#line 28 "/home/fpl/softwares/basilisk/src/bcg.h"
 {







    double un = dt*_stencil_val(__FILE__,__LINE__,uf.y,0,0,0)/(val_fm_y(fm.y,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = _stencil_val(__FILE__,__LINE__,f,0,i,0) + (val_src(src,0,0,0) + val_src(src,0,-1,0))*dt/4. + s*(1. - s*un)*_stencil_val(__FILE__,__LINE__,g.y,0,i,0)*Delta/2.;





    IF (val_fm_x(fm.x,0,i,0) && val_fm_x(fm.x,1,i,0)) {
      double vn = (_stencil_val(__FILE__,__LINE__,uf.x,0,i,0) + _stencil_val(__FILE__,__LINE__,uf.x,1,i,0))/(val_fm_x(fm.x,0,i,0) + val_fm_x(fm.x,1,i,0));
      double fyy = vn < 0. ? _stencil_val(__FILE__,__LINE__,f,1,i,0) - _stencil_val(__FILE__,__LINE__,f,0,i,0) : _stencil_val(__FILE__,__LINE__,f,0,i,0) - _stencil_val(__FILE__,__LINE__,f,-1,i,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/fpl/softwares/basilisk/src/bcg.h"
    _stencil_val(__FILE__,__LINE__,flux.y,0,0,0) = f2*_stencil_val(__FILE__,__LINE__,uf.y,0,0,0);
  } }  }}  end_foreach_face_stencil()
#line 59
 } if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/fpl/softwares/basilisk/src/bcg.h", 28, "dt");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 59

strongif (!is_constant(fm.x) && !is_constant(src)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_src
#define val_src(a,i,j,k) val(a,i,j,k)
#undef fine_src
#define fine_src(a,i,j,k) fine(a,i,j,k)
#undef coarse_src
#define coarse_src(a,i,j,k) coarse(a,i,j,k)
#line 28
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 28
{

#line 28 "/home/fpl/softwares/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.x,0,0,0)/(val_fm_x(fm.x,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val_src(src,0,0,0) + val_src(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (val_fm_y(fm.y,i,0,0) && val_fm_y(fm.y,i,1,0)) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(val_fm_y(fm.y,i,0,0) + val_fm_y(fm.y,i,1,0));
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/fpl/softwares/basilisk/src/bcg.h"
    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 28
{

#line 28 "/home/fpl/softwares/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.y,0,0,0)/(val_fm_y(fm.y,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val_src(src,0,0,0) + val_src(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (val_fm_x(fm.x,0,i,0) && val_fm_x(fm.x,1,i,0)) {
      double vn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(val_fm_x(fm.x,0,i,0) + val_fm_x(fm.x,1,i,0));
      double fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/fpl/softwares/basilisk/src/bcg.h"
    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 59
 end_foreach_face(); }
strongif (is_constant(fm.x) && !is_constant(src)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_src
#define val_src(a,i,j,k) val(a,i,j,k)
#undef fine_src
#define fine_src(a,i,j,k) fine(a,i,j,k)
#undef coarse_src
#define coarse_src(a,i,j,k) coarse(a,i,j,k)
#line 28
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 28
{

#line 28 "/home/fpl/softwares/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.x,0,0,0)/(val_fm_x(fm.x,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val_src(src,0,0,0) + val_src(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (val_fm_y(fm.y,i,0,0) && val_fm_y(fm.y,i,1,0)) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(val_fm_y(fm.y,i,0,0) + val_fm_y(fm.y,i,1,0));
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/fpl/softwares/basilisk/src/bcg.h"
    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 28
{

#line 28 "/home/fpl/softwares/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.y,0,0,0)/(val_fm_y(fm.y,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val_src(src,0,0,0) + val_src(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (val_fm_x(fm.x,0,i,0) && val_fm_x(fm.x,1,i,0)) {
      double vn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(val_fm_x(fm.x,0,i,0) + val_fm_x(fm.x,1,i,0));
      double fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/fpl/softwares/basilisk/src/bcg.h"
    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 59
 end_foreach_face(); }
strongif (!is_constant(fm.x) && is_constant(src)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const double _const_src = _constant[src.i -_NVARMAX];
NOT_UNUSED(_const_src);
#undef val_src
#define val_src(a,i,j,k) _const_src
#undef fine_src
#define fine_src(a,i,j,k) _const_src
#undef coarse_src
#define coarse_src(a,i,j,k) _const_src
#line 28
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 28
{

#line 28 "/home/fpl/softwares/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.x,0,0,0)/(val_fm_x(fm.x,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val_src(src,0,0,0) + val_src(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (val_fm_y(fm.y,i,0,0) && val_fm_y(fm.y,i,1,0)) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(val_fm_y(fm.y,i,0,0) + val_fm_y(fm.y,i,1,0));
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/fpl/softwares/basilisk/src/bcg.h"
    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 28
{

#line 28 "/home/fpl/softwares/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.y,0,0,0)/(val_fm_y(fm.y,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val_src(src,0,0,0) + val_src(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (val_fm_x(fm.x,0,i,0) && val_fm_x(fm.x,1,i,0)) {
      double vn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(val_fm_x(fm.x,0,i,0) + val_fm_x(fm.x,1,i,0));
      double fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/fpl/softwares/basilisk/src/bcg.h"
    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 59
 end_foreach_face(); }
strongif (is_constant(fm.x) && is_constant(src)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_src = _constant[src.i -_NVARMAX];
NOT_UNUSED(_const_src);
#undef val_src
#define val_src(a,i,j,k) _const_src
#undef fine_src
#define fine_src(a,i,j,k) _const_src
#undef coarse_src
#define coarse_src(a,i,j,k) _const_src
#line 28
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 28
{

#line 28 "/home/fpl/softwares/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.x,0,0,0)/(val_fm_x(fm.x,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val_src(src,0,0,0) + val_src(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (val_fm_y(fm.y,i,0,0) && val_fm_y(fm.y,i,1,0)) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(val_fm_y(fm.y,i,0,0) + val_fm_y(fm.y,i,1,0));
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/fpl/softwares/basilisk/src/bcg.h"
    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 28
{

#line 28 "/home/fpl/softwares/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.y,0,0,0)/(val_fm_y(fm.y,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val_src(src,0,0,0) + val_src(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (val_fm_x(fm.x,0,i,0) && val_fm_x(fm.x,1,i,0)) {
      double vn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(val_fm_x(fm.x,0,i,0) + val_fm_x(fm.x,1,i,0));
      double fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/fpl/softwares/basilisk/src/bcg.h"
    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 59
 end_foreach_face(); } }
 delete (((scalar []){g.x,g.y,{-1}})); }






struct Advection {
  scalar * tracers;
  vector u;
  double dt;
  scalar * src;
};

void advection (struct Advection p)
{




  scalar * lsrc = p.src;
  if (!lsrc)
    strongif (p.tracers) for (scalar s = *p.tracers, *_i103 = p.tracers; ((scalar *)&s)->i >= 0; s = *++_i103)
      lsrc = list_append (lsrc, zeroc);
  if (!(list_len(p.tracers) == list_len(lsrc))) qassert ("/home/fpl/softwares/basilisk/src/bcg.h", 84, "list_len(p.tracers) == list_len(lsrc)");

  scalar f, src;
  scalar * _i2 = p.tracers; scalar * _i3 = lsrc; strongif (p.tracers) for (f = *p.tracers, src = *lsrc; ((scalar *)&f)->i >= 0; f = *++_i2, src = *++_i3) {
    vector flux= new_face_vector("flux");
    tracer_fluxes (f, p.u, flux, p.dt, src);

     { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/bcg.h", .line = 91,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 91
foreach_stencil(){

#line 91 "/home/fpl/softwares/basilisk/src/bcg.h"

      {
#line 92

        _stencil_val(__FILE__,__LINE__,f,0,0,0) += p.dt*(_stencil_val(__FILE__,__LINE__,flux.x,0,0,0) - _stencil_val(__FILE__,__LINE__,flux.x,1,0,0))/(Delta*val_cm(cm,0,0,0));
#line 92

        _stencil_val(__FILE__,__LINE__,f,0,0,0) += p.dt*(_stencil_val(__FILE__,__LINE__,flux.y,0,0,0) - _stencil_val(__FILE__,__LINE__,flux.y,0,1,0))/(Delta*val_cm(cm,0,0,0));}; } end_foreach_stencil(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 91
foreach_stencil(){

#line 91 "/home/fpl/softwares/basilisk/src/bcg.h"

      {
#line 92

        _stencil_val(__FILE__,__LINE__,f,0,0,0) += p.dt*(_stencil_val(__FILE__,__LINE__,flux.x,0,0,0) - _stencil_val(__FILE__,__LINE__,flux.x,1,0,0))/(Delta*val_cm(cm,0,0,0));
#line 92

        _stencil_val(__FILE__,__LINE__,f,0,0,0) += p.dt*(_stencil_val(__FILE__,__LINE__,flux.y,0,0,0) - _stencil_val(__FILE__,__LINE__,flux.y,0,1,0))/(Delta*val_cm(cm,0,0,0));}; } end_foreach_stencil(); }  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 93

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 91
foreach(){

#line 91 "/home/fpl/softwares/basilisk/src/bcg.h"

      {
#line 92

        val(f,0,0,0) += p.dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/(Delta*val_cm(cm,0,0,0));
#line 92

        val(f,0,0,0) += p.dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/(Delta*val_cm(cm,0,0,0));}; } end_foreach(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 91
foreach(){

#line 91 "/home/fpl/softwares/basilisk/src/bcg.h"

      {
#line 92

        val(f,0,0,0) += p.dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/(Delta*val_cm(cm,0,0,0));
#line 92

        val(f,0,0,0) += p.dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/(Delta*val_cm(cm,0,0,0));}; } end_foreach(); } }



   delete (((scalar []){flux.x,flux.y,{-1}})); }

  if (!p.src)
    pfree (lsrc,__func__,__FILE__,__LINE__);
}
#line 30 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"



#line 1 "./viscosity.h"
#line 1 "/home/fpl/softwares/basilisk/src/viscosity.h"
#line 1 "./poisson.h"
#line 1 "/home/fpl/softwares/basilisk/src/poisson.h"
#line 32 "/home/fpl/softwares/basilisk/src/poisson.h"
void mg_cycle (scalar * a, scalar * res, scalar * da,
        void (* relax) (scalar * da, scalar * res,
          int depth, void * data),
        void * data,
        int nrelax, int minlevel, int maxlevel)
{




  restriction (res);





  minlevel = min (minlevel, maxlevel);
  for (int l = minlevel; l <= maxlevel; l++) {




    if (l == minlevel)
       { foreach_level_or_leaf (l){

#line 55 "/home/fpl/softwares/basilisk/src/poisson.h"

 strongif (da) for (scalar s = *da, *_i104 = da; ((scalar *)&s)->i >= 0; s = *++_i104)
  
     val(s,0,0,0) = 0.; } end_foreach_level_or_leaf(); }





    else
       { foreach_level (l){

#line 65 "/home/fpl/softwares/basilisk/src/poisson.h"

 strongif (da) for (scalar s = *da, *_i105 = da; ((scalar *)&s)->i >= 0; s = *++_i105)
  
     val(s,0,0,0) = bilinear (point, s); } end_foreach_level(); }





    boundary_level (da, l);
    for (int i = 0; i < nrelax; i++) {
      relax (da, res, l, data);
      boundary_level (da, l);
    }
  }




   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/poisson.h", .line = 84,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 84 "/home/fpl/softwares/basilisk/src/poisson.h"
 {
    scalar s, ds;
    scalar * _i4 = a; scalar * _i5 = da; strongif (a) for (s = *a, ds = *da; ((scalar *)&s)->i >= 0; s = *++_i4, ds = *++_i5)
     
 _stencil_val(__FILE__,__LINE__,s,0,0,0) += _stencil_val(__FILE__,__LINE__,ds,0,0,0);
  } } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 89
foreach(){

#line 84 "/home/fpl/softwares/basilisk/src/poisson.h"
 {
    scalar s, ds;
    scalar * _i4 = a; scalar * _i5 = da; strongif (a) for (s = *a, ds = *da; ((scalar *)&s)->i >= 0; s = *++_i4, ds = *++_i5)
     
 val(s,0,0,0) += val(ds,0,0,0);
  } } end_foreach(); }
}
#line 102 "/home/fpl/softwares/basilisk/src/poisson.h"
int NITERMAX = 100, NITERMIN = 1;
double TOLERANCE = 1e-3;




typedef struct {
  int i;
  double resb, resa;
  double sum;
  int nrelax;
  int minlevel;
} mgstats;
#line 125 "/home/fpl/softwares/basilisk/src/poisson.h"
struct MGSolve {
  scalar * a, * b;
  double (* residual) (scalar * a, scalar * b, scalar * res,
         void * data);
  void (* relax) (scalar * da, scalar * res, int depth,
    void * data);
  void * data;

  int nrelax;
  scalar * res;
  int minlevel;
  double tolerance;
};

mgstats mg_solve (struct MGSolve p)
{





  scalar * da = list_clone (p.a), * res = p.res;
  if (!res)
    res = list_clone (p.b);






  for (int b = 0; b < nboundary; b++)
    strongif (da) for (scalar s = *da, *_i106 = da; ((scalar *)&s)->i >= 0; s = *++_i106)
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b];




  mgstats s = {0};
  double sum = 0.;
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _sum = sum;
{ double sum = _sum; NOT_UNUSED(sum);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/poisson.h", .line = 164,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 164 "/home/fpl/softwares/basilisk/src/poisson.h"

    strongif (p.b) for (scalar s = *p.b, *_i107 = p.b; ((scalar *)&s)->i >= 0; s = *++_i107)
      sum += _stencil_val(__FILE__,__LINE__,s,0,0,0); } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 166

#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:sum)) {

#line 164
foreach (){

#line 164 "/home/fpl/softwares/basilisk/src/poisson.h"

    strongif (p.b) for (scalar s = *p.b, *_i107 = p.b; ((scalar *)&s)->i >= 0; s = *++_i107)
      sum += val(s,0,0,0); } end_foreach();mpi_all_reduce_array (&sum, double, MPI_SUM, 1);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 166
 }
  s.sum = sum;
  s.nrelax = p.nrelax > 0 ? p.nrelax : 4;




  double resb;
  resb = s.resb = s.resa = p.residual (p.a, p.b, res, p.data);






  if (p.tolerance == 0.)
    p.tolerance = TOLERANCE;
  for (s.i = 0;
       s.i < NITERMAX && (s.i < NITERMIN || s.resa > p.tolerance);
       s.i++) {
    mg_cycle (p.a, res, da, p.relax, p.data,
       s.nrelax,
       p.minlevel,
       grid->maxdepth);
    s.resa = p.residual (p.a, p.b, res, p.data);
#line 199 "/home/fpl/softwares/basilisk/src/poisson.h"
    if (s.resa > p.tolerance) {
      if (resb/s.resa < 1.2 && s.nrelax < 100)
 s.nrelax++;
      else if (resb/s.resa > 10 && s.nrelax > 2)
 s.nrelax--;
    }







    resb = s.resa;
  }
  s.minlevel = p.minlevel;




  if (s.resa > p.tolerance) {
    scalar v = p.a[0];
    fprintf (ferr,
      "WARNING: convergence for %s not reached after %d iterations\n"
      "  res: %g sum: %g nrelax: %d\n", _attribute[v.i].name,
      s.i, s.resa, s.sum, s.nrelax), fflush (ferr);
  }




  if (!p.res)
    delete (res), pfree (res,__func__,__FILE__,__LINE__);
  delete (da), pfree (da,__func__,__FILE__,__LINE__);

  return s;
}
#line 258 "/home/fpl/softwares/basilisk/src/poisson.h"
struct Poisson {
  scalar a, b;
   vector alpha;
   scalar lambda;
  double tolerance;
  int nrelax, minlevel;
  scalar * res;
};





static void relax (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct Poisson * p = (struct Poisson *) data;
   vector alpha = p->alpha;
   scalar lambda = p->lambda;
#line 296 "/home/fpl/softwares/basilisk/src/poisson.h"
  scalar c = a;






   { 
strongif (!is_constant(lambda) && !is_constant(alpha.x)) {
#undef val_lambda
#define val_lambda(a,i,j,k) val(a,i,j,k)
#undef fine_lambda
#define fine_lambda(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 303
foreach_level_or_leaf (l){

#line 303 "/home/fpl/softwares/basilisk/src/poisson.h"
 {
    double n = - sq(Delta)*val(b,0,0,0), d = - val_lambda(lambda,0,0,0)*sq(Delta);
    {
#line 305
 {
      n += val_alpha_x(alpha.x,1,0,0)*val(a,1,0,0) + val_alpha_x(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val_alpha_x(alpha.x,1,0,0) + val_alpha_x(alpha.x,0,0,0);
    }
#line 305
 {
      n += val_alpha_y(alpha.y,0,1,0)*val(a,0,1,0) + val_alpha_y(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val_alpha_y(alpha.y,0,1,0) + val_alpha_y(alpha.y,0,0,0);
    }}
#line 319 "/home/fpl/softwares/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); }
strongif (is_constant(lambda) && !is_constant(alpha.x)) {
const double _const_lambda = _constant[lambda.i -_NVARMAX];
NOT_UNUSED(_const_lambda);
#undef val_lambda
#define val_lambda(a,i,j,k) _const_lambda
#undef fine_lambda
#define fine_lambda(a,i,j,k) _const_lambda
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _const_lambda
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 303
foreach_level_or_leaf (l){

#line 303 "/home/fpl/softwares/basilisk/src/poisson.h"
 {
    double n = - sq(Delta)*val(b,0,0,0), d = - val_lambda(lambda,0,0,0)*sq(Delta);
    {
#line 305
 {
      n += val_alpha_x(alpha.x,1,0,0)*val(a,1,0,0) + val_alpha_x(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val_alpha_x(alpha.x,1,0,0) + val_alpha_x(alpha.x,0,0,0);
    }
#line 305
 {
      n += val_alpha_y(alpha.y,0,1,0)*val(a,0,1,0) + val_alpha_y(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val_alpha_y(alpha.y,0,1,0) + val_alpha_y(alpha.y,0,0,0);
    }}
#line 319 "/home/fpl/softwares/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); }
strongif (!is_constant(lambda) && is_constant(alpha.x)) {
#undef val_lambda
#define val_lambda(a,i,j,k) val(a,i,j,k)
#undef fine_lambda
#define fine_lambda(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 303
foreach_level_or_leaf (l){

#line 303 "/home/fpl/softwares/basilisk/src/poisson.h"
 {
    double n = - sq(Delta)*val(b,0,0,0), d = - val_lambda(lambda,0,0,0)*sq(Delta);
    {
#line 305
 {
      n += val_alpha_x(alpha.x,1,0,0)*val(a,1,0,0) + val_alpha_x(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val_alpha_x(alpha.x,1,0,0) + val_alpha_x(alpha.x,0,0,0);
    }
#line 305
 {
      n += val_alpha_y(alpha.y,0,1,0)*val(a,0,1,0) + val_alpha_y(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val_alpha_y(alpha.y,0,1,0) + val_alpha_y(alpha.y,0,0,0);
    }}
#line 319 "/home/fpl/softwares/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); }
strongif (is_constant(lambda) && is_constant(alpha.x)) {
const double _const_lambda = _constant[lambda.i -_NVARMAX];
NOT_UNUSED(_const_lambda);
#undef val_lambda
#define val_lambda(a,i,j,k) _const_lambda
#undef fine_lambda
#define fine_lambda(a,i,j,k) _const_lambda
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _const_lambda
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 303
foreach_level_or_leaf (l){

#line 303 "/home/fpl/softwares/basilisk/src/poisson.h"
 {
    double n = - sq(Delta)*val(b,0,0,0), d = - val_lambda(lambda,0,0,0)*sq(Delta);
    {
#line 305
 {
      n += val_alpha_x(alpha.x,1,0,0)*val(a,1,0,0) + val_alpha_x(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val_alpha_x(alpha.x,1,0,0) + val_alpha_x(alpha.x,0,0,0);
    }
#line 305
 {
      n += val_alpha_y(alpha.y,0,1,0)*val(a,0,1,0) + val_alpha_y(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val_alpha_y(alpha.y,0,1,0) + val_alpha_y(alpha.y,0,0,0);
    }}
#line 319 "/home/fpl/softwares/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); } }
#line 338 "/home/fpl/softwares/basilisk/src/poisson.h"
}






static double residual (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct Poisson * p = (struct Poisson *) data;
   vector alpha = p->alpha;
   scalar lambda = p->lambda;



  double maxres = 0.;


  vector g= new_face_vector("g");
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/poisson.h", .line = 358,
    .each = "foreach_face", .first = _first_call
  };

strongif (!is_constant(alpha.x)) {
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 358
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 358
{

#line 358 "/home/fpl/softwares/basilisk/src/poisson.h"

    _stencil_val(__FILE__,__LINE__,g.x,0,0,0) = val_alpha_x(alpha.x,0,0,0)*((_stencil_val(__FILE__,__LINE__,a,0,0,0) - _stencil_val(__FILE__,__LINE__,a,0 -1,0,0))/Delta); }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 358
{

#line 358 "/home/fpl/softwares/basilisk/src/poisson.h"

    _stencil_val(__FILE__,__LINE__,g.y,0,0,0) = val_alpha_y(alpha.y,0,0,0)*((_stencil_val(__FILE__,__LINE__,a,0,0,0) - _stencil_val(__FILE__,__LINE__,a,0,0 -1,0))/Delta); }  }}  end_foreach_face_stencil()
#line 359
 }
strongif (is_constant(alpha.x)) {
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 358
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 358
{

#line 358 "/home/fpl/softwares/basilisk/src/poisson.h"

    _stencil_val(__FILE__,__LINE__,g.x,0,0,0) = val_alpha_x(alpha.x,0,0,0)*((_stencil_val(__FILE__,__LINE__,a,0,0,0) - _stencil_val(__FILE__,__LINE__,a,0 -1,0,0))/Delta); }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 358
{

#line 358 "/home/fpl/softwares/basilisk/src/poisson.h"

    _stencil_val(__FILE__,__LINE__,g.y,0,0,0) = val_alpha_y(alpha.y,0,0,0)*((_stencil_val(__FILE__,__LINE__,a,0,0,0) - _stencil_val(__FILE__,__LINE__,a,0,0 -1,0))/Delta); }  }}  end_foreach_face_stencil()
#line 359
 }  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 359

strongif (!is_constant(alpha.x)) {
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 358
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 358
{

#line 358 "/home/fpl/softwares/basilisk/src/poisson.h"

    val(g.x,0,0,0) = val_alpha_x(alpha.x,0,0,0)*((val(a,0,0,0) - val(a,0 -1,0,0))/Delta); }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 358
{

#line 358 "/home/fpl/softwares/basilisk/src/poisson.h"

    val(g.y,0,0,0) = val_alpha_y(alpha.y,0,0,0)*((val(a,0,0,0) - val(a,0,0 -1,0))/Delta); }  }}  end_foreach_face_generic()
#line 359
 end_foreach_face(); }
strongif (is_constant(alpha.x)) {
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 358
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 358
{

#line 358 "/home/fpl/softwares/basilisk/src/poisson.h"

    val(g.x,0,0,0) = val_alpha_x(alpha.x,0,0,0)*((val(a,0,0,0) - val(a,0 -1,0,0))/Delta); }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 358
{

#line 358 "/home/fpl/softwares/basilisk/src/poisson.h"

    val(g.y,0,0,0) = val_alpha_y(alpha.y,0,0,0)*((val(a,0,0,0) - val(a,0,0 -1,0))/Delta); }  }}  end_foreach_face_generic()
#line 359
 end_foreach_face(); } }
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _maxres = maxres;
{ double maxres = _maxres; NOT_UNUSED(maxres);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/poisson.h", .line = 360,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(lambda)) {
#undef val_lambda
#define val_lambda(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_lambda
#define fine_lambda(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 360
foreach_stencil(){

#line 360 "/home/fpl/softwares/basilisk/src/poisson.h"
 {
    _stencil_val(__FILE__,__LINE__,res,0,0,0) = _stencil_val(__FILE__,__LINE__,b,0,0,0) - val_lambda(lambda,0,0,0)*_stencil_val(__FILE__,__LINE__,a,0,0,0);
    {
#line 362

      _stencil_val(__FILE__,__LINE__,res,0,0,0) -= (_stencil_val(__FILE__,__LINE__,g.x,1,0,0) - _stencil_val(__FILE__,__LINE__,g.x,0,0,0))/Delta;
#line 362

      _stencil_val(__FILE__,__LINE__,res,0,0,0) -= (_stencil_val(__FILE__,__LINE__,g.y,0,1,0) - _stencil_val(__FILE__,__LINE__,g.y,0,0,0))/Delta;}






    IF (fabs (_stencil_val(__FILE__,__LINE__,res,0,0,0)) > maxres)
      maxres = fabs (_stencil_val(__FILE__,__LINE__,res,0,0,0));
  } } end_foreach_stencil(); }
strongif (is_constant(lambda)) {
const double _const_lambda = _constant[lambda.i -_NVARMAX];
NOT_UNUSED(_const_lambda);
#undef val_lambda
#define val_lambda(a,i,j,k) _const_lambda
#undef fine_lambda
#define fine_lambda(a,i,j,k) _const_lambda
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _const_lambda
#line 360
foreach_stencil(){

#line 360 "/home/fpl/softwares/basilisk/src/poisson.h"
 {
    _stencil_val(__FILE__,__LINE__,res,0,0,0) = _stencil_val(__FILE__,__LINE__,b,0,0,0) - val_lambda(lambda,0,0,0)*_stencil_val(__FILE__,__LINE__,a,0,0,0);
    {
#line 362

      _stencil_val(__FILE__,__LINE__,res,0,0,0) -= (_stencil_val(__FILE__,__LINE__,g.x,1,0,0) - _stencil_val(__FILE__,__LINE__,g.x,0,0,0))/Delta;
#line 362

      _stencil_val(__FILE__,__LINE__,res,0,0,0) -= (_stencil_val(__FILE__,__LINE__,g.y,0,1,0) - _stencil_val(__FILE__,__LINE__,g.y,0,0,0))/Delta;}






    IF (fabs (_stencil_val(__FILE__,__LINE__,res,0,0,0)) > maxres)
      maxres = fabs (_stencil_val(__FILE__,__LINE__,res,0,0,0));
  } } end_foreach_stencil(); }  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 372

#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)) {

#line 360

strongif (!is_constant(lambda)) {
#undef val_lambda
#define val_lambda(a,i,j,k) val(a,i,j,k)
#undef fine_lambda
#define fine_lambda(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) coarse(a,i,j,k)
#line 360
foreach (){

#line 360 "/home/fpl/softwares/basilisk/src/poisson.h"
 {
    val(res,0,0,0) = val(b,0,0,0) - val_lambda(lambda,0,0,0)*val(a,0,0,0);
    {
#line 362

      val(res,0,0,0) -= (val(g.x,1,0,0) - val(g.x,0,0,0))/Delta;
#line 362

      val(res,0,0,0) -= (val(g.y,0,1,0) - val(g.y,0,0,0))/Delta;}






    if (fabs (val(res,0,0,0)) > maxres)
      maxres = fabs (val(res,0,0,0));
  } } end_foreach(); }
strongif (is_constant(lambda)) {
const double _const_lambda = _constant[lambda.i -_NVARMAX];
NOT_UNUSED(_const_lambda);
#undef val_lambda
#define val_lambda(a,i,j,k) _const_lambda
#undef fine_lambda
#define fine_lambda(a,i,j,k) _const_lambda
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _const_lambda
#line 360
foreach (){

#line 360 "/home/fpl/softwares/basilisk/src/poisson.h"
 {
    val(res,0,0,0) = val(b,0,0,0) - val_lambda(lambda,0,0,0)*val(a,0,0,0);
    {
#line 362

      val(res,0,0,0) -= (val(g.x,1,0,0) - val(g.x,0,0,0))/Delta;
#line 362

      val(res,0,0,0) -= (val(g.y,0,1,0) - val(g.y,0,0,0))/Delta;}






    if (fabs (val(res,0,0,0)) > maxres)
      maxres = fabs (val(res,0,0,0));
  } } end_foreach(); }mpi_all_reduce_array (&maxres, double, MPI_MAX, 1);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 372
 }
#line 390 "/home/fpl/softwares/basilisk/src/poisson.h"
  { double _ret =  maxres; delete (((scalar []){g.x,g.y,{-1}}));  return _ret; }
 delete (((scalar []){g.x,g.y,{-1}})); }
#line 402 "/home/fpl/softwares/basilisk/src/poisson.h"
mgstats poisson (struct Poisson p)
{






  if (!p.alpha.x.i)
    p.alpha = unityf;
  if (!p.lambda.i)
    p.lambda = zeroc;




  vector alpha = p.alpha;
  scalar lambda = p.lambda;
  restriction (((scalar []){alpha.x,alpha.y,lambda,{-1}}));





  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;

  scalar a = p.a, b = p.b;
  mgstats s = mg_solve ((struct MGSolve){((scalar []){a,{-1}}), ((scalar []){b,{-1}}), residual, relax,
   &p, p.nrelax, p.res, .minlevel = max(1, p.minlevel)});




  if (p.tolerance)
    TOLERANCE = defaultol;

  return s;
}
#line 460 "/home/fpl/softwares/basilisk/src/poisson.h"
struct Project {
  vector uf;
  scalar p;
  vector alpha;
  double dt;
  int nrelax;
};


mgstats project (struct Project q)
{ trace ("project", "/home/fpl/softwares/basilisk/src/poisson.h", 470);
  vector uf = q.uf;
  scalar p = q.p;
   vector alpha = q.alpha.x.i ? q.alpha : unityf;
  double dt = q.dt ? q.dt : 1.;
  int nrelax = q.nrelax ? q.nrelax : 4;






  scalar div= new_scalar("div");
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dt = dt;
{ double dt = _dt; NOT_UNUSED(dt);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/poisson.h", .line = 483,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 483 "/home/fpl/softwares/basilisk/src/poisson.h"
 {
    _stencil_val(__FILE__,__LINE__,div,0,0,0) = 0.;
    {
#line 485

      _stencil_val(__FILE__,__LINE__,div,0,0,0) += _stencil_val(__FILE__,__LINE__,uf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,uf.x,0,0,0);
#line 485

      _stencil_val(__FILE__,__LINE__,div,0,0,0) += _stencil_val(__FILE__,__LINE__,uf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,uf.y,0,0,0);}
    _stencil_val(__FILE__,__LINE__,div,0,0,0) /= dt*Delta;
  } } end_foreach_stencil(); if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/fpl/softwares/basilisk/src/poisson.h", 483, "dt");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 488
foreach(){

#line 483 "/home/fpl/softwares/basilisk/src/poisson.h"
 {
    val(div,0,0,0) = 0.;
    {
#line 485

      val(div,0,0,0) += val(uf.x,1,0,0) - val(uf.x,0,0,0);
#line 485

      val(div,0,0,0) += val(uf.y,0,1,0) - val(uf.y,0,0,0);}
    val(div,0,0,0) /= dt*Delta;
  } } end_foreach(); }
#line 499 "/home/fpl/softwares/basilisk/src/poisson.h"
  mgstats mgp = poisson ((struct Poisson){p, div, alpha,
    .tolerance = TOLERANCE/sq(dt), .nrelax = nrelax});




   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dt = dt;
{ double dt = _dt; NOT_UNUSED(dt);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/poisson.h", .line = 505,
    .each = "foreach_face", .first = _first_call
  };

strongif (!is_constant(alpha.x)) {
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 505
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 505
{

#line 505 "/home/fpl/softwares/basilisk/src/poisson.h"

    _stencil_val(__FILE__,__LINE__,uf.x,0,0,0) -= dt*val_alpha_x(alpha.x,0,0,0)*((_stencil_val(__FILE__,__LINE__,p,0,0,0) - _stencil_val(__FILE__,__LINE__,p,0 -1,0,0))/Delta); }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 505
{

#line 505 "/home/fpl/softwares/basilisk/src/poisson.h"

    _stencil_val(__FILE__,__LINE__,uf.y,0,0,0) -= dt*val_alpha_y(alpha.y,0,0,0)*((_stencil_val(__FILE__,__LINE__,p,0,0,0) - _stencil_val(__FILE__,__LINE__,p,0,0 -1,0))/Delta); }  }}  end_foreach_face_stencil()
#line 506
 }
strongif (is_constant(alpha.x)) {
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 505
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 505
{

#line 505 "/home/fpl/softwares/basilisk/src/poisson.h"

    _stencil_val(__FILE__,__LINE__,uf.x,0,0,0) -= dt*val_alpha_x(alpha.x,0,0,0)*((_stencil_val(__FILE__,__LINE__,p,0,0,0) - _stencil_val(__FILE__,__LINE__,p,0 -1,0,0))/Delta); }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 505
{

#line 505 "/home/fpl/softwares/basilisk/src/poisson.h"

    _stencil_val(__FILE__,__LINE__,uf.y,0,0,0) -= dt*val_alpha_y(alpha.y,0,0,0)*((_stencil_val(__FILE__,__LINE__,p,0,0,0) - _stencil_val(__FILE__,__LINE__,p,0,0 -1,0))/Delta); }  }}  end_foreach_face_stencil()
#line 506
 } if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/fpl/softwares/basilisk/src/poisson.h", 505, "dt");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 506

strongif (!is_constant(alpha.x)) {
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 505
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 505
{

#line 505 "/home/fpl/softwares/basilisk/src/poisson.h"

    val(uf.x,0,0,0) -= dt*val_alpha_x(alpha.x,0,0,0)*((val(p,0,0,0) - val(p,0 -1,0,0))/Delta); }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 505
{

#line 505 "/home/fpl/softwares/basilisk/src/poisson.h"

    val(uf.y,0,0,0) -= dt*val_alpha_y(alpha.y,0,0,0)*((val(p,0,0,0) - val(p,0,0 -1,0))/Delta); }  }}  end_foreach_face_generic()
#line 506
 end_foreach_face(); }
strongif (is_constant(alpha.x)) {
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 505
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 505
{

#line 505 "/home/fpl/softwares/basilisk/src/poisson.h"

    val(uf.x,0,0,0) -= dt*val_alpha_x(alpha.x,0,0,0)*((val(p,0,0,0) - val(p,0 -1,0,0))/Delta); }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 505
{

#line 505 "/home/fpl/softwares/basilisk/src/poisson.h"

    val(uf.y,0,0,0) -= dt*val_alpha_y(alpha.y,0,0,0)*((val(p,0,0,0) - val(p,0,0 -1,0))/Delta); }  }}  end_foreach_face_generic()
#line 506
 end_foreach_face(); } }

  { mgstats _ret =  mgp; delete (((scalar []){div,{-1}}));  end_trace("project", "/home/fpl/softwares/basilisk/src/poisson.h", 508);  return _ret; }
 delete (((scalar []){div,{-1}}));  end_trace("project", "/home/fpl/softwares/basilisk/src/poisson.h", 509); }
#line 2 "/home/fpl/softwares/basilisk/src/viscosity.h"

struct Viscosity {
  vector u;
  vector mu;
  scalar rho;
  double dt;
  int nrelax;
  scalar * res;
};
#line 25 "/home/fpl/softwares/basilisk/src/viscosity.h"
static void relax_viscosity (scalar * a, scalar * b, int l, void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
   vector mu = p->mu;
   scalar rho = p->rho;
  double dt = p->dt;
  vector u = (*((vector *)&(a[0]))), r = (*((vector *)&(b[0])));




  vector w = u;


   { 
strongif (!is_constant(rho) && !is_constant(mu.x)) {
#undef val_rho
#define val_rho(a,i,j,k) val(a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#line 39
foreach_level_or_leaf (l){

#line 39 "/home/fpl/softwares/basilisk/src/viscosity.h"
 {
    {
#line 40

      val(w.x,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val_mu_y(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)
#line 58 "/home/fpl/softwares/basilisk/src/viscosity.h"
      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).x + dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)

          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)




        ));
#line 40

      val(w.y,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val_mu_x(mu.x,1,0,0)*(val(u.y,1,0,0) +
     (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
     (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)
#line 58 "/home/fpl/softwares/basilisk/src/viscosity.h"
      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).y + dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)

          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)




        ));}
  } } end_foreach_level_or_leaf(); }
strongif (is_constant(rho) && !is_constant(mu.x)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,i,j,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#line 39
foreach_level_or_leaf (l){

#line 39 "/home/fpl/softwares/basilisk/src/viscosity.h"
 {
    {
#line 40

      val(w.x,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val_mu_y(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)
#line 58 "/home/fpl/softwares/basilisk/src/viscosity.h"
      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).x + dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)

          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)




        ));
#line 40

      val(w.y,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val_mu_x(mu.x,1,0,0)*(val(u.y,1,0,0) +
     (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
     (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)
#line 58 "/home/fpl/softwares/basilisk/src/viscosity.h"
      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).y + dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)

          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)




        ));}
  } } end_foreach_level_or_leaf(); }
strongif (!is_constant(rho) && is_constant(mu.x)) {
#undef val_rho
#define val_rho(a,i,j,k) val(a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#line 39
foreach_level_or_leaf (l){

#line 39 "/home/fpl/softwares/basilisk/src/viscosity.h"
 {
    {
#line 40

      val(w.x,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val_mu_y(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)
#line 58 "/home/fpl/softwares/basilisk/src/viscosity.h"
      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).x + dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)

          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)




        ));
#line 40

      val(w.y,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val_mu_x(mu.x,1,0,0)*(val(u.y,1,0,0) +
     (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
     (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)
#line 58 "/home/fpl/softwares/basilisk/src/viscosity.h"
      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).y + dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)

          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)




        ));}
  } } end_foreach_level_or_leaf(); }
strongif (is_constant(rho) && is_constant(mu.x)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,i,j,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
const struct { double x, y; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#line 39
foreach_level_or_leaf (l){

#line 39 "/home/fpl/softwares/basilisk/src/viscosity.h"
 {
    {
#line 40

      val(w.x,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val_mu_y(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)
#line 58 "/home/fpl/softwares/basilisk/src/viscosity.h"
      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).x + dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)

          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)




        ));
#line 40

      val(w.y,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val_mu_x(mu.x,1,0,0)*(val(u.y,1,0,0) +
     (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
     (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)
#line 58 "/home/fpl/softwares/basilisk/src/viscosity.h"
      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).y + dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)

          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)




        ));}
  } } end_foreach_level_or_leaf(); } }
#line 85 "/home/fpl/softwares/basilisk/src/viscosity.h"
}

static double residual_viscosity (scalar * a, scalar * b, scalar * resl,
      void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
   vector mu = p->mu;
   scalar rho = p->rho;
  double dt = p->dt;
  vector u = (*((vector *)&(a[0]))), r = (*((vector *)&(b[0]))), res = (*((vector *)&(resl[0])));
  double maxres = 0.;
#line 104 "/home/fpl/softwares/basilisk/src/viscosity.h"
  boundary_internal ((scalar *)(((vector []){{u.x,u.y},{{-1},{-1}}})), "/home/fpl/softwares/basilisk/src/viscosity.h", 104);

  {
#line 106
 {
    vector taux= new_face_vector("taux");
     { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/viscosity.h", .line = 108,
    .each = "foreach_face", .first = _first_call
  };

strongif (!is_constant(mu.x)) {
#undef val_mu_x
#define val_mu_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 108
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 108
{

#line 108 "/home/fpl/softwares/basilisk/src/viscosity.h"

      _stencil_val(__FILE__,__LINE__,taux.x,0,0,0) = 2.*val_mu_x(mu.x,0,0,0)*(_stencil_val(__FILE__,__LINE__,u.x,0,0,0) - _stencil_val(__FILE__,__LINE__,u.x,-1,0,0))/Delta; }  }}  end_foreach_face_stencil()
#line 109
 }
strongif (is_constant(mu.x)) {
const struct { double x, y; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#line 108
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 108
{

#line 108 "/home/fpl/softwares/basilisk/src/viscosity.h"

      _stencil_val(__FILE__,__LINE__,taux.x,0,0,0) = 2.*val_mu_x(mu.x,0,0,0)*(_stencil_val(__FILE__,__LINE__,u.x,0,0,0) - _stencil_val(__FILE__,__LINE__,u.x,-1,0,0))/Delta; }  }}  end_foreach_face_stencil()
#line 109
 }  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 109

strongif (!is_constant(mu.x)) {
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#line 108
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 108
{

#line 108 "/home/fpl/softwares/basilisk/src/viscosity.h"

      val(taux.x,0,0,0) = 2.*val_mu_x(mu.x,0,0,0)*(val(u.x,0,0,0) - val(u.x,-1,0,0))/Delta; }  }}  end_foreach_face_generic()
#line 109
 end_foreach_face(); }
strongif (is_constant(mu.x)) {
const struct { double x, y; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#line 108
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 108
{

#line 108 "/home/fpl/softwares/basilisk/src/viscosity.h"

      val(taux.x,0,0,0) = 2.*val_mu_x(mu.x,0,0,0)*(val(u.x,0,0,0) - val(u.x,-1,0,0))/Delta; }  }}  end_foreach_face_generic()
#line 109
 end_foreach_face(); } }

       { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/viscosity.h", .line = 111,
    .each = "foreach_face", .first = _first_call
  };

strongif (!is_constant(mu.x)) {
#undef val_mu_x
#define val_mu_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 111
foreach_face_stencil() { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 111
{

#line 111 "/home/fpl/softwares/basilisk/src/viscosity.h"

 _stencil_val(__FILE__,__LINE__,taux.y,0,0,0) = val_mu_y(mu.y,0,0,0)*(_stencil_val(__FILE__,__LINE__,u.x,0,0,0) - _stencil_val(__FILE__,__LINE__,u.x,0,-1,0) +
      (_stencil_val(__FILE__,__LINE__,u.y,1,-1,0) + _stencil_val(__FILE__,__LINE__,u.y,1,0,0))/4. -
      (_stencil_val(__FILE__,__LINE__,u.y,-1,-1,0) + _stencil_val(__FILE__,__LINE__,u.y,-1,0,0))/4.)/Delta; }  }}  end_foreach_face_stencil()
#line 114
 }
strongif (is_constant(mu.x)) {
const struct { double x, y; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#line 111
foreach_face_stencil() { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 111
{

#line 111 "/home/fpl/softwares/basilisk/src/viscosity.h"

 _stencil_val(__FILE__,__LINE__,taux.y,0,0,0) = val_mu_y(mu.y,0,0,0)*(_stencil_val(__FILE__,__LINE__,u.x,0,0,0) - _stencil_val(__FILE__,__LINE__,u.x,0,-1,0) +
      (_stencil_val(__FILE__,__LINE__,u.y,1,-1,0) + _stencil_val(__FILE__,__LINE__,u.y,1,0,0))/4. -
      (_stencil_val(__FILE__,__LINE__,u.y,-1,-1,0) + _stencil_val(__FILE__,__LINE__,u.y,-1,0,0))/4.)/Delta; }  }}  end_foreach_face_stencil()
#line 114
 }  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 114

strongif (!is_constant(mu.x)) {
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#line 111
foreach_face_generic() { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 111
{

#line 111 "/home/fpl/softwares/basilisk/src/viscosity.h"

 val(taux.y,0,0,0) = val_mu_y(mu.y,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,-1,0) +
      (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
      (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 114
 end_foreach_face(); }
strongif (is_constant(mu.x)) {
const struct { double x, y; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#line 111
foreach_face_generic() { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 111
{

#line 111 "/home/fpl/softwares/basilisk/src/viscosity.h"

 val(taux.y,0,0,0) = val_mu_y(mu.y,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,-1,0) +
      (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
      (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 114
 end_foreach_face(); } }







     { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dt = dt;
 double _maxres = maxres;
{ double dt = _dt; NOT_UNUSED(dt);
 double maxres = _maxres; NOT_UNUSED(maxres);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/viscosity.h", .line = 122,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(rho)) {
#undef val_rho
#define val_rho(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 122
foreach_stencil(){

#line 122 "/home/fpl/softwares/basilisk/src/viscosity.h"
 {
      double d = 0.;
      {
#line 124

 d += _stencil_val(__FILE__,__LINE__,taux.x,1,0,0) - _stencil_val(__FILE__,__LINE__,taux.x,0,0,0);
#line 124

 d += _stencil_val(__FILE__,__LINE__,taux.y,0,1,0) - _stencil_val(__FILE__,__LINE__,taux.y,0,0,0);}
      _stencil_val(__FILE__,__LINE__,res.x,0,0,0) = _stencil_val(__FILE__,__LINE__,r.x,0,0,0) - ((coord){1.,1.}).x*_stencil_val(__FILE__,__LINE__,u.x,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      IF (fabs (_stencil_val(__FILE__,__LINE__,res.x,0,0,0)) > maxres)
 maxres = fabs (_stencil_val(__FILE__,__LINE__,res.x,0,0,0));
    } } end_foreach_stencil(); }
strongif (is_constant(rho)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,i,j,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#line 122
foreach_stencil(){

#line 122 "/home/fpl/softwares/basilisk/src/viscosity.h"
 {
      double d = 0.;
      {
#line 124

 d += _stencil_val(__FILE__,__LINE__,taux.x,1,0,0) - _stencil_val(__FILE__,__LINE__,taux.x,0,0,0);
#line 124

 d += _stencil_val(__FILE__,__LINE__,taux.y,0,1,0) - _stencil_val(__FILE__,__LINE__,taux.y,0,0,0);}
      _stencil_val(__FILE__,__LINE__,res.x,0,0,0) = _stencil_val(__FILE__,__LINE__,r.x,0,0,0) - ((coord){1.,1.}).x*_stencil_val(__FILE__,__LINE__,u.x,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      IF (fabs (_stencil_val(__FILE__,__LINE__,res.x,0,0,0)) > maxres)
 maxres = fabs (_stencil_val(__FILE__,__LINE__,res.x,0,0,0));
    } } end_foreach_stencil(); } if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/fpl/softwares/basilisk/src/viscosity.h", 122, "dt");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 129

#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)) {

#line 122

strongif (!is_constant(rho)) {
#undef val_rho
#define val_rho(a,i,j,k) val(a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,i,j,k)
#line 122
foreach (){

#line 122 "/home/fpl/softwares/basilisk/src/viscosity.h"
 {
      double d = 0.;
      {
#line 124

 d += val(taux.x,1,0,0) - val(taux.x,0,0,0);
#line 124

 d += val(taux.y,0,1,0) - val(taux.y,0,0,0);}
      val(res.x,0,0,0) = val(r.x,0,0,0) - ((coord){1.,1.}).x*val(u.x,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      if (fabs (val(res.x,0,0,0)) > maxres)
 maxres = fabs (val(res.x,0,0,0));
    } } end_foreach(); }
strongif (is_constant(rho)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,i,j,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#line 122
foreach (){

#line 122 "/home/fpl/softwares/basilisk/src/viscosity.h"
 {
      double d = 0.;
      {
#line 124

 d += val(taux.x,1,0,0) - val(taux.x,0,0,0);
#line 124

 d += val(taux.y,0,1,0) - val(taux.y,0,0,0);}
      val(res.x,0,0,0) = val(r.x,0,0,0) - ((coord){1.,1.}).x*val(u.x,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      if (fabs (val(res.x,0,0,0)) > maxres)
 maxres = fabs (val(res.x,0,0,0));
    } } end_foreach(); }mpi_all_reduce_array (&maxres, double, MPI_MAX, 1);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 129
 }
   delete (((scalar []){taux.x,taux.y,{-1}})); }
#line 106
 {
    vector taux= new_face_vector("taux");
     { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/viscosity.h", .line = 108,
    .each = "foreach_face", .first = _first_call
  };

strongif (!is_constant(mu.y)) {
#undef val_mu_y
#define val_mu_y(a,j,i,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_mu_x
#define val_mu_x(a,j,i,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 108
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 108
{

#line 108 "/home/fpl/softwares/basilisk/src/viscosity.h"

      _stencil_val(__FILE__,__LINE__,taux.y,0,0,0) = 2.*val_mu_y(mu.y,0,0,0)*(_stencil_val(__FILE__,__LINE__,u.y,0,0,0) - _stencil_val(__FILE__,__LINE__,u.y,-1,0,0))/Delta; }  }}  end_foreach_face_stencil()
#line 109
 }
strongif (is_constant(mu.y)) {
const struct { double x, y; } _const_mu = {_constant[mu.y.i -_NVARMAX], _constant[mu.x.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_y
#define val_mu_y(a,j,i,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_x
#define val_mu_x(a,j,i,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#line 108
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 108
{

#line 108 "/home/fpl/softwares/basilisk/src/viscosity.h"

      _stencil_val(__FILE__,__LINE__,taux.y,0,0,0) = 2.*val_mu_y(mu.y,0,0,0)*(_stencil_val(__FILE__,__LINE__,u.y,0,0,0) - _stencil_val(__FILE__,__LINE__,u.y,-1,0,0))/Delta; }  }}  end_foreach_face_stencil()
#line 109
 }  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 109

strongif (!is_constant(mu.y)) {
#undef val_mu_y
#define val_mu_y(a,j,i,k) val(a,j,i,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,j,i,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,j,i,k)
#undef val_mu_x
#define val_mu_x(a,j,i,k) val(a,j,i,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,j,i,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,j,i,k)
#line 108
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_y()) {
#line 108
{

#line 108 "/home/fpl/softwares/basilisk/src/viscosity.h"

      val(taux.y,0,0,0) = 2.*val_mu_y(mu.y,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 109
 end_foreach_face(); }
strongif (is_constant(mu.y)) {
const struct { double x, y; } _const_mu = {_constant[mu.y.i -_NVARMAX], _constant[mu.x.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_y
#define val_mu_y(a,j,i,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_x
#define val_mu_x(a,j,i,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#line 108
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_y()) {
#line 108
{

#line 108 "/home/fpl/softwares/basilisk/src/viscosity.h"

      val(taux.y,0,0,0) = 2.*val_mu_y(mu.y,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 109
 end_foreach_face(); } }

       { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/viscosity.h", .line = 111,
    .each = "foreach_face", .first = _first_call
  };

strongif (!is_constant(mu.y)) {
#undef val_mu_y
#define val_mu_y(a,j,i,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_mu_x
#define val_mu_x(a,j,i,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 111
foreach_face_stencil() { int jg = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 111
{

#line 111 "/home/fpl/softwares/basilisk/src/viscosity.h"

 _stencil_val(__FILE__,__LINE__,taux.x,0,0,0) = val_mu_x(mu.x,0,0,0)*(_stencil_val(__FILE__,__LINE__,u.y,0,0,0) - _stencil_val(__FILE__,__LINE__,u.y,0,-1,0) +
      (_stencil_val(__FILE__,__LINE__,u.x,1,-1,0) + _stencil_val(__FILE__,__LINE__,u.x,1,0,0))/4. -
      (_stencil_val(__FILE__,__LINE__,u.x,-1,-1,0) + _stencil_val(__FILE__,__LINE__,u.x,-1,0,0))/4.)/Delta; }  }}  end_foreach_face_stencil()
#line 114
 }
strongif (is_constant(mu.y)) {
const struct { double x, y; } _const_mu = {_constant[mu.y.i -_NVARMAX], _constant[mu.x.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_y
#define val_mu_y(a,j,i,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_x
#define val_mu_x(a,j,i,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#line 111
foreach_face_stencil() { int jg = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 111
{

#line 111 "/home/fpl/softwares/basilisk/src/viscosity.h"

 _stencil_val(__FILE__,__LINE__,taux.x,0,0,0) = val_mu_x(mu.x,0,0,0)*(_stencil_val(__FILE__,__LINE__,u.y,0,0,0) - _stencil_val(__FILE__,__LINE__,u.y,0,-1,0) +
      (_stencil_val(__FILE__,__LINE__,u.x,1,-1,0) + _stencil_val(__FILE__,__LINE__,u.x,1,0,0))/4. -
      (_stencil_val(__FILE__,__LINE__,u.x,-1,-1,0) + _stencil_val(__FILE__,__LINE__,u.x,-1,0,0))/4.)/Delta; }  }}  end_foreach_face_stencil()
#line 114
 }  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 114

strongif (!is_constant(mu.y)) {
#undef val_mu_y
#define val_mu_y(a,j,i,k) val(a,j,i,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,j,i,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,j,i,k)
#undef val_mu_x
#define val_mu_x(a,j,i,k) val(a,j,i,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,j,i,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,j,i,k)
#line 111
foreach_face_generic() { int jg = -1; VARIABLES;  strongif (is_face_x()) {
#line 111
{

#line 111 "/home/fpl/softwares/basilisk/src/viscosity.h"

 val(taux.x,0,0,0) = val_mu_x(mu.x,0,0,0)*(val(u.y,0,0,0) - val(u.y,-1,0,0) +
      (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
      (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 114
 end_foreach_face(); }
strongif (is_constant(mu.y)) {
const struct { double x, y; } _const_mu = {_constant[mu.y.i -_NVARMAX], _constant[mu.x.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_y
#define val_mu_y(a,j,i,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_x
#define val_mu_x(a,j,i,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#line 111
foreach_face_generic() { int jg = -1; VARIABLES;  strongif (is_face_x()) {
#line 111
{

#line 111 "/home/fpl/softwares/basilisk/src/viscosity.h"

 val(taux.x,0,0,0) = val_mu_x(mu.x,0,0,0)*(val(u.y,0,0,0) - val(u.y,-1,0,0) +
      (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
      (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 114
 end_foreach_face(); } }







     { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dt = dt;
 double _maxres = maxres;
{ double dt = _dt; NOT_UNUSED(dt);
 double maxres = _maxres; NOT_UNUSED(maxres);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/viscosity.h", .line = 122,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(rho)) {
#undef val_rho
#define val_rho(a,j,i,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 122
foreach_stencil(){

#line 122 "/home/fpl/softwares/basilisk/src/viscosity.h"
 {
      double d = 0.;
      {
#line 124

 d += _stencil_val(__FILE__,__LINE__,taux.y,1,0,0) - _stencil_val(__FILE__,__LINE__,taux.y,0,0,0);
#line 124

 d += _stencil_val(__FILE__,__LINE__,taux.x,0,1,0) - _stencil_val(__FILE__,__LINE__,taux.x,0,0,0);}
      _stencil_val(__FILE__,__LINE__,res.y,0,0,0) = _stencil_val(__FILE__,__LINE__,r.y,0,0,0) - ((coord){1.,1.}).y*_stencil_val(__FILE__,__LINE__,u.y,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      IF (fabs (_stencil_val(__FILE__,__LINE__,res.y,0,0,0)) > maxres)
 maxres = fabs (_stencil_val(__FILE__,__LINE__,res.y,0,0,0));
    } } end_foreach_stencil(); }
strongif (is_constant(rho)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,j,i,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#line 122
foreach_stencil(){

#line 122 "/home/fpl/softwares/basilisk/src/viscosity.h"
 {
      double d = 0.;
      {
#line 124

 d += _stencil_val(__FILE__,__LINE__,taux.y,1,0,0) - _stencil_val(__FILE__,__LINE__,taux.y,0,0,0);
#line 124

 d += _stencil_val(__FILE__,__LINE__,taux.x,0,1,0) - _stencil_val(__FILE__,__LINE__,taux.x,0,0,0);}
      _stencil_val(__FILE__,__LINE__,res.y,0,0,0) = _stencil_val(__FILE__,__LINE__,r.y,0,0,0) - ((coord){1.,1.}).y*_stencil_val(__FILE__,__LINE__,u.y,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      IF (fabs (_stencil_val(__FILE__,__LINE__,res.y,0,0,0)) > maxres)
 maxres = fabs (_stencil_val(__FILE__,__LINE__,res.y,0,0,0));
    } } end_foreach_stencil(); } if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/fpl/softwares/basilisk/src/viscosity.h", 122, "dt");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 129

#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)) {

#line 122

strongif (!is_constant(rho)) {
#undef val_rho
#define val_rho(a,j,i,k) val(a,j,i,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,j,i,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,j,i,k)
#line 122
foreach (){

#line 122 "/home/fpl/softwares/basilisk/src/viscosity.h"
 {
      double d = 0.;
      {
#line 124

 d += val(taux.y,0,1,0) - val(taux.y,0,0,0);
#line 124

 d += val(taux.x,1,0,0) - val(taux.x,0,0,0);}
      val(res.y,0,0,0) = val(r.y,0,0,0) - ((coord){1.,1.}).y*val(u.y,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      if (fabs (val(res.y,0,0,0)) > maxres)
 maxres = fabs (val(res.y,0,0,0));
    } } end_foreach(); }
strongif (is_constant(rho)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,j,i,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#line 122
foreach (){

#line 122 "/home/fpl/softwares/basilisk/src/viscosity.h"
 {
      double d = 0.;
      {
#line 124

 d += val(taux.y,0,1,0) - val(taux.y,0,0,0);
#line 124

 d += val(taux.x,1,0,0) - val(taux.x,0,0,0);}
      val(res.y,0,0,0) = val(r.y,0,0,0) - ((coord){1.,1.}).y*val(u.y,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      if (fabs (val(res.y,0,0,0)) > maxres)
 maxres = fabs (val(res.y,0,0,0));
    } } end_foreach(); }mpi_all_reduce_array (&maxres, double, MPI_MAX, 1);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 129
 }
   delete (((scalar []){taux.x,taux.y,{-1}})); }}
#line 159 "/home/fpl/softwares/basilisk/src/viscosity.h"
  return maxres;
}




mgstats viscosity (struct Viscosity p)
{ trace ("viscosity", "/home/fpl/softwares/basilisk/src/viscosity.h", 166);
  vector u = p.u, r= new_vector("r");
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/viscosity.h", .line = 168,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 168 "/home/fpl/softwares/basilisk/src/viscosity.h"

    {
#line 169

      _stencil_val(__FILE__,__LINE__,r.x,0,0,0) = _stencil_val(__FILE__,__LINE__,u.x,0,0,0);
#line 169

      _stencil_val(__FILE__,__LINE__,r.y,0,0,0) = _stencil_val(__FILE__,__LINE__,u.y,0,0,0);}; } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 170
foreach(){

#line 168 "/home/fpl/softwares/basilisk/src/viscosity.h"

    {
#line 169

      val(r.x,0,0,0) = val(u.x,0,0,0);
#line 169

      val(r.y,0,0,0) = val(u.y,0,0,0);}; } end_foreach(); }

  vector mu = p.mu;
  scalar rho = p.rho;
  restriction (((scalar []){mu.x,mu.y,rho,{-1}}));

  { mgstats _ret =  mg_solve ((struct MGSolve){(scalar *)((vector []){{u.x,u.y},{{-1},{-1}}}), (scalar *)((vector []){{r.x,r.y},{{-1},{-1}}}),
     residual_viscosity, relax_viscosity, &p, p.nrelax, p.res}); delete (((scalar []){r.x,r.y,{-1}}));  end_trace("viscosity", "/home/fpl/softwares/basilisk/src/viscosity.h", 177);  return _ret; }
 delete (((scalar []){r.x,r.y,{-1}}));  end_trace("viscosity", "/home/fpl/softwares/basilisk/src/viscosity.h", 178); }


mgstats viscosity_explicit (struct Viscosity p)
{ trace ("viscosity_explicit", "/home/fpl/softwares/basilisk/src/viscosity.h", 182);
  vector u = p.u, r= new_vector("r");
  mgstats mg = {0};
  mg.resb = residual_viscosity ((scalar *)((vector []){{u.x,u.y},{{-1},{-1}}}), (scalar *)((vector []){{u.x,u.y},{{-1},{-1}}}), (scalar *)((vector []){{r.x,r.y},{{-1},{-1}}}), &p);
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/viscosity.h", .line = 186,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 186 "/home/fpl/softwares/basilisk/src/viscosity.h"

    {
#line 187

      _stencil_val(__FILE__,__LINE__,u.x,0,0,0) += _stencil_val(__FILE__,__LINE__,r.x,0,0,0);
#line 187

      _stencil_val(__FILE__,__LINE__,u.y,0,0,0) += _stencil_val(__FILE__,__LINE__,r.y,0,0,0);}; } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 188
foreach(){

#line 186 "/home/fpl/softwares/basilisk/src/viscosity.h"

    {
#line 187

      val(u.x,0,0,0) += val(r.x,0,0,0);
#line 187

      val(u.y,0,0,0) += val(r.y,0,0,0);}; } end_foreach(); }
  { mgstats _ret =  mg; delete (((scalar []){r.x,r.y,{-1}}));  end_trace("viscosity_explicit", "/home/fpl/softwares/basilisk/src/viscosity.h", 189);  return _ret; }
 delete (((scalar []){r.x,r.y,{-1}}));  end_trace("viscosity_explicit", "/home/fpl/softwares/basilisk/src/viscosity.h", 190); }
#line 34 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
#line 44 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
scalar p= {0};
vector u= {{1},{2}}, g= {{3},{4}};
scalar pf= {5};
vector uf= {{6},{7}};
#line 70 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
 vector mu = {{_NVARMAX + 0},{_NVARMAX + 1}}, a = {{_NVARMAX + 0},{_NVARMAX + 1}}, alpha = {{_NVARMAX + 2},{_NVARMAX + 3}};
 scalar rho = {(_NVARMAX + 4)};
mgstats mgp, mgpf, mgu;
bool stokes = false;
#line 91 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
static void _set_boundary0 (void) { _attribute[p.i].boundary[right] = _boundary0; _attribute[p.i].boundary_homogeneous[right] = _boundary0_homogeneous; _attribute[p.i].dirty = true; } 
static void _set_boundary1 (void) { _attribute[p.i].boundary[left] = _boundary1; _attribute[p.i].boundary_homogeneous[left] = _boundary1_homogeneous; _attribute[p.i].dirty = true; } 
#line 101 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
static void _set_boundary2 (void) { _attribute[p.i].boundary[top] = _boundary2; _attribute[p.i].boundary_homogeneous[top] = _boundary2_homogeneous; _attribute[p.i].dirty = true; } 
static void _set_boundary3 (void) { _attribute[p.i].boundary[bottom] = _boundary3; _attribute[p.i].boundary_homogeneous[bottom] = _boundary3_homogeneous; _attribute[p.i].dirty = true; } 
#line 126 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
static int defaults_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults_0 (const int i, const double t, Event * _ev) { trace ("defaults_0", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 126); 
{

  CFL = 0.8;




  _attribute[p.i].nodump = _attribute[pf.i].nodump = true;




  if (alpha.x.i == unityf.x.i) {
    alpha = fm;
    rho = cm;
  }
  else if (!is_constant(alpha.x)) {
    vector alphav = alpha;
     { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", .line = 145,
    .each = "foreach_face", .first = _first_call
  };

strongif (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 145
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 145
{

#line 145 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

      _stencil_val(__FILE__,__LINE__,alphav.x,0,0,0) = val_fm_x(fm.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 145
{

#line 145 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

      _stencil_val(__FILE__,__LINE__,alphav.y,0,0,0) = val_fm_y(fm.y,0,0,0); }  }}  end_foreach_face_stencil()
#line 146
 }
strongif (is_constant(fm.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 145
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 145
{

#line 145 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

      _stencil_val(__FILE__,__LINE__,alphav.x,0,0,0) = val_fm_x(fm.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 145
{

#line 145 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

      _stencil_val(__FILE__,__LINE__,alphav.y,0,0,0) = val_fm_y(fm.y,0,0,0); }  }}  end_foreach_face_stencil()
#line 146
 }  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 146

strongif (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#line 145
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 145
{

#line 145 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

      val(alphav.x,0,0,0) = val_fm_x(fm.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 145
{

#line 145 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

      val(alphav.y,0,0,0) = val_fm_y(fm.y,0,0,0); }  }}  end_foreach_face_generic()
#line 146
 end_foreach_face(); }
strongif (is_constant(fm.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 145
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 145
{

#line 145 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

      val(alphav.x,0,0,0) = val_fm_x(fm.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 145
{

#line 145 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

      val(alphav.y,0,0,0) = val_fm_y(fm.y,0,0,0); }  }}  end_foreach_face_generic()
#line 146
 end_foreach_face(); } }
  }






  _attribute[uf.x.i].refine = refine_face_solenoidal;
#line 173 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
 end_trace("defaults_0", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 173); } return 0; } 





static int default_display_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int default_display (const int i, const double t, Event * _ev) { trace ("default_display", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 179); 
  display ((struct _display){"squares (color = 'u.x', spread = -1);"}); end_trace("default_display", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 180);  return 0; } 





double dtmax;

static int init_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int init (const int i, const double t, Event * _ev) { trace ("init", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 188); 
{
  trash (((vector []){{uf.x,uf.y},{{-1},{-1}}}));
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", .line = 191,
    .each = "foreach_face", .first = _first_call
  };

strongif (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 191
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 191
{

#line 191 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    _stencil_val(__FILE__,__LINE__,uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*((_stencil_val(__FILE__,__LINE__,u.x,0,0,0) + _stencil_val(__FILE__,__LINE__,u.x,0 -1,0,0))/2.); }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 191
{

#line 191 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    _stencil_val(__FILE__,__LINE__,uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*((_stencil_val(__FILE__,__LINE__,u.y,0,0,0) + _stencil_val(__FILE__,__LINE__,u.y,0,0 -1,0))/2.); }  }}  end_foreach_face_stencil()
#line 192
 }
strongif (is_constant(fm.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 191
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 191
{

#line 191 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    _stencil_val(__FILE__,__LINE__,uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*((_stencil_val(__FILE__,__LINE__,u.x,0,0,0) + _stencil_val(__FILE__,__LINE__,u.x,0 -1,0,0))/2.); }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 191
{

#line 191 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    _stencil_val(__FILE__,__LINE__,uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*((_stencil_val(__FILE__,__LINE__,u.y,0,0,0) + _stencil_val(__FILE__,__LINE__,u.y,0,0 -1,0))/2.); }  }}  end_foreach_face_stencil()
#line 192
 }  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 192

strongif (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#line 191
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 191
{

#line 191 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.); }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 191
{

#line 191 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.); }  }}  end_foreach_face_generic()
#line 192
 end_foreach_face(); }
strongif (is_constant(fm.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 191
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 191
{

#line 191 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.); }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 191
{

#line 191 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.); }  }}  end_foreach_face_generic()
#line 192
 end_foreach_face(); } }




  event ("properties");





  dtmax = DT;
  event ("stability");
 end_trace("init", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 205); } return 0; } 
#line 214 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
static int set_dtmax_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int set_dtmax (const int i, const double t, Event * _ev) { trace ("set_dtmax", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 214);  dtmax = DT; end_trace("set_dtmax", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 214);  return 0; } 

static int stability_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int stability (const int i, const double t, Event * _ev) { trace ("stability", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 216);  {
  dt = dtnext (stokes ? dtmax : timestep (uf, dtmax));
 end_trace("stability", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 218); } return 0; } 







static int vof_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int vof (const int i, const double t, Event * _ev) { trace ("vof", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 226); ; end_trace("vof", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 226);  return 0; } 
static int tracer_advection_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int tracer_advection (const int i, const double t, Event * _ev) { trace ("tracer_advection", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 227); ; end_trace("tracer_advection", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 227);  return 0; } 
static int tracer_diffusion_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int tracer_diffusion (const int i, const double t, Event * _ev) { trace ("tracer_diffusion", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 228); ; end_trace("tracer_diffusion", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 228);  return 0; } 






static int properties_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int properties (const int i, const double t, Event * _ev) { trace ("properties", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 235); ; end_trace("properties", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 235);  return 0; } 
#line 247 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
void prediction ()
{
  vector du;
  {
#line 250
 {
    scalar s = new_scalar("s");
    du.x = s;
  }
#line 250
 {
    scalar s = new_scalar("s");
    du.y = s;
  }}

  if (_attribute[u.x.i].gradient)
     { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", .line = 256,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 256 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

      {
#line 257
 {





   _stencil_val(__FILE__,__LINE__,du.x,0,0,0) = _attribute[u.x.i].gradient (_stencil_val(__FILE__,__LINE__,u.x,-1,0,0), _stencil_val(__FILE__,__LINE__,u.x,0,0,0), _stencil_val(__FILE__,__LINE__,u.x,1,0,0))/Delta;
      }
#line 257
 {





   _stencil_val(__FILE__,__LINE__,du.y,0,0,0) = _attribute[u.y.i].gradient (_stencil_val(__FILE__,__LINE__,u.y,0,-1,0), _stencil_val(__FILE__,__LINE__,u.y,0,0,0), _stencil_val(__FILE__,__LINE__,u.y,0,1,0))/Delta;
      }} } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 264
foreach(){

#line 256 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

      {
#line 257
 {





   val(du.x,0,0,0) = _attribute[u.x.i].gradient (val(u.x,-1,0,0), val(u.x,0,0,0), val(u.x,1,0,0))/Delta;
      }
#line 257
 {





   val(du.y,0,0,0) = _attribute[u.y.i].gradient (val(u.y,0,-1,0), val(u.y,0,0,0), val(u.y,0,1,0))/Delta;
      }} } end_foreach(); }
  else
     { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", .line = 266,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 266 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

      {
#line 267
 {





   _stencil_val(__FILE__,__LINE__,du.x,0,0,0) = (_stencil_val(__FILE__,__LINE__,u.x,1,0,0) - _stencil_val(__FILE__,__LINE__,u.x,-1,0,0))/(2.*Delta);
    }
#line 267
 {





   _stencil_val(__FILE__,__LINE__,du.y,0,0,0) = (_stencil_val(__FILE__,__LINE__,u.y,0,1,0) - _stencil_val(__FILE__,__LINE__,u.y,0,-1,0))/(2.*Delta);
    }} } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 274
foreach(){

#line 266 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

      {
#line 267
 {





   val(du.x,0,0,0) = (val(u.x,1,0,0) - val(u.x,-1,0,0))/(2.*Delta);
    }
#line 267
 {





   val(du.y,0,0,0) = (val(u.y,0,1,0) - val(u.y,0,-1,0))/(2.*Delta);
    }} } end_foreach(); }

  trash (((vector []){{uf.x,uf.y},{{-1},{-1}}}));
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dt = dt;
{ double dt = _dt; NOT_UNUSED(dt);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", .line = 277,
    .each = "foreach_face", .first = _first_call
  };

strongif (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 277
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 277
{

#line 277 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(_stencil_val(__FILE__,__LINE__,u.x,0,0,0) + _stencil_val(__FILE__,__LINE__,u.x,-1,0,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    _stencil_val(__FILE__,__LINE__,uf.x,0,0,0) = _stencil_val(__FILE__,__LINE__,u.x,i,0,0) + (_stencil_val(__FILE__,__LINE__,g.x,0,0,0) + _stencil_val(__FILE__,__LINE__,g.x,-1,0,0))*dt/4. + s*(1. - s*un)*_stencil_val(__FILE__,__LINE__,du.x,i,0,0)*Delta/2.;

    IF (val_fm_y(fm.y,i,0,0) && val_fm_y(fm.y,i,1,0)) {
      double fyy = _stencil_val(__FILE__,__LINE__,u.y,i,0,0) < 0. ? _stencil_val(__FILE__,__LINE__,u.x,i,1,0) - _stencil_val(__FILE__,__LINE__,u.x,i,0,0) : _stencil_val(__FILE__,__LINE__,u.x,i,0,0) - _stencil_val(__FILE__,__LINE__,u.x,i,-1,0);
      _stencil_val(__FILE__,__LINE__,uf.x,0,0,0) -= dt*_stencil_val(__FILE__,__LINE__,u.y,i,0,0)*fyy/(2.*Delta);
    }







    _stencil_val(__FILE__,__LINE__,uf.x,0,0,0) *= val_fm_x(fm.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 277
{

#line 277 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(_stencil_val(__FILE__,__LINE__,u.y,0,0,0) + _stencil_val(__FILE__,__LINE__,u.y,0,-1,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    _stencil_val(__FILE__,__LINE__,uf.y,0,0,0) = _stencil_val(__FILE__,__LINE__,u.y,0,i,0) + (_stencil_val(__FILE__,__LINE__,g.y,0,0,0) + _stencil_val(__FILE__,__LINE__,g.y,0,-1,0))*dt/4. + s*(1. - s*un)*_stencil_val(__FILE__,__LINE__,du.y,0,i,0)*Delta/2.;

    IF (val_fm_x(fm.x,0,i,0) && val_fm_x(fm.x,1,i,0)) {
      double fyy = _stencil_val(__FILE__,__LINE__,u.x,0,i,0) < 0. ? _stencil_val(__FILE__,__LINE__,u.y,1,i,0) - _stencil_val(__FILE__,__LINE__,u.y,0,i,0) : _stencil_val(__FILE__,__LINE__,u.y,0,i,0) - _stencil_val(__FILE__,__LINE__,u.y,-1,i,0);
      _stencil_val(__FILE__,__LINE__,uf.y,0,0,0) -= dt*_stencil_val(__FILE__,__LINE__,u.x,0,i,0)*fyy/(2.*Delta);
    }







    _stencil_val(__FILE__,__LINE__,uf.y,0,0,0) *= val_fm_y(fm.y,0,0,0);
  } }  }}  end_foreach_face_stencil()
#line 294
 }
strongif (is_constant(fm.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 277
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 277
{

#line 277 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(_stencil_val(__FILE__,__LINE__,u.x,0,0,0) + _stencil_val(__FILE__,__LINE__,u.x,-1,0,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    _stencil_val(__FILE__,__LINE__,uf.x,0,0,0) = _stencil_val(__FILE__,__LINE__,u.x,i,0,0) + (_stencil_val(__FILE__,__LINE__,g.x,0,0,0) + _stencil_val(__FILE__,__LINE__,g.x,-1,0,0))*dt/4. + s*(1. - s*un)*_stencil_val(__FILE__,__LINE__,du.x,i,0,0)*Delta/2.;

    IF (val_fm_y(fm.y,i,0,0) && val_fm_y(fm.y,i,1,0)) {
      double fyy = _stencil_val(__FILE__,__LINE__,u.y,i,0,0) < 0. ? _stencil_val(__FILE__,__LINE__,u.x,i,1,0) - _stencil_val(__FILE__,__LINE__,u.x,i,0,0) : _stencil_val(__FILE__,__LINE__,u.x,i,0,0) - _stencil_val(__FILE__,__LINE__,u.x,i,-1,0);
      _stencil_val(__FILE__,__LINE__,uf.x,0,0,0) -= dt*_stencil_val(__FILE__,__LINE__,u.y,i,0,0)*fyy/(2.*Delta);
    }







    _stencil_val(__FILE__,__LINE__,uf.x,0,0,0) *= val_fm_x(fm.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 277
{

#line 277 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(_stencil_val(__FILE__,__LINE__,u.y,0,0,0) + _stencil_val(__FILE__,__LINE__,u.y,0,-1,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    _stencil_val(__FILE__,__LINE__,uf.y,0,0,0) = _stencil_val(__FILE__,__LINE__,u.y,0,i,0) + (_stencil_val(__FILE__,__LINE__,g.y,0,0,0) + _stencil_val(__FILE__,__LINE__,g.y,0,-1,0))*dt/4. + s*(1. - s*un)*_stencil_val(__FILE__,__LINE__,du.y,0,i,0)*Delta/2.;

    IF (val_fm_x(fm.x,0,i,0) && val_fm_x(fm.x,1,i,0)) {
      double fyy = _stencil_val(__FILE__,__LINE__,u.x,0,i,0) < 0. ? _stencil_val(__FILE__,__LINE__,u.y,1,i,0) - _stencil_val(__FILE__,__LINE__,u.y,0,i,0) : _stencil_val(__FILE__,__LINE__,u.y,0,i,0) - _stencil_val(__FILE__,__LINE__,u.y,-1,i,0);
      _stencil_val(__FILE__,__LINE__,uf.y,0,0,0) -= dt*_stencil_val(__FILE__,__LINE__,u.x,0,i,0)*fyy/(2.*Delta);
    }







    _stencil_val(__FILE__,__LINE__,uf.y,0,0,0) *= val_fm_y(fm.y,0,0,0);
  } }  }}  end_foreach_face_stencil()
#line 294
 } if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 277, "dt");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 294

strongif (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#line 277
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 277
{

#line 277 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. + s*(1. - s*un)*val(du.x,i,0,0)*Delta/2.;

    if (val_fm_y(fm.y,i,0,0) && val_fm_y(fm.y,i,1,0)) {
      double fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);
      val(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);
    }







    val(uf.x,0,0,0) *= val_fm_x(fm.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 277
{

#line 277 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. + s*(1. - s*un)*val(du.y,0,i,0)*Delta/2.;

    if (val_fm_x(fm.x,0,i,0) && val_fm_x(fm.x,1,i,0)) {
      double fyy = val(u.x,0,i,0) < 0. ? val(u.y,1,i,0) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,-1,i,0);
      val(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fyy/(2.*Delta);
    }







    val(uf.y,0,0,0) *= val_fm_y(fm.y,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 294
 end_foreach_face(); }
strongif (is_constant(fm.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 277
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 277
{

#line 277 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. + s*(1. - s*un)*val(du.x,i,0,0)*Delta/2.;

    if (val_fm_y(fm.y,i,0,0) && val_fm_y(fm.y,i,1,0)) {
      double fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);
      val(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);
    }







    val(uf.x,0,0,0) *= val_fm_x(fm.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 277
{

#line 277 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. + s*(1. - s*un)*val(du.y,0,i,0)*Delta/2.;

    if (val_fm_x(fm.x,0,i,0) && val_fm_x(fm.x,1,i,0)) {
      double fyy = val(u.x,0,i,0) < 0. ? val(u.y,1,i,0) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,-1,i,0);
      val(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fyy/(2.*Delta);
    }







    val(uf.y,0,0,0) *= val_fm_y(fm.y,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 294
 end_foreach_face(); } }

  delete ((scalar *)((vector []){{du.x,du.y},{{-1},{-1}}}));
}
#line 308 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
static int advection_term_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int advection_term (const int i, const double t, Event * _ev) { trace ("advection_term", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 308); 
{
  if (!stokes) {
    prediction();
    mgpf = project ((struct Project){uf, pf, alpha, dt/2., mgpf.nrelax});
    advection ((struct Advection){(scalar *)((vector []){{u.x,u.y},{{-1},{-1}}}), uf, dt, (scalar *)((vector []){{g.x,g.y},{{-1},{-1}}})});
  }
 end_trace("advection_term", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 315); } return 0; } 







static void correction (double dt)
{
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dt = dt;
{ double dt = _dt; NOT_UNUSED(dt);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", .line = 325,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 325 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    {
#line 326

      _stencil_val(__FILE__,__LINE__,u.x,0,0,0) += dt*_stencil_val(__FILE__,__LINE__,g.x,0,0,0);
#line 326

      _stencil_val(__FILE__,__LINE__,u.y,0,0,0) += dt*_stencil_val(__FILE__,__LINE__,g.y,0,0,0);}; } end_foreach_stencil(); if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 325, "dt");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 327
foreach(){

#line 325 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    {
#line 326

      val(u.x,0,0,0) += dt*val(g.x,0,0,0);
#line 326

      val(u.y,0,0,0) += dt*val(g.y,0,0,0);}; } end_foreach(); }
}
#line 337 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
static int viscous_term_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int viscous_term (const int i, const double t, Event * _ev) { trace ("viscous_term", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 337); 
{
  if (constant(mu.x) != 0.) {
    correction (dt);
    mgu = viscosity ((struct Viscosity){u, mu, rho, dt, mgu.nrelax});
    correction (-dt);
  }




  if (!is_constant(a.x)) {
    vector af = a;
    trash (((vector []){{af.x,af.y},{{-1},{-1}}}));
     { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", .line = 351,
    .each = "foreach_face", .first = _first_call
  };
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 351
{

#line 351 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

      _stencil_val(__FILE__,__LINE__,af.x,0,0,0) = 0.; }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 351
{

#line 351 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

      _stencil_val(__FILE__,__LINE__,af.y,0,0,0) = 0.; }  }}  end_foreach_face_stencil()
#line 352
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 352
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 351
{

#line 351 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

      val(af.x,0,0,0) = 0.; }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 351
{

#line 351 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

      val(af.y,0,0,0) = 0.; }  }}  end_foreach_face_generic()
#line 352
 end_foreach_face(); }
  }
 end_trace("viscous_term", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 354); } return 0; } 
#line 373 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
static int acceleration_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int acceleration (const int i, const double t, Event * _ev) { trace ("acceleration", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 373); 
{
  trash (((vector []){{uf.x,uf.y},{{-1},{-1}}}));
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dt = dt;
{ double dt = _dt; NOT_UNUSED(dt);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", .line = 376,
    .each = "foreach_face", .first = _first_call
  };

strongif (!is_constant(fm.x) && !is_constant(a.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_a_x
#define val_a_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 376
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 376
{

#line 376 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    _stencil_val(__FILE__,__LINE__,uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(((_stencil_val(__FILE__,__LINE__,u.x,0,0,0) + _stencil_val(__FILE__,__LINE__,u.x,0 -1,0,0))/2.) + dt*val_a_x(a.x,0,0,0)); }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 376
{

#line 376 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    _stencil_val(__FILE__,__LINE__,uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(((_stencil_val(__FILE__,__LINE__,u.y,0,0,0) + _stencil_val(__FILE__,__LINE__,u.y,0,0 -1,0))/2.) + dt*val_a_y(a.y,0,0,0)); }  }}  end_foreach_face_stencil()
#line 377
 }
strongif (is_constant(fm.x) && !is_constant(a.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_a_x
#define val_a_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 376
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 376
{

#line 376 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    _stencil_val(__FILE__,__LINE__,uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(((_stencil_val(__FILE__,__LINE__,u.x,0,0,0) + _stencil_val(__FILE__,__LINE__,u.x,0 -1,0,0))/2.) + dt*val_a_x(a.x,0,0,0)); }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 376
{

#line 376 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    _stencil_val(__FILE__,__LINE__,uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(((_stencil_val(__FILE__,__LINE__,u.y,0,0,0) + _stencil_val(__FILE__,__LINE__,u.y,0,0 -1,0))/2.) + dt*val_a_y(a.y,0,0,0)); }  }}  end_foreach_face_stencil()
#line 377
 }
strongif (!is_constant(fm.x) && is_constant(a.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#line 376
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 376
{

#line 376 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    _stencil_val(__FILE__,__LINE__,uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(((_stencil_val(__FILE__,__LINE__,u.x,0,0,0) + _stencil_val(__FILE__,__LINE__,u.x,0 -1,0,0))/2.) + dt*val_a_x(a.x,0,0,0)); }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 376
{

#line 376 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    _stencil_val(__FILE__,__LINE__,uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(((_stencil_val(__FILE__,__LINE__,u.y,0,0,0) + _stencil_val(__FILE__,__LINE__,u.y,0,0 -1,0))/2.) + dt*val_a_y(a.y,0,0,0)); }  }}  end_foreach_face_stencil()
#line 377
 }
strongif (is_constant(fm.x) && is_constant(a.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#line 376
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 376
{

#line 376 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    _stencil_val(__FILE__,__LINE__,uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(((_stencil_val(__FILE__,__LINE__,u.x,0,0,0) + _stencil_val(__FILE__,__LINE__,u.x,0 -1,0,0))/2.) + dt*val_a_x(a.x,0,0,0)); }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 376
{

#line 376 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    _stencil_val(__FILE__,__LINE__,uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(((_stencil_val(__FILE__,__LINE__,u.y,0,0,0) + _stencil_val(__FILE__,__LINE__,u.y,0,0 -1,0))/2.) + dt*val_a_y(a.y,0,0,0)); }  }}  end_foreach_face_stencil()
#line 377
 } if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 376, "dt");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 377

strongif (!is_constant(fm.x) && !is_constant(a.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#line 376
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 376
{

#line 376 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*val_a_x(a.x,0,0,0)); }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 376
{

#line 376 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*val_a_y(a.y,0,0,0)); }  }}  end_foreach_face_generic()
#line 377
 end_foreach_face(); }
strongif (is_constant(fm.x) && !is_constant(a.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#line 376
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 376
{

#line 376 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*val_a_x(a.x,0,0,0)); }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 376
{

#line 376 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*val_a_y(a.y,0,0,0)); }  }}  end_foreach_face_generic()
#line 377
 end_foreach_face(); }
strongif (!is_constant(fm.x) && is_constant(a.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#line 376
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 376
{

#line 376 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*val_a_x(a.x,0,0,0)); }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 376
{

#line 376 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*val_a_y(a.y,0,0,0)); }  }}  end_foreach_face_generic()
#line 377
 end_foreach_face(); }
strongif (is_constant(fm.x) && is_constant(a.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#line 376
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 376
{

#line 376 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*val_a_x(a.x,0,0,0)); }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 376
{

#line 376 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*val_a_y(a.y,0,0,0)); }  }}  end_foreach_face_generic()
#line 377
 end_foreach_face(); } }
 end_trace("acceleration", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 378); } return 0; } 
#line 387 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
void centered_gradient (scalar p, vector g)
{





  vector gf= new_face_vector("gf");
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", .line = 395,
    .each = "foreach_face", .first = _first_call
  };

strongif (!is_constant(fm.x) && !is_constant(a.x) && !is_constant(alpha.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_a_x
#define val_a_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 395
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    _stencil_val(__FILE__,__LINE__,gf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)*(_stencil_val(__FILE__,__LINE__,p,0,0,0) - _stencil_val(__FILE__,__LINE__,p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    _stencil_val(__FILE__,__LINE__,gf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)*(_stencil_val(__FILE__,__LINE__,p,0,0,0) - _stencil_val(__FILE__,__LINE__,p,0,-1,0))/Delta; }  }}  end_foreach_face_stencil()
#line 396
 }
strongif (is_constant(fm.x) && !is_constant(a.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_a_x
#define val_a_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 395
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    _stencil_val(__FILE__,__LINE__,gf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)*(_stencil_val(__FILE__,__LINE__,p,0,0,0) - _stencil_val(__FILE__,__LINE__,p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    _stencil_val(__FILE__,__LINE__,gf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)*(_stencil_val(__FILE__,__LINE__,p,0,0,0) - _stencil_val(__FILE__,__LINE__,p,0,-1,0))/Delta; }  }}  end_foreach_face_stencil()
#line 396
 }
strongif (!is_constant(fm.x) && is_constant(a.x) && !is_constant(alpha.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 395
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    _stencil_val(__FILE__,__LINE__,gf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)*(_stencil_val(__FILE__,__LINE__,p,0,0,0) - _stencil_val(__FILE__,__LINE__,p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    _stencil_val(__FILE__,__LINE__,gf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)*(_stencil_val(__FILE__,__LINE__,p,0,0,0) - _stencil_val(__FILE__,__LINE__,p,0,-1,0))/Delta; }  }}  end_foreach_face_stencil()
#line 396
 }
strongif (is_constant(fm.x) && is_constant(a.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 395
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    _stencil_val(__FILE__,__LINE__,gf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)*(_stencil_val(__FILE__,__LINE__,p,0,0,0) - _stencil_val(__FILE__,__LINE__,p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    _stencil_val(__FILE__,__LINE__,gf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)*(_stencil_val(__FILE__,__LINE__,p,0,0,0) - _stencil_val(__FILE__,__LINE__,p,0,-1,0))/Delta; }  }}  end_foreach_face_stencil()
#line 396
 }
strongif (!is_constant(fm.x) && !is_constant(a.x) && is_constant(alpha.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_a_x
#define val_a_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 395
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    _stencil_val(__FILE__,__LINE__,gf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)*(_stencil_val(__FILE__,__LINE__,p,0,0,0) - _stencil_val(__FILE__,__LINE__,p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    _stencil_val(__FILE__,__LINE__,gf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)*(_stencil_val(__FILE__,__LINE__,p,0,0,0) - _stencil_val(__FILE__,__LINE__,p,0,-1,0))/Delta; }  }}  end_foreach_face_stencil()
#line 396
 }
strongif (is_constant(fm.x) && !is_constant(a.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_a_x
#define val_a_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 395
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    _stencil_val(__FILE__,__LINE__,gf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)*(_stencil_val(__FILE__,__LINE__,p,0,0,0) - _stencil_val(__FILE__,__LINE__,p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    _stencil_val(__FILE__,__LINE__,gf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)*(_stencil_val(__FILE__,__LINE__,p,0,0,0) - _stencil_val(__FILE__,__LINE__,p,0,-1,0))/Delta; }  }}  end_foreach_face_stencil()
#line 396
 }
strongif (!is_constant(fm.x) && is_constant(a.x) && is_constant(alpha.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 395
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    _stencil_val(__FILE__,__LINE__,gf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)*(_stencil_val(__FILE__,__LINE__,p,0,0,0) - _stencil_val(__FILE__,__LINE__,p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    _stencil_val(__FILE__,__LINE__,gf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)*(_stencil_val(__FILE__,__LINE__,p,0,0,0) - _stencil_val(__FILE__,__LINE__,p,0,-1,0))/Delta; }  }}  end_foreach_face_stencil()
#line 396
 }
strongif (is_constant(fm.x) && is_constant(a.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 395
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    _stencil_val(__FILE__,__LINE__,gf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)*(_stencil_val(__FILE__,__LINE__,p,0,0,0) - _stencil_val(__FILE__,__LINE__,p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    _stencil_val(__FILE__,__LINE__,gf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)*(_stencil_val(__FILE__,__LINE__,p,0,0,0) - _stencil_val(__FILE__,__LINE__,p,0,-1,0))/Delta; }  }}  end_foreach_face_stencil()
#line 396
 }  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 396

strongif (!is_constant(fm.x) && !is_constant(a.x) && !is_constant(alpha.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 395
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 396
 end_foreach_face(); }
strongif (is_constant(fm.x) && !is_constant(a.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 395
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 396
 end_foreach_face(); }
strongif (!is_constant(fm.x) && is_constant(a.x) && !is_constant(alpha.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 395
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 396
 end_foreach_face(); }
strongif (is_constant(fm.x) && is_constant(a.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 395
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 396
 end_foreach_face(); }
strongif (!is_constant(fm.x) && !is_constant(a.x) && is_constant(alpha.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 395
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 396
 end_foreach_face(); }
strongif (is_constant(fm.x) && !is_constant(a.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 395
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 396
 end_foreach_face(); }
strongif (!is_constant(fm.x) && is_constant(a.x) && is_constant(alpha.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 395
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 396
 end_foreach_face(); }
strongif (is_constant(fm.x) && is_constant(a.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 395
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 395
{

#line 395 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 396
 end_foreach_face(); } }





  trash (((vector []){{g.x,g.y},{{-1},{-1}}}));
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", .line = 403,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 403
foreach_stencil(){

#line 403 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    {
#line 404

      _stencil_val(__FILE__,__LINE__,g.x,0,0,0) = (_stencil_val(__FILE__,__LINE__,gf.x,0,0,0) + _stencil_val(__FILE__,__LINE__,gf.x,1,0,0))/(val_fm_x(fm.x,0,0,0) + val_fm_x(fm.x,1,0,0) + 0.);
#line 404

      _stencil_val(__FILE__,__LINE__,g.y,0,0,0) = (_stencil_val(__FILE__,__LINE__,gf.y,0,0,0) + _stencil_val(__FILE__,__LINE__,gf.y,0,1,0))/(val_fm_y(fm.y,0,0,0) + val_fm_y(fm.y,0,1,0) + 0.);}; } end_foreach_stencil(); }
strongif (is_constant(fm.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 403
foreach_stencil(){

#line 403 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    {
#line 404

      _stencil_val(__FILE__,__LINE__,g.x,0,0,0) = (_stencil_val(__FILE__,__LINE__,gf.x,0,0,0) + _stencil_val(__FILE__,__LINE__,gf.x,1,0,0))/(val_fm_x(fm.x,0,0,0) + val_fm_x(fm.x,1,0,0) + 0.);
#line 404

      _stencil_val(__FILE__,__LINE__,g.y,0,0,0) = (_stencil_val(__FILE__,__LINE__,gf.y,0,0,0) + _stencil_val(__FILE__,__LINE__,gf.y,0,1,0))/(val_fm_y(fm.y,0,0,0) + val_fm_y(fm.y,0,1,0) + 0.);}; } end_foreach_stencil(); }  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 405

strongif (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#line 403
foreach(){

#line 403 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    {
#line 404

      val(g.x,0,0,0) = (val(gf.x,0,0,0) + val(gf.x,1,0,0))/(val_fm_x(fm.x,0,0,0) + val_fm_x(fm.x,1,0,0) + 0.);
#line 404

      val(g.y,0,0,0) = (val(gf.y,0,0,0) + val(gf.y,0,1,0))/(val_fm_y(fm.y,0,0,0) + val_fm_y(fm.y,0,1,0) + 0.);}; } end_foreach(); }
strongif (is_constant(fm.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 403
foreach(){

#line 403 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

    {
#line 404

      val(g.x,0,0,0) = (val(gf.x,0,0,0) + val(gf.x,1,0,0))/(val_fm_x(fm.x,0,0,0) + val_fm_x(fm.x,1,0,0) + 0.);
#line 404

      val(g.y,0,0,0) = (val(gf.y,0,0,0) + val(gf.y,0,1,0))/(val_fm_y(fm.y,0,0,0) + val_fm_y(fm.y,0,1,0) + 0.);}; } end_foreach(); } }
 delete (((scalar []){gf.x,gf.y,{-1}})); }






static int projection_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int projection (const int i, const double t, Event * _ev) { trace ("projection", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 413); 
{
  mgp = project ((struct Project){uf, p, alpha, dt, mgp.nrelax});
  centered_gradient (p, g);




  correction (dt);
 end_trace("projection", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 422); } return 0; } 





static int end_timestep_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int end_timestep (const int i, const double t, Event * _ev) { trace ("end_timestep", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 428); ; end_trace("end_timestep", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 428);  return 0; } 
#line 438 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
static int adapt_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int adapt (const int i, const double t, Event * _ev) { trace ("adapt", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 438);  {






  event ("properties");
 end_trace("adapt", "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 446); } return 0; } 
#line 19 "60d_plate_advancing_simulation.c"
#line 1 "vtk.h"
#line 1 "/home/fpl/softwares/basilisk/src/vtk.h"
void output_vtk (scalar * list, int n, FILE * fp, bool linear)
{
  fputs ("# vtk DataFile Version 2.0\n"
  "Basilisk\n"
  "ASCII\n"
  "DATASET STRUCTURED_GRID\n", fp);
  fprintf (fp, "DIMENSIONS %d %d 1\n", n, n);
  fprintf (fp, "POINTS %d double\n", n*n);

  if (linear)
    boundary_internal ((scalar *)(list), "/home/fpl/softwares/basilisk/src/vtk.h", 11);

  double fn = n;
  double Delta = L0/fn;
  for (int i = 0; i < n; i++) {
    double x = Delta*i + X0 + Delta/2.;
    for (int j = 0; j < n; j++) {
      double y = Delta*j + Y0 + Delta/2.;
      fprintf (fp, "%g %g 0\n", x, y);
    }
  }
  fprintf (fp, "POINT_DATA %d\n", n*n);
  strongif (list) for (scalar s = *list, *_i108 = list; ((scalar *)&s)->i >= 0; s = *++_i108) {
    fprintf (fp, "SCALARS %s double\n", _attribute[s.i].name);
    fputs ("LOOKUP_TABLE default\n", fp);
    double fn = n;
    double Delta = L0/fn;
    for (int i = 0; i < n; i++) {
      double x = Delta*i + X0 + Delta/2.;
      for (int j = 0; j < n; j++) {
 double y = Delta*j + Y0 + Delta/2., v;
 if (linear)
   v = interpolate ((struct _interpolate){s, x, y});
 else {
   Point point = locate ((struct _locate){x, y});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 35 "/home/fpl/softwares/basilisk/src/vtk.h"

   v = point.level >= 0 ? val(s,0,0,0) : nodata;
 }
 fprintf (fp, "%g\n", v);
      }
    }
  }
  fflush (fp);
}
#line 20 "60d_plate_advancing_simulation.c"
#line 1 "adapt_wavelet_leave_interface.h"
#line 1 "/home/fpl/softwares/basilisk/src/adapt_wavelet_leave_interface.h"
struct Adapt_leave_interface {
  scalar * slist;
  scalar * vol_frac;
  double * max;
  int maxlevel;
  int minlevel;
  int padding;
  scalar * list;
};


astats adapt_wavelet_leave_interface (struct Adapt_leave_interface p)
{
  if (p.list == NULL)
    p.list = all;
  if (is_constant(cm))
    restriction (p.slist);
  else {
    scalar * listr = list_concat (((scalar []){cm,{-1}}), p.slist);
    restriction (listr);
    pfree (listr,__func__,__FILE__,__LINE__);
  }


  astats st = {0, 0};
  scalar * listc = NULL;
  strongif (p.list) for (scalar s = *p.list, *_i109 = p.list; ((scalar *)&s)->i >= 0; s = *++_i109)
    if (!is_constant(s) && _attribute[s.i].restriction != no_restriction)
      listc = list_add (listc, s);


  if (p.minlevel < 1)
    p.minlevel = 1;
  ((Tree *)grid)->refined.n = 0;
  static const int refined = 1 << user, too_fine = 1 << (user + 1);
   { foreach_cell(){

#line 36 "/home/fpl/softwares/basilisk/src/adapt_wavelet_leave_interface.h"
 {
    if (is_active(cell)) {
      static const int too_coarse = 1 << (user + 2);
      if (is_leaf (cell)) {
 if (cell.flags & too_coarse) {
   cell.flags &= ~too_coarse;
   refine_cell (point, listc, refined, &((Tree *)grid)->refined);
   st.nf++;
 }
 continue;
      }
      else {
 if (cell.flags & refined) {

   cell.flags &= ~too_coarse;
   continue;
 }

 bool local = is_local(cell);
 if (!local)
    { foreach_child()
     if (is_local(cell))
       local = true, foreach_child_break(); end_foreach_child(); }
 if (local) {
   int i = 0;
   static const int just_fine = 1 << (user + 3);
   strongif (p.slist) for (scalar s = *p.slist, *_i110 = p.slist; ((scalar *)&s)->i >= 0; s = *++_i110) {
     double max = p.max[i++], sc[1 << 2];
     int c = 0;
      { foreach_child()
       sc[c++] = val(s,0,0,0); end_foreach_child(); }
     _attribute[s.i].prolongation (point, s);
     c = 0;
      { foreach_child() {
       double e = fabs(sc[c] - val(s,0,0,0));
       if (e > max && level < p.maxlevel) {
  cell.flags &= ~too_fine;
  cell.flags |= too_coarse;
       }
       else if ((e <= max/1.5 || level > p.maxlevel) &&
         !(cell.flags & (too_coarse|just_fine))) {
  if (level >= p.minlevel)
    cell.flags |= too_fine;
       }
       else if (!(cell.flags & too_coarse)) {
  cell.flags &= ~too_fine;
  cell.flags |= just_fine;
       }

       strongif (p.vol_frac) for (scalar vf = *p.vol_frac, *_i111 = p.vol_frac; ((scalar *)&vf)->i >= 0; vf = *++_i111) {
                if (val(vf,0,0,0) > 0.0001 && val(vf,0,0,0) < 0.9999 && level < p.maxlevel) {
                  cell.flags |= too_coarse;
                  cell.flags &= ~too_fine;
               cell.flags &= ~just_fine;
                    if (p.padding > 0){
                        { foreach_neighbor(p.padding){
                          cell.flags |= too_coarse;
                          cell.flags &= ~too_fine;
                          cell.flags &= ~just_fine;
                        } end_foreach_neighbor(); }
                    }
                }
              }


       val(s,0,0,0) = sc[c++];
     } end_foreach_child(); }
   }
    { foreach_child() {
     cell.flags &= ~just_fine;
     if (!is_leaf(cell)) {
       cell.flags &= ~too_coarse;
       if (level >= p.maxlevel)
  cell.flags |= too_fine;
     }
     else if (!is_active(cell))
       cell.flags &= ~too_coarse;
   } end_foreach_child(); }
 }
      }
    }
    else
      continue;
  } } end_foreach_cell(); }
  mpi_boundary_refine (listc);



  for (int l = depth(); l >= 0; l--) {
     { foreach_cell(){

#line 125 "/home/fpl/softwares/basilisk/src/adapt_wavelet_leave_interface.h"

      if (!(cell.pid < 0)) {
 if (level == l) {
   if (!is_leaf(cell)) {
     if (cell.flags & refined)

       cell.flags &= ~(refined|too_fine);
     else if (cell.flags & too_fine) {
       if (is_local(cell) && coarsen_cell (point, listc))
  st.nc++;
       cell.flags &= ~too_fine;
     }
   }
   if (cell.flags & too_fine)
     cell.flags &= ~too_fine;
   else if (level > 0 && (aparent(0,0,0).flags & too_fine))
     aparent(0,0,0).flags &= ~too_fine;
   continue;
 }
 else if (is_leaf(cell))
   continue;
      } } end_foreach_cell(); }
    mpi_boundary_coarsen (l, too_fine);
  }
  pfree (listc,__func__,__FILE__,__LINE__);

  mpi_all_reduce (st.nf, MPI_INT, MPI_SUM);
  mpi_all_reduce (st.nc, MPI_INT, MPI_SUM);
  if (st.nc || st.nf)
    mpi_boundary_update (p.list);

  return st;
}
#line 21 "60d_plate_advancing_simulation.c"
#line 1 "contact.h"
#line 1 "/home/fpl/softwares/basilisk/src/contact.h"
#line 11 "/home/fpl/softwares/basilisk/src/contact.h"
coord interface_normal (Point point, scalar c);
#if _call_interface_normal
static coord _interface_normal (Point point, scalar c);
#endif

#line 11


#line 1 "fractions.h"
#line 1 "/home/fpl/softwares/basilisk/src/fractions.h"
#line 12 "/home/fpl/softwares/basilisk/src/fractions.h"
#line 1 "geometry.h"
#line 1 "/home/fpl/softwares/basilisk/src/geometry.h"
#line 28 "/home/fpl/softwares/basilisk/src/geometry.h"
double line_alpha (double c, coord n)
{
  double alpha, n1, n2;

  n1 = fabs (n.x); n2 = fabs (n.y);
  if (n1 > n2)
    swap (double, n1, n2);

  c = clamp (c, 0., 1.);
  double v1 = n1/2.;
  if (c <= v1/n2)
    alpha = sqrt (2.*c*n1*n2);
  else if (c <= 1. - v1/n2)
    alpha = c*n2 + v1;
  else
    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));

  if (n.x < 0.)
    alpha += n.x;
  if (n.y < 0.)
    alpha += n.y;

  return alpha - (n.x + n.y)/2.;
}
#line 133 "/home/fpl/softwares/basilisk/src/geometry.h"
double line_area (double nx, double ny, double alpha)
{
  double a, v, area;

  alpha += (nx + ny)/2.;
  if (nx < 0.) {
    alpha -= nx;
    nx = - nx;
  }
  if (ny < 0.) {
    alpha -= ny;
    ny = - ny;
  }

  if (alpha <= 0.)
    return 0.;

  if (alpha >= nx + ny)
    return 1.;

  if (nx < 1e-10)
    area = alpha/ny;
  else if (ny < 1e-10)
    area = alpha/nx;
  else {
    v = sq(alpha);

    a = alpha - nx;
    if (a > 0.)
      v -= a*a;

    a = alpha - ny;
    if (a > 0.)
      v -= a*a;

    area = v/(2.*nx*ny);
  }

  return clamp (area, 0., 1.);
}
#line 237 "/home/fpl/softwares/basilisk/src/geometry.h"
double rectangle_fraction (coord n, double alpha, coord a, coord b)
{
  coord n1;
  {
#line 240
 {
    alpha -= n.x*(b.x + a.x)/2.;
    n1.x = n.x*(b.x - a.x);
  }
#line 240
 {
    alpha -= n.y*(b.y + a.y)/2.;
    n1.y = n.y*(b.y - a.y);
  }}
  return line_area(n1.x, n1.y, alpha);
}
#line 262 "/home/fpl/softwares/basilisk/src/geometry.h"
int facets (coord n, double alpha, coord p[2])
{
  int i = 0;
  for (double s = -0.5; s <= 0.5; s += 1.)
    {
#line 266

      if (fabs (n.y) > 1e-4 && i < 2) {
 double a = (alpha - s*n.x)/n.y;
 if (a >= -0.5 && a <= 0.5) {
   p[i].x = s;
   p[i++].y = a;
 }
      }
#line 266

      if (fabs (n.x) > 1e-4 && i < 2) {
 double a = (alpha - s*n.y)/n.x;
 if (a >= -0.5 && a <= 0.5) {
   p[i].y = s;
   p[i++].x = a;
 }
      }}
  return i;
}
#line 352 "/home/fpl/softwares/basilisk/src/geometry.h"
double line_length_center (coord m, double alpha, coord * p)
{
  alpha += (m.x + m.y)/2.;

  coord n = m;
  {
#line 357

    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
#line 357

    if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }}

  p->x = p->y = p->z = 0.;

  if (alpha <= 0. || alpha >= n.x + n.y)
    return 0.;

  {
#line 368

    if (n.x < 1e-4) {
      p->x = 0.;
      p->y = (m.y < 0. ? 1. - alpha : alpha) - 0.5;
      return 1.;
    }
#line 368

    if (n.y < 1e-4) {
      p->y = 0.;
      p->x = (m.x < 0. ? 1. - alpha : alpha) - 0.5;
      return 1.;
    }}

  if (alpha >= n.x) {
    p->x += 1.;
    p->y += (alpha - n.x)/n.y;
  }
  else
    p->x += alpha/n.x;

  double ax = p->x, ay = p->y;
  if (alpha >= n.y) {
    p->y += 1.;
    ay -= 1.;
    p->x += (alpha - n.y)/n.x;
    ax -= (alpha - n.y)/n.x;
  }
  else {
    p->y += alpha/n.y;
    ay -= alpha/n.y;
  }

  {
#line 394
 {
    p->x /= 2.;
    p->x = clamp (p->x, 0., 1.);
    if (m.x < 0.)
      p->x = 1. - p->x;
    p->x -= 0.5;
  }
#line 394
 {
    p->y /= 2.;
    p->y = clamp (p->y, 0., 1.);
    if (m.y < 0.)
      p->y = 1. - p->y;
    p->y -= 0.5;
  }}

  return sqrt (ax*ax + ay*ay);
}
#line 482 "/home/fpl/softwares/basilisk/src/geometry.h"
void line_center (coord m, double alpha, double a, coord * p)
{
  alpha += (m.x + m.y)/2.;

  coord n = m;
  {
#line 487

    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
#line 487

    if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }}

  p->z = 0.;
  if (alpha <= 0.) {
    p->x = p->y = -0.5;
    return;
  }

  if (alpha >= n.x + n.y) {
    p->x = p->y = 0.;
    return;
  }

  {
#line 504

    if (n.x < 1e-4) {
      p->x = 0.;
      p->y = sign(m.y)*(a/2. - 0.5);
      return;
    }
#line 504

    if (n.y < 1e-4) {
      p->y = 0.;
      p->x = sign(m.x)*(a/2. - 0.5);
      return;
    }}

  p->x = p->y = cube(alpha);

  {
#line 513
 {
    double b = alpha - n.x;
    if (b > 0.) {
      p->x -= sq(b)*(alpha + 2.*n.x);
      p->y -= cube(b);
    }
  }
#line 513
 {
    double b = alpha - n.y;
    if (b > 0.) {
      p->y -= sq(b)*(alpha + 2.*n.y);
      p->x -= cube(b);
    }
  }}

  {
#line 521
 {
    p->x /= 6.*sq(n.x)*n.y*a;
    p->x = sign(m.x)*(p->x - 0.5);
  }
#line 521
 {
    p->y /= 6.*sq(n.y)*n.x*a;
    p->y = sign(m.y)*(p->y - 0.5);
  }}
}
#line 13 "/home/fpl/softwares/basilisk/src/fractions.h"






#line 1 "myc2d.h"
#line 1 "/home/fpl/softwares/basilisk/src/myc2d.h"





coord mycs (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 7 "/home/fpl/softwares/basilisk/src/myc2d.h"

  int ix;
  double c_t,c_b,c_r,c_l;
  double mx0,my0,mx1,my1,mm1,mm2;


  c_t = val(c,-1,1,0) + val(c,0,1,0) + val(c,1,1,0);
  c_b = val(c,-1,-1,0) + val(c,0,-1,0) + val(c,1,-1,0);
  c_r = val(c,1,-1,0) + val(c,1,0,0) + val(c,1,1,0);
  c_l = val(c,-1,-1,0) + val(c,-1,0,0) + val(c,-1,1,0);



  mx0 = 0.5*(c_l-c_r);
  my0 = 0.5*(c_b-c_t);


  if (fabs(mx0) <= fabs(my0)) {
    my0 = my0 > 0. ? 1. : -1.;
    ix = 1;
  }
  else {
    mx0 = mx0 > 0. ? 1. : -1.;
    ix = 0;
  }


  mm1 = val(c,-1,-1,0) + 2.0*val(c,-1,0,0) + val(c,-1,1,0);
  mm2 = val(c,1,-1,0) + 2.0*val(c,1,0,0) + val(c,1,1,0);
  mx1 = mm1 - mm2 + 1.e-30;
  mm1 = val(c,-1,-1,0) + 2.0*val(c,0,-1,0) + val(c,1,-1,0);
  mm2 = val(c,-1,1,0) + 2.0*val(c,0,1,0) + val(c,1,1,0);
  my1 = mm1 - mm2 + 1.e-30;


  if (ix) {
    mm1 = fabs(my1);
    mm1 = fabs(mx1)/mm1;
    if (mm1 > fabs(mx0)) {
      mx0 = mx1;
      my0 = my1;
    }
  }
  else {
    mm1 = fabs(mx1);
    mm1 = fabs(my1)/mm1;
    if (mm1 > fabs(my0)) {
      mx0 = mx1;
      my0 = my1;
    }
  }



  mm1 = fabs(mx0) + fabs(my0);
  coord n = {mx0/mm1, my0/mm1};

  return n;

#if _call_mycs
}
#define _IN_STENCIL 1

#line 6
static coord _mycs (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 7 "/home/fpl/softwares/basilisk/src/myc2d.h"

  int ix;
  double c_t,c_b,c_r,c_l;
  double mx0,my0,mx1,my1,mm1,mm2;


  c_t = _stencil_val(__FILE__,__LINE__,c,-1,1,0) + _stencil_val(__FILE__,__LINE__,c,0,1,0) + _stencil_val(__FILE__,__LINE__,c,1,1,0);
  c_b = _stencil_val(__FILE__,__LINE__,c,-1,-1,0) + _stencil_val(__FILE__,__LINE__,c,0,-1,0) + _stencil_val(__FILE__,__LINE__,c,1,-1,0);
  c_r = _stencil_val(__FILE__,__LINE__,c,1,-1,0) + _stencil_val(__FILE__,__LINE__,c,1,0,0) + _stencil_val(__FILE__,__LINE__,c,1,1,0);
  c_l = _stencil_val(__FILE__,__LINE__,c,-1,-1,0) + _stencil_val(__FILE__,__LINE__,c,-1,0,0) + _stencil_val(__FILE__,__LINE__,c,-1,1,0);



  mx0 = 0.5*(c_l-c_r);
  my0 = 0.5*(c_b-c_t);


  IF (fabs(mx0) <= fabs(my0)) {
    my0 = my0 > 0. ? 1. : -1.;
    ix = 1;
  }
   {
    mx0 = mx0 > 0. ? 1. : -1.;
    ix = 0;
  }


  mm1 = _stencil_val(__FILE__,__LINE__,c,-1,-1,0) + 2.0*_stencil_val(__FILE__,__LINE__,c,-1,0,0) + _stencil_val(__FILE__,__LINE__,c,-1,1,0);
  mm2 = _stencil_val(__FILE__,__LINE__,c,1,-1,0) + 2.0*_stencil_val(__FILE__,__LINE__,c,1,0,0) + _stencil_val(__FILE__,__LINE__,c,1,1,0);
  mx1 = mm1 - mm2 + 1.e-30;
  mm1 = _stencil_val(__FILE__,__LINE__,c,-1,-1,0) + 2.0*_stencil_val(__FILE__,__LINE__,c,0,-1,0) + _stencil_val(__FILE__,__LINE__,c,1,-1,0);
  mm2 = _stencil_val(__FILE__,__LINE__,c,-1,1,0) + 2.0*_stencil_val(__FILE__,__LINE__,c,0,1,0) + _stencil_val(__FILE__,__LINE__,c,1,1,0);
  my1 = mm1 - mm2 + 1.e-30;


  IF (ix) {
    mm1 = fabs(my1);
    mm1 = fabs(mx1)/mm1;
    IF (mm1 > fabs(mx0)) {
      mx0 = mx1;
      my0 = my1;
    }
  }
   {
    mm1 = fabs(mx1);
    mm1 = fabs(my1)/mm1;
    IF (mm1 > fabs(my0)) {
      mx0 = mx1;
      my0 = my1;
    }
  }



  mm1 = fabs(mx0) + fabs(my0);
  coord n = {mx0/mm1, my0/mm1};

  return n;

#undef _IN_STENCIL

#endif

#line 65
}
#line 20 "/home/fpl/softwares/basilisk/src/fractions.h"
#line 41 "/home/fpl/softwares/basilisk/src/fractions.h"
void fraction_refine (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 42 "/home/fpl/softwares/basilisk/src/fractions.h"






  double cc = val(c,0,0,0);
  if (cc <= 0. || cc >= 1.)
     { foreach_child()
      val(c,0,0,0) = cc; end_foreach_child(); }
  else {




    coord n = mycs (point, c);
    double alpha = line_alpha (cc, n);






     { foreach_child() {
      static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
      coord nc;
      {
#line 68

 nc.x = child.x*n.x;
#line 68

 nc.y = child.y*n.y;}
      val(c,0,0,0) = rectangle_fraction (nc, alpha, a, b);
    } end_foreach_child(); }
  }

#if _call_fraction_refine
}
#define _IN_STENCIL 1
#define mycs _mycs

#line 41
static void _fraction_refine (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 42 "/home/fpl/softwares/basilisk/src/fractions.h"






  double cc = _stencil_val(__FILE__,__LINE__,c,0,0,0);
  IF (cc <= 0. || cc >= 1.)
     { foreach_child()
      _stencil_val(__FILE__,__LINE__,c,0,0,0) = cc; end_foreach_child(); }
   {




    coord n = mycs (point, c);
    double alpha = line_alpha (cc, n);






     { foreach_child() {
      static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
      coord nc;
      {
#line 68

 nc.x = child.x*n.x;
#line 68

 nc.y = child.y*n.y;}
      _stencil_val(__FILE__,__LINE__,c,0,0,0) = rectangle_fraction (nc, alpha, a, b);
    } end_foreach_child(); }
  }

#undef mycs
#undef _IN_STENCIL

#endif

#line 73
}











static void alpha_refine (Point point, scalar alpha)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 86 "/home/fpl/softwares/basilisk/src/fractions.h"

  vector n = _attribute[alpha.i].n;
  double alphac = 2.*val(alpha,0,0,0);
  coord m;
  {
#line 90

    m.x = val(n.x,0,0,0);
#line 90

    m.y = val(n.y,0,0,0);}
   { foreach_child() {
    val(alpha,0,0,0) = alphac;
    {
#line 94

      val(alpha,0,0,0) -= child.x*m.x/2.;
#line 94

      val(alpha,0,0,0) -= child.y*m.y/2.;}
  } end_foreach_child(); }

#if _call_alpha_refine
}
#define _IN_STENCIL 1

#line 85
static void _alpha_refine (Point point, scalar alpha)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 86 "/home/fpl/softwares/basilisk/src/fractions.h"

  vector n = _attribute[alpha.i].n;
  double alphac = 2.*_stencil_val(__FILE__,__LINE__,alpha,0,0,0);
  coord m;
  {
#line 90

    m.x = _stencil_val(__FILE__,__LINE__,n.x,0,0,0);
#line 90

    m.y = _stencil_val(__FILE__,__LINE__,n.y,0,0,0);}
   { foreach_child() {
    _stencil_val(__FILE__,__LINE__,alpha,0,0,0) = alphac;
    {
#line 94

      _stencil_val(__FILE__,__LINE__,alpha,0,0,0) -= child.x*m.x/2.;
#line 94

      _stencil_val(__FILE__,__LINE__,alpha,0,0,0) -= child.y*m.y/2.;}
  } end_foreach_child(); }

#undef _IN_STENCIL

#endif

#line 97
}
#line 121 "/home/fpl/softwares/basilisk/src/fractions.h"
struct Fractions {
  scalar Phi;
  scalar c;
  vector s;
  double val;
};


void fractions (struct Fractions a)
{ trace ("fractions", "/home/fpl/softwares/basilisk/src/fractions.h", 130);
  scalar Phi = a.Phi;
  scalar c = a.c;
  vector s = (a.s).x.i ? (a.s) : new_face_vector("s");
  double val = a.val;
#line 145 "/home/fpl/softwares/basilisk/src/fractions.h"
  vector p;
  p.x = s.y; p.y = s.x;
#line 155 "/home/fpl/softwares/basilisk/src/fractions.h"
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _val = val;
{ double val = _val; NOT_UNUSED(val);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/fractions.h", .line = 155,
    .each = "foreach_face", .first = _first_call
  };
foreach_face_stencil() { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 155
{

#line 155 "/home/fpl/softwares/basilisk/src/fractions.h"
 {





    IF ((_stencil_val(__FILE__,__LINE__,Phi,0,0,0) - val)*(_stencil_val(__FILE__,__LINE__,Phi,1,0,0) - val) < 0.) {






      _stencil_val(__FILE__,__LINE__,p.x,0,0,0) = (_stencil_val(__FILE__,__LINE__,Phi,0,0,0) - val)/(_stencil_val(__FILE__,__LINE__,Phi,0,0,0) - _stencil_val(__FILE__,__LINE__,Phi,1,0,0));
      IF (_stencil_val(__FILE__,__LINE__,Phi,0,0,0) < val)
 _stencil_val(__FILE__,__LINE__,p.x,0,0,0) = 1. - _stencil_val(__FILE__,__LINE__,p.x,0,0,0);
    }
#line 180 "/home/fpl/softwares/basilisk/src/fractions.h"
    
      _stencil_val(__FILE__,__LINE__,p.x,0,0,0) = (_stencil_val(__FILE__,__LINE__,Phi,0,0,0) > val || _stencil_val(__FILE__,__LINE__,Phi,1,0,0) > val);
  } }  }}  { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 155
{

#line 155 "/home/fpl/softwares/basilisk/src/fractions.h"
 {





    IF ((_stencil_val(__FILE__,__LINE__,Phi,0,0,0) - val)*(_stencil_val(__FILE__,__LINE__,Phi,0,1,0) - val) < 0.) {






      _stencil_val(__FILE__,__LINE__,p.y,0,0,0) = (_stencil_val(__FILE__,__LINE__,Phi,0,0,0) - val)/(_stencil_val(__FILE__,__LINE__,Phi,0,0,0) - _stencil_val(__FILE__,__LINE__,Phi,0,1,0));
      IF (_stencil_val(__FILE__,__LINE__,Phi,0,0,0) < val)
 _stencil_val(__FILE__,__LINE__,p.y,0,0,0) = 1. - _stencil_val(__FILE__,__LINE__,p.y,0,0,0);
    }
#line 180 "/home/fpl/softwares/basilisk/src/fractions.h"
    
      _stencil_val(__FILE__,__LINE__,p.y,0,0,0) = (_stencil_val(__FILE__,__LINE__,Phi,0,0,0) > val || _stencil_val(__FILE__,__LINE__,Phi,0,1,0) > val);
  } }  }}  end_foreach_face_stencil()
#line 182
 if (_first_call) {
 if (val != _val)
   reduction_warning ("/home/fpl/softwares/basilisk/src/fractions.h", 155, "val");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 182
foreach_face_generic() { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 155
{

#line 155 "/home/fpl/softwares/basilisk/src/fractions.h"
 {





    if ((val(Phi,0,0,0) - val)*(val(Phi,1,0,0) - val) < 0.) {






      val(p.x,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,1,0,0));
      if (val(Phi,0,0,0) < val)
 val(p.x,0,0,0) = 1. - val(p.x,0,0,0);
    }
#line 180 "/home/fpl/softwares/basilisk/src/fractions.h"
    else
      val(p.x,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,1,0,0) > val);
  } }  }}  { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 155
{

#line 155 "/home/fpl/softwares/basilisk/src/fractions.h"
 {





    if ((val(Phi,0,0,0) - val)*(val(Phi,0,1,0) - val) < 0.) {






      val(p.y,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,0,1,0));
      if (val(Phi,0,0,0) < val)
 val(p.y,0,0,0) = 1. - val(p.y,0,0,0);
    }
#line 180 "/home/fpl/softwares/basilisk/src/fractions.h"
    else
      val(p.y,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,0,1,0) > val);
  } }  }}  end_foreach_face_generic()
#line 182
 end_foreach_face(); }
#line 205 "/home/fpl/softwares/basilisk/src/fractions.h"
  scalar s_z = c;
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _val = val;
{ double val = _val; NOT_UNUSED(val);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/fractions.h", .line = 206,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 206 "/home/fpl/softwares/basilisk/src/fractions.h"


  {
#line 240 "/home/fpl/softwares/basilisk/src/fractions.h"
    coord n;
    double nn = 0.;
    {
#line 242
 {
      n.x = _stencil_val(__FILE__,__LINE__,p.y,0,0,0) - _stencil_val(__FILE__,__LINE__,p.y,1,0,0);
      nn += fabs(n.x);
    }
#line 242
 {
      n.y = _stencil_val(__FILE__,__LINE__,p.x,0,0,0) - _stencil_val(__FILE__,__LINE__,p.x,0,1,0);
      nn += fabs(n.y);
    }}





    IF (nn == 0.)
      _stencil_val(__FILE__,__LINE__,s_z,0,0,0) = _stencil_val(__FILE__,__LINE__,p.x,0,0,0);
     {





      {
#line 259

 n.x /= nn;
#line 259

 n.y /= nn;}






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
#line 269

   IF (_stencil_val(__FILE__,__LINE__,p.x,0,i,0) > 0. && _stencil_val(__FILE__,__LINE__,p.x,0,i,0) < 1.) {
     double a = sign(_stencil_val(__FILE__,__LINE__,Phi,0,i,0) - val)*(_stencil_val(__FILE__,__LINE__,p.x,0,i,0) - 0.5);
     alpha += n.x*a + n.y*(i - 0.5);
     ni++;
   }
#line 269

   IF (_stencil_val(__FILE__,__LINE__,p.y,i,0,0) > 0. && _stencil_val(__FILE__,__LINE__,p.y,i,0,0) < 1.) {
     double a = sign(_stencil_val(__FILE__,__LINE__,Phi,i,0,0) - val)*(_stencil_val(__FILE__,__LINE__,p.y,i,0,0) - 0.5);
     alpha += n.y*a + n.x*(i - 0.5);
     ni++;
   }}
#line 283 "/home/fpl/softwares/basilisk/src/fractions.h"
      IF (ni == 0)
 _stencil_val(__FILE__,__LINE__,s_z,0,0,0) = max (_stencil_val(__FILE__,__LINE__,p.x,0,0,0), _stencil_val(__FILE__,__LINE__,p.y,0,0,0));
      IF (ni != 4)
 _stencil_val(__FILE__,__LINE__,s_z,0,0,0) = line_area (n.x, n.y, alpha/ni);
       {



 _stencil_val(__FILE__,__LINE__,s_z,0,0,0) = 0.;

      }
    }
  } } end_foreach_stencil(); if (_first_call) {
 if (val != _val)
   reduction_warning ("/home/fpl/softwares/basilisk/src/fractions.h", 206, "val");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 295
foreach(){

#line 206 "/home/fpl/softwares/basilisk/src/fractions.h"


  {
#line 240 "/home/fpl/softwares/basilisk/src/fractions.h"
    coord n;
    double nn = 0.;
    {
#line 242
 {
      n.x = val(p.y,0,0,0) - val(p.y,1,0,0);
      nn += fabs(n.x);
    }
#line 242
 {
      n.y = val(p.x,0,0,0) - val(p.x,0,1,0);
      nn += fabs(n.y);
    }}





    if (nn == 0.)
      val(s_z,0,0,0) = val(p.x,0,0,0);
    else {





      {
#line 259

 n.x /= nn;
#line 259

 n.y /= nn;}






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
#line 269

   if (val(p.x,0,i,0) > 0. && val(p.x,0,i,0) < 1.) {
     double a = sign(val(Phi,0,i,0) - val)*(val(p.x,0,i,0) - 0.5);
     alpha += n.x*a + n.y*(i - 0.5);
     ni++;
   }
#line 269

   if (val(p.y,i,0,0) > 0. && val(p.y,i,0,0) < 1.) {
     double a = sign(val(Phi,i,0,0) - val)*(val(p.y,i,0,0) - 0.5);
     alpha += n.y*a + n.x*(i - 0.5);
     ni++;
   }}
#line 283 "/home/fpl/softwares/basilisk/src/fractions.h"
      if (ni == 0)
 val(s_z,0,0,0) = max (val(p.x,0,0,0), val(p.y,0,0,0));
      else if (ni != 4)
 val(s_z,0,0,0) = line_area (n.x, n.y, alpha/ni);
      else {



 val(s_z,0,0,0) = 0.;

      }
    }
  } } end_foreach(); }
#line 347 "/home/fpl/softwares/basilisk/src/fractions.h"
 { strongif (!(a.s).x.i) delete (((scalar []){s.x,s.y,{-1}})); }  end_trace("fractions", "/home/fpl/softwares/basilisk/src/fractions.h", 347); }
#line 391 "/home/fpl/softwares/basilisk/src/fractions.h"
coord youngs_normal (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 392 "/home/fpl/softwares/basilisk/src/fractions.h"

  coord n;
  double nn = 0.;
  if (!(2 == 2)) qassert ("/home/fpl/softwares/basilisk/src/fractions.h", 395, "dimension == 2");
  {
#line 396
 {
    n.x = (val(c,-1,1,0) + 2.*val(c,-1,0,0) + val(c,-1,-1,0) -
    val(c,+1,1,0) - 2.*val(c,+1,0,0) - val(c,+1,-1,0));
    nn += fabs(n.x);
  }
#line 396
 {
    n.y = (val(c,1,-1,0) + 2.*val(c,0,-1,0) + val(c,-1,-1,0) -
    val(c,1,+1,0) - 2.*val(c,0,+1,0) - val(c,-1,+1,0));
    nn += fabs(n.y);
  }}

  if (nn > 0.)
    {
#line 403

      n.x /= nn;
#line 403

      n.y /= nn;}
  else
    n.x = 1.;
  return n;

#if _call_youngs_normal
}
#define _IN_STENCIL 1

#line 391
static coord _youngs_normal (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 392 "/home/fpl/softwares/basilisk/src/fractions.h"

  coord n;
  double nn = 0.;
  IF (!(2 == 2)) _stencil_qassert (__FILE__,__LINE__,"/home/fpl/softwares/basilisk/src/fractions.h", 395, "dimension == 2");
  {
#line 396
 {
    n.x = (_stencil_val(__FILE__,__LINE__,c,-1,1,0) + 2.*_stencil_val(__FILE__,__LINE__,c,-1,0,0) + _stencil_val(__FILE__,__LINE__,c,-1,-1,0) -
    _stencil_val(__FILE__,__LINE__,c,+1,1,0) - 2.*_stencil_val(__FILE__,__LINE__,c,+1,0,0) - _stencil_val(__FILE__,__LINE__,c,+1,-1,0));
    nn += fabs(n.x);
  }
#line 396
 {
    n.y = (_stencil_val(__FILE__,__LINE__,c,1,-1,0) + 2.*_stencil_val(__FILE__,__LINE__,c,0,-1,0) + _stencil_val(__FILE__,__LINE__,c,-1,-1,0) -
    _stencil_val(__FILE__,__LINE__,c,1,+1,0) - 2.*_stencil_val(__FILE__,__LINE__,c,0,+1,0) - _stencil_val(__FILE__,__LINE__,c,-1,+1,0));
    nn += fabs(n.y);
  }}

  IF (nn > 0.)
    {
#line 403

      n.x /= nn;
#line 403

      n.y /= nn;}
  
    n.x = 1.;
  return n;

#undef _IN_STENCIL

#endif

#line 408
}





coord facet_normal (Point point, scalar c, vector s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 415 "/home/fpl/softwares/basilisk/src/fractions.h"

  if (s.x.i >= 0) {
    coord n;
    double nn = 0.;
    {
#line 419
 {
      n.x = val(s.x,0,0,0) - val(s.x,1,0,0);
      nn += fabs(n.x);
    }
#line 419
 {
      n.y = val(s.y,0,0,0) - val(s.y,0,1,0);
      nn += fabs(n.y);
    }}
    if (nn > 0.)
      {
#line 424

 n.x /= nn;
#line 424

 n.y /= nn;}
    else
      {
#line 427

 n.x = 1./2;
#line 427

 n.y = 1./2;}
    return n;
  }
  return interface_normal (point, c);

#if _call_facet_normal
}
#define _IN_STENCIL 1
#define interface_normal _interface_normal

#line 414
static coord _facet_normal (Point point, scalar c, vector s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 415 "/home/fpl/softwares/basilisk/src/fractions.h"

  IF (s.x.i >= 0) {
    coord n;
    double nn = 0.;
    {
#line 419
 {
      n.x = _stencil_val(__FILE__,__LINE__,s.x,0,0,0) - _stencil_val(__FILE__,__LINE__,s.x,1,0,0);
      nn += fabs(n.x);
    }
#line 419
 {
      n.y = _stencil_val(__FILE__,__LINE__,s.y,0,0,0) - _stencil_val(__FILE__,__LINE__,s.y,0,1,0);
      nn += fabs(n.y);
    }}
    IF (nn > 0.)
      {
#line 424

 n.x /= nn;
#line 424

 n.y /= nn;}
    
      {
#line 427

 n.x = 1./2;
#line 427

 n.y = 1./2;}
    return n;
  }
  return interface_normal (point, c);

#undef interface_normal
#undef _IN_STENCIL

#endif

#line 432
}
#line 441 "/home/fpl/softwares/basilisk/src/fractions.h"

void reconstruction (const scalar c, vector n, scalar alpha)
{ trace ("reconstruction", "/home/fpl/softwares/basilisk/src/fractions.h", 443);
   { 
#define interface_normal _interface_normal
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/fractions.h", .line = 444,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 444 "/home/fpl/softwares/basilisk/src/fractions.h"
 {





    IF (_stencil_val(__FILE__,__LINE__,c,0,0,0) <= 0. || _stencil_val(__FILE__,__LINE__,c,0,0,0) >= 1.) {
      _stencil_val(__FILE__,__LINE__,alpha,0,0,0) = 0.;
      {
#line 452

 _stencil_val(__FILE__,__LINE__,n.x,0,0,0) = 0.;
#line 452

 _stencil_val(__FILE__,__LINE__,n.y,0,0,0) = 0.;}
    }
     {






      coord m = interface_normal (point, c);
      {
#line 463

 _stencil_val(__FILE__,__LINE__,n.x,0,0,0) = m.x;
#line 463

 _stencil_val(__FILE__,__LINE__,n.y,0,0,0) = m.y;}
      _stencil_val(__FILE__,__LINE__,alpha,0,0,0) = line_alpha (_stencil_val(__FILE__,__LINE__,c,0,0,0), m);
    }
  } } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#undef interface_normal
#line 467
foreach(){

#line 444 "/home/fpl/softwares/basilisk/src/fractions.h"
 {





    if (val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) {
      val(alpha,0,0,0) = 0.;
      {
#line 452

 val(n.x,0,0,0) = 0.;
#line 452

 val(n.y,0,0,0) = 0.;}
    }
    else {






      coord m = interface_normal (point, c);
      {
#line 463

 val(n.x,0,0,0) = m.x;
#line 463

 val(n.y,0,0,0) = m.y;}
      val(alpha,0,0,0) = line_alpha (val(c,0,0,0), m);
    }
  } } end_foreach(); }
#line 476 "/home/fpl/softwares/basilisk/src/fractions.h"
  {
#line 476

    _attribute[n.x.i].refine = _attribute[n.x.i].prolongation = refine_injection;
#line 476

    _attribute[n.y.i].refine = _attribute[n.y.i].prolongation = refine_injection;}




  _attribute[alpha.i].n = n;
  _attribute[alpha.i].refine = _attribute[alpha.i].prolongation = alpha_refine;

 end_trace("reconstruction", "/home/fpl/softwares/basilisk/src/fractions.h", 485); }
#line 505 "/home/fpl/softwares/basilisk/src/fractions.h"
struct OutputFacets {
  scalar c;
  FILE * fp;
  vector s;
};


void output_facets (struct OutputFacets p)
{ trace ("output_facets", "/home/fpl/softwares/basilisk/src/fractions.h", 513);
  scalar c = p.c;
  vector s = p.s;
  if (!p.fp) p.fp = fout;
  if (!s.x.i) s.x.i = -1;

   { 
#define facet_normal _facet_normal
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/fractions.h", .line = 519,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 519 "/home/fpl/softwares/basilisk/src/fractions.h"

    IF (_stencil_val(__FILE__,__LINE__,c,0,0,0) > 1e-6 && _stencil_val(__FILE__,__LINE__,c,0,0,0) < 1. - 1e-6) {
      coord n = facet_normal (point, c, s);
      double alpha = line_alpha (_stencil_val(__FILE__,__LINE__,c,0,0,0), n);

      coord segment[2];
      IF (facets (n, alpha, segment) == 2)
 _stencil_fprintf (__FILE__,__LINE__,p.fp, "%g %g\n%g %g\n\n",
   x + segment[0].x*Delta, y + segment[0].y*Delta,
   x + segment[1].x*Delta, y + segment[1].y*Delta);
#line 538 "/home/fpl/softwares/basilisk/src/fractions.h"
    } } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#undef facet_normal
#line 538
foreach(){

#line 519 "/home/fpl/softwares/basilisk/src/fractions.h"

    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n = facet_normal (point, c, s);
      double alpha = line_alpha (val(c,0,0,0), n);

      coord segment[2];
      if (facets (n, alpha, segment) == 2)
 fprintf (p.fp, "%g %g\n%g %g\n\n",
   x + segment[0].x*Delta, y + segment[0].y*Delta,
   x + segment[1].x*Delta, y + segment[1].y*Delta);
#line 538 "/home/fpl/softwares/basilisk/src/fractions.h"
    } } end_foreach(); }

  fflush (p.fp);
 end_trace("output_facets", "/home/fpl/softwares/basilisk/src/fractions.h", 541); }








double interface_area (scalar c)
{ trace ("interface_area", "/home/fpl/softwares/basilisk/src/fractions.h", 551);
  double area = 0.;
   { 
#define interface_normal _interface_normal
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _area = area;
{ double area = _area; NOT_UNUSED(area);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/fractions.h", .line = 553,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 553 "/home/fpl/softwares/basilisk/src/fractions.h"

    IF (_stencil_val(__FILE__,__LINE__,c,0,0,0) > 1e-6 && _stencil_val(__FILE__,__LINE__,c,0,0,0) < 1. - 1e-6) {
      coord n = interface_normal (point, c), p;
      double alpha = line_alpha (_stencil_val(__FILE__,__LINE__,c,0,0,0), n);
      area += pow(Delta, 2 - 1)*line_length_center(n,alpha,&p);
    } } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#undef interface_normal
#line 558

#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:area)) {

#line 553
foreach (){

#line 553 "/home/fpl/softwares/basilisk/src/fractions.h"

    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n = interface_normal (point, c), p;
      double alpha = line_alpha (val(c,0,0,0), n);
      area += pow(Delta, 2 - 1)*line_length_center(n,alpha,&p);
    } } end_foreach();mpi_all_reduce_array (&area, double, MPI_SUM, 1);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 558
 }
  { double _ret =  area; end_trace("interface_area", "/home/fpl/softwares/basilisk/src/fractions.h", 559);  return _ret; }
 end_trace("interface_area", "/home/fpl/softwares/basilisk/src/fractions.h", 560); }
#line 15 "/home/fpl/softwares/basilisk/src/contact.h"
#line 1 "curvature.h"
#line 1 "/home/fpl/softwares/basilisk/src/curvature.h"
#line 12 "/home/fpl/softwares/basilisk/src/curvature.h"
static void curvature_restriction (Point point, scalar kappa)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 13 "/home/fpl/softwares/basilisk/src/curvature.h"

  double k = 0., s = 0.;
   { foreach_child()
    if (val(kappa,0,0,0) != nodata)
      k += val(kappa,0,0,0), s++; end_foreach_child(); }
  val(kappa,0,0,0) = s ? k/s : nodata;

#if _call_curvature_restriction
}
#define _IN_STENCIL 1

#line 12
static void _curvature_restriction (Point point, scalar kappa)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 13 "/home/fpl/softwares/basilisk/src/curvature.h"

  double k = 0., s = 0.;
   { foreach_child()
    IF (_stencil_val(__FILE__,__LINE__,kappa,0,0,0) != nodata)
      k += _stencil_val(__FILE__,__LINE__,kappa,0,0,0), s++; end_foreach_child(); }
  _stencil_val(__FILE__,__LINE__,kappa,0,0,0) = s ? k/s : nodata;

#undef _IN_STENCIL

#endif

#line 19
}







static void curvature_prolongation (Point point, scalar kappa)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 28 "/home/fpl/softwares/basilisk/src/curvature.h"

   { foreach_child() {
    double sk = 0., s = 0.;
    for (int i = 0; i <= 1; i++)

      for (int j = 0; j <= 1; j++)




   if (coarse(kappa,child.x*i,child.y*j,child.z*k) != nodata)
     sk += coarse(kappa,child.x*i,child.y*j,child.z*k), s++;
    val(kappa,0,0,0) = s ? sk/s : nodata;
  } end_foreach_child(); }

#if _call_curvature_prolongation
}
#define _IN_STENCIL 1

#line 27
static void _curvature_prolongation (Point point, scalar kappa)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 28 "/home/fpl/softwares/basilisk/src/curvature.h"

   { foreach_child() {
    double sk = 0., s = 0.;
    for (int i = 0; i <= 1; i++)

      for (int j = 0; j <= 1; j++)




   IF (_stencil_coarse(__FILE__,__LINE__,kappa,child.x*i,child.y*j,child.z*k) != nodata)
     sk += _stencil_coarse(__FILE__,__LINE__,kappa,child.x*i,child.y*j,child.z*k), s++;
    _stencil_val(__FILE__,__LINE__,kappa,0,0,0) = s ? sk/s : nodata;
  } end_foreach_child(); }

#undef _IN_STENCIL

#endif

#line 42
}
#line 66 "/home/fpl/softwares/basilisk/src/curvature.h"
#line 1 "heights.h"
#line 1 "/home/fpl/softwares/basilisk/src/heights.h"
#line 29 "/home/fpl/softwares/basilisk/src/heights.h"
static inline double height (double H) {
  return H > 20./2. ? H - 20. : H < -20./2. ? H + 20. : H;
}

static inline int orientation (double H) {
  return fabs(H) > 20./2.;
}
#line 49 "/home/fpl/softwares/basilisk/src/heights.h"
static void half_column (Point point, scalar c, vector h, vector cs, int j)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 50 "/home/fpl/softwares/basilisk/src/heights.h"







  const int complete = -1;

  {
#line 59
 {







    double S = val(c,0,0,0), H = S, ci, a;







    typedef struct { int s; double h; } HState;
    HState state = {0, 0};
    if (j == 1) {




      if (val(h.x,0,0,0) == 300.)
 state.s = complete, state.h = nodata;




      else {
 int s = (val(h.x,0,0,0) + 20./2.)/100.;
 state.h = val(h.x,0,0,0) - 100.*s;
 state.s = s - 1;
      }





      if (state.s != complete)
 S = state.s, H = state.h;
    }
#line 109 "/home/fpl/softwares/basilisk/src/heights.h"
    for (int i = 1; i <= 4; i++) {
      ci = i <= 2 ? val(c,i*j,0,0) : val(cs.x,(i - 2)*j,0,0);
      H += ci;




      if (S > 0. && S < 1.) {
 S = ci;
 if (ci <= 0. || ci >= 1.) {







   H -= i*ci;
   break;
 }
      }
#line 138 "/home/fpl/softwares/basilisk/src/heights.h"
      else if (S >= 1. && ci <= 0.) {
 H = (H - 0.5)*j + (j == -1)*20.;
 S = complete;
 break;
      }
      else if (S <= 0. && ci >= 1.) {
 H = (i + 0.5 - H)*j + (j == 1)*20.;
 S = complete;
 break;
      }
#line 156 "/home/fpl/softwares/basilisk/src/heights.h"
      else if (S == ci && modf(H, &a))
 break;
    }





    if (j == -1) {







      if (S != complete && ((val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) ||
       (S > 0. && S < 1.)))
 val(h.x,0,0,0) = 300.;
      else if (S == complete)
 val(h.x,0,0,0) = H;
      else





 val(h.x,0,0,0) = H + 100.*(1. + (S >= 1.));
    }
    else {
#line 195 "/home/fpl/softwares/basilisk/src/heights.h"
      if (state.s != complete ||
   (S == complete && fabs(height(H)) < fabs(height(state.h))))
 state.s = S, state.h = H;





      if (state.s != complete)
 val(h.x,0,0,0) = nodata;
      else
 val(h.x,0,0,0) = (state.h > 1e10 ? nodata : state.h);
    }
  }
#line 59
 {







    double S = val(c,0,0,0), H = S, ci, a;







    typedef struct { int s; double h; } HState;
    HState state = {0, 0};
    if (j == 1) {




      if (val(h.y,0,0,0) == 300.)
 state.s = complete, state.h = nodata;




      else {
 int s = (val(h.y,0,0,0) + 20./2.)/100.;
 state.h = val(h.y,0,0,0) - 100.*s;
 state.s = s - 1;
      }





      if (state.s != complete)
 S = state.s, H = state.h;
    }
#line 109 "/home/fpl/softwares/basilisk/src/heights.h"
    for (int i = 1; i <= 4; i++) {
      ci = i <= 2 ? val(c,0,i*j,0) : val(cs.y,0,(i - 2)*j,0);
      H += ci;




      if (S > 0. && S < 1.) {
 S = ci;
 if (ci <= 0. || ci >= 1.) {







   H -= i*ci;
   break;
 }
      }
#line 138 "/home/fpl/softwares/basilisk/src/heights.h"
      else if (S >= 1. && ci <= 0.) {
 H = (H - 0.5)*j + (j == -1)*20.;
 S = complete;
 break;
      }
      else if (S <= 0. && ci >= 1.) {
 H = (i + 0.5 - H)*j + (j == 1)*20.;
 S = complete;
 break;
      }
#line 156 "/home/fpl/softwares/basilisk/src/heights.h"
      else if (S == ci && modf(H, &a))
 break;
    }





    if (j == -1) {







      if (S != complete && ((val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) ||
       (S > 0. && S < 1.)))
 val(h.y,0,0,0) = 300.;
      else if (S == complete)
 val(h.y,0,0,0) = H;
      else





 val(h.y,0,0,0) = H + 100.*(1. + (S >= 1.));
    }
    else {
#line 195 "/home/fpl/softwares/basilisk/src/heights.h"
      if (state.s != complete ||
   (S == complete && fabs(height(H)) < fabs(height(state.h))))
 state.s = S, state.h = H;





      if (state.s != complete)
 val(h.y,0,0,0) = nodata;
      else
 val(h.y,0,0,0) = (state.h > 1e10 ? nodata : state.h);
    }
  }}

#if _call_half_column
}
#define _IN_STENCIL 1

#line 49
static void _half_column (Point point, scalar c, vector h, vector cs, int j)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 50 "/home/fpl/softwares/basilisk/src/heights.h"







  const int complete = -1;

  {
#line 59
 {







    double S = _stencil_val(__FILE__,__LINE__,c,0,0,0), H = S, ci, a;







    typedef struct { int s; double h; } HState;
    HState state = {0, 0};
    IF (j == 1) {




      IF (_stencil_val(__FILE__,__LINE__,h.x,0,0,0) == 300.)
 state.s = complete, state.h = nodata;




       {
 int s = (_stencil_val(__FILE__,__LINE__,h.x,0,0,0) + 20./2.)/100.;
 state.h = _stencil_val(__FILE__,__LINE__,h.x,0,0,0) - 100.*s;
 state.s = s - 1;
      }





      IF (state.s != complete)
 S = state.s, H = state.h;
    }
#line 109 "/home/fpl/softwares/basilisk/src/heights.h"
    for (int i = 1; i <= 4; i++) {
      ci = i <= 2 ? _stencil_val(__FILE__,__LINE__,c,i*j,0,0) : _stencil_val(__FILE__,__LINE__,cs.x,(i - 2)*j,0,0);
      H += ci;




      IF (S > 0. && S < 1.) {
 S = ci;
 IF (ci <= 0. || ci >= 1.) {







   H -= i*ci;
   ;
 }
      }
#line 138 "/home/fpl/softwares/basilisk/src/heights.h"
      IF (S >= 1. && ci <= 0.) {
 H = (H - 0.5)*j + (j == -1)*20.;
 S = complete;
 ;
      }
      IF (S <= 0. && ci >= 1.) {
 H = (i + 0.5 - H)*j + (j == 1)*20.;
 S = complete;
 ;
      }
#line 156 "/home/fpl/softwares/basilisk/src/heights.h"
      IF (S == ci && modf(H, &a))
 ;
    }





    IF (j == -1) {







      IF (S != complete && ((_stencil_val(__FILE__,__LINE__,c,0,0,0) <= 0. || _stencil_val(__FILE__,__LINE__,c,0,0,0) >= 1.) ||
       (S > 0. && S < 1.)))
 _stencil_val(__FILE__,__LINE__,h.x,0,0,0) = 300.;
      IF (S == complete)
 _stencil_val(__FILE__,__LINE__,h.x,0,0,0) = H;
      





 _stencil_val(__FILE__,__LINE__,h.x,0,0,0) = H + 100.*(1. + (S >= 1.));
    }
     {
#line 195 "/home/fpl/softwares/basilisk/src/heights.h"
      IF (state.s != complete ||
   (S == complete && fabs(height(H)) < fabs(height(state.h))))
 state.s = S, state.h = H;





      IF (state.s != complete)
 _stencil_val(__FILE__,__LINE__,h.x,0,0,0) = nodata;
      
 _stencil_val(__FILE__,__LINE__,h.x,0,0,0) = (state.h > 1e10 ? nodata : state.h);
    }
  }
#line 59
 {







    double S = _stencil_val(__FILE__,__LINE__,c,0,0,0), H = S, ci, a;







    typedef struct { int s; double h; } HState;
    HState state = {0, 0};
    IF (j == 1) {




      IF (_stencil_val(__FILE__,__LINE__,h.y,0,0,0) == 300.)
 state.s = complete, state.h = nodata;




       {
 int s = (_stencil_val(__FILE__,__LINE__,h.y,0,0,0) + 20./2.)/100.;
 state.h = _stencil_val(__FILE__,__LINE__,h.y,0,0,0) - 100.*s;
 state.s = s - 1;
      }





      IF (state.s != complete)
 S = state.s, H = state.h;
    }
#line 109 "/home/fpl/softwares/basilisk/src/heights.h"
    for (int i = 1; i <= 4; i++) {
      ci = i <= 2 ? _stencil_val(__FILE__,__LINE__,c,0,i*j,0) : _stencil_val(__FILE__,__LINE__,cs.y,0,(i - 2)*j,0);
      H += ci;




      IF (S > 0. && S < 1.) {
 S = ci;
 IF (ci <= 0. || ci >= 1.) {







   H -= i*ci;
   ;
 }
      }
#line 138 "/home/fpl/softwares/basilisk/src/heights.h"
      IF (S >= 1. && ci <= 0.) {
 H = (H - 0.5)*j + (j == -1)*20.;
 S = complete;
 ;
      }
      IF (S <= 0. && ci >= 1.) {
 H = (i + 0.5 - H)*j + (j == 1)*20.;
 S = complete;
 ;
      }
#line 156 "/home/fpl/softwares/basilisk/src/heights.h"
      IF (S == ci && modf(H, &a))
 ;
    }





    IF (j == -1) {







      IF (S != complete && ((_stencil_val(__FILE__,__LINE__,c,0,0,0) <= 0. || _stencil_val(__FILE__,__LINE__,c,0,0,0) >= 1.) ||
       (S > 0. && S < 1.)))
 _stencil_val(__FILE__,__LINE__,h.y,0,0,0) = 300.;
      IF (S == complete)
 _stencil_val(__FILE__,__LINE__,h.y,0,0,0) = H;
      





 _stencil_val(__FILE__,__LINE__,h.y,0,0,0) = H + 100.*(1. + (S >= 1.));
    }
     {
#line 195 "/home/fpl/softwares/basilisk/src/heights.h"
      IF (state.s != complete ||
   (S == complete && fabs(height(H)) < fabs(height(state.h))))
 state.s = S, state.h = H;





      IF (state.s != complete)
 _stencil_val(__FILE__,__LINE__,h.y,0,0,0) = nodata;
      
 _stencil_val(__FILE__,__LINE__,h.y,0,0,0) = (state.h > 1e10 ? nodata : state.h);
    }
  }}

#undef _IN_STENCIL

#endif

#line 209
}
#line 222 "/home/fpl/softwares/basilisk/src/heights.h"
static void column_propagation (vector h)
{
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/heights.h", .line = 224,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 224 "/home/fpl/softwares/basilisk/src/heights.h"

    for (int i = -2; i <= 2; i++)
      {
#line 226

 IF (fabs(height(_stencil_val(__FILE__,__LINE__,h.x,i,0,0))) <= 3.5 &&
     fabs(height(_stencil_val(__FILE__,__LINE__,h.x,i,0,0)) + i) < fabs(height(_stencil_val(__FILE__,__LINE__,h.x,0,0,0))))
   _stencil_val(__FILE__,__LINE__,h.x,0,0,0) = _stencil_val(__FILE__,__LINE__,h.x,i,0,0) + i;
#line 226

 IF (fabs(height(_stencil_val(__FILE__,__LINE__,h.y,0,i,0))) <= 3.5 &&
     fabs(height(_stencil_val(__FILE__,__LINE__,h.y,0,i,0)) + i) < fabs(height(_stencil_val(__FILE__,__LINE__,h.y,0,0,0))))
   _stencil_val(__FILE__,__LINE__,h.y,0,0,0) = _stencil_val(__FILE__,__LINE__,h.y,0,i,0) + i;}; } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 229

#if _OPENMP
  #undef OMP
  #define OMP(x)
#endif
#line 224
foreach (){

#line 224 "/home/fpl/softwares/basilisk/src/heights.h"

    for (int i = -2; i <= 2; i++)
      {
#line 226

 if (fabs(height(val(h.x,i,0,0))) <= 3.5 &&
     fabs(height(val(h.x,i,0,0)) + i) < fabs(height(val(h.x,0,0,0))))
   val(h.x,0,0,0) = val(h.x,i,0,0) + i;
#line 226

 if (fabs(height(val(h.y,0,i,0))) <= 3.5 &&
     fabs(height(val(h.y,0,i,0)) + i) < fabs(height(val(h.y,0,0,0))))
   val(h.y,0,0,0) = val(h.y,0,i,0) + i;}; } end_foreach();
#if _OPENMP
  #undef OMP
  #define OMP(x) _Pragma(#x)
#endif
#line 229
 }
}
#line 288 "/home/fpl/softwares/basilisk/src/heights.h"

#line 288

static void refine_h_x (Point point, scalar h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 290 "/home/fpl/softwares/basilisk/src/heights.h"





  bool complete = true;
   { foreach_child() {
    for (int i = -2; i <= 2; i++)
      if (allocated(i,0,0) &&
   !(!is_leaf(neighbor(i,0,0)) && !neighbor(i,0,0).neighbors && neighbor(i,0,0).pid >= 0) && !(neighbor(i,0,0).pid < 0) &&
   fabs(height(val(h,i,0,0))) <= 3.5 &&
   fabs(height(val(h,i,0,0)) + i) < fabs(height(val(h,0,0,0))))
 val(h,0,0,0) = val(h,i,0,0) + i;
    if (val(h,0,0,0) == nodata)
      complete = false;
  } end_foreach_child(); }
  if (complete)
    return;
#line 316 "/home/fpl/softwares/basilisk/src/heights.h"
  int ori = orientation(val(h,0,0,0));

  for (int i = -1; i <= 1; i++)
    if (val(h,0,i,0) == nodata || orientation(val(h,0,i,0)) != ori)
      return;

  double h0 = (30.*height(val(h,0,0,0)) + height(val(h,0,1,0)) + height(val(h,0,-1,0)))/16.
    + 20.*ori;
  double dh = (height(val(h,0,1,0)) - height(val(h,0,-1,0)))/4.;
   { foreach_child()
    if (val(h,0,0,0) == nodata)
      val(h,0,0,0) = h0 + dh*child.y - child.x/2.; end_foreach_child(); }
#line 351 "/home/fpl/softwares/basilisk/src/heights.h"

#if _call_refine_h_x
}
#define _IN_STENCIL 1

#line 289
static void _refine_h_x (Point point, scalar h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 290 "/home/fpl/softwares/basilisk/src/heights.h"





  bool complete = true;
   { foreach_child() {
    for (int i = -2; i <= 2; i++)
      IF (allocated(i,0,0) &&
   !(!is_leaf(neighbor(i,0,0)) && !neighbor(i,0,0).neighbors && neighbor(i,0,0).pid >= 0) && !(neighbor(i,0,0).pid < 0) &&
   fabs(height(_stencil_val(__FILE__,__LINE__,h,i,0,0))) <= 3.5 &&
   fabs(height(_stencil_val(__FILE__,__LINE__,h,i,0,0)) + i) < fabs(height(_stencil_val(__FILE__,__LINE__,h,0,0,0))))
 _stencil_val(__FILE__,__LINE__,h,0,0,0) = _stencil_val(__FILE__,__LINE__,h,i,0,0) + i;
    IF (_stencil_val(__FILE__,__LINE__,h,0,0,0) == nodata)
      complete = false;
  } end_foreach_child(); }
  IF (complete)
    return;
#line 316 "/home/fpl/softwares/basilisk/src/heights.h"
  int ori = orientation(_stencil_val(__FILE__,__LINE__,h,0,0,0));

  for (int i = -1; i <= 1; i++)
    IF (_stencil_val(__FILE__,__LINE__,h,0,i,0) == nodata || orientation(_stencil_val(__FILE__,__LINE__,h,0,i,0)) != ori)
      return;

  double h0 = (30.*height(_stencil_val(__FILE__,__LINE__,h,0,0,0)) + height(_stencil_val(__FILE__,__LINE__,h,0,1,0)) + height(_stencil_val(__FILE__,__LINE__,h,0,-1,0)))/16.
    + 20.*ori;
  double dh = (height(_stencil_val(__FILE__,__LINE__,h,0,1,0)) - height(_stencil_val(__FILE__,__LINE__,h,0,-1,0)))/4.;
   { foreach_child()
    IF (_stencil_val(__FILE__,__LINE__,h,0,0,0) == nodata)
      _stencil_val(__FILE__,__LINE__,h,0,0,0) = h0 + dh*child.y - child.x/2.; end_foreach_child(); }
#line 351 "/home/fpl/softwares/basilisk/src/heights.h"

#undef _IN_STENCIL

#endif

#line 351
}
#line 288

static void refine_h_y (Point point, scalar h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 290 "/home/fpl/softwares/basilisk/src/heights.h"





  bool complete = true;
   { foreach_child() {
    for (int i = -2; i <= 2; i++)
      if (allocated(0,i,0) &&
   !(!is_leaf(neighbor(0,i,0)) && !neighbor(0,i,0).neighbors && neighbor(0,i,0).pid >= 0) && !(neighbor(0,i,0).pid < 0) &&
   fabs(height(val(h,0,i,0))) <= 3.5 &&
   fabs(height(val(h,0,i,0)) + i) < fabs(height(val(h,0,0,0))))
 val(h,0,0,0) = val(h,0,i,0) + i;
    if (val(h,0,0,0) == nodata)
      complete = false;
  } end_foreach_child(); }
  if (complete)
    return;
#line 316 "/home/fpl/softwares/basilisk/src/heights.h"
  int ori = orientation(val(h,0,0,0));

  for (int i = -1; i <= 1; i++)
    if (val(h,i,0,0) == nodata || orientation(val(h,i,0,0)) != ori)
      return;

  double h0 = (30.*height(val(h,0,0,0)) + height(val(h,1,0,0)) + height(val(h,-1,0,0)))/16.
    + 20.*ori;
  double dh = (height(val(h,1,0,0)) - height(val(h,-1,0,0)))/4.;
   { foreach_child()
    if (val(h,0,0,0) == nodata)
      val(h,0,0,0) = h0 + dh*child.x - child.y/2.; end_foreach_child(); }
#line 351 "/home/fpl/softwares/basilisk/src/heights.h"

#if _call_refine_h_y
}
#define _IN_STENCIL 1

#line 289
static void _refine_h_y (Point point, scalar h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 290 "/home/fpl/softwares/basilisk/src/heights.h"





  bool complete = true;
   { foreach_child() {
    for (int i = -2; i <= 2; i++)
      IF (allocated(0,i,0) &&
   !(!is_leaf(neighbor(0,i,0)) && !neighbor(0,i,0).neighbors && neighbor(0,i,0).pid >= 0) && !(neighbor(0,i,0).pid < 0) &&
   fabs(height(_stencil_val(__FILE__,__LINE__,h,i,0,0))) <= 3.5 &&
   fabs(height(_stencil_val(__FILE__,__LINE__,h,i,0,0)) + i) < fabs(height(_stencil_val(__FILE__,__LINE__,h,0,0,0))))
 _stencil_val(__FILE__,__LINE__,h,0,0,0) = _stencil_val(__FILE__,__LINE__,h,i,0,0) + i;
    IF (_stencil_val(__FILE__,__LINE__,h,0,0,0) == nodata)
      complete = false;
  } end_foreach_child(); }
  IF (complete)
    return;
#line 316 "/home/fpl/softwares/basilisk/src/heights.h"
  int ori = orientation(_stencil_val(__FILE__,__LINE__,h,0,0,0));

  for (int i = -1; i <= 1; i++)
    IF (_stencil_val(__FILE__,__LINE__,h,0,i,0) == nodata || orientation(_stencil_val(__FILE__,__LINE__,h,0,i,0)) != ori)
      return;

  double h0 = (30.*height(_stencil_val(__FILE__,__LINE__,h,0,0,0)) + height(_stencil_val(__FILE__,__LINE__,h,0,1,0)) + height(_stencil_val(__FILE__,__LINE__,h,0,-1,0)))/16.
    + 20.*ori;
  double dh = (height(_stencil_val(__FILE__,__LINE__,h,0,1,0)) - height(_stencil_val(__FILE__,__LINE__,h,0,-1,0)))/4.;
   { foreach_child()
    IF (_stencil_val(__FILE__,__LINE__,h,0,0,0) == nodata)
      _stencil_val(__FILE__,__LINE__,h,0,0,0) = h0 + dh*child.x - child.y/2.; end_foreach_child(); }
#line 351 "/home/fpl/softwares/basilisk/src/heights.h"

#undef _IN_STENCIL

#endif

#line 351
}







void heights (scalar c, vector h)
{ trace ("heights", "/home/fpl/softwares/basilisk/src/heights.h", 360);
  vector s= new_vector("s");
  {
#line 362

    for (int i = 0; i < nboundary; i++)
      _attribute[s.x.i].boundary[i] = _attribute[c.i].boundary[i];
#line 362

    for (int i = 0; i < nboundary; i++)
      _attribute[s.y.i].boundary[i] = _attribute[c.i].boundary[i];}





  restriction (((scalar []){c,{-1}}));
  for (int j = -1; j <= 1; j += 2) {





     { foreach_level(0){

#line 377 "/home/fpl/softwares/basilisk/src/heights.h"

      {
#line 378

        val(h.x,0,0,0) = nodata;
#line 378

        val(h.y,0,0,0) = nodata;}; } end_foreach_level(); }

    for (int l = 1; l <= depth(); l++) {




       { foreach_level (l){

#line 386 "/home/fpl/softwares/basilisk/src/heights.h"

 {
#line 387

   val(s.x,0,0,0) = val(c,2*j,0,0);
#line 387

   val(s.y,0,0,0) = val(c,0,2*j,0);}; } end_foreach_level(); }
#line 398 "/home/fpl/softwares/basilisk/src/heights.h"
       { foreach_level (l - 1){

#line 398 "/home/fpl/softwares/basilisk/src/heights.h"

 {
#line 399
 {
   val(s.x,0,0,0) = val(c,j,0,0);
   val(s.x,j,0,0) = val(c,2*j,0,0);
        }
#line 399
 {
   val(s.y,0,0,0) = val(c,0,j,0);
   val(s.y,0,j,0) = val(c,0,2*j,0);
        }} } end_foreach_level(); }






       { foreach_halo (prolongation, l - 1){

#line 409 "/home/fpl/softwares/basilisk/src/heights.h"

 {
#line 410

   _attribute[c.i].prolongation (point, s.x);
#line 410

   _attribute[c.i].prolongation (point, s.y);}; } end_foreach_halo(); }
      { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, (scalar *)((vector []){{s.x,s.y},{{-1},{-1}}}), l); };





       { foreach_level (l){

#line 418 "/home/fpl/softwares/basilisk/src/heights.h"

        half_column (point, c, h, s, j); } end_foreach_level(); }
    }
  }






  {
#line 428
 {
    _attribute[h.x.i].prolongation = no_data;
    _attribute[h.x.i].restriction = no_restriction;
    _attribute[h.x.i].dirty = true;
  }
#line 428
 {
    _attribute[h.y.i].prolongation = no_data;
    _attribute[h.y.i].restriction = no_restriction;
    _attribute[h.y.i].dirty = true;
  }}




  column_propagation (h);






  {
#line 444

    _attribute[h.x.i].prolongation = refine_h_x;
#line 444

    _attribute[h.y.i].prolongation = refine_h_y;}
 delete (((scalar []){s.x,s.y,{-1}}));  end_trace("heights", "/home/fpl/softwares/basilisk/src/heights.h", 446); }










#line 67 "/home/fpl/softwares/basilisk/src/curvature.h"



#line 69

static double kappa_y (Point point, vector h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 71 "/home/fpl/softwares/basilisk/src/curvature.h"

  int ori = orientation(val(h.y,0,0,0));
  for (int i = -1; i <= 1; i++)
    if (val(h.y,i,0,0) == nodata || orientation(val(h.y,i,0,0)) != ori)
      return nodata;
  double hx = (val(h.y,1,0,0) - val(h.y,-1,0,0))/2.;
  double hxx = (val(h.y,1,0,0) + val(h.y,-1,0,0) - 2.*val(h.y,0,0,0))/Delta;
  return hxx/pow(1. + sq(hx), 3/2.);

#if _call_kappa_y
}
#define _IN_STENCIL 1

#line 70
static double _kappa_y (Point point, vector h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 71 "/home/fpl/softwares/basilisk/src/curvature.h"

  int ori = orientation(_stencil_val(__FILE__,__LINE__,h.y,0,0,0));
  for (int i = -1; i <= 1; i++)
    IF (_stencil_val(__FILE__,__LINE__,h.y,i,0,0) == nodata || orientation(_stencil_val(__FILE__,__LINE__,h.y,i,0,0)) != ori)
      return nodata;
  double hx = (_stencil_val(__FILE__,__LINE__,h.y,1,0,0) - _stencil_val(__FILE__,__LINE__,h.y,-1,0,0))/2.;
  double hxx = (_stencil_val(__FILE__,__LINE__,h.y,1,0,0) + _stencil_val(__FILE__,__LINE__,h.y,-1,0,0) - 2.*_stencil_val(__FILE__,__LINE__,h.y,0,0,0))/Delta;
  return hxx/pow(1. + sq(hx), 3/2.);

#undef _IN_STENCIL

#endif

#line 79
}
#line 69

static double kappa_x (Point point, vector h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 71 "/home/fpl/softwares/basilisk/src/curvature.h"

  int ori = orientation(val(h.x,0,0,0));
  for (int i = -1; i <= 1; i++)
    if (val(h.x,0,i,0) == nodata || orientation(val(h.x,0,i,0)) != ori)
      return nodata;
  double hx = (val(h.x,0,1,0) - val(h.x,0,-1,0))/2.;
  double hxx = (val(h.x,0,1,0) + val(h.x,0,-1,0) - 2.*val(h.x,0,0,0))/Delta;
  return hxx/pow(1. + sq(hx), 3/2.);

#if _call_kappa_x
}
#define _IN_STENCIL 1

#line 70
static double _kappa_x (Point point, vector h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 71 "/home/fpl/softwares/basilisk/src/curvature.h"

  int ori = orientation(_stencil_val(__FILE__,__LINE__,h.x,0,0,0));
  for (int i = -1; i <= 1; i++)
    IF (_stencil_val(__FILE__,__LINE__,h.x,i,0,0) == nodata || orientation(_stencil_val(__FILE__,__LINE__,h.x,i,0,0)) != ori)
      return nodata;
  double hx = (_stencil_val(__FILE__,__LINE__,h.x,1,0,0) - _stencil_val(__FILE__,__LINE__,h.x,-1,0,0))/2.;
  double hxx = (_stencil_val(__FILE__,__LINE__,h.x,1,0,0) + _stencil_val(__FILE__,__LINE__,h.x,-1,0,0) - 2.*_stencil_val(__FILE__,__LINE__,h.x,0,0,0))/Delta;
  return hxx/pow(1. + sq(hx), 3/2.);

#undef _IN_STENCIL

#endif

#line 79
}


#line 81

static coord normal_y (Point point, vector h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 83 "/home/fpl/softwares/basilisk/src/curvature.h"

  coord n = {nodata, nodata, nodata};
  if (val(h.y,0,0,0) == nodata)
    return n;
  int ori = orientation(val(h.y,0,0,0));
  if (val(h.y,-1,0,0) != nodata && orientation(val(h.y,-1,0,0)) == ori) {
    if (val(h.y,1,0,0) != nodata && orientation(val(h.y,1,0,0)) == ori)
      n.x = (val(h.y,-1,0,0) - val(h.y,1,0,0))/2.;
    else
      n.x = val(h.y,-1,0,0) - val(h.y,0,0,0);
  }
  else if (val(h.y,1,0,0) != nodata && orientation(val(h.y,1,0,0)) == ori)
    n.x = val(h.y,0,0,0) - val(h.y,1,0,0);
  else
    return n;
  double nn = (ori ? -1. : 1.)*sqrt(1. + sq(n.x));
  n.x /= nn;
  n.y = 1./nn;
  return n;

#if _call_normal_y
}
#define _IN_STENCIL 1

#line 82
static coord _normal_y (Point point, vector h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 83 "/home/fpl/softwares/basilisk/src/curvature.h"

  coord n = {nodata, nodata, nodata};
  IF (_stencil_val(__FILE__,__LINE__,h.y,0,0,0) == nodata)
    return n;
  int ori = orientation(_stencil_val(__FILE__,__LINE__,h.y,0,0,0));
  IF (_stencil_val(__FILE__,__LINE__,h.y,-1,0,0) != nodata && orientation(_stencil_val(__FILE__,__LINE__,h.y,-1,0,0)) == ori) {
    IF (_stencil_val(__FILE__,__LINE__,h.y,1,0,0) != nodata && orientation(_stencil_val(__FILE__,__LINE__,h.y,1,0,0)) == ori)
      n.x = (_stencil_val(__FILE__,__LINE__,h.y,-1,0,0) - _stencil_val(__FILE__,__LINE__,h.y,1,0,0))/2.;
    
      n.x = _stencil_val(__FILE__,__LINE__,h.y,-1,0,0) - _stencil_val(__FILE__,__LINE__,h.y,0,0,0);
  }
  IF (_stencil_val(__FILE__,__LINE__,h.y,1,0,0) != nodata && orientation(_stencil_val(__FILE__,__LINE__,h.y,1,0,0)) == ori)
    n.x = _stencil_val(__FILE__,__LINE__,h.y,0,0,0) - _stencil_val(__FILE__,__LINE__,h.y,1,0,0);
  
    return n;
  double nn = (ori ? -1. : 1.)*sqrt(1. + sq(n.x));
  n.x /= nn;
  n.y = 1./nn;
  return n;

#undef _IN_STENCIL

#endif

#line 102
}
#line 81

static coord normal_x (Point point, vector h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 83 "/home/fpl/softwares/basilisk/src/curvature.h"

  coord n = {nodata, nodata, nodata};
  if (val(h.x,0,0,0) == nodata)
    return n;
  int ori = orientation(val(h.x,0,0,0));
  if (val(h.x,0,-1,0) != nodata && orientation(val(h.x,0,-1,0)) == ori) {
    if (val(h.x,0,1,0) != nodata && orientation(val(h.x,0,1,0)) == ori)
      n.y = (val(h.x,0,-1,0) - val(h.x,0,1,0))/2.;
    else
      n.y = val(h.x,0,-1,0) - val(h.x,0,0,0);
  }
  else if (val(h.x,0,1,0) != nodata && orientation(val(h.x,0,1,0)) == ori)
    n.y = val(h.x,0,0,0) - val(h.x,0,1,0);
  else
    return n;
  double nn = (ori ? -1. : 1.)*sqrt(1. + sq(n.y));
  n.y /= nn;
  n.x = 1./nn;
  return n;

#if _call_normal_x
}
#define _IN_STENCIL 1

#line 82
static coord _normal_x (Point point, vector h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 83 "/home/fpl/softwares/basilisk/src/curvature.h"

  coord n = {nodata, nodata, nodata};
  IF (_stencil_val(__FILE__,__LINE__,h.x,0,0,0) == nodata)
    return n;
  int ori = orientation(_stencil_val(__FILE__,__LINE__,h.x,0,0,0));
  IF (_stencil_val(__FILE__,__LINE__,h.x,-1,0,0) != nodata && orientation(_stencil_val(__FILE__,__LINE__,h.x,-1,0,0)) == ori) {
    IF (_stencil_val(__FILE__,__LINE__,h.x,1,0,0) != nodata && orientation(_stencil_val(__FILE__,__LINE__,h.x,1,0,0)) == ori)
      n.y = (_stencil_val(__FILE__,__LINE__,h.x,-1,0,0) - _stencil_val(__FILE__,__LINE__,h.x,1,0,0))/2.;
    
      n.y = _stencil_val(__FILE__,__LINE__,h.x,-1,0,0) - _stencil_val(__FILE__,__LINE__,h.x,0,0,0);
  }
  IF (_stencil_val(__FILE__,__LINE__,h.x,1,0,0) != nodata && orientation(_stencil_val(__FILE__,__LINE__,h.x,1,0,0)) == ori)
    n.y = _stencil_val(__FILE__,__LINE__,h.x,0,0,0) - _stencil_val(__FILE__,__LINE__,h.x,1,0,0);
  
    return n;
  double nn = (ori ? -1. : 1.)*sqrt(1. + sq(n.y));
  n.y /= nn;
  n.x = 1./nn;
  return n;

#undef _IN_STENCIL

#endif

#line 102
}
#line 179 "/home/fpl/softwares/basilisk/src/curvature.h"
static double height_curvature (Point point, scalar c, vector h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 180 "/home/fpl/softwares/basilisk/src/curvature.h"







  typedef struct {
    double n;
    double (* kappa) (Point, vector);
  } NormKappa;
  struct { NormKappa x, y, z; } n;
  {
#line 192

    n.x.n = val(c,1,0,0) - val(c,-1,0,0), n.x.kappa = kappa_x;
#line 192

    n.y.n = val(c,0,1,0) - val(c,0,-1,0), n.y.kappa = kappa_y;}
  double (* kappaf) (Point, vector) = NULL; NOT_UNUSED (kappaf);




  if (fabs(n.x.n) < fabs(n.y.n))
    swap (NormKappa, n.x, n.y);
#line 211 "/home/fpl/softwares/basilisk/src/curvature.h"
  double kappa = nodata;
  {
#line 212

    if (kappa == nodata) {
      kappa = n.x.kappa (point, h);
      if (kappa != nodata) {
 kappaf = n.x.kappa;
 if (n.x.n < 0.)
   kappa = - kappa;
      }
    }
#line 212

    if (kappa == nodata) {
      kappa = n.y.kappa (point, h);
      if (kappa != nodata) {
 kappaf = n.y.kappa;
 if (n.y.n < 0.)
   kappa = - kappa;
      }
    }}

  if (kappa != nodata) {




    if (fabs(kappa) > 1./Delta)
      kappa = sign(kappa)/Delta;
#line 247 "/home/fpl/softwares/basilisk/src/curvature.h"
  }

  return kappa;

#if _call_height_curvature
}
#define _IN_STENCIL 1

#line 179
static double _height_curvature (Point point, scalar c, vector h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 180 "/home/fpl/softwares/basilisk/src/curvature.h"







  typedef struct {
    double n;
    double (* kappa) (Point, vector);
  } NormKappa;
  struct { NormKappa x, y, z; } n;
  {
#line 192

    n.x.n = _stencil_val(__FILE__,__LINE__,c,1,0,0) - _stencil_val(__FILE__,__LINE__,c,-1,0,0), n.x.kappa = kappa_x;
#line 192

    n.y.n = _stencil_val(__FILE__,__LINE__,c,0,1,0) - _stencil_val(__FILE__,__LINE__,c,0,-1,0), n.y.kappa = kappa_y;}
  double (* kappaf) (Point, vector) = NULL; NOT_UNUSED (kappaf);




  IF (fabs(n.x.n) < fabs(n.y.n))
    swap (NormKappa, n.x, n.y);
#line 211 "/home/fpl/softwares/basilisk/src/curvature.h"
  double kappa = nodata;
  {
#line 212

    IF (kappa == nodata) {
      kappa = n.x.kappa (point, h);
      IF (kappa != nodata) {
 kappaf = n.x.kappa;
 IF (n.x.n < 0.)
   kappa = - kappa;
      }
    }
#line 212

    IF (kappa == nodata) {
      kappa = n.y.kappa (point, h);
      IF (kappa != nodata) {
 kappaf = n.y.kappa;
 IF (n.y.n < 0.)
   kappa = - kappa;
      }
    }}

  IF (kappa != nodata) {




    IF (fabs(kappa) > 1./Delta)
      kappa = sign(kappa)/Delta;
#line 247 "/home/fpl/softwares/basilisk/src/curvature.h"
  }

  return kappa;

#undef _IN_STENCIL

#endif

#line 250
}






coord height_normal (Point point, scalar c, vector h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 258 "/home/fpl/softwares/basilisk/src/curvature.h"







  typedef struct {
    double n;
    coord (* normal) (Point, vector);
  } NormNormal;
  struct { NormNormal x, y, z; } n;
  {
#line 270

    n.x.n = val(c,1,0,0) - val(c,-1,0,0), n.x.normal = normal_x;
#line 270

    n.y.n = val(c,0,1,0) - val(c,0,-1,0), n.y.normal = normal_y;}




  if (fabs(n.x.n) < fabs(n.y.n))
    swap (NormNormal, n.x, n.y);
#line 288 "/home/fpl/softwares/basilisk/src/curvature.h"
  coord normal = {nodata, nodata, nodata};
  {
#line 289

    if (normal.x == nodata)
      normal = n.x.normal (point, h);
#line 289

    if (normal.y == nodata)
      normal = n.y.normal (point, h);}

  return normal;

#if _call_height_normal
}
#define _IN_STENCIL 1

#line 257
static coord _height_normal (Point point, scalar c, vector h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 258 "/home/fpl/softwares/basilisk/src/curvature.h"







  typedef struct {
    double n;
    coord (* normal) (Point, vector);
  } NormNormal;
  struct { NormNormal x, y, z; } n;
  {
#line 270

    n.x.n = _stencil_val(__FILE__,__LINE__,c,1,0,0) - _stencil_val(__FILE__,__LINE__,c,-1,0,0), n.x.normal = normal_x;
#line 270

    n.y.n = _stencil_val(__FILE__,__LINE__,c,0,1,0) - _stencil_val(__FILE__,__LINE__,c,0,-1,0), n.y.normal = normal_y;}




  IF (fabs(n.x.n) < fabs(n.y.n))
    swap (NormNormal, n.x, n.y);
#line 288 "/home/fpl/softwares/basilisk/src/curvature.h"
  coord normal = {nodata, nodata, nodata};
  {
#line 289

    IF (normal.x == nodata)
      normal = n.x.normal (point, h);
#line 289

    IF (normal.y == nodata)
      normal = n.y.normal (point, h);}

  return normal;

#undef _IN_STENCIL

#endif

#line 294
}
#line 330 "/home/fpl/softwares/basilisk/src/curvature.h"
#line 1 "parabola.h"
#line 1 "/home/fpl/softwares/basilisk/src/parabola.h"
#line 1 "utils.h"
#line 2 "/home/fpl/softwares/basilisk/src/parabola.h"






typedef struct {
  coord o;

  coord m;
  double ** M, rhs[3], a[3];
#line 21 "/home/fpl/softwares/basilisk/src/parabola.h"
} ParabolaFit;

static void parabola_fit_init (ParabolaFit * p, coord o, coord m)
{
  {
#line 25

    p->o.x = o.x;
#line 25

    p->o.y = o.y;}

  {
#line 28

    p->m.x = m.x;
#line 28

    p->m.y = m.y;}
  normalize (&p->m);
  int n = 3;
#line 65 "/home/fpl/softwares/basilisk/src/parabola.h"
  p->M = (double **) matrix_new (n, n, sizeof(double));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      p->M[i][j] = 0.;
    p->rhs[i] = 0.;
  }
}

static void parabola_fit_add (ParabolaFit * p, coord m, double w)
{

  double x1 = m.x - p->o.x, y1 = m.y - p->o.y;
  double x = p->m.y*x1 - p->m.x*y1;
  double y = p->m.x*x1 + p->m.y*y1;
  double x2 = w*x*x, x3 = x2*x, x4 = x3*x;
  p->M[0][0] += x4;
  p->M[1][0] += x3; p->M[1][1] += x2;
  p->M[2][1] += w*x; p->M[2][2] += w;
  p->rhs[0] += x2*y; p->rhs[1] += w*x*y; p->rhs[2] += w*y;
#line 111 "/home/fpl/softwares/basilisk/src/parabola.h"
}

static double parabola_fit_solve (ParabolaFit * p)
{

  p->M[0][1] = p->M[1][0];
  p->M[0][2] = p->M[2][0] = p->M[1][1];
  p->M[1][2] = p->M[2][1];
  double pivmin = matrix_inverse (p->M, 3, 1e-10);
  if (pivmin) {
    p->a[0] = p->M[0][0]*p->rhs[0] + p->M[0][1]*p->rhs[1] + p->M[0][2]*p->rhs[2];
    p->a[1] = p->M[1][0]*p->rhs[0] + p->M[1][1]*p->rhs[1] + p->M[1][2]*p->rhs[2];
  }
  else
    p->a[0] = p->a[1] = 0.;
#line 158 "/home/fpl/softwares/basilisk/src/parabola.h"
  matrix_free (p->M);
  return pivmin;
}

static double parabola_fit_curvature (ParabolaFit * p,
          double kappamax, double * kmax)
{
  double kappa;

  double dnm = 1. + sq(p->a[1]);
  kappa = - 2.*p->a[0]/pow(dnm, 3/2.);
  if (kmax)
    *kmax = fabs (kappa);
#line 190 "/home/fpl/softwares/basilisk/src/parabola.h"
  if (fabs (kappa) > kappamax) {
    if (kmax)
      *kmax = kappamax;
    return kappa > 0. ? kappamax : - kappamax;
  }
  return kappa;
}
#line 331 "/home/fpl/softwares/basilisk/src/curvature.h"






static int independents (coord * p, int n)
{
  if (n < 2)
    return n;
  int ni = 1;
  for (int j = 1; j < n; j++) {
    bool depends = false;
    for (int i = 0; i < j && !depends; i++) {
      double d2 = 0.;
      {
#line 346

 d2 += sq(p[i].x - p[j].x);
#line 346

 d2 += sq(p[i].y - p[j].y);}
      depends = (d2 < sq(0.5));
    }
    ni += !depends;
  }
  return ni;
}






static double height_curvature_fit (Point point, scalar c, vector h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 361 "/home/fpl/softwares/basilisk/src/curvature.h"






  coord ip[2 == 2 ? 6 : 27];
  int n = 0;




  {
#line 373
 {





    int n1 = 0, n2 = 0;

    for (int i = -1; i <= 1; i++)
      if (val(h.y,i,0,0) != nodata) {
 if (orientation(val(h.y,i,0,0))) n1++; else n2++;
      }







    int ori = (n1 > n2);







    for (int i = -1; i <= 1; i++)
      if (val(h.y,i,0,0) != nodata && orientation(val(h.y,i,0,0)) == ori)
 ip[n].x = i, ip[n++].y = height(val(h.y,i,0,0));






  }
#line 373
 {





    int n1 = 0, n2 = 0;

    for (int i = -1; i <= 1; i++)
      if (val(h.x,0,i,0) != nodata) {
 if (orientation(val(h.x,0,i,0))) n1++; else n2++;
      }







    int ori = (n1 > n2);







    for (int i = -1; i <= 1; i++)
      if (val(h.x,0,i,0) != nodata && orientation(val(h.x,0,i,0)) == ori)
 ip[n].y = i, ip[n++].x = height(val(h.x,0,i,0));






  }}





  if (independents (ip, n) < (2 == 2 ? 3 : 9))
    return nodata;





  coord m = mycs (point, c), fc;
  double alpha = line_alpha (val(c,0,0,0), m);
  double area = line_length_center(m,alpha,&fc);
  ParabolaFit fit;
  parabola_fit_init (&fit, fc, m);

  NOT_UNUSED(area);
  parabola_fit_add (&fit, fc, .1);
#line 438 "/home/fpl/softwares/basilisk/src/curvature.h"
  for (int i = 0; i < n; i++)
    parabola_fit_add (&fit, ip[i], 1.);
  parabola_fit_solve (&fit);
  double kappa = parabola_fit_curvature (&fit, 2., NULL)/Delta;



  return kappa;

#if _call_height_curvature_fit
}
#define _IN_STENCIL 1
#define mycs _mycs

#line 360
static double _height_curvature_fit (Point point, scalar c, vector h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 361 "/home/fpl/softwares/basilisk/src/curvature.h"






  coord ip[2 == 2 ? 6 : 27];
  int n = 0;




  {
#line 373
 {





    int n1 = 0, n2 = 0;

    for (int i = -1; i <= 1; i++)
      IF (_stencil_val(__FILE__,__LINE__,h.y,i,0,0) != nodata) {
 IF (orientation(_stencil_val(__FILE__,__LINE__,h.y,i,0,0))) n1++;  n2++;
      }







    int ori = (n1 > n2);







    for (int i = -1; i <= 1; i++)
      IF (_stencil_val(__FILE__,__LINE__,h.y,i,0,0) != nodata && orientation(_stencil_val(__FILE__,__LINE__,h.y,i,0,0)) == ori)
 ip[n].x = i, ip[n++].y = height(_stencil_val(__FILE__,__LINE__,h.y,i,0,0));






  }
#line 373
 {





    int n1 = 0, n2 = 0;

    for (int i = -1; i <= 1; i++)
      IF (_stencil_val(__FILE__,__LINE__,h.x,0,i,0) != nodata) {
 IF (orientation(_stencil_val(__FILE__,__LINE__,h.x,0,i,0))) n1++;  n2++;
      }







    int ori = (n1 > n2);







    for (int i = -1; i <= 1; i++)
      IF (_stencil_val(__FILE__,__LINE__,h.x,0,i,0) != nodata && orientation(_stencil_val(__FILE__,__LINE__,h.x,0,i,0)) == ori)
 ip[n].y = i, ip[n++].x = height(_stencil_val(__FILE__,__LINE__,h.x,0,i,0));






  }}





  IF (independents (ip, n) < (2 == 2 ? 3 : 9))
    return nodata;





  coord m = mycs (point, c), fc;
  double alpha = line_alpha (_stencil_val(__FILE__,__LINE__,c,0,0,0), m);
  double area = line_length_center(m,alpha,&fc);
  ParabolaFit fit;
  parabola_fit_init (&fit, fc, m);

  NOT_UNUSED(area);
  parabola_fit_add (&fit, fc, .1);
#line 438 "/home/fpl/softwares/basilisk/src/curvature.h"
  for (int i = 0; i < n; i++)
    parabola_fit_add (&fit, ip[i], 1.);
  parabola_fit_solve (&fit);
  double kappa = parabola_fit_curvature (&fit, 2., NULL)/Delta;



  return kappa;

#undef mycs
#undef _IN_STENCIL

#endif

#line 446
}






static double centroids_curvature_fit (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 454 "/home/fpl/softwares/basilisk/src/curvature.h"






  coord m = mycs (point, c), fc;
  double alpha = line_alpha (val(c,0,0,0), m);
  line_length_center(m,alpha,&fc);
  ParabolaFit fit;
  parabola_fit_init (&fit, fc, m);





  coord r = {x,y,z};
   { foreach_neighbor(1)
    if (val(c,0,0,0) > 0. && val(c,0,0,0) < 1.) {
      coord m = mycs (point, c), fc;
      double alpha = line_alpha (val(c,0,0,0), m);
      double area = line_length_center(m,alpha,&fc);
      coord rn = {x,y,z};
      {
#line 477

 fc.x += (rn.x - r.x)/Delta;
#line 477

 fc.y += (rn.y - r.y)/Delta;}
      parabola_fit_add (&fit, fc, area);
    } end_foreach_neighbor(); }
  parabola_fit_solve (&fit);
  double kappa = parabola_fit_curvature (&fit, 2., NULL)/Delta;



  return kappa;

#if _call_centroids_curvature_fit
}
#define _IN_STENCIL 1
#define mycs _mycs

#line 453
static double _centroids_curvature_fit (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 454 "/home/fpl/softwares/basilisk/src/curvature.h"






  coord m = mycs (point, c), fc;
  double alpha = line_alpha (_stencil_val(__FILE__,__LINE__,c,0,0,0), m);
  line_length_center(m,alpha,&fc);
  ParabolaFit fit;
  parabola_fit_init (&fit, fc, m);





  coord r = {x,y,z};
   { foreach_neighbor(1)
    IF (_stencil_val(__FILE__,__LINE__,c,0,0,0) > 0. && _stencil_val(__FILE__,__LINE__,c,0,0,0) < 1.) {
      coord m = mycs (point, c), fc;
      double alpha = line_alpha (_stencil_val(__FILE__,__LINE__,c,0,0,0), m);
      double area = line_length_center(m,alpha,&fc);
      coord rn = {x,y,z};
      {
#line 477

 fc.x += (rn.x - r.x)/Delta;
#line 477

 fc.y += (rn.y - r.y)/Delta;}
      parabola_fit_add (&fit, fc, area);
    } end_foreach_neighbor(); }
  parabola_fit_solve (&fit);
  double kappa = parabola_fit_curvature (&fit, 2., NULL)/Delta;



  return kappa;

#undef mycs
#undef _IN_STENCIL

#endif

#line 487
}
#line 500 "/home/fpl/softwares/basilisk/src/curvature.h"
static inline bool interfacial (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 501 "/home/fpl/softwares/basilisk/src/curvature.h"

  if (val(c,0,0,0) >= 1.) {
    for (int i = -1; i <= 1; i += 2)
      {
#line 504

 if (val(c,i,0,0) <= 0.)
   return true;
#line 504

 if (val(c,0,i,0) <= 0.)
   return true;}
  }
  else if (val(c,0,0,0) <= 0.) {
    for (int i = -1; i <= 1; i += 2)
      {
#line 510

 if (val(c,i,0,0) >= 1.)
   return true;
#line 510

 if (val(c,0,i,0) >= 1.)
   return true;}
  }
  else
    return true;
  return false;

#if _call_interfacial
}
#define _IN_STENCIL 1

#line 500
static bool _interfacial (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 501 "/home/fpl/softwares/basilisk/src/curvature.h"

  IF (_stencil_val(__FILE__,__LINE__,c,0,0,0) >= 1.) {
    for (int i = -1; i <= 1; i += 2)
      {
#line 504

 IF (_stencil_val(__FILE__,__LINE__,c,i,0,0) <= 0.)
   return true;
#line 504

 IF (_stencil_val(__FILE__,__LINE__,c,0,i,0) <= 0.)
   return true;}
  }
  IF (_stencil_val(__FILE__,__LINE__,c,0,0,0) <= 0.) {
    for (int i = -1; i <= 1; i += 2)
      {
#line 510

 IF (_stencil_val(__FILE__,__LINE__,c,i,0,0) >= 1.)
   return true;
#line 510

 IF (_stencil_val(__FILE__,__LINE__,c,0,i,0) >= 1.)
   return true;}
  }
  
    return true;
  return false;

#undef _IN_STENCIL

#endif

#line 517
}
#line 530 "/home/fpl/softwares/basilisk/src/curvature.h"
typedef struct {
  int h;
  int f;
  int a;
  int c;
} cstats;

struct Curvature {
  scalar c, kappa;
  double sigma;
  bool add;
};


cstats curvature (struct Curvature p)
{ trace ("curvature", "/home/fpl/softwares/basilisk/src/curvature.h", 545);
  scalar c = p.c, kappa = p.kappa;
  double sigma = p.sigma ? p.sigma : 1.;
  int sh = 0, sf = 0, sa = 0, sc = 0;
  vector ch = _attribute[c.i].height, h = (ch).x.i ? (ch) : new_vector("h");
  if (!ch.x.i)
    heights (c, h);






  _attribute[kappa.i].refine = _attribute[kappa.i].prolongation = curvature_prolongation;
  _attribute[kappa.i].restriction = curvature_restriction;






  scalar k= new_scalar("k");
  scalar_clone (k, kappa);

   { 
#define height_curvature _height_curvature
#define height_curvature_fit _height_curvature_fit
#define interfacial _interfacial
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  int _sh = sh;
 int _sf = sf;
{ int sh = _sh; NOT_UNUSED(sh);
 int sf = _sf; NOT_UNUSED(sf);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/curvature.h", .line = 569,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 569 "/home/fpl/softwares/basilisk/src/curvature.h"
 {




    IF (!interfacial (point, c))
      _stencil_val(__FILE__,__LINE__,k,0,0,0) = nodata;





    IF ((_stencil_val(__FILE__,__LINE__,k,0,0,0) = height_curvature (point, c, h)) != nodata)
      sh++;
    IF ((_stencil_val(__FILE__,__LINE__,k,0,0,0) = height_curvature_fit (point, c, h)) != nodata)
      sf++;
  } } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#undef height_curvature
#undef height_curvature_fit
#undef interfacial
#line 585

#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:sh)  reduction(+:sf)) {

#line 569
foreach(){

#line 569 "/home/fpl/softwares/basilisk/src/curvature.h"
 {




    if (!interfacial (point, c))
      val(k,0,0,0) = nodata;





    else if ((val(k,0,0,0) = height_curvature (point, c, h)) != nodata)
      sh++;
    else if ((val(k,0,0,0) = height_curvature_fit (point, c, h)) != nodata)
      sf++;
  } } end_foreach();mpi_all_reduce_array (&sh, int, MPI_SUM, 1);
mpi_all_reduce_array (&sf, int, MPI_SUM, 1);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 585
 }

   { 
#define centroids_curvature_fit _centroids_curvature_fit
#define interfacial _interfacial
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _sigma = sigma;
 int _sa = sa;
 int _sc = sc;
{ double sigma = _sigma; NOT_UNUSED(sigma);
 int sa = _sa; NOT_UNUSED(sa);
 int sc = _sc; NOT_UNUSED(sc);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/curvature.h", .line = 587,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 587 "/home/fpl/softwares/basilisk/src/curvature.h"
 {





    double kf;
    IF (_stencil_val(__FILE__,__LINE__,k,0,0,0) < nodata)
      kf = _stencil_val(__FILE__,__LINE__,k,0,0,0);
    IF (interfacial (point, c)) {





      double sk = 0., a = 0.;
       { foreach_neighbor(1)
 IF (_stencil_val(__FILE__,__LINE__,k,0,0,0) < nodata)
   sk += _stencil_val(__FILE__,__LINE__,k,0,0,0), a++; end_foreach_neighbor(); }
      IF (a > 0.)
 kf = sk/a, sa++;
      




 kf = centroids_curvature_fit (point, c), sc++;
    }
    
      kf = nodata;




    IF (kf == nodata)
      _stencil_val(__FILE__,__LINE__,kappa,0,0,0) = nodata;
    IF (p.add)
      _stencil_val(__FILE__,__LINE__,kappa,0,0,0) += sigma*kf;
    
      _stencil_val(__FILE__,__LINE__,kappa,0,0,0) = sigma*kf;
  } } end_foreach_stencil(); if (_first_call) {
 if (sigma != _sigma)
   reduction_warning ("/home/fpl/softwares/basilisk/src/curvature.h", 587, "sigma");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#undef centroids_curvature_fit
#undef interfacial
#line 627

#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:sa)  reduction(+:sc)) {

#line 587
foreach (){

#line 587 "/home/fpl/softwares/basilisk/src/curvature.h"
 {





    double kf;
    if (val(k,0,0,0) < nodata)
      kf = val(k,0,0,0);
    else if (interfacial (point, c)) {





      double sk = 0., a = 0.;
       { foreach_neighbor(1)
 if (val(k,0,0,0) < nodata)
   sk += val(k,0,0,0), a++; end_foreach_neighbor(); }
      if (a > 0.)
 kf = sk/a, sa++;
      else




 kf = centroids_curvature_fit (point, c), sc++;
    }
    else
      kf = nodata;




    if (kf == nodata)
      val(kappa,0,0,0) = nodata;
    else if (p.add)
      val(kappa,0,0,0) += sigma*kf;
    else
      val(kappa,0,0,0) = sigma*kf;
  } } end_foreach();mpi_all_reduce_array (&sa, int, MPI_SUM, 1);
mpi_all_reduce_array (&sc, int, MPI_SUM, 1);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 627
 }

  { cstats _ret =  (cstats){sh, sf, sa, sc}; delete (((scalar []){k,{-1}}));  { strongif (!(ch).x.i) delete (((scalar []){h.x,h.y,{-1}})); }  end_trace("curvature", "/home/fpl/softwares/basilisk/src/curvature.h", 629);  return _ret; }
 delete (((scalar []){k,{-1}}));  { strongif (!(ch).x.i) delete (((scalar []){h.x,h.y,{-1}})); }  end_trace("curvature", "/home/fpl/softwares/basilisk/src/curvature.h", 630); }
#line 649 "/home/fpl/softwares/basilisk/src/curvature.h"

#line 649

static double pos_x (Point point, vector h, coord * G, coord * Z)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 651 "/home/fpl/softwares/basilisk/src/curvature.h"

  if (fabs(height(val(h.x,0,0,0))) > 1.)
    return nodata;
  coord o = {x, y, z};
  o.x += height(val(h.x,0,0,0))*Delta;
  double pos = 0.;
  {
#line 657

    pos += (o.x - Z->x)*G->x;
#line 657

    pos += (o.y - Z->y)*G->y;}
  return pos;

#if _call_pos_x
}
#define _IN_STENCIL 1

#line 650
static double _pos_x (Point point, vector h, coord * G, coord * Z)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 651 "/home/fpl/softwares/basilisk/src/curvature.h"

  IF (fabs(height(_stencil_val(__FILE__,__LINE__,h.x,0,0,0))) > 1.)
    return nodata;
  coord o = {x, y, z};
  o.x += height(_stencil_val(__FILE__,__LINE__,h.x,0,0,0))*Delta;
  double pos = 0.;
  {
#line 657

    pos += (o.x - Z->x)*G->x;
#line 657

    pos += (o.y - Z->y)*G->y;}
  return pos;

#undef _IN_STENCIL

#endif

#line 660
}
#line 649

static double pos_y (Point point, vector h, coord * G, coord * Z)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 651 "/home/fpl/softwares/basilisk/src/curvature.h"

  if (fabs(height(val(h.y,0,0,0))) > 1.)
    return nodata;
  coord o = {x, y, z};
  o.y += height(val(h.y,0,0,0))*Delta;
  double pos = 0.;
  {
#line 657

    pos += (o.y - Z->y)*G->y;
#line 657

    pos += (o.x - Z->x)*G->x;}
  return pos;

#if _call_pos_y
}
#define _IN_STENCIL 1

#line 650
static double _pos_y (Point point, vector h, coord * G, coord * Z)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 651 "/home/fpl/softwares/basilisk/src/curvature.h"

  IF (fabs(height(_stencil_val(__FILE__,__LINE__,h.y,0,0,0))) > 1.)
    return nodata;
  coord o = {x, y, z};
  o.y += height(_stencil_val(__FILE__,__LINE__,h.y,0,0,0))*Delta;
  double pos = 0.;
  {
#line 657

    pos += (o.y - Z->y)*G->y;
#line 657

    pos += (o.x - Z->x)*G->x;}
  return pos;

#undef _IN_STENCIL

#endif

#line 660
}







static double height_position (Point point, scalar f, vector h,
          coord * G, coord * Z)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 670 "/home/fpl/softwares/basilisk/src/curvature.h"







  typedef struct {
    double n;
    double (* pos) (Point, vector, coord *, coord *);
  } NormPos;
  struct { NormPos x, y, z; } n;
  {
#line 682

    n.x.n = val(f,1,0,0) - val(f,-1,0,0), n.x.pos = pos_x;
#line 682

    n.y.n = val(f,0,1,0) - val(f,0,-1,0), n.y.pos = pos_y;}




  if (fabs(n.x.n) < fabs(n.y.n))
    swap (NormPos, n.x, n.y);
#line 700 "/home/fpl/softwares/basilisk/src/curvature.h"
  double pos = nodata;
  {
#line 701

    if (pos == nodata)
      pos = n.x.pos (point, h, G, Z);
#line 701

    if (pos == nodata)
      pos = n.y.pos (point, h, G, Z);}

  return pos;

#if _call_height_position
}
#define _IN_STENCIL 1
#define pos_x _pos_x
#define pos_y _pos_y

#line 668
static double _height_position (Point point, scalar f, vector h,
          coord * G, coord * Z)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 670 "/home/fpl/softwares/basilisk/src/curvature.h"







  typedef struct {
    double n;
    double (* pos) (Point, vector, coord *, coord *);
  } NormPos;
  struct { NormPos x, y, z; } n;
  {
#line 682

    n.x.n = _stencil_val(__FILE__,__LINE__,f,1,0,0) - _stencil_val(__FILE__,__LINE__,f,-1,0,0), n.x.pos = pos_x;
#line 682

    n.y.n = _stencil_val(__FILE__,__LINE__,f,0,1,0) - _stencil_val(__FILE__,__LINE__,f,0,-1,0), n.y.pos = pos_y;}




  IF (fabs(n.x.n) < fabs(n.y.n))
    swap (NormPos, n.x, n.y);
#line 700 "/home/fpl/softwares/basilisk/src/curvature.h"
  double pos = nodata;
  {
#line 701

    IF (pos == nodata)
      pos = n.x.pos (point, h, G, Z);
#line 701

    IF (pos == nodata)
      pos = n.y.pos (point, h, G, Z);}

  return pos;

#undef pos_x
#undef pos_y
#undef _IN_STENCIL

#endif

#line 706
}
#line 717 "/home/fpl/softwares/basilisk/src/curvature.h"
struct Position {
  scalar f, pos;
  coord G, Z;
  bool add;
};

void position (struct Position p)
{
  scalar f = p.f, pos = p.pos;
  coord * G = &p.G, * Z = &p.Z;






  _attribute[pos.i].refine = _attribute[pos.i].prolongation = curvature_prolongation;
  _attribute[pos.i].restriction = curvature_restriction;


  vector fh = _attribute[f.i].height, h = (fh).x.i ? (fh) : new_vector("h");
  if (!fh.x.i)
    heights (f, h);
   { 
#define mycs _mycs
#define interfacial _interfacial
#define height_position _height_position
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  coord _G[10] = {0};
 coord _Z[10] = {0};
{ coord * G = _G; NOT_UNUSED(G);
 coord * Z = _Z; NOT_UNUSED(Z);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/curvature.h", .line = 740,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 740 "/home/fpl/softwares/basilisk/src/curvature.h"
 {
    IF (interfacial (point, f)) {
      double hp = height_position (point, f, h, G, Z);
      IF (hp == nodata) {





 coord n = mycs (point, f), o = {x,y,z}, c;
 double alpha = line_alpha (_stencil_val(__FILE__,__LINE__,f,0,0,0), n);
 line_length_center(n,alpha,&c);
 hp = 0.;
 {
#line 753

   hp += (o.x + Delta*c.x - Z->x)*G->x;
#line 753

   hp += (o.y + Delta*c.y - Z->y)*G->y;}
      }
      IF (p.add)
 _stencil_val(__FILE__,__LINE__,pos,0,0,0) += hp;
      
 _stencil_val(__FILE__,__LINE__,pos,0,0,0) = hp;
    }
    
      _stencil_val(__FILE__,__LINE__,pos,0,0,0) = nodata;
  } } end_foreach_stencil(); if (_first_call) {
 for (int i = 0; i < (10*sizeof(coord)); i++)
   if (((char *)_G)[i] != 0) {
     reduction_warning ("/home/fpl/softwares/basilisk/src/curvature.h", 740, "G");
     break; }
 }
 if (_first_call) {
 for (int i = 0; i < (10*sizeof(coord)); i++)
   if (((char *)_Z)[i] != 0) {
     reduction_warning ("/home/fpl/softwares/basilisk/src/curvature.h", 740, "Z");
     break; }
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#undef mycs
#undef interfacial
#undef height_position
#line 763
foreach(){

#line 740 "/home/fpl/softwares/basilisk/src/curvature.h"
 {
    if (interfacial (point, f)) {
      double hp = height_position (point, f, h, G, Z);
      if (hp == nodata) {





 coord n = mycs (point, f), o = {x,y,z}, c;
 double alpha = line_alpha (val(f,0,0,0), n);
 line_length_center(n,alpha,&c);
 hp = 0.;
 {
#line 753

   hp += (o.x + Delta*c.x - Z->x)*G->x;
#line 753

   hp += (o.y + Delta*c.y - Z->y)*G->y;}
      }
      if (p.add)
 val(pos,0,0,0) += hp;
      else
 val(pos,0,0,0) = hp;
    }
    else
      val(pos,0,0,0) = nodata;
  } } end_foreach(); }
 { strongif (!(fh).x.i) delete (((scalar []){h.x,h.y,{-1}})); } }
#line 16 "/home/fpl/softwares/basilisk/src/contact.h"






coord interface_normal (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 23 "/home/fpl/softwares/basilisk/src/contact.h"

  coord n;
  if (!_attribute[c.i].height.x.i || (n = height_normal (point, c, _attribute[c.i].height)).x == nodata)
    n = mycs (point, c);
  return n;

#if _call_interface_normal
}
#define _IN_STENCIL 1
#define mycs _mycs
#define height_normal _height_normal

#line 22
static coord _interface_normal (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 23 "/home/fpl/softwares/basilisk/src/contact.h"

  coord n;
  IF (!_attribute[c.i].height.x.i || (n = height_normal (point, c, _attribute[c.i].height)).x == nodata)
    n = mycs (point, c);
  return n;

#undef mycs
#undef height_normal
#undef _IN_STENCIL

#endif

#line 28
}
#line 43 "/home/fpl/softwares/basilisk/src/contact.h"
extern scalar * interfaces;

static int init_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int init_0 (const int i, const double t, Event * _ev) { trace ("init_0", "/home/fpl/softwares/basilisk/src/contact.h", 45);  {
  strongif (interfaces) for (scalar c = *interfaces, *_i112 = interfaces; ((scalar *)&c)->i >= 0; c = *++_i112)
    if (_attribute[c.i].height.x.i)
      heights (c, _attribute[c.i].height);
 end_trace("init_0", "/home/fpl/softwares/basilisk/src/contact.h", 49); } return 0; } 

static int vof_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int vof_0 (const int i, const double t, Event * _ev) { trace ("vof_0", "/home/fpl/softwares/basilisk/src/contact.h", 51);  {
  strongif (interfaces) for (scalar c = *interfaces, *_i113 = interfaces; ((scalar *)&c)->i >= 0; c = *++_i113)
    if (_attribute[c.i].height.x.i)
      heights (c, _attribute[c.i].height);
 end_trace("vof_0", "/home/fpl/softwares/basilisk/src/contact.h", 55); } return 0; } 
#line 22 "60d_plate_advancing_simulation.c"
#line 1 "tension.h"
#line 1 "/home/fpl/softwares/basilisk/src/tension.h"
#line 15 "/home/fpl/softwares/basilisk/src/tension.h"
#line 1 "iforce.h"
#line 1 "/home/fpl/softwares/basilisk/src/iforce.h"
#line 20 "/home/fpl/softwares/basilisk/src/iforce.h"










static int defaults_1_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults_1 (const int i, const double t, Event * _ev) { trace ("defaults_1", "/home/fpl/softwares/basilisk/src/iforce.h", 30);  {
  if (is_constant(a.x)) {
    a = new_face_vector("a");
     { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/iforce.h", .line = 33,
    .each = "foreach_face", .first = _first_call
  };
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 33
{

#line 33 "/home/fpl/softwares/basilisk/src/iforce.h"

      _stencil_val(__FILE__,__LINE__,a.x,0,0,0) = 0.; }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 33
{

#line 33 "/home/fpl/softwares/basilisk/src/iforce.h"

      _stencil_val(__FILE__,__LINE__,a.y,0,0,0) = 0.; }  }}  end_foreach_face_stencil()
#line 34
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 34
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 33
{

#line 33 "/home/fpl/softwares/basilisk/src/iforce.h"

      val(a.x,0,0,0) = 0.; }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 33
{

#line 33 "/home/fpl/softwares/basilisk/src/iforce.h"

      val(a.y,0,0,0) = 0.; }  }}  end_foreach_face_generic()
#line 34
 end_foreach_face(); }
  }
 end_trace("defaults_1", "/home/fpl/softwares/basilisk/src/iforce.h", 36); } return 0; } 






static int acceleration_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int acceleration_0 (const int i, const double t, Event * _ev) { trace ("acceleration_0", "/home/fpl/softwares/basilisk/src/iforce.h", 43); 
{





  scalar * list = NULL;
  strongif (interfaces) for (scalar f = *interfaces, *_i114 = interfaces; ((scalar *)&f)->i >= 0; f = *++_i114)
    if (_attribute[f.i].phi.i) {
      list = list_add (list, f);






       { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/iforce.h", .line = 60,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 60 "/home/fpl/softwares/basilisk/src/iforce.h"

 _stencil_val(__FILE__,__LINE__,f,0,0,0) = clamp (_stencil_val(__FILE__,__LINE__,f,0,0,0), 0., 1.); } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 61
foreach(){

#line 60 "/home/fpl/softwares/basilisk/src/iforce.h"

 val(f,0,0,0) = clamp (val(f,0,0,0), 0., 1.); } end_foreach(); }
    }
#line 72 "/home/fpl/softwares/basilisk/src/iforce.h"
  strongif (list) for (scalar f = *list, *_i115 = list; ((scalar *)&f)->i >= 0; f = *++_i115) {
    _attribute[f.i].prolongation = _attribute[p.i].prolongation;
    _attribute[f.i].dirty = true;
  }
#line 86 "/home/fpl/softwares/basilisk/src/iforce.h"
  vector ia = a;
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/iforce.h", .line = 87,
    .each = "foreach_face", .first = _first_call
  };

strongif (!is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 87
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 87
{

#line 87 "/home/fpl/softwares/basilisk/src/iforce.h"

    strongif (list) for (scalar f = *list, *_i116 = list; ((scalar *)&f)->i >= 0; f = *++_i116)
      IF (_stencil_val(__FILE__,__LINE__,f,0,0,0) != _stencil_val(__FILE__,__LINE__,f,-1,0,0) && val_fm_x(fm.x,0,0,0) > 0.) {
#line 99 "/home/fpl/softwares/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (_stencil_val(__FILE__,__LINE__,phi,0,0,0) < nodata && _stencil_val(__FILE__,__LINE__,phi,-1,0,0) < nodata) ?
   (_stencil_val(__FILE__,__LINE__,phi,0,0,0) + _stencil_val(__FILE__,__LINE__,phi,-1,0,0))/2. :
   _stencil_val(__FILE__,__LINE__,phi,0,0,0) < nodata ? _stencil_val(__FILE__,__LINE__,phi,0,0,0) :
   _stencil_val(__FILE__,__LINE__,phi,-1,0,0) < nodata ? _stencil_val(__FILE__,__LINE__,phi,-1,0,0) :
   0.;

 _stencil_val(__FILE__,__LINE__,ia.x,0,0,0) += val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*phif*(_stencil_val(__FILE__,__LINE__,f,0,0,0) - _stencil_val(__FILE__,__LINE__,f,-1,0,0))/Delta;
      } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 87
{

#line 87 "/home/fpl/softwares/basilisk/src/iforce.h"

    strongif (list) for (scalar f = *list, *_i116 = list; ((scalar *)&f)->i >= 0; f = *++_i116)
      IF (_stencil_val(__FILE__,__LINE__,f,0,0,0) != _stencil_val(__FILE__,__LINE__,f,0,-1,0) && val_fm_y(fm.y,0,0,0) > 0.) {
#line 99 "/home/fpl/softwares/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (_stencil_val(__FILE__,__LINE__,phi,0,0,0) < nodata && _stencil_val(__FILE__,__LINE__,phi,0,-1,0) < nodata) ?
   (_stencil_val(__FILE__,__LINE__,phi,0,0,0) + _stencil_val(__FILE__,__LINE__,phi,0,-1,0))/2. :
   _stencil_val(__FILE__,__LINE__,phi,0,0,0) < nodata ? _stencil_val(__FILE__,__LINE__,phi,0,0,0) :
   _stencil_val(__FILE__,__LINE__,phi,0,-1,0) < nodata ? _stencil_val(__FILE__,__LINE__,phi,0,-1,0) :
   0.;

 _stencil_val(__FILE__,__LINE__,ia.y,0,0,0) += val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*phif*(_stencil_val(__FILE__,__LINE__,f,0,0,0) - _stencil_val(__FILE__,__LINE__,f,0,-1,0))/Delta;
      } }  }}  end_foreach_face_stencil()
#line 108
 }
strongif (is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 87
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 87
{

#line 87 "/home/fpl/softwares/basilisk/src/iforce.h"

    strongif (list) for (scalar f = *list, *_i116 = list; ((scalar *)&f)->i >= 0; f = *++_i116)
      IF (_stencil_val(__FILE__,__LINE__,f,0,0,0) != _stencil_val(__FILE__,__LINE__,f,-1,0,0) && val_fm_x(fm.x,0,0,0) > 0.) {
#line 99 "/home/fpl/softwares/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (_stencil_val(__FILE__,__LINE__,phi,0,0,0) < nodata && _stencil_val(__FILE__,__LINE__,phi,-1,0,0) < nodata) ?
   (_stencil_val(__FILE__,__LINE__,phi,0,0,0) + _stencil_val(__FILE__,__LINE__,phi,-1,0,0))/2. :
   _stencil_val(__FILE__,__LINE__,phi,0,0,0) < nodata ? _stencil_val(__FILE__,__LINE__,phi,0,0,0) :
   _stencil_val(__FILE__,__LINE__,phi,-1,0,0) < nodata ? _stencil_val(__FILE__,__LINE__,phi,-1,0,0) :
   0.;

 _stencil_val(__FILE__,__LINE__,ia.x,0,0,0) += val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*phif*(_stencil_val(__FILE__,__LINE__,f,0,0,0) - _stencil_val(__FILE__,__LINE__,f,-1,0,0))/Delta;
      } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 87
{

#line 87 "/home/fpl/softwares/basilisk/src/iforce.h"

    strongif (list) for (scalar f = *list, *_i116 = list; ((scalar *)&f)->i >= 0; f = *++_i116)
      IF (_stencil_val(__FILE__,__LINE__,f,0,0,0) != _stencil_val(__FILE__,__LINE__,f,0,-1,0) && val_fm_y(fm.y,0,0,0) > 0.) {
#line 99 "/home/fpl/softwares/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (_stencil_val(__FILE__,__LINE__,phi,0,0,0) < nodata && _stencil_val(__FILE__,__LINE__,phi,0,-1,0) < nodata) ?
   (_stencil_val(__FILE__,__LINE__,phi,0,0,0) + _stencil_val(__FILE__,__LINE__,phi,0,-1,0))/2. :
   _stencil_val(__FILE__,__LINE__,phi,0,0,0) < nodata ? _stencil_val(__FILE__,__LINE__,phi,0,0,0) :
   _stencil_val(__FILE__,__LINE__,phi,0,-1,0) < nodata ? _stencil_val(__FILE__,__LINE__,phi,0,-1,0) :
   0.;

 _stencil_val(__FILE__,__LINE__,ia.y,0,0,0) += val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*phif*(_stencil_val(__FILE__,__LINE__,f,0,0,0) - _stencil_val(__FILE__,__LINE__,f,0,-1,0))/Delta;
      } }  }}  end_foreach_face_stencil()
#line 108
 }
strongif (!is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 87
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 87
{

#line 87 "/home/fpl/softwares/basilisk/src/iforce.h"

    strongif (list) for (scalar f = *list, *_i116 = list; ((scalar *)&f)->i >= 0; f = *++_i116)
      IF (_stencil_val(__FILE__,__LINE__,f,0,0,0) != _stencil_val(__FILE__,__LINE__,f,-1,0,0) && val_fm_x(fm.x,0,0,0) > 0.) {
#line 99 "/home/fpl/softwares/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (_stencil_val(__FILE__,__LINE__,phi,0,0,0) < nodata && _stencil_val(__FILE__,__LINE__,phi,-1,0,0) < nodata) ?
   (_stencil_val(__FILE__,__LINE__,phi,0,0,0) + _stencil_val(__FILE__,__LINE__,phi,-1,0,0))/2. :
   _stencil_val(__FILE__,__LINE__,phi,0,0,0) < nodata ? _stencil_val(__FILE__,__LINE__,phi,0,0,0) :
   _stencil_val(__FILE__,__LINE__,phi,-1,0,0) < nodata ? _stencil_val(__FILE__,__LINE__,phi,-1,0,0) :
   0.;

 _stencil_val(__FILE__,__LINE__,ia.x,0,0,0) += val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*phif*(_stencil_val(__FILE__,__LINE__,f,0,0,0) - _stencil_val(__FILE__,__LINE__,f,-1,0,0))/Delta;
      } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 87
{

#line 87 "/home/fpl/softwares/basilisk/src/iforce.h"

    strongif (list) for (scalar f = *list, *_i116 = list; ((scalar *)&f)->i >= 0; f = *++_i116)
      IF (_stencil_val(__FILE__,__LINE__,f,0,0,0) != _stencil_val(__FILE__,__LINE__,f,0,-1,0) && val_fm_y(fm.y,0,0,0) > 0.) {
#line 99 "/home/fpl/softwares/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (_stencil_val(__FILE__,__LINE__,phi,0,0,0) < nodata && _stencil_val(__FILE__,__LINE__,phi,0,-1,0) < nodata) ?
   (_stencil_val(__FILE__,__LINE__,phi,0,0,0) + _stencil_val(__FILE__,__LINE__,phi,0,-1,0))/2. :
   _stencil_val(__FILE__,__LINE__,phi,0,0,0) < nodata ? _stencil_val(__FILE__,__LINE__,phi,0,0,0) :
   _stencil_val(__FILE__,__LINE__,phi,0,-1,0) < nodata ? _stencil_val(__FILE__,__LINE__,phi,0,-1,0) :
   0.;

 _stencil_val(__FILE__,__LINE__,ia.y,0,0,0) += val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*phif*(_stencil_val(__FILE__,__LINE__,f,0,0,0) - _stencil_val(__FILE__,__LINE__,f,0,-1,0))/Delta;
      } }  }}  end_foreach_face_stencil()
#line 108
 }
strongif (is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 87
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 87
{

#line 87 "/home/fpl/softwares/basilisk/src/iforce.h"

    strongif (list) for (scalar f = *list, *_i116 = list; ((scalar *)&f)->i >= 0; f = *++_i116)
      IF (_stencil_val(__FILE__,__LINE__,f,0,0,0) != _stencil_val(__FILE__,__LINE__,f,-1,0,0) && val_fm_x(fm.x,0,0,0) > 0.) {
#line 99 "/home/fpl/softwares/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (_stencil_val(__FILE__,__LINE__,phi,0,0,0) < nodata && _stencil_val(__FILE__,__LINE__,phi,-1,0,0) < nodata) ?
   (_stencil_val(__FILE__,__LINE__,phi,0,0,0) + _stencil_val(__FILE__,__LINE__,phi,-1,0,0))/2. :
   _stencil_val(__FILE__,__LINE__,phi,0,0,0) < nodata ? _stencil_val(__FILE__,__LINE__,phi,0,0,0) :
   _stencil_val(__FILE__,__LINE__,phi,-1,0,0) < nodata ? _stencil_val(__FILE__,__LINE__,phi,-1,0,0) :
   0.;

 _stencil_val(__FILE__,__LINE__,ia.x,0,0,0) += val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*phif*(_stencil_val(__FILE__,__LINE__,f,0,0,0) - _stencil_val(__FILE__,__LINE__,f,-1,0,0))/Delta;
      } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 87
{

#line 87 "/home/fpl/softwares/basilisk/src/iforce.h"

    strongif (list) for (scalar f = *list, *_i116 = list; ((scalar *)&f)->i >= 0; f = *++_i116)
      IF (_stencil_val(__FILE__,__LINE__,f,0,0,0) != _stencil_val(__FILE__,__LINE__,f,0,-1,0) && val_fm_y(fm.y,0,0,0) > 0.) {
#line 99 "/home/fpl/softwares/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (_stencil_val(__FILE__,__LINE__,phi,0,0,0) < nodata && _stencil_val(__FILE__,__LINE__,phi,0,-1,0) < nodata) ?
   (_stencil_val(__FILE__,__LINE__,phi,0,0,0) + _stencil_val(__FILE__,__LINE__,phi,0,-1,0))/2. :
   _stencil_val(__FILE__,__LINE__,phi,0,0,0) < nodata ? _stencil_val(__FILE__,__LINE__,phi,0,0,0) :
   _stencil_val(__FILE__,__LINE__,phi,0,-1,0) < nodata ? _stencil_val(__FILE__,__LINE__,phi,0,-1,0) :
   0.;

 _stencil_val(__FILE__,__LINE__,ia.y,0,0,0) += val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*phif*(_stencil_val(__FILE__,__LINE__,f,0,0,0) - _stencil_val(__FILE__,__LINE__,f,0,-1,0))/Delta;
      } }  }}  end_foreach_face_stencil()
#line 108
 }  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 108

strongif (!is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 87
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 87
{

#line 87 "/home/fpl/softwares/basilisk/src/iforce.h"

    strongif (list) for (scalar f = *list, *_i116 = list; ((scalar *)&f)->i >= 0; f = *++_i116)
      if (val(f,0,0,0) != val(f,-1,0,0) && val_fm_x(fm.x,0,0,0) > 0.) {
#line 99 "/home/fpl/softwares/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < nodata && val(phi,-1,0,0) < nodata) ?
   (val(phi,0,0,0) + val(phi,-1,0,0))/2. :
   val(phi,0,0,0) < nodata ? val(phi,0,0,0) :
   val(phi,-1,0,0) < nodata ? val(phi,-1,0,0) :
   0.;

 val(ia.x,0,0,0) += val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*phif*(val(f,0,0,0) - val(f,-1,0,0))/Delta;
      } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 87
{

#line 87 "/home/fpl/softwares/basilisk/src/iforce.h"

    strongif (list) for (scalar f = *list, *_i116 = list; ((scalar *)&f)->i >= 0; f = *++_i116)
      if (val(f,0,0,0) != val(f,0,-1,0) && val_fm_y(fm.y,0,0,0) > 0.) {
#line 99 "/home/fpl/softwares/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < nodata && val(phi,0,-1,0) < nodata) ?
   (val(phi,0,0,0) + val(phi,0,-1,0))/2. :
   val(phi,0,0,0) < nodata ? val(phi,0,0,0) :
   val(phi,0,-1,0) < nodata ? val(phi,0,-1,0) :
   0.;

 val(ia.y,0,0,0) += val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*phif*(val(f,0,0,0) - val(f,0,-1,0))/Delta;
      } }  }}  end_foreach_face_generic()
#line 108
 end_foreach_face(); }
strongif (is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 87
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 87
{

#line 87 "/home/fpl/softwares/basilisk/src/iforce.h"

    strongif (list) for (scalar f = *list, *_i116 = list; ((scalar *)&f)->i >= 0; f = *++_i116)
      if (val(f,0,0,0) != val(f,-1,0,0) && val_fm_x(fm.x,0,0,0) > 0.) {
#line 99 "/home/fpl/softwares/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < nodata && val(phi,-1,0,0) < nodata) ?
   (val(phi,0,0,0) + val(phi,-1,0,0))/2. :
   val(phi,0,0,0) < nodata ? val(phi,0,0,0) :
   val(phi,-1,0,0) < nodata ? val(phi,-1,0,0) :
   0.;

 val(ia.x,0,0,0) += val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*phif*(val(f,0,0,0) - val(f,-1,0,0))/Delta;
      } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 87
{

#line 87 "/home/fpl/softwares/basilisk/src/iforce.h"

    strongif (list) for (scalar f = *list, *_i116 = list; ((scalar *)&f)->i >= 0; f = *++_i116)
      if (val(f,0,0,0) != val(f,0,-1,0) && val_fm_y(fm.y,0,0,0) > 0.) {
#line 99 "/home/fpl/softwares/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < nodata && val(phi,0,-1,0) < nodata) ?
   (val(phi,0,0,0) + val(phi,0,-1,0))/2. :
   val(phi,0,0,0) < nodata ? val(phi,0,0,0) :
   val(phi,0,-1,0) < nodata ? val(phi,0,-1,0) :
   0.;

 val(ia.y,0,0,0) += val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*phif*(val(f,0,0,0) - val(f,0,-1,0))/Delta;
      } }  }}  end_foreach_face_generic()
#line 108
 end_foreach_face(); }
strongif (!is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 87
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 87
{

#line 87 "/home/fpl/softwares/basilisk/src/iforce.h"

    strongif (list) for (scalar f = *list, *_i116 = list; ((scalar *)&f)->i >= 0; f = *++_i116)
      if (val(f,0,0,0) != val(f,-1,0,0) && val_fm_x(fm.x,0,0,0) > 0.) {
#line 99 "/home/fpl/softwares/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < nodata && val(phi,-1,0,0) < nodata) ?
   (val(phi,0,0,0) + val(phi,-1,0,0))/2. :
   val(phi,0,0,0) < nodata ? val(phi,0,0,0) :
   val(phi,-1,0,0) < nodata ? val(phi,-1,0,0) :
   0.;

 val(ia.x,0,0,0) += val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*phif*(val(f,0,0,0) - val(f,-1,0,0))/Delta;
      } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 87
{

#line 87 "/home/fpl/softwares/basilisk/src/iforce.h"

    strongif (list) for (scalar f = *list, *_i116 = list; ((scalar *)&f)->i >= 0; f = *++_i116)
      if (val(f,0,0,0) != val(f,0,-1,0) && val_fm_y(fm.y,0,0,0) > 0.) {
#line 99 "/home/fpl/softwares/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < nodata && val(phi,0,-1,0) < nodata) ?
   (val(phi,0,0,0) + val(phi,0,-1,0))/2. :
   val(phi,0,0,0) < nodata ? val(phi,0,0,0) :
   val(phi,0,-1,0) < nodata ? val(phi,0,-1,0) :
   0.;

 val(ia.y,0,0,0) += val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*phif*(val(f,0,0,0) - val(f,0,-1,0))/Delta;
      } }  }}  end_foreach_face_generic()
#line 108
 end_foreach_face(); }
strongif (is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 87
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 87
{

#line 87 "/home/fpl/softwares/basilisk/src/iforce.h"

    strongif (list) for (scalar f = *list, *_i116 = list; ((scalar *)&f)->i >= 0; f = *++_i116)
      if (val(f,0,0,0) != val(f,-1,0,0) && val_fm_x(fm.x,0,0,0) > 0.) {
#line 99 "/home/fpl/softwares/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < nodata && val(phi,-1,0,0) < nodata) ?
   (val(phi,0,0,0) + val(phi,-1,0,0))/2. :
   val(phi,0,0,0) < nodata ? val(phi,0,0,0) :
   val(phi,-1,0,0) < nodata ? val(phi,-1,0,0) :
   0.;

 val(ia.x,0,0,0) += val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*phif*(val(f,0,0,0) - val(f,-1,0,0))/Delta;
      } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 87
{

#line 87 "/home/fpl/softwares/basilisk/src/iforce.h"

    strongif (list) for (scalar f = *list, *_i116 = list; ((scalar *)&f)->i >= 0; f = *++_i116)
      if (val(f,0,0,0) != val(f,0,-1,0) && val_fm_y(fm.y,0,0,0) > 0.) {
#line 99 "/home/fpl/softwares/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < nodata && val(phi,0,-1,0) < nodata) ?
   (val(phi,0,0,0) + val(phi,0,-1,0))/2. :
   val(phi,0,0,0) < nodata ? val(phi,0,0,0) :
   val(phi,0,-1,0) < nodata ? val(phi,0,-1,0) :
   0.;

 val(ia.y,0,0,0) += val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*phif*(val(f,0,0,0) - val(f,0,-1,0))/Delta;
      } }  }}  end_foreach_face_generic()
#line 108
 end_foreach_face(); } }






  strongif (list) for (scalar f = *list, *_i117 = list; ((scalar *)&f)->i >= 0; f = *++_i117) {
    _attribute[f.i].prolongation = fraction_refine;
    _attribute[f.i].dirty = true;
  }






  strongif (list) for (scalar f = *list, *_i118 = list; ((scalar *)&f)->i >= 0; f = *++_i118) {
    scalar phi = _attribute[f.i].phi;
    delete (((scalar []){phi,{-1}}));
    _attribute[f.i].phi.i = 0;
  }
  pfree (list,__func__,__FILE__,__LINE__);
 end_trace("acceleration_0", "/home/fpl/softwares/basilisk/src/iforce.h", 131); } return 0; } 
#line 16 "/home/fpl/softwares/basilisk/src/tension.h"








#line 36 "/home/fpl/softwares/basilisk/src/tension.h"
static int stability_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int stability_0 (const int i, const double t, Event * _ev) { trace ("stability_0", "/home/fpl/softwares/basilisk/src/tension.h", 36); 
{





  double amin = HUGE, amax = -HUGE, dmin = HUGE;
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _amin = amin;
 double _amax = amax;
 double _dmin = dmin;
{ double amin = _amin; NOT_UNUSED(amin);
 double amax = _amax; NOT_UNUSED(amax);
 double dmin = _dmin; NOT_UNUSED(dmin);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/tension.h", .line = 44,
    .each = "foreach_face", .first = _first_call
  };

strongif (!is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 44
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 44
{

#line 44 "/home/fpl/softwares/basilisk/src/tension.h"

    IF (val_fm_x(fm.x,0,0,0) > 0.) {
      IF (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) > amax) amax = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
      IF (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) < amin) amin = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
      IF (Delta < dmin) dmin = Delta;
    } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 44
{

#line 44 "/home/fpl/softwares/basilisk/src/tension.h"

    IF (val_fm_y(fm.y,0,0,0) > 0.) {
      IF (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) > amax) amax = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
      IF (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) < amin) amin = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
      IF (Delta < dmin) dmin = Delta;
    } }  }}  end_foreach_face_stencil()
#line 49
 }
strongif (is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 44
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 44
{

#line 44 "/home/fpl/softwares/basilisk/src/tension.h"

    IF (val_fm_x(fm.x,0,0,0) > 0.) {
      IF (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) > amax) amax = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
      IF (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) < amin) amin = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
      IF (Delta < dmin) dmin = Delta;
    } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 44
{

#line 44 "/home/fpl/softwares/basilisk/src/tension.h"

    IF (val_fm_y(fm.y,0,0,0) > 0.) {
      IF (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) > amax) amax = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
      IF (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) < amin) amin = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
      IF (Delta < dmin) dmin = Delta;
    } }  }}  end_foreach_face_stencil()
#line 49
 }
strongif (!is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 44
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 44
{

#line 44 "/home/fpl/softwares/basilisk/src/tension.h"

    IF (val_fm_x(fm.x,0,0,0) > 0.) {
      IF (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) > amax) amax = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
      IF (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) < amin) amin = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
      IF (Delta < dmin) dmin = Delta;
    } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 44
{

#line 44 "/home/fpl/softwares/basilisk/src/tension.h"

    IF (val_fm_y(fm.y,0,0,0) > 0.) {
      IF (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) > amax) amax = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
      IF (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) < amin) amin = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
      IF (Delta < dmin) dmin = Delta;
    } }  }}  end_foreach_face_stencil()
#line 49
 }
strongif (is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 44
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 44
{

#line 44 "/home/fpl/softwares/basilisk/src/tension.h"

    IF (val_fm_x(fm.x,0,0,0) > 0.) {
      IF (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) > amax) amax = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
      IF (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) < amin) amin = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
      IF (Delta < dmin) dmin = Delta;
    } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 44
{

#line 44 "/home/fpl/softwares/basilisk/src/tension.h"

    IF (val_fm_y(fm.y,0,0,0) > 0.) {
      IF (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) > amax) amax = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
      IF (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) < amin) amin = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
      IF (Delta < dmin) dmin = Delta;
    } }  }}  end_foreach_face_stencil()
#line 49
 }  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 49

#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(min:amin)  reduction(max:amax)  reduction(min:dmin)) {

#line 44

strongif (!is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 44
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 44
{

#line 44 "/home/fpl/softwares/basilisk/src/tension.h"

    if (val_fm_x(fm.x,0,0,0) > 0.) {
      if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) > amax) amax = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
      if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) < amin) amin = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
      if (Delta < dmin) dmin = Delta;
    } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 44
{

#line 44 "/home/fpl/softwares/basilisk/src/tension.h"

    if (val_fm_y(fm.y,0,0,0) > 0.) {
      if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) > amax) amax = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
      if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) < amin) amin = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
      if (Delta < dmin) dmin = Delta;
    } }  }}  end_foreach_face_generic()
#line 49
 end_foreach_face(); }
strongif (is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 44
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 44
{

#line 44 "/home/fpl/softwares/basilisk/src/tension.h"

    if (val_fm_x(fm.x,0,0,0) > 0.) {
      if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) > amax) amax = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
      if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) < amin) amin = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
      if (Delta < dmin) dmin = Delta;
    } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 44
{

#line 44 "/home/fpl/softwares/basilisk/src/tension.h"

    if (val_fm_y(fm.y,0,0,0) > 0.) {
      if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) > amax) amax = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
      if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) < amin) amin = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
      if (Delta < dmin) dmin = Delta;
    } }  }}  end_foreach_face_generic()
#line 49
 end_foreach_face(); }
strongif (!is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 44
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 44
{

#line 44 "/home/fpl/softwares/basilisk/src/tension.h"

    if (val_fm_x(fm.x,0,0,0) > 0.) {
      if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) > amax) amax = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
      if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) < amin) amin = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
      if (Delta < dmin) dmin = Delta;
    } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 44
{

#line 44 "/home/fpl/softwares/basilisk/src/tension.h"

    if (val_fm_y(fm.y,0,0,0) > 0.) {
      if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) > amax) amax = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
      if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) < amin) amin = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
      if (Delta < dmin) dmin = Delta;
    } }  }}  end_foreach_face_generic()
#line 49
 end_foreach_face(); }
strongif (is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 44
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 44
{

#line 44 "/home/fpl/softwares/basilisk/src/tension.h"

    if (val_fm_x(fm.x,0,0,0) > 0.) {
      if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) > amax) amax = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
      if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) < amin) amin = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
      if (Delta < dmin) dmin = Delta;
    } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 44
{

#line 44 "/home/fpl/softwares/basilisk/src/tension.h"

    if (val_fm_y(fm.y,0,0,0) > 0.) {
      if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) > amax) amax = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
      if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) < amin) amin = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
      if (Delta < dmin) dmin = Delta;
    } }  }}  end_foreach_face_generic()
#line 49
 end_foreach_face(); }mpi_all_reduce_array (&amin, double, MPI_MIN, 1);
mpi_all_reduce_array (&amax, double, MPI_MAX, 1);
mpi_all_reduce_array (&dmin, double, MPI_MIN, 1);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 49
 }
  double rhom = (1./amin + 1./amax)/2.;





  double sigma = 0.;
  strongif (interfaces) for (scalar c = *interfaces, *_i119 = interfaces; ((scalar *)&c)->i >= 0; c = *++_i119)
    sigma += _attribute[c.i].sigma;
  if (sigma) {
    double dt = sqrt (rhom*cube(dmin)/(pi*sigma));
    if (dt < dtmax)
      dtmax = dt;
  }
 end_trace("stability_0", "/home/fpl/softwares/basilisk/src/tension.h", 64); } return 0; } 







static int acceleration_1_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int acceleration_1 (const int i, const double t, Event * _ev) { trace ("acceleration_1", "/home/fpl/softwares/basilisk/src/tension.h", 72); 
{




  strongif (interfaces) for (scalar f = *interfaces, *_i120 = interfaces; ((scalar *)&f)->i >= 0; f = *++_i120)
    if (_attribute[f.i].sigma) {





      scalar phi = _attribute[f.i].phi;
      if (phi.i)
 curvature ((struct Curvature){f, phi, _attribute[f.i].sigma, .add = true});
      else {
 phi = new_scalar("phi");
 curvature ((struct Curvature){f, phi, _attribute[f.i].sigma, .add = false});
 _attribute[f.i].phi = phi;
      }
    }
 end_trace("acceleration_1", "/home/fpl/softwares/basilisk/src/tension.h", 94); } return 0; } 
#line 23 "60d_plate_advancing_simulation.c"
#line 1 "two-phase.h"
#line 1 "/home/fpl/softwares/basilisk/src/two-phase.h"
#line 13 "/home/fpl/softwares/basilisk/src/two-phase.h"
#line 1 "vof.h"
#line 1 "/home/fpl/softwares/basilisk/src/vof.h"
#line 27 "/home/fpl/softwares/basilisk/src/vof.h"




#line 44 "/home/fpl/softwares/basilisk/src/vof.h"
extern scalar * interfaces;
extern vector uf;
extern double dt;








#line 54

static double vof_concentration_gradient_x (Point point, scalar c, scalar t)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 56 "/home/fpl/softwares/basilisk/src/vof.h"

  static const double cmin = 0.5;
  double cl = val(c,-1,0,0), cc = val(c,0,0,0), cr = val(c,1,0,0);
  if (_attribute[t.i].inverse)
    cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
  if (cc >= cmin && _attribute[t.i].gradient != zero) {
    if (cr >= cmin) {
      if (cl >= cmin)
 return _attribute[t.i].gradient (val(t,-1,0,0)/cl, val(t,0,0,0)/cc, val(t,1,0,0)/cr)/Delta;
      else
 return (val(t,1,0,0)/cr - val(t,0,0,0)/cc)/Delta;
    }
    else if (cl >= cmin)
      return (val(t,0,0,0)/cc - val(t,-1,0,0)/cl)/Delta;
  }
  return 0.;

#if _call_vof_concentration_gradient_x
}
#define _IN_STENCIL 1

#line 55
static double _vof_concentration_gradient_x (Point point, scalar c, scalar t)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 56 "/home/fpl/softwares/basilisk/src/vof.h"

  static const double cmin = 0.5;
  double cl = _stencil_val(__FILE__,__LINE__,c,-1,0,0), cc = _stencil_val(__FILE__,__LINE__,c,0,0,0), cr = _stencil_val(__FILE__,__LINE__,c,1,0,0);
  IF (_attribute[t.i].inverse)
    cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
  IF (cc >= cmin && _attribute[t.i].gradient != zero) {
    IF (cr >= cmin) {
      IF (cl >= cmin)
 return _attribute[t.i].gradient (_stencil_val(__FILE__,__LINE__,t,-1,0,0)/cl, _stencil_val(__FILE__,__LINE__,t,0,0,0)/cc, _stencil_val(__FILE__,__LINE__,t,1,0,0)/cr)/Delta;
      
 return (_stencil_val(__FILE__,__LINE__,t,1,0,0)/cr - _stencil_val(__FILE__,__LINE__,t,0,0,0)/cc)/Delta;
    }
    IF (cl >= cmin)
      return (_stencil_val(__FILE__,__LINE__,t,0,0,0)/cc - _stencil_val(__FILE__,__LINE__,t,-1,0,0)/cl)/Delta;
  }
  return 0.;

#undef _IN_STENCIL

#endif

#line 72
}
#line 54

static double vof_concentration_gradient_y (Point point, scalar c, scalar t)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 56 "/home/fpl/softwares/basilisk/src/vof.h"

  static const double cmin = 0.5;
  double cl = val(c,0,-1,0), cc = val(c,0,0,0), cr = val(c,0,1,0);
  if (_attribute[t.i].inverse)
    cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
  if (cc >= cmin && _attribute[t.i].gradient != zero) {
    if (cr >= cmin) {
      if (cl >= cmin)
 return _attribute[t.i].gradient (val(t,0,-1,0)/cl, val(t,0,0,0)/cc, val(t,0,1,0)/cr)/Delta;
      else
 return (val(t,0,1,0)/cr - val(t,0,0,0)/cc)/Delta;
    }
    else if (cl >= cmin)
      return (val(t,0,0,0)/cc - val(t,0,-1,0)/cl)/Delta;
  }
  return 0.;

#if _call_vof_concentration_gradient_y
}
#define _IN_STENCIL 1

#line 55
static double _vof_concentration_gradient_y (Point point, scalar c, scalar t)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 56 "/home/fpl/softwares/basilisk/src/vof.h"

  static const double cmin = 0.5;
  double cl = _stencil_val(__FILE__,__LINE__,c,-1,0,0), cc = _stencil_val(__FILE__,__LINE__,c,0,0,0), cr = _stencil_val(__FILE__,__LINE__,c,1,0,0);
  IF (_attribute[t.i].inverse)
    cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
  IF (cc >= cmin && _attribute[t.i].gradient != zero) {
    IF (cr >= cmin) {
      IF (cl >= cmin)
 return _attribute[t.i].gradient (_stencil_val(__FILE__,__LINE__,t,-1,0,0)/cl, _stencil_val(__FILE__,__LINE__,t,0,0,0)/cc, _stencil_val(__FILE__,__LINE__,t,1,0,0)/cr)/Delta;
      
 return (_stencil_val(__FILE__,__LINE__,t,1,0,0)/cr - _stencil_val(__FILE__,__LINE__,t,0,0,0)/cc)/Delta;
    }
    IF (cl >= cmin)
      return (_stencil_val(__FILE__,__LINE__,t,0,0,0)/cc - _stencil_val(__FILE__,__LINE__,t,-1,0,0)/cl)/Delta;
  }
  return 0.;

#undef _IN_STENCIL

#endif

#line 72
}






static void vof_concentration_refine (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 80 "/home/fpl/softwares/basilisk/src/vof.h"

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 80

  scalar f = _attribute[s.i].c;
  if (val_cm(cm,0,0,0) == 0. || (!_attribute[s.i].inverse && val(f,0,0,0) <= 0.) || (_attribute[s.i].inverse && val(f,0,0,0) >= 1.))
     { foreach_child()
      val(s,0,0,0) = 0.; end_foreach_child(); }
  else {
    coord g;
    {
#line 87

      g.x = Delta*vof_concentration_gradient_x (point, f, s);
#line 87

      g.y = Delta*vof_concentration_gradient_y (point, f, s);}
    double sc = _attribute[s.i].inverse ? val(s,0,0,0)/(1. - val(f,0,0,0)) : val(s,0,0,0)/val(f,0,0,0), cmc = 4.*val_cm(cm,0,0,0);
     { foreach_child() {
      val(s,0,0,0) = sc;
      {
#line 92

 val(s,0,0,0) += child.x*g.x*val_cm(cm,-child.x,0,0)/cmc;
#line 92

 val(s,0,0,0) += child.y*g.y*val_cm(cm,0,-child.y,0)/cmc;}
      val(s,0,0,0) *= _attribute[s.i].inverse ? 1. - val(f,0,0,0) : val(f,0,0,0);
    } end_foreach_child(); }
  }
 }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 80

  scalar f = _attribute[s.i].c;
  if (val_cm(cm,0,0,0) == 0. || (!_attribute[s.i].inverse && val(f,0,0,0) <= 0.) || (_attribute[s.i].inverse && val(f,0,0,0) >= 1.))
     { foreach_child()
      val(s,0,0,0) = 0.; end_foreach_child(); }
  else {
    coord g;
    {
#line 87

      g.x = Delta*vof_concentration_gradient_x (point, f, s);
#line 87

      g.y = Delta*vof_concentration_gradient_y (point, f, s);}
    double sc = _attribute[s.i].inverse ? val(s,0,0,0)/(1. - val(f,0,0,0)) : val(s,0,0,0)/val(f,0,0,0), cmc = 4.*val_cm(cm,0,0,0);
     { foreach_child() {
      val(s,0,0,0) = sc;
      {
#line 92

 val(s,0,0,0) += child.x*g.x*val_cm(cm,-child.x,0,0)/cmc;
#line 92

 val(s,0,0,0) += child.y*g.y*val_cm(cm,0,-child.y,0)/cmc;}
      val(s,0,0,0) *= _attribute[s.i].inverse ? 1. - val(f,0,0,0) : val(f,0,0,0);
    } end_foreach_child(); }
  }
 }
#if _call_vof_concentration_refine
}
#define _IN_STENCIL 1
#define vof_concentration_gradient_x _vof_concentration_gradient_x
#define vof_concentration_gradient_y _vof_concentration_gradient_y

#line 79
static void _vof_concentration_refine (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 80 "/home/fpl/softwares/basilisk/src/vof.h"

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 80

  scalar f = _attribute[s.i].c;
  IF (val_cm(cm,0,0,0) == 0. || (!_attribute[s.i].inverse && _stencil_val(__FILE__,__LINE__,f,0,0,0) <= 0.) || (_attribute[s.i].inverse && _stencil_val(__FILE__,__LINE__,f,0,0,0) >= 1.))
     { foreach_child()
      _stencil_val(__FILE__,__LINE__,s,0,0,0) = 0.; end_foreach_child(); }
   {
    coord g;
    {
#line 87

      g.x = Delta*vof_concentration_gradient_x (point, f, s);
#line 87

      g.y = Delta*vof_concentration_gradient_y (point, f, s);}
    double sc = _attribute[s.i].inverse ? _stencil_val(__FILE__,__LINE__,s,0,0,0)/(1. - _stencil_val(__FILE__,__LINE__,f,0,0,0)) : _stencil_val(__FILE__,__LINE__,s,0,0,0)/_stencil_val(__FILE__,__LINE__,f,0,0,0), cmc = 4.*val_cm(cm,0,0,0);
     { foreach_child() {
      _stencil_val(__FILE__,__LINE__,s,0,0,0) = sc;
      {
#line 92

 _stencil_val(__FILE__,__LINE__,s,0,0,0) += child.x*g.x*val_cm(cm,-child.x,0,0)/cmc;
#line 92

 _stencil_val(__FILE__,__LINE__,s,0,0,0) += child.y*g.y*val_cm(cm,0,-child.y,0)/cmc;}
      _stencil_val(__FILE__,__LINE__,s,0,0,0) *= _attribute[s.i].inverse ? 1. - _stencil_val(__FILE__,__LINE__,f,0,0,0) : _stencil_val(__FILE__,__LINE__,f,0,0,0);
    } end_foreach_child(); }
  }
 }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 80

  scalar f = _attribute[s.i].c;
  IF (val_cm(cm,0,0,0) == 0. || (!_attribute[s.i].inverse && _stencil_val(__FILE__,__LINE__,f,0,0,0) <= 0.) || (_attribute[s.i].inverse && _stencil_val(__FILE__,__LINE__,f,0,0,0) >= 1.))
     { foreach_child()
      _stencil_val(__FILE__,__LINE__,s,0,0,0) = 0.; end_foreach_child(); }
   {
    coord g;
    {
#line 87

      g.x = Delta*vof_concentration_gradient_x (point, f, s);
#line 87

      g.y = Delta*vof_concentration_gradient_y (point, f, s);}
    double sc = _attribute[s.i].inverse ? _stencil_val(__FILE__,__LINE__,s,0,0,0)/(1. - _stencil_val(__FILE__,__LINE__,f,0,0,0)) : _stencil_val(__FILE__,__LINE__,s,0,0,0)/_stencil_val(__FILE__,__LINE__,f,0,0,0), cmc = 4.*val_cm(cm,0,0,0);
     { foreach_child() {
      _stencil_val(__FILE__,__LINE__,s,0,0,0) = sc;
      {
#line 92

 _stencil_val(__FILE__,__LINE__,s,0,0,0) += child.x*g.x*val_cm(cm,-child.x,0,0)/cmc;
#line 92

 _stencil_val(__FILE__,__LINE__,s,0,0,0) += child.y*g.y*val_cm(cm,0,-child.y,0)/cmc;}
      _stencil_val(__FILE__,__LINE__,s,0,0,0) *= _attribute[s.i].inverse ? 1. - _stencil_val(__FILE__,__LINE__,f,0,0,0) : _stencil_val(__FILE__,__LINE__,f,0,0,0);
    } end_foreach_child(); }
  }
 }
#undef vof_concentration_gradient_x
#undef vof_concentration_gradient_y
#undef _IN_STENCIL

#endif

#line 97
}





static int defaults_2_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults_2 (const int i, const double t, Event * _ev) { trace ("defaults_2", "/home/fpl/softwares/basilisk/src/vof.h", 103); 
{
  strongif (interfaces) for (scalar c = *interfaces, *_i121 = interfaces; ((scalar *)&c)->i >= 0; c = *++_i121) {
    _attribute[c.i].refine = _attribute[c.i].prolongation = fraction_refine;
    _attribute[c.i].dirty = true;
    scalar * tracers = _attribute[c.i].tracers;
    strongif (tracers) for (scalar t = *tracers, *_i122 = tracers; ((scalar *)&t)->i >= 0; t = *++_i122) {
      _attribute[t.i].restriction = restriction_volume_average;
      _attribute[t.i].refine = _attribute[t.i].prolongation = vof_concentration_refine;
      _attribute[t.i].dirty = true;
      _attribute[t.i].c = c;
    }
  }
 end_trace("defaults_2", "/home/fpl/softwares/basilisk/src/vof.h", 116); } return 0; } 






static int defaults_3_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults_3 (const int i, const double t, Event * _ev) { trace ("defaults_3", "/home/fpl/softwares/basilisk/src/vof.h", 123); 
{
  strongif (interfaces) for (scalar c = *interfaces, *_i123 = interfaces; ((scalar *)&c)->i >= 0; c = *++_i123) {
    scalar * tracers = _attribute[c.i].tracers;
    strongif (tracers) for (scalar t = *tracers, *_i124 = tracers; ((scalar *)&t)->i >= 0; t = *++_i124)
      _attribute[t.i].depends = list_add (_attribute[t.i].depends, c);
  }
 end_trace("defaults_3", "/home/fpl/softwares/basilisk/src/vof.h", 130); } return 0; } 





static int stability_1_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int stability_1 (const int i, const double t, Event * _ev) { trace ("stability_1", "/home/fpl/softwares/basilisk/src/vof.h", 136);  {
  if (CFL > 0.5)
    CFL = 0.5;
 end_trace("stability_1", "/home/fpl/softwares/basilisk/src/vof.h", 139); } return 0; } 
#line 153 "/home/fpl/softwares/basilisk/src/vof.h"

#line 153

static void sweep_x (scalar c, scalar cc, scalar * tcl)
{
  vector n= new_vector("n");
  scalar alpha= new_scalar("alpha"), flux= new_scalar("flux");
  double cfl = 0.;
#line 167 "/home/fpl/softwares/basilisk/src/vof.h"
  scalar * tracers = _attribute[c.i].tracers, * gfl = NULL, * tfluxl = NULL;
  if (tracers) {
    strongif (tracers) for (scalar t = *tracers, *_i125 = tracers; ((scalar *)&t)->i >= 0; t = *++_i125) {
      scalar gf = new_scalar("gf"), flux = new_scalar("flux");
      gfl = list_append (gfl, gf);
      tfluxl = list_append (tfluxl, flux);
    }




     { 
#define vof_concentration_gradient_x _vof_concentration_gradient_x
#define vof_concentration_gradient_y _vof_concentration_gradient_y
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/vof.h", .line = 178,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 178 "/home/fpl/softwares/basilisk/src/vof.h"
 {
      scalar t, gf;
      scalar * _i6 = tracers; scalar * _i7 = gfl; strongif (tracers) for (t = *tracers, gf = *gfl; ((scalar *)&t)->i >= 0; t = *++_i6, gf = *++_i7)
 _stencil_val(__FILE__,__LINE__,gf,0,0,0) = vof_concentration_gradient_x (point, c, t);
    } } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#undef vof_concentration_gradient_x
#undef vof_concentration_gradient_y
#line 182
foreach(){

#line 178 "/home/fpl/softwares/basilisk/src/vof.h"
 {
      scalar t, gf;
      scalar * _i6 = tracers; scalar * _i7 = gfl; strongif (tracers) for (t = *tracers, gf = *gfl; ((scalar *)&t)->i >= 0; t = *++_i6, gf = *++_i7)
 val(gf,0,0,0) = vof_concentration_gradient_x (point, c, t);
    } } end_foreach(); }
  }






  reconstruction (c, n, alpha);
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dt = dt;
 double _cfl = cfl;
{ double dt = _dt; NOT_UNUSED(dt);
 double cfl = _cfl; NOT_UNUSED(cfl);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/vof.h", .line = 191,
    .each = "foreach_face", .first = _first_call
  };

strongif (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 191
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 191
{

#line 191 "/home/fpl/softwares/basilisk/src/vof.h"
 {






    double un = _stencil_val(__FILE__,__LINE__,uf.x,0,0,0)*dt/(Delta*val_fm_x(fm.x,0,0,0) + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    IF (un*val_fm_x(fm.x,0,0,0)*s/(val_cm(cm,0,0,0) + 0.) > cfl)
      cfl = un*val_fm_x(fm.x,0,0,0)*s/(val_cm(cm,0,0,0) + 0.);
#line 221 "/home/fpl/softwares/basilisk/src/vof.h"
    double cf;
    IF (_stencil_val(__FILE__,__LINE__,c,i,0,0) <= 0. || _stencil_val(__FILE__,__LINE__,c,i,0,0) >= 1.)
      cf = _stencil_val(__FILE__,__LINE__,c,i,0,0);
    
      cf = rectangle_fraction ((coord){-s*_stencil_val(__FILE__,__LINE__,n.x,i,0,0), _stencil_val(__FILE__,__LINE__,n.y,i,0,0), _val_higher_dimension(n.x,i,0,0)}, _stencil_val(__FILE__,__LINE__,alpha,i,0,0),
          (coord){-0.5, -0.5, -0.5},
          (coord){s*un - 0.5, 0.5, 0.5});





    _stencil_val(__FILE__,__LINE__,flux,0,0,0) = cf*_stencil_val(__FILE__,__LINE__,uf.x,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; strongif (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = _stencil_val(__FILE__,__LINE__,c,i,0,0);
      IF (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      IF (ci > 1e-10) {
 double ff = _stencil_val(__FILE__,__LINE__,t,i,0,0)/ci + s*min(1., 1. - s*un)*_stencil_val(__FILE__,__LINE__,gf,i,0,0)*Delta/2.;
 _stencil_val(__FILE__,__LINE__,tflux,0,0,0) = ff*cf1*_stencil_val(__FILE__,__LINE__,uf.x,0,0,0);
      }
      
 _stencil_val(__FILE__,__LINE__,tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_stencil()
#line 252
 }
strongif (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 191
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 191
{

#line 191 "/home/fpl/softwares/basilisk/src/vof.h"
 {






    double un = _stencil_val(__FILE__,__LINE__,uf.x,0,0,0)*dt/(Delta*val_fm_x(fm.x,0,0,0) + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    IF (un*val_fm_x(fm.x,0,0,0)*s/(val_cm(cm,0,0,0) + 0.) > cfl)
      cfl = un*val_fm_x(fm.x,0,0,0)*s/(val_cm(cm,0,0,0) + 0.);
#line 221 "/home/fpl/softwares/basilisk/src/vof.h"
    double cf;
    IF (_stencil_val(__FILE__,__LINE__,c,i,0,0) <= 0. || _stencil_val(__FILE__,__LINE__,c,i,0,0) >= 1.)
      cf = _stencil_val(__FILE__,__LINE__,c,i,0,0);
    
      cf = rectangle_fraction ((coord){-s*_stencil_val(__FILE__,__LINE__,n.x,i,0,0), _stencil_val(__FILE__,__LINE__,n.y,i,0,0), _val_higher_dimension(n.x,i,0,0)}, _stencil_val(__FILE__,__LINE__,alpha,i,0,0),
          (coord){-0.5, -0.5, -0.5},
          (coord){s*un - 0.5, 0.5, 0.5});





    _stencil_val(__FILE__,__LINE__,flux,0,0,0) = cf*_stencil_val(__FILE__,__LINE__,uf.x,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; strongif (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = _stencil_val(__FILE__,__LINE__,c,i,0,0);
      IF (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      IF (ci > 1e-10) {
 double ff = _stencil_val(__FILE__,__LINE__,t,i,0,0)/ci + s*min(1., 1. - s*un)*_stencil_val(__FILE__,__LINE__,gf,i,0,0)*Delta/2.;
 _stencil_val(__FILE__,__LINE__,tflux,0,0,0) = ff*cf1*_stencil_val(__FILE__,__LINE__,uf.x,0,0,0);
      }
      
 _stencil_val(__FILE__,__LINE__,tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_stencil()
#line 252
 }
strongif (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 191
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 191
{

#line 191 "/home/fpl/softwares/basilisk/src/vof.h"
 {






    double un = _stencil_val(__FILE__,__LINE__,uf.x,0,0,0)*dt/(Delta*val_fm_x(fm.x,0,0,0) + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    IF (un*val_fm_x(fm.x,0,0,0)*s/(val_cm(cm,0,0,0) + 0.) > cfl)
      cfl = un*val_fm_x(fm.x,0,0,0)*s/(val_cm(cm,0,0,0) + 0.);
#line 221 "/home/fpl/softwares/basilisk/src/vof.h"
    double cf;
    IF (_stencil_val(__FILE__,__LINE__,c,i,0,0) <= 0. || _stencil_val(__FILE__,__LINE__,c,i,0,0) >= 1.)
      cf = _stencil_val(__FILE__,__LINE__,c,i,0,0);
    
      cf = rectangle_fraction ((coord){-s*_stencil_val(__FILE__,__LINE__,n.x,i,0,0), _stencil_val(__FILE__,__LINE__,n.y,i,0,0), _val_higher_dimension(n.x,i,0,0)}, _stencil_val(__FILE__,__LINE__,alpha,i,0,0),
          (coord){-0.5, -0.5, -0.5},
          (coord){s*un - 0.5, 0.5, 0.5});





    _stencil_val(__FILE__,__LINE__,flux,0,0,0) = cf*_stencil_val(__FILE__,__LINE__,uf.x,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; strongif (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = _stencil_val(__FILE__,__LINE__,c,i,0,0);
      IF (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      IF (ci > 1e-10) {
 double ff = _stencil_val(__FILE__,__LINE__,t,i,0,0)/ci + s*min(1., 1. - s*un)*_stencil_val(__FILE__,__LINE__,gf,i,0,0)*Delta/2.;
 _stencil_val(__FILE__,__LINE__,tflux,0,0,0) = ff*cf1*_stencil_val(__FILE__,__LINE__,uf.x,0,0,0);
      }
      
 _stencil_val(__FILE__,__LINE__,tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_stencil()
#line 252
 }
strongif (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 191
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 191
{

#line 191 "/home/fpl/softwares/basilisk/src/vof.h"
 {






    double un = _stencil_val(__FILE__,__LINE__,uf.x,0,0,0)*dt/(Delta*val_fm_x(fm.x,0,0,0) + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    IF (un*val_fm_x(fm.x,0,0,0)*s/(val_cm(cm,0,0,0) + 0.) > cfl)
      cfl = un*val_fm_x(fm.x,0,0,0)*s/(val_cm(cm,0,0,0) + 0.);
#line 221 "/home/fpl/softwares/basilisk/src/vof.h"
    double cf;
    IF (_stencil_val(__FILE__,__LINE__,c,i,0,0) <= 0. || _stencil_val(__FILE__,__LINE__,c,i,0,0) >= 1.)
      cf = _stencil_val(__FILE__,__LINE__,c,i,0,0);
    
      cf = rectangle_fraction ((coord){-s*_stencil_val(__FILE__,__LINE__,n.x,i,0,0), _stencil_val(__FILE__,__LINE__,n.y,i,0,0), _val_higher_dimension(n.x,i,0,0)}, _stencil_val(__FILE__,__LINE__,alpha,i,0,0),
          (coord){-0.5, -0.5, -0.5},
          (coord){s*un - 0.5, 0.5, 0.5});





    _stencil_val(__FILE__,__LINE__,flux,0,0,0) = cf*_stencil_val(__FILE__,__LINE__,uf.x,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; strongif (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = _stencil_val(__FILE__,__LINE__,c,i,0,0);
      IF (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      IF (ci > 1e-10) {
 double ff = _stencil_val(__FILE__,__LINE__,t,i,0,0)/ci + s*min(1., 1. - s*un)*_stencil_val(__FILE__,__LINE__,gf,i,0,0)*Delta/2.;
 _stencil_val(__FILE__,__LINE__,tflux,0,0,0) = ff*cf1*_stencil_val(__FILE__,__LINE__,uf.x,0,0,0);
      }
      
 _stencil_val(__FILE__,__LINE__,tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_stencil()
#line 252
 } if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/fpl/softwares/basilisk/src/vof.h", 191, "dt");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 252

#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel   reduction (max:cfl)) {

#line 191

strongif (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 191
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 191
{

#line 191 "/home/fpl/softwares/basilisk/src/vof.h"
 {






    double un = val(uf.x,0,0,0)*dt/(Delta*val_fm_x(fm.x,0,0,0) + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    if (un*val_fm_x(fm.x,0,0,0)*s/(val_cm(cm,0,0,0) + 0.) > cfl)
      cfl = un*val_fm_x(fm.x,0,0,0)*s/(val_cm(cm,0,0,0) + 0.);
#line 221 "/home/fpl/softwares/basilisk/src/vof.h"
    double cf;
    if (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.)
      cf = val(c,i,0,0);
    else
      cf = rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), _val_higher_dimension(n.x,i,0,0)}, val(alpha,i,0,0),
          (coord){-0.5, -0.5, -0.5},
          (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.x,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; strongif (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = val(c,i,0,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,i,0,0)/ci + s*min(1., 1. - s*un)*val(gf,i,0,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.x,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_generic()
#line 252
 end_foreach_face(); }
strongif (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 191
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 191
{

#line 191 "/home/fpl/softwares/basilisk/src/vof.h"
 {






    double un = val(uf.x,0,0,0)*dt/(Delta*val_fm_x(fm.x,0,0,0) + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    if (un*val_fm_x(fm.x,0,0,0)*s/(val_cm(cm,0,0,0) + 0.) > cfl)
      cfl = un*val_fm_x(fm.x,0,0,0)*s/(val_cm(cm,0,0,0) + 0.);
#line 221 "/home/fpl/softwares/basilisk/src/vof.h"
    double cf;
    if (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.)
      cf = val(c,i,0,0);
    else
      cf = rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), _val_higher_dimension(n.x,i,0,0)}, val(alpha,i,0,0),
          (coord){-0.5, -0.5, -0.5},
          (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.x,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; strongif (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = val(c,i,0,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,i,0,0)/ci + s*min(1., 1. - s*un)*val(gf,i,0,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.x,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_generic()
#line 252
 end_foreach_face(); }
strongif (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 191
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 191
{

#line 191 "/home/fpl/softwares/basilisk/src/vof.h"
 {






    double un = val(uf.x,0,0,0)*dt/(Delta*val_fm_x(fm.x,0,0,0) + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    if (un*val_fm_x(fm.x,0,0,0)*s/(val_cm(cm,0,0,0) + 0.) > cfl)
      cfl = un*val_fm_x(fm.x,0,0,0)*s/(val_cm(cm,0,0,0) + 0.);
#line 221 "/home/fpl/softwares/basilisk/src/vof.h"
    double cf;
    if (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.)
      cf = val(c,i,0,0);
    else
      cf = rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), _val_higher_dimension(n.x,i,0,0)}, val(alpha,i,0,0),
          (coord){-0.5, -0.5, -0.5},
          (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.x,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; strongif (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = val(c,i,0,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,i,0,0)/ci + s*min(1., 1. - s*un)*val(gf,i,0,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.x,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_generic()
#line 252
 end_foreach_face(); }
strongif (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 191
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 191
{

#line 191 "/home/fpl/softwares/basilisk/src/vof.h"
 {






    double un = val(uf.x,0,0,0)*dt/(Delta*val_fm_x(fm.x,0,0,0) + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    if (un*val_fm_x(fm.x,0,0,0)*s/(val_cm(cm,0,0,0) + 0.) > cfl)
      cfl = un*val_fm_x(fm.x,0,0,0)*s/(val_cm(cm,0,0,0) + 0.);
#line 221 "/home/fpl/softwares/basilisk/src/vof.h"
    double cf;
    if (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.)
      cf = val(c,i,0,0);
    else
      cf = rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), _val_higher_dimension(n.x,i,0,0)}, val(alpha,i,0,0),
          (coord){-0.5, -0.5, -0.5},
          (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.x,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; strongif (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = val(c,i,0,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,i,0,0)/ci + s*min(1., 1. - s*un)*val(gf,i,0,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.x,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_generic()
#line 252
 end_foreach_face(); }mpi_all_reduce_array (&cfl, double, MPI_MAX, 1);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 252
 }
  delete (gfl); pfree (gfl,__func__,__FILE__,__LINE__);




  if (cfl > 0.5 + 1e-6)
    fprintf (ferr,
      "WARNING: CFL must be <= 0.5 for VOF (cfl - 0.5 = %g)\n",
      cfl - 0.5), fflush (ferr);
#line 281 "/home/fpl/softwares/basilisk/src/vof.h"
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dt = dt;
{ double dt = _dt; NOT_UNUSED(dt);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/vof.h", .line = 281,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 281
foreach_stencil(){

#line 281 "/home/fpl/softwares/basilisk/src/vof.h"
 {
    _stencil_val(__FILE__,__LINE__,c,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,flux,0,0,0) - _stencil_val(__FILE__,__LINE__,flux,1,0,0) + _stencil_val(__FILE__,__LINE__,cc,0,0,0)*(_stencil_val(__FILE__,__LINE__,uf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,uf.x,0,0,0)))/(val_cm(cm,0,0,0)*Delta);
    scalar t, tc, tflux;
    scalar * _i11 = tracers; scalar * _i12 = tcl; scalar * _i13 = tfluxl; strongif (tracers) for (t = *tracers, tc = *tcl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i11, tc = *++_i12, tflux = *++_i13)
      _stencil_val(__FILE__,__LINE__,t,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,tflux,0,0,0) - _stencil_val(__FILE__,__LINE__,tflux,1,0,0) + _stencil_val(__FILE__,__LINE__,tc,0,0,0)*(_stencil_val(__FILE__,__LINE__,uf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,uf.x,0,0,0)))/(val_cm(cm,0,0,0)*Delta);
  } } end_foreach_stencil(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 281
foreach_stencil(){

#line 281 "/home/fpl/softwares/basilisk/src/vof.h"
 {
    _stencil_val(__FILE__,__LINE__,c,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,flux,0,0,0) - _stencil_val(__FILE__,__LINE__,flux,1,0,0) + _stencil_val(__FILE__,__LINE__,cc,0,0,0)*(_stencil_val(__FILE__,__LINE__,uf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,uf.x,0,0,0)))/(val_cm(cm,0,0,0)*Delta);
    scalar t, tc, tflux;
    scalar * _i11 = tracers; scalar * _i12 = tcl; scalar * _i13 = tfluxl; strongif (tracers) for (t = *tracers, tc = *tcl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i11, tc = *++_i12, tflux = *++_i13)
      _stencil_val(__FILE__,__LINE__,t,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,tflux,0,0,0) - _stencil_val(__FILE__,__LINE__,tflux,1,0,0) + _stencil_val(__FILE__,__LINE__,tc,0,0,0)*(_stencil_val(__FILE__,__LINE__,uf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,uf.x,0,0,0)))/(val_cm(cm,0,0,0)*Delta);
  } } end_foreach_stencil(); } if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/fpl/softwares/basilisk/src/vof.h", 281, "dt");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 286

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 281
foreach(){

#line 281 "/home/fpl/softwares/basilisk/src/vof.h"
 {
    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,1,0,0) + val(cc,0,0,0)*(val(uf.x,1,0,0) - val(uf.x,0,0,0)))/(val_cm(cm,0,0,0)*Delta);
    scalar t, tc, tflux;
    scalar * _i11 = tracers; scalar * _i12 = tcl; scalar * _i13 = tfluxl; strongif (tracers) for (t = *tracers, tc = *tcl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i11, tc = *++_i12, tflux = *++_i13)
      val(t,0,0,0) += dt*(val(tflux,0,0,0) - val(tflux,1,0,0) + val(tc,0,0,0)*(val(uf.x,1,0,0) - val(uf.x,0,0,0)))/(val_cm(cm,0,0,0)*Delta);
  } } end_foreach(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 281
foreach(){

#line 281 "/home/fpl/softwares/basilisk/src/vof.h"
 {
    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,1,0,0) + val(cc,0,0,0)*(val(uf.x,1,0,0) - val(uf.x,0,0,0)))/(val_cm(cm,0,0,0)*Delta);
    scalar t, tc, tflux;
    scalar * _i11 = tracers; scalar * _i12 = tcl; scalar * _i13 = tfluxl; strongif (tracers) for (t = *tracers, tc = *tcl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i11, tc = *++_i12, tflux = *++_i13)
      val(t,0,0,0) += dt*(val(tflux,0,0,0) - val(tflux,1,0,0) + val(tc,0,0,0)*(val(uf.x,1,0,0) - val(uf.x,0,0,0)))/(val_cm(cm,0,0,0)*Delta);
  } } end_foreach(); } }
#line 304 "/home/fpl/softwares/basilisk/src/vof.h"
  delete (tfluxl); pfree (tfluxl,__func__,__FILE__,__LINE__);
 delete (((scalar []){flux,alpha,n.x,n.y,{-1}})); }
#line 153

static void sweep_y (scalar c, scalar cc, scalar * tcl)
{
  vector n= new_vector("n");
  scalar alpha= new_scalar("alpha"), flux= new_scalar("flux");
  double cfl = 0.;
#line 167 "/home/fpl/softwares/basilisk/src/vof.h"
  scalar * tracers = _attribute[c.i].tracers, * gfl = NULL, * tfluxl = NULL;
  if (tracers) {
    strongif (tracers) for (scalar t = *tracers, *_i125 = tracers; ((scalar *)&t)->i >= 0; t = *++_i125) {
      scalar gf = new_scalar("gf"), flux = new_scalar("flux");
      gfl = list_append (gfl, gf);
      tfluxl = list_append (tfluxl, flux);
    }




     { 
#define vof_concentration_gradient_y _vof_concentration_gradient_y
#define vof_concentration_gradient_x _vof_concentration_gradient_x
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/vof.h", .line = 178,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 178 "/home/fpl/softwares/basilisk/src/vof.h"
 {
      scalar t, gf;
      scalar * _i6 = tracers; scalar * _i7 = gfl; strongif (tracers) for (t = *tracers, gf = *gfl; ((scalar *)&t)->i >= 0; t = *++_i6, gf = *++_i7)
 _stencil_val(__FILE__,__LINE__,gf,0,0,0) = vof_concentration_gradient_y (point, c, t);
    } } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#undef vof_concentration_gradient_y
#undef vof_concentration_gradient_x
#line 182
foreach(){

#line 178 "/home/fpl/softwares/basilisk/src/vof.h"
 {
      scalar t, gf;
      scalar * _i6 = tracers; scalar * _i7 = gfl; strongif (tracers) for (t = *tracers, gf = *gfl; ((scalar *)&t)->i >= 0; t = *++_i6, gf = *++_i7)
 val(gf,0,0,0) = vof_concentration_gradient_y (point, c, t);
    } } end_foreach(); }
  }






  reconstruction (c, n, alpha);
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dt = dt;
 double _cfl = cfl;
{ double dt = _dt; NOT_UNUSED(dt);
 double cfl = _cfl; NOT_UNUSED(cfl);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/vof.h", .line = 191,
    .each = "foreach_face", .first = _first_call
  };

strongif (!is_constant(fm.y) && !is_constant(cm)) {
#undef val_fm_y
#define val_fm_y(a,j,i,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,j,i,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_cm
#define val_cm(a,j,i,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 191
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 191
{

#line 191 "/home/fpl/softwares/basilisk/src/vof.h"
 {






    double un = _stencil_val(__FILE__,__LINE__,uf.y,0,0,0)*dt/(Delta*val_fm_y(fm.y,0,0,0) + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    IF (un*val_fm_y(fm.y,0,0,0)*s/(val_cm(cm,0,0,0) + 0.) > cfl)
      cfl = un*val_fm_y(fm.y,0,0,0)*s/(val_cm(cm,0,0,0) + 0.);
#line 221 "/home/fpl/softwares/basilisk/src/vof.h"
    double cf;
    IF (_stencil_val(__FILE__,__LINE__,c,i,0,0) <= 0. || _stencil_val(__FILE__,__LINE__,c,i,0,0) >= 1.)
      cf = _stencil_val(__FILE__,__LINE__,c,i,0,0);
    
      cf = rectangle_fraction ((coord){-s*_stencil_val(__FILE__,__LINE__,n.y,i,0,0), _stencil_val(__FILE__,__LINE__,n.x,i,0,0), _val_higher_dimension(n.y,i,0,0)}, _stencil_val(__FILE__,__LINE__,alpha,i,0,0),
          (coord){-0.5, -0.5, -0.5},
          (coord){s*un - 0.5, 0.5, 0.5});





    _stencil_val(__FILE__,__LINE__,flux,0,0,0) = cf*_stencil_val(__FILE__,__LINE__,uf.y,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; strongif (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = _stencil_val(__FILE__,__LINE__,c,i,0,0);
      IF (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      IF (ci > 1e-10) {
 double ff = _stencil_val(__FILE__,__LINE__,t,i,0,0)/ci + s*min(1., 1. - s*un)*_stencil_val(__FILE__,__LINE__,gf,i,0,0)*Delta/2.;
 _stencil_val(__FILE__,__LINE__,tflux,0,0,0) = ff*cf1*_stencil_val(__FILE__,__LINE__,uf.y,0,0,0);
      }
      
 _stencil_val(__FILE__,__LINE__,tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_stencil()
#line 252
 }
strongif (is_constant(fm.y) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.y.i -_NVARMAX], _constant[fm.x.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_y
#define val_fm_y(a,j,i,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_x
#define val_fm_x(a,j,i,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_cm
#define val_cm(a,j,i,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 191
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 191
{

#line 191 "/home/fpl/softwares/basilisk/src/vof.h"
 {






    double un = _stencil_val(__FILE__,__LINE__,uf.y,0,0,0)*dt/(Delta*val_fm_y(fm.y,0,0,0) + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    IF (un*val_fm_y(fm.y,0,0,0)*s/(val_cm(cm,0,0,0) + 0.) > cfl)
      cfl = un*val_fm_y(fm.y,0,0,0)*s/(val_cm(cm,0,0,0) + 0.);
#line 221 "/home/fpl/softwares/basilisk/src/vof.h"
    double cf;
    IF (_stencil_val(__FILE__,__LINE__,c,i,0,0) <= 0. || _stencil_val(__FILE__,__LINE__,c,i,0,0) >= 1.)
      cf = _stencil_val(__FILE__,__LINE__,c,i,0,0);
    
      cf = rectangle_fraction ((coord){-s*_stencil_val(__FILE__,__LINE__,n.y,i,0,0), _stencil_val(__FILE__,__LINE__,n.x,i,0,0), _val_higher_dimension(n.y,i,0,0)}, _stencil_val(__FILE__,__LINE__,alpha,i,0,0),
          (coord){-0.5, -0.5, -0.5},
          (coord){s*un - 0.5, 0.5, 0.5});





    _stencil_val(__FILE__,__LINE__,flux,0,0,0) = cf*_stencil_val(__FILE__,__LINE__,uf.y,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; strongif (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = _stencil_val(__FILE__,__LINE__,c,i,0,0);
      IF (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      IF (ci > 1e-10) {
 double ff = _stencil_val(__FILE__,__LINE__,t,i,0,0)/ci + s*min(1., 1. - s*un)*_stencil_val(__FILE__,__LINE__,gf,i,0,0)*Delta/2.;
 _stencil_val(__FILE__,__LINE__,tflux,0,0,0) = ff*cf1*_stencil_val(__FILE__,__LINE__,uf.y,0,0,0);
      }
      
 _stencil_val(__FILE__,__LINE__,tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_stencil()
#line 252
 }
strongif (!is_constant(fm.y) && is_constant(cm)) {
#undef val_fm_y
#define val_fm_y(a,j,i,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,j,i,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,j,i,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 191
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 191
{

#line 191 "/home/fpl/softwares/basilisk/src/vof.h"
 {






    double un = _stencil_val(__FILE__,__LINE__,uf.y,0,0,0)*dt/(Delta*val_fm_y(fm.y,0,0,0) + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    IF (un*val_fm_y(fm.y,0,0,0)*s/(val_cm(cm,0,0,0) + 0.) > cfl)
      cfl = un*val_fm_y(fm.y,0,0,0)*s/(val_cm(cm,0,0,0) + 0.);
#line 221 "/home/fpl/softwares/basilisk/src/vof.h"
    double cf;
    IF (_stencil_val(__FILE__,__LINE__,c,i,0,0) <= 0. || _stencil_val(__FILE__,__LINE__,c,i,0,0) >= 1.)
      cf = _stencil_val(__FILE__,__LINE__,c,i,0,0);
    
      cf = rectangle_fraction ((coord){-s*_stencil_val(__FILE__,__LINE__,n.y,i,0,0), _stencil_val(__FILE__,__LINE__,n.x,i,0,0), _val_higher_dimension(n.y,i,0,0)}, _stencil_val(__FILE__,__LINE__,alpha,i,0,0),
          (coord){-0.5, -0.5, -0.5},
          (coord){s*un - 0.5, 0.5, 0.5});





    _stencil_val(__FILE__,__LINE__,flux,0,0,0) = cf*_stencil_val(__FILE__,__LINE__,uf.y,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; strongif (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = _stencil_val(__FILE__,__LINE__,c,i,0,0);
      IF (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      IF (ci > 1e-10) {
 double ff = _stencil_val(__FILE__,__LINE__,t,i,0,0)/ci + s*min(1., 1. - s*un)*_stencil_val(__FILE__,__LINE__,gf,i,0,0)*Delta/2.;
 _stencil_val(__FILE__,__LINE__,tflux,0,0,0) = ff*cf1*_stencil_val(__FILE__,__LINE__,uf.y,0,0,0);
      }
      
 _stencil_val(__FILE__,__LINE__,tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_stencil()
#line 252
 }
strongif (is_constant(fm.y) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.y.i -_NVARMAX], _constant[fm.x.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_y
#define val_fm_y(a,j,i,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_x
#define val_fm_x(a,j,i,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,j,i,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 191
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 191
{

#line 191 "/home/fpl/softwares/basilisk/src/vof.h"
 {






    double un = _stencil_val(__FILE__,__LINE__,uf.y,0,0,0)*dt/(Delta*val_fm_y(fm.y,0,0,0) + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    IF (un*val_fm_y(fm.y,0,0,0)*s/(val_cm(cm,0,0,0) + 0.) > cfl)
      cfl = un*val_fm_y(fm.y,0,0,0)*s/(val_cm(cm,0,0,0) + 0.);
#line 221 "/home/fpl/softwares/basilisk/src/vof.h"
    double cf;
    IF (_stencil_val(__FILE__,__LINE__,c,i,0,0) <= 0. || _stencil_val(__FILE__,__LINE__,c,i,0,0) >= 1.)
      cf = _stencil_val(__FILE__,__LINE__,c,i,0,0);
    
      cf = rectangle_fraction ((coord){-s*_stencil_val(__FILE__,__LINE__,n.y,i,0,0), _stencil_val(__FILE__,__LINE__,n.x,i,0,0), _val_higher_dimension(n.y,i,0,0)}, _stencil_val(__FILE__,__LINE__,alpha,i,0,0),
          (coord){-0.5, -0.5, -0.5},
          (coord){s*un - 0.5, 0.5, 0.5});





    _stencil_val(__FILE__,__LINE__,flux,0,0,0) = cf*_stencil_val(__FILE__,__LINE__,uf.y,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; strongif (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = _stencil_val(__FILE__,__LINE__,c,i,0,0);
      IF (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      IF (ci > 1e-10) {
 double ff = _stencil_val(__FILE__,__LINE__,t,i,0,0)/ci + s*min(1., 1. - s*un)*_stencil_val(__FILE__,__LINE__,gf,i,0,0)*Delta/2.;
 _stencil_val(__FILE__,__LINE__,tflux,0,0,0) = ff*cf1*_stencil_val(__FILE__,__LINE__,uf.y,0,0,0);
      }
      
 _stencil_val(__FILE__,__LINE__,tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_stencil()
#line 252
 } if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/fpl/softwares/basilisk/src/vof.h", 191, "dt");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 252

#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel   reduction (max:cfl)) {

#line 191

strongif (!is_constant(fm.y) && !is_constant(cm)) {
#undef val_fm_y
#define val_fm_y(a,j,i,k) val(a,j,i,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,j,i,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,j,i,k)
#undef val_fm_x
#define val_fm_x(a,j,i,k) val(a,j,i,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,j,i,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,j,i,k)
#undef val_cm
#define val_cm(a,j,i,k) val(a,j,i,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,j,i,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,j,i,k)
#line 191
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_y()) {
#line 191
{

#line 191 "/home/fpl/softwares/basilisk/src/vof.h"
 {






    double un = val(uf.y,0,0,0)*dt/(Delta*val_fm_y(fm.y,0,0,0) + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    if (un*val_fm_y(fm.y,0,0,0)*s/(val_cm(cm,0,0,0) + 0.) > cfl)
      cfl = un*val_fm_y(fm.y,0,0,0)*s/(val_cm(cm,0,0,0) + 0.);
#line 221 "/home/fpl/softwares/basilisk/src/vof.h"
    double cf;
    if (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.)
      cf = val(c,0,i,0);
    else
      cf = rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.x,0,i,0), _val_higher_dimension(n.y,i,0,0)}, val(alpha,0,i,0),
          (coord){-0.5, -0.5, -0.5},
          (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.y,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; strongif (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = val(c,0,i,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,i,0)/ci + s*min(1., 1. - s*un)*val(gf,0,i,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.y,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_generic()
#line 252
 end_foreach_face(); }
strongif (is_constant(fm.y) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.y.i -_NVARMAX], _constant[fm.x.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_y
#define val_fm_y(a,j,i,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_x
#define val_fm_x(a,j,i,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_cm
#define val_cm(a,j,i,k) val(a,j,i,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,j,i,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,j,i,k)
#line 191
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_y()) {
#line 191
{

#line 191 "/home/fpl/softwares/basilisk/src/vof.h"
 {






    double un = val(uf.y,0,0,0)*dt/(Delta*val_fm_y(fm.y,0,0,0) + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    if (un*val_fm_y(fm.y,0,0,0)*s/(val_cm(cm,0,0,0) + 0.) > cfl)
      cfl = un*val_fm_y(fm.y,0,0,0)*s/(val_cm(cm,0,0,0) + 0.);
#line 221 "/home/fpl/softwares/basilisk/src/vof.h"
    double cf;
    if (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.)
      cf = val(c,0,i,0);
    else
      cf = rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.x,0,i,0), _val_higher_dimension(n.y,i,0,0)}, val(alpha,0,i,0),
          (coord){-0.5, -0.5, -0.5},
          (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.y,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; strongif (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = val(c,0,i,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,i,0)/ci + s*min(1., 1. - s*un)*val(gf,0,i,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.y,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_generic()
#line 252
 end_foreach_face(); }
strongif (!is_constant(fm.y) && is_constant(cm)) {
#undef val_fm_y
#define val_fm_y(a,j,i,k) val(a,j,i,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,j,i,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,j,i,k)
#undef val_fm_x
#define val_fm_x(a,j,i,k) val(a,j,i,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,j,i,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,j,i,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,j,i,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 191
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_y()) {
#line 191
{

#line 191 "/home/fpl/softwares/basilisk/src/vof.h"
 {






    double un = val(uf.y,0,0,0)*dt/(Delta*val_fm_y(fm.y,0,0,0) + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    if (un*val_fm_y(fm.y,0,0,0)*s/(val_cm(cm,0,0,0) + 0.) > cfl)
      cfl = un*val_fm_y(fm.y,0,0,0)*s/(val_cm(cm,0,0,0) + 0.);
#line 221 "/home/fpl/softwares/basilisk/src/vof.h"
    double cf;
    if (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.)
      cf = val(c,0,i,0);
    else
      cf = rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.x,0,i,0), _val_higher_dimension(n.y,i,0,0)}, val(alpha,0,i,0),
          (coord){-0.5, -0.5, -0.5},
          (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.y,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; strongif (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = val(c,0,i,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,i,0)/ci + s*min(1., 1. - s*un)*val(gf,0,i,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.y,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_generic()
#line 252
 end_foreach_face(); }
strongif (is_constant(fm.y) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.y.i -_NVARMAX], _constant[fm.x.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_y
#define val_fm_y(a,j,i,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_x
#define val_fm_x(a,j,i,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,j,i,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 191
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_y()) {
#line 191
{

#line 191 "/home/fpl/softwares/basilisk/src/vof.h"
 {






    double un = val(uf.y,0,0,0)*dt/(Delta*val_fm_y(fm.y,0,0,0) + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    if (un*val_fm_y(fm.y,0,0,0)*s/(val_cm(cm,0,0,0) + 0.) > cfl)
      cfl = un*val_fm_y(fm.y,0,0,0)*s/(val_cm(cm,0,0,0) + 0.);
#line 221 "/home/fpl/softwares/basilisk/src/vof.h"
    double cf;
    if (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.)
      cf = val(c,0,i,0);
    else
      cf = rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.x,0,i,0), _val_higher_dimension(n.y,i,0,0)}, val(alpha,0,i,0),
          (coord){-0.5, -0.5, -0.5},
          (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.y,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; strongif (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = val(c,0,i,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,i,0)/ci + s*min(1., 1. - s*un)*val(gf,0,i,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.y,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_generic()
#line 252
 end_foreach_face(); }mpi_all_reduce_array (&cfl, double, MPI_MAX, 1);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 252
 }
  delete (gfl); pfree (gfl,__func__,__FILE__,__LINE__);




  if (cfl > 0.5 + 1e-6)
    fprintf (ferr,
      "WARNING: CFL must be <= 0.5 for VOF (cfl - 0.5 = %g)\n",
      cfl - 0.5), fflush (ferr);
#line 281 "/home/fpl/softwares/basilisk/src/vof.h"
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dt = dt;
{ double dt = _dt; NOT_UNUSED(dt);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/vof.h", .line = 281,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,j,i,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 281
foreach_stencil(){

#line 281 "/home/fpl/softwares/basilisk/src/vof.h"
 {
    _stencil_val(__FILE__,__LINE__,c,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,flux,0,0,0) - _stencil_val(__FILE__,__LINE__,flux,1,0,0) + _stencil_val(__FILE__,__LINE__,cc,0,0,0)*(_stencil_val(__FILE__,__LINE__,uf.y,1,0,0) - _stencil_val(__FILE__,__LINE__,uf.y,0,0,0)))/(val_cm(cm,0,0,0)*Delta);
    scalar t, tc, tflux;
    scalar * _i11 = tracers; scalar * _i12 = tcl; scalar * _i13 = tfluxl; strongif (tracers) for (t = *tracers, tc = *tcl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i11, tc = *++_i12, tflux = *++_i13)
      _stencil_val(__FILE__,__LINE__,t,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,tflux,0,0,0) - _stencil_val(__FILE__,__LINE__,tflux,1,0,0) + _stencil_val(__FILE__,__LINE__,tc,0,0,0)*(_stencil_val(__FILE__,__LINE__,uf.y,1,0,0) - _stencil_val(__FILE__,__LINE__,uf.y,0,0,0)))/(val_cm(cm,0,0,0)*Delta);
  } } end_foreach_stencil(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,j,i,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 281
foreach_stencil(){

#line 281 "/home/fpl/softwares/basilisk/src/vof.h"
 {
    _stencil_val(__FILE__,__LINE__,c,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,flux,0,0,0) - _stencil_val(__FILE__,__LINE__,flux,1,0,0) + _stencil_val(__FILE__,__LINE__,cc,0,0,0)*(_stencil_val(__FILE__,__LINE__,uf.y,1,0,0) - _stencil_val(__FILE__,__LINE__,uf.y,0,0,0)))/(val_cm(cm,0,0,0)*Delta);
    scalar t, tc, tflux;
    scalar * _i11 = tracers; scalar * _i12 = tcl; scalar * _i13 = tfluxl; strongif (tracers) for (t = *tracers, tc = *tcl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i11, tc = *++_i12, tflux = *++_i13)
      _stencil_val(__FILE__,__LINE__,t,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,tflux,0,0,0) - _stencil_val(__FILE__,__LINE__,tflux,1,0,0) + _stencil_val(__FILE__,__LINE__,tc,0,0,0)*(_stencil_val(__FILE__,__LINE__,uf.y,1,0,0) - _stencil_val(__FILE__,__LINE__,uf.y,0,0,0)))/(val_cm(cm,0,0,0)*Delta);
  } } end_foreach_stencil(); } if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/fpl/softwares/basilisk/src/vof.h", 281, "dt");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 286

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,j,i,k) val(a,j,i,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,j,i,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,j,i,k)
#line 281
foreach(){

#line 281 "/home/fpl/softwares/basilisk/src/vof.h"
 {
    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,0,1,0) + val(cc,0,0,0)*(val(uf.y,0,1,0) - val(uf.y,0,0,0)))/(val_cm(cm,0,0,0)*Delta);
    scalar t, tc, tflux;
    scalar * _i11 = tracers; scalar * _i12 = tcl; scalar * _i13 = tfluxl; strongif (tracers) for (t = *tracers, tc = *tcl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i11, tc = *++_i12, tflux = *++_i13)
      val(t,0,0,0) += dt*(val(tflux,0,0,0) - val(tflux,0,1,0) + val(tc,0,0,0)*(val(uf.y,0,1,0) - val(uf.y,0,0,0)))/(val_cm(cm,0,0,0)*Delta);
  } } end_foreach(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,j,i,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 281
foreach(){

#line 281 "/home/fpl/softwares/basilisk/src/vof.h"
 {
    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,0,1,0) + val(cc,0,0,0)*(val(uf.y,0,1,0) - val(uf.y,0,0,0)))/(val_cm(cm,0,0,0)*Delta);
    scalar t, tc, tflux;
    scalar * _i11 = tracers; scalar * _i12 = tcl; scalar * _i13 = tfluxl; strongif (tracers) for (t = *tracers, tc = *tcl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i11, tc = *++_i12, tflux = *++_i13)
      val(t,0,0,0) += dt*(val(tflux,0,0,0) - val(tflux,0,1,0) + val(tc,0,0,0)*(val(uf.y,0,1,0) - val(uf.y,0,0,0)))/(val_cm(cm,0,0,0)*Delta);
  } } end_foreach(); } }
#line 304 "/home/fpl/softwares/basilisk/src/vof.h"
  delete (tfluxl); pfree (tfluxl,__func__,__FILE__,__LINE__);
 delete (((scalar []){flux,alpha,n.x,n.y,{-1}})); }






void vof_advection (scalar * interfaces, int i)
{
  strongif (interfaces) for (scalar c = *interfaces, *_i126 = interfaces; ((scalar *)&c)->i >= 0; c = *++_i126) {
#line 324 "/home/fpl/softwares/basilisk/src/vof.h"
    scalar cc= new_scalar("cc"), * tcl = NULL, * tracers = _attribute[c.i].tracers;
    strongif (tracers) for (scalar t = *tracers, *_i127 = tracers; ((scalar *)&t)->i >= 0; t = *++_i127) {
      scalar tc = new_scalar("tc");
      tcl = list_append (tcl, tc);

      if (_attribute[t.i].refine != vof_concentration_refine) {
 _attribute[t.i].refine = _attribute[t.i].prolongation = vof_concentration_refine;
 _attribute[t.i].restriction = restriction_volume_average;
 _attribute[t.i].dirty = true;
 _attribute[t.i].c = c;
      }

    }
     { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/vof.h", .line = 337,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 337 "/home/fpl/softwares/basilisk/src/vof.h"
 {
      _stencil_val(__FILE__,__LINE__,cc,0,0,0) = (_stencil_val(__FILE__,__LINE__,c,0,0,0) > 0.5);
      scalar t, tc;
      scalar * _i14 = tracers; scalar * _i15 = tcl; strongif (tracers) for (t = *tracers, tc = *tcl; ((scalar *)&t)->i >= 0; t = *++_i14, tc = *++_i15) {
 IF (_attribute[t.i].inverse)
   _stencil_val(__FILE__,__LINE__,tc,0,0,0) = _stencil_val(__FILE__,__LINE__,c,0,0,0) < 0.5 ? _stencil_val(__FILE__,__LINE__,t,0,0,0)/(1. - _stencil_val(__FILE__,__LINE__,c,0,0,0)) : 0.;
 
   _stencil_val(__FILE__,__LINE__,tc,0,0,0) = _stencil_val(__FILE__,__LINE__,c,0,0,0) > 0.5 ? _stencil_val(__FILE__,__LINE__,t,0,0,0)/_stencil_val(__FILE__,__LINE__,c,0,0,0) : 0.;
      }
    } } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 346
foreach(){

#line 337 "/home/fpl/softwares/basilisk/src/vof.h"
 {
      val(cc,0,0,0) = (val(c,0,0,0) > 0.5);
      scalar t, tc;
      scalar * _i14 = tracers; scalar * _i15 = tcl; strongif (tracers) for (t = *tracers, tc = *tcl; ((scalar *)&t)->i >= 0; t = *++_i14, tc = *++_i15) {
 if (_attribute[t.i].inverse)
   val(tc,0,0,0) = val(c,0,0,0) < 0.5 ? val(t,0,0,0)/(1. - val(c,0,0,0)) : 0.;
 else
   val(tc,0,0,0) = val(c,0,0,0) > 0.5 ? val(t,0,0,0)/val(c,0,0,0) : 0.;
      }
    } } end_foreach(); }






    void (* sweep[2]) (scalar, scalar, scalar *);
    int d = 0;
    {
#line 355

      sweep[d++] = sweep_x;
#line 355

      sweep[d++] = sweep_y;}
    for (d = 0; d < 2; d++)
      sweep[(i + d) % 2] (c, cc, tcl);
    delete (tcl), pfree (tcl,__func__,__FILE__,__LINE__);
   delete (((scalar []){cc,{-1}})); }
}

static int vof_1_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int vof_1 (const int i, const double t, Event * _ev) { trace ("vof_1", "/home/fpl/softwares/basilisk/src/vof.h", 363); 
  vof_advection (interfaces, i); end_trace("vof_1", "/home/fpl/softwares/basilisk/src/vof.h", 364);  return 0; } 
#line 14 "/home/fpl/softwares/basilisk/src/two-phase.h"

scalar f= {8}, * interfaces = ((scalar []){{8},{-1}});
double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0.;





vector alphav= {{9},{10}};
scalar rhov= {11};

static int defaults_4_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults_4 (const int i, const double t, Event * _ev) { trace ("defaults_4", "/home/fpl/softwares/basilisk/src/two-phase.h", 25);  {
  alpha = alphav;
  rho = rhov;





  if (mu1 || mu2)
    mu = new_face_vector("mu");




  display ((struct _display){"draw_vof (c = 'f');"});
 end_trace("defaults_4", "/home/fpl/softwares/basilisk/src/two-phase.h", 40); } return 0; } 
#line 64 "/home/fpl/softwares/basilisk/src/two-phase.h"
static int tracer_advection_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int tracer_advection_0 (const int i, const double t, Event * _ev) { trace ("tracer_advection_0", "/home/fpl/softwares/basilisk/src/two-phase.h", 64); 
{
#line 90 "/home/fpl/softwares/basilisk/src/two-phase.h"
  _attribute[f.i].prolongation = refine_bilinear;
  _attribute[f.i].dirty = true;

 end_trace("tracer_advection_0", "/home/fpl/softwares/basilisk/src/two-phase.h", 93); } return 0; } 

static int properties_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int properties_0 (const int i, const double t, Event * _ev) { trace ("properties_0", "/home/fpl/softwares/basilisk/src/two-phase.h", 95); 
{
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _rho1 = rho1;
 double _rho2 = rho2;
 double _mu1 = mu1;
 double _mu2 = mu2;
{ double rho1 = _rho1; NOT_UNUSED(rho1);
 double rho2 = _rho2; NOT_UNUSED(rho2);
 double mu1 = _mu1; NOT_UNUSED(mu1);
 double mu2 = _mu2; NOT_UNUSED(mu2);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/two-phase.h", .line = 97,
    .each = "foreach_face", .first = _first_call
  };

strongif (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 97
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 97
{

#line 97 "/home/fpl/softwares/basilisk/src/two-phase.h"
 {
    double ff = (_stencil_val(__FILE__,__LINE__,f,0,0,0) + _stencil_val(__FILE__,__LINE__,f,-1,0,0))/2.;
    _stencil_val(__FILE__,__LINE__,alphav.x,0,0,0) = val_fm_x(fm.x,0,0,0)/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    IF (mu1 || mu2) {
      vector muv = mu;
      _stencil_val(__FILE__,__LINE__,muv.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 97
{

#line 97 "/home/fpl/softwares/basilisk/src/two-phase.h"
 {
    double ff = (_stencil_val(__FILE__,__LINE__,f,0,0,0) + _stencil_val(__FILE__,__LINE__,f,0,-1,0))/2.;
    _stencil_val(__FILE__,__LINE__,alphav.y,0,0,0) = val_fm_y(fm.y,0,0,0)/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    IF (mu1 || mu2) {
      vector muv = mu;
      _stencil_val(__FILE__,__LINE__,muv.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  } }  }}  end_foreach_face_stencil()
#line 104
 }
strongif (is_constant(fm.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 97
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 97
{

#line 97 "/home/fpl/softwares/basilisk/src/two-phase.h"
 {
    double ff = (_stencil_val(__FILE__,__LINE__,f,0,0,0) + _stencil_val(__FILE__,__LINE__,f,-1,0,0))/2.;
    _stencil_val(__FILE__,__LINE__,alphav.x,0,0,0) = val_fm_x(fm.x,0,0,0)/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    IF (mu1 || mu2) {
      vector muv = mu;
      _stencil_val(__FILE__,__LINE__,muv.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 97
{

#line 97 "/home/fpl/softwares/basilisk/src/two-phase.h"
 {
    double ff = (_stencil_val(__FILE__,__LINE__,f,0,0,0) + _stencil_val(__FILE__,__LINE__,f,0,-1,0))/2.;
    _stencil_val(__FILE__,__LINE__,alphav.y,0,0,0) = val_fm_y(fm.y,0,0,0)/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    IF (mu1 || mu2) {
      vector muv = mu;
      _stencil_val(__FILE__,__LINE__,muv.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  } }  }}  end_foreach_face_stencil()
#line 104
 } if (_first_call) {
 if (rho1 != _rho1)
   reduction_warning ("/home/fpl/softwares/basilisk/src/two-phase.h", 97, "rho1");
 }
 if (_first_call) {
 if (rho2 != _rho2)
   reduction_warning ("/home/fpl/softwares/basilisk/src/two-phase.h", 97, "rho2");
 }
 if (_first_call) {
 if (mu1 != _mu1)
   reduction_warning ("/home/fpl/softwares/basilisk/src/two-phase.h", 97, "mu1");
 }
 if (_first_call) {
 if (mu2 != _mu2)
   reduction_warning ("/home/fpl/softwares/basilisk/src/two-phase.h", 97, "mu2");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 104

strongif (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#line 97
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 97
{

#line 97 "/home/fpl/softwares/basilisk/src/two-phase.h"
 {
    double ff = (val(f,0,0,0) + val(f,-1,0,0))/2.;
    val(alphav.x,0,0,0) = val_fm_x(fm.x,0,0,0)/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    if (mu1 || mu2) {
      vector muv = mu;
      val(muv.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 97
{

#line 97 "/home/fpl/softwares/basilisk/src/two-phase.h"
 {
    double ff = (val(f,0,0,0) + val(f,0,-1,0))/2.;
    val(alphav.y,0,0,0) = val_fm_y(fm.y,0,0,0)/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    if (mu1 || mu2) {
      vector muv = mu;
      val(muv.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  } }  }}  end_foreach_face_generic()
#line 104
 end_foreach_face(); }
strongif (is_constant(fm.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 97
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 97
{

#line 97 "/home/fpl/softwares/basilisk/src/two-phase.h"
 {
    double ff = (val(f,0,0,0) + val(f,-1,0,0))/2.;
    val(alphav.x,0,0,0) = val_fm_x(fm.x,0,0,0)/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    if (mu1 || mu2) {
      vector muv = mu;
      val(muv.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 97
{

#line 97 "/home/fpl/softwares/basilisk/src/two-phase.h"
 {
    double ff = (val(f,0,0,0) + val(f,0,-1,0))/2.;
    val(alphav.y,0,0,0) = val_fm_y(fm.y,0,0,0)/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    if (mu1 || mu2) {
      vector muv = mu;
      val(muv.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  } }  }}  end_foreach_face_generic()
#line 104
 end_foreach_face(); } }

   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _rho1 = rho1;
 double _rho2 = rho2;
{ double rho1 = _rho1; NOT_UNUSED(rho1);
 double rho2 = _rho2; NOT_UNUSED(rho2);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/fpl/softwares/basilisk/src/two-phase.h", .line = 106,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 106
foreach_stencil(){

#line 106 "/home/fpl/softwares/basilisk/src/two-phase.h"

    _stencil_val(__FILE__,__LINE__,rhov,0,0,0) = val_cm(cm,0,0,0)*(clamp(_stencil_val(__FILE__,__LINE__,f,0,0,0),0.,1.)*(rho1 - rho2) + rho2); } end_foreach_stencil(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 106
foreach_stencil(){

#line 106 "/home/fpl/softwares/basilisk/src/two-phase.h"

    _stencil_val(__FILE__,__LINE__,rhov,0,0,0) = val_cm(cm,0,0,0)*(clamp(_stencil_val(__FILE__,__LINE__,f,0,0,0),0.,1.)*(rho1 - rho2) + rho2); } end_foreach_stencil(); } if (_first_call) {
 if (rho1 != _rho1)
   reduction_warning ("/home/fpl/softwares/basilisk/src/two-phase.h", 106, "rho1");
 }
 if (_first_call) {
 if (rho2 != _rho2)
   reduction_warning ("/home/fpl/softwares/basilisk/src/two-phase.h", 106, "rho2");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 107

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 106
foreach(){

#line 106 "/home/fpl/softwares/basilisk/src/two-phase.h"

    val(rhov,0,0,0) = val_cm(cm,0,0,0)*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2); } end_foreach(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 106
foreach(){

#line 106 "/home/fpl/softwares/basilisk/src/two-phase.h"

    val(rhov,0,0,0) = val_cm(cm,0,0,0)*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2); } end_foreach(); } }


  _attribute[f.i].prolongation = fraction_refine;
  _attribute[f.i].dirty = true;

 end_trace("properties_0", "/home/fpl/softwares/basilisk/src/two-phase.h", 113); } return 0; } 
#line 24 "60d_plate_advancing_simulation.c"



int maxlevel = 9;
char name_vtk[100];
double U0;
double H0;
#line 60 "60d_plate_advancing_simulation.c"
double h0;

vector h= {{12},{13}};
double theta0 ;


static void _set_boundary4 (void) { _attribute[uf.x.i].boundary[left] = _boundary4; _attribute[uf.x.i].boundary_homogeneous[left] = _boundary4_homogeneous; _attribute[uf.x.i].dirty = true; } 
static void _set_boundary5 (void) { _attribute[uf.x.i].boundary[right] = _boundary5; _attribute[uf.x.i].boundary_homogeneous[right] = _boundary5_homogeneous; _attribute[uf.x.i].dirty = true; } 
static void _set_boundary6 (void) { _attribute[uf.x.i].boundary[top] = _boundary6; _attribute[uf.x.i].boundary_homogeneous[top] = _boundary6_homogeneous; _attribute[uf.x.i].dirty = true; } 
static void _set_boundary7 (void) { _attribute[uf.x.i].boundary[bottom] = _boundary7; _attribute[uf.x.i].boundary_homogeneous[bottom] = _boundary7_homogeneous; _attribute[uf.x.i].dirty = true; } 

int padding=6;
int main ()
{ _init_solver();
        L0 = 0.015;
        U0 = -0.001 ;
 origin ((struct _origin){-L0/2, -L0/2});
 N = 128;

        _attribute[f.i].sigma = 0.021;
        _attribute[f.i].height = h;
        ;

        theta0 = 120*pi/180.0;
        _attribute[h.y.i].boundary[top] = _boundary8; _attribute[h.y.i].boundary_homogeneous[top] = _boundary8_homogeneous; _attribute[h.y.i].dirty = true;
        _attribute[h.y.i].boundary[bottom] = _boundary9; _attribute[h.y.i].boundary_homogeneous[bottom] = _boundary9_homogeneous; _attribute[h.y.i].dirty = true;


        rho2 = 940;
        mu2 = 0.0088;
 rho1 = 1.2;
 mu1 = 0.0000181;
        run();

 free_solver(); }


static int init_1_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (t = 0);   *ip = i; *tp = t;   return ret; } static int init_1 (const int i, const double t, Event * _ev) { trace ("init_1", "60d_plate_advancing_simulation.c", 97); 
{






        do { scalar phi= new_vertex_scalar("phi");  { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _theta0 = theta0;
{ double theta0 = _theta0; NOT_UNUSED(theta0);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "60d_plate_advancing_simulation.c", .line = 105,
    .each = "foreach_vertex", .first = _first_call
  };
foreach_vertex_stencil(){

#line 105 "60d_plate_advancing_simulation.c"
 _stencil_val(__FILE__,__LINE__,phi,0,0,0) = ( x + 1.51e-3/(tan(theta0)*exp((-y+ 0.0075)/1.51e-3))); } end_foreach_vertex_stencil(); if (_first_call) {
 if (theta0 != _theta0)
   reduction_warning ("60d_plate_advancing_simulation.c", 105, "theta0");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 105
foreach_vertex(){

#line 105 "60d_plate_advancing_simulation.c"
 val(phi,0,0,0) = ( x + 1.51e-3/(tan(theta0)*exp((-y+ 0.0075)/1.51e-3))); } end_foreach_vertex(); } fractions ((struct Fractions){phi, f});  delete (((scalar []){phi,{-1}})); } while(0);
        boundary_internal ((scalar *)(((scalar []){f,{-1}})), "60d_plate_advancing_simulation.c", 106);
 end_trace("init_1", "60d_plate_advancing_simulation.c", 107); } return 0; } 




static int acceleration_2_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int acceleration_2 (const int i, const double t, Event * _ev) { trace ("acceleration_2", "60d_plate_advancing_simulation.c", 112); 
{
        vector av = a;
         { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "60d_plate_advancing_simulation.c", .line = 115,
    .each = "foreach_face", .first = _first_call
  };
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 115
{

#line 115 "60d_plate_advancing_simulation.c"

                _stencil_val(__FILE__,__LINE__,av.x,0,0,0) = -9.81; }  }}  end_foreach_face_stencil()
#line 116
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 116
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 115
{

#line 115 "60d_plate_advancing_simulation.c"

                val(av.x,0,0,0) = -9.81; }  }}  end_foreach_face_generic()
#line 116
 end_foreach_face(); }
 end_trace("acceleration_2", "60d_plate_advancing_simulation.c", 117); } return 0; } 


static void _set_boundary10 (void) { _attribute[u.x.i].boundary[left] = _boundary10; _attribute[u.x.i].boundary_homogeneous[left] = _boundary10_homogeneous; _attribute[u.x.i].dirty = true; } 
static void _set_boundary11 (void) { _attribute[u.y.i].boundary[left] = _boundary11; _attribute[u.y.i].boundary_homogeneous[left] = _boundary11_homogeneous; _attribute[u.y.i].dirty = true; } 


static void _set_boundary12 (void) { _attribute[u.x.i].boundary[right] = _boundary12; _attribute[u.x.i].boundary_homogeneous[right] = _boundary12_homogeneous; _attribute[u.x.i].dirty = true; } 
static void _set_boundary13 (void) { _attribute[u.y.i].boundary[right] = _boundary13; _attribute[u.y.i].boundary_homogeneous[right] = _boundary13_homogeneous; _attribute[u.y.i].dirty = true; } 


static void _set_boundary14 (void) { _attribute[u.x.i].boundary[top] = _boundary14; _attribute[u.x.i].boundary_homogeneous[top] = _boundary14_homogeneous; _attribute[u.x.i].dirty = true; } 
static void _set_boundary15 (void) { _attribute[u.y.i].boundary[top] = _boundary15; _attribute[u.y.i].boundary_homogeneous[top] = _boundary15_homogeneous; _attribute[u.y.i].dirty = true; } 

static void _set_boundary16 (void) { _attribute[u.x.i].boundary[bottom] = _boundary16; _attribute[u.x.i].boundary_homogeneous[bottom] = _boundary16_homogeneous; _attribute[u.x.i].dirty = true; } 
static void _set_boundary17 (void) { _attribute[u.y.i].boundary[bottom] = _boundary17; _attribute[u.y.i].boundary_homogeneous[bottom] = _boundary17_homogeneous; _attribute[u.y.i].dirty = true; } 




static int logfile_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i+=50);   *ip = i; *tp = t;   return ret; } static int logfile (const int i, const double t, Event * _ev) { trace ("logfile", "60d_plate_advancing_simulation.c", 137); 
        fprintf (ferr, "%d %g\n", i, t); end_trace("logfile", "60d_plate_advancing_simulation.c", 138);  return 0; } 







static int dumpfile_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (t = 0.);   *ip = i; *tp = t;   return ret; } static int dumpfile_expr1 (int * ip, double * tp, Event * _ev) {   int i = *ip; double t = *tp;   int ret = (t += 0.01);   *ip = i; *tp = t;   return ret; } static int dumpfile_expr2 (int * ip, double * tp, Event * _ev) {   int i = *ip; double t = *tp;   int ret = (t <= 10);   *ip = i; *tp = t;   return ret; } static int dumpfile (const int i, const double t, Event * _ev) { trace ("dumpfile", "60d_plate_advancing_simulation.c", 146); 
{
    char name[80];
    sprintf(name,"dump-%g",t*100);
    scalar pid= new_scalar("pid");
     { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "60d_plate_advancing_simulation.c", .line = 151,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 151 "60d_plate_advancing_simulation.c"

    {
        _stencil_val(__FILE__,__LINE__,pid,0,0,0) = fmod(pid()*(npe() + 37),npe());
    } } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 154
foreach(){

#line 151 "60d_plate_advancing_simulation.c"

    {
        val(pid,0,0,0) = fmod(pid()*(npe() + 37),npe());
    } } end_foreach(); }
    boundary_internal ((scalar *)(((scalar []){pid,{-1}})), "60d_plate_advancing_simulation.c", 155);
    dump ((struct Dump){name});




        sprintf (name, "bview-%g",t*100);
        dump ((struct Dump){name});

 delete (((scalar []){pid,{-1}}));  end_trace("dumpfile", "60d_plate_advancing_simulation.c", 164); } return 0; } 
#line 192 "60d_plate_advancing_simulation.c"
static int adapt_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i += 5);   *ip = i; *tp = t;   return ret; } static int adapt_0 (const int i, const double t, Event * _ev) { trace ("adapt_0", "60d_plate_advancing_simulation.c", 192);  {
  adapt_wavelet ((struct Adapt){(scalar*)((scalar []){f,u.x,u.y,{-1}}), (double[]){0.0001,0.05,0.05},maxlevel , maxlevel-3 });
 end_trace("adapt_0", "60d_plate_advancing_simulation.c", 194); } return 0; } 
#line 91 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
static double _boundary0 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 90 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

strongif (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 91
return  neumann ((val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); }
strongif (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 91
return  neumann ((val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); }
strongif (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 91
return  neumann ((val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); }
strongif (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 91
return  neumann ((val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); }
strongif (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 91
return  neumann ((val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); }
strongif (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 91
return  neumann ((val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); }
strongif (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 91
return  neumann ((val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); }
strongif (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 91
return  neumann ((val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); } return 0.; } static double _boundary0_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 90 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

strongif (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 91
return  neumann_homogeneous(); }
strongif (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 91
return  neumann_homogeneous(); }
strongif (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 91
return  neumann_homogeneous(); }
strongif (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 91
return  neumann_homogeneous(); }
strongif (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 91
return  neumann_homogeneous(); }
strongif (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 91
return  neumann_homogeneous(); }
strongif (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 91
return  neumann_homogeneous(); }
strongif (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 91
return  neumann_homogeneous(); } return 0.; }
#line 92 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
static double _boundary1 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 91 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

strongif (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann (- (val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0))); }
strongif (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann (- (val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0))); }
strongif (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann (- (val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0))); }
strongif (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann (- (val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0))); }
strongif (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 92
return  neumann (- (val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0))); }
strongif (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 92
return  neumann (- (val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0))); }
strongif (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 92
return  neumann (- (val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0))); }
strongif (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 92
return  neumann (- (val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0))); } return 0.; } static double _boundary1_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 91 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

strongif (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann_homogeneous(); }
strongif (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann_homogeneous(); }
strongif (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann_homogeneous(); }
strongif (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann_homogeneous(); }
strongif (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 92
return  neumann_homogeneous(); }
strongif (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 92
return  neumann_homogeneous(); }
strongif (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 92
return  neumann_homogeneous(); }
strongif (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 92
return  neumann_homogeneous(); } return 0.; }
#line 101 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
static double _boundary2 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 100 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

strongif (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 101
return  neumann ((val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); }
strongif (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 101
return  neumann ((val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); }
strongif (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 101
return  neumann ((val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); }
strongif (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 101
return  neumann ((val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); }
strongif (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 101
return  neumann ((val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); }
strongif (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 101
return  neumann ((val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); }
strongif (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 101
return  neumann ((val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); }
strongif (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 101
return  neumann ((val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); } return 0.; } static double _boundary2_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 100 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

strongif (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 101
return  neumann_homogeneous(); }
strongif (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 101
return  neumann_homogeneous(); }
strongif (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 101
return  neumann_homogeneous(); }
strongif (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 101
return  neumann_homogeneous(); }
strongif (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 101
return  neumann_homogeneous(); }
strongif (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 101
return  neumann_homogeneous(); }
strongif (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 101
return  neumann_homogeneous(); }
strongif (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 101
return  neumann_homogeneous(); } return 0.; }
#line 102 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"
static double _boundary3 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 101 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

strongif (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 102
return  neumann (- (val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0))); }
strongif (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 102
return  neumann (- (val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0))); }
strongif (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 102
return  neumann (- (val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0))); }
strongif (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 102
return  neumann (- (val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0))); }
strongif (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 102
return  neumann (- (val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0))); }
strongif (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 102
return  neumann (- (val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0))); }
strongif (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 102
return  neumann (- (val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0))); }
strongif (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 102
return  neumann (- (val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0))); } return 0.; } static double _boundary3_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 101 "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h"

strongif (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 102
return  neumann_homogeneous(); }
strongif (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 102
return  neumann_homogeneous(); }
strongif (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 102
return  neumann_homogeneous(); }
strongif (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 102
return  neumann_homogeneous(); }
strongif (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 102
return  neumann_homogeneous(); }
strongif (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 102
return  neumann_homogeneous(); }
strongif (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 102
return  neumann_homogeneous(); }
strongif (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 102
return  neumann_homogeneous(); } return 0.; }
#line 66 "60d_plate_advancing_simulation.c"
static double _boundary4 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 65 "60d_plate_advancing_simulation.c"
return  0.; return 0.; } static double _boundary4_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 65 "60d_plate_advancing_simulation.c"
return  0.; return 0.; }
#line 67 "60d_plate_advancing_simulation.c"
static double _boundary5 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 66 "60d_plate_advancing_simulation.c"
return  0.; return 0.; } static double _boundary5_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 66 "60d_plate_advancing_simulation.c"
return  0.; return 0.; }
#line 68 "60d_plate_advancing_simulation.c"
static double _boundary6 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 67 "60d_plate_advancing_simulation.c"
return  0.; return 0.; } static double _boundary6_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 67 "60d_plate_advancing_simulation.c"
return  0.; return 0.; }
#line 69 "60d_plate_advancing_simulation.c"
static double _boundary7 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 68 "60d_plate_advancing_simulation.c"
return  0.; return 0.; } static double _boundary7_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 68 "60d_plate_advancing_simulation.c"
return  0.; return 0.; }
#line 84 "60d_plate_advancing_simulation.c"
static double _boundary8 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 83 "60d_plate_advancing_simulation.c"
return  (val(_s,0,0,0) == nodata ? nodata : val(_s,0,0,0) + (orientation(val(_s,0,0,0)) ? -1. : 1.)/tan(theta0)); return 0.; } static double _boundary8_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 83 "60d_plate_advancing_simulation.c"
return  (val(_s,0,0,0) == nodata ? nodata : val(_s,0,0,0) + (orientation(val(_s,0,0,0)) ? -1. : 1.)/tan(theta0)); return 0.; }
#line 85 "60d_plate_advancing_simulation.c"
static double _boundary9 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 84 "60d_plate_advancing_simulation.c"
return  (val(_s,0,0,0) == nodata ? nodata : val(_s,0,0,0) + (orientation(val(_s,0,0,0)) ? -1. : 1.)/tan(pi/2)); return 0.; } static double _boundary9_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 84 "60d_plate_advancing_simulation.c"
return  (val(_s,0,0,0) == nodata ? nodata : val(_s,0,0,0) + (orientation(val(_s,0,0,0)) ? -1. : 1.)/tan(pi/2)); return 0.; }
#line 120 "60d_plate_advancing_simulation.c"
static double _boundary10 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 119 "60d_plate_advancing_simulation.c"
return  dirichlet(0.); return 0.; } static double _boundary10_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 119 "60d_plate_advancing_simulation.c"
return  dirichlet_homogeneous(); return 0.; }
#line 121 "60d_plate_advancing_simulation.c"
static double _boundary11 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 120 "60d_plate_advancing_simulation.c"
return  dirichlet(0.0); return 0.; } static double _boundary11_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 120 "60d_plate_advancing_simulation.c"
return  dirichlet_homogeneous(); return 0.; }
#line 124 "60d_plate_advancing_simulation.c"
static double _boundary12 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 123 "60d_plate_advancing_simulation.c"
return  dirichlet(0.); return 0.; } static double _boundary12_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 123 "60d_plate_advancing_simulation.c"
return  dirichlet_homogeneous(); return 0.; }
#line 125 "60d_plate_advancing_simulation.c"
static double _boundary13 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 124 "60d_plate_advancing_simulation.c"
return  dirichlet(0.); return 0.; } static double _boundary13_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 124 "60d_plate_advancing_simulation.c"
return  dirichlet_homogeneous(); return 0.; }
#line 128 "60d_plate_advancing_simulation.c"
static double _boundary14 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 127 "60d_plate_advancing_simulation.c"
return  dirichlet(0.); return 0.; } static double _boundary14_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 127 "60d_plate_advancing_simulation.c"
return  dirichlet_homogeneous(); return 0.; }
#line 129 "60d_plate_advancing_simulation.c"
static double _boundary15 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 128 "60d_plate_advancing_simulation.c"
return  dirichlet(0.0001); return 0.; } static double _boundary15_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 128 "60d_plate_advancing_simulation.c"
return  dirichlet_homogeneous(); return 0.; }
#line 131 "60d_plate_advancing_simulation.c"
static double _boundary16 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 130 "60d_plate_advancing_simulation.c"
return  dirichlet(0.); return 0.; } static double _boundary16_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 130 "60d_plate_advancing_simulation.c"
return  dirichlet_homogeneous(); return 0.; }
#line 132 "60d_plate_advancing_simulation.c"
static double _boundary17 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 131 "60d_plate_advancing_simulation.c"
return  dirichlet(0.0); return 0.; } static double _boundary17_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  strongif (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  strongif (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 131 "60d_plate_advancing_simulation.c"
return  dirichlet_homogeneous(); return 0.; }
size_t datasize = 14*sizeof (double);
static int defaults (const int i, const double t, Event * _ev);
static int defaults_expr0 (int * ip, double * tp, Event * _ev);
static int cleanup (const int i, const double t, Event * _ev);
static int cleanup_expr0 (int * ip, double * tp, Event * _ev);
static int defaults_0 (const int i, const double t, Event * _ev);
static int defaults_0_expr0 (int * ip, double * tp, Event * _ev);
static int default_display (const int i, const double t, Event * _ev);
static int default_display_expr0 (int * ip, double * tp, Event * _ev);
static int init (const int i, const double t, Event * _ev);
static int init_expr0 (int * ip, double * tp, Event * _ev);
static int set_dtmax (const int i, const double t, Event * _ev);
static int set_dtmax_expr0 (int * ip, double * tp, Event * _ev);
static int stability (const int i, const double t, Event * _ev);
static int stability_expr0 (int * ip, double * tp, Event * _ev);
static int vof (const int i, const double t, Event * _ev);
static int vof_expr0 (int * ip, double * tp, Event * _ev);
static int tracer_advection (const int i, const double t, Event * _ev);
static int tracer_advection_expr0 (int * ip, double * tp, Event * _ev);
static int tracer_diffusion (const int i, const double t, Event * _ev);
static int tracer_diffusion_expr0 (int * ip, double * tp, Event * _ev);
static int properties (const int i, const double t, Event * _ev);
static int properties_expr0 (int * ip, double * tp, Event * _ev);
static int advection_term (const int i, const double t, Event * _ev);
static int advection_term_expr0 (int * ip, double * tp, Event * _ev);
static int viscous_term (const int i, const double t, Event * _ev);
static int viscous_term_expr0 (int * ip, double * tp, Event * _ev);
static int acceleration (const int i, const double t, Event * _ev);
static int acceleration_expr0 (int * ip, double * tp, Event * _ev);
static int projection (const int i, const double t, Event * _ev);
static int projection_expr0 (int * ip, double * tp, Event * _ev);
static int end_timestep (const int i, const double t, Event * _ev);
static int end_timestep_expr0 (int * ip, double * tp, Event * _ev);
static int adapt (const int i, const double t, Event * _ev);
static int adapt_expr0 (int * ip, double * tp, Event * _ev);
static int init_0 (const int i, const double t, Event * _ev);
static int init_0_expr0 (int * ip, double * tp, Event * _ev);
static int vof_0 (const int i, const double t, Event * _ev);
static int vof_0_expr0 (int * ip, double * tp, Event * _ev);
static int defaults_1 (const int i, const double t, Event * _ev);
static int defaults_1_expr0 (int * ip, double * tp, Event * _ev);
static int acceleration_0 (const int i, const double t, Event * _ev);
static int acceleration_0_expr0 (int * ip, double * tp, Event * _ev);
static int stability_0 (const int i, const double t, Event * _ev);
static int stability_0_expr0 (int * ip, double * tp, Event * _ev);
static int acceleration_1 (const int i, const double t, Event * _ev);
static int acceleration_1_expr0 (int * ip, double * tp, Event * _ev);
static int defaults_2 (const int i, const double t, Event * _ev);
static int defaults_2_expr0 (int * ip, double * tp, Event * _ev);
static int defaults_3 (const int i, const double t, Event * _ev);
static int defaults_3_expr0 (int * ip, double * tp, Event * _ev);
static int stability_1 (const int i, const double t, Event * _ev);
static int stability_1_expr0 (int * ip, double * tp, Event * _ev);
static int vof_1 (const int i, const double t, Event * _ev);
static int vof_1_expr0 (int * ip, double * tp, Event * _ev);
static int defaults_4 (const int i, const double t, Event * _ev);
static int defaults_4_expr0 (int * ip, double * tp, Event * _ev);
static int tracer_advection_0 (const int i, const double t, Event * _ev);
static int tracer_advection_0_expr0 (int * ip, double * tp, Event * _ev);
static int properties_0 (const int i, const double t, Event * _ev);
static int properties_0_expr0 (int * ip, double * tp, Event * _ev);
static int init_1 (const int i, const double t, Event * _ev);
static int init_1_expr0 (int * ip, double * tp, Event * _ev);
static int acceleration_2 (const int i, const double t, Event * _ev);
static int acceleration_2_expr0 (int * ip, double * tp, Event * _ev);
static int logfile (const int i, const double t, Event * _ev);
static int logfile_expr0 (int * ip, double * tp, Event * _ev);
static int dumpfile (const int i, const double t, Event * _ev);
static int dumpfile_expr0 (int * ip, double * tp, Event * _ev);
static int dumpfile_expr1 (int * ip, double * tp, Event * _ev);
static int dumpfile_expr2 (int * ip, double * tp, Event * _ev);
static int adapt_0 (const int i, const double t, Event * _ev);
static int adapt_0_expr0 (int * ip, double * tp, Event * _ev);
static void _set_boundary0 (void);
static void _set_boundary1 (void);
static void _set_boundary2 (void);
static void _set_boundary3 (void);
static void _set_boundary4 (void);
static void _set_boundary5 (void);
static void _set_boundary6 (void);
static void _set_boundary7 (void);
static void _set_boundary10 (void);
static void _set_boundary11 (void);
static void _set_boundary12 (void);
static void _set_boundary13 (void);
static void _set_boundary14 (void);
static void _set_boundary15 (void);
static void _set_boundary16 (void);
static void _set_boundary17 (void);
void _init_solver (void) {
  void init_solver();
  init_solver();
  Events = (Event *) pmalloc (sizeof (Event), __func__, __FILE__, __LINE__);
  Events[0].last = 1;
  event_register ((Event){ 0, 1, defaults, {defaults_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/run.h", 42, "defaults"});
  event_register ((Event){ 0, 1, defaults_0, {defaults_0_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 126, "defaults"});
  event_register ((Event){ 0, 1, defaults_1, {defaults_1_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/iforce.h", 30, "defaults"});
  event_register ((Event){ 0, 1, defaults_2, {defaults_2_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/vof.h", 103, "defaults"});
  event_register ((Event){ 0, 1, defaults_3, {defaults_3_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/vof.h", 123, "defaults"});
  event_register ((Event){ 0, 1, defaults_4, {defaults_4_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/two-phase.h", 25, "defaults"});
  event_register ((Event){ 0, 1, default_display, {default_display_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 179, "default_display"});
  event_register ((Event){ 0, 1, init, {init_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 188, "init"});
  event_register ((Event){ 0, 1, init_0, {init_0_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/contact.h", 45, "init"});
  event_register ((Event){ 0, 1, init_1, {init_1_expr0}, ((int *)0), ((double *)0),
    "60d_plate_advancing_simulation.c", 97, "init"});
  event_register ((Event){ 0, 1, logfile, {logfile_expr0}, ((int *)0), ((double *)0),
    "60d_plate_advancing_simulation.c", 137, "logfile"});
  event_register ((Event){ 0, 3, dumpfile, {dumpfile_expr0, dumpfile_expr1, dumpfile_expr2}, ((int *)0), ((double *)0),
    "60d_plate_advancing_simulation.c", 146, "dumpfile"});
  event_register ((Event){ 0, 1, cleanup, {cleanup_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/run.h", 50, "cleanup"});
  event_register ((Event){ 0, 1, set_dtmax, {set_dtmax_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 214, "set_dtmax"});
  event_register ((Event){ 0, 1, stability, {stability_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 216, "stability"});
  event_register ((Event){ 0, 1, stability_0, {stability_0_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/tension.h", 36, "stability"});
  event_register ((Event){ 0, 1, stability_1, {stability_1_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/vof.h", 136, "stability"});
  event_register ((Event){ 0, 1, vof, {vof_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 226, "vof"});
  event_register ((Event){ 0, 1, vof_0, {vof_0_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/contact.h", 51, "vof"});
  event_register ((Event){ 0, 1, vof_1, {vof_1_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/vof.h", 363, "vof"});
  event_register ((Event){ 0, 1, tracer_advection, {tracer_advection_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 227, "tracer_advection"});
  event_register ((Event){ 0, 1, tracer_advection_0, {tracer_advection_0_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/two-phase.h", 64, "tracer_advection"});
  event_register ((Event){ 0, 1, tracer_diffusion, {tracer_diffusion_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 228, "tracer_diffusion"});
  event_register ((Event){ 0, 1, properties, {properties_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 235, "properties"});
  event_register ((Event){ 0, 1, properties_0, {properties_0_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/two-phase.h", 95, "properties"});
  event_register ((Event){ 0, 1, advection_term, {advection_term_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 308, "advection_term"});
  event_register ((Event){ 0, 1, viscous_term, {viscous_term_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 337, "viscous_term"});
  event_register ((Event){ 0, 1, acceleration, {acceleration_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 373, "acceleration"});
  event_register ((Event){ 0, 1, acceleration_0, {acceleration_0_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/iforce.h", 43, "acceleration"});
  event_register ((Event){ 0, 1, acceleration_1, {acceleration_1_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/tension.h", 72, "acceleration"});
  event_register ((Event){ 0, 1, acceleration_2, {acceleration_2_expr0}, ((int *)0), ((double *)0),
    "60d_plate_advancing_simulation.c", 112, "acceleration"});
  event_register ((Event){ 0, 1, projection, {projection_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 413, "projection"});
  event_register ((Event){ 0, 1, end_timestep, {end_timestep_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 428, "end_timestep"});
  event_register ((Event){ 0, 1, adapt, {adapt_expr0}, ((int *)0), ((double *)0),
    "/home/fpl/softwares/basilisk/src/navier-stokes/centered.h", 438, "adapt"});
  event_register ((Event){ 0, 1, adapt_0, {adapt_0_expr0}, ((int *)0), ((double *)0),
    "60d_plate_advancing_simulation.c", 192, "adapt"});
  void allocate_globals (int);
  allocate_globals (14);
  set_fpe();
  quadtree_methods();
  init_vector ((vector){{12},{13}}, "h");
  init_scalar ((scalar){11}, "rhov");
  init_face_vector ((vector){{9},{10}}, "alphav");
  init_scalar ((scalar){8}, "f");
  init_face_vector ((vector){{6},{7}}, "uf");
  init_scalar ((scalar){5}, "pf");
  init_vector ((vector){{3},{4}}, "g");
  init_vector ((vector){{1},{2}}, "u");
  init_scalar ((scalar){0}, "p");
  init_const_scalar ((scalar){_NVARMAX+5}, "zeroc",  0.);
  init_const_scalar ((scalar){_NVARMAX+4}, "unity",  1.);
  init_const_vector ((vector){{_NVARMAX+2},{_NVARMAX+3}}, "unityf", (double []) {1.,1.,1.});
  init_const_vector ((vector){{_NVARMAX+0},{_NVARMAX+1}}, "zerof", (double []) {0.,0.,0.});
  _set_boundary0();
  _set_boundary1();
  _set_boundary2();
  _set_boundary3();
  _set_boundary4();
  _set_boundary5();
  _set_boundary6();
  _set_boundary7();
  _set_boundary10();
  _set_boundary11();
  _set_boundary12();
  _set_boundary13();
  _set_boundary14();
  _set_boundary15();
  _set_boundary16();
  _set_boundary17();
}
