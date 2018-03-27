#include "simdocs.h"

#define OMPMINOPS       50000

static int gk_exit_on_error = 1;

GK_MKALLOC(gk_i,   int)
GK_MKALLOC(gk_i32, int32_t)
GK_MKALLOC(gk_f,   float)
GK_MKALLOC(gk_fkv, gk_fkv_t)
GK_MKALLOC(gk_dkv, gk_dkv_t)
GK_MKBLAS(gk_f, float, float)


/*************************************************************************
* This function prints an error message and exits 
**************************************************************************/
void errexit(char *f_str,...)
{
  va_list argp;

  va_start(argp, f_str);
  vfprintf(stderr, f_str, argp);
  va_end(argp);

  fprintf(stderr,"\n");
  fflush(stderr);

  if (gk_exit_on_error)
    exit(-2);
}


/*************************************************************************
* This function prints an error message and raises a signum signal
**************************************************************************/
void gk_errexit(int signum, char *f_str,...)
{
  va_list argp;

  va_start(argp, f_str);
  vfprintf(stderr, f_str, argp);
  va_end(argp);

  fprintf(stderr,"\n");
  fflush(stderr);

  if (gk_exit_on_error)
    raise(signum);
}


/*************************************************************************/
/*! This function is my wrapper around malloc that provides the following
    enhancements over malloc:
    * It always allocates one byte of memory, even if 0 bytes are requested.
      This is to ensure that checks of returned values do not lead to NULL
      due to 0 bytes requested.
    * It zeros-out the memory that is allocated. This is for a quick init
      of the underlying datastructures.
*/
/**************************************************************************/
void *gk_malloc(size_t nbytes, char *msg)
{
  void *ptr=NULL;

  if (nbytes == 0)
    nbytes++;  /* This was added to make all the mallocs to actually allocate some memory */

#ifdef USE_DLMALLOC
#ifdef GKMSPACE
  if (gk_mspace == 0)
    gk_mspace = create_mspace(0, 0);

  if (gk_mspace == NULL) {
    gk_errexit(SIGMEM, "***Memory allocation failed for creating gk_mspace.");
    return NULL;
  }

  ptr = (void *)mspace_malloc(gk_mspace, nbytes);
#else
  ptr = (void *)dlmalloc(nbytes);
#endif
#else
  ptr = (void *)malloc(nbytes);
#endif

  if (ptr == NULL) {
    fprintf(stderr, "   Current memory used:  %10"PRId64" bytes\n", (int64_t)gk_GetCurMemoryUsed());
    fprintf(stderr, "   Maximum memory used:  %10"PRId64" bytes\n", (int64_t)gk_GetMaxMemoryUsed());

    gk_errexit(SIGMEM, "***Memory allocation failed for %s. Requested size: %"PRId64" bytes", 
        msg, (int64_t)nbytes);
    return NULL;
  }

  /* zero-out the allocated space */
  memset(ptr, 0, nbytes);

  return ptr;
}


/*************************************************************************
* This function is my wrapper around realloc
**************************************************************************/
void *gk_realloc(void *oldptr, size_t nbytes, char *msg)
{
  void *ptr=NULL;

  nbytes++;  /* This was added to make all the mallocs to actually allocate some memory */

  if (nbytes == 0) {
    gk_free((void **)&oldptr, LTERM);
    return NULL;
  }

#ifdef USE_DLMALLOC
#ifdef GKMSPACE
  if (gk_mspace == 0)
    gk_mspace = create_mspace(0, 0);

  if (gk_mspace == NULL) {
    gk_errexit(SIGMEM, "***Memory allocation failed for creating gk_mspace.");
    return NULL;
  }

  ptr = (void *)mspace_realloc(gk_mspace, oldptr, nbytes);
#else
  ptr = (void *)dlrealloc(oldptr, nbytes);
#endif
#else
  ptr = (void *)realloc(oldptr, nbytes);
#endif

  if (ptr == NULL) {
    fprintf(stderr, "   Maximum memory used:              %10ju bytes\n", gk_GetMaxMemoryUsed());
    fprintf(stderr, "   Current memory used:              %10ju bytes\n", gk_GetCurMemoryUsed());

    gk_errexit(SIGMEM, "***Memory re-allocation failed for %s. Requested size: %zd bytes", msg, nbytes);
    return NULL;
  }

  return ptr;
}


/*************************************************************************
* This function is my wrapper around free, allows multiple pointers    
**************************************************************************/
void gk_free(void **ptr1,...)
{
  va_list plist;
  void **ptr;

  if (*ptr1 != NULL)
#ifdef USE_DLMALLOC
#ifdef GKMSPACE
    mspace_free(gk_mspace, *ptr1);
#else
    dlfree(*ptr1);
#endif
#else
    free(*ptr1);
#endif
  *ptr1 = NULL;

  va_start(plist, ptr1);

  while ((ptr = va_arg(plist, void **)) != LTERM) {
    if (*ptr != NULL)
#ifdef USE_DLMALLOC
#ifdef GKMSPACE
      mspace_free(gk_mspace, *ptr);
#else
      dlfree(*ptr);
#endif
#else
      free(*ptr);
#endif
    *ptr = NULL;
  }

  va_end(plist);
}            


/*************************************************************************/
/*! Sorts an array of gk_fkv_t in decreasing order */
/*************************************************************************/
void gk_fkvsortd(size_t n, gk_fkv_t *base)
{
#define fkey_gt(a, b) ((a)->key > (b)->key)
  GK_MKQSORT(gk_fkv_t, base, n, fkey_gt);
#undef fkey_gt
}


/******************************************************************************/
/*! This function puts the 'topk' largest values in the beginning of the array */
/*******************************************************************************/
int gk_dfkvkselect(size_t n, int topk, gk_fkv_t *cand)
{
  int i, j, lo, hi, mid;
  gk_fkv_t stmp;
  float pivot;

  if (n <= topk)
    return n; /* return if the array has fewer elements than we want */

  for (lo=0, hi=n-1; lo < hi;) {
    mid = lo + ((hi-lo) >> 1);

    /* select the median */
    if (cand[lo].key < cand[mid].key)
      mid = lo;
    if (cand[hi].key > cand[mid].key)
      mid = hi;
    else 
      goto jump_over;
    if (cand[lo].key < cand[mid].key)
      mid = lo;

jump_over:
    QSSWAP(cand[mid], cand[hi], stmp);
    pivot = cand[hi].key;

    /* the partitioning algorithm */
    for (i=lo-1, j=lo; j<hi; j++) {
      if (cand[j].key >= pivot) {
        i++;
        QSSWAP(cand[i], cand[j], stmp);
      }
    }
    i++;
    QSSWAP(cand[i], cand[hi], stmp);


    if (i > topk) 
      hi = i-1;
    else if (i < topk)
      lo = i+1;
    else
      break;
  }

/*
  if (cand[lo].key < cand[hi].key)
    printf("Hmm Error: %d %d %d %f %f\n", i, lo, hi, cand[lo].key, cand[hi].key);


  for (i=topk; i<n; i++) {
    for (j=0; j<topk; j++)
      if (cand[i].key > cand[j].key)
        printf("Hmm Error: %d %d %f %f %d %d\n", i, j, cand[i].key, cand[j].key, lo, hi);
  }
*/

  return topk;
}


/************************************************************************/
/*! \brief Duplicates a string

This function is a replacement for C's standard <em>strdup()</em> function.
The key differences between the two are that gk_strdup():
  - uses the dynamic memory allocation routines of \e GKlib. 
  - it correctly handles NULL input strings.

The string that is returned must be freed by gk_free().

\param orgstr is the string that will be duplicated.
\returns A pointer to the newly created string.
\sa gk_free()
*/
/*************************************************************************/
char *gk_strdup(char *orgstr)
{
  int len;
  char *str=NULL;

  if (orgstr != NULL) {
    len = strlen(orgstr)+1;
    str = gk_malloc(len*sizeof(char), "gk_strdup: str");
    strcpy(str, orgstr);
  }

  return str;
}


/*************************************************************************
* This function opens a file
**************************************************************************/
FILE *gk_fopen(char *fname, char *mode, const char *msg)
{
  FILE *fp;
  char errmsg[8192];

  fp = fopen(fname, mode);
  if (fp != NULL)
    return fp;

  sprintf(errmsg,"file: %s, mode: %s, [%s]", fname, mode, msg);
  perror(errmsg);
  errexit("Failed on gk_fopen()\n");

  return NULL;
}


/*************************************************************************
* This function closes a file
**************************************************************************/
void gk_fclose(FILE *fp)
{
  fclose(fp);
}


/*************************************************************************
* This function checks if a file exists
**************************************************************************/
int gk_fexists(char *fname)
{
  struct stat status;

  if (stat(fname, &status) == -1)
    return 0;

  return S_ISREG(status.st_mode);
}


/*************************************************************************/
/*! This function gets some basic statistics about the file. 
    \param fname is the name of the file
    \param r_nlines is the number of lines in the file. If it is NULL,
           this information is not returned.
    \param r_ntokens is the number of tokens in the file. If it is NULL,
           this information is not returned.
    \param r_max_nlntokens is the maximum number of tokens in any line
           in the file. If it is NULL this information is not returned.
    \param r_nbytes is the number of bytes in the file. If it is NULL,
           this information is not returned.
*/
/*************************************************************************/
void gk_getfilestats(char *fname, gk_idx_t *r_nlines, gk_idx_t *r_ntokens, 
        gk_idx_t *r_max_nlntokens, gk_idx_t *r_nbytes)
{
  gk_idx_t nlines=0, ntokens=0, max_nlntokens=0, nbytes=0;
  gk_idx_t oldntokens=0, intoken=0;
  char buffer[2049], *cptr;
  size_t nread;
  FILE *fpin;

  fpin = gk_fopen(fname, "r", "gk_GetFileStats");

  while (!feof(fpin)) {
    nread = fread(buffer, sizeof(char), 2048, fpin);
    nbytes += nread;

    buffer[nread] = '\0';  /* There is space for this one */
    for (cptr=buffer; *cptr!='\0'; cptr++) {
      if (*cptr == '\n') {
        nlines++;
        ntokens += intoken;
        intoken = 0;
        if (max_nlntokens < ntokens-oldntokens)
          max_nlntokens = ntokens-oldntokens;
        oldntokens = ntokens;
      }
      else if (*cptr == ' ' || *cptr == '\t') {
        ntokens += intoken;
        intoken = 0;
      }
      else {
        intoken = 1;
      }
    }
  }
  ntokens += intoken;
  if (max_nlntokens < ntokens-oldntokens)
    max_nlntokens = ntokens-oldntokens;

  gk_fclose(fpin);

  if (r_nlines != NULL)
    *r_nlines  = nlines;
  if (r_ntokens != NULL)
    *r_ntokens = ntokens;
  if (r_max_nlntokens != NULL)
    *r_max_nlntokens = max_nlntokens;
  if (r_nbytes != NULL)
    *r_nbytes  = nbytes;
}


/*************************************************************************/
/*! This function is the GKlib implementation of glibc's getline()
    function.
    \returns -1 if the EOF has been reached, otherwise it returns the 
             number of bytes read.
*/
/*************************************************************************/
gk_idx_t gk_getline(char **lineptr, size_t *n, FILE *stream)
{
#ifdef __GNUC__
  return getline(lineptr, n, stream);
#else
  size_t i;
  int ch;

  if (feof(stream))
    return -1;  

  /* Initial memory allocation if *lineptr is NULL */
  if (*lineptr == NULL || *n == 0) {
    *n = 1024;
    *lineptr = gk_malloc((*n)*sizeof(char), "gk_getline: lineptr");
  }

  /* get into the main loop */
  i = 0;
  while ((ch = getc(stream)) != EOF) {
    (*lineptr)[i++] = (char)ch;

    /* reallocate memory if reached at the end of the buffer. The +1 is for '\0' */
    if (i+1 == *n) { 
      *n = 2*(*n);
      *lineptr = gk_realloc(*lineptr, (*n)*sizeof(char), "gk_getline: lineptr");
    }
      
    if (ch == '\n')
      break;
  }
  (*lineptr)[i] = '\0';

  return (i == 0 ? -1 : i);
#endif
}


/*************************************************************************
* This function returns the CPU seconds
**************************************************************************/
double gk_WClockSeconds(void)
{
#ifdef __GNUC__
  struct timeval ctime;

  gettimeofday(&ctime, NULL);

  return (double)ctime.tv_sec + (double).000001*ctime.tv_usec;
#else
  return (double)time(NULL);
#endif
}


/*************************************************************************
* This function returns the current ammount of dynamically allocated
* memory that is used by the system
**************************************************************************/
uintmax_t gk_GetCurMemoryUsed()
{
  size_t cused=0;
#ifdef USE_DLMALLOC
  struct mallinfo meminfo;

#ifdef GKMSPACE
  if (gk_mspace != 0) {
    meminfo = mspace_mallinfo(gk_mspace);
    cused = meminfo.uordblks;
  }
#else
  meminfo = dlmallinfo();
  cused = meminfo.uordblks;
#endif
#endif

  return (uintmax_t) cused;
}


/*************************************************************************
* This function returns the maximum ammount of dynamically allocated 
* memory that was used by the system
**************************************************************************/
uintmax_t gk_GetMaxMemoryUsed()
{
  size_t mused=0;
#ifdef USE_DLMALLOC
  struct mallinfo meminfo;

#ifdef GKMSPACE
  if (gk_mspace != 0) {
    meminfo = mspace_mallinfo(gk_mspace);
    mused = meminfo.usmblks;
  }
#else
  meminfo = dlmallinfo();
  mused = meminfo.usmblks;
#endif
#endif

  return (uintmax_t) mused;
}


/*************************************************************************/
/*! Initializes the matrix 
    \param mat is the matrix to be initialized.
*/
/*************************************************************************/
void gk_csr_Init(gk_csr_t *mat)
{
  mat->nrows = mat->ncols = -1;
  mat->rowptr = mat->rowind = mat->colptr = mat->colind = NULL;
  mat->rowids = mat->colids = NULL;
  mat->rowval = mat->colval = NULL;
  mat->rnorms = mat->cnorms = NULL;
}


/*************************************************************************/
/*! Allocate memory for a CSR matrix and initializes it 
    \returns the allocated matrix. The various fields are set to NULL.
*/
/**************************************************************************/
gk_csr_t *gk_csr_Create()
{
  gk_csr_t *mat;

  mat = (gk_csr_t *)gk_malloc(sizeof(gk_csr_t), "gk_csr_Create: mat");

  gk_csr_Init(mat);

  return mat;
}


/*************************************************************************/
/*! Compacts the column-space of the matrix by removing empty columns.
    As a result of the compaction, the column numbers are renumbered. 
    The compaction operation is done in place and only affects the row-based
    representation of the matrix.
   
    \param mat the matrix whose empty columns will be removed.
*/
/**************************************************************************/
void gk_csr_CompactColumns(gk_csr_t *mat)
{
  int i, nrows, ncols, nncols;
  int *rowptr, *rowind, *collen;

  nrows  = mat->nrows;
  ncols  = mat->ncols;
  rowptr = mat->rowptr;
  rowind = mat->rowind;

  collen = gk_ismalloc(ncols, 0, "gk_csr_CompactColumns: collen");

  for (i=0; i<rowptr[nrows]; i++) 
    collen[rowind[i]]++;

  for (nncols=0, i=0; i<ncols; i++) {
    if (collen[i] > 0) 
      collen[i] = nncols++;
  }

  for (i=0; i<rowptr[nrows]; i++) 
    rowind[i] = collen[rowind[i]];

  mat->ncols = nncols;

  gk_free((void **)&collen, LTERM);
}


/*************************************************************************/
/*! Normalizes the rows/columns of the matrix to be unit 
    length.
    \param mat the matrix itself,
    \param what indicates what will be normalized and is obtained by
           specifying GK_CSR_ROW, GK_CSR_COL, GK_CSR_ROW|GK_CSR_COL. 
    \param norm indicates what norm is to normalize to, 1: 1-norm, 2: 2-norm
*/
/**************************************************************************/
void gk_csr_Normalize(gk_csr_t *mat, int what, int norm)
{
  int i, j, n;
  int *ptr;
  float *val, sum;

  if (what&GK_CSR_ROW && mat->rowval) {
    n   = mat->nrows;
    ptr = mat->rowptr;
    val = mat->rowval;

    for (i=0; i<n; i++) {
      for (sum=0.0, j=ptr[i]; j<ptr[i+1]; j++){
	if (norm == 2)
	  sum += val[j]*val[j];
	else if (norm == 1)
	  sum += val[j]; /* assume val[j] > 0 */ 
      }
      if (sum > 0) {
	if (norm == 2)
	  sum=1.0/sqrt(sum); 
	else if (norm == 1)
	  sum=1.0/sum; 
        for (j=ptr[i]; j<ptr[i+1]; j++)
          val[j] *= sum;
	
      }
    }
  }

  if (what&GK_CSR_COL && mat->colval) {
    n   = mat->ncols;
    ptr = mat->colptr;
    val = mat->colval;

    for (i=0; i<n; i++) {
      for (sum=0.0, j=ptr[i]; j<ptr[i+1]; j++)
	if (norm == 2)
	  sum += val[j]*val[j];
	else if (norm == 1)
	  sum += val[j]; 
      if (sum > 0) {
	if (norm == 2)
	  sum=1.0/sqrt(sum); 
	else if (norm == 1)
	  sum=1.0/sum; 
        for (j=ptr[i]; j<ptr[i+1]; j++)
          val[j] *= sum;
      }
    }
  }
}


/*************************************************************************/
/*! Computes the squared of the norms of the rows/columns

    \param mat the matrix itself,
    \param what is either GK_CSR_ROW or GK_CSR_COL indicating which 
           squared norms to compute.

    \note If the rowval/colval arrays are NULL, the matrix is assumed
          to be binary and the norms are computed accordingly.
*/
/**************************************************************************/
void gk_csr_ComputeSquaredNorms(gk_csr_t *mat, int what)
{
  ssize_t i;
  int n;
  int *ptr;
  float *val, *norms;

  switch (what) {
    case GK_CSR_ROW:
      n   = mat->nrows;
      ptr = mat->rowptr;
      val = mat->rowval;

      if (mat->rnorms) gk_free((void **)&mat->rnorms, LTERM);

      norms = mat->rnorms = gk_fsmalloc(n, 0, "gk_csr_ComputeSums: norms");
      break;
    case GK_CSR_COL:
      n   = mat->ncols;
      ptr = mat->colptr;
      val = mat->colval;

      if (mat->cnorms) gk_free((void **)&mat->cnorms, LTERM);

      norms = mat->cnorms = gk_fsmalloc(n, 0, "gk_csr_ComputeSums: norms");
      break;
    default:
      gk_errexit(SIGERR, "Invalid norm type of %d.\n", what);
      return;
  }

  if (val) {
    for (i=0; i<n; i++) 
      norms[i] = gk_fdot(ptr[i+1]-ptr[i], val+ptr[i], 1, val+ptr[i], 1);
  }
  else {
    for (i=0; i<n; i++) 
      norms[i] = ptr[i+1]-ptr[i];
  }
}


/**************************************************************************/
/*! Reads a CSR matrix from the supplied file and stores it the matrix's 
    forward structure.
    \param filename is the file that stores the data.
    \param format is either GK_CSR_FMT_CLUTO or GK_CSR_FMT_CSR specifying the
           type of the input format. The difference between the two formats
           is that GK_CSR_FMT_CLUTO requires the header line, whereas 
           GK_CSR_FMT_CSR does not. 
    \param readvals is either 1 or 0, indicating if the CSR file contains
           values or it does not. It only applies when GK_CSR_FMT_CSR is
           used.
    \param numbering is either 1 or 0, indicating if the numbering of the 
           indices start from 1 or 0, respectively. If they start from 1, 
           they are automatically decreamented during input so that they
           will start from 0. It only applies when GK_CSR_FMT_CSR is
           used.
    \returns the matrix that was read.
*/
/**************************************************************************/
gk_csr_t *gk_csr_Read(char *filename, int format, int readvals, int numbering)
{
  gk_idx_t i, k, nrows, ncols, nnz;
  size_t lnlen;
  int *rowptr, *rowind, ival;
  float *rowval, fval;
  char *line=NULL, *head, *tail;
  FILE *fpin;
  gk_csr_t *mat=NULL;


  if (!gk_fexists(filename)) 
    gk_errexit(SIGERR, "File %s does not exist!\n", filename);

  if (format == GK_CSR_FMT_CLUTO) {
    fpin = gk_fopen(filename, "r", "gk_csr_Read: fpin");
    if (gk_getline(&line, &lnlen, fpin) <= 0)
      gk_errexit(SIGERR, "Premature end of input file: file:%s\n", filename);
    if (sscanf(line, "%"SCNGKIDX" %"SCNGKIDX" %"SCNGKIDX, &nrows, &ncols, &nnz) != 3)
      gk_errexit(SIGERR, "Header line must contain 3 integers.\n");
    readvals = 1;
    numbering = 1;
  }
  else {
    gk_getfilestats(filename, &nrows, &nnz, NULL, NULL);

    if (readvals && nnz%2 == 1)
      gk_errexit(SIGERR, "Error: The number of numbers (%d) in the input file is not even.\n", (int)nnz);
    if (readvals)
      nnz = nnz/2;
    fpin = gk_fopen(filename, "r", "gk_csr_Read: fpin");
  }

  mat = gk_csr_Create();

  mat->nrows = nrows;

  rowptr = mat->rowptr = gk_imalloc(nrows+1, "gk_csr_Read: rowptr");
  rowind = mat->rowind = gk_imalloc(nnz, "gk_csr_Read: rowind");
  rowval = mat->rowval = gk_fsmalloc(nnz, 1.0, "gk_csr_Read: rowval");

  /*----------------------------------------------------------------------
   * Read the sparse matrix file
   *---------------------------------------------------------------------*/
  numbering = (numbering ? - 1 : 0);
  for (ncols=0, rowptr[0]=0, k=0, i=0; i<nrows; i++) {
    if (gk_getline(&line, &lnlen, fpin) == -1)
      gk_errexit(SIGERR, "Premature end of input file: file while reading row %d\n", i);

    /* Parse the string and get the arguments */
    head = line;
    while (1) {
      ival = strtol(head, &tail, 0);
      if (tail == head) 
        break;
      head = tail;
      
      if ((rowind[k] = ival + numbering) < 0)
        gk_errexit(SIGERR, "Error: Invalid column number %d at row %d.\n", ival, i);

      ncols = gk_max(rowind[k], ncols);

      if (readvals) {
        fval = strtof(head, &tail);
        if (tail == head)
          gk_errexit(SIGERR, "Value could not be found for column! Row:%d, NNZ:%d\n", i, k);
        head = tail;

        rowval[k] = fval;
      }
      k++;
    }
    rowptr[i+1] = k;
  }
  mat->ncols = ncols+1;

  if (k != nnz)
    gk_errexit(SIGERR, "gk_csr_Read: Something wrong with the number of nonzeros in the input file. NNZ=%d, ActualNNZ=%d\n", nnz, k);

  gk_fclose(fpin);

  gk_free((void **)&line, LTERM);

  return mat;
}


/*************************************************************************/
/*! Frees only the memory allocated for the matrix's different fields and
    sets them to NULL.
    \param mat is the matrix whose contents will be freed.
*/    
/*************************************************************************/
void gk_csr_FreeContents(gk_csr_t *mat)
{
  gk_free((void *)&mat->rowptr, &mat->rowind, &mat->rowval, &mat->rowids,
          &mat->colptr, &mat->colind, &mat->colval, &mat->colids, 
          &mat->rnorms, &mat->cnorms,
          LTERM);
}


/*************************************************************************/
/*! Frees all the memory allocated for matrix.
    \param mat is the matrix to be freed.
*/
/*************************************************************************/
void gk_csr_Free(gk_csr_t **mat)
{
  if (*mat == NULL)
    return;
  gk_csr_FreeContents(*mat);
  gk_free((void **)mat, LTERM);
}


/*************************************************************************/
/*! Creates a row/column index from the column/row data.
    \param mat the matrix itself,
    \param what is either GK_CSR_ROW or GK_CSR_COL indicating which index
           will be created.
*/
/**************************************************************************/
void gk_csr_CreateIndex(gk_csr_t *mat, int what)
{
  /* 'f' stands for forward, 'r' stands for reverse */
  int i, j, k, nf, nr;
  int *fptr, *find, *rptr, *rind;
  float *fval, *rval;

  switch (what) {
    case GK_CSR_COL:
      nf   = mat->nrows;
      fptr = mat->rowptr;
      find = mat->rowind;
      fval = mat->rowval;

      if (mat->colptr) gk_free((void **)&mat->colptr, LTERM);
      if (mat->colind) gk_free((void **)&mat->colind, LTERM);
      if (mat->colval) gk_free((void **)&mat->colval, LTERM);

      nr   = mat->ncols;
      rptr = mat->colptr = gk_ismalloc(nr+1, 0, "gk_csr_CreateIndex: rptr");
      rind = mat->colind = gk_imalloc(fptr[nf], "gk_csr_CreateIndex: rind");
      rval = mat->colval = (fval ? gk_fmalloc(fptr[nf], "gk_csr_CreateIndex: rval") : NULL);
      break;
    case GK_CSR_ROW:
      nf   = mat->ncols;
      fptr = mat->colptr;
      find = mat->colind;
      fval = mat->colval;

      if (mat->rowptr) gk_free((void **)&mat->rowptr, LTERM);
      if (mat->rowind) gk_free((void **)&mat->rowind, LTERM);
      if (mat->rowval) gk_free((void **)&mat->rowval, LTERM);

      nr   = mat->nrows;
      rptr = mat->rowptr = gk_ismalloc(nr+1, 0, "gk_csr_CreateIndex: rptr");
      rind = mat->rowind = gk_imalloc(fptr[nf], "gk_csr_CreateIndex: rind");
      rval = mat->rowval = (fval ? gk_fmalloc(fptr[nf], "gk_csr_CreateIndex: rval") : NULL);
      break;
    default:
      gk_errexit(SIGERR, "Invalid index type of %d.\n", what);
      return;
  }


  for (i=0; i<nf; i++) {
    for (j=fptr[i]; j<fptr[i+1]; j++)
      rptr[find[j]]++;
  }
  MAKECSR(i, nr, rptr);
  
  if (rptr[nr] > 6*nr) {
    for (i=0; i<nf; i++) {
      for (j=fptr[i]; j<fptr[i+1]; j++) 
        rind[rptr[find[j]]++] = i;
    }
    SHIFTCSR(i, nr, rptr);

    if (fval) {
      for (i=0; i<nf; i++) {
        for (j=fptr[i]; j<fptr[i+1]; j++) 
          rval[rptr[find[j]]++] = fval[j];
      }
      SHIFTCSR(i, nr, rptr);
    }
  }
  else {
    if (fval) {
      for (i=0; i<nf; i++) {
        for (j=fptr[i]; j<fptr[i+1]; j++) {
          k = find[j];
          rind[rptr[k]]   = i;
          rval[rptr[k]++] = fval[j];
        }
      }
    }
    else {
      for (i=0; i<nf; i++) {
        for (j=fptr[i]; j<fptr[i+1]; j++) 
          rind[rptr[find[j]]++] = i;
      }
    }
    SHIFTCSR(i, nr, rptr);
  }
}


/*************************************************************************/
/*! Returns a submatrix containint a set of consecutive rows.
    \param mat is the original matrix.
    \param rstart is the starting row.
    \param nrows is the number of rows from rstart to extract.
    \returns the row structure of the newly created submatrix.
*/
/**************************************************************************/
gk_csr_t *gk_csr_ExtractSubmatrix(gk_csr_t *mat, int rstart, int nrows)
{
  int i;
  gk_csr_t *nmat;

  if (rstart+nrows > mat->nrows)
    return NULL;

  nmat = gk_csr_Create();

  nmat->nrows  = nrows;
  nmat->ncols  = mat->ncols;

  /* copy the row structure */
  if (mat->rowptr)
    nmat->rowptr = gk_icopy(nrows+1, mat->rowptr+rstart, 
                            gk_imalloc(nrows+1, "gk_csr_ExtractSubmatrix: rowptr"));
  for (i=nrows; i>=0; i--)
    nmat->rowptr[i] -= nmat->rowptr[0];
  ASSERT(nmat->rowptr[0] == 0);

  if (mat->rowids)
    nmat->rowids = gk_icopy(nrows, mat->rowids+rstart, 
                            gk_imalloc(nrows, "gk_csr_ExtractSubmatrix: rowids"));
  if (mat->rnorms)
    nmat->rnorms = gk_fcopy(nrows, mat->rnorms+rstart, 
                            gk_fmalloc(nrows, "gk_csr_ExtractSubmatrix: rnorms"));

  if (mat->rsums)
    nmat->rsums = gk_fcopy(nrows, mat->rsums+rstart, 
                            gk_fmalloc(nrows, "gk_csr_ExtractSubmatrix: rsums"));

  ASSERT(nmat->rowptr[nrows] == mat->rowptr[rstart+nrows]-mat->rowptr[rstart]);
  if (mat->rowind)
    nmat->rowind = gk_icopy(mat->rowptr[rstart+nrows]-mat->rowptr[rstart], 
                            mat->rowind+mat->rowptr[rstart], 
                            gk_imalloc(mat->rowptr[rstart+nrows]-mat->rowptr[rstart],
                                       "gk_csr_ExtractSubmatrix: rowind"));
  if (mat->rowval)
    nmat->rowval = gk_fcopy(mat->rowptr[rstart+nrows]-mat->rowptr[rstart], 
                            mat->rowval+mat->rowptr[rstart], 
                            gk_fmalloc(mat->rowptr[rstart+nrows]-mat->rowptr[rstart],
                                       "gk_csr_ExtractSubmatrix: rowval"));

  return nmat;
}


/*************************************************************************/
/*! Finds the n most similar rows (neighbors) to the query using cosine
    similarity.

    \param mat the matrix itself
    \param nqterms is the number of columns in the query
    \param qind is the list of query columns
    \param qval is the list of correspodning query weights
    \param simtype is the type of similarity and is one of GK_CSR_COS,
           GK_CSR_JAC, GK_CSR_MIN
    \param nsim is the maximum number of requested most similar rows.
           If -1 is provided, then everything is returned.
    \param minsim is the minimum similarity of the requested most 
           similar rows.
    \param hits is the result set. This array should be at least
           of length nsim.
    \param i_marker is an array of size equal to the number of rows
           whose values are initialized to -1. If NULL is provided
           then this array is allocated and freed internally.
    \param i_cand is an array of size equal to the number of rows.
           If NULL is provided then this array is allocated and freed 
           internally.
    \returns the number of identified most similar rows, which can be
             smaller than the requested number of nnbrs in those cases
             in which there are no sufficiently many neighbors.
*/
/**************************************************************************/
int gk_csr_GetSimilarRows(gk_csr_t *mat, int nqterms, int *qind, float *qval, 
        int simtype, int nsim, float minsim, gk_fkv_t *hits, int *i_marker,
        gk_fkv_t *i_cand)
{
  int i, ii, j, k, nrows, ncand;
  int *colptr, *colind, *marker;
  float *colval, *rnorms, mynorm, *rsums, mysum;
  gk_fkv_t *cand;

  if (nqterms == 0)
    return 0;

  nrows  = mat->nrows;
  colptr = mat->colptr;
  colind = mat->colind;
  colval = mat->colval;

  marker = (i_marker ? i_marker : gk_ismalloc(nrows, -1, "gk_csr_SimilarRows: marker"));
  cand   = (i_cand   ? i_cand   : gk_fkvmalloc(nrows, "gk_csr_SimilarRows: cand"));

  switch (simtype) {
    case GK_CSR_COS:
      for (ncand=0, ii=0; ii<nqterms; ii++) {
        i = qind[ii];
        for (j=colptr[i]; j<colptr[i+1]; j++) {
          k = colind[j];
          if (marker[k] == -1) {
            cand[ncand].val = k;
            cand[ncand].key = 0;
            marker[k]       = ncand++;
          }
          cand[marker[k]].key += colval[j]*qval[ii];
        }
      }
      break;

    case GK_CSR_JAC:
      for (mynorm=0.0, ncand=0, ii=0; ii<nqterms; ii++) {
        mynorm += qval[ii]*qval[ii];

        i = qind[ii];
        for (j=colptr[i]; j<colptr[i+1]; j++) {
          k = colind[j];
          if (marker[k] == -1) {
            cand[ncand].val = k;
            cand[ncand].key = 0;
            marker[k]       = ncand++;
          }
          cand[marker[k]].key += colval[j]*qval[ii];
        }
      }

      rnorms = mat->rnorms;

      for (i=0; i<ncand; i++)
        cand[i].key = cand[i].key/(rnorms[cand[i].val]+mynorm-cand[i].key);
      break;

    case GK_CSR_MIN:
      for (ncand=0, ii=0; ii<nqterms; ii++) {
        i = qind[ii];
        for (j=colptr[i]; j<colptr[i+1]; j++) {
          k = colind[j];
          if (marker[k] == -1) {
            cand[ncand].val = k;
            cand[ncand].key = 0;
            marker[k]       = ncand++;
          }
          cand[marker[k]].key += gk_min(colval[j], qval[ii]);
        }
      }

      rsums = mat->rsums;
      mysum = gk_fsum(nqterms, qval, 1);

      for (i=0; i<ncand; i++)
        cand[i].key = cand[i].key/(rsums[cand[i].val]+mysum-cand[i].key);
      break;

    default:
      gk_errexit(SIGERR, "Unknown similarity measure %d\n", simtype);
      return -1;
  }

  /* go and prune the hits that are bellow minsim */
  for (j=0, i=0; i<ncand; i++) {
    marker[cand[i].val] = -1;
    if (cand[i].key >= minsim) 
      cand[j++] = cand[i];
  }
  ncand = j;

  if (nsim == -1 || nsim >= ncand) {
    nsim = ncand;
  }
  else {
    nsim = gk_min(nsim, ncand);
    gk_dfkvkselect(ncand, nsim, cand);
  }

  gk_fkvsortd(nsim, cand);
  gk_fkvcopy(nsim, cand, hits);

  if (i_marker == NULL)
    gk_free((void **)&marker, LTERM);
  if (i_cand == NULL)
    gk_free((void **)&cand, LTERM);

  return nsim;
}

