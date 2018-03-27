/*!
\file  
\brief Parsing of command-line arguments
 
This file parses the command line arguments

\date   Started 11/27/09
\author George
\version\verbatim $Id: all_cmdline.c 9504 2011-03-04 00:36:28Z karypis $ \endverbatim
*/


#include "simdocs.h"


/*-------------------------------------------------------------------
 * Command-line options 
 *-------------------------------------------------------------------*/
static struct gk_option long_options[] = {
  {"nnbrs",             1,      0,      CMD_NNBRS},
  {"minsim",            1,      0,      CMD_MINSIM},
  {"verbosity",         1,      0,      CMD_VERBOSITY},

  {"help",              0,      0,      CMD_HELP},
  {0,                   0,      0,      0}
};


/*-------------------------------------------------------------------
 * Mini help
 *-------------------------------------------------------------------*/
static char helpstr[][100] =
{
" ",
"Usage: sd infstem [outfile]",
" ",
" Options",
"  -nnbrs=int",
"     Specifies the maximum number of nearest neighbors.",
"     Default value is 100.",
" ",
"  -minsim=float",
"     The minimum allowed similarity between neighbors. ",
"     Default value is .25.",
" ",
"  -verbosity=int",
"     Specifies the level of debugging information to be displayed.",
"     Default value is 0.",
" ",
"  -help",
"     Prints this message.",
""
};

 

/*************************************************************************/
/*! This is the entry point of the command-line argument parser */
/*************************************************************************/
void cmdline_parse(params_t *params, int argc, char *argv[])
{
  gk_idx_t i, j, k;
  int type=0;
  int c, option_index;

  /* print the command line */
  for (i=0; i<argc; i++)
    printf("%s ", argv[i]);
  printf("\n");

  /* initialize the params data structure */
  params->nnbrs     = 100;
  params->minsim    = 0.25;
  params->verbosity = -1;


  /* Parse the command line arguments  */
  while ((c = gk_getopt_long_only(argc, argv, "", long_options, &option_index)) != -1) {
    switch (c) {
      case CMD_NNBRS:
        if (gk_optarg) {
          if ((params->nnbrs = atoi(gk_optarg)) < 1)
            errexit("The -nnbrs must be greater than 1.\n");
        }
        break;

      case CMD_MINSIM:
        if (gk_optarg) {
          params->minsim = atof(gk_optarg);
          if (params->minsim < 0.0 )
            errexit("The -minsim must be non-negative.\n");
        }
        break;

      case CMD_VERBOSITY:
        if (gk_optarg) {
          params->verbosity = atoi(gk_optarg);
          if (params->verbosity < 0) 
            errexit("The -verbosity must be non-negative.\n");
        }
        break;

      case CMD_HELP:
        for (i=0; strlen(helpstr[i]) > 0; i++)
          printf("%s\n", helpstr[i]);
        exit(EXIT_SUCCESS);
        break;
      default:
        printf("Illegal command-line option(s)\nUse %s -help for a summary of the options.\n", argv[0]);
        exit(EXIT_FAILURE);
    }
  }

  /* Get the input/output file info */
  if (argc-gk_optind == 0) {
    printf("Missing input/output file info.\n  Use %s -help for a summary of the options.\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  params->infstem = gk_strdup(argv[gk_optind++]);
  params->outfile = (gk_optind < argc ? gk_strdup(argv[gk_optind++]) : NULL);
}


