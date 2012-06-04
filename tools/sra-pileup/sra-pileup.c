/*===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*
*/

#include "sra-pileup.vers.h"

#include <kapp/main.h>
#include <kapp/args.h>

#include <klib/text.h>
#include <klib/out.h>
#include <klib/rc.h>
#include <klib/log.h>
#include <klib/vector.h>
#include <klib/printf.h>

#include <kfs/file.h>
#include <kfs/directory.h>
#include <kfs/buffile.h>
#include <kfs/bzip.h>
#include <kfs/gzip.h>

#include <sra/srapath.h>
#include <vdb/manager.h>
#include <align/iterator.h>
#include <align/reference.h>
#include <align/manager.h>

#include <os-native.h>
#include <strtol.h>
#include <sysalloc.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

typedef uint8_t align_tab_select;
enum { primary_ats = 1, secondary_ats = 2, evidence_ats = 4 };


#define OPTION_REF     "aligned-region"
#define ALIAS_REF      "r"

#define OPTION_MINMAPQ "minmapq"
#define ALIAS_MINMAPQ  "q"

#define OPTION_OUTF    "outfile"
#define ALIAS_OUTF     "o"

#define OPTION_TABLE   "table"
#define ALIAS_TABLE    "t"

#define OPTION_DUPS    "duplicates"
#define ALIAS_DUPS     "d"

#define OPTION_MODE    "mode"
#define ALIAS_MODE     "m"

#define OPTION_NOQUAL  "noqual"
#define ALIAS_NOQUAL   "n"

#define OPTION_GZIP    "gzip"
#define ALIAS_GZIP     "g"

#define OPTION_BZIP    "bzip2"
#define ALIAS_BZIP     "b"

enum
{
    sra_pileup_samtools = 0,
    sra_pileup_counters = 1,
    sra_pileup_detect = 2
};


static const char * ref_usage[] = { "Filter by position on genome.",
                                    "Name can either be file specific name",
                                    "(ex: \"chr1\" or \"1\").",
                                    "\"from\" and \"to\" are 1-based coordinates",
                                    NULL };

static const char * minmapq_usage[] = { "Minimum mapq-value, ", 
                                        "alignments with lower mapq",
                                        "will be ignored (default=0)", NULL };

static const char * outf_usage[] = { "Output will be written to this file",
                                     "instead of std-out", NULL };

static const char * table_usage[] = { "What table to use (p/s/e)", 
                                      "p..primary alignments, ",
                                      "s..secondary alignments", 
                                      "e..evidence alignments", 
                                      "(default=p)", NULL };

static const char * dups_usage[] = { "Don't ignore dups (0/1)", NULL };

static const char * mode_usage[] = { "Output-format: 0...samtools, 1...just counters",
                                     "(default=0)", NULL };

static const char * noqual_usage[] = { "Omit qualities in output", NULL };

static const char * gzip_usage[] = { "Compress output using gzip", NULL };

static const char * bzip_usage[] = { "Compress output using bzip2", NULL };

OptDef MyOptions[] =
{
    /*name,           alias,         hfkt, usage-help,    maxcount, needs value, required */
    { OPTION_REF,     ALIAS_REF,     NULL, ref_usage,     0,        true,        false },
    { OPTION_MINMAPQ, ALIAS_MINMAPQ, NULL, minmapq_usage, 1,        true,        false },
    { OPTION_OUTF,    ALIAS_OUTF,    NULL, outf_usage,    1,        true,        false },
    { OPTION_TABLE,   ALIAS_TABLE,   NULL, table_usage,   1,        true,        false },
    { OPTION_DUPS,    ALIAS_DUPS,    NULL, dups_usage,    1,        true,        false },
    { OPTION_MODE,    ALIAS_MODE,    NULL, mode_usage,    1,        true,        false },
    { OPTION_NOQUAL,  ALIAS_NOQUAL,  NULL, noqual_usage,  1,        false,       false },
    { OPTION_GZIP,    ALIAS_GZIP,    NULL, gzip_usage,    1,        false,       false },
    { OPTION_BZIP,    ALIAS_BZIP,    NULL, bzip_usage,    1,        false,       false }
};


/* =========================================================================================== */

static rc_t get_str_option( const Args *args, const char *name, const char ** res )
{
    uint32_t count;
    rc_t rc = ArgsOptionCount( args, name, &count );
    *res = NULL;
    if ( rc != 0 )
    {
        LOGERR( klogInt, rc, "ArgsOptionCount() failed" );
    }
    else
    {
        if ( count > 0 )
        {
            rc = ArgsOptionValue( args, name, 0, res );
            if ( rc != 0 )
            {
                LOGERR( klogInt, rc, "ArgsOptionValue() failed" );
            }
        }
    }
    return rc;
}


static rc_t get_uint32_option( const Args *args, const char *name,
                               uint32_t *res, const uint32_t def )
{
    const char * s;
    rc_t rc = get_str_option( args, name, &s );
    if ( rc == 0 && s != NULL )
        *res = atoi( s );
    else
        *res = def;
    return rc;
}


static rc_t get_bool_option( const Args *args, const char *name, bool *res, const bool def )
{
    uint32_t count;
    rc_t rc = ArgsOptionCount( args, name, &count );
    if ( rc == 0 && count > 0 )
    {
        *res = true;
    }
    else
    {
        *res = def;
    }
    return rc;
}

/* =========================================================================================== */

typedef struct pileup_options
{
    bool ignore_dups;
    bool omit_qualities;
    bool gzip_output;
    bool bzip_output;
    uint32_t minmapq;
    uint32_t output_mode;
    align_tab_select tab_select;
    const char * output_file;
} pileup_options;


static rc_t get_pileup_options( Args * args, pileup_options *opts )
{
    rc_t rc = get_uint32_option( args, OPTION_MINMAPQ, &opts->minmapq, 0 );

    if ( rc == 0 )
    {
        rc = get_str_option( args, OPTION_OUTF, &opts->output_file );
    }

    if ( rc == 0 )
    {
         rc = get_uint32_option( args, OPTION_MODE, &opts->output_mode, sra_pileup_samtools );
    }

    if ( rc == 0 )
    {
        rc = get_bool_option( args, OPTION_DUPS, &opts->ignore_dups, true );
    }

    if ( rc == 0 )
    {
        rc = get_bool_option( args, OPTION_NOQUAL, &opts->omit_qualities, false );
    }

    if ( rc == 0 )
    {
        rc = get_bool_option( args, OPTION_GZIP, &opts->gzip_output, false );
    }

    if ( rc == 0 )
    {
        rc = get_bool_option( args, OPTION_BZIP, &opts->bzip_output, false );
    }

    if ( rc == 0 )
    {
        const char * table2use = NULL;
        rc = get_str_option( args, OPTION_TABLE, &table2use );
        opts->tab_select = primary_ats;
        if ( rc == 0 && table2use != NULL )
        {
            size_t l = string_size ( table2use );
            opts->tab_select = 0;
            if ( ( string_chr ( table2use, l, 'p' ) != NULL )||
                 ( string_chr ( table2use, l, 'P' ) != NULL ) )
            { opts->tab_select |= primary_ats; };

            if ( ( string_chr ( table2use, l, 's' ) != NULL )||
                 ( string_chr ( table2use, l, 'S' ) != NULL ) )
            { opts->tab_select |= secondary_ats; };
        }
    }


    return rc;
}

/* GLOBAL VARIABLES */
struct {
    KWrtWriter org_writer;
    void* org_data;
    KFile* kfile;
    uint64_t pos;
} g_out_writer = { NULL };

pileup_options g_options;

const char UsageDefaultName[] = "sra-pileup";

rc_t CC UsageSummary ( const char * progname )
{
    return KOutMsg ("\n"
                    "Usage:\n"
                    "  %s <path> [options]\n"
                    "\n", progname);
}


rc_t CC Usage ( const Args * args )
{
    const char * progname = UsageDefaultName;
    const char * fullpath = UsageDefaultName;
    rc_t rc;

    if ( args == NULL )
        rc = RC ( rcApp, rcArgv, rcAccessing, rcSelf, rcNull );
    else
        rc = ArgsProgram ( args, &fullpath, &progname );

    if ( rc )
        progname = fullpath = UsageDefaultName;

    UsageSummary ( progname );
    KOutMsg ( "Options:\n" );
    HelpOptionLine ( ALIAS_REF, OPTION_REF, "name[:from-to]", ref_usage );
    HelpOptionLine ( ALIAS_MINMAPQ, OPTION_MINMAPQ, "min. mapq", minmapq_usage );
    HelpOptionLine ( ALIAS_OUTF, OPTION_OUTF, "output-file", outf_usage );
    HelpOptionLine ( ALIAS_TABLE, OPTION_TABLE, "table", table_usage );
    HelpOptionLine ( ALIAS_DUPS, OPTION_DUPS, "duplicates", dups_usage );
    HelpOptionLine ( ALIAS_MODE, OPTION_MODE, "output-modes", mode_usage );
    HelpOptionsStandard ();
    HelpVersion ( fullpath, KAppVersion() );
    return rc;
}


/* Version  EXTERN
 *  return 4-part version code: 0xMMmmrrrr, where
 *      MM = major release
 *      mm = minor release
 *    rrrr = bug-fix release
 */
ver_t CC KAppVersion ( void )
{
    return SRA_PILEUP_VERS;
}

/***************************************
    N (0x4E)  n (0x6E)  <--> 0x0
    A (0x41)  a (0x61)  <--> 0x1
    C (0x43)  c (0x63)  <--> 0x2
    M (0x4D)  m (0x6D)  <--> 0x3
    G (0x47)  g (0x67)  <--> 0x4
    R (0x52)  r (0x72)  <--> 0x5
    S (0x53)  s (0x73)  <--> 0x6
    V (0x56)  v (0x76)  <--> 0x7
    T (0x54)  t (0x74)  <--> 0x8
    W (0x57)  w (0x77)  <--> 0x9
    Y (0x59)  y (0x79)  <--> 0xA
    H (0x48)  h (0x68)  <--> 0xB
    K (0x4B)  k (0x6B)  <--> 0xC
    D (0x44)  d (0x64)  <--> 0xD
    B (0x42)  b (0x62)  <--> 0xE
    N (0x4E)  n (0x6E)  <--> 0xF
***************************************/

INSDC_4na_bin ascii_2_4na_tab[] = 
{
/*         0x0  0x01 0x02 0x03 0x04 0x05 0x06 0x07 0x08 0x09 0x0A 0x0B 0x0C 0x0D 0x0E 0x0F */
/* 0x00 */ 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
/* 0x10 */ 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
/* 0x20 */ 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
/* 0x30 */ 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
/* 0x40 */ 0,   0x1, 0xE, 0x2, 0xD, 0,   0,   0x4, 0xB, 0,   0,   0xC, 0,   0x3, 0,   0,
/* 0x50 */ 0,   0,   0x5, 0x6, 0x8, 0,   0x7, 0x9, 0,   0xA, 0,   0,   0,   0,   0,   0,
/* 0x60 */ 0,   0x1, 0xE, 0x2, 0xD, 0,   0,   0x4, 0xB, 0,   0,   0xC, 0,   0x3, 0,   0,
/* 0x70 */ 0,   0,   0x5, 0x6, 0x8, 0,   0x7, 0x9, 0,   0xA, 0,   0,   0,   0,   0,   0,
/* 0x80 */ 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
/* 0x90 */ 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
/* 0xA0 */ 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
/* 0xB0 */ 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
/* 0xC0 */ 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
/* 0xD0 */ 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
/* 0xE0 */ 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
/* 0xF0 */ 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0
};


char _4na_2_ascii_tab[] =
{
/*  0x0  0x01 0x02 0x03 0x04 0x05 0x06 0x07 0x08 0x09 0x0A 0x0B 0x0C 0x0D 0x0E 0x0F */
    'N', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N',
    'n', 'a', 'c', 'm', 'g', 'r', 's', 'v', 't', 'w', 'y', 'h', 'k', 'd', 'b', 'n'
};


static char _4na_to_ascii( INSDC_4na_bin c, bool reverse )
{
    uint8_t mask = ( reverse ? 0x10 : 0 );
    return _4na_2_ascii_tab[ ( ( c | mask ) & 0x1F ) ];
}

/*
static char * dup_2_ascii( const INSDC_4na_bin * b, size_t len, bool reverse )
{
    char * res = malloc( len + 1 );
    if ( res != NULL )
    {
        uint32_t i;
        for ( i = 0; i < len; ++i )
            res[ i ] = _4na_to_ascii( b[ i ], reverse );
        res[ i ] = 0;
    }
    return res;
}
*/

/* =========================================================================================== */

typedef struct dyn_string
{
    char * data;
    size_t allocated;
    size_t data_len;
} dyn_string;


static rc_t allocated_dyn_string ( dyn_string *self, size_t size )
{
    rc_t rc = 0;
    self->data_len = 0;
    self->data = malloc( size );
    if ( self->data != NULL )
    {
        self->allocated = size;
    }
    else
    {
        self->allocated = 0;
        rc = RC( rcApp, rcNoTarg, rcConstructing, rcMemory, rcExhausted );
    }
    return rc;
}

static void free_dyn_string ( dyn_string *self )
{
    free( self->data );
    self->data = NULL;
    self->allocated = 0;
    self->data_len = 0;
}

static void reset_dyn_string( dyn_string *self )
{
    self->data_len = 0;
}

static rc_t expand_dyn_string( dyn_string *self, size_t new_size )
{
    rc_t rc = 0;
    if ( new_size > self->allocated )
    {
        self->data = realloc ( self->data, new_size );
        if ( self->data != NULL )
        {
            self->allocated = new_size;
        }
        else
        {
            self->allocated = 0;
            self->data_len = 0;
            rc = RC( rcApp, rcNoTarg, rcConstructing, rcMemory, rcExhausted );
        }
    }
    return rc;
}


static rc_t add_char_2_dyn_string( dyn_string *self, const char c )
{
    /* does nothing if self->data_len + 2 < self->allocated */
    rc_t rc = expand_dyn_string( self, self->data_len + 2 );
    if ( rc == 0 )
    {
        self->data[ self->data_len++ ] = c;
        self->data[ self->data_len ] = 0;
    }
    return rc;
}


static rc_t add_string_2_dyn_string( dyn_string *self, const char * s )
{
    rc_t rc;
    size_t size = string_size ( s );
    /* does nothing if self->data_len + size + 1 < self->allocated */
    rc = expand_dyn_string( self, self->data_len + size + 1 );
    if ( rc == 0 )
    {
        string_copy ( &(self->data[ self->data_len ]), self->allocated, s, size );
        self->data_len += size;
        self->data[ self->data_len ] = 0;
    }
    return rc;
}


static rc_t print_2_dyn_string( dyn_string * self, const char *fmt, ... )
{
    rc_t rc = 0;
    bool not_enough;

    do
    {
        size_t num_writ;
        va_list args;
        va_start ( args, fmt );
        rc = string_vprintf ( &(self->data[ self->data_len ]), 
                              self->allocated - ( self->data_len + 1 ),
                              &num_writ,
                              fmt,
                              args );
        va_end ( args );

        if ( rc == 0 )
        {
            self->data_len += num_writ;
            self->data[ self->data_len ] = 0;
        }
        not_enough = ( GetRCState( rc ) == rcInsufficient );
        if ( not_enough )
        {
            rc = expand_dyn_string( self, self->allocated + ( num_writ * 2 ) );
        }
    } while ( not_enough && rc == 0 );
    return rc;
}


/* =========================================================================================== */


static rc_t CC BufferedWriter ( void* self, const char* buffer, size_t bufsize, size_t* num_writ )
{
    rc_t rc = 0;

    assert( buffer != NULL );
    assert( num_writ != NULL );

    do {
        rc = KFileWrite( g_out_writer.kfile, g_out_writer.pos, buffer, bufsize, num_writ );
        if ( rc == 0 )
        {
            buffer += *num_writ;
            bufsize -= *num_writ;
            g_out_writer.pos += *num_writ;
        }
    } while ( rc == 0 && bufsize > 0 );
    return rc;
}


static rc_t set_stdout_to( bool gzip, bool bzip2, const char * filename, size_t bufsize )
{
    rc_t rc = 0;
    if ( gzip && bzip2 )
    {
        rc = RC( rcApp, rcFile, rcConstructing, rcParam, rcAmbiguous );
    }
    else
    {
        KDirectory *dir;
        rc = KDirectoryNativeDir( &dir );
        if ( rc != 0 )
            LOGERR( klogInt, rc, "KDirectoryNativeDir() failed" );
        else
        {
            KFile *of;
            rc = KDirectoryCreateFile ( dir, &of, false, 0664, kcmInit, "%s", filename );
            if ( rc == 0 )
            {
                KFile *buf;
                if ( gzip )
                {
                    KFile *gz;
                    rc = KFileMakeGzipForWrite( &gz, of );
                    if ( rc == 0 )
                    {
                        KFileRelease( of );
                        of = gz;
                    }
                }
                if ( bzip2 )
                {
                    KFile *bz;
                    rc = KFileMakeBzip2ForWrite( &bz, of );
                    if ( rc == 0 )
                    {
                        KFileRelease( of );
                        of = bz;
                    }
                }

                rc = KBufFileMakeWrite( &buf, of, false, bufsize );
                if ( rc == 0 )
                {
                    g_out_writer.kfile = buf;
                    g_out_writer.org_writer = KOutWriterGet();
                    g_out_writer.org_data = KOutDataGet();
                    rc = KOutHandlerSet( BufferedWriter, &g_out_writer );
                    if ( rc != 0 )
                        LOGERR( klogInt, rc, "KOutHandlerSet() failed" );
                }
                KFileRelease( of );
            }
            KDirectoryRelease( dir );
        }
    }
    return rc;
}


static void release_stdout_redirection( void )
{
    KFileRelease( g_out_writer.kfile );
    if( g_out_writer.org_writer != NULL )
    {
        KOutHandlerSet( g_out_writer.org_writer, g_out_writer.org_data );
    }
    g_out_writer.org_writer = NULL;
}


static rc_t CC write_to_FILE( void *f, const char *buffer, size_t bytes, size_t *num_writ )
{
    * num_writ = fwrite ( buffer, 1, bytes, f );
    if ( * num_writ != bytes )
        return RC( rcExe, rcFile, rcWriting, rcTransfer, rcIncomplete );
    return 0;
}


/* =========================================================================================== */


static int cmp_pchar( const char * a, const char * b )
{
    int res = 0;
    if ( ( a != NULL )&&( b != NULL ) )
    {
        size_t len_a = string_size( a );
        size_t len_b = string_size( b );
        res = string_cmp ( a, len_a, b, len_b, ( len_a < len_b ) ? len_b : len_a );
    }
    return res;
}


/* =========================================================================================== */


typedef struct range
{
    uint32_t start;
    uint32_t end;
} range;


static range * make_range( const uint64_t start, const uint64_t end )
{
    range *res = calloc( sizeof *res, 1 );
    if ( res != NULL )
    {
        res->start = start;
        res->end = end;
    }
    return res;
}


static int cmp_range( const range * a, const range * b )
{

    int res = ( a->start - b->start );
    if ( res == 0 )
        res = ( a->end - b->end );
    return res;
}


static bool range_overlapp( const range * a, const range * b )
{
    return ( !( ( b->end < a->start ) || ( b->start > a->end ) ) );
}


/* =========================================================================================== */

typedef struct reference_ranges
{
    BSTNode node;
    const char * name;
    Vector ranges;
} reference_ranges;


static reference_ranges * make_reference_ranges( const char *name )
{
    reference_ranges *res = calloc( sizeof *res, 1 );
    if ( res != NULL )
    {
        res->name = string_dup_measure ( name, NULL );
        VectorInit ( &res->ranges, 0, 5 );
    }
    return res;
}


static int CC cmp_range_wrapper( const void *item, const void *n )
{   return cmp_range( item, n ); }

static rc_t add_reference_range( reference_ranges * self, const uint64_t start, const uint64_t end )
{
    rc_t rc = 0;
    range *r = make_range( start, end );
    if ( r == NULL )
        rc = RC( rcApp, rcNoTarg, rcConstructing, rcMemory, rcExhausted );
    else
    {
        rc = VectorInsert ( &self->ranges, r, NULL, cmp_range_wrapper );
        if ( rc != 0 )
            free( r );
    }
    return rc;
}


#define RR_NAME  1
#define RR_START 2
#define RR_END   3


static void put_c( char *s, size_t size, size_t *dst, char c )
{
    if ( *dst < ( size - 1 ) )
        s[ *dst ] = c;
    (*dst)++;
}

static void finish_txt( char *s, size_t size, size_t *dst )
{
    if ( *dst > size )
        s[ size - 1 ] = 0;
    else
        s[ *dst ] = 0;
    *dst = 0;
}

static uint64_t finish_num( char *s, size_t size, size_t *dst )
{
    uint64_t res = 0;
    char *endp;
    finish_txt( s, size, dst );
    res = strtou64( s, &endp, 10 );
    return res;
}


/* s = refname:1000-2000 */
static void parse_definition( const char *s, char * name, size_t len,
                              uint64_t *start, uint64_t *end )
{
    size_t n = string_size( s );

    *start = 0;
    *end   = 0;
    name[ 0 ] = 0;
    if ( n > 0 )
    {
        size_t i, st, dst = 0;
        char tmp[ 32 ];
        st = RR_NAME;
        for ( i = 0; i < n; ++i )
        {
            char c = s[ i ];
            switch( st )
            {
                case RR_NAME  : if ( c == ':' )
                                {
                                    finish_txt( name, len, &dst );
                                    st = RR_START;
                                }
                                else
                                {
                                    put_c( name, len, &dst, c );
                                }
                                break;

                case RR_START : if ( c == '-' )
                                {
                                    *start = finish_num( tmp, sizeof tmp, &dst );
                                    st = RR_END;
                                }
                                else if ( ( c >= '0' )&&( c <= '9' ) )
                                {
                                    put_c( tmp, sizeof tmp, &dst, c );
                                }
                                break;

                case RR_END   : if ( ( c >= '0' )&&( c <= '9' ) )
                                {
                                    put_c( tmp, sizeof tmp, &dst, c );
                                }
                                break;
            }
        }
        switch( st )
        {
            case RR_NAME  : finish_txt( name, len, &dst );
                            break;

            case RR_START : *start = finish_num( tmp, sizeof tmp, &dst );
                            break;

            case RR_END   : *end = finish_num( tmp, sizeof tmp, &dst );
                            break;
        }
    }
}


static void CC release_range_wrapper( void * item, void * data )
{    free( item ); }


static void free_reference_ranges( reference_ranges * self )
{
    free( (void*)self->name );
    VectorWhack ( &self->ranges, release_range_wrapper, NULL );
    free( self );
}


static void check_reference_ranges( reference_ranges * self )
{
    uint32_t n = VectorLength( &self->ranges );
    uint32_t i = 0;
    range *a = NULL;
    while ( i < n )
    {
        range *b = VectorGet ( &self->ranges, i );
        bool remove = false;
        if ( a != NULL )
        {
            remove = range_overlapp( a, b );
            if ( remove )
            {
                range *r;
                a->end = b->end;
                VectorRemove ( &self->ranges, i, (void**)&r );
                free( r );
                n--;
            }
        }
        if ( !remove )
        {
            a = b;
            ++i;
        }
    }
}



/* =========================================================================================== */


static int CC reference_vs_pchar_wrapper( const void *item, const BSTNode *n )
{
    const reference_ranges * r = ( const reference_ranges * )n;
    return cmp_pchar( (const char *)item, r->name );
}

static reference_ranges * find_reference_ranges( BSTree * tree, const char * name )
{
    return ( reference_ranges * ) BSTreeFind ( tree, name, reference_vs_pchar_wrapper );
}

static int CC ref_vs_ref_wrapper( const BSTNode *item, const BSTNode *n )
{
   const reference_ranges * a = ( const reference_ranges * )item;
   const reference_ranges * b = ( const reference_ranges * )n;
   return cmp_pchar( a->name, b->name );
}

static rc_t add_refrange( BSTree * tree, const char * name, const uint64_t start, const uint64_t end )
{
    rc_t rc;

    reference_ranges * ref = find_reference_ranges( tree, name );
    if ( ref == NULL )
    {
        ref = make_reference_ranges( name );
        if ( ref == NULL )
            rc = RC( rcApp, rcNoTarg, rcConstructing, rcMemory, rcExhausted );
        else
            rc = add_reference_range( ref, start, end );
        if ( rc == 0 )
            rc = BSTreeInsert ( tree, (BSTNode *)ref, ref_vs_ref_wrapper );
        if ( rc != 0 )
            free_reference_ranges( ref );
    }
    else
    {
        rc = add_reference_range( ref, start, end );
    }
    return rc;
}


static rc_t parse_and_add( BSTree * tree, const char * s )
{
    uint64_t start, end;
    char name[ 64 ];
    parse_definition( s, name, sizeof name, &start, &end );
    if ( name[ 0 ] == 0 )
        return RC( rcApp, rcNoTarg, rcConstructing, rcMemory, rcExhausted );
    else
        return add_refrange( tree, name, start, end );
}


static void CC check_refrange_wrapper( BSTNode *n, void *data )
{   check_reference_ranges( ( reference_ranges * ) n ); }

static void check_refranges( BSTree * tree )
{   BSTreeForEach ( tree, false, check_refrange_wrapper, NULL ); }

static void CC release_ref_wrapper( BSTNode *n, void * data )
{    free_reference_ranges( ( reference_ranges * ) n ); }

static void free_refranges( BSTree * tree )
{    BSTreeWhack ( tree, release_ref_wrapper, NULL ); }

static rc_t init_refranges( BSTree * tree, Args * args )
{
    uint32_t count;
    rc_t rc;

    BSTreeInit( tree );
    rc = ArgsOptionCount( args, OPTION_REF, &count );
    if ( rc != 0 )
        LOGERR( klogInt, rc, "ArgsOptionCount() failed" );
    else
    {
        uint32_t i;
        for ( i = 0; i < count && rc == 0; ++i )
        {
            const char * s;
            rc = ArgsOptionValue( args, OPTION_REF, i, &s );
            if ( rc != 0 )
                LOGERR( klogInt, rc, "ArgsOptionValue() failed" );
            else
                rc = parse_and_add( tree, s );
        }
    }
    return rc;
}


static void CC count_refrange_wrapper( BSTNode *n, void *data )
{   
    reference_ranges * rr = ( reference_ranges * ) n;
    uint32_t * count = ( uint32_t * ) data;
    *count += VectorLength( &(rr->ranges) );
}

static uint32_t count_refranges( BSTree * tree )
{
    uint32_t res = 0;
    BSTreeForEach ( tree, false, count_refrange_wrapper, &res );
    return res;
}


typedef struct foreach_refrange_func
{
    rc_t ( CC * on_range ) ( const char * name, uint32_t start, uint32_t end, void *data );
    const char * name;
    void * data;
    rc_t rc;
} foreach_refrange_func;


static void CC foreach_range_vector_wrapper( void *item, void *data )
{
    range * r = ( range * ) item;
    foreach_refrange_func * func = ( foreach_refrange_func * )data;

    if ( func->rc == 0 )
    {
        func->rc = func->on_range( func->name, r->start, r->end, func->data );
    }
}


static void CC foreach_refrange_tree_wrapper( BSTNode *n, void *data )
{   
    reference_ranges * rr = ( reference_ranges * ) n;
    foreach_refrange_func * func = ( foreach_refrange_func * )data;

    if ( func->rc == 0 )
    {
        func->name = rr->name;
        VectorForEach ( &(rr->ranges), false, foreach_range_vector_wrapper, data );
    }
}


static rc_t foreach_refrange( BSTree * tree,
    rc_t ( CC * on_range ) ( const char * name, uint32_t start, uint32_t end, void *data ), 
    void *data )
{
    foreach_refrange_func func;

    func.on_range = on_range;
    func.data = data;
    func.rc = 0;
    BSTreeForEach ( tree, false, foreach_refrange_tree_wrapper, &func );
    return func.rc;
}


/* =========================================================================================== */

typedef struct tool_rec tool_rec;
struct tool_rec
{
    /* orientation towards reference ( false...in ref-orientation / true...reverse) */
    bool reverse;
    /* ptr to quality... */
    uint8_t * quality;
};


static rc_t read_base_and_len( struct VCursor const *curs,
                               const char * name,
                               int64_t row_id,
                               const void ** base,
                               uint32_t * len )
{
    uint32_t column_idx;
    rc_t rc = VCursorGetColumnIdx ( curs, &column_idx, name );
    if ( rc != 0 )
        LOGERR( klogInt, rc, "VCursorGetColumnIdx() failed" );
    else
    {
        uint32_t elem_bits, boff, len_intern;
        const void * ptr;
        rc = VCursorCellDataDirect ( curs, row_id, column_idx, 
                                     &elem_bits, &ptr, &boff, &len_intern );
        if ( rc != 0 )
            LOGERR( klogInt, rc, "VCursorCellDataDirect() failed" );
        else
        {
            if ( len != NULL ) *len = len_intern;
            if ( base != NULL ) *base = ptr;
        }
    }
    return rc;
}


static rc_t CC populate_tooldata( void *obj, const PlacementRecord *placement,
        struct VCursor const *curs,
        INSDC_coord_zero ref_window_start, INSDC_coord_len ref_window_len,
        void *data )
{
    tool_rec * rec = ( tool_rec * ) obj;
    rc_t rc = 0;

    rec->quality = NULL;
    if ( g_options.ignore_dups > 0 )
    {
        const uint32_t * samflags;
        uint32_t samflags_len;
        rc = read_base_and_len( curs, "SAM_FLAGS", placement->id, (const void **)&samflags, &samflags_len );
        if ( rc == 0 )
        {
            if ( ( *samflags & 1024 ) == 1024 )
            {
                rc = RC( rcAlign, rcType, rcAccessing, rcId, rcIgnored );
            }
        }
    }

    if ( rc == 0 )
    {
        const bool * orientation;
        rc = read_base_and_len( curs, "REF_ORIENTATION", placement->id, (const void **)&orientation, NULL );
        if ( rc == 0 )
        {
            rec->reverse = *orientation;
        }
    }

    if ( rc == 0 && g_options.omit_qualities == 0 )
    {
        const uint8_t * quality;
        uint32_t quality_len;

        rc = read_base_and_len( curs, "CLIPPED_QUALITY", placement->id, (const void **)&quality, &quality_len );
        if ( rc == 0 )
        {
            rec->quality = ( uint8_t * )rec;
            rec->quality += sizeof ( * rec );
            memcpy( rec->quality, quality, quality_len );
        }
    }
    return rc;
}


static size_t CC alloc_size( struct VCursor const *curs, int64_t row_id, void *data )
{
    tool_rec * rec;
    size_t res = ( sizeof *rec );

    if ( g_options.omit_qualities == 0 )
    {
        uint32_t q_len;
        rc_t rc = read_base_and_len( curs, "CLIPPED_QUALITY", row_id, NULL, &q_len );
        if ( rc == 0 )
        {
            res += q_len;
        }
    }
    return res;
}


static rc_t walk_ref_position( ReferenceIterator *ref_iter,
                               const PlacementRecord *rec,
                               dyn_string *line,
                               char * qual )
{
    rc_t rc = 0;
    INSDC_coord_zero seq_pos;
    int32_t state = ReferenceIteratorState ( ref_iter, &seq_pos );
    tool_rec *xrec = ( tool_rec * ) PlacementRecordCast ( rec, placementRecordExtension1 );
    bool reverse = xrec->reverse;

    if ( g_options.omit_qualities == 0 )
    {
        *qual = xrec->quality[ seq_pos ];
    }

    if ( ( state & align_iter_invalid ) == align_iter_invalid )
    {
        return add_char_2_dyn_string( line, '?' );
    }

    if ( ( state & align_iter_first ) == align_iter_first )
    {
        char s[ 3 ];
        int32_t c = rec->mapq + 33;
        if ( c > '~' ) { c = '~'; }
        if ( c < 33 ) { c = 33; }
        s[ 0 ] = '^';
        s[ 1 ] = c;
        s[ 2 ] = 0;
        rc = add_string_2_dyn_string( line, s );
    }

    if ( ( state & align_iter_skip ) == align_iter_skip )
    {
        if ( rc == 0 )
        {
            rc = add_char_2_dyn_string( line, '*' );
        }
        if ( g_options.omit_qualities == 0 )
        {
            *qual = xrec->quality[ seq_pos + 1 ];
        }
    }
    else
    {
        if ( rc == 0 )
        {
            if ( ( state & align_iter_match ) == align_iter_match )
            {
                rc = add_char_2_dyn_string( line, ( reverse ? ',' : '.' ) );
            }
            else
            {
    /*                **ptr = _4na_to_ascii( state & 0x0F, reverse ); */
                rc = add_char_2_dyn_string( line, ( reverse ? tolower( state & 0xFF ) : ( state & 0xFF ) ) );
            }
        }
    }

    if ( ( state & align_iter_insert ) == align_iter_insert )
    {
        const INSDC_4na_bin *bases;
        uint32_t i;
        uint32_t n = ReferenceIteratorBasesInserted ( ref_iter, &bases );
        
        rc = print_2_dyn_string( line, "+%u", n );
        for ( i = 0; i < n && rc == 0; ++i )
        {
/*                    (*ptr)[ i ] = _4na_to_ascii( bases[ i ], reverse ); */
            rc = add_char_2_dyn_string( line, ( reverse ? tolower( bases[ i ] ) : bases[ i ] ) );
        }
    }

    if ( ( state & align_iter_delete ) == align_iter_delete )
    {
        const INSDC_4na_bin *bases;
        INSDC_coord_zero ref_pos;
        uint32_t n = ReferenceIteratorBasesDeleted ( ref_iter, &ref_pos, &bases );
        if ( bases != NULL )
        {
            uint32_t i;
            rc = print_2_dyn_string( line, "-%u", n );
            for ( i = 0; i < n && rc == 0; ++i )
            {
                rc = add_char_2_dyn_string( line, _4na_to_ascii( bases[ i ], reverse ) );
            }
            free( (void *) bases );
        }
    }

    if ( ( ( state & align_iter_last ) == align_iter_last )&& ( rc == 0 ) )
    {
        rc = add_char_2_dyn_string( line, '$' );
    }
    return rc;
}


static rc_t walk_position( ReferenceIterator *ref_iter,
                           const char * refname,
                           dyn_string *line,
                           dyn_string *qualities,
                           bool skip_empty )
{
    INSDC_coord_zero pos;
    uint32_t depth;
    INSDC_4na_bin base;

    rc_t rc = ReferenceIteratorPosition ( ref_iter, &pos, &depth, &base );
    if ( rc != 0 )
    {
        if ( GetRCState( rc ) != rcDone )
        {
            LOGERR( klogInt, rc, "ReferenceIteratorNextPos() failed" );
        }
    }
    else if ( ( depth > 0 )||( !skip_empty ) )
    {
        rc = expand_dyn_string( line, ( 5 * depth ) + 100 );
        if ( rc == 0 )
        {
            rc = expand_dyn_string( qualities, depth + 100 );
            if ( rc == 0 )
            {
                rc_t rc1 = 0;
                char c = _4na_to_ascii( base, false );

                reset_dyn_string( line );
                rc = print_2_dyn_string( line, "%s\t%u\t%c\t%u\t", refname, pos + 1, c, depth );
                if ( rc == 0 )
                {
                    if ( depth > 0 )
                    {
                        const PlacementRecord *rec;
                        rc1 = ReferenceIteratorNextPlacement ( ref_iter, &rec );
                        if ( rc1 == 0 )
                        {
                            uint32_t i = 0;

                            /* this is the 3rd level of walking the reference-iterator: 
                               visiting each aligned base on this reference-position */
                            while ( rc1 == 0 )
                            {
                                rc1 = walk_ref_position( ref_iter, rec, line, &( qualities->data[ i++ ] ) );
                                if ( rc1 == 0 )
                                {
                                    rc1 = ReferenceIteratorNextPlacement ( ref_iter, &rec );
                                }
                            }

                            if ( g_options.omit_qualities == 0 )
                            {
                                rc1 = add_char_2_dyn_string( line, '\t' );
                                for ( i = 0; i < depth && rc1 == 0; ++i )
                                {
                                    rc1 = add_char_2_dyn_string( line, qualities->data[ i ] + 33 );
                                }
                            }
                        }
                    }
                    if ( GetRCState( rc1 ) == rcDone ) { rc = 0; } else { rc = rc1; }
                    if ( rc == 0 )
                    {
                        /* only one OUTMSG(()) per line... */
                        OUTMSG(( "%s\n", line->data ));
                    }
                }
            }
        }
    } 
    return rc;
}


static rc_t walk_reference( ReferenceIterator *ref_iter,
                            const char * refname,
                            bool skip_empty )
{
    dyn_string line;
    rc_t rc = allocated_dyn_string ( &line, 4096 );
    if ( rc == 0 )
    {
        dyn_string qualities;
        rc = allocated_dyn_string ( &qualities, 4096 );
        if ( rc == 0 )
        {
            while ( rc == 0 )
            {
                rc = Quitting ();
                if ( rc == 0 )
                {
                    /* this is the 2nd level of walking the reference-iterator: 
                       visiting each position (that has alignments) on this reference */
                    rc = ReferenceIteratorNextPos ( ref_iter, skip_empty );
                    if ( rc != 0 )
                    {
                        if ( GetRCState( rc ) != rcDone )
                            LOGERR( klogInt, rc, "ReferenceIteratorNextPos() failed" );
                    }
                    else
                    {
                        rc = walk_position( ref_iter, refname, &line, &qualities, skip_empty );
                        if ( GetRCState( rc ) == rcDone ) { rc = 0; }
                    }
                }
            }
            free_dyn_string ( &qualities );
        }
        free_dyn_string ( &line );
    }
    return rc;
}


static rc_t walk_counter_position( ReferenceIterator *ref_iter,
                           const char * refname,
                           bool skip_empty )
{
    INSDC_coord_zero pos;
    uint32_t depth;
    INSDC_4na_bin base;

    rc_t rc = ReferenceIteratorPosition ( ref_iter, &pos, &depth, &base );
    if ( rc != 0 )
    {
        if ( GetRCState( rc ) != rcDone )
        {
            LOGERR( klogInt, rc, "ReferenceIteratorNextPos() failed" );
        }
    }
    else if ( ( depth > 0 )||( !skip_empty ) )
    {
        rc_t rc1 = 0;
        const PlacementRecord *rec;
        char c = _4na_to_ascii( base, false );
        uint32_t matches = 0;
        OUTMSG(( "%s\t%u\t%c\t%u", refname, pos + 1, c, depth ));

        rc1 = ReferenceIteratorNextPlacement ( ref_iter, &rec );
        while ( rc1 == 0 )
        {
            INSDC_coord_zero seq_pos;
            int32_t state = ReferenceIteratorState ( ref_iter, &seq_pos );
            if ( ( state & align_iter_match ) == align_iter_match )
            {
                matches++;
            }
            rc1 = ReferenceIteratorNextPlacement ( ref_iter, &rec );
        }

        if ( GetRCState( rc1 ) == rcDone ) { rc = 0; } else { rc = rc1; }
        if ( rc == 0 )
        {
            OUTMSG(( "\t%u\n", matches ));
        }
    } 
    return rc;
}


static rc_t walk_just_counters( ReferenceIterator *ref_iter,
                                const char * refname,
                                bool skip_empty )
{
    rc_t rc = 0;
    while ( rc == 0 )
    {
        rc = Quitting ();
        if ( rc == 0 )
        {
            /* this is the 2nd level of walking the reference-iterator: 
               visiting each position (that has alignments) on this reference */
            rc = ReferenceIteratorNextPos ( ref_iter, skip_empty );
            if ( rc != 0 )
            {
                if ( GetRCState( rc ) != rcDone )
                    LOGERR( klogInt, rc, "ReferenceIteratorNextPos() failed" );
            }
            else
            {
                rc = walk_counter_position( ref_iter, refname, skip_empty );
                if ( GetRCState( rc ) == rcDone ) { rc = 0; }
            }
        }
    }
    return 0;
}


enum { fsm_INIT = 0, fsm_DATA, fsm_GAP1, fsm_GAP2 };

typedef struct fsm_context fsm_context;
struct fsm_context
{
    uint32_t state;
    const char * refname;
    INSDC_coord_zero start;
    INSDC_coord_zero end;
};


static void fsm_initialize( fsm_context * ctx, const char * refname )
{
    ctx->state = fsm_INIT;
    ctx->refname = refname;
    ctx->start = 0;
    ctx->end = 0;
}

static void fsm_finalize( fsm_context * ctx, INSDC_coord_zero pos )
{
    switch( ctx->state )
    {
        case fsm_DATA : ;
        case fsm_GAP1 : OUTMSG(( "%s:%u-%u\n", ctx->refname, ctx->start, pos )); break;
    }
}

/* transition into state 'fsm_DATA' */
static void fsm_data( fsm_context * ctx, INSDC_coord_zero pos )
{
    switch( ctx->state )
    {
        case fsm_INIT : ;
        case fsm_GAP2 : ctx->start = pos; break;
    }
    ctx->state = fsm_DATA;
}

/* transition into state 'fsm_GAP1' */
static void fsm_gap1( fsm_context * ctx, INSDC_coord_zero pos )
{
    if ( ctx->state == fsm_DATA )
    {
        ctx->end = pos;
    }
    ctx->state = fsm_GAP1;
}

/* transition into state 'fsm_GAP2' */
static void fsm_gap2( fsm_context * ctx )
{
    if ( ctx->state == fsm_GAP1 )
    {
        OUTMSG(( "%s:%u-%u\n", ctx->refname, ctx->start, ctx->end ));
    }
    ctx->state = fsm_GAP2;
}


static void fsm_run( fsm_context * ctx, uint32_t depth, INSDC_coord_zero pos, uint32_t maxgap )
{
    switch( ctx->state )
    {
        case fsm_INIT : if ( depth > 0 )
                            fsm_data( ctx, pos );
                        else
                            fsm_gap2( ctx );
                        break;

        case fsm_DATA : if ( depth == 0 )
                            fsm_gap1( ctx, pos );
                        break;

        case fsm_GAP1 : if ( ( pos - ctx->end ) > maxgap )
                            fsm_gap2( ctx );
                        if ( depth > 0 )
                            fsm_data( ctx, pos );
                        break;

        case fsm_GAP2 : if ( depth > 0 )
                            fsm_data( ctx, pos );
                        break;
    }

}

#if 0
static void fsm_show( fsm_context * ctx, INSDC_coord_zero pos )
{
    switch( ctx->state )
    {
        case fsm_INIT : OUTMSG(( "[%u].INIT\n", pos )); break;
        case fsm_DATA : OUTMSG(( "[%u].DATA\n", pos )); break;
        case fsm_GAP1 : OUTMSG(( "[%u].GAP1\n", pos )); break;
        case fsm_GAP2 : OUTMSG(( "[%u].GAP2\n", pos )); break;
    }
}
#endif

static rc_t walk_and_detect( ReferenceIterator *ref_iter,
                             const char * refname, uint32_t maxgap )
{
    rc_t rc = 0;
    INSDC_coord_zero pos;

    /* here we have a little FSM with these states: INIT, DATA, GAP1, GAP2 */
    fsm_context ctx;
    fsm_initialize( &ctx, refname );
    while ( rc == 0 )
    {
        rc = Quitting ();
        if ( rc == 0 )
        {
            /* this is the 2nd level of walking the reference-iterator: 
               visiting each position (that has alignments) on this reference */
            rc = ReferenceIteratorNextPos ( ref_iter, true );
            if ( rc != 0 )
            {
                if ( GetRCState( rc ) != rcDone )
                    LOGERR( klogInt, rc, "ReferenceIteratorNextPos() failed" );
            }
            else
            {
                uint32_t depth;
                rc = ReferenceIteratorPosition ( ref_iter, &pos, &depth, NULL );
                if ( rc == 0 )
                {
                    fsm_run( &ctx, depth, pos, maxgap );
                    /* fsm_show( &ctx, pos ); */
                }
                if ( GetRCState( rc ) == rcDone ) { rc = 0; }
            }
        }
    }
    fsm_finalize( &ctx, pos );
    return 0;
}


static rc_t walk_ref_iter( ReferenceIterator *ref_iter, bool skip_empty )
{
    rc_t rc = 0;
    while( rc == 0 )
    {
        /* this is the 1st level of walking the reference-iterator: 
           visiting each (requested) reference */
        struct ReferenceObj const * refobj;
        rc = ReferenceIteratorNextReference( ref_iter, &refobj );
        if ( rc == 0 )
        {
            const char * refname = NULL;
            rc = ReferenceObj_Name( refobj, &refname );
            if ( rc == 0 )
            {
                switch( g_options.output_mode )
                {
                    case sra_pileup_samtools : rc = walk_reference( ref_iter, refname, skip_empty );
                                               break;

                    case sra_pileup_counters : rc = walk_just_counters( ref_iter, refname, skip_empty );
                                               break;

                    case sra_pileup_detect : rc = walk_and_detect( ref_iter, refname, 200 );
                                             break;

                    default : OUTMSG(( "unknown output-mode '%u'\n", g_options.output_mode ));
                              break;

                }
                if ( GetRCState( rc ) == rcDone ) { rc = 0; }
            }
            else
                LOGERR( klogInt, rc, "ReferenceObj_Name() failed" );
        }
        else
        {
            if ( GetRCState( rc ) != rcDone )
                LOGERR( klogInt, rc, "ReferenceIteratorNextReference() failed" );
        }
    }
    if ( GetRCState( rc ) == rcDone ) { rc = 0; }
    return rc;
}

/* =========================================================================================== */


typedef struct prepare_ctx
{
    ReferenceIterator *ref_iter;
    const VDatabase *db;
    const ReferenceList *reflist;
    uint32_t minmapq;
    align_tab_select select;
} prepare_ctx;


static rc_t add_quality_and_orientation( const VTable *tbl, const VCursor ** cursor )
{
    rc_t rc = VTableCreateCursorRead ( tbl, cursor );
    if ( rc != 0 )
        LOGERR( klogInt, rc, "VTableCreateCursorRead() failed" );

    if ( rc == 0 && g_options.omit_qualities == 0 )
    {
        uint32_t quality_idx;
        rc = VCursorAddColumn ( *cursor, &quality_idx, "CLIPPED_QUALITY" );
        if ( rc != 0 )
            LOGERR( klogInt, rc, "VCursorAddColumn(QUALITY) failed" );
    }

    if ( rc == 0 )
    {
        uint32_t ref_orientation_idx;
        rc = VCursorAddColumn ( *cursor, &ref_orientation_idx, "REF_ORIENTATION" );
        if ( rc != 0 )
            LOGERR( klogInt, rc, "VCursorAddColumn(REF_ORIENTATION) failed" );
    }

    if ( rc == 0 )
    {
        uint32_t sam_flags_idx;
        rc = VCursorAddColumn ( *cursor, &sam_flags_idx, "SAM_FLAGS" );
        if ( rc != 0 )
            LOGERR( klogInt, rc, "VCursorAddColumn(SAM_FLAGS) failed" );
    }
    return rc;
}


static rc_t prepare_prim_cursor( const VDatabase *db, const VCursor ** cursor )
{
    const VTable *tbl;
    rc_t rc = VDatabaseOpenTableRead ( db, &tbl, "PRIMARY_ALIGNMENT" );
    if ( rc != 0 )
        LOGERR( klogInt, rc, "VDatabaseOpenTableRead(PRIMARY_ALIGNMENT) failed" );
    else
    {
        rc = add_quality_and_orientation( tbl, cursor );
        VTableRelease ( tbl );
    }
    return rc;
}


static rc_t prepare_sec_cursor( const VDatabase *db, const VCursor ** cursor )
{
    const VTable *tbl;
    rc_t rc = VDatabaseOpenTableRead ( db, &tbl, "SECONDARY_ALIGNMENT" );
    if ( rc != 0 )
        LOGERR( klogInt, rc, "VDatabaseOpenTableRead(SECONDARY_ALIGNMENT) failed" );
    else
    {
        rc = add_quality_and_orientation( tbl, cursor );
        VTableRelease ( tbl );
    }
    return rc;
}


static rc_t prepare_section( prepare_ctx * ctx, const ReferenceObj *refobj,
                             const char * name, uint32_t start, uint32_t end )
{
    INSDC_coord_len len;
    rc_t rc = ReferenceObj_SeqLength( refobj, &len );
    if ( rc != 0 )
        LOGERR( klogInt, rc, "ReferenceObj_SeqLength() failed" );
    else
    {
        if ( start == 0 ) start = 1;
        if ( ( end == 0 )||( end > len + 1 ) )
        {
            end = ( len - start ) + 1;
        }
        /* depending on ctx->select prepare primary, secondary or both... */
        if ( ( ctx->select & primary_ats ) == primary_ats )
        {
            const VCursor * align_cursor = NULL;
            rc = prepare_prim_cursor( ctx->db, &align_cursor );
            if ( rc == 0 )
            {
/*                OUTMSG(( "prepare primary: <%s> %u..%u\n", name, start, end )); */
                rc = ReferenceIteratorAddPlacements ( ctx->ref_iter,        /* the outer ref-iter */
                                                      refobj,               /* the ref-obj for this chromosome */
                                                      start - 1,            /* start ( zero-based ) */
                                                      end - start + 1,      /* length */
                                                      NULL,                 /* ref-cursor */
                                                      align_cursor,         /* align-cursor */
                                                      primary_align_ids );  /* which id's */
                if ( rc != 0 )
                    LOGERR( klogInt, rc, "ReferenceIteratorAddPlacements() failed" );
                VCursorRelease( align_cursor );
            }
        }
        if ( ( rc == 0 ) && ( ( ctx->select & secondary_ats ) == secondary_ats ) )
        {
            const VCursor * align_cursor = NULL;
            rc = prepare_sec_cursor( ctx->db, &align_cursor );
            if ( rc == 0 )
            {
/*                OUTMSG(( "prepare secondary: <%s> %u..%u\n", name, start, end )); */
                rc = ReferenceIteratorAddPlacements ( ctx->ref_iter,        /* the outer ref-iter */
                                                      refobj,               /* the ref-obj for this chromosome */
                                                      start - 1,            /* start ( zero-based ) */
                                                      end - start + 1,      /* length */
                                                      NULL,                 /* ref-cursor */
                                                      align_cursor,         /* align-cursor */
                                                      secondary_align_ids );  /* which id's */
                if ( rc != 0 )
                    LOGERR( klogInt, rc, "ReferenceIteratorAddPlacements() failed" );
            }
            VCursorRelease( align_cursor );
        }
    }
    return rc;
}


static rc_t prepare_whole_file( prepare_ctx * ctx )
{
    uint32_t count;
    rc_t rc = ReferenceList_Count( ctx->reflist, &count );
    if ( rc != 0 )
        LOGERR( klogInt, rc, "ReferenceList_Count() failed" );
    else
    {
        uint32_t idx;
        for ( idx = 0; idx < count && rc == 0; ++idx )
        {
            const ReferenceObj *refobj;
            rc = ReferenceList_Get( ctx->reflist, &refobj, idx );
            if ( rc != 0 )
                LOGERR( klogInt, rc, "ReferenceList_Get() failed" );
            else
            {
                const char *name;
                rc = ReferenceObj_Name( refobj, &name );
                if ( rc != 0 )
                    LOGERR( klogInt, rc, "ReferenceObj_Name() failed" );
                else
                {
                    rc = prepare_section( ctx, refobj, name, 0, 0 );
                }
            }
        }
    }

    return rc;
}


static rc_t CC prepare_range( const char * name, uint32_t start, uint32_t end, void *data )
{
    prepare_ctx * ctx = ( prepare_ctx * )data;
    const ReferenceObj* refobj;
    rc_t rc = ReferenceList_Find( ctx->reflist, &refobj, name, string_size( name ) );
    if ( rc != 0 )
        LOGERR( klogInt, rc, "ReferenceList_Find() failed" );
    else
    {
        rc = prepare_section( ctx, refobj, name, start, end );
    }
    return rc;
}


static rc_t prepare_pileup_file( Args * args,
                                 ReferenceIterator *ref_iter,
                                 const VDBManager *vdb_mgr,
                                 const char * filename,
                                 BSTree * tree )
{
    const VDatabase *db;
    rc_t rc = VDBManagerOpenDBRead ( vdb_mgr, &db, NULL, "%s", filename );
    if ( rc != 0 )
        LOGERR( klogInt, rc, "VDBManagerOpenDBRead() failed" );
    else
    {
        const ReferenceList *reflist;
        rc = ReferenceList_MakeDatabase( &reflist, db, ereferencelist_4na | ereferencelist_usePrimaryIds,
                                         1024, NULL, 0 );
        if ( rc != 0 )
            LOGERR( klogInt, rc, "ReferenceList_MakeDatabase() failed" );
        else
        {
            prepare_ctx ctx;

            ctx.ref_iter = ref_iter;
            ctx.db = db;
            ctx.reflist = reflist;
            ctx.minmapq = g_options.minmapq;
            ctx.select = g_options.tab_select;
            if ( count_refranges( tree ) == 0 )
            {
                /* the user has not specified a reference-range : use the whole file... */
                rc = prepare_whole_file( &ctx );
            }
            else
            {
                /* pick only the requested ranges... */
                rc = foreach_refrange( tree, prepare_range, &ctx );
            }
            ReferenceList_Release( reflist );
        }
        VDatabaseRelease ( db );
    }
    return rc;
}


static bool is_this_a_filesystem_path( const char * path )
{
    bool res = false;
    size_t i, n = string_size ( path );
    for ( i = 0; i < n && !res; ++i )
    {
        char c = path[ i ];
        res = ( c == '.' || c == '/' || c == '\\' );
    }
    return res;
}


static char *translate_accession( SRAPath *my_sra_path,
                           const char *accession,
                           const size_t bufsize )
{
    rc_t rc;
    char * res = calloc( 1, bufsize );
    if ( res == NULL ) return NULL;

    rc = SRAPathFind( my_sra_path, accession, res, bufsize );
    if ( GetRCState( rc ) == rcNotFound )
    {
        free( res );
        return NULL;
    }
    else if ( GetRCState( rc ) == rcInsufficient )
    {
        /* bufsize was insufficient ---> repeat */
        free( res );
        return translate_accession( my_sra_path, accession, bufsize * 2 );
    }
    else if ( rc != 0 )
    {
        free( res );
        return NULL;
    }
    return res;
}


static rc_t resolve_accession( const KDirectory *my_dir, char ** path )
{
    SRAPath *my_sra_path;
    rc_t rc = 0;

    if ( strchr ( *path, '/' ) != NULL )
        return 0;

    rc = SRAPathMake( &my_sra_path, my_dir );
    if ( rc != 0 )
    {
        if ( GetRCState ( rc ) != rcNotFound || GetRCTarget ( rc ) != rcDylib )
        {
            if ( rc != 0 )
                LOGERR( klogInt, rc, "SRAPathMake() failed" );
        }
        else
            rc = 0;
    }
    else
    {
        if ( !SRAPathTest( my_sra_path, *path ) )
        {
            char *buf = translate_accession( my_sra_path, *path, 64 );
            if ( buf != NULL )
            {
                free( (char*)(*path) );
                *path = buf;
            }
        }
        SRAPathRelease( my_sra_path );
    }
    return rc;
}


static rc_t pileup_main( Args * args )
{
    const AlignMgr *almgr = NULL;
    ReferenceIterator *ref_iter = NULL;
    KDirectory *dir = NULL;
    const VDBManager *vdb_mgr = NULL;
    uint32_t count;

    /* (1) make the align-manager ( necessary to make a ReferenceIterator... ) */
    rc_t rc = AlignMgrMakeRead ( &almgr );
    if ( rc != 0 )
        LOGERR( klogInt, rc, "AlignMgrMake() failed" );

    /* (2) make the reference-iterator */
    if ( rc == 0 )
    {
        PlacementRecordExtendFuncs cb_block;

        cb_block.data = ( void * )almgr;
        cb_block.destroy = NULL;
        cb_block.populate = populate_tooldata;
        cb_block.alloc_size = alloc_size;
        cb_block.fixed_size = 0;

        rc = AlignMgrMakeReferenceIterator ( almgr, &ref_iter, &cb_block, g_options.minmapq );
        if ( rc != 0 )
            LOGERR( klogInt, rc, "AlignMgrMakeReferenceIterator() failed" );
    }

    /* (3) make a k-directory ( necessary to make a vdb-manager ) */
    if ( rc == 0 )
    {
        rc = KDirectoryNativeDir( &dir );
        if ( rc != 0 )
            LOGERR( klogInt, rc, "KDirectoryNativeDir() failed" );
    }

    /* (4) make a vdb-manager */
    if ( rc == 0 )
    {
        rc = VDBManagerMakeRead ( &vdb_mgr, dir );
        if ( rc != 0 )
            LOGERR( klogInt, rc, "VDBManagerMakeRead() failed" );
    }

    /* (5) get the number of cmd-line-parameters loop through the given input-filenames */
    if ( rc == 0 )
    {
        rc = ArgsParamCount( args, &count );
        if ( rc != 0 )
            LOGERR( klogInt, rc, "ArgsParamCount() failed" );
    }

    /* (6) loop through the given input-filenames and load the ref-iter with it's input */
    if ( rc == 0 )
    {
        if ( count > 0 )
        {
            BSTree tree;
            rc = init_refranges( &tree, args );
            if ( rc == 0 )
            {
                uint32_t idx;
                check_refranges( &tree ); /* sanitize input... */
                for ( idx = 0; idx < count && rc == 0; ++idx )
                {
                    const char *filename = NULL;
                    rc = ArgsParamValue( args, idx, &filename );
                    if ( rc != 0 )
                        LOGERR( klogInt, rc, "ArgsParamvalue() failed" );
                    else
                    {
                        char * resolved = string_dup_measure ( filename, NULL );
                        if ( resolved != NULL )
                        {
                            if ( !is_this_a_filesystem_path( resolved ) )
                            {
                                rc = resolve_accession( dir, &resolved );
                            }
                            if ( rc == 0 )
                            {
                                rc = prepare_pileup_file( args, ref_iter, vdb_mgr, resolved, &tree );
                            }
                            free( resolved );
                        }
                        else
                            LOGERR( klogInt, rc, "string_dup_measure() failed" );
                    }
                }
                free_refranges( &tree );
            }
        }
        else
            rc = UsageSummary ( UsageDefaultName );
    }

    /* (7) walk the "loaded" ref-iterator ===> perform the pileup */
    if ( rc == 0 )
    {
        /* ============================================== */
        rc = walk_ref_iter( ref_iter, true );
        /* ============================================== */
    }

    if ( vdb_mgr != NULL ) VDBManagerRelease( vdb_mgr );
    if ( dir != NULL ) KDirectoryRelease( dir );
    if ( ref_iter != NULL ) ReferenceIteratorRelease( ref_iter );
    if ( almgr != NULL ) AlignMgrRelease ( almgr );
    return rc;
}

/* =========================================================================================== */

rc_t CC KMain( int argc, char *argv [] )
{
    rc_t rc = KOutHandlerSet( write_to_FILE, stdout );
    if ( rc != 0 )
        LOGERR( klogInt, rc, "KOutHandlerSet() failed" );
    else
    {
        Args * args;

        KLogHandlerSetStdErr();
        rc = ArgsMakeAndHandle( &args, argc, argv, 1,
            MyOptions, sizeof MyOptions / sizeof MyOptions [ 0 ] );
        if ( rc == 0 )
        {
            rc = get_pileup_options( args, &g_options );
            if ( rc == 0 )
            {
                if ( g_options.output_file != NULL )
                {
                    rc = set_stdout_to( g_options.gzip_output,
                                        g_options.bzip_output,
                                        g_options.output_file,
                                        32 * 1024 );
                }

                if ( rc == 0 )
                {
                    /* ======================= */
                    rc = pileup_main( args );
                    /* ======================= */
                }

                if ( g_options.output_file != NULL )
                    release_stdout_redirection();
            }
            ArgsWhack( args );
        }
    }
    return rc;
}
