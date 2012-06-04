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

#include <kapp/extern.h>
#include <sysalloc.h>

#include <klib/report.h>
#include <klib/rc.h>
#include <klib/container.h>
#include <klib/text.h>
#include <klib/printf.h>
#include <klib/vector.h>
#include <klib/log.h>
#include <klib/out.h>
#include <klib/status.h>
#include <klib/debug.h>
#include <kapp/main.h>
#include <kapp/args.h>

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#define USE_INLINING 0

KAPP_EXTERN bool CC is_valid_name (const char * string)
{
    /* we do not allow leading - or empty names */
    if ((*string == '\0') || (*string == '-'))
        return false;
    /* we do not allow ASCII control or ASCII space
     * but we do not disallow punctuation.  Use at your own risk */
    for ( ; *string != '\0'; ++string)
        if (isspace (*string) || iscntrl (*string))
            return false;
    return true;
}


/* ==========
 * ParamValue
 *   the value for an option is a NUL terminated ASCII / UTF-8 string
 */
typedef char ParamValue [1];

/*
 * Make
 *   allocated memory for a ASCII/UTF-8 string including one extra for NUL
 *   the value passed in must be a NUL termintaed string as well
 */
static
rc_t CC ParamValueMake (ParamValue ** pself, const char * value, size_t value_size)
{
    size_t alloc_size;

    assert (pself);
    assert (value);
    assert (value_size);

    alloc_size = sizeof (ParamValue) + value_size;
    *pself = malloc (alloc_size);

    if (*pself == NULL)
    {
        fprintf (stderr, "Error allocating memory for option parameter %s\n",
                 value);
        return RC (rcRuntime, rcArgv, rcConstructing, rcMemory, rcExhausted);
    }

    string_copy (**pself, alloc_size, value, value_size);
    return 0;
}

/*
 * Whack
 *   undo the Make.  That is free the memory
 */
static
void CC ParamValueWhack (void * self)
{
    assert (self); /* not an absolute requirement but if NULL we got programming error */

    free (self);
}

static
void CC ParamValueVectorWhack (void * self, void * ignored)
{
    assert (self);

    ParamValueWhack(self);
}


/* ==========
 * Option
 *  this is the primary node ofr an option
 *  It contains the name of the long option that is also the key to pulling
 *  out values.  The storage of these nodes is in a BSTree.
 */
typedef struct Option
{
    BSTNode     n;              /* BSTree node storage */

    Vector      values;         /* Vector of set values */
    uint32_t    count;          /* count of times option seen on the command line */
    uint32_t    max_count;      /* if non-zero, how many times is it legal to be seen */
    bool        needs_value;    /* does this option take/require a value? */
    bool        required;       
    size_t      size;           /* name length (size if UTF-8) */
    char        name [1];       /* key value The 1 will be the NUL */
} Option;

static
rc_t CC OptionMake (Option ** pself, const char * name, size_t size, uint32_t max_count, bool needs_value, bool required)
{
    Option *   self;

    assert (pself);
    assert (name);
    assert ((needs_value == true)||(needs_value == false)); /* not really but lets be rigorous */
    assert ((required == true)||(required == false)); /* not really but lets be rigorous */

    self = malloc (sizeof (*self) + size);
    if (self == NULL)
    {
        *pself = NULL;
        return RC (rcRuntime, rcArgv, rcConstructing, rcMemory, rcExhausted);
    }

    if ((self->needs_value = needs_value) != false)
        VectorInit (&self->values,0,4);
    else
        memset (&self->values, sizeof(self->values), 0);
    self->required = required;

    self->count = 0;
    self->max_count = max_count;
    self->size = size;
    string_copy (self->name, size+1, name, size);
    *pself = self;
    return 0;
}

static
void CC OptionWhack (Option * self)
{
    assert (self);

    if (self->needs_value)
        VectorWhack (&self->values, ParamValueVectorWhack, NULL);

    free (self);
}

static
void CC OptionTreeWhack (BSTNode * node, void * ignored)
{
    OptionWhack ((Option*)node);
}
#ifdef OptionGetName_needed
static
const char * CC OptionGetName (const Option * self, size_t * size)
{
    assert (self);
    if (size)
        *size = self->size;
    return self->name;
}
#endif
/*
 * NeedsValue
 *  return bool, does this option require a value
 */
static
bool CC OptionNeedsValue (const Option * self)
{
    assert (self);

    return self->needs_value;
}

/*
 * GetCount
 *  return the count of values seen so far.
 */
static
uint32_t CC OptionGetCount (const Option *self)
{
    assert (self);

    return self->count;
}

/* 
 * GetValue
 *    returns the address of the Nth value (not a copy of the value)
 */
static
rc_t CC OptionGetValue (const Option * self, uint32_t number, const char ** value)
{
    /* SKH -- not sure why this was here. 
       const char * pc; */
    uint32_t count;

    assert (self);
    assert (value);

    count = OptionGetCount(self);

    *value = VectorGet (&self->values, number);
    /* SKH -- this was checking pc, which is uninitialized */
    if (*value == NULL)
        return RC (rcRuntime, rcArgv, rcAccessing, rcIndex, rcExcessive);

    return 0;
}

/*
 * add a value to this node.  If a value isn't needed this is just incrementing the count
 */
static
rc_t CC OptionAddValue (Option * self, const char * value, size_t size)
{
    rc_t rc = 0;
    ParamValue * pvalue;

    assert (self);

    ++self->count;
    if (self->max_count && (self->count > self->max_count))
    {
        --self->count;
        fprintf (stderr, "Too many occurances of %s option\n", self->name);
        return RC (rcRuntime, rcArgv, rcInserting, rcData, rcExcessive);
    }
    if (self->needs_value)
    {
        assert (value);     /* gotta have a value */
        assert (size);      /* value can't be a NUL string */

        rc = ParamValueMake (&pvalue, value, size);
        if (rc == 0)
        {
            /* NOTE: effectively going in as a char ** though we will 
             * pull it out as a char* with the same value */
            rc = VectorAppend (&self->values, NULL, pvalue);
            if (rc)
            {
                fprintf (stderr, "error capturing parameter %s\n", self->name);
                ParamValueWhack (pvalue);
            }
        }
    }
    return rc;
}

static
int CC OptionCmp (const void * item, const BSTNode * n)
{
    const char * name;
    const Option * option;
    size_t l;

    assert (item);
    assert (n);

    name = item;
    option = (Option*)n;
    l = string_size (name);
    return string_cmp (name, l, option->name, option->size, (uint32_t)(l + option->size) );
}

static
int CC OptionSort (const BSTNode * item, const BSTNode * n)
{
    const Option * l = (Option*)item;
    const Option * r = (Option*)n;
    return string_cmp (l->name, l->size, r->name, r->size, (uint32_t)(l->size + r->size) );
}

/* ==========
 */
typedef struct OptAlias
{
    BSTNode     n;              /* BSTree Node storage */
    Option *   option;         /* full name node for which this is an alias */
    size_t      size;
    char        name[1];        /* utf8 name */
} OptAlias;

static
rc_t CC OptAliasMake (OptAlias ** pself, const char * name, size_t size,
                      Option * option)
{
    OptAlias * self;

    assert (pself);
    assert (name);
    assert (size);

    self = malloc (sizeof (*self) + size);
    if (self == NULL)
    {
        *pself = NULL;
        return RC (rcRuntime, rcArgv, rcConstructing, rcMemory, rcExhausted);
    }
    self->option = option;
    self->size = size;
    string_copy (self->name, size + 1, name, size);
    *pself = self;
    return 0;
}

static
void CC OptAliasWhack (const OptAlias * self)
{
    assert (self);
    free ((void*)self);
}

static
void CC OptAliasTreeWhack (BSTNode * node, void * ignored)
{
    assert (node); /* catch programming errors; freeing NULL wouldn't hurt */

    OptAliasWhack ((OptAlias*)node);
}

#if 0
static
const char * CC OptAliasName (const OptAlias * self, size_t * size)
{
    assert (self);

    if (size)
        *size = self->size;
    return self->name;
}
#endif
static
Option * CC OptAliasOption (const OptAlias *self)
{
    assert (self != NULL);

    return self->option;
}

static
int CC OptAliasCmp (const void * item, const BSTNode * n)
{
    const char * name;
    const OptAlias * option;
    size_t l;

    assert (item);
    assert (n);

    name = item;
    option = (OptAlias*)n;
    l = string_size (name);
    return string_cmp (name, l, option->name, option->size, (uint32_t)( l + option->size ) );
}

static
int CC OptAliasSort (const BSTNode * item, const BSTNode * n)
{
    const OptAlias * l;
    const OptAlias * r;

    assert (item);
    assert (n);

    l = (OptAlias*)item;
    r = (OptAlias*)n;
    return string_cmp (l->name, l->size, r->name, r->size, (uint32_t)( l->size + r->size ) );
}

#if NOT_USED_YET
static
rc_t CC OptDefMakeCopy (OptDef ** pself, OptDef * original)
{
    OptDef * self;
    size_t nsize;
    size_t asize;
    size_t hsize;

    nsize = string_size (original->name);
    asize = original->aliases ? string_size (original->aliases) : 0;
    hsize = original->help ? string_size (original->help) : 0;

    self = malloc (sizeof (*self) + nsize + 1 + asize + 1 + hsize + 1);
    if (self == NULL)
    {
        rc_t rc;
        /* assuming DebugMsg is stderr equivalent */
        rc = RC (rcRuntime, rcArgv, rcConstructing, rcMemory, rcExhausted);
        LOGERR (klogFatal, rc, "error creating help for option");
        return rc;
    }
    self->name = (char*)(self+1);
    string_copy (self->name, nsize + 1, original->name, nsize);
    self->aliases = self->name + nsize + 1;
    if (original->aliases)
    {
        string_copy (self->aliases, asize + 1, original->aliases, asize);
    }
    else
        self->aliases[0] = '\0';
    self->help = self->aliases + asize + 1;
    if (original->help)
    {
        string_copy (self->help, hsize + 1, original->help, asize);
    }
    else
        self->help[0] = '\0';
    self->max_count = original->max_count;
    self->needs_value = original->needs_value;
    *pself = self;
    return 0;
}
static
void CC OptDefCopyVectorWhack (void * self, void * ignored)
{
    free (self);
}

#endif

#if NOT_USED_YET

typedef struct HelpGroup
{
    rc_t ( CC * header_fmt) (Args * args, const char * header);
    Vector options;
    const char header[1];
} HelpGroup;


static
rc_t CC HelpGroupMake (HelpGroup ** pself, const char * name)
{
    HelpGroup * self;
    size_t size;

    size = string_size (name);
    self = malloc (sizeof (*self) + size);
    if (self == NULL)
    {
        fprintf (stderr, "Error allocating help group structure %s\n", name);
        return RC (rcRuntime, rcArgv, rcConstructing, rcMemory, rcExhausted);
    }
    string_copy (self->name, size+1, name, size);
    VectorInit (&self->optdefs, 0, 16);

    *pself = self;
    return 0;
}


static
rc_t CC HelpGroupAddOptDef (HelpGroup * self, OptDef * option)
{
    OptDef * opt;
    rc_t rc;

    rc = OptDefCopy (&opt, option);
    if (rc)
        return rc;

    rc = VectorAppend (&self->optdefs, NULL, opt);
    if (rc)
    {
        fprintf (stderr, "Error appending option for help\n");
        OptDefCopyVectorWhack (opt, NULL);
        return rc;
    }
    return 0;
}

static
void CC HelpGroupVectorWhack (void * item, void * ignored)
{
    HelpGroup * self = item;

    assert (self);
    VectorWhack (&self->optdefs, OptDefCopyVectorWhack, NULL);
    free (self);
}
#endif

/* ==========
 */
struct Args
{
    BSTree names;
    BSTree aliases;
    Vector argv;
    Vector params;
    Vector help;
#if NOT_USED_YET
    HelpGroup * def_help;
#endif
};

KAPP_EXTERN rc_t CC ArgsMake (Args ** pself)
{
    rc_t rc;
    Args * self;

    assert (pself);

    self = malloc (sizeof (*self));
    if (self == NULL)
    {
        rc = RC (rcRuntime, rcArgv, rcConstructing, rcMemory, rcExhausted);
    }
    else
    {
#if NOT_USED_YET
        HelpGroup * help;
#endif
        BSTreeInit (&self->names);
        BSTreeInit (&self->aliases);
        VectorInit (&self->argv,0,8);
        VectorInit (&self->params,0,8);
#if NOT_USED_YET
        VectorInit (&self->help,0,16);
        rc = HelpGroupMake (&help, "NCBI Options");
        if (rc)
        {
            ArgsWhack (self);
        }
        else
        {
            self->def_help = help;
            rc = VectorAppend (&self->help, NULL, help);
        }
#else
        rc = 0;
#endif
    }
    *pself = (rc == 0) ? self : NULL;
    return rc;
}

KAPP_EXTERN rc_t CC ArgsWhack (Args * self)
{
    if (self)
    {
        BSTreeWhack (&self->names, OptionTreeWhack, NULL);
        BSTreeWhack (&self->aliases, OptAliasTreeWhack, NULL);
        VectorWhack (&self->argv, ParamValueVectorWhack, NULL);
        VectorWhack (&self->params, NULL, NULL);
#if NOT_USED_YET
        VectorWhack (&self->help, HelpGroupVectorWhack, NULL);
#endif
        free (self);
    }
    return 0;
}

static
rc_t CC ArgsAddOption (Args * self, const OptDef * option)
{
    rc_t rc = 0;
    Option * node;
    Option * prev;
    const char * name;
    size_t size; /* will be used for the help tree side */

    if (self == NULL)
    {
        fprintf (stderr, "Error calling %s with NULL first parameter\n", __func__);
        return RC (rcRuntime, rcArgv, rcConstructing, rcSelf, rcNull);
    }
    if (option == NULL)
    {
        fprintf (stderr, "Error calling %s with NULL second parameter\n", __func__);
        return RC (rcRuntime, rcArgv, rcConstructing, rcParam, rcNull);
    }
    name = option->name;
    if (! is_valid_name (name))
    {
        fprintf (stderr, "Error using illegal option name %s\n", name);
        return RC (rcRuntime, rcArgv, rcConstructing, rcName, rcInvalid);
    }
    size = string_size (name);
    rc = OptionMake (&node, name, size, option->max_count, option->needs_value, option->required);
    if (rc)
    {
        fprintf (stderr, "Error initializing structures to parase parameters\n");
        return rc;
    }
    size ++;

    prev = NULL;
    rc = BSTreeInsertUnique (&self->names, &node->n, (BSTNode**)&prev, OptionSort);
    if (rc)
    {
        if (GetRCState(rc) == rcBusy)
        {
            rc = RC (rcRuntime, rcArgv, rcConstructing, rcName, rcBusy);
            fprintf (stderr, "duplicate option name %s\n", name);
        }
        else
            fprintf (stderr, "unknown error inserting %s\n", name);

        OptionWhack (node);
        return rc;

    }
    if (option->aliases)     /* any single character aliases? */
    {
        const char * startc;
        const char * endc;
        size_t incr;

        for ((startc = option->aliases),(endc = startc + string_size (startc));
             startc < endc; startc += incr)
        {
            OptAlias * snode;
            OptAlias * sprev;
            uint32_t c;
            char alias_name [8]; /* 4 should be enough + 1 for NUL */

            incr = utf8_utf32 (&c, startc, endc);
            if (incr < 0)
            {
                rc = RC (rcRuntime, rcArgv, rcConstructing, rcName, rcCorrupt);
                fprintf (stderr, "Error parsing alias string %s from %s for %s\n",
                         startc, option->aliases, name);
                break;
            }
            if (incr > 4)
            {
                rc = RC (rcRuntime, rcArgv, rcConstructing, rcName, rcCorrupt);
                fprintf (stderr, "Error parsing UTF-8 string %s\n", startc);
                break;
            }
            string_copy (alias_name, sizeof (alias_name), startc, incr);
            if (! is_valid_name (alias_name))
            {
                rc = RC (rcRuntime, rcArgv, rcConstructing, rcName, rcInvalid);
                fprintf (stderr, "Error using invalid alias name %s", alias_name);
                break;
            }
            size += incr + 1;
            rc = OptAliasMake (&snode, alias_name, incr, node);
            if (rc)
            {
                fprintf (stderr, "Error creating structure for alias %s for parameter %s\n", 
                         alias_name, node->name);
                break;
            }
            rc = BSTreeInsertUnique (&self->aliases, &snode->n, (BSTNode**)&sprev, OptAliasSort);
            if (rc)
            {
                if (GetRCState(rc) == rcExists)
                {
                    rc = RC (rcRuntime, rcArgv, rcConstructing, rcName, rcExists);
                    fprintf (stderr, "duplicate alias name %s\n", alias_name);
                }
                else
                    fprintf (stderr, "%d unknown error inserting alias %s\n", rc, alias_name);

                OptAliasWhack (snode);
                break;
            }
        }
    }
    if (rc == 0)
    {
    }
    return rc;
}

KAPP_EXTERN rc_t CC ArgsAddOptionArray (Args * self, const OptDef * option, uint32_t count /*, 
                                                                                             rc_t (*header_fmt)(Args * args, const char * header),
                                                                                             const char * header */)
{
    rc_t rc;
#if NOT_USED_YET
    HelpGroup * hg;

    rc = HelpGroupMake (&hg, header, header_fmt, option, count);
    if (rc == 0)
    {

        rc = VectorAppend (&self->help, NULL, hg);
        if (rc == 0)
        {

            /* count might be zero for help only sections */
            for (rc = 0; count > 0; --count, ++option)
            {
                rc = ArgsAddOption (self, option);
                if (rc)
                    break;
            }
            if (rc == 0)
                return 0;
        }
        else
            HelpGroupVectorWhack (hg, NULL);
    }
#else
    for (rc = 0; (rc == 0) && (count > 0); --count, ++option)
    {
        rc = ArgsAddOption (self, option);
    }
#endif
    return rc;
}

/*
 */
KAPP_EXTERN rc_t CC next_arg (const Args * self, int * pix, int max, const char ** pvalue)
{
    int ix;

    if (*pix >= max)
        return RC (rcApp, rcArgv, rcAccessing, rcString, rcNotFound);

    ix = *pix;
    ix++;
    *pvalue = (const char *) VectorGet (&self->argv, ix);
    *pix = ix;
    return 0;
}

typedef struct ArgsReqCheckData
{
    Option * missing_option;
} ArgsReqCheckData;

static
bool CC ArgsCheckRequired (BSTNode * n, void * _data)
{
#if 0
/* Need to rethink this. If help is asked for this shouldn't fail */
    ArgsReqCheckData * data;
    Option * opt;

    data = _data;
    opt = (Option*)n;

    if (opt->required && ! opt->count)
    {
        data->missing_option = opt;
        return true;
    }
#endif
    return false;
}

KAPP_EXTERN rc_t CC ArgsParse (Args * self, int argc, char *argv[])
{
    rc_t rc = 0;
    int ix;
    size_t jx;
    Option * node;
    const char * parg;
    char * equal_sign;
    ParamValue * arg;
    char name [32];
    const char * value = NULL;
    size_t value_len;
    bool needs_value;

    /* first capture original argument list and store in our format */
    for (ix = 0; ix < argc; ++ix)
    {
        size_t len;

        parg = argv[ix];
        len = strlen (parg);

        rc = ParamValueMake (&arg, parg, len);
        if (rc == 0)
            rc = VectorAppend (&self->argv, NULL, arg);
        if (rc)
            break;
    }
    
    if (rc)
        return rc;

    for (ix = 1; ix < argc; ++ix)
    {
        parg = (const char *)VectorGet (&self->argv, ix);
        if (parg[0] != '-')
        {
            /* we can do this because it is already a (const char *)
             * and (ParamValue *) after the first loop */
            rc = VectorAppend (&self->params, NULL, parg);
            if (rc)
                break;
        }
        else
        {
            node = NULL;
            if (parg[1] == '-')
            {
                size_t nlen=string_copy (name, sizeof (name), parg+2, string_size (parg+2));
                equal_sign = string_chr (name, nlen, '=');
                if (equal_sign)
                    *equal_sign = '\0';

                node = (Option*)BSTreeFind (&self->names, name, OptionCmp);
                if (node == NULL)
                {
                    rc = RC (rcApp, rcArgv, rcParsing, rcParam, rcUnknown);
                    fprintf (stderr, "Unknown argument %s\n", name);
                    break;
                }
                else
                {
                    needs_value = OptionNeedsValue (node);
                    if (needs_value)
                    {
                        if (equal_sign != NULL)
                            value = parg + 2 + ((equal_sign+1) - name);
                        else if( ix + 1 >= argc ) {
                            rc = RC (rcRuntime, rcArgv, rcParsing, rcParam, rcExhausted);
                            fprintf (stderr, "Option '%s' is missing a value\n", node->name);
                        }
                        else
                            rc = next_arg (self, &ix, argc, &value);

                        if (rc)
                            break;

                        assert (value != NULL);

                        value_len = string_size (value);

                        rc = OptionAddValue (node, value, value_len);
                        if (rc)
                        {
                            break;
                        }
                    }
                    else
                    {
                        rc = OptionAddValue (node, NULL, 0);
                        if (rc)
                        {
                            break;
                        }
                    }
                }
            }
            else
            {
                const char * end;
                uint32_t name_len;

                end = parg + string_size (parg);
                for (jx = 1; parg[jx]; jx += name_len)
                {
                    OptAlias *alias;
                    uint32_t c;

                    name_len = utf8_utf32 ( &c, parg+jx, end); 
                    string_copy (name, sizeof (name), parg+jx, name_len);

                    alias = (OptAlias*)BSTreeFind (&self->aliases, name, OptAliasCmp);
                    if (alias == NULL)
                    {
                        rc = RC (rcApp, rcArgv, rcParsing, rcParam, rcUnknown);
                        fprintf (stderr, "Unknown argument %s\n", name);
                    }
                    else
                    {
                        node = OptAliasOption (alias);
                        if (OptionNeedsValue(node))
                        {
                            if (parg[jx+name_len] == '=')
                            {
                                ++jx;
                                if (parg[jx+name_len])
                                    value = parg + jx + name_len;
                                else
                                {
                                    rc = RC (rcRuntime, rcArgv, rcParsing, rcParam, rcExhausted);
                                    fprintf (stderr, "Value missing with alias followed by =\n");
                                }
                            }
                            else if(parg[jx+name_len])
                            {
                                value = parg + jx + name_len;
                            }
                            else if( ix + 1 >= argc ) {
                                rc = RC (rcRuntime, rcArgv, rcParsing, rcParam, rcExhausted);
                                fprintf (stderr, "Option '%s' is missing a value\n", node->name);
                            }
                            else
                                value = argv [++ix];
                            if (rc == 0)
                                rc = OptionAddValue (node, value, string_size (value));
                            break;
                        }
                        else
                        {
                            rc = OptionAddValue (node, NULL, 0);
                            if (rc)
                                break;
                        }
                    }
                }
            }
        }
        if (rc)
            break;
    }

    /* check for required parameters */
    if (rc == 0)
    {
        bool b;
        ArgsReqCheckData r;


        b = BSTreeDoUntil (&self->names, false, ArgsCheckRequired, &r);
        if (b)
        {
            fprintf (stderr, "Error missing required parameter %.*s\n",
                     (int)r.missing_option->size, r.missing_option->name);
            rc = RC (rcRuntime, rcArgv, rcParsing, rcParam, rcNotFound);
        }
    }

    if (rc) {
        ReportSilence();
    }

    return rc;
}

KAPP_EXTERN rc_t CC ArgsOptionCount (const Args * self, const char * option_name, uint32_t * count)
{
    rc_t rc;

    if (self == NULL)
        return RC (rcRuntime, rcArgv, rcAccessing, rcSelf, rcNull);
    else if (count == NULL)
        return RC (rcRuntime, rcArgv, rcAccessing, rcParam, rcNull);
    else
    {
        const Option * node;

        node = (const Option*)BSTreeFind (&self->names, option_name, OptionCmp);
        if (node == NULL)
        {
            rc = RC (rcRuntime, rcArgv, rcAccessing, rcName, rcNotFound);
            PLOGERR (klogWarn, (klogWarn, rc, "Option name not found '%s'", option_name));
            return rc;
        }

        *count = OptionGetCount (node);
        return 0;
    }
}

KAPP_EXTERN rc_t CC ArgsOptionValue (const Args * self, const char * option_name, uint32_t iteration,
                                     const char ** value_string)
{
    const Option * node;

    if (self == NULL)
        return RC (rcRuntime, rcArgv, rcAccessing, rcSelf, rcNull);

    if ((option_name == NULL) || (value_string == NULL))
        return RC (rcRuntime, rcArgv, rcAccessing, rcParam, rcNull);

    *value_string = NULL;

    node = (const Option*)BSTreeFind (&self->names, option_name, OptionCmp);
    if (node == NULL)
        return RC (rcRuntime, rcArgv, rcAccessing, rcName, rcNotFound);
    else
        return  OptionGetValue (node, iteration, value_string);
}

KAPP_EXTERN rc_t CC ArgsParamCount (const Args * self, uint32_t * count)
{
    if (self == NULL)
        return RC (rcRuntime, rcArgv, rcAccessing, rcSelf, rcNull);
    else if (count == NULL)
        return RC (rcRuntime, rcArgv, rcAccessing, rcParam, rcNull);

    *count = VectorLength (&self->params);
    return 0;
}

KAPP_EXTERN rc_t CC ArgsParamValue (const Args * self, uint32_t iteration, const char ** value_string)
{
    if (self == NULL)
        return RC (rcRuntime, rcArgv, rcAccessing, rcSelf, rcNull);

    if (value_string == NULL)
        return RC (rcRuntime, rcArgv, rcAccessing, rcParam, rcNull);

    if (iteration >= VectorLength (&self->params))
    {
        *value_string = NULL;
        return RC (rcRuntime, rcArgv, rcAccessing, rcParam, rcExcessive);
    }

    *value_string = (const char*) VectorGet (&self->params, iteration);
    return 0;
}

KAPP_EXTERN rc_t CC ArgsArgvCount (const Args * self, uint32_t * count)
{
    if (self == NULL)
        return RC (rcRuntime, rcArgv, rcAccessing, rcSelf, rcNull);
    else if (count == NULL)
        return RC (rcRuntime, rcArgv, rcAccessing, rcArgv, rcNull);

    *count = VectorLength (&self->argv);
    return 0;
}
#ifdef ArgsArgc
#undef ArgsArgc
#endif

KAPP_EXTERN rc_t CC ArgsArgc (const Args * self, uint32_t * count)
{
    return ArgsArgvCount (self, count);
}

KAPP_EXTERN rc_t CC ArgsArgvValue (const Args * self, uint32_t iteration, const char ** value_string)
{
    if (self == NULL)
        return RC (rcRuntime, rcArgv, rcAccessing, rcSelf, rcNull);

    if (value_string == NULL)
        return RC (rcRuntime, rcArgv, rcAccessing, rcArgv, rcNull);

    if (iteration >= VectorLength (&self->argv))
    {
        *value_string = NULL;
        return RC (rcRuntime, rcArgv, rcAccessing, rcArgv, rcExcessive);
    }

    *value_string = (const char*) VectorGet (&self->argv, iteration);
    return 0;
}

#define USAGE_MAX_SIZE 81
static
const char * verbose_usage[] = 
{ "Increase the verbosity level of the program.",
  "Use multiple times for more verbosity.", NULL };
static
const char * debug_usage[] = 
{ "Turn on debug output for module.",
  "All flags if not specified.", NULL };
static
const char * help_usage[] = 
{ "Output a brief explantion for the program.", NULL };
static
const char * report_usage[] = 
{ "Control program execution environment report generation (if implemented).",
    "One of (never|error|always). Default is error", NULL };
static
const char * version_usage[] = 
{ "Display the version of the program then quit.", NULL };
static 
char log0 [USAGE_MAX_SIZE];
static 
char log1 [USAGE_MAX_SIZE];
static
const char * log_usage[] = 
{ "Logging level as number or enum string.", log0, log1, NULL };

static
void CC gen_log_usage (const char ** _buffers)
{
    static const char div[] = "|";

    /* these have to be mutable for this to work */
    char ** buffers = (char **) _buffers;

    char * pc;
    char * p0;
    char * p1;
    size_t rem;
    size_t used;
    KLogLevel level;

    rc_t rc;
    char buffv[USAGE_MAX_SIZE];


    p0 = pc = buffers[1];
    p1 = pc = buffers[2];

    *p0 = *p1 = '\0';

    rem = USAGE_MAX_SIZE; /* makes an assumption */
    pc = buffv;
    for (level = klogLevelMin; level <= klogLevelMax; ++level)
    {
        rc = KLogLevelExplain (level, pc, rem, &used);
        if (rc || used == 0)
        {
            p0 = NULL;
            pc = NULL;
            break;
        }
        pc += used;
        rem -= used;
        strcpy (pc, div);
        pc += sizeof div - 1;
        rem -= sizeof div - 1;
    }
    if (p0)
    {
        pc -= sizeof div - 1;
        rem += sizeof div - 1;
        *pc = '\0';
        
        rc = string_printf (p0, USAGE_MAX_SIZE, & used,
                            "One of (%s) or (%u-%u)",
                            buffv, klogLevelMin, klogLevelMax);
        if (used == 0)
            p0 = NULL;
    }
    rc = KLogLevelExplain (KLogLevelGet(), buffv, sizeof buffv, &used);
    if (rc || used == 0)
        p1 = NULL;
    else
    {
        buffv[used] = '\0';
        rc = string_printf (p1, USAGE_MAX_SIZE, & used,
                         "Current/default is %s",
                         buffv);
        if (used == 0)
            p1 = NULL;
    }
    if (p0 == NULL)
    {
        p0 = p1;
        p1 = NULL;
    }
}



OptDef StandardOptions[]  =
{
    {
        OPTION_HELP,    ALIAS_HELP,    NULL, 
        help_usage,
        OPT_UNLIM, false, false
    },
    {
        OPTION_VERSION, ALIAS_VERSION, NULL,
        version_usage,
        OPT_UNLIM, false, false
    },
    {
        OPTION_LOG_LEVEL, ALIAS_LOG_LEVEL, gen_log_usage,
        log_usage,
        OPT_UNLIM, true, false
    },
    {
        OPTION_VERBOSE, ALIAS_VERBOSE, NULL,
        verbose_usage,
        OPT_UNLIM, false, false
    },
    {
        OPTION_DEBUG, ALIAS_DEBUG, NULL,
        debug_usage, 
        OPT_UNLIM, true, false
    },
    {   /* OPTION_REPORT is used in klib/report.c */
        OPTION_REPORT, NULL, NULL,
        report_usage, 
        OPT_UNLIM, true, false
    }
};

KAPP_EXTERN rc_t CC ArgsAddStandardOptions(Args * self)
{
    return ArgsAddOptionArray (self, StandardOptions,
                               sizeof (StandardOptions) / sizeof (OptDef)
                               /*, NULL, NULL */ );
}

KAPP_EXTERN rc_t CC ArgsMakeStandardOptions (Args** pself)
{
    Args * self;
    rc_t rc;

    rc = ArgsMake (&self);
    if (rc == 0)
    {
        rc = ArgsAddStandardOptions (self);
        if (rc)
            ArgsWhack (self);
    }
    *pself = (rc == 0) ? self : NULL;
    return rc;
}

KAPP_EXTERN rc_t CC ArgsHandleHelp (Args * self)
{
    uint32_t count;
    rc_t rc;

    rc = ArgsOptionCount (self, OPTION_HELP, &count);
    if (rc == 0)
    {
        if (count)
        {
            /* this is a call into the main program code */
            Usage(self);
            ArgsWhack (self);
            exit (0);
        }
    }
    return rc;
}


KAPP_EXTERN rc_t CC ArgsHandleVersion (Args * self)
{
    uint32_t count;
    rc_t rc;

    rc = ArgsOptionCount (self, OPTION_VERSION, &count);
    if (rc == 0)
    {
        if (count)
        {
            const char * progname = UsageDefaultName;
            const char * fullpath = UsageDefaultName;

            if (self)
                rc = ArgsProgram (self, &fullpath, &progname);

            HelpVersion (fullpath, KAppVersion());

            ArgsWhack (self);
            exit (0);
        }
    }
    return rc;
}


KAPP_EXTERN rc_t CC ArgsHandleLogLevel (const Args * self)
{
    uint32_t count;
    rc_t rc;

    rc = ArgsOptionCount (self, OPTION_LOG_LEVEL, &count);
    if (rc == 0)
    {
        if (count == 0)
        {
        }
        else
        {
            const char * value;
            uint32_t ix;

            for (ix = 0; ix < count; ++ix)
            {
                rc = ArgsOptionValue (self, OPTION_LOG_LEVEL,
                                      ix, &value);
                if (rc == 0)
                    rc = LogLevelSet (value);
                if (rc)
                    break;
            }
        }
    }
    return rc;
}

KAPP_EXTERN rc_t CC ArgsHandleStatusLevel (const Args * self)
{
    uint32_t count;
    rc_t rc;

    rc = ArgsOptionCount (self, OPTION_VERBOSE, &count);
    if (rc == 0) {
        KStsLevelSet (count);
    }
    return rc;
}
#if _DEBUGGING
rc_t CC ArgsHandleDebug (const Args * self)
{
    uint32_t count;
    rc_t rc;

    rc = ArgsOptionCount (self, OPTION_DEBUG, &count);
    if (rc == 0)
    {
        if (count == 0)
        {
        }
        else
        {
            const char * value;
            uint32_t ix;

            for (ix = 0; ix < count; ++ix)
            {
                rc = ArgsOptionValue (self, OPTION_DEBUG,
                                      ix, &value);
                if (rc == 0)
                    rc = KDbgSetString (value);
                if (rc)
                    break;
            }
        }
    }
    return rc;
}
#endif

KAPP_EXTERN rc_t CC ArgsHandleStandardOptions (Args * self)
{
    rc_t rc;

    do
    {
        rc = ArgsHandleHelp (self);
        if (rc)
            break;

        rc = ArgsHandleVersion (self);
        if (rc)
            break;

        rc = ArgsHandleLogLevel (self);
        if (rc)
            break;

        rc = ArgsHandleStatusLevel (self);
        if (rc)
            break;

#if _DEBUGGING
        rc = ArgsHandleDebug (self);
#endif
    } while (0); /* not a real loop */
    return rc;
}

KAPP_EXTERN rc_t CC ArgsMakeAndHandle (Args ** pself, int argc, char ** argv, uint32_t table_count, ...)
{
    rc_t rc;
    Args * self;

    *pself = NULL;
    rc = ArgsMakeStandardOptions (&self);
    if (rc == 0)
    {
        do
        {
            if (table_count)
            {
                va_list ap;

                va_start (ap, table_count);

                while (table_count-- && (rc == 0))
                {
                    OptDef * options;
                    uint32_t opt_count;

                    options = va_arg (ap, OptDef *);
                    opt_count = va_arg (ap, uint32_t);

                    rc = ArgsAddOptionArray (self, options, opt_count /* , NULL, NULL */);
                }
                va_end (ap);
                if (rc)
                    break;
            }

            rc = ArgsParse (self, argc, argv);
            if (rc)
                break;

            rc = ArgsHandleStandardOptions (self);
            if (rc)
                break;

            *pself = self;

            return 0;

        } while (0);
    
        ArgsWhack (self);
    }
    return rc;
}




/* NOTE:
 * This needs to move into a unix/win32 seperated file
 * and quite probably outside of args.
 */

KAPP_EXTERN const char * CC trim_path (const char * full_name)
{
    const char * name;

    name = strrchr (full_name, '/');
    if (name == NULL)
        name = full_name;
    else
        ++name; /* skip past '/' */
    return name;
}


KAPP_EXTERN rc_t CC ArgsProgram (const Args * args, const char ** fullpath, const char ** progname)
{
    const char * defaultname = UsageDefaultName;
    const char * f;
    rc_t rc;

    rc = ArgsArgvValue (args, 0, &f);
    if (rc == 0)
    {
        if (fullpath)
            *fullpath = f;
        if (progname)
            *progname = trim_path (f);
    }
    else
    {
        f = defaultname;
        
        if (fullpath != NULL)
        {
            if (*fullpath == NULL)
                *fullpath = f;
            else
                f = *fullpath;
        }
        if (progname)
        {
            if (*progname == NULL)
                *progname = trim_path (f);
        }
    }
    return rc;
}

void CC HelpVersion (const char * fullpath, ver_t version)
{
    OUTMSG (("\n%s : %.3V\n\n", fullpath, version));
}


static void print_indented( const size_t first_indent, const size_t indent,
                            const size_t max_line_len, const char ** msgs )
{
    const char * msg;
    size_t line_len;

    if ( *msgs == NULL )
    {
        OUTMSG(( "\n" ));
        return;
    }

    if ( first_indent < indent )
    {
        OUTMSG(( "%*s", indent - first_indent, " " ));
        line_len = indent;
    }
    else
    {
        OUTMSG(( "  " ));
        line_len = first_indent + 2;
    }
    while ( ( msg = *msgs++ ) != NULL )
    {
        while ( msg != NULL )
        {
            const char * space = strchr( msg, ' ' );
            if ( space != NULL )
            {
                /* space found, can we print the word on the current line? */
                int wordlen = ( space - msg );
                if ( ( line_len + wordlen + 1 ) < max_line_len )
                {
                    if ( wordlen > 1 )
                        OUTMSG(( "%.*s", wordlen + 1, msg )); /* yes */
                }
                else
                {
                    OUTMSG(( "\n%*s%.*s", indent, " ", wordlen + 1, msg )); /* no: new line */
                    line_len = indent;
                }
                line_len += ( wordlen + 1 );
                msg += ( wordlen + 1 );
            }
            else
            {
                /* no space found, can we print the string on the current line? */
                size_t remainder = string_size( msg );
                if ( line_len + remainder < max_line_len )
                {
                    OUTMSG(( "%s ", msg )); /* yes */
                }
                else
                {
                    OUTMSG(( "\n%*s%s ", indent, " ", msg )); /* no: new line */
                    line_len = indent;
                }
                line_len += remainder;
                msg = NULL; /* we are done with one source-line... */
            }
        }
    }
    OUTMSG(( "\n" ));
}

void CC HelpOptionLine(const char * alias, const char * option, const char * param, const char ** msgs)
{
    int n, msgc;
/*    const char * msg; */
#define INDENT 2
#define MSG_INDENT 35
#define MSG_MAXLEN 80

    bool has_alias = (alias != NULL && alias[0] != '\0');
    bool has_opt = (option != NULL && option[0] != '\0');

    if( has_alias || has_opt ) {
        OUTMSG(("%*s%n", INDENT, " ", &n));
        msgc = n;
        if( has_alias ) {
            OUTMSG(("-%s%n", alias, &n));
            msgc += n;
        }
        if( has_alias && has_opt ) {
            OUTMSG(("|"));
            msgc++;
        }
        if( has_opt ) {
            OUTMSG(("--%s%n", option, &n));
            msgc += n;
        }
        if( param != NULL) {
            OUTMSG((" <%s>%n", param, &n));
            msgc += n;
        }
        print_indented( msgc, MSG_INDENT, MSG_MAXLEN, msgs );
    }
}

void CC HelpParamLine (const char * param, const char * const * msgs)
{
    int msgc;
    const char * msg;

    msg = *msgs++;

    if (param)
    {
        OUTMSG (("%*s%s%n", INDENT, " ", param, &msgc));
	if (msg == NULL)
	    OUTMSG (("\n"));
	else
	{
	    OUTMSG (("%*s%s\n", MSG_INDENT-msgc, " ", msg));
	}
    }
    if (msg != NULL)
	while ((msg = *msgs++) != NULL)
	    OUTMSG (("%*s%s\n", MSG_INDENT, " ", msg));
}

void CC HelpOptionsStandard (void)
{
    HelpOptionLine (ALIAS_HELP1, OPTION_HELP, NULL, help_usage);

    HelpOptionLine (ALIAS_VERSION, OPTION_VERSION, NULL, version_usage);

    gen_log_usage (log_usage);
    HelpOptionLine (ALIAS_LOG_LEVEL, OPTION_LOG_LEVEL, "level", log_usage);

    HelpOptionLine (ALIAS_VERBOSE, OPTION_VERBOSE, NULL, verbose_usage);
    HelpOptionLine (NULL, OPTION_REPORT, "type", report_usage);
#if _DEBUGGING
    HelpOptionLine (ALIAS_DEBUG, OPTION_DEBUG, "Module[-Flag]", debug_usage); 
#endif
}


rc_t CC MiniUsage (const Args * args)
{
    KWrtWriter w;
    void * d;
    const char * progname;
    rc_t rc;

    w = KOutWriterGet();
    d = KOutDataGet();

    rc = ArgsProgram (args, NULL, &progname);
    if (rc)
        progname = UsageDefaultName;
    KOutHandlerSetStdErr();
    UsageSummary (progname);
    KOutMsg ("\nUse option --help for more information.\n\n");

    KOutHandlerSet (w,d);

    return rc;
}
