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
#include <align/extern.h>

#include <klib/rc.h>
#include <klib/vector.h>
#include <klib/refcount.h>
#include <klib/sort.h>
#include <klib/text.h>
#include <klib/out.h>
#include <insdc/insdc.h>
#include <align/manager.h>
#include <align/iterator.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>


typedef struct window
{
    INSDC_coord_zero first;
    INSDC_coord_zero last;
    INSDC_coord_len len;
} window;


typedef struct pi_entry
{
    DLNode n;                       /* to have it in a DLList */
    PlacementIterator *pi;          /* the placement-iterator we have added */
    window nxt_avail;               /* the next available position of the placement-iterator */
    window w;                       /* the window of the placement-iterator */
} pi_entry;


typedef struct pi_ref
{
    DLNode n;                       /* to have it in a DLList */
    const char * name;              /* the name of the reference it referes to */
    DLList pi_entries;              /* it has a DLList of pi_entry-struct's */
    pi_entry * current_entry;       /* the current entry to walk pi_entries, one !!record!! at a time */
    INSDC_coord_zero first_pos;     /* where do we have to start? */
    INSDC_coord_zero last_pos;      /* where do we have to stop? */
    bool first_pos_initialized;     /* is the start-position initialized */
    struct ReferenceObj const *refobj;
} pi_ref;


struct PlacementSetIterator
{
    KRefcount refcount;
    struct AlignMgr const *amgr;    /* the alignment-manager... ( right now: we store it, but that's it )*/
    DLList pi_refs;                 /* a list of references we have to iterate over... */
    pi_ref * current_pi_ref;        /* what is the current pi_ref, we are handling ? */
};


LIB_EXPORT rc_t CC AlignMgrMakePlacementSetIterator ( struct AlignMgr const *self,
    PlacementSetIterator **iter )
{
    rc_t rc = 0;
    if ( self == NULL )
        rc = RC( rcAlign, rcIterator, rcConstructing, rcSelf, rcNull );
    else
    {
        if ( iter == NULL  )
            rc = RC( rcAlign, rcIterator, rcConstructing, rcParam, rcNull );
        else
        {
            PlacementSetIterator * psi = calloc( sizeof * psi, 1 );
            if ( psi == NULL )
                rc = RC( rcAlign, rcIterator, rcConstructing, rcMemory, rcExhausted );
            else
            {
                rc = AlignMgrAddRef ( self );
                if ( rc == 0 )
                {
                    KRefcountInit( &psi->refcount, 1, "PlacementSetIterator", "Make", "align" );
                    psi->amgr = self;
                    psi->current_pi_ref = NULL;          /* we don't know that yet */
                    DLListInit( &psi->pi_refs );
                }
            }

            if ( rc == 0 )
                *iter = psi;
            else
                free( psi );
        }
    }

    return rc;
}


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


typedef struct pi_ref_cb_ctx
{
    const char * name;
    pi_ref *res;
} pi_ref_cb_ctx;

static bool CC find_pi_ref_callback( DLNode *n, void *data )
{
    pi_ref_cb_ctx *ctx = ( pi_ref_cb_ctx * )data;
    pi_ref * pr = ( pi_ref * ) n;
    if ( cmp_pchar( ctx->name, pr->name ) == 0 )
    {
        ctx->res = pr;
        return true;
    }
    else
    {
        return false;
    }
}


static pi_ref * find_pi_ref( const DLList * list, const char * name )
{
    pi_ref_cb_ctx ctx;
    ctx.res = NULL;
    ctx.name = name;
    DLListDoUntil ( list, false, find_pi_ref_callback, &ctx );
    return ctx.res;
}


static rc_t make_pi_ref( pi_ref ** pr, const char * name )
{
    rc_t rc = 0;
    *pr = calloc( 1, sizeof ** pr );
    if ( *pr == NULL )
        rc = RC( rcAlign, rcIterator, rcConstructing, rcMemory, rcExhausted );
    else
    {
        (*pr)->name = name;
        DLListInit( &( (*pr)->pi_entries ) );
    }
    return rc;
}


static rc_t add_to_pi_ref( pi_ref * pr, window * w, PlacementIterator *pi )
{
    rc_t rc = 0;
    bool added = false;
    pi_entry * pie = calloc( 1, sizeof *pie );
    if ( pie == NULL )
        rc = RC( rcAlign, rcIterator, rcConstructing, rcMemory, rcExhausted );
    else
    {
        pie->pi = pi;       /* store the placement-iterator in it's entry-struct */
        pie->w.first = w->first;
        pie->w.last = w->last;
        pie->w.len = w->len;

        rc = PlacementIteratorNextAvailPos ( pi, &(pie->nxt_avail.first), &(pie->nxt_avail.len) );
        if ( ( rc == 0 ) || ( GetRCState( rc ) == rcDone ) )
        {
            pie->nxt_avail.last = pie->nxt_avail.first + pie->nxt_avail.len - 1;
            if ( pie->nxt_avail.last >= w->first && pie->nxt_avail.first <= w->last )
            {
                /* finally add the iterator to our list */
                DLListPushTail ( &pr->pi_entries, ( DLNode * )pie );
                added = true;
            }
            rc = 0;
        }
        if ( !added )
        {
            PlacementIteratorRelease ( pi );
            free( pie );
        }
    }
    return rc;
}


LIB_EXPORT rc_t CC PlacementSetIteratorAddPlacementIterator ( PlacementSetIterator *self,
    PlacementIterator *pi )
{
    rc_t rc = 0;
    if ( self == NULL )
        rc = RC( rcAlign, rcIterator, rcConstructing, rcSelf, rcNull );
    else
    {
        if ( pi == NULL  )
            rc = RC( rcAlign, rcIterator, rcConstructing, rcParam, rcNull );
        else
        {
            const char * name;      /* what reference are we aligning against */
            window w;               /* where does the pi start/end, against said reference */

            /* to find the name of the reference used, important for adding the iterator */
            rc = PlacementIteratorRefWindow ( pi, &name, &(w.first), &(w.len) );
            if ( rc == 0 )
            {
                pi_ref * pr = find_pi_ref( &self->pi_refs, name );
                w.last = w.first + w.len - 1;
                if ( pr == NULL )
                {
                    /* we do not have a pi_ref yet with this name: make one! */
                    rc = make_pi_ref( &pr, name );
                    if ( rc == 0 )
                    {
                        DLListPushTail ( &self->pi_refs, ( DLNode * )pr );
                    }
                }
                if ( rc == 0 )
                {
                    /* add the placement-iterator to the newly-made or existing pi_ref! */
                    rc = add_to_pi_ref( pr, &w, pi );
                }
                if ( rc == 0 )
                {
                    if ( pr->first_pos_initialized )
                    {
                        if ( w.first < pr->first_pos )
                            pr->first_pos = w.first;
                        if ( w.last > pr->last_pos )
                            pr->last_pos = w.last;
                    }
                    else
                    {
                        pr->first_pos = w.first;
                        pr->last_pos = w.last;
                        pr->first_pos_initialized = true;
                    }
                }
            }
        }
    }
    return rc;
}


LIB_EXPORT rc_t CC PlacementSetIteratorAddRef ( const PlacementSetIterator *cself )
{
    rc_t rc = 0;
    if ( cself == NULL )
        rc = RC( rcAlign, rcIterator, rcAttaching, rcSelf, rcNull );
    else
    {
        if ( KRefcountAdd( &cself->refcount, "PlacementSetIterator" ) != krefOkay )
        {
            rc = RC( rcAlign, rcIterator, rcAttaching, rcError, rcUnexpected );
        }
    }
    return rc;
}


static void CC pi_entry_whacker( DLNode *n, void *data )
{
    pi_entry * pie = ( pi_entry * )n;
    PlacementIteratorRelease ( pie->pi );
    free( pie );
}

static void CC pi_ref_whacker( DLNode *n, void *data )
{
    pi_ref * pr = ( pi_ref * )n;
    DLListWhack ( &pr->pi_entries, pi_entry_whacker, NULL );
    free( pr );
}

LIB_EXPORT rc_t CC PlacementSetIteratorRelease ( const PlacementSetIterator *cself )
{
    rc_t rc = 0;
    if ( cself == NULL )
        rc = RC( rcAlign, rcIterator, rcReleasing, rcSelf, rcNull );
    else
    {
        if ( KRefcountDrop( &cself->refcount, "PlacementSetIterator" ) == krefWhack )
        {
            PlacementSetIterator * self = ( PlacementSetIterator * ) cself;
            /* release the DLList of pi-ref's and the pi's in it... */
            DLListWhack ( &self->pi_refs, pi_ref_whacker, NULL );
            AlignMgrRelease ( self->amgr );
            free( self );
        }
    }
    return rc;
}


LIB_EXPORT rc_t CC PlacementSetIteratorNextReference ( PlacementSetIterator *self,
    INSDC_coord_zero *first_pos, INSDC_coord_zero *last_pos,
    struct ReferenceObj const ** refobj )
{
    rc_t rc = 0;
    if ( refobj != NULL ) { *refobj = NULL; }
    if ( first_pos != NULL ) { *first_pos = 0; }
    if ( last_pos != NULL ) { *last_pos = 0; }

    if ( self == NULL )
        rc = RC( rcAlign, rcIterator, rcReleasing, rcSelf, rcNull );
    else
    {
        if ( self->current_pi_ref != NULL )
        {
            pi_ref_whacker( (DLNode *)self->current_pi_ref, NULL );
        }
        self->current_pi_ref = ( pi_ref * )DLListPopHead ( &self->pi_refs );
        if ( self->current_pi_ref == NULL )
        {
            rc = RC( rcAlign, rcIterator, rcAccessing, rcOffset, rcDone ); 
        }
        else
        {
            /* if the caller wants to know the ref-obj... */
            if ( refobj != NULL )
            {
                pi_entry * pie = ( pi_entry * )DLListHead( &(self->current_pi_ref->pi_entries) );
                if ( pie != NULL )
                {
                    rc = PlacementIteratorRefObj( pie->pi, refobj );
                }
                else
                {
                    rc = RC( rcAlign, rcIterator, rcAccessing, rcOffset, rcDone ); 
                }
            }
            /* if the caller wants to know the starting position on this reference...*/
            if ( first_pos != NULL )
            {
                *first_pos = self->current_pi_ref->first_pos;
            }
            if ( last_pos != NULL )
            {
                *last_pos = self->current_pi_ref->last_pos;
            }
            /* start with the first entry when looping through the placement-records... */
            self->current_pi_ref->current_entry = NULL;
        }
    }
    return rc;
}



typedef struct pi_ref_nxt_avail_pos_ctx
{
    uint32_t count;
    INSDC_coord_zero min_pos;
    INSDC_coord_len min_len;
    bool min_pos_initialized;
    rc_t rc;
} pi_ref_nxt_avail_pos_ctx;

static void CC nxt_avail_pos_cb( DLNode * n, void * data )
{
    pi_ref_nxt_avail_pos_ctx * ctx = ( pi_ref_nxt_avail_pos_ctx * ) data;
    if ( ctx->rc == 0 )
    {
        pi_entry * pie = ( pi_entry * )n;
        rc_t rc = PlacementIteratorNextAvailPos ( pie->pi, &(pie->nxt_avail.first), &(pie->nxt_avail.len) );
        if ( rc == 0 )
        {
            pie->nxt_avail.last = pie->nxt_avail.first + pie->nxt_avail.len - 1;
            ( ctx->count )++;
            if ( ctx->min_pos_initialized )
            {
                if ( pie->nxt_avail.first < ctx->min_pos )
                {
                    ctx->min_pos = pie->nxt_avail.first;
                    ctx->min_len = pie->nxt_avail.len;
                }
            }
            else
            {
                ctx->min_pos = pie->nxt_avail.first;
                ctx->min_len = pie->nxt_avail.len;
                ctx->min_pos_initialized = true;
            }
        }
        else
        {
            if ( GetRCState( rc ) != rcDone )
                ctx->rc = rc;
        }
    }
}


LIB_EXPORT rc_t CC PlacementSetIteratorNextAvailPos ( const PlacementSetIterator *cself,
    INSDC_coord_zero *pos, INSDC_coord_len *len )
{
    rc_t rc = 0;
    if ( cself == NULL )
        rc = RC( rcAlign, rcIterator, rcAccessing, rcSelf, rcNull );
    else
    {
        if ( pos == NULL )
            rc = RC( rcAlign, rcIterator, rcAccessing, rcParam, rcNull );
        else
        {
            PlacementSetIterator *self = ( PlacementSetIterator * )cself;
            if ( self->current_pi_ref == NULL )
            {
                rc = RC( rcAlign, rcIterator, rcAccessing, rcOffset, rcDone );
            }
            else
            {
                /* loop through all the pi_entry int the current_pi_ref */
                pi_ref_nxt_avail_pos_ctx ctx;
                ctx.count = 0;
                ctx.rc = 0;
                ctx.min_pos = 0;
                ctx.min_len = 0;
                ctx.min_pos_initialized = false;
                DLListForEach ( &(self->current_pi_ref->pi_entries), false, nxt_avail_pos_cb, &ctx );
                rc = ctx.rc;
                if ( ctx.count == 0 )
                {
                    rc = RC( rcAlign, rcIterator, rcAccessing, rcOffset, rcDone );
                }
                else
                {
                    *pos = ctx.min_pos;
                    if ( len != NULL )
                    {
                        *len = ctx.min_len;
                    }
                }
            } 
        }
    }
    return rc;
}


static void unlink_all_before( DLList * list, INSDC_coord_zero pos )
{
    pi_entry *res = ( pi_entry * )DLListHead( list );
    while ( res != NULL )
    {
        pi_entry *nxt = ( pi_entry * )DLNodeNext( ( DLNode * )res );
        if ( res->w.last < pos )
        {
            DLListUnlink ( list, ( DLNode * )res );
            pi_entry_whacker( ( DLNode * )res, NULL );
        }
        res = nxt;
    }
}


LIB_EXPORT rc_t CC PlacementSetIteratorNextRecordAt ( PlacementSetIterator *self,
    INSDC_coord_zero pos, const PlacementRecord **rec )
{
    rc_t rc = 0;
    pi_ref * pr;
    bool done;

    if ( rec == NULL )
        return RC( rcAlign, rcIterator, rcAccessing, rcParam, rcNull );
    *rec = NULL;
    if ( self == NULL )
        return RC( rcAlign, rcIterator, rcAccessing, rcSelf, rcNull );
    if ( self->current_pi_ref == NULL )
        /* no more reference to iterator over! the iterator is done! */
        return RC( rcAlign, rcIterator, rcAccessing, rcOffset, rcDone );

    pr = self->current_pi_ref;
    done = false;
    do
    {
        if ( pr->current_entry == NULL )
        {
            unlink_all_before( &(pr->pi_entries), pos );
            pr->current_entry = ( pi_entry * )DLListHead( &(pr->pi_entries) );
        }
        done = ( pr->current_entry == NULL );
        rc = ( done ? RC( rcAlign, rcIterator, rcAccessing, rcOffset, rcDone ) : 0 );

        if ( rc == 0 )
        {
            rc = PlacementIteratorNextRecordAt ( pr->current_entry->pi, pos, rec );
            done = ( GetRCState( rc ) != rcDone );
            if ( !done )
            {
                pr->current_entry = ( pi_entry * )DLNodeNext( ( DLNode * )pr->current_entry );
                done = ( pr->current_entry == NULL );
                rc = ( done ? RC( rcAlign, rcIterator, rcAccessing, rcOffset, rcDone ) : 0 );
            }
        }
    } while ( !done );
    return rc;
}


LIB_EXPORT rc_t CC PlacementSetIteratorNextIdAt ( PlacementSetIterator *self,
    INSDC_coord_zero pos, int64_t *row_id, INSDC_coord_len *len )
{
    rc_t rc = 0;
    if ( self == NULL )
        rc = RC( rcAlign, rcIterator, rcAccessing, rcSelf, rcNull );
    else
    {

    }
    return rc;
}
