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
#include <klib/container.h>
#include <klib/refcount.h>
#include <klib/sort.h>
#include <klib/text.h>
#include <klib/out.h>
#include <insdc/insdc.h>
#include <align/reference.h>
#include <align/iterator.h>
#include <align/manager.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>


struct ReferenceIterator
{
    KRefcount refcount;
    struct AlignMgr const *amgr;            /* the alignment-manager... */

    DLList records;                         /* has list of records... */
    int32_t min_mapq;                       /* has a minimum mapq-value... */
    PlacementRecordExtendFuncs ext_func;    /* has a struct with record-extension-functions from client*/
    PlacementRecordExtendFuncs int_func;    /* has a struct with record-extension-functions for itself*/

    uint32_t depth;                         /* how many records are in the list */
    INSDC_coord_zero current_pos;           /* what is the current ref-position on the current ref. */
    INSDC_coord_zero last_pos;              /* what is the last ref-position on the current ref. */
    INSDC_coord_zero nxt_avail_pos;         /* what is the next available ref-position on the current ref. */
    PlacementRecord *current_rec;           /* the current-record at the current position */
    bool last_rec_reached;                  /* do we have reached the last record at a given position ? */
    bool need_init;                         /* do we need to init for the first next()-call */
    PlacementSetIterator * pl_set_iter;     /* holds a list of placement-iterators */
    struct ReferenceObj const * refobj;     /* cached result of ReferenceIteratorNextReference(...) */
};


LIB_EXPORT void CC RefIterRecordDestroy ( void *obj, void *data )
{
    /* nothing to do, because there are no sub-allocations etc. ... */
}


LIB_EXPORT rc_t CC RefIterRecordPopulate ( void *obj,
    const PlacementRecord *placement, struct VCursor const *curs,
    INSDC_coord_zero ref_window_start, INSDC_coord_len ref_window_len, void *data )
{
    /* read the data required to build a Alignment-Iterator,
       then create the Alignment-Iterator into the already allocated memory */
    return AlignIteratorRecordPopulate ( obj, placement, curs, ref_window_start, ref_window_len, data );
}


LIB_EXPORT size_t CC RefIterRecordSize ( struct VCursor const *curs,
    int64_t row_id, void *data )
{
    /* discover the size of the ref-iter-part to be allocated... */
    return AlignIteratorRecordSize ( curs, row_id, data );
}


static void CC RefIterDestroyRecPart( void *obj, void *data )
{
    AlignmentIterator *iter = ( AlignmentIterator * )obj;
    if ( iter != NULL )
        AlignmentIteratorRelease( iter );
}

LIB_EXPORT rc_t CC AlignMgrMakeReferenceIterator ( struct AlignMgr const *self,
    ReferenceIterator **iter, const PlacementRecordExtendFuncs *ext_1, int32_t min_mapq )
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
            ReferenceIterator * refi = calloc( sizeof * refi, 1 );
            if ( refi == NULL )
                rc = RC( rcAlign, rcIterator, rcConstructing, rcMemory, rcExhausted );
            else
            {
                KRefcountInit( &refi->refcount, 1, "ReferenceIterator", "Make", "align" );
                refi->min_mapq = min_mapq;
                if ( ext_1 != NULL )
                {
                    refi->ext_func.data = ext_1->data;
                    refi->ext_func.destroy = ext_1->destroy;
                    refi->ext_func.populate = ext_1->populate;
                    refi->ext_func.alloc_size = ext_1->alloc_size;
                    refi->ext_func.fixed_size = ext_1->fixed_size;
                }

                refi->int_func.data = ( void * )self;     /* have it point to the AlignMgr... */
                refi->int_func.destroy = RefIterDestroyRecPart;
                refi->int_func.populate = RefIterRecordPopulate;
                refi->int_func.alloc_size = RefIterRecordSize; 

                DLListInit( &(refi->records) );
                rc = AlignMgrMakePlacementSetIterator ( self, &refi->pl_set_iter );
                refi->need_init = true;
            }

            if ( rc == 0 )
                rc = AlignMgrAddRef ( self );

            if ( rc == 0 )
            {
                refi->amgr = self;
                *iter = refi;
            }
            else
                free( refi );
        }
    }
    return rc;
}


LIB_EXPORT rc_t CC ReferenceIteratorAddRef ( const ReferenceIterator *self )
{
    rc_t rc = 0;
    if ( self == NULL )
        rc = RC( rcAlign, rcIterator, rcAttaching, rcSelf, rcNull );
    else
    {
        if ( KRefcountAdd( &self->refcount, "ReferenceIterator" ) != krefOkay )
        {
            rc = RC( rcAlign, rcIterator, rcAttaching, rcError, rcUnexpected );
        }
    }
    return rc;
}


static void CC whack_the_placement_record( DLNode *n, void *data )
{
    PlacementRecord * rec = ( PlacementRecord * )n;
    PlacementRecordWhack ( rec );
}


static void clear_recordlist( DLList * list )
{
    DLListWhack ( list, whack_the_placement_record, NULL );
}

LIB_EXPORT rc_t CC ReferenceIteratorRelease ( const ReferenceIterator *cself )
{
    rc_t rc = 0;
    if ( cself == NULL )
        rc = RC( rcAlign, rcIterator, rcReleasing, rcSelf, rcNull );
    else
    {
        if ( KRefcountDrop( &cself->refcount, "ReferenceIterator" ) == krefWhack )
        {
            ReferenceIterator * self = ( ReferenceIterator * ) cself;
            /* we 'own' the records! - we have to destroy them, if some are left in here */
            clear_recordlist( &self->records );
            rc = PlacementSetIteratorRelease ( self->pl_set_iter );
            AlignMgrRelease ( self->amgr );
            free( self );
        }
    }
    return rc;
}


LIB_EXPORT rc_t CC ReferenceIteratorAddPlacementIterator( ReferenceIterator *self,
    PlacementIterator *pi )
{
    rc_t rc = 0;
    if ( self == NULL )
        rc = RC( rcAlign, rcIterator, rcConstructing, rcSelf, rcNull );
    else
    {
        if ( pi == NULL )
            rc = RC( rcAlign, rcIterator, rcConstructing, rcParam, rcNull );
        else
        {
            rc = PlacementSetIteratorAddPlacementIterator ( self->pl_set_iter, pi );
        }
    }
    return rc;
}


#define ALIGN_COL_COUNT 4

static const char * align_cols[ ALIGN_COL_COUNT ] = 
{ "(I32)CLIPPED_REF_OFFSET",
  "(bool)CLIPPED_HAS_REF_OFFSET",
  "(bool)CLIPPED_HAS_MISMATCH",
  "(INSDC:dna:text)CLIPPED_READ" };


static rc_t prepare_align_cursor( struct VCursor const *align )
{
    rc_t rc = 0;
    uint32_t i, throw_away_idx;

    for ( i = 0; i < ALIGN_COL_COUNT && rc == 0; ++i )
        rc =VCursorAddColumn ( align, &throw_away_idx, "%s", align_cols[ i ] );
    return rc;
}


LIB_EXPORT rc_t CC ReferenceIteratorAddPlacements( ReferenceIterator *self,
     struct ReferenceObj const *ref_obj, INSDC_coord_zero ref_pos, INSDC_coord_len ref_len,
     struct VCursor const *ref, struct VCursor const *align, align_id_src ids )
{
    rc_t rc = 0;
    if ( self == NULL )
        rc = RC( rcAlign, rcIterator, rcConstructing, rcSelf, rcNull );
    else
    {
        if ( ref_obj == NULL )
            rc = RC( rcAlign, rcIterator, rcConstructing, rcParam, rcNull );
        else
        {
            if ( align != NULL )
                rc = prepare_align_cursor( align );

            if ( rc == 0 )
            {
                PlacementIterator *pi;

                rc = ReferenceObj_MakePlacementIterator ( ref_obj, &pi, ref_pos, ref_len, self->min_mapq,
                        ref, align, ids, &self->int_func, &self->ext_func );
                if ( rc == 0 )
                {
                    rc = PlacementSetIteratorAddPlacementIterator ( self->pl_set_iter, pi );
                }
            }
            ReferenceObj_Release( ref_obj );
        }
    }
    return rc;
}


static rc_t fill_recordlist( ReferenceIterator *self, INSDC_coord_zero pos )
{
    rc_t rc = 0;
    while ( rc == 0 )
    {
        const PlacementRecord *rec;
        /* from the placement-set-iterator into our list... */
        rc = PlacementSetIteratorNextRecordAt ( self->pl_set_iter, pos, &rec );
        if ( rc == 0 )
        {
            if ( rec->pos == pos )
            {
                self->depth++;
                DLListPushTail ( &self->records, (DLNode *)rec );
/*                OUTMSG(( " add [%lu]%u(%u)", rec->id, rec->pos, rec->len )); */
            }
            else
            {
                PlacementRecordWhack ( rec );
            }
        }
    }
    if ( GetRCState( rc ) == rcDone ) rc = 0;
    return rc;
}


static uint32_t remove_invalid_records( DLList * list, INSDC_coord_zero pos )
{
    uint32_t res = 0;
    PlacementRecord *rec = ( PlacementRecord * )DLListHead( list );
    while ( rec != NULL )
    {
        PlacementRecord *nxt = ( PlacementRecord * )DLNodeNext( ( DLNode * )rec );
        bool remove = ( ( rec->pos + rec->len ) <= pos );
        if ( !remove )
        {
            AlignmentIterator * al_iter = PlacementRecordCast ( rec, placementRecordExtension0 );
            int32_t state = AlignmentIteratorState ( al_iter, NULL );
            remove = ( ( state & align_iter_invalid ) == align_iter_invalid );
        }
        if ( remove )
        {
            DLListUnlink ( list, ( DLNode * )rec );
            PlacementRecordWhack ( rec );
        }
        else
        {
            res++;
        }
        rec = nxt;
    }
    return res;
}


static void inc_alignment_iterators( ReferenceIterator *self, INSDC_coord_zero pos )
{
    PlacementRecord *rec = ( PlacementRecord * )DLListHead( &self->records );
    while ( rec != NULL )
    {
        AlignmentIterator * al_iter = PlacementRecordCast ( rec, placementRecordExtension0 );
        if ( rec->pos <= pos && al_iter != NULL )
        {
            AlignmentIteratorNext ( al_iter );
        }
        rec = ( PlacementRecord * )DLNodeNext( ( DLNode * )rec );
    }
}


LIB_EXPORT rc_t CC ReferenceIteratorNextReference ( ReferenceIterator *self,
    struct ReferenceObj const ** refobj )
{
    rc_t rc = 0;
    if ( self == NULL )
        rc = RC( rcAlign, rcIterator, rcAccessing, rcSelf, rcNull );
    else
    {
        struct ReferenceObj const * robj;
        rc = PlacementSetIteratorNextReference ( self->pl_set_iter,
            &self->current_pos, &self->last_pos, &robj );
        clear_recordlist( &self->records );
        if ( rc == 0 )
        {
            /* cache the returned refobj in order to get to reference-bases later... */
            self->refobj = robj;
            self->need_init = true;
        }
        else
        {
            self->refobj = NULL;
        }
        if ( refobj != NULL )
        {
            *refobj = self->refobj;
        }
    }
    return rc;
}


static rc_t first_ref_iter_nxt_pos( ReferenceIterator *self, bool skip_empty )
{
    rc_t rc = 0;
    uint32_t diff;
    self->need_init = false;

    do
    {
        rc = PlacementSetIteratorNextAvailPos ( self->pl_set_iter, &self->nxt_avail_pos, NULL );
        if ( ( rc == 0 ) && ( self->nxt_avail_pos <= self->current_pos ) )
        {
            rc = fill_recordlist( self, self->nxt_avail_pos );
            diff = ( self->nxt_avail_pos - diff );
        }
    } while ( ( self->nxt_avail_pos < self->current_pos ) && ( rc == 0 ) && diff != 0 );

    /* jump over gaps, if requested ... */
    if ( skip_empty && ( self->current_pos < self->nxt_avail_pos ) && self->depth == 0 )
    {
        self->current_pos = self->nxt_avail_pos;
    }

    if ( GetRCState( rc ) == rcDone ) rc = 0;
    self->last_rec_reached = false;
    self->depth = remove_invalid_records( &self->records, self->current_pos );
    return rc;
}

LIB_EXPORT rc_t CC ReferenceIteratorNextPos ( ReferenceIterator *self, bool skip_empty )
{
    rc_t rc = 0;
    if ( self == NULL )
        rc = RC( rcAlign, rcIterator, rcAccessing, rcSelf, rcNull );
    else
    {
        if ( self->need_init )
        {
            rc = first_ref_iter_nxt_pos( self, skip_empty );
        }
        else
        {
            /* increment the current position */
            self->current_pos++;

            if ( self->current_pos <= self->last_pos )
            {
                /* jump over gaps, if requested ... */
                if ( self->depth == 0 && skip_empty )
                {
                    self->current_pos = self->nxt_avail_pos;
                }

                /* increment the internal alignment-iterator of every placement-record */
                inc_alignment_iterators( self, self->current_pos );

                /* loop through the list to look if we have to remove records,
                   that do end before this new position */
                self->depth = remove_invalid_records( &self->records, self->current_pos );

                rc = fill_recordlist( self, self->current_pos );
                if ( rc == 0 )
                {
                    /* set our sights to the next position... */
                    rc = PlacementSetIteratorNextAvailPos ( self->pl_set_iter, &self->nxt_avail_pos, NULL );
                    if ( GetRCState( rc ) == rcDone )
                    {
                        if ( self->depth > 0 )
                        {
                            rc = 0;
                        }
                        else if ( !skip_empty )
                        {
                            if ( self->current_pos <= self->last_pos ) rc = 0;
                        }
                    }
                }
                self->last_rec_reached = false;
            }
            else
            {
                rc = RC( rcAlign, rcIterator, rcAccessing, rcOffset, rcDone );
                self->last_rec_reached = true;
                clear_recordlist( &self->records );
            }
        }
        /* load the current record with the first record out of our list */
        /* self->current_rec = ( PlacementRecord * )DLListHead( &self->records ); */
        self->current_rec = NULL;
    }
    return rc;
}


LIB_EXPORT rc_t CC ReferenceIteratorPosition ( const ReferenceIterator *self,
    INSDC_coord_zero *pos, uint32_t * depth, INSDC_4na_bin * base )
{
    rc_t rc = 0;
    if ( self == NULL )
        rc = RC( rcAlign, rcIterator, rcAccessing, rcSelf, rcNull );
    else
    {
        /* return position, many records our record-list, and the base at this position */
        if ( pos != NULL )
        {
            *pos = self->current_pos;
        }

        if ( depth != NULL )
        {
            *depth = self->depth;
        }

        if ( base != NULL )
        {
            uint32_t written;
            *base = 0;
            /* problem! how to get the base if depth == 0 */
            if ( self->current_rec != NULL )
            {
                rc = ReferenceObj_Read( self->current_rec->ref, self->current_pos, 1, base, &written );
            }
            else if ( self->refobj != NULL )
            {
                rc = ReferenceObj_Read( self->refobj, self->current_pos, 1, base, &written );
            }
        }
    }
    return rc;
}


LIB_EXPORT rc_t CC ReferenceIteratorNextPlacement ( ReferenceIterator *self,
    const PlacementRecord **rec )
{
    rc_t rc = 0;
    if ( self == NULL )
        rc = RC( rcAlign, rcIterator, rcAccessing, rcSelf, rcNull );
    else
    {
        if ( self->current_rec == NULL )
        {
            if ( !self->last_rec_reached )
            {
                self->current_rec = ( PlacementRecord * )DLListHead( &self->records );
            }
        }
        else
        {
            self->current_rec = ( PlacementRecord * )DLNodeNext( ( DLNode * )self->current_rec );
        }

        if ( rec != NULL ) *rec = self->current_rec;

        if ( self->current_rec == NULL )
        {
            rc = RC( rcAlign, rcIterator, rcAccessing, rcOffset, rcDone );
            self->last_rec_reached = true;
        }
    }
    return rc;
}


LIB_EXPORT int32_t CC ReferenceIteratorState ( const ReferenceIterator *self,
    INSDC_coord_zero *seq_pos )
{
    int32_t res = align_iter_invalid;
    if ( seq_pos != NULL )
    {
        *seq_pos = 0;
    }
    if ( self != NULL )
    {
        /* PlacementRecordCast returns NULL if self->current_rec is NULL */
        AlignmentIterator * al_iter = PlacementRecordCast ( self->current_rec, placementRecordExtension0 );
        if ( al_iter != NULL )
            res = AlignmentIteratorState ( al_iter, seq_pos );
    }
    return res;
}


LIB_EXPORT uint32_t CC ReferenceIteratorBasesInserted ( const ReferenceIterator *self,
    const INSDC_4na_bin **bases )
{
    uint32_t res = align_iter_invalid;
    if ( bases != NULL )
    {
        *bases = NULL;
    }
    if ( self != NULL )
    {
        /* PlacementRecordCast returns NULL if self->current_rec is NULL */
        AlignmentIterator * al_iter = PlacementRecordCast ( self->current_rec, placementRecordExtension0 );
        if ( al_iter != NULL )
            res = AlignmentIteratorBasesInserted( al_iter, bases );
    }
    return res;
}


LIB_EXPORT uint32_t CC ReferenceIteratorBasesDeleted ( const ReferenceIterator *self,
    INSDC_coord_zero *pos, const INSDC_4na_bin **bases )
{
    uint32_t res = align_iter_invalid;
    if ( bases != NULL )
    {
        *bases = NULL;
    }
    if ( self != NULL )
    {
        /* PlacementRecordCast returns NULL if self->current_rec is NULL */
        AlignmentIterator * al_iter = PlacementRecordCast ( self->current_rec, placementRecordExtension0 );
        if ( al_iter != NULL )
        {
            INSDC_coord_zero temp_pos;
            res = AlignmentIteratorBasesDeleted( al_iter, &temp_pos );
            if ( ( res & align_iter_invalid ) != align_iter_invalid )
            {
                if ( pos != NULL ) { *pos = temp_pos; }
                /* where to get the reference-bases from ? PlacementRecord.ref ! */
                if ( res > 0 && bases != NULL )
                {
                    uint8_t * buffer = malloc( res );
                    if ( buffer != NULL )
                    {
                        INSDC_coord_len written;
                        rc_t rc = ReferenceObj_Read( self->current_rec->ref, temp_pos, res, buffer, &written );
                        if ( rc == 0 )
                        {
                            *bases = buffer;
                        }
                    }
                }
            }
        }
    }
    return res;
}
