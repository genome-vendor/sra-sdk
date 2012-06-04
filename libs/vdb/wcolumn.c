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

#include <vdb/extern.h>

#include "klib/symbol.h"
#include "column-priv.h"
#include "dbmgr-priv.h"
#include "schema-priv.h"
#include "schema-expr.h"
#include "schema-parse.h"
#include "cursor-priv.h"
#include "prod-priv.h"
#include "blob-priv.h"
#include "page-map.h"

#include <vdb/manager.h>
#include <vdb/cursor.h>
#include <kdb/column.h>
#include <klib/log.h>
#include <klib/rc.h>
#include <bitstr.h>
#include <sysalloc.h>

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define DFLT_TRIGGER ( 2 * 1024 * 1024 )
#define USE_RLE 1

/*--------------------------------------------------------------------------
 * VDBMemBuffer
 *  manages a linked list of VDBMem objects
 */
typedef struct VDBMemBuffer VDBMemBuffer;
struct VDBMemBuffer
{
    bitsz_t start;      /* offset in bits from the first memory page start to the first unwritten element */
    bitsz_t end;        /* bit offset just past last unwritten element */
    bitsz_t marker;     /* a cursor of sorts that should be in the page referred to by mem */
    /* NOTE: marker is not always referring to an element and the user of this object
     * is required to track this */
    const VDBMem *mem;  /* a memory page that should contain the marker below */
};

/* these are compile time constants */
#define VDBMemBitsInPage()      ((bitsz_t)(sizeof ((VDBMem*)0)->page * 8))
#define VDBMemBitsInPageMask()  (VDBMemBitsInPage()-1)

/* InitHead
 * initialize buffer to point to left edge of first element
 * buffer may be empty
 *
 * Parameters are expected to be good, no production run-time verification applied.
 */
static
void VDBMemBufferInitHead ( VDBMemBuffer *self, void *mem,
    bitsz_t start, bitsz_t elem_size, uint64_t elem_count )
{
    /* during development these are validated: element size can't ever be zero
     * the memory page must be an actual page or if its not then there can't be
     * any elements or offset with in the NULL page */
    assert (self != NULL);
    assert (elem_size != 0);
    assert ((mem != NULL) || ((mem == NULL) && (start == 0) && (elem_count == 0)));

    self -> mem = ( const void* ) mem;
    self -> start = self -> marker = start;     /* cursor at first element */
    self -> end = start + elem_size * elem_count; /* just past last element */
}

/* InitTail
 *  initialize buffer to point to last bit
 *  buffer must not be empty
 *
 * Parameters are expected to be good, no verification applied.
 */
static
void VDBMemBufferInitTail ( VDBMemBuffer *self, void *mem,
    bitsz_t start, bitsz_t elem_size, uint64_t elem_count )
{
    /* during development these are validated: element size can't ever be zero
     * the memory page must be an actual page and there must be valid elements */
    assert (self != NULL);
    assert (elem_size != 0);
    assert ( mem != NULL );
    assert (elem_count != 0);
    assert (start < VDBMemBitsInPage());

    self -> mem = ( const void* ) mem;
    self -> start = start;
    self -> end = start + elem_size * elem_count;
    assert ( self -> end > start );
    self -> marker = self -> end - 1;   /* cursor at last element's last bit */
    /* NOTE: with marker not referring to an element the first "advance" must
     * take that off by one bit into account rather than seek the exact number of bits
     * in elements that we'd want to move the marker */
}

/* Advance
 *  sets marker
 */
static
bool VDBMemBufferAdvance ( VDBMemBuffer *self, bitsz_t bits, bool reverse )
{
    if ( bits != 0 )
    {
        bitsz_t marker;
        const VDBMem *mem;
        uint64_t cur, targ;

        if ( ! reverse )                        /* forward */
            marker = self -> marker + bits;
        else                                    /* backward */
        {
            if ( self -> marker < bits )        /* can't back up past beginning of buffer list */
                return false;
            marker = self -> marker - bits;
            if ( marker < self -> start )       /* can't back up past stored begining of valid data */
                return false;
        }

        /* anything beyond buffer is an illegal seek */
        if ( marker >= self -> end )
            return false;

        /* convert markers to zero-based page ids */
        cur = self -> marker / ( sizeof mem -> page * 8 );
        targ = marker / ( sizeof mem -> page * 8 );
        if ( cur != targ )      /* change pages if we have to */
        {
            /* rewind loop */
            for ( mem = self -> mem; cur > targ; -- cur )
            {
                mem = ( const VDBMem* ) DLNodePrev ( & mem -> n );
                assert ( mem != NULL );
            }

            /* forward loop */
            for ( ; cur < targ; ++ cur )
            {
                mem = ( const VDBMem* ) DLNodeNext ( & mem -> n );
                assert ( mem != NULL );
            }

            self -> mem = mem;
        }

        /* take new marker */
        self -> marker = marker;
    }

    return true;
}

/* Access
 *  return pointer to current element (by marker/cursor)
 *  along with bit offset and number of BITS available
 */
static
const void *VDBMemBufferAccess ( VDBMemBuffer const *self, bitsz_t *boff, bitsz_t *avail )
{
    const VDBMem *mem = self -> mem;
    bitsz_t bits, marker = self -> marker;

    /* if nothing was ever added to the buffer list */
    assert (marker < self->end);
    if ( mem == NULL )
    {
        if ( boff != NULL )
            * boff = 0;
        if ( avail != NULL )
            * avail = 0;
        return NULL;
    }

    /* bits available in page past marker */
    bits =  VDBMemBitsInPage() - ( marker & VDBMemBitsInPageMask()) ;

    /* limit to end of buffer if not full page */
    if ( marker + bits > self -> end )
        bits = self -> end - marker;

    /* return bit offset */
    if ( boff != NULL )
        * boff = marker & 7;

    /* return bytes available */
    if ( avail != NULL )
        * avail = bits;

    /* return byte pointer */
/*     return & mem -> page [ ( marker >> 3 ) & ( sizeof mem -> page - 1 ) ]; */
    return & mem->page[(marker & VDBMemBitsInPageMask()) >> 3];
}



/*--------------------------------------------------------------------------
 * VColumn
 */

static 
void WColumnDestroy (WColumn * self)
{
#if PROD_REFCOUNT && ! PROD_ALL_IN_CURSOR
    PROD_TRACK_REFCOUNT (VProductionRelease,self->out);
    VProductionRelease ( self->out, NULL);
    PROD_TRACK_REFCOUNT (VProductionRelease,self->val);
    VProductionRelease (self->val, NULL);
#endif
}
/* Whack
 *  perform read-only cleanup
 */
void CC VColumnWhack ( void *item, void *data )
{
    VColumn *self = item;
    VCursor *curs = data;

    /* remove from cursor */
    if ( curs != NULL )
    {
        VectorSwap ( & curs -> row, self -> ord, NULL, & item );
        VCursorCacheSwap ( & curs -> col, & self -> scol -> cid, NULL, & item );
    }

    if ( ! self -> read_only )
    {
        WColumn *wself = (WColumn *)self;

        if ( wself -> page != NULL )
        {
            TRACK_BLOB (VBlobRelease, wself->page);
            (void)VBlobRelease ( wself -> page );
        }

        KDataBufferWhack ( & wself -> dflt );
        DLListWhack ( & wself -> data, VDBMemRelease, wself -> mgr );
        DLListWhack ( & wself -> rowmap, VDBMemRelease, wself -> mgr );
        VDBManagerSever ( wself -> mgr );
        WColumnDestroy (wself);

    }

    VColumnDestroy ( self );
    free ( self );
}


/* Make - PRIVATE
 *  make a write column
 */
rc_t WColumnMake ( VColumn **colp, const VSchema *schema,
    const SColumn *scol, const SExpression *blob_limit, VDBManager *mgr )
{
    rc_t rc;
    WColumn *col;

    assert ( colp != NULL );
    assert ( schema != NULL );
    assert ( scol != NULL );
    assert ( mgr != NULL );

    col = malloc ( sizeof * col );
    if ( col == NULL )
        rc = RC ( rcVDB, rcColumn, rcConstructing, rcMemory, rcExhausted );
    else
    {
        memset ( col, 0, sizeof * col );
        rc = VColumnInit ( & col -> dad, schema, scol );
        if ( rc == 0 )
        {
            col -> mgr = VDBManagerAttach ( mgr );
            DLListInit ( & col -> data );
            DLListInit ( & col -> rowmap );

            if ( scol -> limit != NULL )
                blob_limit = scol -> limit;

            if ( blob_limit == NULL )
            {
#ifdef DFLT_TRIGGER
                col -> trigger = DFLT_TRIGGER;
#else
                -- col -> trigger;
#endif
            }
            else
            {
                /* evaluate column blob limit */
                uint64_t trigger;
                rc = eval_uint64_expr ( schema, blob_limit, & trigger );
                col -> trigger = ( size_t ) trigger;
            }

            if ( rc == 0 )
            {
                * colp = & col -> dad;
                return 0;
            }
        }

        free ( col );
    }

    * colp = NULL;
    return rc;
}


/* IdRange
 *  returns id range for column or page
 */
rc_t VColumnIdRange ( const VColumn *vcol, int64_t *first, int64_t *last )
{
    rc_t rc;
    const WColumn *self = ( const WColumn* ) vcol;

    assert ( self != NULL );
    assert ( first != NULL && last != NULL );

    if ( self -> dad . in != NULL )
        return VColumnIdRangeRead ( & self -> dad, first, last );

    if ( self -> val == NULL )
        rc = RC ( rcVDB, rcColumn, rcAccessing, rcRange, rcNotOpen );
    else
    {
        /* a little silly, but set max range in 64-bit
           without complaints from 32-bit compilers */
        * first = 1;
        * first <<= 63;
        * last = ~ * first;

        /* now intersect this range with all physical sources */
        rc = VProductionColumnIdRange ( self -> val, first, last );
        if ( rc == 0 )
            return 0;
    }

    * first = * last = 0;

    return rc;
}

/* SetDefault
 *  capture default row data
 */
rc_t WColumnSetDefault ( VColumn *vcol,
    bitsz_t elem_bits, const void *buffer, bitsz_t boff, uint64_t len )
{
    rc_t rc;
    bitsz_t elem_size;
    WColumn *self = ( WColumn* ) vcol;

    assert ( elem_bits != 0 );
    assert ( buffer != NULL || ( boff == 0 && len == 0 ) );

    /* test "compatibility" of elem_bits
       this is used to interpret "len" */
    elem_size = VTypedescSizeof ( & self -> dad . desc );
    if ( elem_bits < elem_size && elem_size % elem_bits != 0 )
        return RC ( rcVDB, rcColumn, rcUpdating, rcType, rcInconsistent );
    if ( elem_bits > elem_size && elem_bits % elem_size != 0 )
        return RC ( rcVDB, rcColumn, rcUpdating, rcType, rcInconsistent );

    /* allow NULL setting */
    if ( buffer == NULL )
    {
        KDataBufferWhack ( & self -> dflt );
        memset ( & self -> dflt, 0, sizeof self -> dflt );
        self -> have_dflt = true;
        return 0;
    }

    /* set the element size */
    rc = KDataBufferCast ( & self -> dflt, & self -> dflt, elem_bits, false );
    if ( rc != 0 )
        return rc;

    /* set the length */
    rc = KDataBufferResize ( & self -> dflt, len );
    if ( rc != 0 )
    {
        assert ( KDataBufferWritable ( & self -> dflt ) );
        return rc;
    }

    /* copy in the row */
    bitcpy ( self -> dflt . base, 0, buffer, boff, len * elem_bits );
    self -> have_dflt = true;
    return 0;
}

/* OpenRow
 *  update state
 *
 *  "const_row_id" [ IN, CONST ] - id of row being opened. useful
 *  only on initial open when no other rows are buffered.
 */
void CC WColumnOpenRow ( void *item, void *const_row_id )
{
    WColumn *self = item;
    int64_t row_id = * ( const int64_t* ) const_row_id;

    assert ( ! self -> row_written );
    if ( self -> start_id != self -> end_id )
    {
        assert ( row_id == self -> end_id );
    }
    else
    {
        /* capture row id */
        self -> start_id = self -> end_id = self -> cutoff_id = row_id;
        assert ( self -> data_off == 0 );
        assert ( self -> rowmap_off == 0 );
        assert ( self -> num_rows == 0 );
        assert ( self -> num_elems == 0 );
        assert ( self -> row_len == 0 );
    }
}


/* Write
 */
rc_t WColumnWrite ( VColumn *cself,
    bitsz_t elem_bits, const void *buffer, bitsz_t boff, uint64_t len )
{
    rc_t rc;
    bitsz_t elem_size;
    WColumn *self = ( WColumn* ) cself;

    DLList data;
    VDBMem *mem;
    uint8_t *dst;
    bitsz_t num_bits, num_writ, to_write, doff;

    assert ( elem_bits != 0 );
    assert ( buffer != NULL || ( boff == 0 && len == 0 ) );

    /* test "compatibility" of elem_bits
       this is used to interpret "len" */
    elem_size = VTypedescSizeof ( & self -> dad . desc );
    if ( elem_bits < elem_size && elem_size % elem_bits != 0 )
        return RC ( rcVDB, rcColumn, rcUpdating, rcType, rcInconsistent );
    if ( elem_bits > elem_size && elem_bits % elem_size != 0 )
        return RC ( rcVDB, rcColumn, rcUpdating, rcType, rcInconsistent );

    /* allow empty row */
    if ( len == 0 )
    {
        self -> row_written = true;
        return 0;
    }

    /* disallow any further modifications */
    if ( self -> row_committed )
        return RC ( rcVDB, rcColumn, rcUpdating, rcColumn, rcBusy );

    /* the number of bits to write */
    num_bits = ( bitsz_t ) elem_bits * len;

    /* prepare data row pointer */
    doff = self -> data_off + ( self -> num_elems + self -> row_len ) * elem_size;
    doff &= sizeof mem -> page * 8 - 1;

    /* prepare whack lists */
    DLListInit ( & data );

    /* prepare initial destination pointer */
    mem = ( VDBMem* ) DLListTail ( & self -> data );
    dst = ( mem == NULL ) ? NULL : mem -> page;

    /* append row data */
    for ( rc = 0, num_writ = 0; num_writ < num_bits; doff += to_write, num_writ += to_write )
    {
        bitsz_t avail;

        /* allocate more buffer space */
        if ( ( doff & ( sizeof mem -> page * 8 - 1 ) ) == 0 )
        {
            rc = VDBManagerMakeMem ( self -> mgr, & mem );
            if ( rc != 0 )
                break;
            DLListPushTail ( & data, & mem -> n );
            dst = mem -> page;
            doff = 0;
        }

        /* decide on the bits to write */
        avail = ( sizeof mem -> page * 8 ) - doff;
        to_write = num_bits - num_writ;
        if ( avail < to_write )
            to_write = avail;

        /* copy in the data */
        bitcpy ( dst, doff, buffer, boff + num_writ, to_write );
    }

    /* if all were written, accept changes */
    if ( rc == 0 )
    {
        self -> row_len += num_bits / elem_size;
        DLListAppendList ( & self -> data, & data );
        self -> row_written = true;
    }

    /* drop any uncommitted buffers */
    DLListWhack ( & data, VDBMemRelease, self -> mgr );

    return rc;
}

/* RowDefaults
 *  if a row has not been written but has a default value,
 *  that value is written to the row. if no default exists,
 *  an error is generated.
 *
 *  "rc" [ OUT, DEFAULT ZERO ] - preset to 0
 *
 *  returns true if any error occurs ( i.e. "*rc != 0" )
 */
bool CC WColumnRowDefaults ( void *item, void *data )
{
    WColumn *self = item;
    rc_t *rc = data;

    /* nothing to do if row written */
    if ( self -> row_written )
        return false;

    /* error if no default value */
    if ( ! self -> have_dflt )
    {
        * rc = RC ( rcVDB, rcColumn, rcClosing, rcRow, rcIncomplete );
        (void)PLOGERR(klogErr, (klogErr, *rc, "Column: $(col)", "col=%.*s", self->dad.scol->name->name.size, self->dad.scol->name->name.addr));
        return true;
    }

    /* detect NULL row as default */
    if ( self -> dflt . elem_bits == 0 )
    {
        /* need to clip here, i.e. need to cut off blobs */
        * rc = RC ( rcVDB, rcColumn, rcClosing, rcRow, rcNull );
        (void)PLOGERR(klogWarn, (klogWarn, *rc, "Column: $(col)", "col=%.*s", self->dad.scol->name->name.size, self->dad.scol->name->name.addr));
        return false;
    }

    /* write default data */
    * rc = WColumnWrite ( & self -> dad, self -> dflt . elem_bits,
        self -> dflt . base, 0, self -> dflt . elem_count );

    return ( * rc != 0 ) ? true : false;
}

/* CommitRow
 *  closes the row to further writes and accepts
 *  all data written so far as complete. if the accumulated
 *  page data trigger a flush, the flush parameter is set.
 *
 *  "end_id" [ IN/OUT ] - used to calculate the minimum
 *  end_id for pages. if the column decides that it has too
 *  much data in its buffer and wants a cutoff < current
 *  value, it can lower the id.
 *
 *  returns true if there was a memory error.
 */
bool CC WColumnCommitRow ( void *item, void *data )
{
    WColumn *self = item;
    int64_t *end_id = data;

    VDBMem *mem;
    bitsz_t elem_bits;
    WColumnRowMap *rm;
    size_t rmoff, cur_size;

    /* if no data were written and that's okay, ignore */
    if ( ! self -> row_written )
    {
        self -> row_committed = true;
        return false;
    }

    /* likely to need this later */
    elem_bits = VTypedescSizeof ( & self -> dad . desc );

    /* last buffer */
    mem = ( VDBMem* ) DLListTail ( & self -> rowmap );

    /* detect a prior row */
    if ( self -> num_rows != 0 )
    {
        /* byte offset to rowmap entry */
        rmoff = ( self -> num_rows + self -> rowmap_off - 1 ) * sizeof * rm;

        assert ( mem != NULL );

        /* point to prior row entry */
        rm = ( void* ) & mem -> page [ rmoff & ( sizeof mem -> page - 1 ) ];

        /* RLE comparison starts with lengths */
        if ( USE_RLE && rm -> len == self -> row_len )
        {
            /* two consecutive rows with same length */
            if ( self -> row_len == 0 )
            {
                /* when the length is zero, there's nothing much to do */
                ++ rm -> cnt;
                self -> row_committed = true;
                return false;
            }
            else
            {
                VDBMemBuffer cur, prior;
                bitsz_t cnt, to_cmp, num_bits = ( bitsz_t ) elem_bits * self -> row_len;

                /* create first view onto buffer at tail */
                VDBMemBufferInitTail ( & cur, DLListTail ( & self -> data ),
                    self -> data_off, elem_bits, self -> num_elems + self -> row_len );

                /* back up to beginning of uncommitted row */
                if ( ! VDBMemBufferAdvance ( & cur, num_bits - 1, true ) )
                    return true;

                /* "cur" points at the current row. prepare "prior" to point to
                   the row before that which is the last committed row, which we
                   know is the same length as the current row */
                prior = cur;
                assert ( self -> num_elems >= self -> row_len );
                if ( ! VDBMemBufferAdvance ( & prior, num_bits, true ) )
                    return true;

                /* comparison loop */
                for ( cnt = to_cmp = 0; cnt < num_bits; cnt += to_cmp )
                {
                    const void *p, *c;
                    bitsz_t coff, cavail, poff, pavail;

                    VDBMemBufferAdvance ( & prior, to_cmp, false );
                    p = VDBMemBufferAccess ( & prior, & poff, & pavail );
                    c = VDBMemBufferAccess ( & cur, & coff, & cavail );
                    if ( p == NULL || c == NULL )
                        return true;

                    to_cmp = num_bits - cnt;
                    if ( to_cmp > pavail )
                        to_cmp = pavail;
                    if ( to_cmp > cavail )
                        to_cmp = cavail;

                    /* compare bits */
                    if ( bitcmp ( p, poff, c, coff, to_cmp ) != 0 )
                        break;

                    /* advance buffers */
                    VDBMemBufferAdvance ( & cur, to_cmp, false );
                }

                /* if the rows were identical */
                if ( cnt == num_bits )
                {
                    /* just bump the repeat count */
                    ++ rm -> cnt;

                    if ((self->cutoff_id != self->start_id) && (self->cutoff_id + 1 == *end_id))
                        self->cutoff_id = *end_id;
#if 1
                    /* drop any extra buffers used during write */
                    for ( mem = ( VDBMem* ) DLNodeNext ( & prior . mem -> n );
                          mem != NULL;
                          mem = ( VDBMem* ) DLNodeNext ( & prior . mem -> n ) )
                    {
                        DLListUnlink ( & self -> data, & mem -> n );
                        VDBMemRelease ( & mem -> n, self -> mgr );
                    }
#else
                    while ( DLListTail ( & self -> data ) != prior . mem -> n )
                    {
                        mem = ( VDBMem * ) DLListPopTail ( & self -> data );
                        VDBMemRelease ( & mem -> n, self -> mgr );
                    }
#endif                    

                    /* consider the row committed */
                    self -> row_len = 0;
                    self -> row_committed = true;
                    goto check_size;
                }
            }
        }
    }

    /* offset to next rowmap entry */
    rmoff = ( self -> num_rows + self -> rowmap_off ) * sizeof * rm;

    /* check for space */
    if ( mem == NULL ||( rmoff & ( sizeof mem -> page - 1 ) ) == 0 )
    {
        rc_t rc = VDBManagerMakeMem ( self -> mgr, & mem );
        if ( rc != 0 )
            return true;
        DLListPushTail ( & self -> rowmap, & mem -> n );
        rmoff = 0;
    }

    /* mem might not be initialized! */
    /* new rowmap entry */
    rm = ( void* ) & mem -> page [ rmoff & ( sizeof mem -> page - 1 ) ];
#if _DEBUGGING
    rm -> start_id = self -> end_id;
#endif
    rm -> len = self -> row_len;
    rm -> cnt = 1;

    self -> num_elems += self -> row_len;
    ++ self -> num_rows;
    self -> row_len = 0;
    self -> row_committed = true;
check_size:
    /* current size in bytes */
    cur_size = ( (size_t)self -> num_elems * elem_bits + 7 ) >> 3;

    /* if the buffer is large enough */
    if ( cur_size >= self -> trigger )
    {

        /* if size just crossed the trigger boundary and 
         * cutoff_id has not been advanced yet */
        if ( self -> cutoff_id == self -> start_id )
        {
            self -> cutoff_id = * end_id;
        }

        /* or perhaps the buffer is too large */
        else if ( ( cur_size + cur_size ) >= self -> trigger * 3 )
        {
            /* set to min of current end or our cutoff */
            if ( self -> cutoff_id < * end_id )
            {
                * end_id = self -> cutoff_id;
            }
        }
    }

    /* if the row range is too great ( 3G rows ) */
    else if ( ( self -> end_id - self -> start_id ) >= 0xC0000000 )
    {
        /* if row range has just crossed the boundary and 
         * cutoff_id has not been advanced yet */
        if ( self -> cutoff_id == self -> start_id )
        {
            self -> cutoff_id = * end_id;
        }

        /* set to min of current end or our cutoff */
        else if ( self -> cutoff_id < * end_id )
        {
            * end_id = self -> cutoff_id;
        }
    }

    return false;
}

/* CloseRow
 *  discards uncommitted data
 *  update state
 */
void CC WColumnCloseRow ( void *item, void *ignore )
{
    WColumn *self = item;

    if ( self -> row_len != 0 )
    {
        /* discard any extra buffers used */
        VDBMem *mem;
        bitsz_t elem_bits = VTypedescSizeof ( & self -> dad . desc );
        bitsz_t boff = self -> data_off + self -> num_elems * elem_bits;
        bitsz_t origin = boff + self -> row_len * elem_bits;

        /* set origin to the start of last buffer */
        origin = ( origin - 1 ) & ~ ( bitsz_t ) ( sizeof mem -> page * 8 - 1 );

        if ((boff == 0) && (origin == 0)) /* we got nothing */
        {
            /* delete all but it should only be one if any */
            while ((mem = (VDBMem*)DLListPopTail (&self->data)) != NULL)
                VDBMemRelease ( & mem -> n, self -> mgr );
        }
        else
        {
            /* not reached, not tested, not understood */
            for ( ; origin >= boff ; origin -= sizeof mem -> page * 8 )
            {
                mem = ( VDBMem* ) DLListPopTail ( & self -> data );
                assert ( mem != NULL );
                VDBMemRelease ( & mem -> n, self -> mgr );
            }
        }
    }

    if ( self -> row_committed )
        ++ self -> end_id;

    self -> row_len = 0;
    self -> row_written = false;
    self -> row_committed = false;
}

/* BufferPage
 *  captures page range
 *
 *  "end_id" [ IN, CONST ] - half-closed id of buffered range end
 *  column should capture this information for creating page
 *  id range either on demand, or pre-prepared.
 *
 *  returns true if there was a memory error.
 */
static
rc_t WColumnMakePage ( WColumn *self, int64_t end_id,
   bitsz_t elem_bits, uint64_t num_elems, uint64_t num_rows )
{
    rc_t rc;

    if ( self -> page != NULL )
    {
        TRACK_BLOB (VBlobRelease, self->page);
        ( void ) VBlobRelease ( self -> page );
        self -> page = NULL;
    }

    rc = VBlobNew ( & self -> page, self -> start_id, end_id - 1, self -> dad . scol -> name -> name . addr );
    TRACK_BLOB (VBlobNew, self->page);

    if ( rc == 0 )
    {
        VBlob *vblob = self -> page;
        rc = KDataBufferMake ( & vblob -> data, elem_bits, num_elems );
        if ( rc == 0 )
        {
            rc = PageMapNew ( & vblob -> pm, num_rows );
            if ( rc == 0 )
                return 0;
        }

        TRACK_BLOB (VBlobRelease, vblob);
        (void)VBlobRelease ( vblob );
        self -> page = NULL;
    }
    return rc;
}

#define ROWMAPS_PER_PAGE (sizeof(((const VDBMem *)0)->page) / sizeof(WColumnRowMap))
#define ROWMAP_BITS (sizeof(WColumnRowMap) * 8)
/* false return is good, true return is failure of some sort */
bool CC WColumnBufferPage ( void *item, void *const_end_id )
{
    rc_t rc;
    WColumn *self = item;
    int64_t id, end_id = * ( const int64_t* ) const_end_id;

    VDBMem *mem;
    bool seek_ok;
    const WColumnRowMap *rm = NULL;
    VDBMemBuffer rowmap, data;
    uint64_t num_elems, num_rows;
    uint64_t rowmap_off;
    bitsz_t cnt, num_bits, to_cpy, elem_bits;
    bool advance = true;

    /* initialize buffer onto rowmap */
    VDBMemBufferInitHead ( & rowmap, DLListHead ( & self -> rowmap ),
        self -> rowmap_off * ROWMAP_BITS, ROWMAP_BITS, self -> num_rows );

    /* walk row map to count rows and elements through end_id not the whole list */
    num_elems = num_rows = 0;
    for ( id = self -> start_id;
          id < end_id;
          id += rm -> cnt )
    {
        rm = VDBMemBufferAccess ( & rowmap, NULL, NULL );
        if ( rm == NULL )
            return true;
        if (advance)
            advance  = VDBMemBufferAdvance ( & rowmap, ROWMAP_BITS, false );
        else
        {
            /* assert (0); */
            return true;
        }
        num_elems += rm -> len;
        assert (rm->start_id == id);
        ++ num_rows;
    }

    /* allocate page */
    elem_bits = VTypedescSizeof ( & self -> dad . desc );
    rc = WColumnMakePage ( self, end_id, elem_bits, num_elems, num_rows );
    if ( rc != 0 )
        return true;

    /* initialize data buffer */
    VDBMemBufferInitHead ( & data, DLListHead ( & self -> data ),
          self -> data_off, elem_bits, self -> num_elems + self -> row_len );

    /* copy data */
    num_bits = KDataBufferBits ( & self -> page -> data );
    for ( seek_ok = true, cnt = 0; cnt < num_bits; cnt += to_cpy )
    {
        bitsz_t avail, soff;
        const void *src = VDBMemBufferAccess ( & data, & soff, & avail );
        if ( src == NULL )
            return true;

        to_cpy = num_bits - cnt;
        if ( to_cpy > avail )
            to_cpy = avail;

        bitcpy ( self -> page -> data . base, cnt, src, soff, to_cpy );

        seek_ok = VDBMemBufferAdvance ( & data, to_cpy, false );
    }

    /* if last seek failed, "to_cpy" represents the bits remaining */
    if ( seek_ok )
        to_cpy = 0;

    /* capture left edge of data */
    self -> data_off = data . marker + to_cpy;

    /* if keeping the last row, back off by that amount:
     * we could only "over shoot" by one in the loop above */
    if ( id > end_id )
    {
        self -> data_off -= elem_bits * rm -> len;
    }
    /* drop data pages */
#if 1
    for ( ; ; ) {
        bitsz_t to_drop = 0;

        if (self->data_off >= sizeof(mem->page) * 8)
            to_drop = sizeof(mem->page) * 8;
        else if (self->data_off == data.end)
            to_drop = data.end;
        if (to_drop) {
            mem = ( VDBMem* ) DLListPopHead ( & self -> data );
            VDBMemRelease ( & mem -> n, self -> mgr );
            self->data_off -= to_drop;
            data.end -= to_drop;
        }
        else break;
    }
#else
    for ( ; self -> data_off >= sizeof mem -> page * 8; self -> data_off -= sizeof mem -> page * 8 )
    {
        mem = ( VDBMem* ) DLListPopHead ( & self -> data );
        VDBMemRelease ( & mem -> n, self -> mgr );
    }
#endif

    /* re-initialize rowmap buffer */
    VDBMemBufferInitHead ( & rowmap, DLListHead ( & self -> rowmap ),
        self -> rowmap_off * ROWMAP_BITS, ROWMAP_BITS, self -> num_rows );

    /* copy row map */
    for ( ; num_rows > 1; -- num_rows )
    {
        rm = VDBMemBufferAccess ( & rowmap, NULL, NULL );
        assert ( rm != NULL );
        rc = PageMapAppendSomeRows ( self -> page -> pm, rm -> len, rm -> cnt );
        if ( rc != 0 )
            return true;
        VDBMemBufferAdvance ( & rowmap, ROWMAP_BITS, false );
    }
    for ( seek_ok = true; num_rows > 0; -- num_rows )
    {
        uint64_t row_count;
        
        rm = VDBMemBufferAccess ( & rowmap, NULL, NULL );
        assert ( rm != NULL );
        
        row_count = rm -> cnt - ( id - end_id );
        rc = PageMapAppendSomeRows ( self -> page -> pm, rm -> len, row_count);
        if ( rc != 0 )
            return true;
#if _DEBUGGING
        if (id - end_id)
        {
            ((WColumnRowMap *)rm)->start_id += row_count;
        }
#endif
        seek_ok = VDBMemBufferAdvance ( & rowmap, ROWMAP_BITS, false );
    }

    /* capture left edge of rowmap:  that is capture the number of bits to the new
     * offset to the first unwritten rowmap that might be a different page than the
     * current rowmap offset */
    rowmap_off = rowmap . marker / ROWMAP_BITS;

    /* if last advance failed, incorporate bits */
    if ( ! seek_ok )
        ++rowmap_off;

    /* if last rowmap entry needs to be preserved */
    if ( id > end_id )
    {
        /* back up to preserve row */
        --rowmap_off;
        /* adjust row count */
        ( ( WColumnRowMap* ) rm ) -> cnt = id - end_id;
    }

    /* drop rowmap pages and fix up other fields
     * this is less complex if we are clearing out all pages */
    if (rowmap_off - self->rowmap_off == self->num_rows)
    {
        /* delete all */
        while ((mem = (VDBMem*)DLListPopHead (&self->rowmap)) != NULL)
            VDBMemRelease (&mem->n, self->mgr);

        self->num_elems = 0;
        self->num_rows = 0;
        self->rowmap_off = 0;
    }
    else
    {
        uint64_t rowmaps_this_page;

        /* page must be exactly divisible by sizeof ( WColumnRowMap ) */
        assert ( VDBMemBitsInPage () % sizeof ( WColumnRowMap ) == 0 );

        /* Kind of ugly:  The first page might have less that the
         * maximum number of rowmaps due to its having an offset
         * to the first current rowmap.  When that page is dropped
         * the number of rows needs to be reduced by that partial
         * number of rowmaps.  But this seemed preferable to
         * a more awkward check for a first partial page followed by
         * a loop over full pages. */
        for ( rowmaps_this_page = ROWMAPS_PER_PAGE - self->rowmap_off;
              rowmap_off >= ROWMAPS_PER_PAGE;
              rowmaps_this_page = ROWMAPS_PER_PAGE )
        {
            mem = ( VDBMem* ) DLListPopHead ( & self -> rowmap );
            assert ( mem != NULL );
            VDBMemRelease ( & mem -> n, self -> mgr );
            rowmap_off -= ROWMAPS_PER_PAGE;
            assert ( self->num_rows > rowmaps_this_page );
            self->num_rows -= rowmaps_this_page;
            self->rowmap_off = 0;     /* old rowmap_off is no longer on the same page */
        }
        /* new number of rows needs to lose the rows on this page that are done;
         * that might be just the new rowmap_off or it might be new - old
         * this expression assumes new is page decremented in the loop and
         * the old is dropped in the loop or that the loop is not run
         * and the new and old are both from the same page
         */
        assert ( self->num_rows > (rowmap_off - self->rowmap_off));
        self->num_rows -= (rowmap_off - self->rowmap_off);
        self->rowmap_off = rowmap_off;

        /* reinitialize buffer onto rowmap */
        VDBMemBufferInitHead ( & rowmap, DLListHead ( & self -> rowmap ),
                               rowmap_off * ROWMAP_BITS, ROWMAP_BITS, self -> num_rows );

        /* walk row map to count actual elements */
        for ( num_elems = 0, num_rows = 0; num_rows != self -> num_rows; ++num_rows )
        {
            rm = VDBMemBufferAccess ( & rowmap, NULL, NULL );
            if ( rm == NULL )
                return true;
            VDBMemBufferAdvance ( & rowmap, ROWMAP_BITS, false );
            num_elems += rm->len;
        }
        self->num_elems = num_elems;
    }

    self->start_id = self->cutoff_id = end_id;
    return false;
}

/* ReadBlob
 *  reads an input blob
 *  called as a result of commit page which reads the validation production
 */
rc_t WColumnReadBlob ( WColumn *self, VBlob **vblob, int64_t id )
{
    if ( self -> page == NULL )
        return RC ( rcVDB, rcColumn, rcReading, rcBuffer, rcNotFound );
    if ( id < self -> page -> start_id || id > self -> page -> stop_id )
        return RC ( rcVDB, rcColumn, rcReading, rcRow, rcNotFound );

    * vblob = self -> page;
    (void)VBlobAddRef ( self -> page );
    TRACK_BLOB(VBlobAddRef, self->page);

    return 0;
}

/* DropPage
 *  drops any page buffers created
 */
void CC WColumnDropPage ( void *item, void *ignore )
{
    WColumn *self = item;
    if ( self -> page != NULL )
    {
        TRACK_BLOB(VBlobRelease,self->page);
        (void)VBlobRelease ( self -> page );
        self -> page = NULL;
    }
}
