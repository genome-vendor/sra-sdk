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

#ifndef _h_cursor_priv_
#define _h_cursor_priv_

#ifndef _h_vdb_extern_
#include <vdb/extern.h>
#endif

#ifndef _h_klib_vector_
#include <klib/vector.h>
#endif

#ifndef _h_klib_container_
#include <klib/container.h>
#endif

#ifndef _h_klib_refcount_
#include <klib/refcount.h>
#endif

#ifndef KONST
#define KONST
#endif

#ifndef SKONST
#define SKONST const
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include "blob-priv.h"


#define VCURSOR_WRITE_MODES_SUPPORTED 0

/* temporary - until the full kproc for Windows is operational */
#if WINDOWS
#define VCURSOR_FLUSH_THREAD 0
#else
#define VCURSOR_FLUSH_THREAD 1
#endif

/*--------------------------------------------------------------------------
 * forwards
 */
struct KLock;
struct KCondition;
struct KThread;
struct KDlset;
struct VTable;
struct VCtxId;
struct VSchema;
struct SColumn;
struct VColumn;
struct VPhysical;


/*--------------------------------------------------------------------------
 * VCursorCache
 *  a dual-level vector
 */
typedef struct VCursorCache VCursorCache;
struct VCursorCache
{
    Vector cache;
};


/* Init
 */
#define VCursorCacheInit( self, start, block ) \
    VectorInit ( & ( self ) -> cache, start, block )

/* Whack
 */
void VCursorCacheWhack ( VCursorCache *self,
    void ( CC * whack ) ( void *item, void *data ), void *data );

/* Get
 *  retrieve object by cid
 */
void *VCursorCacheGet ( const VCursorCache *self, struct VCtxId const *cid );

/* Set
 *  set object by cid
 */
rc_t VCursorCacheSet ( VCursorCache *self,
    struct VCtxId const *cid, const void *item );

/* Swap
 *  swap object by cid
 */
rc_t VCursorCacheSwap ( VCursorCache *self,
    struct VCtxId const *cid, const void *item, void **prior );


/*--------------------------------------------------------------------------
 * VCursor
 */
enum
{
    vcConstruct,
    vcFailed,
    vcReady,
    vcRowOpen,
    vcRowCommitted,
    vcPageCommit
};

enum
{
    vfReady,
    vfBusy,
    vfFgErr,
    vfBgErr,
    vfExit
};

#define LAST_BLOB_CACHE_SIZE 32

struct VCursor
{
    /* row id */
    int64_t row_id;

    /* half-closed page range */
    int64_t start_id, end_id;

    /* starting id for flush */
    volatile int64_t flush_id;

    /* attached reference to table */
    struct VTable KONST *tbl;

    /* cursor-specific schema and table */
    struct VSchema SKONST *schema;
    struct STable SKONST *stbl;

    /* background flush thread objects */
    struct KThread *flush_thread;
    struct KLock *flush_lock;
    struct KCondition *flush_cond;

    /* user data */
    void *user;
    void ( CC * user_whack ) ( void *data );

    /* external named cursor parameters */    
    BSTree named_params;

    /* read-only blob cache */
    BSTree blob_cache;
    DLList blob_lru;
    size_t blob_cache_capacity;
    size_t blob_cache_contents;
    /* last blob cache */
    VBlob *last_blob_cache[LAST_BLOB_CACHE_SIZE]; /** last blob to be cached per given col_idx, limiting col_idx <=32 **/


    /* external row of VColumn* by ord ( owned ) */
    Vector row;

    /* column objects by cid ( not-owned ) */
    VCursorCache col;

    /* physical columns by cid ( owned ) */
    VCursorCache phys;

    /* productions by cid ( not-owned ) */
    VCursorCache prod;

    /* intermediate productions ( owned ) */
    Vector owned;

    /* trigger productions ( not-owned ) */
    Vector trig;

    KRefcount refcount;

    volatile uint32_t flush_cnt;

    /* foreground state */
    uint8_t state;

    /* flush_state */
    volatile uint8_t flush_state;

    bool read_only;

    /* support for sradb-v1 API */
    bool permit_add_column;
    bool permit_post_open_add;
};


/* Make
 */
rc_t VCursorMake ( struct VCursor **cursp, struct VTable const *tbl );

/* Whack
 * Destroy
 */
rc_t VCursorWhack ( struct VCursor *self );
rc_t VCursorDestroy ( struct VCursor *self );

/* SupplementSchema
 *  scan table for physical column names
 *  create transparent yet incomplete (untyped) columns for unknown names
 *  create incomplete (untyped) physical columns for forwarded names
 *  repeat process on static columns, except create complete (fully typed) objects
 */
rc_t VCursorSupplementSchema ( struct VCursor const *self );

/* MakeColumn
 */
rc_t VCursorMakeColumn ( struct VCursor *self,
    struct VColumn **col, struct SColumn const *scol );

/* SetRowIdRead - PRIVATE
 *  seek to given row id
 *
 *  "row_id" [ IN ] - row id to select
 */
rc_t VCursorSetRowIdRead ( struct VCursor *self, int64_t row_id );

/* Open
 */
rc_t VCursorOpenRead ( struct VCursor *self, struct KDlset const *libs );

/* ListReadableColumns
 *  performs an insert of '*' to cursor
 *  attempts to resolve all read rules
 *  records all SColumns that successfully resolved
 *  populates BTree with VColumnRef objects
 */
rc_t VCursorListReadableColumns ( struct VCursor *self, BSTree *columns );

/* ListWritableColumns
 *  walks list of SPhysicals and trigger SProductions
 *  attempts to resolve all write rules
 *  records any SColumn that can be reached
 *  populates BTree with VColumnRef objects
 */
rc_t VCursorListWritableColumns ( struct VCursor *self, BSTree *columns );

/* PostOpenAdd
 *  handle opening of a column after the cursor is opened
 */
rc_t VCursorPostOpenAdd ( struct VCursor *self, struct VColumn *col );
rc_t VCursorPostOpenAddRead ( struct VCursor *self, struct VColumn *col );

/* OpenRowRead
 * CloseRowRead
 */
rc_t VCursorOpenRowRead ( struct VCursor *self );
rc_t VCursorCloseRowRead ( struct VCursor *self );


#ifdef __cplusplus
}
#endif

#endif /* _h_cursor_priv_ */
