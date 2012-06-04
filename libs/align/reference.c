/*===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was readten as part of
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
#include <klib/vector.h>
#include <klib/out.h>
#include <insdc/insdc.h>
#include <vdb/manager.h>
#include <vdb/database.h>
#include <vdb/cursor.h>
#include <vdb/table.h>
#include <vdb/vdb-priv.h>
#include <align/iterator.h>
#include <align/reference.h>
#include <align/refseq-mgr.h>
#include <os-native.h>
#include <sysalloc.h>

#include "reader-cmn.h"
#include "reference-cmn.h"
#include "debug.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <assert.h>

enum EReferenceList_ColNames {
    ereflst_cn_READ_dna,
    ereflst_cn_READ_4na,
    ereflst_cn_SEQ_LEN,
    ereflst_cn_PRIMARY_ALIGNMENT_IDS,
    ereflst_cn_SECONDARY_ALIGNMENT_IDS,
    ereflst_cn_EVIDENCE_ALIGNMENT_IDS,
    ereflst_cn_OVERLAP_REF_POS,
    ereflst_cn_OVERLAP_REF_LEN
};

const TableReaderColumn ReferenceList_cols[] =
{
    {0, "(INSDC:dna:text)READ", {NULL}, 0, 0},
    {0, "(INSDC:4na:bin)READ", {NULL}, 0, ercol_Skip},
    {0, "SEQ_LEN", {NULL}, 0, 0},
    {0, "PRIMARY_ALIGNMENT_IDS", {NULL}, 0, ercol_Skip},
    {0, "SECONDARY_ALIGNMENT_IDS", {NULL}, 0, ercol_Skip},
    {0, "EVIDENCE_ALIGNMENT_IDS", {NULL}, 0, ercol_Skip},
    {0, "OVERLAP_REF_POS", {NULL}, 0, ercol_Optional},
    {0, "OVERLAP_REF_LEN", {NULL}, 0, ercol_Optional},
    {0, NULL, {NULL}, 0, 0}
};

enum EPlacementIterator_ColNames {
    eplacementiter_cn_REF_POS,
    eplacementiter_cn_REF_LEN,
    eplacementiter_cn_MAPQ
};

const TableReaderColumn PlacementIterator_cols[] =
{
    {0, "REF_POS", {NULL}, 0, 0},
    {0, "REF_LEN", {NULL}, 0, 0},
    {0, "MAPQ", {NULL}, 0, 0},
    {0, NULL, {NULL}, 0, 0}
};

struct ReferenceList {
    KRefcount refcount;
    const RefSeqMgr* refseqmgr;
    const VCursor* cursor;
    BSTree name_tree;
    BSTree seqid_tree;
    uint32_t options;
    size_t cache;
    uint32_t max_seq_len;
    uint32_t nodes_qty;
    uint32_t nodes_max_qty;
    const TableReader* reader;
    TableReaderColumn reader_cols[sizeof(ReferenceList_cols)/sizeof(ReferenceList_cols[0])];
    const TableReader* iter;
    TableReaderColumn iter_cols[sizeof(PlacementIterator_cols)/sizeof(PlacementIterator_cols[0])];
    /* last are children using realloc!! */
    ReferenceObj* nodes[2];
};

struct ReferenceObj {
    /* we use this in 2 lists so we need to adjust results of Find, etc calls for name tree!! */
    BSTNode by_seqid; /* primary key */
    BSTNode by_name; /* struct addr will be by_name[-1]; */
    ReferenceList* mgr;
    uint32_t id;
    uint32_t bin;
    char* name;
    char* seqid;
    bool circular;
    bool read_present;
    int64_t start_rowid;
    int64_t end_rowid;
    INSDC_coord_len seq_len;
};

static
int CC ReferenceObj_CmpSeqId(const void *item, const BSTNode *n)
{
    return strcasecmp((const char*)item, ((const ReferenceObj*)n)->seqid);
}

static
int CC ReferenceObj_SortSeqId(const BSTNode *item, const BSTNode *n)
{
    return ReferenceObj_CmpSeqId(((const ReferenceObj*)item)->seqid, n);
}

static
int CC ReferenceObj_CmpName(const void *item, const BSTNode *n)
{
    return strcasecmp((const char*)item, ((const ReferenceObj*)&n[-1])->name);
}

static
int CC ReferenceObj_SortName(const BSTNode *item, const BSTNode *n)
{
    return ReferenceObj_CmpName(((const ReferenceObj*)&item[-1])->name, n);
}

static
rc_t ReferenceObj_Alloc(ReferenceObj** self, const char* seqid, size_t seqid_sz, const char* name, size_t name_sz)
{
    rc_t rc = 0;
    ReferenceObj* obj;

    if( self == NULL || seqid == NULL || seqid_sz == 0 || name == NULL || name_sz == 0 ) {
        rc = RC(rcAlign, rcIndex, rcConstructing, rcParam, rcNull);
    } else if( (obj = calloc(1,  sizeof(*obj) + seqid_sz + 1 + name_sz + 1)) == NULL ) {
        rc = RC(rcAlign, rcIndex, rcConstructing, rcMemory, rcExhausted);
    } else {
        obj->seqid = (char*)&obj[1];
        obj->name = obj->seqid;
        obj->name += seqid_sz + 1;
        memcpy(obj->seqid, seqid, seqid_sz);
        obj->seqid[seqid_sz] = '\0';
        memcpy(obj->name, name, name_sz);
        obj->name[name_sz] = '\0';
        *self = obj;
    }
    return rc;
}

LIB_EXPORT rc_t CC ReferenceList_MakeCursor(const ReferenceList** cself, const VCursor* cursor, uint32_t options,
                                            const char* filt_name, const uint32_t numbins)
{
    rc_t rc = 0;
    ReferenceList* self = NULL;

    if( cself == NULL || cursor == NULL ) {
        rc = RC(rcAlign, rcType, rcConstructing, rcParam, rcNull);
    } else if( (self = calloc(1, sizeof(*self))) == NULL ) {
        rc = RC(rcAlign, rcType, rcConstructing, rcMemory, rcExhausted);
    } else {
        const TableReader* tmp = NULL;
        uint32_t bin_size;
        int bin_num = -1;
        TableReaderColumn h[] =
        {
            {0, "NAME", {NULL}, 0, 0},       /*0*/
            {0, "SEQ_ID", {NULL}, 0, 0},     /*1*/
            {0, "SEQ_LEN", {NULL}, 0, 0},    /*2*/
            {0, "CIRCULAR", {NULL}, 0, 0},   /*3*/
            {0, "MAX_SEQ_LEN", {NULL}, 0, 0},/*4*/
            {0, "SEQ_START", {NULL}, 0, 0},  /*5*/
            {0, "CMP_READ", {NULL}, 0, 0},   /*6*/
            {0, NULL, {NULL}, 0, 0}
        };
        KRefcountInit(&self->refcount, 1, "ReferenceList", "Make", "align" );
        BSTreeInit(&self->name_tree);
        BSTreeInit(&self->seqid_tree);
        self->options = options;
        self->nodes_max_qty = sizeof(self->nodes) / sizeof(self->nodes[0]);

        if( (rc = VCursorAddRef(self->cursor = cursor)) == 0 &&
            (rc = TableReader_MakeCursor(&tmp, cursor, h)) == 0 ) {
                int64_t start, tbl_start,tbl_stop;
                uint64_t count;
                const KIndex* iname = NULL;
                bool only_one = false;

                /* index is optional */
                rc_t rctmp = TableReader_OpenIndex(tmp, "i_name", &iname);
                ALIGN_DBGERRP("index '%s' was not found", rctmp, "i_name");

                rc = TableReader_IdRange(tmp, &start, &count);
                assert(rc == 0);
                tbl_start = start;
		tbl_stop  = start + count -1;
                if(numbins > 0) {
                    bin_size = (count + numbins -1) / numbins;
                } else {
                    bin_size = 0;
                }
                if(iname && filt_name) {
                    if(bin_size == 0) {
                        only_one = true;
                    }
                    if(strncmp(filt_name, "START_ROW:", 10)) {
                        if(bin_size == 0 || strncmp(filt_name, "BIN_NUM:", 8)) {
                            rc = KIndexFindText(iname, filt_name, &start, &count, NULL, NULL);
                            if(rc == 0 && bin_size > 0) {
                                /** change start to the beginning of the bin **/
                                bin_num = (start - tbl_start) / bin_size;
                            }
                        } else {
                            bin_num = atoi(filt_name + 8);
                        }
                        if(bin_size > 0) {
                            start = tbl_start + bin_size * bin_num;
                            rc = TableReader_ReadRow(tmp, start);
                            if(rc == 0) {
                                int64_t r_start;
                                char name[4096];
                                if(h[0].len < sizeof(name)){
                                    memcpy(name, h[0].base.str, h[0].len);
                                    name[h[0].len] = '\0';
                                    rc = KIndexFindText(iname, name, &r_start, &count, NULL, NULL);
                                    if(rc == 0 && start > r_start) { /*** move start to the beginning of the fully contained sequence **/
                                        start = r_start + count;
                                    }
                                } else {
                                    rc = RC(rcAlign, rcType, rcConstructing, rcName, rcTooLong);
                                }
                            }
                        }
                    } else {
                        int64_t req_start = atoi(filt_name + 10);
                        if(req_start >= start && req_start < start + count){
                            int64_t delta = req_start - start;
                            start = req_start;
                            count -= delta;
                        } else {
                            rc = RC(rcAlign, rcType, rcConstructing, rcId, rcOutofrange);
                        }
                    }
                }
                if( rc == 0 ) {
                    ReferenceObj* node = NULL;
                    uint32_t last_name_len = 0;
                    bool read_determination_done = false;

                    while(rc == 0  && start <= tbl_stop) {
                        if(bin_num < 0 && count == 0) {
                            /*** normal loop without binning ***/
                            break;
                        }
                        if( (rc = TableReader_ReadRow(tmp, start)) == 0 ) {
                            if( node == NULL || last_name_len != h[0].len || strncmp(h[0].base.str, node->name, h[0].len) != 0 ) {
                                uint32_t cur_bin = (bin_size > 0) ? (start-tbl_start) / bin_size : 0;
                                if( only_one && self->nodes_qty == 1 ) {
                                    break;
                                }
                                if( bin_num >= 0 && cur_bin != bin_num ) {
                                    break;
                                }
                                if( node == NULL && h[4].len > 0 ) {
                                    self->max_seq_len = h[4].base.u32[0];
                                }
                                if( self->nodes_qty == self->nodes_max_qty ) {
                                    ReferenceList* tmp = realloc(self, sizeof(*self) + sizeof(node) * self->nodes_max_qty);
                                    if( tmp == NULL ) {
                                        rc = RC(rcAlign, rcType, rcConstructing, rcMemory, rcExhausted);
                                    } else {
                                        self = tmp;
                                        self->nodes_max_qty += sizeof(self->nodes) / sizeof(self->nodes[0]);
                                    }
                                }
                                if( rc == 0 &&
                                    (rc = ReferenceObj_Alloc(&node, h[1].base.str, h[1].len, h[0].base.str, h[0].len)) == 0 ) {
                                        node->id = self->nodes_qty;
                                        self->nodes[self->nodes_qty++] = node;
                                        last_name_len = h[0].len;
                                        node->circular = h[3].len ? h[3].base.buul[0] : false;
                                        node->start_rowid = start;
                                        node->seq_len = 0;
                                        node->bin = cur_bin;
                                        read_determination_done = false;
                                        if( (rc = BSTreeInsertUnique(&self->seqid_tree, &node->by_seqid, NULL, ReferenceObj_SortSeqId)) == 0 ) {
                                            rc = BSTreeInsertUnique(&self->name_tree, &node->by_name, NULL, ReferenceObj_SortName);
                                        }
                                }
                            }
                            if( rc == 0 ) {
                                if( h[6].len > 0 ) {/** CMP_READ > 0 -- truly local ***/
                                    node->read_present = true; 
                                    read_determination_done = true;
                                } else if( h[5].base.coord1[0] != 0 ) { /*** truly remote ***/
                                    node->read_present = false;
                                    read_determination_done = true;
                                } /*** else still not sure **/
                                if( read_determination_done && iname != NULL ) {
                                    /* scroll to last row for this reference projecting the seq_len */
                                    int64_t r_start;
                                    uint64_t r_count;
                                    if( KIndexFindText(iname, node->name, &r_start, &r_count, NULL, NULL) == 0 ) {
                                        assert(node->start_rowid == r_start);
                                        /* not last ref row */
                                        if( start != r_start + r_count - 1 ) {
                                            /* we need to pickup last row SEQ_LEN for this reference
                                            so we step back 2 rows in table from this ref end row
                                            and also skip rows already scanned for read presence */
                                            r_count -= (start - r_start) + 2;
                                            node->seq_len += h[2].base.coord_len[0] * r_count;
                                            start += r_count;
                                            count -= r_count;
                                        }
                                    }
                                }
                                node->seq_len += h[2].base.coord_len[0];
                                node->end_rowid = start;
                            }
                        } else if( GetRCState(rc) == rcNotFound && GetRCObject(rc) == rcRow ) {
                            rc = 0;
                        }
                        start++;
                        count--;
                    }
                    for(start = 0; rc == 0 && start < self->nodes_qty; start++) {
                        self->nodes[start]->mgr = self;
                    }
                    if( rc == 0 && self->max_seq_len == 0 ) {
                        rc = RC(rcAlign, rcType, rcConstructing, rcData, rcCorrupt);
                    }
                }
                KIndexRelease(iname);
        }
        TableReader_Whack(tmp);
    }
    if( rc == 0 ) {
        *cself = self;
        ALIGN_DBG("created 0x%p with cursor 0x%p", self, cursor);
    } else {
        *cself = NULL;
        ReferenceList_Release(self);
        ALIGN_DBGERRP("failed for cursor 0x%p", rc, cursor);
    }
    return rc;
}

LIB_EXPORT rc_t CC ReferenceList_MakeTable(const ReferenceList** cself, const VTable* table, uint32_t options,
                                           size_t cache, const char* filt_name, const uint32_t numbins)
{
    rc_t rc = 0;
    const VCursor* curs;

    if( table == NULL ) {
        rc = RC(rcAlign, rcType, rcConstructing, rcParam, rcNull);
    } else if( (rc = VTableCreateCachedCursorRead(table, &curs, cache)) == 0 &&
               (rc = VCursorPermitPostOpenAdd(curs)) == 0 ) {
        if( (rc = ReferenceList_MakeCursor(cself, curs, options, filt_name, numbins)) == 0 ) {
            ((ReferenceList*)(*cself))->cache = cache;
        }
        VCursorRelease(curs);
    }
    ALIGN_DBGERRP("failed for table 0x%p", rc, table);
    return rc;
}

LIB_EXPORT rc_t CC ReferenceList_MakeDatabase(const ReferenceList** cself, const VDatabase* db, uint32_t options,
                                              size_t cache, const char* name, const uint32_t numbins)
{
    rc_t rc = 0;
    const VTable* tbl = NULL;
    const char* nm = "REFERENCE";
    /*const char* nm = (options & ereferencelist_useEvidence) ? "EVIDENCE_INTERVAL" : "REFERENCE";*/

    if( db == NULL ) {
        rc = RC(rcAlign, rcType, rcConstructing, rcParam, rcNull);
    } else if( (rc = VDatabaseOpenTableRead(db, &tbl, nm)) == 0 ) {
        rc = ReferenceList_MakeTable(cself, tbl, options, cache, name, numbins);
        VTableRelease(tbl);
    }
    ALIGN_DBGERRP("failed for database 0x%p", rc, db);
    return rc;
}

LIB_EXPORT rc_t CC ReferenceList_MakePath(const ReferenceList** cself, const VDBManager* vmgr, const char* dbpath,
                                          uint32_t options, size_t cache, const char* name, const uint32_t numbins)
{
    rc_t rc = 0;
    const VDatabase* db = NULL;

    if( vmgr == NULL || dbpath == NULL ) {
        rc = RC(rcAlign, rcType, rcConstructing, rcParam, rcNull);
    } else if( (rc = VDBManagerOpenDBRead(vmgr, &db, NULL, "%s", dbpath)) == 0 ) {
        rc = ReferenceList_MakeDatabase(cself, db, options, cache, name, numbins);
        VDatabaseRelease(db);
    }
    ALIGN_DBGERRP("failed for database %s", rc, dbpath);
    return rc;
}

static
rc_t ReferenceList_OpenCursor(ReferenceList* self)
{
    rc_t rc = 0;

    assert(self != NULL);

    memcpy(self->reader_cols, ReferenceList_cols, sizeof(ReferenceList_cols));
    if( self->options & ereferencelist_4na ) {
        self->reader_cols[ereflst_cn_READ_dna].flags |= ercol_Skip;
        self->reader_cols[ereflst_cn_READ_4na].flags &= ~ercol_Skip;
    }
    if( self->options & ereferencelist_usePrimaryIds ) {
        self->reader_cols[ereflst_cn_PRIMARY_ALIGNMENT_IDS].flags &= ~ercol_Skip;
    }
    if( self->options & ereferencelist_useSecondaryIds ) {
        self->reader_cols[ereflst_cn_SECONDARY_ALIGNMENT_IDS].flags &= ~ercol_Skip;
    }
    if( self->options & ereferencelist_useEvidenceIds ) {
        self->reader_cols[ereflst_cn_EVIDENCE_ALIGNMENT_IDS].flags &= ~ercol_Skip;
    }
    if( !(self->options &
        (ereferencelist_usePrimaryIds | ereferencelist_useSecondaryIds | ereferencelist_useEvidenceIds )) ) {
        self->reader_cols[ereflst_cn_OVERLAP_REF_POS].flags |= ercol_Skip;
        self->reader_cols[ereflst_cn_OVERLAP_REF_LEN].flags |= ercol_Skip;
    }
    rc = TableReader_MakeCursor(&self->reader, self->cursor, self->reader_cols);
    ALIGN_DBGERR(rc);
    return rc;
}

static
rc_t ReferenceList_OpenCursor2(ReferenceList* self, align_id_src ids)
{
    rc_t rc = 0;
    const VDatabase* db = NULL;
    const VTable* vtbl = NULL;

    assert(self != NULL);

    if( ids != primary_align_ids && ids != secondary_align_ids && ids != evidence_align_ids ) {
        rc = RC(rcAlign, rcType, rcAccessing, rcParam, rcOutofrange);
    } else if( ids == primary_align_ids && !(self->options & ereferencelist_usePrimaryIds) ) {
        rc = RC(rcAlign, rcType, rcAccessing, rcParam, rcInvalid);
    } else if( ids == secondary_align_ids && !(self->options & ereferencelist_useSecondaryIds) ) {
        rc = RC(rcAlign, rcType, rcAccessing, rcParam, rcInvalid);
    } else if( ids == evidence_align_ids && !(self->options & ereferencelist_useEvidenceIds) ) {
        rc = RC(rcAlign, rcType, rcAccessing, rcParam, rcInvalid);
    } else if( (rc = VCursorOpenParentRead(self->cursor, &vtbl)) == 0 &&
               (rc = VTableOpenParentRead(vtbl, &db)) == 0 &&
               (rc = VTableRelease(vtbl)) == 0 &&
               (rc = VDatabaseOpenTableRead(db, &vtbl, ids == primary_align_ids ? "PRIMARY_ALIGNMENT" :
               (ids == secondary_align_ids ? "SECONDARY_ALIGNMENT" : "EVIDENCE_ALIGNMENT"))) == 0 ) {
        memcpy(self->iter_cols, PlacementIterator_cols, sizeof(PlacementIterator_cols));
        rc = TableReader_Make(&self->iter, vtbl, self->iter_cols, self->cache);
    }
    VDatabaseRelease(db);
    VTableRelease(vtbl);
    ALIGN_DBGERR(rc);
    return rc;
}

static
rc_t ReferenceList_RefSeqMgr(const ReferenceList* cself, const RefSeqMgr** rmgr)
{
    rc_t rc = 0;

    assert(rmgr != NULL);

    if( cself->refseqmgr == NULL ) {
        const VTable* vtbl = NULL;
        const VDBManager* vmgr;
        if( (rc = VCursorOpenParentRead(cself->cursor, &vtbl)) == 0 &&
            (rc = VTableOpenManagerRead(vtbl, &vmgr)) == 0 ) {
            rc = RefSeqMgr_Make(&((ReferenceList*)cself)->refseqmgr, vmgr,
                                (cself->options & ereferencelist_4na) ? errefseq_4NA : 0, cself->cache, 2);
            VDBManagerRelease(vmgr);
        }
        VTableRelease(vtbl);
    }
    *rmgr = cself->refseqmgr;
    ALIGN_DBGERR(rc);
    return rc;
}

LIB_EXPORT rc_t CC ReferenceList_AddRef(const ReferenceList *cself)
{
    rc_t rc = 0;
    if( cself != NULL ) {
        if( KRefcountAdd(&cself->refcount, "ReferenceList") != krefOkay ) {
            rc = RC(rcAlign, rcType, rcAttaching, rcError, rcUnexpected);
        }
    }
    return rc;
}

LIB_EXPORT void CC ReferenceList_Release(const ReferenceList* cself)
{
    if( cself != NULL ) {
        if( KRefcountDrop(&cself->refcount, "ReferenceList") == krefWhack ) {
            ReferenceList* self = (ReferenceList*)cself;
            TableReader_Whack(self->reader);
            TableReader_Whack(cself->iter);
            RefSeqMgr_Release(self->refseqmgr);
            while( self->nodes_qty-- > 0 ) {
                free(self->nodes[self->nodes_qty]);
            }
            VCursorRelease(cself->cursor);
            KRefcountWhack(&self->refcount, "ReferenceList");
            free(self);
        }
    }
}

LIB_EXPORT rc_t CC ReferenceList_Count(const ReferenceList* cself, uint32_t* count)
{
    rc_t rc = 0;
    if( cself == NULL || count == NULL ) {
        rc = RC(rcAlign, rcType, rcAccessing, rcParam, rcNull);
    } else {
        *count = cself->nodes_qty;
    }
    ALIGN_DBGERR(rc);
    return rc;
}

LIB_EXPORT rc_t CC ReferenceList_Find(const ReferenceList* cself, const ReferenceObj** obj, const char* key, size_t key_sz)
{
    rc_t rc = 0;
    char buf[4096], *b = buf;

    if( cself == NULL || obj == NULL || key == NULL ) {
        rc = RC(rcAlign, rcType, rcSearching, rcParam, rcNull);
    } else if( key_sz >= sizeof(buf) && (b = malloc(key_sz + 1)) == NULL ) {
        rc = RC(rcAlign, rcType, rcSearching, rcMemory, rcExhausted);
    } else {
        memcpy(b, key, key_sz);
        b[key_sz] = '\0';
        *obj = (ReferenceObj*)BSTreeFind(&cself->seqid_tree, b, ReferenceObj_CmpSeqId);
        if( *obj == NULL ) {
            const BSTNode* n = BSTreeFind(&cself->name_tree, b, ReferenceObj_CmpName);
            if( n != NULL ) {
                *obj = (ReferenceObj*)&n[-1];
            }
        }
        if( *obj == NULL ) {
            rc = RC(rcAlign, rcType, rcSearching, rcItem, rcNotFound);
        } else if( (rc = ReferenceList_AddRef(cself)) != 0 ) {
            *obj = NULL;
        }
        if( b != buf ) {
            free(b);
        }
    }
    ALIGN_DBGERR(rc);
    return rc;
}

LIB_EXPORT rc_t CC ReferenceList_Get(const ReferenceList* cself, const ReferenceObj** obj, uint32_t idx)
{
    rc_t rc = 0;
    if( cself == NULL || obj == NULL || idx >= cself->nodes_qty ) {
        rc = RC(rcAlign, rcType, rcRetrieving, rcParam, rcInvalid);
    } else {
        if( (rc = ReferenceList_AddRef(cself)) == 0 ) {
            *obj = cself->nodes[idx];
        } else {
            *obj = NULL;
        }
    }
    ALIGN_DBGERR(rc);
    return rc;
}

LIB_EXPORT rc_t ReferenceObj_AddRef(const ReferenceObj *cself)
{
    if( cself == NULL ) {
        return RC(rcAlign, rcType, rcAccessing, rcParam, rcInvalid);
    } else {
        return ReferenceList_AddRef( cself->mgr );
    }
}

LIB_EXPORT void ReferenceObj_Release(const ReferenceObj *cself)
{
    ReferenceList_Release(cself ? cself->mgr : NULL);
}

LIB_EXPORT rc_t CC ReferenceObj_Idx(const ReferenceObj* cself, uint32_t* idx)
{
    rc_t rc = 0;
    if( cself == NULL || idx == NULL ) {
        rc = RC(rcAlign, rcType, rcAccessing, rcParam, rcInvalid);
    } else {
        *idx = cself->id;
    }
    ALIGN_DBGERR(rc);
    return rc;
}

LIB_EXPORT rc_t CC ReferenceObj_IdRange(const ReferenceObj* cself, int64_t* start, int64_t* stop)
{
    rc_t rc = 0;
    if( cself == NULL || (start == NULL && stop == NULL) ) {
        rc = RC(rcAlign, rcType, rcAccessing, rcParam, rcInvalid);
    } else {
        if( start != NULL ) {
            *start = cself->start_rowid;
        }
        if( stop != NULL ) {
            *stop = cself->end_rowid;
        }
    }
    ALIGN_DBGERR(rc);
    return rc;
}
LIB_EXPORT rc_t CC ReferenceObj_Bin(const ReferenceObj* cself, uint32_t* bin)
{
    rc_t rc = 0;
    if( cself == NULL || bin == NULL ) {
        rc = RC(rcAlign, rcType, rcAccessing, rcParam, rcInvalid);
    } else {
        *bin = cself->bin;
    }
    ALIGN_DBGERR(rc);
    return rc;

}
LIB_EXPORT rc_t CC ReferenceObj_SeqId(const ReferenceObj* cself, const char** seqid)
{
    rc_t rc = 0;
    if( cself == NULL || seqid == NULL ) {
        rc = RC(rcAlign, rcType, rcAccessing, rcParam, rcInvalid);
    } else {
        *seqid = cself->seqid;
    }
    ALIGN_DBGERR(rc);
    return rc;
}

LIB_EXPORT rc_t CC ReferenceObj_Name(const ReferenceObj* cself, const char** name)
{
    rc_t rc = 0;
    if( cself == NULL || name == NULL ) {
        rc = RC(rcAlign, rcType, rcAccessing, rcParam, rcInvalid);
    } else {
        *name = cself->name;
    }
    ALIGN_DBGERR(rc);
    return rc;
}

LIB_EXPORT rc_t CC ReferenceObj_SeqLength(const ReferenceObj* cself, INSDC_coord_len* len)
{
    rc_t rc = 0;
    if( cself == NULL || len == NULL ) {
        rc = RC(rcAlign, rcType, rcAccessing, rcParam, rcInvalid);
    } else {
        *len = cself->seq_len;
    }
    ALIGN_DBGERR(rc);
    return rc;
}

LIB_EXPORT rc_t CC ReferenceObj_Circular(const ReferenceObj* cself, bool* circular)
{
    rc_t rc = 0;
    if( cself == NULL || circular == NULL ) {
        rc = RC(rcAlign, rcType, rcAccessing, rcParam, rcInvalid);
    } else {
        *circular = cself->circular;
    }
    ALIGN_DBGERR(rc);
    return rc;
}

LIB_EXPORT rc_t CC ReferenceObj_External(const ReferenceObj* cself, bool* external, char** path)
{
    rc_t rc = 0;

    if( cself == NULL || external == NULL ) {
        rc = RC(rcAlign, rcType, rcAccessing, rcParam, rcInvalid);
    } else {
        const RefSeqMgr* rmgr;
        *external = !cself->read_present;
        if( path != NULL && !cself->read_present && (rc = ReferenceList_RefSeqMgr(cself->mgr, &rmgr)) == 0 ) {
            *path = NULL;
            rc = RefSeqMgr_Exists(rmgr, cself->seqid, strlen(cself->seqid), path);
            if( GetRCObject(rc) == rcTable && GetRCState(rc) == rcNotFound ) {
                rc = 0;
            }
        }
    }
    ALIGN_DBGERR(rc);
    return rc;
}

LIB_EXPORT rc_t CC ReferenceObj_Read(const ReferenceObj* cself, INSDC_coord_zero offset, INSDC_coord_len len,
                                     uint8_t* buffer, INSDC_coord_len* written)
{
    rc_t rc = 0;

    if( cself == NULL || buffer == NULL || written == NULL ) {
        rc = RC(rcAlign, rcType, rcAccessing, rcParam, rcInvalid);
    } else if( (rc = ReferenceSeq_ReOffset(cself->circular, cself->seq_len, &offset)) == 0 ) {
        if( cself->mgr->reader != NULL || (rc = ReferenceList_OpenCursor(cself->mgr)) == 0 ) {

            int cid = (cself->mgr->options & ereferencelist_4na) ? ereflst_cn_READ_4na : ereflst_cn_READ_dna;
            INSDC_coord_len q = 0;
            *written = 0;
            do {
                int64_t rowid = cself->start_rowid + offset / cself->mgr->max_seq_len;
                INSDC_coord_zero s = offset % cself->mgr->max_seq_len;
                if( (rc = TableReader_ReadRow(cself->mgr->reader, rowid)) == 0 ) {
                    q = cself->mgr->reader_cols[ereflst_cn_SEQ_LEN].base.coord_len[0] - s;
                    if( q > len ) {
                        q = len;
                    }
                    memcpy(&buffer[*written], &cself->mgr->reader_cols[cid].base.str[s], q);
                    *written += q;
                    offset += q;
                    len -= q;
                }
                /* SEQ_LEN < MAX_SEQ_LEN is last row unless it is CIRCULAR */
                if( cself->mgr->reader_cols[ereflst_cn_SEQ_LEN].base.coord_len[0] < cself->mgr->max_seq_len ) {
                    if( !cself->circular ) {
                        break;
                    }
                    offset = 0;
                }
            } while(rc == 0 && q > 0 && len > 0 );
        }
    }
    ALIGN_DBGERR(rc);
    return rc;
}

typedef struct PlacementRecExtensionInfo PlacementRecExtensionInfo;
struct PlacementRecExtensionInfo
{
    /* data, destructor and size for extension 1 */
    void * data;
    void ( CC * destroy ) ( void *obj, void *data );
    size_t size;
};


LIB_EXPORT void * CC PlacementRecordCast ( const PlacementRecord *self, uint32_t ext )
{
    void * res = NULL;
    if ( self != NULL )
    {
        uint8_t * ptr = ( uint8_t * ) self;
        PlacementRecExtensionInfo * ext_info;
        switch( ext )
        {
            case placementRecordExtension0 : ptr += ( ( sizeof * self ) + ( 2 * ( sizeof * ext_info ) ) );
                                             break;

            case placementRecordExtension1 : ext_info = ( PlacementRecExtensionInfo * )( ptr + sizeof ( *self ) );
                                             ptr += ( ( sizeof * self ) + ( 2 * ( sizeof * ext_info ) ) + ext_info->size );
                                             break;
            default : ptr = NULL;
                      break;
        }
        res = ( void * )ptr;
    }
    return res;
}


LIB_EXPORT void CC PlacementRecordWhack( const PlacementRecord *cself )
{
    if ( cself != NULL ) 
    {
        PlacementRecord * self = ( PlacementRecord * )cself;
        PlacementRecExtensionInfo * ext_info;
        uint8_t * ptr = ( uint8_t * )self;
        ptr += sizeof( *self );
        ext_info = ( PlacementRecExtensionInfo * ) ptr;

        /* destroy from the outer callback-block beginning */
        if ( ext_info[ 1 ].destroy != NULL )
        {
            void *obj = PlacementRecordCast ( self, placementRecordExtension1 );
            ext_info[ 1 ].destroy( obj, ext_info[ 1 ].data );
        }

        if ( ext_info[ 0 ].destroy != NULL )
        {
            void *obj = PlacementRecordCast ( self, placementRecordExtension0 );
            ext_info[ 0 ].destroy( obj, ext_info[ 0 ].data );
        }
        /* now deallocate ( or put back into pool ) */
        free( self );
    }
}


struct PlacementIterator {
    const ReferenceObj* obj;
    INSDC_coord_zero ref_pos;
    INSDC_coord_zero ref_pos_end;
    INSDC_coord_zero ref_pos_overlaped;
    int32_t min_mapq;
    /* own reader in case of ref cursor based construction */
    const TableReader* ref_reader;
    TableReaderColumn* ref_cols;
    TableReaderColumn ref_cols_own[sizeof(ReferenceList_cols)/sizeof(ReferenceList_cols[0])];
    /* own reader in case of align cursor based construction */
    const TableReader* align_reader;
    TableReaderColumn* align_cols;
    TableReaderColumn align_cols_own[sizeof(PlacementIterator_cols)/sizeof(PlacementIterator_cols[0])];
    /* current reference table row */
    int64_t row;
    const TableReaderColumn* ids_col;
    Vector ids;
    /* PlacementRecord c-tor params */
    PlacementRecordExtendFuncs ext_0;
    PlacementRecordExtendFuncs ext_1;

    const VCursor* align_curs;
};


LIB_EXPORT rc_t CC ReferenceObj_MakePlacementIterator ( const ReferenceObj* cself,
    PlacementIterator **iter, INSDC_coord_zero ref_pos, INSDC_coord_len ref_len, int32_t min_mapq,
    struct VCursor const *ref_cur, struct VCursor const *align_cur, align_id_src ids,
    const PlacementRecordExtendFuncs *ext_0, const PlacementRecordExtendFuncs *ext_1 )
{
    rc_t rc = 0;
    PlacementIterator* o = NULL;

    if ( cself == NULL || iter == NULL )
    {
        rc = RC(rcAlign, rcType, rcAccessing, rcParam, rcInvalid);
    }
    else if ( ( rc = ReferenceSeq_ReOffset( cself->circular, cself->seq_len, &ref_pos ) ) != 0 )
    {
    }
    else if ( ( o = calloc( 1, sizeof( *o ) ) ) == NULL )
    {
        rc = RC( rcAlign, rcType, rcAccessing, rcMemory, rcExhausted );
    }
    else if ( ( rc = ReferenceList_AddRef( cself->mgr ) ) == 0 )
    {
        ReferenceList* mgr = cself->mgr;
        o->obj = cself;
        o->min_mapq = min_mapq;

        if ( ext_0 != NULL )
        {
            o->ext_0.data = ext_0->data;
            o->ext_0.destroy = ext_0->destroy;
            o->ext_0.populate = ext_0->populate;
            o->ext_0.alloc_size = ext_0->alloc_size;
            o->ext_0.fixed_size = ext_0->fixed_size;
        }

        if ( ext_1 != NULL )
        {
            o->ext_1.data = ext_1->data;
            o->ext_1.destroy = ext_1->destroy;
            o->ext_1.populate = ext_1->populate;
            o->ext_1.alloc_size = ext_1->alloc_size;
            o->ext_1.fixed_size = ext_1->fixed_size;
        }

        if ( ref_cur == NULL )
        {
            if ( mgr->reader != NULL || ( rc = ReferenceList_OpenCursor( mgr ) ) == 0 )
            {
                o->ref_reader = mgr->reader;
                o->ref_cols = mgr->reader_cols;
            }
        }
        else
        {
            memcpy( o->ref_cols_own, ReferenceList_cols, sizeof( o->ref_cols_own ) );
            o->ref_cols = o->ref_cols_own;
            rc = TableReader_MakeCursor( &o->ref_reader, ref_cur, o->ref_cols_own );
        }

        if ( align_cur == NULL )
        {
            if ( mgr->iter != NULL || ( rc = ReferenceList_OpenCursor2( mgr, ids ) ) == 0 )
            {
                o->align_reader = mgr->iter;
                o->align_cols = mgr->iter_cols;
            }
        }
        else
        {
            memcpy( o->align_cols_own, PlacementIterator_cols, sizeof( o->align_cols_own ) );
            o->align_cols = o->align_cols_own;
            o->align_curs = align_cur;
            rc = TableReader_MakeCursor( &o->align_reader, align_cur, o->align_cols );
        }

        if ( rc == 0 )
        {
            int64_t row = cself->start_rowid + ref_pos / mgr->max_seq_len;
            /* get effective starting offset based on overlap
               from alignments which started before the requested pos */
            if ( ( rc = TableReader_ReadRow( o->ref_reader, row ) ) == 0 )
            {
                if ( o->ref_cols[ereflst_cn_OVERLAP_REF_LEN].idx != 0 &&
                     o->ref_cols[ereflst_cn_OVERLAP_REF_LEN].base.coord_len[0] < (ref_pos % mgr->max_seq_len) )
                {
                    /* position beyond know overlap */
                    o->ref_pos_overlaped = ref_pos;
                }
                else if ( o->ref_cols[ereflst_cn_OVERLAP_REF_POS].idx != 0 )
                {
                    /* start with offset where earliest alignment begin */
                    o->ref_pos_overlaped = o->ref_cols[ereflst_cn_OVERLAP_REF_POS].base.coord0[0];
                }
                else
                {
                    /* default is step back 10 rows */
                    o->ref_pos_overlaped = 10 * mgr->max_seq_len;
                    if ( cself->circular )
                    {
                        if ( o->ref_pos_overlaped > cself->seq_len )
                        {
                            /* go back no more than one full length */
                            o->ref_pos_overlaped = cself->seq_len;
                        }
                        o->ref_pos_overlaped = ref_pos - o->ref_pos_overlaped; /* could become negative */
                    }
                    else
                    {
                        o->ref_pos_overlaped = ref_pos < o->ref_pos_overlaped ? 0 : (ref_pos - o->ref_pos_overlaped);
                    }
                }
                o->ref_pos = ref_pos;
                if ( !cself->circular && ref_pos + ref_len >= cself->seq_len )
                {
                    o->ref_pos_end = cself->seq_len - ref_pos - 1;
                }
                else
                {
                    o->ref_pos_end = ref_pos + ref_len - 1;
                }
                o->row = cself->start_rowid - 1 + o->ref_pos_overlaped / mgr->max_seq_len;
                VectorInit( &o->ids, 0, 100 );

                o->ids_col = &o->ref_cols[ids == primary_align_ids ? ereflst_cn_PRIMARY_ALIGNMENT_IDS :
                        ( ids == secondary_align_ids ? ereflst_cn_SECONDARY_ALIGNMENT_IDS : ereflst_cn_EVIDENCE_ALIGNMENT_IDS ) ];
            }
        }
    }

    if ( rc == 0 )
    {
        *iter = o;
        ALIGN_DBG( "iter for %s:%s opened 0x%p", cself->seqid, cself->name, o );
    }
    else
    {
        *iter = NULL;
        PlacementIteratorRelease(o);
        ALIGN_DBGERRP( "iter for %s:%s", rc, cself->seqid, cself->name );
    }
    return rc;
}

LIB_EXPORT rc_t CC PlacementIteratorAddRef ( const PlacementIterator *cself )
{
    return ReferenceList_AddRef(cself ? cself->obj->mgr : NULL);
}


static void CC PlacementIterator_whack_recs( void *item, void *data )
{
    PlacementRecordWhack( ( PlacementRecord * ) item );
}

LIB_EXPORT rc_t CC PlacementIteratorRelease ( const PlacementIterator *cself )
{
    if( cself != NULL ) {
        PlacementIterator* self = (PlacementIterator*)cself;

        VectorWhack( &self->ids, PlacementIterator_whack_recs, NULL );

        if( self->ref_reader != self->obj->mgr->reader ) {
            TableReader_Whack(self->ref_reader);
        }
        if( self->align_reader != self->obj->mgr->iter ) {
            TableReader_Whack(self->align_reader);
        }
        ReferenceList_Release(self->obj->mgr);
        free(self);
    }
    return 0;
}

LIB_EXPORT rc_t CC PlacementIteratorRefWindow(const PlacementIterator *self,
                                              const char **idstr, INSDC_coord_zero* pos, INSDC_coord_len* len)
{
    rc_t rc = 0;

    if( self == NULL || (idstr == NULL || pos == NULL || len == NULL) ) {
        rc = RC(rcAlign, rcType, rcAccessing, rcParam, rcInvalid);
    } else {
        if( idstr != NULL ) {
            *idstr = self->obj->seqid;
        }
        if( pos != NULL ) {
            *pos = self->ref_pos;
        }
        if( len != NULL ) {
            *len = self->ref_pos_end - self->ref_pos + 1;
        }
    }
    ALIGN_DBGERR(rc);
    return rc;
}


LIB_EXPORT rc_t CC PlacementIteratorRefObj( const PlacementIterator * self,
                                            struct ReferenceObj const ** refobj )
{
    rc_t rc = 0;

    if( self == NULL || refobj == NULL ) {
        rc = RC( rcAlign, rcType, rcAccessing, rcParam, rcInvalid );
    } else {
        *refobj = self->obj;
    }
    ALIGN_DBGERR( rc );
    return rc;
}


#if _DEBUGGING
static
void CC PlacementRecordVector_dump( void *item, void *data )
{
    const PlacementRecord* i = (const PlacementRecord*)item;
    ALIGN_DBG(" {%u, %u, %li}", i->pos, i->len, i->id);
}
#endif

static
int CC PlacementRecordVector_cmp(const void** left, const void** right, void* data)
{
    const PlacementRecord* l = *((const PlacementRecord**)left);
    const PlacementRecord* r = *((const PlacementRecord**)right);

    /* order by pos asc, len desc, id asc */
    int64_t d = (int64_t)(l->pos) - (int64_t)(r->pos);
    if( d == 0 ) {
        d = (int64_t)(r->len) - (int64_t)(l->len);
        if( d == 0 ) {
            d = (int64_t)(l->id) - (int64_t)(r->id);
        }
    }
    return d < 0 ? 1 : (d > 0 ? -1 : 0);
}


static rc_t allocate_populate_rec( const PlacementIterator *cself,
                                   PlacementRecord **rec,
                                   struct VCursor const *curs,
                                   int64_t id,
                                   INSDC_coord_zero apos,
                                   INSDC_coord_len alen )
{
    rc_t rc = 0;

    /* first check mapq, maybe we can skip the whole allocation... */
    if ( cself->align_cols[eplacementiter_cn_MAPQ].base.i32[ 0 ] < cself->min_mapq )
    {
        rc = RC( rcAlign, rcType, rcAccessing, rcId, rcIgnored );
    }
    else
    {
        PlacementRecExtensionInfo * ext_info;

        /* use callback or fixed size to discover the size of portions 0 and 1 */
        size_t size0 = ( ( cself->ext_0.alloc_size != NULL ) ?
                    cself->ext_0.alloc_size( curs, id, cself->ext_0.data ) : cself->ext_0.fixed_size );
        size_t size1 = ( ( cself->ext_1.alloc_size != NULL ) ?
                    cself->ext_1.alloc_size( curs, id, cself->ext_1.data ) : cself->ext_1.fixed_size );
        /* allocate the record ( or take it from a pool ) */
        *rec = calloc( 1, ( sizeof **rec ) + ( 2 * ( sizeof *ext_info ) ) + size0 + size1 );
        if ( *rec == NULL )
        {
            rc = RC( rcAlign, rcType, rcAccessing, rcMemory, rcExhausted );
        }
        else
        {
            uint8_t * ptr = ( uint8_t * )( * rec );
            ptr += sizeof ( **rec );
            ext_info = ( PlacementRecExtensionInfo * )ptr;

            /* prepopulate the core-record : */
            (*rec)->id  = id;               /* the row-id */
            (*rec)->ref = cself->obj;       /* the ReferenceObj it refers to */
            (*rec)->pos = apos;             /* the positon on the reference */
            (*rec)->len = alen;             /* the length on the reference */
            (*rec)->mapq = cself->align_cols[eplacementiter_cn_MAPQ].base.i32[ 0 ]; /* mapq */

            ext_info[ 0 ].data = cself->ext_0.data;          /* the opt. context ptr. */
            ext_info[ 0 ].destroy = cself->ext_0.destroy;    /* the opt. destructor */
            ext_info[ 0 ].size = size0;                      /* discovered size from above */

            ext_info[ 1 ].data = cself->ext_1.data;          /* the opt. context ptr. */
            ext_info[ 1 ].destroy = cself->ext_1.destroy;    /* the opt. destructor */
            ext_info[ 1 ].size = size1;                      /* discovered size from above */

            /* pass the record now to the opt. populate-callbacks */
            if ( cself->ext_0.populate != NULL )
            {
                void * obj = PlacementRecordCast ( *rec, placementRecordExtension0 );
                rc = cself->ext_0.populate( obj, *rec, curs,
                                            cself->ref_pos,
                                            cself->ref_pos_end - cself->ref_pos + 1,
                                            cself->ext_0.data );
            }

            if ( rc == 0 && cself->ext_1.populate != NULL )
            {
                void * obj = PlacementRecordCast ( *rec, placementRecordExtension1 );
                rc = cself->ext_1.populate( obj, *rec, curs, 
                                            cself->ref_pos,
                                            cself->ref_pos_end - cself->ref_pos + 1,
                                            cself->ext_1.data );
                if ( rc != 0 )
                {
                    /* destroy ext_0 */
                    if ( cself->ext_0.destroy != NULL )
                    {
                        void *obj = PlacementRecordCast ( *rec, placementRecordExtension0 );
                        cself->ext_0.destroy( obj, cself->ext_0.data );
                    }
                }
            }

            if ( rc != 0 )
            {
                /* free */
                free( *rec );
                *rec = NULL;
                /* let the caller continue with the next alignment */
                rc = RC( rcAlign, rcType, rcAccessing, rcId, rcIgnored );
            }
        }
    }
    return rc;
}

LIB_EXPORT rc_t CC PlacementIteratorNextAvailPos(const PlacementIterator *cself,
                                                 INSDC_coord_zero *pos, INSDC_coord_len *len)
{
    rc_t rc = 0;

    if( cself == NULL || (pos == NULL && len == NULL) ) {
        rc = RC(rcAlign, rcType, rcSelecting, rcParam, rcInvalid);
    } else {
        while( rc == 0 && VectorLength(&cself->ids) == 0 ) {
            PlacementIterator* self = (PlacementIterator*)cself;
            /* read ids */
            /*ALIGN_DBG("ref row: %li-%li-%li", self->obj->start_rowid, self->row + 1, self->obj->end_rowid);*/
            if( ++self->row > self->obj->end_rowid ) {
                rc = RC(rcAlign, rcType, rcSelecting, rcRange, rcDone);
            } else if( (rc = TableReader_ReadRow(cself->ref_reader, self->row)) == 0 ) {
                uint32_t i;
                /* fill out vector */
                /*ALIGN_DBG("align rows: %u", cself->ids_col->len);*/
                for ( i = 0; rc == 0 && i < cself->ids_col->len; i++ )
                {
                    if ( ( rc = TableReader_ReadRow( cself->align_reader, cself->ids_col->base.i64[i] ) ) == 0 )
                    {
                        INSDC_coord_zero apos = cself->align_cols[ eplacementiter_cn_REF_POS ].base.coord0[ 0 ];
                        INSDC_coord_len alen = self->align_cols[ eplacementiter_cn_REF_LEN ].base.coord_len[ 0 ];
                        /*ALIGN_DBG("align row: {%li, %u, %u}", cself->ids_col->base.i64[i], apos, alen);*/
			if(apos > cself->ref_pos_end){ 
				/*** we are beyond right edge of the window **/
				rc = RC(rcAlign, rcType, rcSelecting, rcRange, rcDone);
			} else if ( apos + alen < cself->ref_pos){
				 /* skip overlaps which do not reach original starting pos */
			} else {
                            PlacementRecord *rec;
                            rc = allocate_populate_rec( cself, &rec, cself->align_curs,
                                                        cself->ids_col->base.i64[i], apos, alen );
                            if ( rc == 0 )
                            {
                                /*ALIGN_DBG("align %p: {%li, %u, %u} - added[%u]", rec, cself->ids_col->base.i64[i],
                                    apos, alen, VectorLength(&cself->ids));*/
                                if( ( rc = VectorAppend( &self->ids, NULL, rec ) ) != 0 )
                                {
                                    PlacementRecordWhack(rec);
                                }
                            }
                            else
                            {
                                if ( GetRCState( rc ) == rcIgnored )
                                {
                                    rc = 0; /* do not break the loop if a record is filtered out! */
                                }
                            }
                        }

                    }
                }
                if( (rc == 0 || GetRCState(rc) == rcDone) && VectorLength(&cself->ids) > 0 ) {
                    VectorReorder(&self->ids, PlacementRecordVector_cmp, NULL);
#if _DEBUGGING
                    ALIGN_DBG("REFERENCE row %li %u recs order by pos asc, len desc, id asc",
                                cself->row, VectorLength(&cself->ids));
                    VectorForEach(&self->ids, true, PlacementRecordVector_dump, NULL);
#endif
                }
            }
        }
        if( rc == 0 || GetRCState(rc) == rcDone) {
            if( VectorLength(&cself->ids) > 0 ) {
                PlacementRecord* r = VectorLast(&cself->ids);
                if( pos != NULL ) {
                    *pos = r->pos;
                }
                if( len != NULL ) {
                    *len = r->len;
                }
                if ( r->pos > cself->ref_pos_end ) {
                    rc = RC(rcAlign, rcType, rcSelecting, rcRange, rcDone);
                }
            }
        }
    }
    if( GetRCState(rc) != rcDone ) {
        ALIGN_DBGERR(rc);
    }
    return rc;
}

LIB_EXPORT rc_t CC PlacementIteratorNextRecordAt(PlacementIterator *cself,
                                                 INSDC_coord_zero pos, const PlacementRecord **rec )
{
    rc_t rc = 0;
    if( cself == NULL || rec == NULL ) {
        rc = RC(rcAlign, rcType, rcAccessing, rcParam, rcInvalid);
    } else {
        uint32_t vlen = VectorLength(&cself->ids); 
        *rec = NULL;
        if( vlen > 0 ) {
            PlacementRecord* r = VectorLast(&cself->ids);
            if( r->pos == pos ) {
                VectorRemove(&((PlacementIterator*)cself)->ids, vlen - 1, (void**)rec);
            }
        }
    }
    if( rc == 0 && *rec == NULL ) {
        rc = RC(rcAlign, rcType, rcSelecting, rcOffset, rcDone);
    } else {
        ALIGN_DBGERR(rc);
    }
    return rc;
}

LIB_EXPORT rc_t CC PlacementIteratorNextIdAt(PlacementIterator *cself,
                                             INSDC_coord_zero pos, int64_t *row_id, INSDC_coord_len *len )
{
    rc_t rc = 0;
    const PlacementRecord* r = NULL;

    if( cself == NULL || row_id == NULL ) {
        rc = RC(rcAlign, rcType, rcAccessing, rcParam, rcInvalid);
    } else if( (rc = PlacementIteratorNextRecordAt(cself, pos, &r)) == 0 ) {
        *row_id = r->id;
        if( len != NULL ) {
            *len = r->len;
        }
        PlacementRecordWhack(r);
    }
    if( GetRCState(rc) != rcDone ) {
        ALIGN_DBGERR(rc);
    }
    return rc;
}
