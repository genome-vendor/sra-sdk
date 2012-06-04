/*==============================================================================
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
*/
#include <align/extern.h>

#include <kapp/main.h>
#include <klib/log.h>
#include <klib/rc.h>
#include <klib/sort.h>
#include <klib/data-buffer.h>
#include <klib/container.h>
#include <klib/checksum.h>
#include <kfs/mmap.h>
#include <kfs/file.h>
#include <kdb/manager.h>
#include <vdb/database.h>
#include <vdb/table.h>
#include <vdb/cursor.h>
#include <vdb/manager.h>
#include <vdb/vdb-priv.h>
#include <sra/sradb.h>

#include <align/writer-reference.h>
#include <align/writer-refseq.h>
#include <align/refseq-mgr.h>
#include "refseq-mgr-priv.h"
#include "writer-ref.h"
#include "reader-cmn.h"
#include "reference-cmn.h"
#include "debug.h"
#include <os-native.h>
#include <sysalloc.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

struct ReferenceMgr {
    const TableWriterRef* writer;
    const KDirectory* dir;
    BSTree tree;
    uint32_t options;
    const RefSeqMgr* rmgr;
    VDatabase* db;
    size_t cache;
    uint32_t num_open_max;
    uint32_t num_open;
    uint64_t usage;
    uint32_t max_seq_len;
    /* last seq end */
    int64_t ref_rowid;
    int32_t compress_buf[64 * 1024];
    /* must be last element of struct! */
    uint8_t seq_buf[TableWriterRefSeq_MAX_SEQ_LEN];
};

struct ReferenceSeq {
    BSTNode dad;
    const ReferenceMgr* mgr;
    bool local;
    char* id;
    char* accession;
    bool circular;
    union {
        struct {
            const KMMap* map;
            uint64_t file_size;
            const uint8_t* fasta_start;
            const uint8_t* fasta_end;
            uint64_t line_sz;
            uint8_t md5[16];
        } local;
        struct {
            const RefSeq* o;
            uint64_t usage;
        } refseq;
    } u;
    /* ref table position */
    int64_t start_rowid;
    /* total reference length */
    INSDC_coord_len seq_len;
};

static
int CC ReferenceSeq_Cmp(const void *item, const BSTNode *n)
{
    return strcasecmp((const char*)item, ((const ReferenceSeq*)n)->id);
}

static
int CC ReferenceSeq_Sort(const BSTNode *item, const BSTNode *n)
{
    return ReferenceSeq_Cmp(((const ReferenceSeq*)item)->id, n);
}

static
void CC ReferenceSeq_Unused( BSTNode *node, void *data )
{
    ReferenceSeq* n = (ReferenceSeq*)node;
    ReferenceSeq** d = (ReferenceSeq**)data;

    if( !n->local && n->u.refseq.o != NULL ) {
        if( *d == NULL || (*d)->u.refseq.usage > n->u.refseq.usage ) {
            *d = n;
        }
    }
}

static
rc_t ReferenceSeq_Alloc(ReferenceSeq** self, const char* id, size_t id_sz, const char* accession, size_t accession_sz)
{
    rc_t rc = 0;
    ReferenceSeq* obj;
    uint32_t sz = id_sz + 1 + accession_sz + 1;

    if( self == NULL || id == NULL || id_sz == 0 || accession == NULL || accession_sz == 0 ) {
        rc = RC(rcAlign, rcIndex, rcConstructing, rcParam, rcNull);
    } else if( (obj = calloc(1,  sizeof(*obj) + sz + (sizeof(size_t) - (sz % sizeof(size_t))))) == NULL ) {
        rc = RC(rcAlign, rcIndex, rcConstructing, rcMemory, rcExhausted);
    } else {
        obj->id = (char*)&obj[1];
        obj->accession = obj->id;
        obj->accession += id_sz + 1;
        memcpy(obj->id, id, id_sz);
        memcpy(obj->accession, accession, accession_sz);
        *self = obj;
    }
    return rc;
}

static
void CC ReferenceSeq_Whack(BSTNode *n, void *data)
{
    ReferenceSeq* self = (ReferenceSeq*)n;
    if( self->local ) {
        KMMapRelease(self->u.local.map);
        if( self->u.local.line_sz == 0 ) {
            free((uint8_t*)self->u.local.fasta_start);
        }
    } else {
        RefSeq_Release(self->u.refseq.o);
    }
    free(self);
}

struct OpenConfigFile_ctx {
    char const *name;
    KDirectory const *dir;
    KFile const **kfp;
    rc_t rc;
};

static
bool OpenConfigFile(char const server[], char const volume[], void *Ctx)
{
    struct OpenConfigFile_ctx *ctx = Ctx;
    KDirectory const *dir;
    
    if( volume == NULL ) {
        ctx->rc = KDirectoryOpenDirRead(ctx->dir, &dir, false, "%s", server);
    } else {
        ctx->rc = KDirectoryOpenDirRead(ctx->dir, &dir, false, "%s/%s", server, volume);
    }
    if (ctx->rc == 0) {
        ctx->rc = KDirectoryOpenFileRead(dir, ctx->kfp, ctx->name);
        KDirectoryRelease(dir);
        if (ctx->rc == 0) {
            return true;
        }
    }
    return false;
}

static
rc_t FindAndOpenConfigFile(const RefSeqMgr* rmgr, KDirectory const *dir, KFile const **kfp, char const conf[])
{
    rc_t rc = KDirectoryOpenFileRead(dir, kfp, conf);
    
    if(rc) {
        struct OpenConfigFile_ctx ctx;

        ctx.name = conf;
        ctx.dir = dir;
        ctx.kfp = kfp;
        ctx.rc = 0;
        
        rc = RefSeqMgr_ForEachVolume(rmgr, OpenConfigFile, &ctx);
        if (rc == 0 && *kfp == NULL) {
            rc = RC(rcAlign, rcIndex, rcConstructing, rcFile, rcNotFound);
        }
    }
    return rc;
}

static
rc_t ReferenceMgr_Conf(const RefSeqMgr* rmgr, BSTree* tree, KDirectory const *dir, char const conf[])
{
    rc_t rc = 0;
    const KFile* kf = NULL;

    assert(tree != NULL);
    assert(dir != NULL);

    if( conf != NULL &&
        (rc = FindAndOpenConfigFile(rmgr, dir, &kf, conf)) == 0 ) {
        const KMMap* mm;
        size_t page_sz;
        const char* str;

        if( (rc = KMMapMakeRead(&mm, kf)) == 0 ) {
            if( (rc = KMMapSize(mm, &page_sz)) == 0 && page_sz > 0 &&
                (rc = KMMapAddrRead(mm, (const void**)&str)) == 0 ) {
                const char* end = str + page_sz;
                while( rc == 0 && str < end ) {
                    const char* id, *accession;
                    size_t id_sz = 0, accession_sz = 0;

                    while( str < end && (str[0] == '\n' || str[0] == '\r') ) {
                        str++;
                    }
                    if( str >= end ) {
                        break;
                    }
                    id = str;
                    while( str < end && !isspace(str[0]) ) {
                        str++; id_sz++;
                    }
                    while( str < end && (str[0] == ' ' || str[0] == '\t') ) {
                        str++;
                    }
                    accession = str;
                    while( str < end && !isspace(str[0]) ) {
                        str++; accession_sz++;
                    }
                    if( id_sz == 0 || accession_sz == 0 ) {
                        rc = RC(rcAlign, rcIndex, rcConstructing, rcFileFormat, rcInvalid);
                        break;
                    }
                    /* skip to eol */
                    while( str < end && str[0] != '\n' && str[0] != '\r' ) {
                        str++;
                    }
                    if( RefSeqMgr_Exists(rmgr, accession, accession_sz, NULL) == 0 ) {
                        ReferenceSeq* tmp;
                        if( (rc = ReferenceSeq_Alloc(&tmp, id, id_sz, accession, accession_sz)) == 0 ) {
                            if( (rc = BSTreeInsertUnique(tree, &tmp->dad, NULL, ReferenceSeq_Sort)) != 0 ) {
                                ReferenceSeq_Whack(&tmp->dad, NULL);
                                break;
                            }
                            tmp->local = false;
                            ALIGN_DBG("RefSeq translation added '%s' -> '%s'", tmp->id, tmp->accession);
                        }
                    } else {
                        /* skips unknown references, later may be picked up as local files */
                        (void)PLOGMSG(klogWarn, (klogWarn, "RefSeq table '$(acc)' for '$(id)' was not found",
                                                           "acc=%.*s,id=%.*s", accession_sz, accession, id_sz, id));
                    }
                }
            }
            KMMapRelease(mm);
        }
        KFileRelease(kf);
    }
    return rc;
}

LIB_EXPORT rc_t CC ReferenceMgr_Make(const ReferenceMgr** cself, VDatabase* db, const VDBManager* vmgr,
                                     const uint32_t options, const char* conf, const char* path,
                                     uint32_t max_seq_len, size_t cache, uint32_t num_open)
{
    rc_t rc = 0;
    KDirectory* dir = NULL;
    ReferenceMgr* obj = NULL;
    uint32_t wopt = 0;
    
    wopt |= (options & ewrefmgr_co_allREADs) ? ewref_co_SaveRead : 0;
    wopt |= (options & ewrefmgr_co_Coverage) ? ewref_co_Coverage : 0;

    if( max_seq_len == 0 ) {
        max_seq_len = sizeof(obj->seq_buf);
    }
    if( cself == NULL ) {
        rc = RC(rcAlign, rcIndex, rcConstructing, rcParam, rcNull);
    } else if( (obj = calloc(1, sizeof(*obj) + max_seq_len - sizeof(obj->seq_buf))) == NULL ) {
        rc = RC(rcAlign, rcIndex, rcConstructing, rcMemory, rcExhausted);

    } else if( (rc = KDirectoryNativeDir(&dir)) == 0 ) {
        BSTreeInit(&obj->tree);
        obj->options = options;
        obj->cache = cache;
        obj->num_open_max = num_open;
        obj->max_seq_len = max_seq_len;
        if( rc == 0 && db != NULL && (rc = VDatabaseAddRef(db)) == 0 ) {
            obj->db = db;
        }
        if( rc == 0 ) {
            if( path == NULL ) {
                if( (rc = KDirectoryAddRef(dir)) == 0 ) {
                    obj->dir = dir;
                }
            } else {
                rc = KDirectoryOpenDirRead(dir, &obj->dir, false, path);
            }
            if( rc == 0 ) {
                rc = RefSeqMgr_Make(&obj->rmgr, vmgr, 0, cache, num_open);
            }
        }
        if (rc == 0 && (rc = ReferenceMgr_Conf(obj->rmgr, &obj->tree, obj->dir, conf)) != 0 ) {
            (void)PLOGERR(klogErr, (klogErr, rc, "failed to open configuration $(file)", "file=%s/%s", path ? path : ".", conf));
        }
    }
    KDirectoryRelease(dir);
    if( rc == 0 ) {
        *cself = obj;
        ALIGN_DBG("conf %s, local path '%s'", conf ? conf : "", path ? path : "");
    } else {
        ReferenceMgr_Release(obj, false, NULL, false);
        ALIGN_DBGERR(rc);
    }
    return rc;
}

#define ID_CHUNK_SZ 32
typedef struct TChunk_struct {
    struct TChunk_struct* next;
    int64_t id[ID_CHUNK_SZ];
} TChunk;

typedef struct {
    TChunk* head;
    TChunk* tail;
    uint8_t tail_qty; /* id count within tail */
    uint32_t id_qty; /* id count overall */
    uint32_t mismatches;
    uint32_t indels;
    uint8_t* hilo; /* array of max_seq_len */
} TCover;

static
void ReferenceMgr_TCoverRelease(TCover* c)
{
    TChunk* h = c->head;
    while( h != NULL ) {
        TChunk* n = h->next;
        free(h);
        h = n;
    }
    c->head = NULL;
}

static
rc_t ReferenceMgr_ReCover(const ReferenceMgr* cself, uint64_t rows)
{
    rc_t rc = 0;
    int i;
    uint64_t new_rows = 0;
    const TableWriterRefCoverage* cover_writer = NULL;

    TableReaderColumn acols[] =
    {
        {0, "REF_ID", {NULL}, 0, 0},
        {0, "REF_START", {NULL}, 0, 0},
        {0, "REF_LEN", {NULL}, 0, 0},
        {0, "(INSDC:dna:text)MISMATCH", {NULL}, 0, 0},
        {0, "REF_OFFSET", {NULL}, 0, 0},
        {0, "HAS_REF_OFFSET", {NULL}, 0, 0},
        {0, NULL, {NULL}, 0, 0}
    };
    const int64_t** al_ref_id = &acols[0].base.i64;
    const INSDC_coord_zero** al_ref_start = &acols[1].base.coord0;
    const INSDC_coord_len** al_ref_len = &acols[2].base.coord_len;
    const TableReaderColumn* mismatch_len =  &acols[3];
    const TableReaderColumn* ref_offset = &acols[4];
    const TableReaderColumn* has_ref_offset = &acols[5];

    /* TODO TO DO TBD change to work properly with mutiple reads in alignment table */

    if( (rc = TableWriterRefCoverage_Make(&cover_writer, cself->db, 0)) == 0 ) {
        /* order is important see ReferenceSeqCoverage struct */
        const char* tbls[] = {"PRIMARY_ALIGNMENT", "SECONDARY_ALIGNMENT", "EVIDENCE_ALIGNMENT"};
#define tbls_qty (sizeof(tbls)/sizeof(tbls[0]))
        rc_t rc1 = 0;
        int64_t ref_from = 1, rr;
        struct {
            int64_t* v;
            size_t sz;
            uint64_t q;
        } ids[tbls_qty];

        ALIGN_DBG("sizeof(TCover) %u, sizeof(TChunk) %u", sizeof(TCover), sizeof(TChunk));
        memset(ids, 0, sizeof(ids));
        while( rc == 0 && ref_from <= rows ) {
            uint32_t i;
            uint64_t ref_qty = rows;
            TCover* data = NULL;
            uint8_t* hilo;
            /* allocate mem for up to ref_qty of reference coverage data up to 4Gb */
#if ARCH_BITS == 32
const uint32_t RECOVER_MAX_MEM = ( 1UL << 30 ); /* 1 GB */
#else
const uint64_t RECOVER_MAX_MEM = ( 4UL << 30 ); /* 4 GB */
#endif
            do {
                data = NULL;
                if( (ref_qty * (sizeof(*data) + cself->max_seq_len)) > RECOVER_MAX_MEM ||
                    (data = calloc(ref_qty, sizeof(*data) + cself->max_seq_len)) == NULL ) {
                    ref_qty = ref_qty / 10 * 9;
                }
                if( ref_qty < 1 ) {
                    rc = RC(rcAlign, rcTable, rcCommitting, rcMemory, rcExhausted);
                }
                hilo = (uint8_t*)&data[ref_qty];
            } while( rc == 0 && data == NULL );
            /* grep through tables for coverage data */
            ALIGN_DBG("covering REFERENCE rowid range [%ld:%ld]", ref_from, ref_from + ref_qty - 1);
            for(i = 0; rc == 0 && i < tbls_qty; i++) {
                const VTable* table = NULL;
                const TableReader* reader = NULL;
                int64_t al_from;
                uint64_t al_qty;

                ALIGN_DBG("covering REFERENCE with %s", tbls[i]);
                if( (rc = VDatabaseOpenTableRead(cself->db, &table, tbls[i])) != 0 ) {
                    if( GetRCObject(rc) == rcTable && GetRCState(rc) == rcNotFound ) {
                        ALIGN_DBG("table %s was not found, ignored", tbls[i]);
                        rc = 0;
                        continue;
                    } else {
                        break;
                    }
                }
                if( (rc = TableReader_Make(&reader, table, acols, cself->cache)) == 0 &&
                    (rc = TableReader_IdRange(reader, &al_from, &al_qty)) == 0 ) {
                    int64_t al_rowid;

                    for(al_rowid = al_from; rc == 0 && al_rowid < al_from + al_qty; al_rowid++) {
                        int64_t ref_r, al_ref_id_end;
                        uint64_t j, seq_start = 0, refseq_start;

                        if( (rc = TableReader_ReadRow(reader, al_rowid)) != 0 ) {
                            break;
                        }
                        /* alignment can run across multiple reference chunks */
                        al_ref_id_end = **al_ref_id + (**al_ref_start + **al_ref_len) / cself->max_seq_len;
                        refseq_start = **al_ref_start;
                        ALIGN_DBG("al row %ld has REF_ID [%ld,%ld], REF_START: %i, REF_LEN %u",
                                   al_rowid, **al_ref_id, al_ref_id_end, **al_ref_start, **al_ref_len);
                        for(ref_r = **al_ref_id; ref_r <= al_ref_id_end; ref_r++) {
                            uint64_t refseq_len = cself->max_seq_len - refseq_start;
                            if( seq_start + refseq_len > **al_ref_len ) {
                                refseq_len = **al_ref_len - seq_start;
                            }
                            ALIGN_DBG("covered ref_id %ld [%ld,%ld]", ref_r, refseq_start, refseq_start + refseq_len);
                            if( ref_r >= ref_from && ref_r < ref_from + ref_qty ) {
                                uint64_t k = ref_r - ref_from;
                                ALIGN_DBG("%ld is a match for %ld[%lu]", ref_r, ref_from + k, k);
                                if( data[k].tail_qty == ID_CHUNK_SZ || data[k].head == NULL ) {
                                    TChunk* x = malloc(sizeof(*(data[k].tail)));
                                    while( x == NULL ) {
                                        /* release last ref cover_writer record and retry */
                                        ReferenceMgr_TCoverRelease(&data[--ref_qty]);
                                        ALIGN_DBG("downsize covering REFERENCE rowid range [%ld:%ld] from %s",
                                                  ref_from, ref_from + ref_qty - 1, tbls[i]);
                                        if( ref_qty < 1 ) {
                                            rc = RC(rcAlign, rcTable, rcCommitting, rcMemory, rcExhausted);
                                            break;
                                        } else if( ref_r >= ref_from + ref_qty ) {
                                            break;
                                        }
                                    }
                                    if( rc != 0 || x == NULL ) {
                                        break;
                                    }
                                    if( data[k].head == NULL ) {
                                        data[k].head = x;
                                    } else {
                                        data[k].tail->next = x;
                                    }
                                    x->next = NULL;
                                    data[k].tail = x;
                                    data[k].tail_qty = 0;
                                    data[k].hilo = &hilo[cself->max_seq_len * k];
                                }
                                /* 2 bits designated src table */
                                data[k].tail->id[data[k].tail_qty] = i;
                                data[k].tail->id[data[k].tail_qty] <<= 62;
                                data[k].tail->id[data[k].tail_qty] = al_rowid;
                                data[k].tail_qty++;
                                data[k].id_qty++;
                                if( ref_r == **al_ref_id ) {
                                    /* write those to 1st chunk only */
                                    data[k].mismatches += mismatch_len->len;
                                    data[k].indels += ref_offset->len;
                                    if( ref_offset->len > 0 && has_ref_offset->base.buul[0] && ref_offset->base.i32[0] < 0 ) {
                                        data[k].indels--;
                                    }
                                }
                                for(j = refseq_start; j < refseq_start + refseq_len; j++) {
                                    if( data[k].hilo[j] < UINT8_MAX ) {
                                        data[k].hilo[j]++;
                                    }
                                }
                            }
                            seq_start += refseq_len;
                            refseq_start = 0;
                        }
                        rc = rc ? rc : Quitting();
                    }
                }
                TableReader_Whack(reader);
                VTableRelease(table);
            }
            /* prep and write coverage data */
            for(rr = 0; rc == 0 && rr < ref_qty; rr++) {
                TChunk* x;
                ReferenceSeqCoverage c;

                ALIGN_DBGF(("ref rowid %ld:", ref_from + rr));
                memset(&c, 0, sizeof(c));
                x = data[rr].head;
                while( rc == 0 && rr < ref_qty && x != NULL ) {
                    uint32_t q = data[rr].tail == x ? data[rr].tail_qty : ID_CHUNK_SZ;
                    for(i = 0; i < q; i++) {
                        int t = (x->id[i] >> 62) & 0x3;
#if _ARCH_BITS == 32
                        x->id[i] &= 0x3FFFFFFFFFFFFFULL;
#else
                        x->id[i] &= 0x3FFFFFFFFFFFFFUL;
#endif
                        assert(t >= 0 && t < tbls_qty);
                        ALIGN_DBGF((" %lu[%c][%i]", x->id[i], t));
                        while( ids[t].sz <= ids[t].q ) {
                            int64_t* n = realloc(ids[t].v, (ids[t].sz += ID_CHUNK_SZ) * sizeof(*ids[t].v));
                            if( n == NULL ) {
                                ReferenceMgr_TCoverRelease(&data[--ref_qty]);
                                ALIGN_DBG("downsize covering REFERENCE rowid range [%ld:%ld]", ref_from, ref_from + ref_qty - 1);
                                if( ref_qty < 1 ) {
                                    rc = RC(rcAlign, rcTable, rcCommitting, rcMemory, rcExhausted);
                                    break;
                                } else if( rr >= ref_qty ) {
                                    break;
                                }
                            } else {
                                ids[t].v = n;
                            }
                        }
                        if( rc != 0 || rr >= ref_qty ) {
                            break;
                        }
                        ids[t].v[ids[t].q++] = x->id[i];
                    }
                    x = x->next;
                }
                ALIGN_DBGF(("\n"));
                for(i = 0; rc == 0 && i < tbls_qty; i++) {
                    c.ids[i].buffer = ids[i].v;
                    c.ids[i].elements = ids[i].q;
                    ALIGN_DBG(" %i qty %lu,", i, c.ids[i].elements);
                }
                if( rc == 0 ) {
                    c.mismatches = data[rr].mismatches;
                    c.indels = data[rr].indels;
                    if( data[rr].hilo != NULL ) {
                        memset(&c.low, 0xFF, sizeof(c.low));
                        for(i = 0; i < cself->max_seq_len; i++) {
                            if( c.high < data[rr].hilo[i] ) {
                                c.high = data[rr].hilo[i];
                            }
                            if( c.low > data[rr].hilo[i] ) {
                                c.low = data[rr].hilo[i];
                            }
                        }
                    }
                    rc = TableWriterRefCoverage_Write(cover_writer, ref_from + rr, &c);
                }
                ReferenceMgr_TCoverRelease(&data[rr]);
            }
            for(; rr < ref_qty; rr++) {
                /* clean up in case of rc != 0 */
                ReferenceMgr_TCoverRelease(&data[rr]);
            }
            free(data);
            ref_from += ref_qty;
            ref_qty = rows - ref_from + 1;
        }
        for(i = 0; i < tbls_qty; i++) {
            free(ids[i].v);
        }
        rc1 = TableWriterRefCoverage_Whack(cover_writer, rc == 0, &new_rows);
        rc = rc ? rc : rc1;
        if( rc == 0 && rows != new_rows ) {
            rc = RC(rcAlign, rcTable, rcCommitting, rcData, rcInconsistent);
        }
    }
    return rc;
}

LIB_EXPORT rc_t CC ReferenceMgr_Release(const ReferenceMgr* cself, bool commit, uint64_t* rows, bool build_coverage)
{
    rc_t rc = 0;
    if( cself != NULL ) {
        uint64_t rr, *r = rows ? rows : &rr;
        ReferenceMgr* self = (ReferenceMgr*)cself;

        rc = TableWriterRef_Whack(self->writer, commit, r);
        BSTreeWhack(&self->tree, ReferenceSeq_Whack, NULL);
        KDirectoryRelease(self->dir);
        if( rc == 0 && build_coverage && commit ) {
            VTable* t = NULL;
            rc = ReferenceMgr_ReCover(cself, *r);
            if( VDatabaseOpenTableUpdate(self->db, &t, "SECONDARY_ALIGNMENT") == 0 ) {
                VTableDropColumn(t, "TMP_MISMATCH");
            }
            VTableRelease(t);
        }
        VDatabaseRelease(self->db);
        RefSeqMgr_Release(self->rmgr);
        free(self);
    }
    return rc;
}

static
rc_t ReferenceMgr_ImportFastaFile(const ReferenceMgr* cself, const KFile* kf, const char* id, ReferenceSeq** seq)
{
    rc_t rc = 0;
    const KMMap* map = NULL;
    const uint8_t* addr;
    size_t file_size;
    MD5State md5;

    assert(cself != NULL);
    assert(kf != NULL);
    if( seq != NULL ) {
        assert(id != NULL);
    }

    if( (rc = KMMapMakeRead(&map, kf)) == 0 &&
        (rc = KMMapAddrRead(map, (const void**)&addr)) == 0 && (rc = KMMapSize(map, &file_size)) == 0  ) {
        ReferenceMgr* self = (ReferenceMgr*)cself;
        const uint8_t* end = &addr[file_size];
        ReferenceSeq* obj = NULL;
        while( rc == 0 && addr < end ) {
            size_t line_sz, line_no, bad_line_1st;
            /* skip empty space between sequences */
            while( addr < end && isspace(addr[0])) {
                addr++;
            }
            if( addr < end && addr[0] == '>' ) {
                const uint8_t* nl, *sp = addr;
                while( sp < end && !isspace((++sp)[0]) );
                nl = sp;
                while( nl[0] != '\n' && nl[0] != '\r' && nl < end ) { nl++; }
                while( (nl[0] == '\n' || nl[0] == '\r') && nl < end ) { nl++; }
                if( sp >= end || nl >= end ) {
                    rc = RC(rcAlign, rcFile, rcReading, rcFormat, rcInvalid);
                } else if( id != NULL && obj != NULL) {
                    /* call with id means only one seq per this file is expected */
                    rc = RC(rcAlign, rcFile, rcReading, rcItem, rcUnexpected);
                } else {
                    const char* id2 = id;
                    size_t id_sz;
                    if( id == NULL ) {
                        id2 = (const char*)(&addr[1]);
                        id_sz = sp - addr - 1;
                    } else {
                        id_sz = strlen(id);
                    }
                    if( (rc = ReferenceSeq_Alloc(&obj, id2, id_sz, (const char*)&addr[1], sp - addr - 1)) == 0 &&
                        (rc = KMMapAddRef(map)) == 0 ) {
                        obj->local = true;
                        obj->mgr = cself;
                        obj->u.local.file_size = file_size;
                        obj->u.local.map = map;
                        obj->u.local.fasta_start = nl;
                        obj->u.local.fasta_end = nl + 1;
                        MD5StateInit(&md5);
                    }
                }
                addr = nl;
            } else if( addr < end ) {
                rc = RC(rcAlign, rcFile, rcReading, rcFormat, rcInvalid);
            }
            line_sz = line_no = bad_line_1st = 0;
            while( rc == 0 && addr < end ) {
                /* read, validate and md5 the sequence */
                if( addr[0] == '\n' || addr[0] == '\r' ) {
                    addr++;
                    line_no++;
                    if( obj->u.local.line_sz == 0 ) {
                        /* remember 1st line len */
                        obj->u.local.line_sz = line_sz;
                    } else if( line_sz != obj->u.local.line_sz ) {
                        if( bad_line_1st == 0 ) {
                            bad_line_1st = line_no;
                        }
                    }
                    line_sz = 0;
                } else if( addr[0] == '>' ) {
                    break;
                } else if( strchr(INSDC_4na_accept_CHARSET, addr[0]) == NULL ) {
                    rc = RC(rcAlign, rcFile, rcReading, rcFormat, rcInvalid);
                } else {
                    INSDC_4na_bin c = addr[0];
                    c = (c == '.') ? 'N' : toupper(c);
                    line_sz++;
                    MD5StateAppend(&md5, &c, 1);
                    obj->seq_len++;
                    obj->u.local.fasta_end = ++addr;
                }
            }
            if( rc == 0 && obj != NULL ) {
                ReferenceSeq* existing = NULL;
                MD5StateFinish(&md5, obj->u.local.md5);
                if( obj->u.local.line_sz == 0 ) {
                    /* if line_sz still not set, set to last */
                    obj->u.local.line_sz = line_sz;
                }
                if( obj->seq_len == 0 ) {
                    rc = RC(rcAlign, rcFile, rcReading, rcItem, rcIncomplete);
                } else if( (rc = BSTreeInsertUnique(&self->tree, &obj->dad, (BSTNode**)&existing, ReferenceSeq_Sort)) == 0 ) {
                    ALIGN_DBG("RefSeq local fasta added '%s' -> '%s'", obj->id, obj->accession);
                    if( bad_line_1st != 0 && bad_line_1st < line_no ) {
                        uint8_t* x = malloc(obj->seq_len);
                        /* need to make copy of malformed fasta to mem */
                        ALIGN_DBG("making memcpy since format is broken for '%s'", obj->id);
                        if( x == NULL ) {
                            rc = RC(rcAlign, rcFile, rcReading, rcMemory, rcExhausted);
                        } else {
                            INSDC_coord_len l;
                            obj->u.local.line_sz = 0;
                            if( (rc = ReferenceSeq_Read(obj, 0, obj->seq_len, x, &l)) == 0 ) {
                                assert( obj->seq_len == l );
                                obj->u.local.fasta_start = x;
                                obj->u.local.fasta_end = obj->u.local.fasta_start + obj->seq_len;
                            } else {
                                free(x);
                            }
                        }
                    }
                } else if( rc != 0 ) {
                    if( existing != NULL ) {
                        const uint8_t* rmd5;
                        INSDC_coord_len len;
                        
                        rc = 0;
                        if( existing->local ) {
                            if( obj->seq_len != existing->seq_len ) {
                                rc = RC(rcAlign, rcFile, rcReading, rcItem, rcDuplicate);
                                ALIGN_DBG("%s length mismatch %u <> %u", obj->id, obj->seq_len, existing->seq_len);
                            } else if( memcmp(obj->u.local.md5, existing->u.local.md5, sizeof(obj->u.local.md5)) != 0 ) {
                                rc = RC(rcAlign, rcFile, rcReading, rcItem, rcDuplicate);
                                ALIGN_DBG("%s local md5 mismatch", obj->id);
                            }
                        } else {
                            const RefSeq* r = existing->u.refseq.o;
                            if( r == NULL ) {
                                rc = RefSeqMgr_GetSeq(cself->rmgr, &r, existing->accession, strlen(existing->accession));
                            }
                            if( rc != 0 ) {
                                ALIGN_DBGERRP("%s refseq opening", rc, existing->id);
                            } else if( (rc = RefSeq_SeqLength(r, &len)) == 0 && obj->seq_len != len ) {
                                ALIGN_DBG("%s refseq and %s length mismatch", existing->id, obj->id);
                                rc = RC(rcAlign, rcFile, rcReading, rcItem, rcDuplicate);
                            } else if( (rc = RefSeq_MD5(r, &rmd5)) == 0 && memcmp(obj->u.local.md5, rmd5, sizeof(obj->u.local.md5)) != 0 ) {
                                ALIGN_DBGF(("%s refseq and %s md5 mismatch '", existing->id, obj->id));
                                for(len = 0; len < sizeof(obj->u.local.md5); len++) {
                                    ALIGN_DBGF(("%02hx", obj->u.local.md5[len]));
                                }
                                ALIGN_DBGF(("' <> '"));
                                for(len = 0; len < sizeof(obj->u.local.md5); len++) {
                                    ALIGN_DBGF(("%02hx", rmd5[len]));
                                }
                                ALIGN_DBGF(("'\n"));
                                rc = RC(rcAlign, rcFile, rcReading, rcItem, rcDuplicate);
                            }
                            if( r != existing->u.refseq.o ) {
                                RefSeq_Release(r);
                            }
                        }
                        ReferenceSeq_Whack(&obj->dad, NULL);
                        if( rc == 0 ) {
                            ALIGN_DBG("note: %s full duplicate", obj->id);
                            obj = existing;
                        }
                    } else {
                        ReferenceSeq_Whack(&obj->dad, NULL);
                    }
                }
                if( rc == 0 && seq != NULL ) {
                    *seq = obj;
                }
            }
        }
    }
    KMMapRelease(map);
    ALIGN_DBGERR(rc);
    return rc;
}

static
rc_t ReferenceSeq_ReadDirect(ReferenceSeq* self, INSDC_coord_zero offset, INSDC_coord_len len, bool read_circular,
                             uint8_t* buffer, INSDC_coord_len* written)
{
    rc_t rc = 0;

    *written = 0;
    if( len == 0 ) {
    } else if( (rc = ReferenceSeq_ReOffset(self->circular || read_circular, self->seq_len, &offset)) != 0 ) {
    } else if( self->local ) {
        /* translate offset into file map dimensions */
        const uint8_t* start = self->u.local.fasta_start + offset;
        start += self->u.local.line_sz ? offset / self->u.local.line_sz : 0; /* add \n on each line */
        assert(start < self->u.local.fasta_end);
        do {
            while( len > *written && start < self->u.local.fasta_end ) {
                while( start < self->u.local.fasta_end && isspace(*start) ) {
                    start++;
                }
                while( len > *written && start < self->u.local.fasta_end && !isspace(*start) ) {
                    buffer[*written] = *start;
                    *written = *written + 1;
                    start++;
                }
            }
            if( read_circular && start >= self->u.local.fasta_end ) {
                start = self->u.local.fasta_start;
            }
        } while( len > *written && start < self->u.local.fasta_end );
    } else {
        /* we need to trim len to actual length of seq */
        if( !read_circular && (offset + len) >= self->seq_len ) {
            len = self->seq_len - offset;
        }
        if( rc == 0 && (rc = RefSeq_Read(self->u.refseq.o, offset, len, buffer, written)) == 0 ) {
            self->u.refseq.usage = ++((ReferenceMgr*)self->mgr)->usage;
        }
    }
    ALIGN_DBGERR(rc);
    return rc;
}

static
rc_t CC ReferenceMgr_Find(const ReferenceMgr* cself, const char* id, ReferenceSeq** seq, const KFile** kf)
{
    rc_t rc = 0;

    assert(cself != NULL && id != NULL && seq != NULL && kf != NULL);

    *seq = (ReferenceSeq*)BSTreeFind(&cself->tree, id, ReferenceSeq_Cmp);
    *kf = NULL;
    if( *seq == NULL ) {
        /* ask efSeqMgr */
        size_t id_len = strlen(id);
        if( RefSeqMgr_Exists(cself->rmgr, id, id_len, NULL) == 0 ) {
            ReferenceSeq* obj;
            if( (rc = ReferenceSeq_Alloc(&obj, id, id_len, id, id_len)) == 0 ) {
                if( (rc = BSTreeInsertUnique((BSTree*)&cself->tree, &obj->dad, NULL, ReferenceSeq_Sort)) == 0 ) {
                    *seq = obj;
                    ALIGN_DBG("RefSeq added on-fly '%s' -> '%s'", (*seq)->id, (*seq)->accession);
                } else {
                    ReferenceSeq_Whack(&obj->dad, NULL);
                }
            }
        } else {
            /* try local file */
            if( KDirectoryOpenFileRead(cself->dir, kf, "%s.fasta", id) == 0 ) {
                ALIGN_DBG("RefSeq added from file %s.fasta", id);
            } else if( KDirectoryOpenFileRead(cself->dir, kf, "%s.fa", id) == 0 ) {
                ALIGN_DBG("RefSeq added from file %s.fa", id);
            } else {
                rc = RC(rcAlign, rcFile, rcOpening, rcFile, rcNotFound);
            }
        }
    }
    ALIGN_DBGERR(rc);
    return rc;
}

LIB_EXPORT rc_t CC ReferenceMgr_GetSeq(const ReferenceMgr* cself, const ReferenceSeq** seq, const char* id)
{
    rc_t rc = 0;
    ReferenceSeq* obj;
    
    if( cself == NULL || seq == NULL || id == NULL ) {
        rc = RC(rcAlign, rcFile, rcConstructing, rcParam, rcNull);
    } else {
        ReferenceMgr* mgr = (ReferenceMgr*)cself;
        const KFile* kf;
        
        *seq = NULL;
        rc = ReferenceMgr_Find(cself, id, &obj, &kf);
        if( rc == 0 && obj == NULL ) {
            /* it is local file */
            rc = ReferenceMgr_ImportFastaFile(cself, kf, id, &obj);
            KFileRelease(kf);
        }
        if( rc == 0 && !obj->local && obj->u.refseq.o == NULL ) {
            if( cself->num_open_max > 0 && cself->num_open >= cself->num_open_max ) {
                ReferenceSeq* old = NULL;
                BSTreeForEach(&cself->tree, false, ReferenceSeq_Unused, &old);
                if( old != NULL ) {
                    RefSeq_Release(old->u.refseq.o);
                    old->u.refseq.o = NULL;
                    mgr->num_open--;
                }
            }
            rc = RefSeqMgr_GetSeq(cself->rmgr, &obj->u.refseq.o, obj->accession, strlen(obj->accession));
            obj->mgr = cself;
            if( rc == 0 &&
                (rc = RefSeq_Circular(obj->u.refseq.o, &obj->circular)) == 0 &&
                (rc = RefSeq_SeqLength(obj->u.refseq.o, &obj->seq_len)) == 0 ) {
                mgr->num_open++;
            }
        }
        if( rc == 0 && obj->start_rowid == 0 ) {
            /* append to the whole thing to REFERENCE table since we encounter it for the 1st time */
            TableWriterRefData data;
            INSDC_coord_len len = 0;
            INSDC_coord_zero offset = 0;
            obj->start_rowid = mgr->ref_rowid + 1;
            data.name.buffer = obj->id;
            data.name.elements = strlen(obj->id);
            data.read.buffer = mgr->seq_buf;
            data.seq_id.buffer = obj->accession;
            data.seq_id.elements = strlen(obj->accession);
            data.force_READ_write = obj->local || (cself->options & ewrefmgr_co_allREADs);
            data.circular = obj->circular;
            if (cself->writer == NULL) {
                uint32_t wopt = 0;

                wopt |= (mgr->options & ewrefmgr_co_allREADs) ? ewref_co_SaveRead : 0;
                wopt |= (mgr->options & ewrefmgr_co_Coverage) ? ewref_co_Coverage : 0;
                if( (rc = TableWriterRef_Make(&mgr->writer, mgr->db, wopt)) == 0 ) {
                    TableWriterData mlen;
                    mlen.buffer = &mgr->max_seq_len;
                    mlen.elements = 1;
                    rc = TableWriterRef_WriteDefaultData(mgr->writer, ewrefseq_cn_MAX_SEQ_LEN, &mlen);
                }
            }
            if (rc == 0) {
                do {
                    if( (rc = ReferenceSeq_ReadDirect(obj, offset, mgr->max_seq_len, false,
                                                      mgr->seq_buf, &len)) == 0 && len > 0 ) {
                        data.read.elements = len;
                        rc = TableWriterRef_Write(cself->writer, &data, NULL);
                        offset += len;
                        mgr->ref_rowid++;
                    }
                } while( rc == 0 && len > 0 && offset < obj->seq_len );
            }
        }
    }
    if( rc == 0 ) {
        *seq = obj;
    } else {
        ALIGN_DBGERR(rc);
    }
    return rc;
}

LIB_EXPORT rc_t CC ReferenceMgr_Verify(const ReferenceMgr* cself, const char* id, INSDC_coord_len seq_len, const uint8_t* md5)
{
    rc_t rc = 0;
    ReferenceSeq* rseq;
    const KFile* kf;
    
    if( cself == NULL || id == NULL ) {
        rc = RC(rcAlign, rcFile, rcConstructing, rcParam, rcNull);
    } else if( (rc = ReferenceMgr_Find(cself, id, &rseq, &kf)) == 0 ) {
        if( rseq == NULL ) {
            uint64_t size = 0;
            if( (rc = KFileSize(kf, &size)) == 0 && size == 0 ) {
                rc = RC(rcAlign, rcTable, rcValidating, rcSize, rcEmpty);
            }
            KFileRelease(kf);
        } else if( rseq->local ) {
            if( rseq->seq_len != seq_len ) {
                rc = RC(rcAlign, rcFile, rcValidating, rcSize, rcUnequal);
                ALIGN_DBGERRP("%s->%s SEQ_LEN verification", rc, id, rseq->accession);
            }
            if( rc == 0 && md5 != NULL && memcmp(md5, rseq->u.local.md5, sizeof(rseq->u.local.md5)) != 0 ) {
                unsigned i;
                rc = RC(rcAlign, rcTable, rcValidating, rcChecksum, rcUnequal);
                ALIGN_DBGERRP("%s->%s MD5 verification", rc, id, rseq->accession);
                ALIGN_DBGF((" local '"));
                for(i = 0; i < sizeof(rseq->u.local.md5); i++) {
                    ALIGN_DBGF(("%02hx", rseq->u.local.md5[i]));
                }
                ALIGN_DBGF(("'  <> match '"));
                for(i = 0; i < sizeof(rseq->u.local.md5); i++) {
                    ALIGN_DBGF(("%02hx", md5[i]));
                }
                ALIGN_DBGF(("'\n"));
            } else {
                ALIGN_DBG("%s->%s MD5 verification ok", id, rseq->accession);
            }
        } else {
            const RefSeq* tmp = rseq->u.refseq.o;
            INSDC_coord_len o_len;
            const uint8_t* o_md5;

            if( tmp == NULL ) {
                if( (rc = RefSeqMgr_GetSeq(cself->rmgr, &tmp, rseq->accession, strlen(rseq->accession))) != 0 ||
                    (rc = RefSeq_SeqLength(tmp, &o_len)) != 0 ) {
                    ALIGN_DBGERRP("%s->%s verification", rc, id, rseq->accession);
                }
            } else {
                o_len = rseq->seq_len;
            }
            if( rc == 0 ) {
                if( seq_len != 0 && o_len != seq_len ) {
                    rc = RC(rcAlign, rcTable, rcValidating, rcSize, rcUnequal);
                    ALIGN_DBGERRP("%s->%s SEQ_LEN verification", rc, id, rseq->accession);
                } else {
                    ALIGN_DBG("%s->%s SEQ_LEN verification ok", id, rseq->accession);
                }
            }
            if( rc == 0 && (rc = RefSeq_MD5(tmp, &o_md5)) == 0 ) {
                if( md5 != NULL && o_md5 != NULL && memcmp(md5, o_md5, sizeof(rseq->u.local.md5)) != 0 ) {
                    rc = RC(rcAlign, rcTable, rcValidating, rcChecksum, rcUnequal);
                    ALIGN_DBGERRP("%s->%s MD5 verification", rc, id, rseq->accession);
                } else {
                    ALIGN_DBG("%s->%s MD5 verification ok", id, rseq->accession);
                }
            }
            if( tmp != rseq->u.refseq.o ) {
                RefSeq_Release(tmp);
            }
        }
    }
    if( rc == 0 ) {
        ALIGN_DBG("%s verification ok", id);
    } else {
        ALIGN_DBGERRP("%s verification", rc, id);
    }
    return rc;
}

LIB_EXPORT rc_t CC ReferenceMgr_FastaPath(const ReferenceMgr* cself, const char* fasta_path)
{
    rc_t rc = 0;
    KDirectory* dir;

    if( cself == NULL || fasta_path == NULL ) {
        rc = RC(rcAlign, rcFile, rcConstructing, rcParam, rcNull);
    } else if( (rc = KDirectoryNativeDir(&dir)) == 0 ) {
        const KFile* kf;
        if( (rc = KDirectoryOpenFileRead(dir, &kf, "%s", fasta_path)) == 0 ) {
            rc = ReferenceMgr_FastaFile(cself, kf);
            KFileRelease(kf);
        }
        KDirectoryRelease(dir);
    }
    ALIGN_DBGERRP("from file %s", rc, fasta_path);
    return rc;
}

LIB_EXPORT rc_t CC ReferenceMgr_FastaFile(const ReferenceMgr* cself, const KFile* file)
{
    if( cself == NULL || file == NULL ) {
        return RC(rcAlign, rcFile, rcConstructing, rcParam, rcNull);
    }
    return ReferenceMgr_ImportFastaFile(cself, file, NULL, NULL);
}

static
rc_t cigar2offset(uint32_t options, const void* cigar, uint32_t cigar_len, int32_t* out, uint32_t out_sz, uint32_t out_used,
                  INSDC_coord_len* seq_len, INSDC_coord_len* ref_len, INSDC_coord_len* max_ref_len)
{
    rc_t rc = 0;
    uint32_t i, last_match = 0, op_len, last_soft_len = 0;
    unsigned char op;
    static char const op_txt[] = "MIDNSHP=XB";

#if _DEBUGGING
    ALIGN_C_DBGF(("%s:%u: ", __func__, __LINE__));
    if( options & ewrefmgr_cmp_Binary ) {
        ALIGN_C_DBGF(("bin cigar '"));
        for(i = 0; i < cigar_len; i++) {
            const uint32_t* c = cigar;
            op = c[i] & 0x0F;
            op_len = c[i] >> 4;
            ALIGN_C_DBGF(("%u%c", op_len, op < 10 ? op_txt[op] : '?'));
        }
    } else {
        ALIGN_C_DBGF(("cigar '%.*s", cigar_len, cigar));
    }
    ALIGN_C_DBGF(("'[%u]\n", cigar_len));
#endif

    *seq_len = 0;
    *ref_len = 0;
    *max_ref_len = 0;
    memset(out, 0, out_used * sizeof(*out) * 2);
    for(i = 0; rc == 0 && i < cigar_len; i++) {
        if( options & ewrefmgr_cmp_Binary ) {
            const uint32_t* c = cigar;
            op = c[i] & 0x0F;
            if( op < 10 ) {
                op = op_txt[op];
            }
            op_len = c[i] >> 4;
        } else {
            const unsigned char* c = cigar;
            op_len = 0;
            while( c[i] != '\0' && isdigit(c[i]) ) {
                op_len *= 10;
                op_len += c[i++] - '0';
            }
            /* avoid intersecting with 4-bit binary above */
            op = c[i] <= 0x0F ? 0xFF : c[i];
        }
        if (op == 'P') {
            continue;
        }
        if( op != 'S' && op != 'I' ) {
            last_soft_len = 0;
        }
        switch(op) {
            case 'M':
            case '=':
            case 'X':
                *seq_len += op_len;
                *ref_len += op_len;
                *max_ref_len += op_len;
                last_match = *seq_len;
                break;
            case 'S':
            case 'I':
                last_soft_len += op_len;
                *seq_len += op_len;
                *max_ref_len += op_len;
                if( last_match < out_sz ) {
                    out[last_match] += -op_len;
                } else {
                    rc = RC(rcAlign, rcFile, rcProcessing, rcData, rcInconsistent);
                }
                break;
            case 'B':
                /* Complete Genomics CIGAR style specific:
                   overlap between consecutive reads
                   ex: sequence 6 bases: ACACTG, reference 2 bases: ACTG,
                   cigar will be: 2M2B2M
                   no need to move sequence position
                */
                if( last_match < out_sz ) {
                    out[last_match] -= op_len;
                } else {
                    rc = RC(rcAlign, rcFile, rcProcessing, rcData, rcInconsistent);
                }
                if( *ref_len < op_len ) {
                    *ref_len = 0;
                } else {
                    *ref_len -= op_len;
                }
                *max_ref_len += op_len;
                break;
            case 'D':
            case 'N':
                if( last_match < out_sz ) {
                    out[last_match] += op_len;
                } else {
                    rc = RC(rcAlign, rcFile, rcProcessing, rcData, rcInconsistent);
                }
                *ref_len += op_len;
                *max_ref_len += op_len;
                break;
            case 'H':
            default:
                rc = RC(rcAlign, rcFile, rcProcessing, rcData, rcUnrecognized);
                break;
        }
    }
    if( !(options & ewrefmgr_cmp_Exact) && last_soft_len > 0 ) {
        ALIGN_C_DBG("trailing soft clip %u @ %u", last_soft_len, last_match);
        out[last_match] += last_soft_len;
    }
#if _DEBUGGING
    ALIGN_C_DBGF(("%s:%u: offsets: ", __func__, __LINE__));
    last_match = last_match > *seq_len ? last_match : *seq_len;
    for(i = 0; i <= last_match; i++) {
        if( i == *seq_len ) {
            ALIGN_C_DBGF((" | "));
        }
        ALIGN_C_DBGF(("%i%c", out[i], ((i + 1) % 5) ? ',' : ' '));
    }
    ALIGN_C_DBGF((" cigar derived seq_len %u, ref_len %u, max_ref_len %u\n", *seq_len, *ref_len, *max_ref_len));
#endif
    return rc;
}

LIB_EXPORT rc_t CC ReferenceMgr_Compress(const ReferenceMgr* cself, uint32_t options,
                                         const char* id, INSDC_coord_zero offset,
                                         const char* seq, INSDC_coord_len seq_len,
                                         const void* cigar, uint32_t cigar_len,
                                         INSDC_coord_zero allele_offset, const char* allele, INSDC_coord_len allele_len,
                                         const void* allele_cigar, uint32_t allele_cigar_len,
                                         TableWriterAlgnData* data)
{
    rc_t rc = 0;
    const ReferenceSeq* refseq;

    if( cself == NULL || id == NULL ) {
        rc = RC(rcAlign, rcFile, rcProcessing, rcParam, rcNull);
    } else if( (rc = ReferenceMgr_GetSeq(cself, &refseq, id)) == 0 ) {
        rc = ReferenceSeq_Compress(refseq, options, offset, seq, seq_len, cigar, cigar_len,
                                   allele_offset, allele, allele_len, allele_cigar, allele_cigar_len, data);
        ReferenceSeq_Release(refseq);
    }
    ALIGN_C_DBGERR(rc);
    return rc;
}

LIB_EXPORT rc_t CC ReferenceSeq_Compress(const ReferenceSeq* cself, uint32_t options, INSDC_coord_zero offset,
                                         const char* seq, INSDC_coord_len seq_len,
                                         const void* cigar, uint32_t cigar_len,
                                         INSDC_coord_zero allele_offset, const char* allele, INSDC_coord_len allele_len,
                                         const void* allele_cigar, uint32_t allele_cigar_len,
                                         TableWriterAlgnData* data)
{
    rc_t rc = 0;

    if( cself == NULL || seq == NULL || cigar == NULL || cigar_len == 0 || data == NULL ||
        (!(allele == NULL && allele_len == 0 && allele_cigar == NULL && allele_cigar_len == 0) &&
         !(allele != NULL && allele_cigar != NULL && allele_cigar_len != 0)) ) {
        rc = RC(rcAlign, rcFile, rcProcessing, rcParam, rcInvalid);
    } else if( seq_len > sizeof(cself->mgr->compress_buf) / sizeof(cself->mgr->compress_buf[0]) ) {
        rc = RC(rcAlign, rcFile, rcProcessing, rcBuffer, rcInsufficient);
    } else if( (rc = ReferenceSeq_ReOffset(cself->circular, cself->seq_len, &offset)) == 0 ) {
        INSDC_coord_len i, seq_pos = 0, allele_ref_end = 0, ref_len = 0, rl = 0, max_rl = 0;
        INSDC_coord_zero* read_start = &((INSDC_coord_zero*)(data->read_start.buffer))[data->ploidy];
        INSDC_coord_len* read_len = &((INSDC_coord_len*)(data->read_len.buffer))[data->ploidy];
        bool* has_ref_offset, *has_mismatch;
        int32_t* ref_offset;
        uint8_t* mismatch, ref_buf[64 * 1024];
#if _DEBUGGING
        uint64_t i_ref_offset_elements, i_mismatch_elements;
        char x[4096];
#endif

        if( data->ploidy == 0 ) {
            data->has_ref_offset.elements = seq_len;
            data->ref_offset.elements = 0;
            data->has_mismatch.elements = seq_len;
            data->mismatch.elements = 0;
            *read_start = 0;
        } else {
            data->has_ref_offset.elements += seq_len;
            data->has_mismatch.elements += seq_len;
            *read_start = read_start[-1] + read_len[-1];
        }
        *read_len = seq_len;
        has_ref_offset = &((bool*)(data->has_ref_offset.buffer))[*read_start];
        ref_offset = (int32_t*)(data->ref_offset.buffer);
        has_mismatch = &((bool*)(data->has_mismatch.buffer))[*read_start];
        mismatch = (uint8_t*)(data->mismatch.buffer);

#if _DEBUGGING
        i_ref_offset_elements = data->ref_offset.elements;
        i_mismatch_elements = data->mismatch.elements;
        ALIGN_C_DBG("align%s '%.*s'[%u] to '%s:%s' at %i", (options & ewrefmgr_cmp_Exact) ? " EXACT" : "",
                    seq_len, seq, seq_len, cself->id, cself->accession, offset);
#endif
        if( allele != NULL ) {
            /* determine length of reference for subst by allele */
            ALIGN_C_DBG("apply allele %.*s[%u] at %i w/cigar below",
                        allele_len, allele, allele_len, allele_offset);
            rc = cigar2offset(options, allele_cigar, allele_cigar_len,
                    (int32_t*)(cself->mgr->compress_buf), sizeof(cself->mgr->compress_buf) / sizeof(cself->mgr->compress_buf[0]),
                    allele_len, &seq_pos, &allele_ref_end, &max_rl);
            /* where allele ends on reference */
            allele_ref_end += allele_offset;
        }
        if( rc == 0 ) {
            rc = cigar2offset(options, cigar, cigar_len,
                             (int32_t*)(cself->mgr->compress_buf), sizeof(cself->mgr->compress_buf) / sizeof(cself->mgr->compress_buf[0]),
                             seq_len, &seq_pos, &rl, &max_rl);
        }
        if( allele != NULL ) {
            if( allele_offset + allele_ref_end < offset || allele_offset >= offset + rl ) {
                (void)PLOGMSG(klogWarn, (klogWarn,
                    "allele $(a) offset $(ao) $(ac) is not within referenced region in $(id) at offset $(ro) $(rc)",
                    "a=%.*s,ao=%i,ac=%.*s,id=%s,ro=%i,rc=%.*s",
                    allele_len, allele, allele_offset, (options & ewrefmgr_cmp_Binary) ? 0 : allele_cigar_len, allele_cigar,
                    cself->accession, offset, (options & ewrefmgr_cmp_Binary) ? 0 : cigar_len, cigar));
                allele = NULL;
            }
        }
        if( rc == 0 ) {
            if( !(options & ewrefmgr_cmp_Exact) ) {
                /* do not write leading deletes just move reference offset */
                if( cself->mgr->compress_buf[0] > 0 ) {
                    offset += cself->mgr->compress_buf[0];
                    rl -= cself->mgr->compress_buf[0];
                    max_rl -= cself->mgr->compress_buf[0];
                    ((int32_t*)(cself->mgr->compress_buf))[0] = 0;
                    ALIGN_C_DBG("adjusted offset %i", offset);
                }
            }
            ref_len = rl;
            if( (offset + max_rl) > cself->seq_len && !cself->circular ) {
                max_rl = cself->seq_len - offset;
                if( max_rl < rl ) {
                    /* ref_len used for compression cannot be shorter than ref_len derived from cigar,
                       if there is a shortage it will fail later here */
                    max_rl = rl;
                }
                ALIGN_C_DBG("max_ref_len truncated to %u cause it goes beyond refseq length %lu at offset %i",
                             max_rl, cself->seq_len, offset);
            }
            ALIGN_C_DBG("chosen REF_LEN %u, ref len for match %u", ref_len, max_rl);

            assert(seq_len == seq_pos);
            if( max_rl > sizeof(ref_buf) ) {
                rc = RC(rcAlign, rcFile, rcProcessing, rcBuffer, rcInsufficient);
            } else {
                int64_t ref_pos;

                if( allele != NULL ) {
                    /* subst allele in reference */
                    if( allele_offset < offset ) {
                        /* move allele start inside referenced chunk */
                        if( allele_len < offset - allele_offset ) {
                            /* shouldnt happen, there is a check + warning above */
                            rc = RC(rcAlign, rcFile, rcProcessing, rcSize, rcInsufficient);
                        } else {
                            allele += offset - allele_offset;
                            allele_len -= offset - allele_offset;
                            rl = 0;
                        }
                    } else {
                        /* fetch portion of reference which comes before allele */
                        rl = allele_offset - offset;
                        rc = ReferenceSeq_ReadDirect((ReferenceSeq*)cself, offset, rl, true, ref_buf, &i);
                        if( rc == 0 && rl != i ) {
                            /* here we need to test it otherwise excessive portion of allele could be fetch in next if */
                            rc = RC(rcAlign, rcFile, rcProcessing, rcRange, rcExcessive);
                        }
                    }
                    if( rc == 0 && allele_len < (max_rl - rl) ) {
                        memcpy(&ref_buf[rl], allele, allele_len);
                        rl += allele_len;
                        /* append tail of actual reference */
                        rc = ReferenceSeq_ReadDirect((ReferenceSeq*)cself, allele_ref_end, max_rl - rl, true, &ref_buf[rl], &i);
                        rl += i;
                    } else if( rc == 0 ) {
                        /* allele is longer than needed */
                        memcpy(&ref_buf[rl], allele, max_rl - rl);
                        rl = max_rl;
                    }
                } else {
                    rc = ReferenceSeq_ReadDirect((ReferenceSeq*)cself, offset, max_rl, true, ref_buf, &rl);
                }
                if( rc != 0 || max_rl != rl ) {
                    rc = rc ? rc : RC(rcAlign, rcFile, rcProcessing, rcRange, rcExcessive);
                    ALIGN_C_DBGERRP("refseq is shorter: at offset %i need %u bases", rc, offset, max_rl);
                } else {
                    for(seq_pos = 0, ref_pos = 0; seq_pos < seq_len; seq_pos++, ref_pos++) {
                        has_ref_offset[seq_pos] = cself->mgr->compress_buf[seq_pos] != 0;
                        if( has_ref_offset[seq_pos] ) {
                            ref_offset[data->ref_offset.elements++] = cself->mgr->compress_buf[seq_pos];
                            ref_pos += cself->mgr->compress_buf[seq_pos];
                        }
                        if( (ref_pos < 0) || (ref_pos >= max_rl) || 
                            ((toupper(ref_buf[ref_pos]) != toupper(seq[seq_pos])) && (seq[seq_pos] != '=')) ) {
                            has_mismatch[seq_pos] = true;
                            mismatch[data->mismatch.elements++] = seq[seq_pos];
                        } else {
                            /* match! */
                            has_mismatch[seq_pos] = false;
                        }
                    }
                }
            }
        }
#if _DEBUGGING
        if(rc == 0 ) {
            INSDC_coord_len j;
            memset(x, '-', sizeof(x) - 1);
            x[sizeof(x) - 2] = '\0';

            ALIGN_C_DBG("ref: %.*s [%u]", max_rl, ref_buf, max_rl);
            ALIGN_C_DBGF(("%s:%u: ref: ", __func__, __LINE__));
            for(seq_pos = 0, j = 0, rl = 0, i = 0; seq_pos < seq_len; seq_pos++, j++) {
                if( has_ref_offset[seq_pos] ) {
                    if( ref_offset[i_ref_offset_elements + rl] > 0 ) {
                        ALIGN_C_DBGF(("%.*s", (uint32_t)(ref_offset[i_ref_offset_elements + rl]), &ref_buf[j]));
                    } else {
                        i = -ref_offset[i_ref_offset_elements + rl];
                    }
                    j += ref_offset[i_ref_offset_elements + rl];
                    rl++;
                }
                ALIGN_C_DBGF(("%c", (j < 0 || j >= max_rl) ? '-' : (i > 0) ? tolower(ref_buf[j]) : ref_buf[j]));
                if( i > 0  ) {
                    i--;
                }
            }
            ALIGN_C_DBGF(("\n%s:%u: seq: ", __func__, __LINE__));
            for(i = 0, j = 0; i < seq_len; i++) {
                if( has_ref_offset[i] && ref_offset[i_ref_offset_elements + j++] > 0 ) {
                    ALIGN_C_DBGF(("%.*s", ref_offset[i_ref_offset_elements + j - 1], x));
                }
                ALIGN_C_DBGF(("%c", seq[i]));
            }
            ALIGN_C_DBGF((" [%u]\n", seq_len));
            ALIGN_C_DBGF(("%s:%u: hro: ", __func__, __LINE__));
            for(i = 0, j = 0; i < seq_len; i++) {
                if( has_ref_offset[i] && ref_offset[i_ref_offset_elements + j++] > 0 ) {
                    ALIGN_C_DBGF(("%.*s", ref_offset[i_ref_offset_elements + j - 1], x));
                }
                ALIGN_C_DBGF(("%c", has_ref_offset[i] + '0'));
            }
            ALIGN_C_DBGF((", ro:"));
            for(i = i_ref_offset_elements; i < data->ref_offset.elements; i++) {
                ALIGN_C_DBGF((" %i,", ref_offset[i]));
            }
            ALIGN_C_DBGF(("[%u]\n", data->ref_offset.elements - i_ref_offset_elements));
            ALIGN_C_DBGF(("%s:%u: hmm: ", __func__, __LINE__));
            for(i = 0, j = 0; i < seq_len; i++) {
                if( has_ref_offset[i] && ref_offset[i_ref_offset_elements + j++] > 0 ) {
                    ALIGN_C_DBGF(("%.*s", ref_offset[i_ref_offset_elements + j - 1], x));
                }
                ALIGN_C_DBGF(("%c", has_mismatch[i] + '0'));
            }
            ALIGN_C_DBGF((", mm: '%.*s'[%u]\n", (int)(data->mismatch.elements - i_mismatch_elements),
                &mismatch[i_mismatch_elements], data->mismatch.elements - i_mismatch_elements));
        }
#endif
        if( rc == 0 && allele == NULL ) {
            if( data->ploidy == 0 ) {
                data->ref_1st_row_id = cself->start_rowid;
                data->effective_offset = offset;
                data->ref_len = ref_len;
                ALIGN_C_DBGF(("%s:%u: reference 1st ROW_ID %li OFFSET %i REF_LEN %u",
                    __func__, __LINE__, data->ref_1st_row_id, data->effective_offset, data->ref_len));
                if( data->ref_id.buffer ) {
                    ((int64_t*)(data->ref_id.buffer))[0] = cself->start_rowid + offset / cself->mgr->max_seq_len;
                    data->ref_id.elements = 1;
                    ALIGN_C_DBGF((" REF_ID %li", ((int64_t*)(data->ref_id.buffer))[0]));
                }
                if( data->ref_start.buffer ) {
                    ((INSDC_coord_zero*)(data->ref_start.buffer))[0] = offset % cself->mgr->max_seq_len;
                    data->ref_start.elements = 1;
                    ALIGN_C_DBGF((" REF_START %i", ((INSDC_coord_zero*)(data->ref_start.buffer))[0]));
                }
                if( data->global_ref_start.buffer ) {
                    ((uint64_t*)(data->global_ref_start.buffer))[0] = (cself->start_rowid - 1) * cself->mgr->max_seq_len + offset;
                    data->global_ref_start.elements = 1;
                    ALIGN_C_DBGF((" GLOBAL_REF_START %lu", ((uint64_t*)(data->global_ref_start.buffer))[0]));
                }
                ALIGN_C_DBGF(("\n"));
            } else {
                if( data->ref_1st_row_id != cself->start_rowid || data->effective_offset != offset ) {
                    rc = RC(rcAlign, rcFile, rcProcessing, rcData, rcInconsistent);
                    (void)PLOGERR(klogErr, (klogErr, rc,
                        "all reads in alignment record must align to same refseq at same location $(r1)@$(o1) <> $(r2):$(a2)@$(o2)",
                        "r1=%li,o1=%i,r2=%s,a2=%s,o2=%i", data->ref_1st_row_id, data->effective_offset, cself->id, cself->accession, offset));
                } else if( data->ref_len != ref_len ) {
                    rc = RC(rcAlign, rcFile, rcProcessing, rcData, rcInconsistent);
                    (void)PLOGERR(klogErr, (klogErr, rc,
                        "all reads in alignment record must have same size projection on refseq $(rl1) <> $(rl2) $(r):$(a)@$(o)",
                        "rl1=%u,rl2=%u,r=%s,a=%s,o=%i", data->ref_len, ref_len, cself->id, cself->accession, offset));
                }
            }
        }
        if( rc == 0 ) {
            data->ploidy++;
            data->read_start.elements = data->ploidy;
            data->read_len.elements = data->ploidy;
        }
    }
    ALIGN_C_DBGERR(rc);
    return rc;
}

LIB_EXPORT rc_t CC ReferenceSeq_Read(const ReferenceSeq* cself, INSDC_coord_zero offset, INSDC_coord_len len,
                                     uint8_t* buffer, INSDC_coord_len* ref_len)
{
    rc_t rc = 0;

    if( cself == NULL || buffer == NULL || ref_len == NULL ) {
        rc = RC(rcAlign, rcFile, rcReading, rcParam, rcNull);
    } else {
        rc = ReferenceSeq_ReadDirect((ReferenceSeq*)cself, offset, len, true, buffer, ref_len);
    }
    ALIGN_DBGERR(rc);
    return rc;
}

LIB_EXPORT rc_t CC ReferenceSeq_Get1stRow(const ReferenceSeq* cself, int64_t* row_id)
{
    rc_t rc = 0;

    if( cself == NULL || row_id == NULL ) {
        rc = RC(rcAlign, rcFile, rcReading, rcParam, rcNull);
    } else {
        *row_id = cself->start_rowid;
    }
    return rc;
}

LIB_EXPORT rc_t CC ReferenceSeq_AddCoverage(const ReferenceSeq* cself, INSDC_coord_zero offset, const ReferenceSeqCoverage* data)
{
    rc_t rc = 0;

    if( cself == NULL || data == NULL) {
        rc = RC(rcAlign, rcFile, rcReading, rcParam, rcNull);
    } else if( !(cself->mgr->options & ewrefmgr_co_Coverage) ) {
        rc = RC( rcAlign, rcType, rcWriting, rcData, rcUnexpected);
        ALIGN_DBGERRP("coverage %s", rc, "data");
    } else if( (rc = ReferenceSeq_ReOffset(cself->circular, cself->seq_len, &offset)) == 0 ) {
        rc = TableWriterRef_WriteCoverage(cself->mgr->writer, cself->start_rowid, offset, data);
    }
    ALIGN_DBGERR(rc);
    return rc;
}

LIB_EXPORT rc_t CC ReferenceSeq_Release(const ReferenceSeq *cself)
{
    return 0;
}
