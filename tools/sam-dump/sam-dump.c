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
#include <klib/report.h> /* ReportInit */
#include <klib/container.h>
#include <klib/log.h>
#include <klib/out.h>
#include <klib/status.h>
#include <klib/rc.h>
#include <klib/vector.h>
#include <klib/printf.h>
#include <kfs/file.h>
#include <kfs/buffile.h>
#include <kfs/gzip.h>
#include <kfs/bzip.h>
#include <kdb/meta.h>
#include <kdb/namelist.h>
#include <kapp/main.h>
#include <kapp/args.h>
#include <insdc/insdc.h>
#include <insdc/sra.h>
#include <vdb/report.h>
#include <vdb/manager.h>
#include <vdb/database.h>
#include <vdb/table.h>
#include <vdb/cursor.h>
#include <vdb/vdb-priv.h>
#include <vdb/schema.h>
#include <vdb/dependencies.h>
#include <sra/sraschema.h>
#include <sra/srapath.h>
#include <align/dna-reverse-cmpl.h>
#include <align/iterator.h>
#include <align/reference.h>

#include <kfs/directory.h>
#include <os-native.h>
#include <sysalloc.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strtol.h>
#include <ctype.h>
#include <assert.h>

#include "debug.h"
#include "sam-dump.vers.h"

#if _ARCH_BITS == 64
#define USE_MATE_CACHE 1
#else
#define USE_MATE_CACHE 0
#endif

#define CURSOR_CACHE (32 * 1024 * 1024)

typedef struct TAlignedRegion_struct {
    char name[1024];
    struct {
        INSDC_coord_zero from;
        INSDC_coord_zero to;
    } r[10240];
    int rq;
    INSDC_coord_zero max_to;
} TAlignedRegion;

typedef struct TMatepairDistance_struct {
    uint64_t from;
    uint64_t to;
} TMatepairDistance;

struct {
    const char* accession;
    const char* path;
    bool use_seqid;
    bool unaligned;
    bool long_cigar;
    bool reheader;
    bool noheader;
    bool hide_identical;
    bool fasta;
    bool fastq;
    bool spot_group_in_name;
    bool reverse_unaligned;
    bool only_primaries;
    int cg_style;
    const char* name_prefix;
    /* region filter data */
    TAlignedRegion* region;
    uint32_t region_qty;
    /* distance filter data */
    bool mp_dist_unknown;
    TMatepairDistance* mp_dist;
    uint32_t mp_dist_qty;
    uint32_t test_rows;
    const ReferenceList* ref_list;
    /* mate info cache */
    int64_t mate_row_gap_cachable;
    uint32_t comments_qty;
    const char** comments;
} param;

typedef union UData_union {
    const void* v;
    const uint32_t* u32;
    const int32_t* i32;
    const int64_t* i64;
    const uint64_t* u64;
    const uint8_t* u8;
    const char* str;
    INSDC_coord_one* coord1;
    INSDC_coord_zero* coord0;
    INSDC_coord_len* coord_len;
    INSDC_coord_val* coord_val;
    INSDC_SRA_xread_type* read_type;
    INSDC_SRA_read_filter* read_filter;
} UData;

typedef struct SCol_struct {
    const char* name;
    uint32_t idx;
    UData base;
    uint32_t len;
    bool optional;
} SCol;

typedef struct STable_struct {
    const char* name;
    const VTable* vtbl;
} STable;

typedef struct SCursCache_struct {
    KVector* cache;
    uint32_t sam_flags;
    INSDC_coord_zero pnext;
    int32_t tlen;
    const ReferenceObj* ref;
    /* cache stats */
    uint64_t projected;
    uint64_t added;
    uint64_t hit;
    uint64_t bad;
} SCursCache;

typedef struct SCurs_struct {
    const STable* tbl;
    const VCursor* vcurs;
    SCursCache* cache;
    SCursCache cache_local;
    uint64_t col_reads_qty;
} SCurs;

enum ealg_col {
    alg_SEQ_NAME = 0,
    alg_MAPQ,
    alg_CIGAR,
    alg_READ,
    alg_SAM_QUALITY,
    alg_SPOT_GROUP,
    alg_SEQ_SPOT_GROUP,
    alg_SEQ_READ_ID,
    alg_EDIT_DISTANCE,
    alg_MATE_ALIGN_ID,
    alg_MATE_REF_NAME,
    alg_SAM_FLAGS,
    alg_SEQ_SPOT_ID,
    alg_MATE_REF_POS,
    alg_REGION_FILTER,
    alg_REF_NAME = alg_REGION_FILTER,
    alg_REF_SEQ_ID,
    alg_REF_POS,
    alg_REF_LEN,
    alg_DISTANCE_FILTER,
    alg_TEMPLATE_LEN = alg_DISTANCE_FILTER,
    alg_CG_GC,
    alg_CG_GS,
    alg_CG_GQ
};

SCol g_alg_col_tmpl[] = {
    {"SEQ_NAME", 0, {NULL}, 0, false},
    {"MAPQ", 0, {NULL}, 0, false},
    {"?CIGAR column name?", 0, {NULL}, 0, false},
    {"?READ column name?", 0, {NULL}, 0, false},
    {"SAM_QUALITY", 0, {NULL}, 0, false},
    {"SPOT_GROUP", 0, {NULL}, 0, true},
    {"SEQ_SPOT_GROUP", 0, {NULL}, 0, true},
    {"SEQ_READ_ID", 0, {NULL}, 0, true},
    {"EDIT_DISTANCE", 0, {NULL}, 0, false},
    /* start cols used as standalone in DumpUnaligned */
    /* MATE_ALIGN_ID col must preceeed
       MATE_REF_NAME, MATE_REF_POS, SAM_FLAGS, TEMPLATE_LEN for cache to work */
    {"MATE_ALIGN_ID", 0, {NULL}, 0, false},
    {"?MATE_REF_NAME column name?", 0, {NULL}, 0, false},
    {"SAM_FLAGS", 0, {NULL}, 0, false},
    {"SEQ_SPOT_ID", 0, {NULL}, 0, true},
    {"MATE_REF_POS", 0, {NULL}, 0, false},
    /* these are read before any other for filtering so they must be last */
    {"REF_NAME", 0, {NULL}, 0, false},
    {"REF_SEQ_ID", 0, {NULL}, 0, false},
    {"REF_POS", 0, {NULL}, 0, false},
    /* end cols used as standalone in DumpUnaligned */
    {"REF_LEN", 0, {NULL}, 0, false},
    {"TEMPLATE_LEN", 0, {NULL}, 0, false},
    {NULL, 0, {NULL}, 0, false}, /* GS tag */
    {NULL, 0, {NULL}, 0, false}, /* GC tag */
    {NULL, 0, {NULL}, 0, false}, /* GQ tag */
    {NULL, 0, {NULL}, 0, false}
};

static
rc_t RefSeqPrint(void)
{
    rc_t rc = 0;
    uint32_t i, count = 0;
    
    rc = ReferenceList_Count(param.ref_list, &count);
    for(i = 0; rc == 0 && i < count; i++) {
        const ReferenceObj* obj;
        if( (rc = ReferenceList_Get(param.ref_list, &obj, i)) == 0 ) {
            const char* name = NULL;
            const char* seqid = NULL;
            INSDC_coord_len len;
            if( (rc = ReferenceObj_SeqId(obj, &seqid)) == 0 && (rc = ReferenceObj_Name(obj, &name)) == 0 &&
                (rc = ReferenceObj_SeqLength(obj, &len)) == 0 ) {
                const char* nm;
                if( param.use_seqid && seqid != NULL && seqid[0] != '\0' ) {
                    nm = seqid;
                } else {
                    nm = name;
                }
                OUTMSG(("@SQ\tSN:%s", nm));
                if( nm != seqid && seqid != NULL && seqid[0] != '\0' && strcmp(nm, seqid) != 0 ) {
                    OUTMSG(("\tAS:%s", seqid));
                }
                OUTMSG(("\tLN:%u\n", len));
            }
            ReferenceObj_Release(obj);
        }
    }
    return rc;
}

#if USE_MATE_CACHE
static 
rc_t Cache_Init(SCursCache* c)
{
    if( c != NULL ) {
        memset(c, 0, sizeof(*c));
        return KVectorMake(&c->cache);
    }
    return 0;
}

static 
void Cache_Close(const char* name, SCursCache* c)
{
    if( c != NULL ) {
        KVectorRelease(c->cache);
        if( c->added > 0 ) {
            SAM_DUMP_DBG(2, ("%s cache stats: projected %lu added of those %lu; "
                             "hits %lu of those broken %lu;\n",
                             name, c->projected, c->added, c->hit, c->bad));
        }
    }
    memset(c, 0, sizeof(*c));
}

static
rc_t Cache_Add(uint64_t key, const SCurs* curs, const SCol* cols)
{
    /* compact values for mate record to cache as:
        pos_delta - 32bit, ref_proj - 11bit, flags - 11bit, rnext_idx - 10bit = 64bit
    */
    rc_t rc = 0;
    const ReferenceObj* r = NULL;
    uint32_t rid = 0;
    uint64_t val = 0;
    int64_t mate_id = cols[alg_MATE_ALIGN_ID].len > 0 ? cols[alg_MATE_ALIGN_ID].base.i64[0] : 0;

    if( (rc = ReferenceList_Find(param.ref_list, &r, cols[alg_REF_NAME].base.str, cols[alg_REF_NAME].len)) == 0 ) {
        rc = ReferenceObj_Idx(r, &rid);
    }
#if _DEBUGGING
    {{
        const char* rname = NULL;
        curs->cache->projected++;
        ReferenceObj_Name(r, &rname);
        SAM_DUMP_DBG(10, ("to cache row %li for mate %li: %u,%s[%hu],%u,%u,%i",
            key, mate_id, cols[alg_SAM_FLAGS].base.u32[0], rname, rid,
            cols[alg_REF_POS].len ? cols[alg_REF_POS].base.coord0[0] : 0,
            cols[alg_MATE_REF_POS].len ? cols[alg_MATE_REF_POS].base.coord0[0] : 0,
            cols[alg_TEMPLATE_LEN].len > 0 ? cols[alg_TEMPLATE_LEN].base.i32[0] : 0));
    }}
#endif
    if( rc == 0 && !(rid & 0xFC00) ) {
        int64_t pos_delta64;
        int32_t pos_delta32;

        if( mate_id != 0 ) {
            const ReferenceObj* rm = NULL;
            uint32_t rm_id;

            if( (rc = ReferenceList_Find(param.ref_list, &rm, cols[alg_MATE_REF_NAME].base.str, cols[alg_MATE_REF_NAME].len)) == 0 ) {
                rc = ReferenceObj_Idx(rm, &rm_id);
            }
            assert(rm != NULL);
            if( rc == 0 && rid != rm_id ) {
                const char* rm_name = NULL;
                ReferenceObj_Name(rm, &rm_name);
                mate_id = 0;
                SAM_DUMP_DBG(10, (" mate ref differ: %s[%hu]!", rm_name, rm_id));
            }
            ReferenceObj_Release(rm);
        }
        if( mate_id != 0 ) {
            pos_delta64 = cols[alg_MATE_REF_POS].base.coord0[0] - cols[alg_REF_POS].base.coord0[0];
        } else {
            pos_delta64 = cols[alg_REF_POS].base.coord0[0];
        }
        pos_delta32 = pos_delta64;
        if( pos_delta64 == pos_delta32 ) {
            int64_t ref_proj;
            if( mate_id == 0 ) {
                ref_proj = 0; /* indicates full value in delta */
            } else if( cols[alg_TEMPLATE_LEN].base.i32[0] < 0 ) {
                assert(pos_delta32 <= 0);
                ref_proj = pos_delta32 - cols[alg_TEMPLATE_LEN].base.i32[0];
            } else {
                assert(pos_delta32 >= 0);
                ref_proj = cols[alg_TEMPLATE_LEN].base.i32[0] - pos_delta32;
            }
            if( !(ref_proj & 0xFFFFF800) ) {
                val = (pos_delta64 << 32) | (ref_proj << 21) | (cols[alg_SAM_FLAGS].base.u32[0] << 10) | rid;
                rc = KVectorSetU64(curs->cache->cache, key, val);
            }
        }
    }
    ReferenceObj_Release(r);

#if _DEBUGGING
    if( val == 0 ) {
        SAM_DUMP_DBG(10, (" --> out of range\n"));
    } else {
        SAM_DUMP_DBG(10, (" --> %016lX\n", val));
        curs->cache->added++;
    }
#endif
    return rc;
}

static
rc_t Cache_Get(const SCurs* curs, uint64_t key, uint64_t* val)
{
    rc_t rc;
    if( (rc = KVectorGetU64(curs->cache->cache, key, val)) == 0 ) {
        uint32_t id = (*val & 0x2FF);
#if _DEBUGGING
        curs->cache->hit++;
#endif
        KVectorUnset(curs->cache->cache, key);
        if( (rc = ReferenceList_Get(param.ref_list, &curs->cache->ref, id)) != 0 ) {
            *val = 0;
            curs->cache->ref = NULL;
            rc = RC(rcExe, rcNoTarg, rcSearching, rcItem, rcNotFound);
#if _DEBUGGING
            curs->cache->bad++;
#endif
        } else {
            SAM_DUMP_DBG(10, ("from cache row %li %016lX", key, *val));
        }
    }
    return rc;
}

static
void Cache_Unpack(uint64_t val, int64_t mate_id, const SCurs* curs, SCol* cols)
{
    int32_t pos_delta = (val & 0xFFFFFFFF00000000) >> 32;
    uint32_t ref_proj = (val & 0x00000000FFE00000) >> 21;
    uint32_t flags = (val & 0x00000000001FFC00) >> 10;

    if( mate_id != 0 ) {
        /* adjust flags for mate record */
        curs->cache->sam_flags = (flags & 0x1) |
                                 (flags & ((flags & 0xC) ? 0x0 : 0x2)) |
                                 ((flags & 0x8) >> 1) |
                                 ((flags & 0x4) << 1) |
                                 ((flags & 0x20) >> 1) |
                                 ((flags & 0x10) << 1) |
                                 ((flags & 0x40) ? 0x80 : 0x40) |
                                 (flags & 0x700);
    } else {
        /* preserve flags as if original records is restored */
        curs->cache->sam_flags = flags;
    }
    cols[alg_SAM_FLAGS].base.u32 = &curs->cache->sam_flags;
    cols[alg_SAM_FLAGS].len = sizeof(curs->cache->sam_flags);

    if( param.use_seqid ) {
        ReferenceObj_SeqId(curs->cache->ref, &cols[alg_MATE_REF_NAME].base.str);
    } else {
        ReferenceObj_Name(curs->cache->ref, &cols[alg_MATE_REF_NAME].base.str);
    }
    cols[alg_MATE_REF_NAME].len = strlen(cols[alg_MATE_REF_NAME].base.str);
    if( ref_proj == 0 ) {
        curs->cache->pnext = pos_delta;
        curs->cache->tlen = 0;
    } else if( pos_delta > 0 ) {
        curs->cache->pnext = (cols[alg_REF_POS].len > 0 ? cols[alg_REF_POS].base.coord0[0] : 0) - pos_delta;
        curs->cache->tlen = - (ref_proj + pos_delta);
    } else {
        curs->cache->pnext = (cols[alg_REF_POS].len > 0 ? cols[alg_REF_POS].base.coord0[0] : 0) - pos_delta;
        curs->cache->tlen = ref_proj - pos_delta;
    }
    cols[alg_MATE_REF_POS].base.coord0 = &curs->cache->pnext;
    cols[alg_MATE_REF_POS].len = sizeof(curs->cache->pnext);
    cols[alg_TEMPLATE_LEN].base.i32 = &curs->cache->tlen;
    cols[alg_TEMPLATE_LEN].len = sizeof(curs->cache->tlen);
    {{
        uint32_t id;
        ReferenceObj_Idx(curs->cache->ref, &id);
        SAM_DUMP_DBG(10, (" --> mate %li: %u,%s[%hu],%u,%i\n",
            mate_id, curs->cache->sam_flags, cols[alg_MATE_REF_NAME].base.str,
            id, curs->cache->pnext, curs->cache->tlen));
    }}
}
#endif /* USE_MATE_CACHE */

static
rc_t OpenVTable(const VDatabase* db, STable* tbl, const char *name, bool optional)
{
    rc_t rc = VDatabaseOpenTableRead(db, &tbl->vtbl, name);
    if( GetRCState(rc) == rcNotFound && optional ) {
        rc = 0;
        tbl->vtbl = NULL;
    }
    tbl->name = name;
    return rc;
}

static
rc_t Cursor_Open(const STable* tbl, SCurs* curs, SCol* cols, SCursCache* cache)
{
    rc_t rc = 0;

    curs->vcurs = NULL;
    if(tbl == NULL || tbl->vtbl == NULL) {
        return 0;
    }
    rc = VTableCreateCachedCursorRead(tbl->vtbl, &curs->vcurs, CURSOR_CACHE);
    while(rc == 0 && cols->name != NULL) {
        if( (rc = VCursorAddColumn(curs->vcurs, &cols->idx, cols->name)) != 0 ) {
            if( GetRCObject(rc) == rcColumn &&
                (GetRCState(rc) == rcExists || (GetRCState(rc) == rcNotFound && cols->optional)) ) {
                rc = 0;
            } else {
                PLOGERR(klogErr, (klogErr, rc, "table $(t) column $(c)", "t=%s,c=%s", tbl->name, cols->name));
            }
        }
        cols++;
    }
    if( rc == 0 ) {
        if( (rc = VCursorOpen(curs->vcurs)) == 0 ) {
#if USE_MATE_CACHE
            if( cache != NULL ) {
                curs->cache = cache;
            } else {
                curs->cache = &curs->cache_local;
                rc = Cache_Init(&curs->cache_local);
            }
#endif /* USE_MATE_CACHE */
            curs->tbl = tbl;
        }
    }
    if( rc != 0 ) {
        VCursorRelease(curs->vcurs);
        curs->vcurs = NULL;
        if( rc != KLogLastErrorCode() ) {
            PLOGERR(klogErr, (klogErr, rc, "table $(t)", "t=%s", tbl->name));
        }
    } else {
        SAM_DUMP_DBG(2, ("%s: table %s\n", __func__, curs->tbl->name));
    }
    return rc;
}

static 
void Cursor_Close(SCurs* curs)
{
    if( curs != NULL && curs->vcurs != NULL ) {
        SAM_DUMP_DBG(2, ("%s: table %s, columns rows read %lu\n", __func__, curs->tbl->name, curs->col_reads_qty));
        VCursorRelease(curs->vcurs);
#if USE_MATE_CACHE
        if( curs->cache == &curs->cache_local ) {
            Cache_Close(curs->tbl->name, curs->cache);
        }
#endif /* USE_MATE_CACHE */
        memset(curs, 0, sizeof(*curs));
    }
}

static
rc_t Cursor_Read(const SCurs* curs, int64_t row_id, SCol* cols)
{
    rc_t rc = 0;
    SCol* c = NULL;

    if( (rc = VCursorCloseRow(curs->vcurs)) == 0 &&
        (rc = VCursorSetRowId(curs->vcurs, row_id)) == 0 &&
        (rc = VCursorOpenRow(curs->vcurs)) == 0 ) {
        c = cols;
        while( rc == 0 && c->name != NULL ) {
            if( c->idx == 0 ||
                (rc = VCursorCellData(curs->vcurs, c->idx, NULL, &c->base.v, NULL, &c->len)) == 0 ) {
#if _DEBUGGING
                ((SCurs*)curs)->col_reads_qty += (c->idx != 0 ? 1 : 0);
#endif
                c++;
            }
        }
    }
    if( rc != 0 ) {
        PLOGERR(klogWarn, (klogWarn, rc, "reading $(t) row $(r), column $(c)", "t=%s,r=%li,c=%s",
            curs->tbl->name, row_id, c ? c->name : "<none>"));
    }
    return rc;
}

struct {
    KWrtWriter writer;
    void* data;
    KFile* kfile;
    uint64_t pos;
} g_out_writer = {NULL};

static
rc_t CC BufferedWriter(void* self, const char* buffer, size_t bufsize, size_t* num_writ)
{
    rc_t rc = 0;

    assert(buffer != NULL);
    assert(num_writ != NULL);

    do {
        if( (rc = KFileWrite(g_out_writer.kfile, g_out_writer.pos, buffer, bufsize, num_writ)) == 0 ) {
            buffer += *num_writ;
            bufsize -= *num_writ;
            g_out_writer.pos += *num_writ;
        }
    } while(rc == 0 && bufsize > 0);
    return rc;
}

static
rc_t BufferedWriterMake(bool gzip, bool bzip2)
{
    rc_t rc = 0;

    if( gzip && bzip2 ) {
        rc = RC(rcApp, rcFile, rcConstructing, rcParam, rcAmbiguous);
    } else if( g_out_writer.writer != NULL ) {
        rc = RC(rcApp, rcFile, rcConstructing, rcParam, rcAmbiguous);
    } else if( (rc = KFileMakeStdOut(&g_out_writer.kfile)) == 0 ) {
        g_out_writer.pos = 0;
        if( gzip ) {
            KFile* gz;
            if( (rc = KFileMakeGzipForWrite(&gz, g_out_writer.kfile)) == 0 ) {
                KFileRelease(g_out_writer.kfile);
                g_out_writer.kfile = gz;
            }
        } else if( bzip2 ) {
            KFile* bz;
            if( (rc = KFileMakeBzip2ForWrite(&bz, g_out_writer.kfile)) == 0 ) {
                KFileRelease(g_out_writer.kfile);
                g_out_writer.kfile = bz;
            }
        }
        if( rc == 0 ) {
            KFile* buf;
            if( (rc = KBufFileMakeWrite(&buf, g_out_writer.kfile, false, 128 * 1024)) == 0 ) {
                KFileRelease(g_out_writer.kfile);
                g_out_writer.kfile = buf;
                g_out_writer.writer = KOutWriterGet();
                g_out_writer.data = KOutDataGet();
                rc = KOutHandlerSet(BufferedWriter, &g_out_writer);
            }
        }
    }
    return rc;
}

static
void BufferedWriterRelease( bool flush )
{
    if( flush ) {
        /* avoid flushing buffered data after failure */
        KFileRelease(g_out_writer.kfile);
    }
    if( g_out_writer.writer != NULL ) {
        KOutHandlerSet(g_out_writer.writer, g_out_writer.data);
    }
    g_out_writer.writer = NULL;
}

typedef struct ReadGroup {
    BSTNode node;
    char name[1024];
} ReadGroup;

static
int CC ReadGroup_sort( const BSTNode *item, const BSTNode *node )
{
    return strcmp(((const ReadGroup*)item)->name, ((const ReadGroup*)node)->name);
}

static
void CC ReadGroup_dump( BSTNode *n, void *data )
{
    const ReadGroup* g = (ReadGroup*)n;
    if( g->name[0] != '\0' && strcasecmp(g->name, "default") ) {
        OUTMSG(("@RG\tID:%s\n", g->name));
    }
}

static
rc_t CC DumpReadGroupsScan(const STable* tbl)
{
    rc_t rc = 0;
    SCurs curs;
    BSTree tree;

    SCol cols[] = {
        {"SPOT_GROUP", 0, {NULL}, 0, false},
        {NULL, 0, {NULL}, 0, false}
    };

    BSTreeInit(&tree);
    memset(&curs, 0, sizeof(curs));
    if( (rc = Cursor_Open(tbl, &curs, cols, NULL)) == 0 ) {
        int64_t start;
        uint64_t count;

        if( (rc = VCursorIdRange(curs.vcurs, 0, &start, &count)) == 0 ) {
            ReadGroup* node = NULL;
            uint32_t last_len = 0;

            while( count > 0 && rc == 0 ) {
                if( (rc = Cursor_Read(&curs, start, cols)) == 0 && cols[0].len != 0 ) {
                    if( node == NULL || last_len != cols[0].len || strncmp(cols[0].base.str, node->name, cols[0].len) != 0 ) {
                        node = malloc(sizeof(*node));
                        if( node == NULL ) {
                            rc = RC(rcExe, rcNode, rcConstructing, rcMemory, rcExhausted);
                        } else if( cols[0].len > sizeof(node->name) ) {
                            rc = RC(rcExe, rcString, rcCopying, rcBuffer, rcInsufficient);
                        } else {
                            last_len = cols[0].len;
                            strncpy(node->name, cols[0].base.str, last_len);
                            node->name[last_len] = '\0';
                            rc = BSTreeInsertUnique(&tree, &node->node, NULL, ReadGroup_sort);
                            if (GetRCState(rc) == rcExists) {
                                free(node);
                                rc = 0;
                            }
                        }
                    }
                } else if( GetRCState(rc) == rcNotFound && GetRCObject(rc) == rcRow ) {
                    rc = 0;
                }
                start++;
                count--;
            }
        }
        Cursor_Close(&curs);
    }
    if( rc == 0 ) {
        BSTreeForEach(&tree, false, ReadGroup_dump, NULL);
    } else if( rc != KLogLastErrorCode() ) {
        PLOGERR(klogErr, (klogErr, rc, "$(f)", "f=%s", __func__));
    }
    BSTreeWhack(&tree, NULL, NULL);
    return rc;
}

rc_t CC DumpReadGroups(const STable* tbl)
{
    rc_t rc = 0;
    const KMetadata* m;

    /* try getting list from stats meta */
    if( (rc = VTableOpenMetadataRead(tbl->vtbl, &m)) == 0 ) {
        const KMDataNode* n;
        if( (rc = KMetadataOpenNodeRead(m, &n, "/STATS/SPOT_GROUP")) == 0 ) {
            KNamelist* names;
            if( (rc = KMDataNodeListChild(n, &names)) == 0 ) {
                uint32_t i, q;
                if( (rc = KNamelistCount(names, &q)) == 0 ) {
                    for(i = 0; rc == 0 && i < q; i++) {
                        const char* nm;
                        if( (rc = KNamelistGet(names, i, &nm)) == 0 ) {
                            /* hack so printed using same func */
                            ReadGroup_dump((BSTNode*)&nm[-sizeof(BSTNode)], NULL);
                        }
                    }
                }
                KNamelistRelease(names);
            }
            KMDataNodeRelease(n);
        }
        KMetadataRelease(m);
    }
    if( GetRCState(rc) == rcNotFound ) {
        rc = DumpReadGroupsScan(tbl);
    } else if( rc != 0 && rc != KLogLastErrorCode() ) {
        PLOGERR(klogErr, (klogErr, rc, "$(f)", "f=%s", __func__));
    }
    return rc;
}

enum eseq_col {
    seq_READ = 0,
    seq_QUALITY,
    seq_SPOT_GROUP,
    seq_READ_START,
    seq_READ_LEN,
    seq_READ_TYPE,
    seq_READ_FILTER,
    seq_NAME,
    seq_PRIMARY_ALIGNMENT_ID
};

SCol g_cseq_col[] = {
    {"READ", 0, {NULL}, 0, false},
    {"(INSDC:quality:text:phred_33)QUALITY", 0, {NULL}, 0, false},
    {"SPOT_GROUP", 0, {NULL}, 0, true},
    {"READ_START", 0, {NULL}, 0, true},
    {"READ_LEN", 0, {NULL}, 0, true},
    {"READ_TYPE", 0, {NULL}, 0, true},
    {"READ_FILTER", 0, {NULL}, 0, true},
    {"NAME", 0, {NULL}, 0, true},
    /* must be last in list to avoid reading all columns */
    {"PRIMARY_ALIGNMENT_ID", 0, {NULL}, 0, true},
    {NULL, 0, {NULL}, 0, false}
};
SCurs g_cseq;

static
rc_t SeqCursorOpen(const STable* tbl)
{
    memset(&g_cseq, 0, sizeof(g_cseq));
    return param.unaligned ? Cursor_Open(tbl, &g_cseq, g_cseq_col, NULL) : 0;
}

static
void SeqCursorClose(void)
{
    Cursor_Close(&g_cseq);
}

static
rc_t Cursor_ReadAlign(const SCurs* curs, int64_t row_id, SCol* cols, uint32_t idx)
{
    rc_t rc = 0;
    SCol* c = NULL;
    SCol* mate_id = NULL;
#if USE_MATE_CACHE
    uint64_t cache_val = 0;
    bool cache_miss = false;
#endif /* USE_MATE_CACHE */

    if( (rc = VCursorCloseRow(curs->vcurs)) == 0 &&
        (rc = VCursorSetRowId(curs->vcurs, row_id)) == 0 &&
        (rc = VCursorOpenRow(curs->vcurs)) == 0 ) {
        for(; rc == 0 && cols[idx].name != NULL; idx++) {
            c = &cols[idx];
            if( c->idx != 0 ) {
#if USE_MATE_CACHE
                if( mate_id != NULL && curs->cache != NULL && !cache_miss &&
                    (idx == alg_SAM_FLAGS || idx == alg_MATE_REF_NAME || idx == alg_MATE_REF_POS || idx == alg_TEMPLATE_LEN) &&
                    mate_id->idx && mate_id->len > 0 && mate_id->base.i64[0] > 0 ) {
                    if( cache_val != 0 ) {
                        continue;
                    }
                    if( (rc = Cache_Get(curs, mate_id->base.u64[0], &cache_val)) == 0 ) {
                        continue;
                    } else if( !(GetRCObject(rc) == rcItem && GetRCState(rc) == rcNotFound) ) {
                        break;
                    } else {
                        /* avoid multiple lookups in cache */
                        cache_miss = true;
                    }
                }
#endif /* USE_MATE_CACHE */
                if( (rc = VCursorCellData(curs->vcurs, c->idx, NULL, &c->base.v, NULL, &c->len)) == 0 ) {
                    if( idx == alg_MATE_ALIGN_ID ) {
                        mate_id = &cols[alg_MATE_ALIGN_ID];
                    }
#if _DEBUGGING
                    ((SCurs*)curs)->col_reads_qty++;
#endif
                }
            }
        }
    }
    if( rc != 0 ) {
        PLOGERR(klogWarn, (klogWarn, rc, "reading $(t) row $(r), column $(c)", "t=%s,r=%li,c=%s",
            curs->tbl->name, row_id, c ? c->name : "<none>"));
#if USE_MATE_CACHE
    } else if( curs->cache == NULL ) {
    } else if( cache_val == 0 ) {
        /* this row is not from cache */
        int64_t mate_align_id = (mate_id != NULL && mate_id->len > 0) ? mate_id->base.i64[0] : 0;
        if( cols[0].idx != 0 && /* we have full cursor which means we are reading alignment table */
            /* no mate -> unaligned (not secondary!) */
            ( (mate_align_id == 0 && param.unaligned && curs->cache != &curs->cache_local) ||
            /* 2nd mate with higher rowid and more than set gap rows away */
              (mate_align_id != 0 && mate_align_id > row_id && (mate_align_id - row_id) > param.mate_row_gap_cachable) ) ) {

          rc = Cache_Add(row_id, curs, cols);
        }
    } else {
        Cache_Unpack(cache_val, row_id, curs, cols);
#endif /* USE_MATE_CACHE */
    }
    return rc;
}

static
void DumpName(const char* name, size_t name_len,
              const char spot_group_sep, const char* spot_group, size_t spot_group_len)
{
    size_t nm;
    if( param.name_prefix != NULL ) {
        OUTMSG(("%s.", param.name_prefix));
    }
    BufferedWriter(NULL, name, name_len, &nm);
    if( param.spot_group_in_name && spot_group_len > 0 ) {
        BufferedWriter(NULL, &spot_group_sep, 1, &nm);
        BufferedWriter(NULL, spot_group, spot_group_len, &nm);
    }
}

static
void DumpUnalignedFastX(const SCol cols[], uint32_t read_id, INSDC_coord_zero readStart, INSDC_coord_len readLen)
{
    size_t nm;
    /* fast[AQ] represnted in SAM fields:
       [@|>]QNAME unaligned
       SEQ
       +
       QUAL
    */
    BufferedWriter(NULL, param.fastq ? "@" : ">", 1, &nm);
    /* QNAME: [PFX.]SEQUENCE:NAME[#SPOT_GROUP] */
    DumpName(cols[seq_NAME].base.str, cols[seq_NAME].len, '#', cols[seq_SPOT_GROUP].base.str, cols[seq_SPOT_GROUP].len);
    if( read_id > 0 ) {
        OUTMSG(("/%u", read_id));
    }
    BufferedWriter(NULL, " unaligned\n", 11, &nm);
    /* SEQ: SEQUENCE.READ */
    BufferedWriter(NULL, &cols[seq_READ].base.str[readStart], readLen, &nm);
    if( param.fastq ) {
        /* QUAL: SEQUENCE.QUALITY */
        BufferedWriter(NULL, "\n+\n", 3, &nm);
        BufferedWriter(NULL, &cols[seq_QUALITY].base.str[readStart], readLen, &nm);
    }
    BufferedWriter(NULL, "\n", 1, &nm);
}

static
void DumpAlignedFastX(const SCol cols[], uint32_t read_id, bool primary)
{
    size_t nm;

    /* fast[AQ] represnted in SAM fields:
       [@|>]QNAME primary|secondary ref=RNAME pos=POS mapq=MAPQ
       SEQ
       +
       QUAL
    */
    BufferedWriter(NULL, param.fastq ? "@" : ">", 1, &nm);
    /* QNAME: [PFX.]SEQ_NAME[#SPOT_GROUP] */
    nm = cols[alg_SPOT_GROUP].len ? alg_SPOT_GROUP : alg_SEQ_SPOT_GROUP;
    DumpName(cols[alg_SEQ_NAME].base.str, cols[alg_SEQ_NAME].len, '#', cols[nm].base.str, cols[nm].len);
    if( read_id > 0 ) {
        OUTMSG(("/%u", read_id));
    }
    if( primary ) {
        BufferedWriter(NULL, " primary", 8, &nm);
    } else {
        BufferedWriter(NULL, " secondary", 10, &nm);
    }
    /* RNAME: REF_NAME or REF_SEQ_ID */
    BufferedWriter(NULL, " ref=", 5, &nm);
    if( param.use_seqid ) {
        BufferedWriter(NULL, cols[alg_REF_SEQ_ID].base.str, cols[alg_REF_SEQ_ID].len, &nm);
    } else {
        BufferedWriter(NULL, cols[alg_REF_NAME].base.str, cols[alg_REF_NAME].len, &nm);
    }
    /* POS: REF_POS, MAPQ: MAPQ */
    OUTMSG((" pos=%u mapq=%i\n", cols[alg_REF_POS].base.coord0[0] + 1, cols[alg_MAPQ].base.i32[0]));
    
    /* SEQ: READ */
    BufferedWriter(NULL, cols[alg_READ].base.str, cols[alg_READ].len, &nm);
    if( param.fastq ) {
        /* QUAL: SAM_QUALITY */
        BufferedWriter(NULL, "\n+\n", 3, &nm);
        BufferedWriter(NULL, cols[alg_SAM_QUALITY].base.str, cols[alg_SAM_QUALITY].len, &nm);
    }
    BufferedWriter(NULL, "\n", 1, &nm);
}

static
void DumpUnalignedSAM(const SCol cols[], uint32_t flags, INSDC_coord_zero readStart, INSDC_coord_len readLen,
                      const char* rnext, uint32_t rnext_len, INSDC_coord_zero pnext)
{
    unsigned i;
    size_t nm;
    /* QNAME: [PFX.]NAME[.SPOT_GROUP] */
    DumpName(cols[seq_NAME].base.str, cols[seq_NAME].len, '.', cols[seq_SPOT_GROUP].base.str, cols[seq_SPOT_GROUP].len);

    /* all these fields are const text for now */
    OUTMSG(("\t%u\t*\t0\t0\t*\t%.*s\t%u\t0\t",
            flags, rnext_len ? rnext_len : 1, rnext_len ? rnext : "*", pnext));
    /* SEQ: SEQUENCE.READ */
    if(flags & 0x10) {
        for(i = 0; i < readLen; i++) {
            char base;
            DNAReverseCompliment(&cols[seq_READ].base.str[readStart + readLen - 1 - i], &base, 1);
            BufferedWriter(NULL, &base, 1, &nm);
        }
    } else {
        BufferedWriter(NULL, &cols[seq_READ].base.str[readStart], readLen, &nm);
    }
    BufferedWriter(NULL, "\t", 1, &nm);
    /* QUAL: SEQUENCE.QUALITY */
    if(flags & 0x10) {
        for(i = 0; i < readLen; i++) {
            BufferedWriter(NULL, &cols[seq_QUALITY].base.str[readStart + readLen - 1 - i], 1, &nm);
        }
    } else {
        BufferedWriter(NULL, &cols[seq_QUALITY].base.str[readStart], readLen, &nm);
    }
    /* optional fields: */
    if( cols[seq_SPOT_GROUP].len > 0 ) {
        /* read group */
        BufferedWriter(NULL, "\tRG:Z:", 6, &nm);
        BufferedWriter(NULL, cols[seq_SPOT_GROUP].base.str, cols[seq_SPOT_GROUP].len, &nm);
    }
    BufferedWriter(NULL, "\n", 1, &nm);
}

static
void DumpAlignedSAM(const SCol cols[])
{
    size_t nm;
    uint32_t flags;

    /* QNAME: [SPOT_GROUP.]SEQ_NAME */
    nm = cols[alg_SPOT_GROUP].len ? alg_SPOT_GROUP : alg_SEQ_SPOT_GROUP;
    DumpName(cols[alg_SEQ_NAME].base.str, cols[alg_SEQ_NAME].len, '.', cols[nm].base.str, cols[nm].len);
    /* FLAG: SAM_FLAGS */
    flags = cols[alg_SAM_FLAGS].base.u32[0];
    if( !param.unaligned /** not going to dump unaligned **/
        && (flags&1)       /** but we have sequenced multiple fragments **/
        && !(flags&2) ) {   /** and not all of them align **/
        /*** remove flags talking about multiple reads **/
        flags &= ~0xC9; /* turn off 0x001 0x008 0x040 0x080 */
    }
    OUTMSG(("\t%u\t", flags));
    /* RNAME: REF_NAME or REF_SEQ_ID */
    if( param.use_seqid ) {
        BufferedWriter(NULL, cols[alg_REF_SEQ_ID].base.str, cols[alg_REF_SEQ_ID].len, &nm);
        BufferedWriter(NULL, "\t", 1, &nm);
    } else {
        BufferedWriter(NULL, cols[alg_REF_NAME].base.str, cols[alg_REF_NAME].len, &nm);
        BufferedWriter(NULL, "\t", 1, &nm);
    }
    /* POS: REF_POS */
    OUTMSG(("%u\t", cols[alg_REF_POS].base.coord0[0] + 1));
    /* MAPQ: MAPQ */
    OUTMSG(("%i\t", cols[alg_MAPQ].base.i32[0]));
    /* CIGAR: CIGAR_* */
    BufferedWriter(NULL, cols[alg_CIGAR].base.str, cols[alg_CIGAR].len, &nm);
    BufferedWriter(NULL, "\t", 1, &nm);
    
    /* RNEXT: MATE_REF_NAME or '*' */
    /* PNEXT: MATE_REF_POS or 0 */
    if( cols[alg_MATE_REF_NAME].len ) {
        if( cols[alg_MATE_REF_NAME].len == cols[alg_REF_NAME].len &&
            memcmp(cols[alg_MATE_REF_NAME].base.str, cols[alg_REF_NAME].base.str, cols[alg_MATE_REF_NAME].len) == 0 ) {
            BufferedWriter(NULL, "=\t", 2, &nm);
        } else {
            BufferedWriter(NULL, cols[alg_MATE_REF_NAME].base.str, cols[alg_MATE_REF_NAME].len, &nm);
            BufferedWriter(NULL, "\t", 1, &nm);
        }
        OUTMSG(("%u\t", cols[alg_MATE_REF_POS].base.coord0[0] + 1));
    } else {
        BufferedWriter(NULL, "*\t0\t", 4, &nm);
    }
    /* TLEN: TEMPLATE_LEN */
    OUTMSG(("%i\t", cols[alg_TEMPLATE_LEN].base.i32[0]));
    /* SEQ: READ */
    BufferedWriter(NULL, cols[alg_READ].base.str, cols[alg_READ].len, &nm);
    BufferedWriter(NULL, "\t", 1, &nm);
    /* QUAL: SAM_QUALITY */
    BufferedWriter(NULL, cols[alg_SAM_QUALITY].base.str, cols[alg_SAM_QUALITY].len, &nm);

    /* optional fields: */
    if( cols[alg_SPOT_GROUP].len > 0 ) {
        /* read group */
        BufferedWriter(NULL, "\tRG:Z:", 6, &nm);
        BufferedWriter(NULL, cols[alg_SPOT_GROUP].base.str, cols[alg_SPOT_GROUP].len, &nm);
    } else if( cols[alg_SEQ_SPOT_GROUP].len > 0 ) {
        /* backward compatibility */
        BufferedWriter(NULL, "\tRG:Z:", 6, &nm);
        BufferedWriter(NULL, cols[alg_SEQ_SPOT_GROUP].base.str, cols[alg_SEQ_SPOT_GROUP].len, &nm);
    }
    if (param.cg_style && cols[alg_CG_GC].len) {
        OUTMSG(("\t%.*s", cols[alg_CG_GC].len, cols[alg_CG_GC].base.str));
        OUTMSG(("\t%.*s", cols[alg_CG_GS].len, cols[alg_CG_GS].base.str));
        OUTMSG(("\t%.*s", cols[alg_CG_GQ].len, cols[alg_CG_GQ].base.str));
    }
    /* edit distance */
    OUTMSG(("\tNM:i:%i\n", cols[alg_EDIT_DISTANCE].len ? cols[alg_EDIT_DISTANCE].base.i32[0] : 0));
}

static
rc_t DumpUnalignedReads(const SCol calg_col[], int64_t row_id, uint64_t* rcount)
{
    rc_t rc = 0;
    uint32_t i, nreads = 0;

    if( calg_col != NULL ) {
        if( calg_col[alg_SAM_FLAGS].base.u32[0] & 0x02 ) {
            /* skip all aligned rows by flag */
            return rc;
        }
        /* get primary alignments only */
        if( (rc = Cursor_Read(&g_cseq, calg_col[alg_SEQ_SPOT_ID].base.i64[0], &g_cseq_col[seq_PRIMARY_ALIGNMENT_ID])) == 0 ) {
            for(i = 0; i < g_cseq_col[seq_PRIMARY_ALIGNMENT_ID].len; i++) {
                if( g_cseq_col[seq_PRIMARY_ALIGNMENT_ID].base.i64[i] != 0 ) {
                    if( g_cseq_col[seq_PRIMARY_ALIGNMENT_ID].base.i64[i] < row_id ) {
                        /* unaligned were printed with 1st aligment */
                        return rc;
                    }
                } else {
                    nreads++;
                }
            }
            if( nreads == g_cseq_col[seq_PRIMARY_ALIGNMENT_ID].len ) {
                /* skip all aligned rows by actual data, if flag above is not set properly */
                return rc;
            }
            row_id = calg_col[alg_SEQ_SPOT_ID].base.i64[0];
        }
    }
    if( rc == 0 && (rc = Cursor_Read(&g_cseq, row_id, g_cseq_col)) == 0 ) {
        nreads = g_cseq_col[seq_READ_LEN].idx != 0 ? g_cseq_col[seq_READ_LEN].len : 1;

        for(i = 0; i < nreads; i++) {
            INSDC_coord_zero readStart;
            INSDC_coord_len readLen;

            if( g_cseq_col[seq_PRIMARY_ALIGNMENT_ID].idx != 0 && g_cseq_col[seq_PRIMARY_ALIGNMENT_ID].base.i64[i] != 0 ) {
                continue;
            }
            if( g_cseq_col[seq_READ_TYPE].idx != 0 && !(g_cseq_col[seq_READ_TYPE].base.read_type[i] & SRA_READ_TYPE_BIOLOGICAL) ) {
                continue;
            }
            readLen = g_cseq_col[seq_READ_LEN].idx ? g_cseq_col[seq_READ_LEN].base.coord_len[i] : g_cseq_col[seq_READ].len;
            if( readLen == 0 ) {
                continue;
            }
            readStart = g_cseq_col[seq_READ_START].idx ? g_cseq_col[seq_READ_START].base.coord0[i] : 0;
            if( param.fasta || param.fastq) {
                DumpUnalignedFastX(g_cseq_col, nreads > 1 ? i + 1 : 0, readStart, readLen);
            } else {
                uint32_t cflags = 0x4;
                if( param.reverse_unaligned ) {
                    if( g_cseq_col[seq_READ_TYPE].base.read_type[i] & SRA_READ_TYPE_REVERSE ) {
                        cflags |= 0x10;
                    }
                    if( g_cseq_col[seq_READ_TYPE].base.read_type[i == nreads - 1 ? 0 : (i + 1)] & SRA_READ_TYPE_REVERSE ) {
                        cflags |= 0x20;
                    }
                }
                if( g_cseq_col[seq_READ_FILTER].idx != 0 ) {
                    if( g_cseq_col[seq_READ_FILTER].base.read_filter[i] == SRA_READ_FILTER_REJECT ) {
                        cflags |= 0x200;
                    } else if( g_cseq_col[seq_READ_FILTER].base.read_filter[i] == SRA_READ_FILTER_CRITERIA ) {
                        cflags |= 0x400;
                    }
                }
                if( calg_col == NULL ) {
                    DumpUnalignedSAM(g_cseq_col, cflags |
                                    (nreads > 1 ? (0x1 | 0x8 | (i == 0 ? 0x40 : 0x00) | (i == nreads - 1 ? 0x80 : 0x00)): 0x00),
                                    readStart, readLen, NULL, 0, 0);
                } else {
                    int c = param.use_seqid ? alg_REF_SEQ_ID : alg_REF_NAME;
                    uint16_t flags = cflags | 0x1 |
                                     ((calg_col[alg_SAM_FLAGS].base.u32[0] & 0x10) << 1) |
                                     ((calg_col[alg_SAM_FLAGS].base.u32[0] & 0x40) ? 0x80 : 0x40);
                    DumpUnalignedSAM(g_cseq_col, flags, readStart, readLen,
                        calg_col[c].base.str, calg_col[c].len, calg_col[alg_REF_POS].base.coord0[0] + 1);
                }
            }
            *rcount = *rcount + 1;
        }
    }
    return rc;
}

static
bool AlignRegionFilter(const SCol* cols)
{
    if( cols[alg_REF_NAME].len != 0 || cols[alg_REF_SEQ_ID].len != 0 ) {
        uint32_t i, j, k;

        assert(cols[alg_REF_POS].len == cols[alg_REF_LEN].len);

        for(i = 0; i < param.region_qty; i++) {
            for(j = 0; j < cols[alg_REF_POS].len; j++) {
                for(k = 0; k < param.region[i].rq; k++) {
                    if( !( cols[alg_REF_POS].base.coord0[j] + cols[alg_REF_LEN].base.coord_len[j] < param.region[i].r[k].from ||
                           cols[alg_REF_POS].base.coord0[j] > param.region[i].r[k].to ) )
                    {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

static
bool AlignDistanceFilter(const SCol* cols)
{
    if( param.mp_dist_qty != 0 || param.mp_dist_unknown ) {
        if( cols[alg_TEMPLATE_LEN].len == 0 && param.mp_dist_unknown ) {
            return true;
        } else {
            uint32_t i, j;
            for(i = 0; i < param.mp_dist_qty; i++) {
                for(j = 0; j < cols[alg_TEMPLATE_LEN].len; j++) {
                    if( (cols[alg_TEMPLATE_LEN].base.i32[j] == 0 && param.mp_dist_unknown) ||
                        (param.mp_dist[i].from <= cols[alg_TEMPLATE_LEN].base.i32[j] &&
                         cols[alg_TEMPLATE_LEN].base.i32[j] <= param.mp_dist[i].to) ) {
                        return true;
                    }
                }
            }
        }
        return false;
    }
    return true;
}

typedef struct cgOp_s {
    uint16_t length;
    uint8_t type; /* 0: match, 1: insert, 2: delete */
    char code;
} cgOp;

/* gap contains the indices of the wobbles in op
 * gap[0] is between read 1 and 2; it is the 'B'
 * gap[1] is between read 2 and 3; it is an 'N' and likely 0N
 * gap[2] is between read 3 and 4; it is an 'N'
 */

static rc_t CIGAR_to_CG_Ops(cgOp op[], unsigned const maxOps, unsigned *opCnt,
                            unsigned gap[3],
                            char const cigar[], unsigned const ciglen)
{
    unsigned i;
    unsigned ops = 0;
    unsigned rbc = 0; /* read base count */
    unsigned gapno;
    int st = 0;

    *opCnt = 0;
    for (i = 0; i < ciglen; ++ops) {
        char opChar;
        int opLen;
        int n;
        int const check = sscanf(cigar + i, "%i%c%n", &opLen, &opChar, &n);
        
        if (ops + 1 >= maxOps)
            return RC(rcExe, rcData, rcReading, rcBuffer, rcInsufficient);
        assert(check == 2);
        i += n;
        
        op[ops].length = opLen;
        op[ops].code = opChar;
        switch (opChar) {
        case 'M':
        case '=':
        case 'X':
            op[ops].type = 0;
            break;
        case 'I':
        case 'S':
            op[ops].type = 1;
            break;
        case 'D':
            op[ops].type = 2;
            break;
        default:
            return RC(rcExe, rcData, rcReading, rcConstraint, rcViolated);
        }
        if (st == 0) {
            if (op[ops].type != 2)
                rbc += opLen;
            if (rbc == 5 || rbc >= 10) {
                gap[0] = ops + 1;
                st = 1;
            }
        }
    }
    if (st == 0)
        return RC(rcExe, rcData, rcReading, rcFormat, rcNotFound); /* CG pattern not found */
    if (rbc == 5) {
        /* pattern 5,10,10,10 */
        /* make it look like 10,10,10,5 pattern by reversing order of operations */
        unsigned j = gap[0];

        if (op[j].type == 1) {
            /* our inserts make the op to right shorter, since we're reversing
             * we need to undo it and make the op to left shorter */ 
            op[j + 1].length += op[j].length;
            op[j - 1].length -= op[j].length;
        }

        for (i = 0, j = ops - 1; i < j; ++i, --j) {
            cgOp temp = op[i];
            op[i] = op[j];
            op[j] = temp;
        }
        st = -1;
    }
    for (rbc = i = 0; i < ops; ++i) {
        unsigned const nrbc = (op[i].type == 2 ? 0 : op[i].length) + rbc;
        
        if (nrbc > 10) {
            if (ops + 1 >= maxOps)
                return RC(rcExe, rcData, rcReading, rcBuffer, rcInsufficient);
            memmove(&op[i + 1], &op[i], (ops - i) * sizeof(op[0])); ++ops;
            op[i + 1].length -= (op[i].length = 10 - rbc);
            rbc = 0;
        }
        else if (nrbc == 10) {
            rbc = 0;
        }
        else {
            rbc = nrbc;
        }
    }
    for (gapno = 3, rbc = i = 0; i < ops - 1; ++i) {
        if (op[i].type == 2)
            continue;
        
        rbc += op[i].length;
        if (rbc < 10) {
            continue;
        }
        assert(rbc == 10);
        rbc = 0;
        if (gapno == 0)
            return RC(rcExe, rcData, rcReading, rcConstraint, rcViolated);
        gap[--gapno] = i + 1;
        if (gapno == 0) {
            if (st == -1) {
                ++i;
                if (op[i].type != 1) {
                    if (ops + 1 >= maxOps)
                        return RC(rcExe, rcData, rcReading, rcBuffer, rcInsufficient);
                    memmove(&op[i + 1], &op[i], (ops - i) * sizeof(op[0])); ++ops;
                    op[i].length = 0;
                    op[i].type = 1;
                }
                op[i].code = 'B';
                
                /* our inserts make the op to right shorter; B's don't do this */
                op[i + 1].length += op[i].length;
            }
            else {
                --gap[0];
                if (op[i].type != 1) {
                    if (ops + 1 >= maxOps)
                        return RC(rcExe, rcData, rcReading, rcBuffer, rcInsufficient);
                    memmove(&op[i + 1], &op[i], (ops - i) * sizeof(op[0])); ++ops;
                    op[i].length = 0;
                    op[i].type = 1;
                }
                op[i].code = 'B';
                
                /* our inserts make the op to left shorter; B's don't do this */
                op[i - 1].length += op[i].length;
            }
            break;
        }
        if (op[i + 1].type != 2) {
            if (ops + 1 >= maxOps)
                return RC(rcExe, rcData, rcReading, rcBuffer, rcInsufficient);
            memmove(&op[i + 1], &op[i], (ops - i) * sizeof(op[0])); ++ops;
            ++i;
            op[i].length = 0;
            op[i].type = 2;
            op[i].code = 'N';
        }
        else {
            ++i;
            op[i].code = 'N';
        }
    }
    if (gapno != 0)
        return RC(rcExe, rcData, rcReading, rcConstraint, rcViolated);
    if (st == -1) {
        /* pattern 5,10,10,10 */
        unsigned j;
        
        for (i = 0, j = ops - 1; i < j; ++i, --j) {
            cgOp temp = op[i];
            op[i] = op[j];
            op[j] = temp;
        }
        st = -1;
        for (i = 0; i < 3; ++i) {
            gap[i] = (ops - 1) - gap[i];
        }
    }
    
    *opCnt = ops;
    return 0;
}

static rc_t GenerateGCData(SCol cols[])
{
    static char newCIGAR[35 * 11];
    rc_t rc = 0;
    
    memset(&cols[alg_CG_GC], 0, sizeof(cols[alg_CG_GC]));
    memset(&cols[alg_CG_GS], 0, sizeof(cols[alg_CG_GS]));
    memset(&cols[alg_CG_GQ], 0, sizeof(cols[alg_CG_GQ]));
    
    if (cols[alg_READ].len == 35 && cols[alg_SAM_QUALITY].len == 35) {
        unsigned gap[3];
        cgOp cigOp[35];
        unsigned opCnt;
        unsigned i;
        unsigned j;
        size_t sz;
        
        rc = CIGAR_to_CG_Ops(cigOp, sizeof(cigOp)/sizeof(cigOp[0]), &opCnt, gap, cols[alg_CIGAR].base.str, cols[alg_CIGAR].len);
        if (rc) return 0;
        
        if (param.cg_style == 1) {
            static char newSeq[35];
            static char newQual[35];
            static char GC[40];
            static char GS[40];
            static char GQ[40];
            unsigned const B_len = cigOp[gap[0]].length;
            unsigned const B_at = gap[0] < gap[2] ? 5 : 30;
            
            memcpy(newSeq, cols[alg_READ].base.v, 35);
            memcpy(newQual, cols[alg_SAM_QUALITY].base.v, 35);
            
            cols[alg_CG_GC].base.v = GC;
            cols[alg_CG_GS].base.v = GS;
            cols[alg_CG_GQ].base.v = GQ;
            
            cols[alg_READ].base.v = newSeq;
            cols[alg_READ].len = 35 - B_len;
            
            cols[alg_SAM_QUALITY].base.v = newQual;
            cols[alg_SAM_QUALITY].len = 35 - B_len;
            
            if (gap[0] < gap[2]) {
                string_printf(GC, sizeof(GC), &sz, "GC:Z:%uS%uG%uS", 5 - B_len, B_len, 30 - B_len);
                cols[alg_CG_GC].len = sz;
                
                string_printf(GS, sizeof(GS), &sz, "GS:Z:%.*s", 2 * B_len, &newSeq[5 - B_len]);
                cols[alg_CG_GS].len = sz;
                
                string_printf(GQ, sizeof(GQ), &sz, "GQ:Z:%.*s", 2 * B_len, &newQual[5 - B_len]);
                cols[alg_CG_GQ].len = sz;
                
                cigOp[gap[0] - 1].length -= B_len;
            }
            else {
                string_printf(GC, sizeof(GC), &sz, "GC:Z:%uS%uG%uS", 30 - B_len, B_len, 5 - B_len);
                cols[alg_CG_GC].len = sz;
                
                string_printf(GS, sizeof(GS), &sz, "GS:Z:%.*s", 2 * B_len, &newSeq[30 - B_len]);
                cols[alg_CG_GS].len = sz;
                
                string_printf(GQ, sizeof(GQ), &sz, "GQ:Z:%.*s", 2 * B_len, &newQual[30 - B_len]);
                cols[alg_CG_GQ].len = sz;
                
                cigOp[gap[0] + 1].length -= B_len;
            }
            for (i = B_at; i < B_at + B_len; ++i) {
                int const Lq = newQual[i - B_len];
                int const Rq = newQual[i];
                
                if (Lq < Rq) {
                    newSeq[i - B_len] = newSeq[i];
                    newQual[i - B_len] = Rq;
                }
                else {
                    newSeq[i] = newSeq[i - B_len];
                    newQual[i] = Lq;
                }
            }
            memmove(&newSeq [B_at], &newSeq [B_at + B_len], 35 - B_at - B_len);
            memmove(&newQual[B_at], &newQual[B_at + B_len], 35 - B_at - B_len);
            
            /* remove B; eg 5M2I10MxN10MyN10M -> 5M10MxN10MyN10M -> 15MxN10MyN10M */
            memmove(&cigOp[gap[0]], &cigOp[gap[0] + 1], (opCnt - (gap[0] + 1)) * sizeof(cigOp[0])); --opCnt;
            /* remove zero length ops */
            for (j = i = 0; i < opCnt; ++i, ++j) {
                if (cigOp[i].length == 0) {
                    ++j;
                    --opCnt;
                }
                cigOp[i] = cigOp[j];
            }
            /* merge adjacent ops */
            for (i = opCnt; i > 1; ) {
                --i;
                if (cigOp[i - 1].code == cigOp[i].code) {
                    cigOp[i - 1].length += cigOp[i].length;
                    memmove(&cigOp[i], &cigOp[i + 1], (opCnt - 1 - i) * sizeof(cigOp[0]));
                    --opCnt;
                }
            }
        }
        else {
            /* print CG CIGAR in its full glory */
        }

        for (i = j = 0; i < opCnt && rc == 0; ++i) {
            rc = string_printf(&newCIGAR[j], sizeof(newCIGAR) - j, &sz, "%u%c", cigOp[i].length, cigOp[i].code);
            j += sz;
        }
        cols[alg_CIGAR].base.v = newCIGAR;
        cols[alg_CIGAR].len = j;
    }
    return rc;
}

static
rc_t DumpAlignedRowList(const SCurs* curs, SCol* cols, const SCol* ids,
                        bool primary, uint64_t* rcount)
{
    rc_t rc = 0;
    uint32_t i;

    for(i = 0; rc == 0 && i < ids->len; i++) {
        if( (rc = Cursor_Read(curs, ids->base.i64[i], &cols[alg_REGION_FILTER])) == 0 ) {
            if( AlignRegionFilter(cols) && AlignDistanceFilter(cols) ) {
                if( (rc = Cursor_ReadAlign(curs, ids->base.i64[i], cols, 0)) == 0 ) {
                    *rcount = *rcount + 1;
                    if( param.fasta || param.fastq) {
                        DumpAlignedFastX(cols, cols[alg_SEQ_READ_ID].base.coord1[0], primary);
                    } else {
                        if (param.cg_style)
                            rc = GenerateGCData(cols);
                        if (rc == 0)
                            DumpAlignedSAM(cols);
                    }
                    if( rc == 0 && primary && param.unaligned ) {
                        rc = DumpUnalignedReads(cols, ids->base.i64[i], rcount);
                    }
                }
            }
        }
        rc = rc ? rc : Quitting();
    }
    return rc;
}

static
rc_t DumpAligned(const STable* talgP, const STable* talgS, const STable* tref, SCursCache* seq_cache)
{
    rc_t rc = 0;
    SCurs calgP, calgS;
    SCol calg_colP[sizeof(g_alg_col_tmpl) / sizeof(g_alg_col_tmpl[0])];
    SCol calg_colS[sizeof(g_alg_col_tmpl) / sizeof(g_alg_col_tmpl[0])];

    memset(&calgP, 0, sizeof(calgP));
    memset(&calgS, 0, sizeof(calgS));
    memcpy(calg_colP, g_alg_col_tmpl, sizeof(calg_colP));
    memcpy(calg_colS, g_alg_col_tmpl, sizeof(calg_colS));
    
    /* copy columns for secondary */
    if( (rc = Cursor_Open(talgP, &calgP, calg_colP, seq_cache)) == 0 &&
        (rc = Cursor_Open(talgS, &calgS, calg_colS, NULL)) == 0 ) {

        int64_t start = 0;
        uint64_t count = 0;

        if( param.region_qty == 0 ) {
            SAM_DUMP_DBG(2, ("%s PRIMARY_ALIGNMENTs\n", param.accession));
            if( rc == 0 && (rc = VCursorIdRange(calgP.vcurs, 0, &start, &count)) == 0 ) {
                uint64_t rcount = 0;
                if( param.test_rows != 0 && count > param.test_rows ) {
                    count = param.test_rows;
                }
                SAM_DUMP_DBG(2, ("range from %ld qty %lu\n", start, count));
                while( count > 0 && rc == 0 ) {
                    if( (rc = Cursor_Read(&calgP, start, &calg_colP[alg_DISTANCE_FILTER])) == 0 ) {
                        if( AlignDistanceFilter(calg_colP) ) {
                            if( (rc = Cursor_ReadAlign(&calgP, start, calg_colP, 0)) == 0 ) {
                                ++rcount;
                                if( param.fasta || param.fastq ) {
                                    DumpAlignedFastX(calg_colP, calg_colP[alg_SEQ_READ_ID].base.coord1[0], true);
                                } else {
                                    if (param.cg_style)
                                        rc = GenerateGCData(calg_colP);
                                    if (rc == 0)
                                        DumpAlignedSAM(calg_colP);
                                }
                            } else if( GetRCState(rc) == rcNotFound && GetRCObject(rc) == rcRow ) {
                                rc = 0;
                            }
                        }
                    } else if( GetRCState(rc) == rcNotFound && GetRCObject(rc) == rcRow ) {
                        rc = 0;
                    }
                    start++;
                    count--;
                    rc = rc ? rc : Quitting();
                }
                PLOGMSG(klogInfo, (klogInfo, "$(a): $(c) primary sequences", "a=%s,c=%lu", param.accession, rcount));
            }
            /* close early to release cache */
            Cursor_Close(&calgP);

            if( !param.only_primaries && calgS.vcurs != NULL && rc == 0 ) {
                SAM_DUMP_DBG(2, ("%s SECONDARY_ALIGNMENTs\n", param.accession));
                if( (rc = VCursorIdRange(calgS.vcurs, 0, &start, &count)) == 0 ) {
                    uint64_t rcount = 0;
                    if( param.test_rows != 0 && count > param.test_rows ) {
                        count = param.test_rows;
                    }
                    SAM_DUMP_DBG(2, ("range from %ld qty %lu\n", start, count));
                    while( count > 0 && rc == 0 ) {
                        if( (rc = Cursor_Read(&calgS, start, &calg_colS[alg_DISTANCE_FILTER])) == 0 ) {
                            if( AlignDistanceFilter(calg_colS) ) {
                                if( (rc = Cursor_ReadAlign(&calgS, start, calg_colS, 0)) == 0 ) {
                                    ++rcount;
                                    if( param.fasta || param.fastq ) {
                                        DumpAlignedFastX(calg_colS, calg_colS[alg_SEQ_READ_ID].base.coord1[0], false);
                                    } else {
                                        if (param.cg_style)
                                            rc = GenerateGCData(calg_colS);
                                        if (rc == 0)
                                            DumpAlignedSAM(calg_colS);
                                    }
                                }
                            }
                        }
                        if( GetRCState(rc) == rcNotFound && GetRCObject(rc) == rcRow ) {
                            rc = 0;
                        }
                        start++;
                        count--;
                        rc = rc ? rc : Quitting();
                    }
                    PLOGMSG(klogInfo, (klogInfo, "$(a): $(c) secondary sequences", "a=%s,c=%lu", param.accession, rcount));
                }
            }
        } else {
            /* use index to set REF_NAME ranges */
            uint32_t r;
            const KIndex* iname = NULL;
            SCurs cref;

            SCol cref_col[] = {
                {"MAX_SEQ_LEN", 0, {NULL}, 0, false},
                {"PRIMARY_ALIGNMENT_IDS", 0, {NULL}, 0, false},
                {"SECONDARY_ALIGNMENT_IDS", 0, {NULL}, 0, false},
                {NULL, 0, {NULL}, 0, false}
            };
            enum eref_col {
                ref_MAX_SEQ_LEN = 0,
                ref_PRIMARY_ALIGNMENT_IDS,
                ref_SECONDARY_ALIGNMENT_IDS
            };

            memset(&cref, 0, sizeof(cref));
            if( (rc = VTableOpenIndexRead(tref->vtbl, &iname, "i_name")) == 0 &&
                (rc = Cursor_Open(tref, &cref, cref_col, NULL)) == 0 ) {
                for(r = 0; rc == 0 && r < param.region_qty; r++ ) {
                    if( (rc = KIndexFindText(iname, param.region[r].name, &start, &count, NULL, NULL)) == 0 ) {
                        bool skip_initial = true;
                        uint64_t cur_pos = 0, rcount = 0;
                        uint32_t max_seq_len = 0;

                        SAM_DUMP_DBG(2, ("REFERENCE %s index range is [%lu:%lu]\n", param.region[r].name, start, start + count - 1));
                        while( count > 0 && rc == 0 ) {
                            if( (rc = Cursor_Read(&cref, start, cref_col)) == 0 ) {
                                if( skip_initial ) {
                                    /* scroll to row with 1st region offset - 1 so algnmts tails in the range are not lost */
                                    uint64_t inc = param.region[r].r[0].from / cref_col[ref_MAX_SEQ_LEN].base.u32[0];
                                    max_seq_len = cref_col[ref_MAX_SEQ_LEN].base.u32[0];
                                    skip_initial = false;
                                    inc = inc ? inc - 1 : 0;
                                    if( start + inc != start ) {
                                        start += inc;
                                        count -= inc;
                                        cur_pos = max_seq_len * inc;
                                        continue;
                                    }
                                } else if( cur_pos > param.region[r].max_to ) {
                                    break;
                                }
                                /*SAM_DUMP_DBG(2, ("row %s index range is [%lu:%lu] pos %lu\n",
                                    param.region[r].name, start, start + count - 1, cur_pos));*/
                                if( (rc = DumpAlignedRowList(&calgP, calg_colP, &cref_col[ref_PRIMARY_ALIGNMENT_IDS],
                                                             true, &rcount)) == 0 ) {
                                    if( !param.only_primaries && calgS.vcurs != NULL ) {
                                        rc = DumpAlignedRowList(&calgS, calg_colS,
                                                                &cref_col[ref_SECONDARY_ALIGNMENT_IDS], false, &rcount);
                                    }
                                }
                            } else if( GetRCState(rc) == rcNotFound && GetRCObject(rc) == rcRow ) {
                                rc = 0;
                            }
                            start++;
                            count--;
                            cur_pos += max_seq_len;
                            rc = rc ? rc : Quitting();
                            if( param.test_rows != 0 && rcount > param.test_rows ) {
                                break;
                            }
                        }
                        PLOGMSG(klogInfo, (klogInfo, "$(a): $(c) region sequences", "a=%s,c=%lu", param.accession, rcount));
                    } else if( GetRCState(rc) == rcNotFound ) {
                        PLOGMSG(klogWarn, (klogWarn, "REFERENCE $(r) not present in data", "r=%s", param.region[r].name));
                        rc = 0;
                    }
                }
            }
            Cursor_Close(&cref);
            KIndexRelease(iname);
            Cursor_Close(&calgP);
            VTableRelease(talgP->vtbl);
        }
    }
    Cursor_Close(&calgS);
    return rc;
}

static
rc_t DumpUnaligned(const STable* talg)
{
    rc_t rc = 0;
    int64_t start = 0;
    uint64_t count = 0;

    SCurs calg;
    SCol calg_col[sizeof(g_alg_col_tmpl) / sizeof(g_alg_col_tmpl[0])];

    if( talg != NULL ) {
        memcpy(calg_col, g_alg_col_tmpl, sizeof(calg_col));
        /* cut off columns which will not bu used */
        calg_col[alg_REF_POS + 1].name = NULL;
    }
    memset(&calg, 0, sizeof(calg));
    if( (rc = VCursorIdRange(g_cseq.vcurs, 0, &start, &count)) == 0 ) {
        uint64_t rcount = 0;
        if( param.test_rows != 0 && count > param.test_rows ) {
            count = param.test_rows;
        }
        SAM_DUMP_DBG(2, ("%s SEQUENCE table range from %ld qty %lu\n", param.accession, start, count));
        while( count > 0 && rc == 0 ) {
            uint32_t i;

            if( g_cseq_col[seq_PRIMARY_ALIGNMENT_ID].idx == 0 ) {
                rc = DumpUnalignedReads(NULL, start, &rcount);
            } else {
                /* avoid reading whole sequence cursor data unnecessarily */
                if( (rc = Cursor_Read(&g_cseq, start, &g_cseq_col[seq_PRIMARY_ALIGNMENT_ID])) == 0 ) {
                    /* find if its completely unaligned */
                    int64_t min_prim_id = 0;
                    bool has_unaligned = false;
                    for(i = 0; rc == 0 && i < g_cseq_col[seq_PRIMARY_ALIGNMENT_ID].len; i++) {
                        int64_t x = g_cseq_col[seq_PRIMARY_ALIGNMENT_ID].base.i64[i];
                        has_unaligned |= x == 0;
                        if( (min_prim_id == 0 && x != 0) || min_prim_id < x ) {
                            min_prim_id = x;
                        }
                    }
                    if( min_prim_id == 0 ) {
                        rc = DumpUnalignedReads(NULL, start, &rcount);
                    } else if( has_unaligned ) {
                        if( calg.vcurs == NULL ) {
                            rc = Cursor_Open(talg, &calg, &calg_col[alg_MATE_ALIGN_ID], g_cseq.cache);
                        }
                        if( rc == 0 ) {
#if USE_MATE_CACHE
                            uint64_t val;
                            if( (rc = Cache_Get(&calg, min_prim_id, &val)) == 0 ) {
                                calg_col[alg_REF_POS].len = 0;
                                Cache_Unpack(val, 0, &calg, calg_col);
                                calg_col[alg_SEQ_SPOT_ID].base.i64 = &start;
                                calg_col[alg_SEQ_SPOT_ID].len = 1;
                                memcpy(&calg_col[alg_REF_NAME], &calg_col[alg_MATE_REF_NAME], sizeof(SCol));
                                memcpy(&calg_col[alg_REF_SEQ_ID], &calg_col[alg_MATE_REF_NAME], sizeof(SCol));
                                memcpy(&calg_col[alg_REF_POS], &calg_col[alg_MATE_REF_POS], sizeof(SCol));
                            } else if( !(GetRCState(rc) == rcNotFound && GetRCObject(rc) == rcItem) ) {
                                break;
                            } else {
#endif /* USE_MATE_CACHE */
                                rc = Cursor_ReadAlign(&calg, min_prim_id, calg_col, alg_MATE_ALIGN_ID);
#if USE_MATE_CACHE
                            }
#endif /* USE_MATE_CACHE */
                            rc = rc ? rc : DumpUnalignedReads(calg_col, min_prim_id, &rcount);
                        }
                    }
                } else if( GetRCState(rc) == rcNotFound && GetRCObject(rc) == rcRow ) {
                    rc = 0;
                }
            }
            start++;
            count--;
            rc = rc ? rc : Quitting();
        }
        PLOGMSG(klogInfo, (klogInfo, "$(a): $(c) unaligned sequences", "a=%s,c=%lu", param.accession, rcount));
    }
    Cursor_Close(&calg);
    return rc;
}

static
rc_t Dump(uint32_t idx, bool multi_run, bool* error_reported)
{
    rc_t rc = 0;
    const VDBManager* mgr = NULL;
    uint32_t i;

    assert(error_reported);

    if( (rc = VDBManagerMakeRead(&mgr, NULL)) == 0 ) {

        const VDatabase* db;

        ReportSetVDBManager(mgr);

        if( (rc = VDBManagerOpenDBRead(mgr, &db, NULL, param.path)) == 0 ) {
            STable seq, algP, algS, ref;

            memset(&seq, 0, sizeof(seq));
            memset(&algP, 0, sizeof(algP));
            memset(&algS, 0, sizeof(algS));
            memset(&ref, 0, sizeof(ref));

            ReportResetDatabase(param.path, db);

            if( (rc = OpenVTable(db, &ref, "REFERENCE", true)) == 0 &&
                (rc = OpenVTable(db, &seq, "SEQUENCE", false)) == 0 &&
                (rc = OpenVTable(db, &algP, "PRIMARY_ALIGNMENT", true)) == 0 &&
                ((!param.only_primaries && (rc = OpenVTable(db, &algS, "SECONDARY_ALIGNMENT", true)) == 0) ||
                 param.only_primaries) ) {

                if( !param.noheader && !param.reheader && !multi_run ) {
                    /* grab header from db meta node */
                    const KMetadata* m;
                    if( (rc = VDatabaseOpenMetadataRead(db, &m)) == 0 ) {
                        const KMDataNode* n;
                        if( (rc = KMetadataOpenNodeRead(m, &n, "BAM_HEADER")) == 0 ) {
                            size_t offset = 0, num_read, remaining = ~0;
                            char buffer[40960];
                            while(rc == 0 && remaining > 0 ) {
                                if( (rc = KMDataNodeRead(n, offset, buffer, sizeof(buffer),
                                                         &num_read, &remaining)) == 0 ) {
                                    OUTMSG(("%.*s", ( uint32_t ) num_read, buffer));
                                    offset += num_read;
                                }
                            }
                            if( rc == 0 && buffer[num_read - 1] != '\n' ) {
                                OUTMSG(("\n"));
                            }
                            KMDataNodeRelease(n);
                        } else if( GetRCState(rc) == rcNotFound ) {
                            param.reheader = true;
                            rc = 0;
                        }
                        KMetadataRelease(m);
                    }
                }
                if( ref.vtbl != NULL ) {
                    rc = ReferenceList_MakeTable(&param.ref_list, ref.vtbl, 0, CURSOR_CACHE, NULL, 0);
                } else if( algP.vtbl == NULL ) {
                    /* all is unaligned - operate as on a single SEQUENCE table below */
                    param.unaligned = true;
                }
                if( rc == 0 && !param.noheader && (param.reheader || multi_run) ) {
                    if( !multi_run || idx == 0 ) {
                        OUTMSG(("@HD\tVN:1.3\n"));
                    }
                    if( !multi_run ) {
                        if( param.ref_list != NULL ) {
                            rc = RefSeqPrint();
                        }
                        if( rc == 0 ) {
                            rc = DumpReadGroups(&seq);
                        }
                    }
                }
                if( rc == 0 && !param.noheader && idx == 0) {
                    for(i = 0; rc == 0 && i < param.comments_qty; i++) {
                        OUTMSG(("@CO\t%s\n", param.comments[i]));
                    }
                }
                if( rc == 0 && (rc = SeqCursorOpen(&seq)) == 0 ) {
                    if( algP.vtbl != NULL ) {
                        rc = DumpAligned(&algP, &algS, &ref, g_cseq.cache);
                    }
                    if( rc == 0 && param.unaligned && param.region_qty == 0 ) {
                        rc = DumpUnaligned(algP.vtbl != NULL ? &algP : NULL);
                    }
                    SeqCursorClose();
                }
                ReferenceList_Release(param.ref_list);
            }
            if( UIError(rc, db, NULL) ) {
                UIDatabaseLOGError(rc, db, true);
                *error_reported = true;
            }
            VTableRelease(algS.vtbl);
            VTableRelease(algP.vtbl);
            VTableRelease(ref.vtbl);
            VTableRelease(seq.vtbl);
            VDatabaseRelease(db);
        } else {
            STable seq;
            VSchema* schema = NULL;

            seq.name = "SEQUENCE";
            param.unaligned = true;
UseLegacy:
            if( (rc = VDBManagerOpenTableRead(mgr, &seq.vtbl, schema, param.path)) == 0 ) {
                ReportResetTable(param.path, seq.vtbl);

                if( !param.noheader && idx == 0) {
                    OUTMSG(("@HD\tVN:1.3\n"));
                    rc = DumpReadGroups(&seq);
                    for(i = 0; rc == 0 && i < param.comments_qty; i++) {
                        OUTMSG(("@CO\t%s\n", param.comments[i]));
                    }
                }
                if( rc == 0 && (rc = SeqCursorOpen(&seq)) == 0 ) {
                    rc = DumpUnaligned(NULL);
                    SeqCursorClose();
                }
                VTableRelease(seq.vtbl);
            } else if (UIError(rc, NULL, NULL)) {
                UITableLOGError(rc, NULL, true);
                *error_reported = true;
            }
            
            if( rc != 0 && schema == NULL ) {
                if( (rc = VDBManagerMakeSRASchema(mgr, &schema)) == 0 ) {
                    goto UseLegacy;
                }
            }
            VSchemaRelease(schema);
        }
    }
    VDBManagerRelease(mgr);
    return rc;
}

ver_t CC KAppVersion( void )
{
    return SAM_DUMP_VERS;
}

const char* unaligned_usage[] = {"Output unaligned reads", NULL};
const char* primaryonly_usage[] = {"Output only primary alignments", NULL};
const char* cigartype_usage[] = {"Output long version of CIGAR", NULL};
const char* cigarCG_usage[] = {"Output CG version of CIGAR", NULL};
const char* header_usage[] = {"Always reconstruct header", NULL};
const char* noheader_usage[] = {"Do not output headers", NULL};
const char* region_usage[] = {"Filter by position on genome.",
                              "Name can either be file specific name (ex: \"chr1\" or \"1\").",
                              "\"from\" and \"to\" are 1-based coordinates", NULL};
const char* distance_usage[] = {"Filter by distance between matepairs.",
                                "Use \"unknown\" to find matepairs split between the references.",
                                "Use from-to to limit matepair distance on the same reference", NULL};
const char* seq_id_usage[] = {"Print reference SEQ_ID in RNAME instead of NAME", NULL};
const char* identicalbases_usage[] = {"Output '=' if base is identical to reference", NULL};
const char* gzip_usage[] = {"Compress output using gzip", NULL};
const char* bzip2_usage[] = {"Compress output using bzip2", NULL};
const char* qname_usage[] = {"Add .SPOT_GROUP to QNAME", NULL};
const char* fasta_usage[] = {"Produce Fasta formatted output", NULL};
const char* fastq_usage[] = {"Produce FastQ formatted output", NULL};
const char* prefix_usage[] = {"Prefix QNAME: prefix.QNAME", NULL};
const char* reverse_usage[] = {"Reverse unaligned reads according to read type", NULL};
const char* comment_usage[] = {"Add comment to header. Use multiple times for several lines. Use quotes", NULL};
const char* CG_usage[] = {"Apply CG fixups and output CG-specific columns", NULL};

const char* usage_params[] =
{
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    "'text'",
    "name[:from-to]",
    "from-to|unknown",
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    "prefix",
    NULL,
    NULL,
    NULL,
    NULL
};

enum eArgs {
    earg_unaligned = 0,
    earg_prim_only,
    earg_cigartype,
    earg_cigarCG,
    earg_header,
    earg_noheader,
    earg_comment,
    earg_region,
    earg_distance,
    earg_seq_id,
    earg_identicalbases,
    earg_gzip,
    earg_bzip2,
    earg_qname,
    earg_fastq,
    earg_fasta,
    earg_prefix,
    earg_reverse,
    earg_test_rows,
    earg_mate_row_gap_cachable,
    earg_cg_style
};

OptDef DumpArgs[] =
{
    {"unaligned", "u", NULL, unaligned_usage, 0, false, false},
    {"primary", "1", NULL, primaryonly_usage, 0, false, false},
    {"cigar-long", "c", NULL, cigartype_usage, 0, false, false},
    {"cigar-CG", NULL, NULL, cigarCG_usage, 0, false, false},
    {"header", "r", NULL, header_usage, 0, false, false},
    {"no-header", "n", NULL, noheader_usage, 0, false, false},
    {"header-comment", NULL, NULL, comment_usage, 0, true, false},
    {"aligned-region", NULL, NULL, region_usage, 0, true, false},
    {"matepair-distance", NULL, NULL, distance_usage, 0, true, false},
    {"seqid", "s", NULL, seq_id_usage, 0, false, false},
    {"hide-identical", "=", NULL, identicalbases_usage, 0, false, false},
    {"gzip", NULL, NULL, gzip_usage, 0, false, false},
    {"bzip2", NULL, NULL, bzip2_usage, 0, false, false},
    {"spot-group", "g", NULL, qname_usage, 0, false, false},
#ifdef NCBI
    {"fastq", NULL, NULL, fastq_usage, 0, false, false},
    {"fasta", NULL, NULL, fasta_usage, 0, false, false},
#else
    {"fastq", NULL, NULL, NULL, 0, false, false},
    {"fasta", NULL, NULL, NULL, 0, false, false},
#endif
    {"prefix", "p", NULL, prefix_usage, 0, true, false},
    {"reverse", NULL, NULL, reverse_usage, 0, false, false},
    {"test-rows", NULL, NULL, NULL, 0, true, false},
    {"mate-cache-row-gap", NULL, NULL, NULL, 0, true, false},
    {"CG", NULL, NULL, CG_usage, 0, false, false}
};

const char UsageDefaultName[] = "sam-dump";

rc_t CC UsageSummary (const char * progname)
{
    return KOutMsg ( "Usage:\n"
        "\t%s [options] path-to-run[ path-to-run ...]\n\n", progname );
}


rc_t CC Usage( const Args* args )
{
    const char * progname = UsageDefaultName;
    const char * fullpath = UsageDefaultName;
    rc_t rc;
    unsigned i;

    rc = ArgsProgram(args, &fullpath, &progname);

    UsageSummary(progname);

    OUTMSG (("Options:\n"));
    for(i = 0; i < sizeof(DumpArgs)/sizeof(DumpArgs[0]); i++ ) {
        if( DumpArgs[i].help != NULL ) {
            HelpOptionLine(DumpArgs[i].aliases, DumpArgs[i].name, usage_params[i], DumpArgs[i].help);
        }
    }
    OUTMSG (("\n"));
    HelpOptionsStandard();

    HelpVersion(fullpath, KAppVersion());

    return rc;
}

rc_t ResolvePath(const char* accession, const char** path)
{
    rc_t rc = 0;
    static char tblpath[4096];
    static SRAPath* pmgr = NULL;

    if( accession == NULL && path == NULL ) {
        SRAPathRelease(pmgr);
    } else if( pmgr != NULL || (rc = SRAPathMake(&pmgr, NULL)) == 0 ||
               (GetRCState(rc) == rcNotFound && GetRCTarget(rc) == rcDylib) ) {
        *path = tblpath;
        tblpath[0] = '\0';
        rc = 0;
        do {
            if( pmgr != NULL && !SRAPathTest(pmgr, accession) ) {
                /* try to resolve the path using mgr */
                if ( (rc = SRAPathFind(pmgr, accession, tblpath, sizeof(tblpath))) == 0 ) {
                    break;
                }
            }
            if( strlen(accession) >= sizeof(tblpath) ) {
                rc = RC(rcExe, rcPath, rcResolving, rcBuffer, rcInsufficient);
            } else {
                strcpy(tblpath, accession);
                rc = 0;
            }
        } while(false);
    }
    return rc;
}

rc_t CC KMain( int argc, char* argv[] )
{
    rc_t rc = 0;
    Args* args;
    const char* errmsg = "stop";

    bool error_reported = false;

    memset(&g_out_writer, 0, sizeof(g_out_writer));
    KOutHandlerSetStdOut();
    KStsHandlerSetStdErr();
    KLogHandlerSetStdErr();
    ( void ) KDbgHandlerSetStdErr();

    ReportBuildDate (  __DATE__ );

    if( (rc = ArgsMakeAndHandle(&args, argc, argv, 1, DumpArgs, sizeof(DumpArgs)/sizeof(DumpArgs[0]))) == 0 ) {
        uint32_t pcount, count[sizeof(DumpArgs)/sizeof(DumpArgs[0])];

        memset(&param, 0, sizeof(param));
        memset(&count, 0, sizeof(count));

        if( (rc = ArgsParamCount(args, &pcount)) != 0 || pcount < 1 ) {
            errmsg = "";
            rc = argc < 2 ? 0 : RC(rcExe, rcArgv, rcParsing, rcParam, rcInsufficient);
            MiniUsage(args);
        } else if( (rc = ArgsOptionCount(args, DumpArgs[earg_unaligned].name, &count[earg_unaligned])) != 0 ) {
            errmsg = DumpArgs[earg_unaligned].name;
        } else if( (rc = ArgsOptionCount(args, DumpArgs[earg_prim_only].name, &count[earg_prim_only])) != 0 ) {
            errmsg = DumpArgs[earg_prim_only].name;
        } else if( (rc = ArgsOptionCount(args, DumpArgs[earg_cigartype].name, &count[earg_cigartype])) != 0 ) {
            errmsg = DumpArgs[earg_cigartype].name;
        } else if( (rc = ArgsOptionCount(args, DumpArgs[earg_cigarCG].name, &count[earg_cigarCG])) != 0 || count[earg_cigarCG] > 1 ) {
            errmsg = DumpArgs[earg_cigartype].name;
        } else if( (rc = ArgsOptionCount(args, DumpArgs[earg_header].name, &count[earg_header])) != 0 ) {
            errmsg = DumpArgs[earg_header].name;
        } else if( (rc = ArgsOptionCount(args, DumpArgs[earg_noheader].name, &count[earg_noheader])) != 0 ) {
            errmsg = DumpArgs[earg_noheader].name;
        } else if( (rc = ArgsOptionCount(args, DumpArgs[earg_comment].name, &param.comments_qty)) != 0 ) {
            errmsg = DumpArgs[earg_comment].name;
        } else if( (rc = ArgsOptionCount(args, DumpArgs[earg_region].name, &count[earg_region])) != 0 ) {
            errmsg = DumpArgs[earg_region].name;
        } else if( (rc = ArgsOptionCount(args, DumpArgs[earg_distance].name, &count[earg_distance])) != 0 ) {
            errmsg = DumpArgs[earg_distance].name;
        } else if( (rc = ArgsOptionCount(args, DumpArgs[earg_seq_id].name, &count[earg_seq_id])) != 0 ) {
            errmsg = DumpArgs[earg_seq_id].name;
        } else if( (rc = ArgsOptionCount(args, DumpArgs[earg_identicalbases].name, &count[earg_identicalbases])) != 0 ) {
            errmsg = DumpArgs[earg_identicalbases].name;
        } else if( (rc = ArgsOptionCount(args, DumpArgs[earg_gzip].name, &count[earg_gzip])) != 0 ) {
            errmsg = DumpArgs[earg_gzip].name;
        } else if( (rc = ArgsOptionCount(args, DumpArgs[earg_bzip2].name, &count[earg_bzip2])) != 0 ) {
            errmsg = DumpArgs[earg_bzip2].name;
        } else if( (rc = ArgsOptionCount(args, DumpArgs[earg_qname].name, &count[earg_qname])) != 0 ) {
            errmsg = DumpArgs[earg_qname].name;
        } else if( (rc = ArgsOptionCount(args, DumpArgs[earg_fastq].name, &count[earg_fastq])) != 0 ) {
            errmsg = DumpArgs[earg_fastq].name;
        } else if( (rc = ArgsOptionCount(args, DumpArgs[earg_fasta].name, &count[earg_fasta])) != 0 ) {
            errmsg = DumpArgs[earg_fasta].name;
        } else if( (rc = ArgsOptionCount(args, DumpArgs[earg_reverse].name, &count[earg_reverse])) != 0 ) {
            errmsg = DumpArgs[earg_reverse].name;
        } else if( (rc = ArgsOptionCount(args, DumpArgs[earg_prefix].name, &count[earg_prefix])) != 0 || count[earg_prefix] > 1) {
            rc = rc ? rc : RC(rcExe, rcArgv, rcParsing, rcParam, rcExcessive);
            errmsg = DumpArgs[earg_prefix].name;
        } else if( count[earg_prefix] > 0 && (rc = ArgsOptionValue(args, DumpArgs[earg_prefix].name, 0, &param.name_prefix)) != 0 ) {
            errmsg = DumpArgs[earg_prefix].name;
        } else if( (rc = ArgsOptionCount(args, DumpArgs[earg_test_rows].name, &count[earg_test_rows])) != 0 || count[earg_test_rows] > 1 ) {
            errmsg = DumpArgs[earg_test_rows].name;
        } else if( (rc = ArgsOptionCount(args, DumpArgs[earg_mate_row_gap_cachable].name, &count[earg_mate_row_gap_cachable])) != 0 || count[earg_mate_row_gap_cachable] > 1 ) {
            errmsg = DumpArgs[earg_mate_row_gap_cachable].name;
        } else if( (rc = ArgsOptionCount(args, DumpArgs[earg_cg_style].name, &count[earg_cg_style])) != 0 || count[earg_cg_style] > 1 ) {
            errmsg = DumpArgs[earg_cg_style].name;
        } else {
            uint32_t p, i;
            const char* arg;

            if( param.comments_qty > 0 ) {
                param.comments = calloc(param.comments_qty, sizeof(*param.comments));
                if( param.comments == NULL ) {
                    rc = RC(rcExe, rcArgv, rcProcessing, rcMemory, rcExhausted);
                }
                for(i = 0; rc == 0 && i < param.comments_qty; i++) {
                    if( (rc = ArgsOptionValue(args, DumpArgs[earg_comment].name, i, &arg)) == 0 ) {
                        param.comments[i] = arg;
                    }
                }
            }

            if( count[earg_test_rows] > 0 ) {
                if( (rc = ArgsOptionValue(args, DumpArgs[earg_test_rows].name, 0, &arg)) == 0 ) {
                    param.test_rows = strtou32(arg, NULL, 10);
                    if( param.test_rows > 0 ) {
                        SAM_DUMP_DBG(2, ("output only 1st %u rows from each table\n", param.test_rows));
                    }
                }
            }
            if( count[earg_mate_row_gap_cachable] > 0 ) {
                if( (rc = ArgsOptionValue(args, DumpArgs[earg_mate_row_gap_cachable].name, 0, &arg)) == 0 ) {
                    param.mate_row_gap_cachable = strtou32(arg, NULL, 10);
                }
            } else {
                param.mate_row_gap_cachable = 10000;
            }
            if( param.mate_row_gap_cachable > 0 ) {
                SAM_DUMP_DBG(2, ("--%s: mate cache disabled for mates with row distance less than %li\n",
                                 DumpArgs[earg_mate_row_gap_cachable].name, param.mate_row_gap_cachable));
            }
            for(p = 0; rc == 0 && p < count[earg_region]; p++) {
                /* name[:[from][-[to]]] 1-based!!! */
                TAlignedRegion r;

                errmsg = DumpArgs[earg_region].name;
                if( (rc = ArgsOptionValue(args, DumpArgs[earg_region].name, p, &arg)) == 0 ) {
                    const char* c = strchr(arg, ':');
                    if( c == NULL ) {
                        strncpy(r.name, arg, sizeof(r.name));
                        r.rq = 0;
                    } else {
                        INSDC_coord_zero* v;

                        r.r[0].from = (c - arg) > sizeof(r.name) ? sizeof(r.name) : c - arg;
                        strncpy(r.name, arg, r.r[0].from);
                        r.name[r.r[0].from] = '\0';
                        r.rq = 1;
                        r.r[0].from = -1;
                        r.r[0].to = -1;
                        r.max_to = 0;
                        v = &r.r[0].from;
                        while(rc == 0 && *++c != '\0') {
                            if( *c == '-' ) {
                                v = &r.r[0].to;
                            } else if( *c == '+' ) {
                                if( *v != 0 ) {
                                    rc = RC(rcExe, rcArgv, rcProcessing, rcParam, rcOutofrange);
                                }
                            } else if( !isdigit(*c) ) {
                                rc = RC(rcExe, rcArgv, rcProcessing, rcParam, rcOutofrange);
                            } else {
                                if( *v == -1 ) {
                                    *v = 0;
                                }
                                *v = *v * 10 + (*c - '0');
                            }
                        }
                        /* convert to 0-based offset */
                        if( r.r[0].from > 0 ) {
                            r.r[0].from--;
                        } else if( r.r[0].to > 0 ) {
                            r.r[0].from = 0;
                        }
                        if(r.r[0].to > 0 ) {
                            r.r[0].to--;
                        } else if( r.r[0].from >= 0 && r.r[0].to < 0 ) {
                            r.r[0].to = r.r[0].from;
                        }
                        if( r.r[0].from < 0 && r.r[0].to < 0 ) {
                            r.rq = 0;
                        } else if( r.r[0].from > r.r[0].to ) {
                            uint64_t x = r.r[0].from;
                            r.r[0].from = r.r[0].to;
                            r.r[0].to = x;
                        }
                    }
                    if( rc == 0 ) {
                        TAlignedRegion* x = NULL;
                        for(i = 0; i < param.region_qty; i++) {
                            if( strcmp(param.region[i].name, r.name) == 0 ) {
                                x = &param.region[i];
                                break;
                            }
                        }
                        if( x == NULL ) {
                            if( (x = realloc(param.region, sizeof(*param.region) * ++param.region_qty)) == NULL ) {
                                rc = RC(rcExe, rcArgv, rcProcessing, rcMemory, rcExhausted);
                            } else {
                                param.region = x;
                                memcpy(&param.region[param.region_qty - 1], &r, sizeof(r));
                            }
                        } else {
                            int32_t k = x->rq;
                            for(i = 0; i < x->rq; i++) {
                                /* sort by from asc */
                                if( r.r[0].from <= x->r[i].from ) {
                                    k = i;
                                    break;
                                }
                            }
                            if( k >= 0 ) {
                                /* insert at k position */
                                if( x->rq >= sizeof(x->r) / sizeof(x->r[0]) ) {
                                    rc = RC(rcExe, rcArgv, rcProcessing, rcBuffer, rcInsufficient);
                                } else {
                                    memmove(&x->r[k + 1], &x->r[k], sizeof(x->r[0]) * (x->rq - k));
                                    x->r[k].from = r.r[0].from;
                                    x->r[k].to = r.r[0].to;
                                    x->rq++;
                                }
                            }
                        }
                    }
                }
            }
            for(p = 0; p < param.region_qty; p++) {
                SAM_DUMP_DBG(2, ("filter by %s\n", param.region[p].name));
                if( param.region[p].rq == 0 ) {
                    param.region[p].rq = 1;
                    param.region[p].r[0].from = 0;
                    param.region[p].r[0].to = 0x7FFFFFFF;
                }
                for(i = 0; i < param.region[p].rq; i++) {
                    SAM_DUMP_DBG(2, ("   range: [%u:%u]\n", param.region[p].r[i].from, param.region[p].r[i].to));
                    if( param.region[p].max_to < param.region[p].r[i].to ) {
                        param.region[p].max_to = param.region[p].r[i].to;
                    }
                }
            }
            for(p = 0; rc == 0 && p < count[earg_distance]; p++) {
                /* from[-to] | [from]-to | unknown */
                errmsg = DumpArgs[earg_distance].name;
                if( (rc = ArgsOptionValue(args, DumpArgs[earg_distance].name, p, &arg)) != 0 ) {
                } else if( strcasecmp(arg, "unknown") == 0 ) {
                    param.mp_dist_unknown = true;
                } else {
                    TMatepairDistance* p;
                    if( (p = realloc(param.mp_dist, sizeof(*param.mp_dist) * ++param.mp_dist_qty)) == NULL ) {
                        rc = RC(rcExe, rcArgv, rcProcessing, rcMemory, rcExhausted);
                    } else {
                        uint64_t* v;
                        param.mp_dist = p;
                        p = &param.mp_dist[param.mp_dist_qty - 1];
                        p->from = 0;
                        p->to = 0;
                        v = &p->from;
                        while(rc == 0 && *arg != '\0') {
                            if( *arg == '-' ) {
                                v = &p->to;
                            } else if( *arg == '+' ) {
                                if( *v != 0 ) {
                                    rc = RC(rcExe, rcArgv, rcProcessing, rcParam, rcOutofrange);
                                }
                            } else if( !isdigit(*arg) ) {
                                rc = RC(rcExe, rcArgv, rcProcessing, rcParam, rcOutofrange);
                            } else {
                                *v = *v * 10 + (*arg - '0');
                            }
                            arg++;
                        }
                        if( p->from > p->to && p->to != 0 ) {
                            uint64_t x = p->from;
                            p->from = p->to;
                            p->to = x;
                        }
                        if( p->from == 0 && p->to == 0 ) {
                            param.mp_dist_qty--;
                        }
                    }
                }
            }
            for(p = 0; p < param.mp_dist_qty; p++) {
                if( param.mp_dist_unknown ) {
                    SAM_DUMP_DBG(2, ("distance 'unknown'\n"));
                }
                SAM_DUMP_DBG(2, ("distance [%lu-%lu]\n", param.mp_dist[p].from, param.mp_dist[p].to));
            }
            param.unaligned = count[earg_unaligned] > 0;
            param.only_primaries = count[earg_prim_only] > 0;
            param.long_cigar = count[earg_cigartype] > 0;
            param.reheader = count[earg_header] > 0;
            param.use_seqid = (count[earg_seq_id] > 0) || (pcount > 1);
            param.hide_identical = count[earg_identicalbases] > 0;
            param.fasta = count[earg_fasta] > 0;
            param.fastq = count[earg_fastq] > 0;
            param.reverse_unaligned = count[earg_reverse] > 0;
            param.spot_group_in_name = (count[earg_qname] > 0) || (pcount > 1);
            param.noheader = (count[earg_noheader] > 0) || param.fasta || param.fastq;
            
            if (count[earg_cigarCG] == 0 && count[earg_cg_style] == 0) {
                param.cg_style = 0;
            }
            else if (count[earg_cigarCG] == 0) {
                param.cg_style = 1;
            }
            else if (count[earg_cg_style] == 0) {
                param.cg_style = 2;
            }
            else {
                rc = RC(rcExe, rcArgv, rcProcessing, rcParam, rcInconsistent);
                errmsg = "cigar-CG and CG are mutually exclusive";
                param.cg_style = 0;
            }
            
            if( rc == 0 ) {
                /* prep alignment table cursor columns */
                g_alg_col_tmpl[alg_MATE_REF_NAME].name = param.use_seqid ? "MATE_REF_SEQ_ID" : "MATE_REF_NAME";
                if( param.fasta || param.fastq ) {
                    g_alg_col_tmpl[alg_READ].name = "RAW_READ";
                } else {
                    g_alg_col_tmpl[alg_READ].name = param.hide_identical? "MISMATCH_READ" : "READ";
                }
                g_alg_col_tmpl[alg_CIGAR].name = param.long_cigar ? "CIGAR_LONG" : "CIGAR_SHORT";
            }

            if( rc == 0 ) {
                rc = BufferedWriterMake(count[earg_gzip] > 0, count[earg_bzip2] > 0);
            }

            for(p = 0; rc == 0 && p < pcount; p++) {
                char* arg;
                if( (rc = ArgsParamValue(args, p, (const char**)&arg)) == 0 ) {
                    int i;

                    ReportResetObject(arg);

                    /* remove trailing /\ */
                    for(i = strlen(arg) - 1; i >= 0; i--) {
                        if( arg[i] != '/' && arg[i] != '\\' ) {
                            break;
                        }
                        arg[i] = '\0';
                    }
                    if( (rc = ResolvePath(arg, &param.path)) == 0 ) {
                        /* use last path element as accession */
                        param.accession = param.path;
                        for(i = strlen(param.path) - 1; i >= 0; i--) {
                            if( param.path[i] == '/' || param.path[i] == '\\' ) {
                                param.accession = &param.path[i + 1];
                                break;
                            }
                        }
                        rc = Dump(p, pcount > 1, &error_reported);
                    }
                }
            }
            BufferedWriterRelease(rc == 0);
            ResolvePath(NULL, NULL);
        }
        ArgsWhack(args);
    }
    if( rc != 0 && !error_reported ) {
        if( errmsg[0] ) {
            LOGERR(klogErr, rc, errmsg);
        } else if( KLogLastErrorCode() != rc ) {
            LOGERR(klogErr, rc, "stop");
        }
    }
    {
        /* Report execution environment if necessary */
        rc_t rc2 = ReportFinalize(rc);
        if (rc == 0)
        {   rc = rc2; }
    }
    return rc;
}
