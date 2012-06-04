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
#include <klib/text.h>
#include <klib/printf.h>
#include <klib/log.h>
#include <kdb/manager.h>
#include <kfg/config.h>
#include <insdc/insdc.h>
#include <vdb/database.h>
#include <vdb/vdb-priv.h>
#include <vdb/cursor.h>
#include <align/refseq-mgr.h>
#include <os-native.h>
#include <sysalloc.h>

#include "refseq-mgr-priv.h"
#include "debug.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <assert.h>

struct RefSeqMgr {
    BSTree tree;
    KConfig *kfg;
    const VDBManager* vmgr;
    uint32_t reader_options;
    size_t cache;
    uint32_t num_open_max;
    uint32_t num_open;
    uint64_t usage;
    RefSeq* last_used;
    BSTree vdbs;
};

static
rc_t RefSeqMgr_KfgReadStr(const KConfig* kfg, const char* path, char* value, size_t value_sz)
{
    rc_t rc = 0;
    const KConfigNode *node;

    if ( (rc = KConfigOpenNodeRead(kfg, &node, path)) == 0 ) {
        size_t num_read, remaining;
        if( (rc = KConfigNodeRead(node, 0, value, value_sz - 1, &num_read, &remaining)) == 0 ) {
            if( remaining != 0 ) {
                rc = RC(rcAlign, rcIndex, rcConstructing, rcString, rcTooLong);
            } else {
                value[num_read] = '\0';
            }
        }
        KConfigNodeRelease(node);
    } else if( GetRCState(rc) == rcNotFound ) {
        rc = 0;
        value[0] = '\0';
    }
    return rc;
}

rc_t RefSeqMgr_ForEachVolume(const RefSeqMgr* cself, RefSeqMgr_ForEachVolume_callback cb, void *data)
{
    rc_t rc = 0;
    char servers[4096];
    char volumes[4096];

    if( cself == NULL || cb == NULL ) {
        rc = RC(rcAlign, rcType, rcConstructing, rcParam, rcNull);
    } else if( cb(".", NULL, data) ) {
        /* found in local dir */
    } else if( (rc = RefSeqMgr_KfgReadStr(cself->kfg, "refseq/paths", servers, sizeof(servers))) != 0 ) {
        ALIGN_DBGERRP("%s", rc, "RefSeqMgr_KfgReadStr(paths)");
    } else {
        bool found = false;
        if( servers[0] != '\0' ) {
            char *srv_sep;
            char *srv_rem = servers;
            do {
                char const* server = srv_rem;

                srv_sep = strchr(server, ':');
                if(srv_sep) {
                    srv_rem = srv_sep + 1;
                    *srv_sep = 0;
                }
                if( cb(server, NULL, data) ) {
                    found = true;
                    break;
                }
            } while(srv_sep);
        }
        if( !found ) {
            if( (rc = RefSeqMgr_KfgReadStr(cself->kfg, "refseq/servers", servers, sizeof(servers))) != 0 ||
                (rc = RefSeqMgr_KfgReadStr(cself->kfg, "refseq/volumes", volumes, sizeof(volumes))) != 0 ) {
                ALIGN_DBGERRP("%s", rc, "RefSeqMgr_KfgReadStr(servers/volumes)");
            } else {
                /* servers and volumes are deprecated and optional */
                if( servers[0] != '\0' || volumes[0] != '\0' ) {
                    char *srv_sep;
                    char *srv_rem = servers;
                    do {
                        char vol[4096];
                        char const *server = srv_rem;
                        char *vol_rem = vol;
                        char *vol_sep;

                        strcpy(vol, volumes);
                        srv_sep = strchr(server, ':');
                        if(srv_sep) {
                            srv_rem = srv_sep + 1;
                            *srv_sep = 0;
                        }
                        do {
                            char const *volume = vol_rem;

                            vol_sep = strchr(volume, ':');
                            if(vol_sep) {
                                vol_rem = vol_sep + 1;
                                *vol_sep = 0;
                            }
                            found = cb(server, volume, data);
                        } while(!found && vol_sep);
                    } while(!found && srv_sep);
                }
            }
        }
    }
    return rc;
}

struct RefSeq {
    BSTNode dad;
    const RefSeqMgr* mgr;
    const TableReaderRefSeq* reader;
    uint64_t usage;
    uint32_t accession_sz;
    /* must be last member of struct */
    char accession[1];
};

typedef struct RefSeqMgr_Db_struct {
    BSTNode dad;
    char key[4096];
    const VDatabase* db;
} RefSeqMgr_Db;

struct FindTable_ctx {
    RefSeqMgr* self;
    const KDBManager* kmgr;
    const VTable** tbl;
    char const* name;
    uint32_t name_sz;
    bool found;
    rc_t rc;
    char** path; /* optional */
};

static
int CC RefSeqMgr_DbSort(const BSTNode* item, const BSTNode* node)
{
    return strcmp(((const RefSeqMgr_Db*)item)->key, ((const RefSeqMgr_Db*)node)->key);
}

static
int CC RefSeqMgr_FindDb(const void *item, const BSTNode *node)
{
    return strcmp((const char*)item, ((const RefSeqMgr_Db*)node)->key);
}

static
bool FindTable(char const server[], char const volume[], void *data)
{
    struct FindTable_ctx *ctx = data;
    int pathType;
    const char* vol_sep = "/";
    
    if( volume == NULL ) {
        volume = "";
        vol_sep = "";
    }
    ALIGN_DBG("try table '%s%s%s/%.*s'", server, vol_sep, volume, ctx->name_sz, ctx->name);
    pathType = KDBManagerPathType(ctx->kmgr, "%s%s%s/%.*s", server, vol_sep, volume, ctx->name_sz, ctx->name);
    pathType &= ~kptAlias;

    if( pathType == kptTable ) {
        ALIGN_DBG("found table '%s%s%s/%.*s'", server, vol_sep, volume, ctx->name_sz, ctx->name);
        ctx->found = true;
        if (ctx->tbl) {
            ctx->rc = VDBManagerOpenTableRead(ctx->self->vmgr, ctx->tbl, NULL, "%s%s%s/%.*s", server, vol_sep, volume, ctx->name_sz, ctx->name);
        }
    } else if( pathType == kptNotFound ) {
        /* can be kar */
        if( ctx->name_sz > 7 ) {
            uint32_t i = 0;

            /* check for pattern '\w{4}\d{2}[\.\d]+' */
            while( isalpha(ctx->name[i]) && i < 4 ) {
                i++;
            }
            if( i == 4 && isdigit(ctx->name[i]) && isdigit(ctx->name[++i]) ) {
                while( ++i < ctx->name_sz ) {
                    if( !isdigit(ctx->name[i]) && ctx->name[i] != '.' ) {
                        break;
                    }
                }
            }
            if( i == ctx->name_sz ) {
                char key[4096];
                size_t n;

                if( string_printf(key, sizeof(key), &n, "%s%s%s/%.*s", server, vol_sep, volume, 6, ctx->name) != 0 ) {
                    LOGMSG(klogWarn, "kar-vdb lookup buffer");
                } else {
                    RefSeqMgr_Db* db = NULL;
                    ALIGN_DBG("retry as '%s/%.*s'", key, ctx->name_sz, ctx->name);
                    db = (RefSeqMgr_Db*)BSTreeFind(&ctx->self->vdbs, key, RefSeqMgr_FindDb);
                    if( db == NULL ) {
                        const VDatabase* vdb = NULL;
                        if( VDBManagerOpenDBRead(ctx->self->vmgr, &vdb, NULL, key) == 0 ) {
                            db = malloc(sizeof(*db));
                            if( db == NULL ) {
                                ctx->rc = RC(rcAlign, rcDatabase, rcOpening, rcMemory, rcExhausted);
                            } else {
                                strcpy(db->key, key);
                                db->db = vdb;
                                ctx->rc = BSTreeInsertUnique(&ctx->self->vdbs, &db->dad, NULL, RefSeqMgr_DbSort);
                            }
                        }
                        if( ctx->rc != 0 ) {
                            VDatabaseRelease(vdb);
                            free(db);
                            db = NULL;
                        }
                    }
                    if( db != NULL ) {
                        const VTable* tbl = NULL;

                        if( VDatabaseOpenTableRead(db->db, &tbl, "%.*s", ctx->name_sz, ctx->name) == 0 ) {
                            ctx->found = true;
                            if(ctx->tbl) {
                                *(ctx->tbl) = tbl;
                            } else {
                                VTableRelease(tbl);
                            }
                            if( ctx->path != NULL ) {

                            }
                        }
                    }
                }
            }
        }
    }
    if( ctx->found && ctx->path != NULL ) {
        size_t path_sz = strlen(server) + strlen(vol_sep) + strlen(volume) + 1 + ctx->name_sz + 1;
        char* path = malloc(path_sz);
        if( path == NULL ) {
            ctx->rc = RC(rcAlign, rcPath, rcConstructing, rcMemory, rcExhausted);
        } else {
            if( (ctx->rc = string_printf(path, path_sz, &path_sz, "%s%s%s/%.*s", server, vol_sep, volume, ctx->name_sz, ctx->name)) == 0 ) {
                *(ctx->path) = path;
            } else {
                free(path);
            }
        }
    }
    return ctx->found;
}

static
rc_t RefSeqMgr_FindTable(const RefSeqMgr* cself, char const accession[], uint32_t accession_sz, VTable const **tbl, char** path)
{
    rc_t rc;
    struct FindTable_ctx ctx;
    
    if( (rc = VDBManagerOpenKDBManagerRead(cself->vmgr, &ctx.kmgr)) != 0 ) {
        ALIGN_DBGERRP("%s", rc, "VDBManagerOpenKDBManagerRead");
    } else {
        ctx.self = (RefSeqMgr*)cself;
        ctx.name = accession;
        ctx.name_sz = accession_sz;
        ctx.found = false;
        ctx.tbl = tbl;
        ctx.rc = 0;
        ctx.path = path;
        
        rc = RefSeqMgr_ForEachVolume(cself, FindTable, &ctx);
        if(rc == 0 && ctx.rc == 0 && !ctx.found) {
            rc = RC(rcAlign, rcTable, rcOpening, rcTable, rcNotFound);
        }
        KDBManagerRelease(ctx.kmgr);
    }
    return rc ? rc : ctx.rc;
}

LIB_EXPORT rc_t CC RefSeqMgr_Make(const RefSeqMgr** cself, const VDBManager* vmgr, uint32_t reader_options, size_t cache, uint32_t keep_open_num)
{
    rc_t rc = 0;
    RefSeqMgr* obj = NULL;
    
    if( cself == NULL || vmgr == NULL ) {
        rc = RC(rcAlign, rcIndex, rcConstructing, rcParam, rcNull);
    } else if( (obj = calloc(1, sizeof(*obj))) == NULL ) {
        rc = RC(rcAlign, rcIndex, rcConstructing, rcMemory, rcExhausted);
    } else if( (rc = KConfigMake(&obj->kfg, NULL)) == 0 &&
               (rc = VDBManagerAddRef(vmgr)) == 0 ) {
        BSTreeInit(&obj->tree);
        BSTreeInit(&obj->vdbs);
        obj->vmgr = vmgr;
        obj->cache = cache;
        obj->num_open_max = keep_open_num;
        obj->reader_options = reader_options;
    }
    if( rc == 0 ) {
        *cself = obj;
        ALIGN_DBG("max_open: %u", obj->num_open_max);
    } else {
        RefSeqMgr_Release(obj);
        ALIGN_DBGERR(rc);
    }
    return rc;
}

static
void CC RefSeqMgr_RefSeqWhack(BSTNode *n, void *data)
{
    RefSeq* self = (RefSeq*)n;
    TableReaderRefSeq_Whack(self->reader);
    free(self);
}

static
void CC RefSeqMgr_VDbRelease(BSTNode *n, void *data)
{
    RefSeqMgr_Db* self = (RefSeqMgr_Db*)n;
    VDatabaseRelease(self->db);
    free(self);
}

LIB_EXPORT rc_t CC RefSeqMgr_Release(const RefSeqMgr* cself)
{
    if( cself != NULL ) {
        RefSeqMgr* self = (RefSeqMgr*)cself;
        BSTreeWhack(&self->tree, RefSeqMgr_RefSeqWhack, NULL);
        BSTreeWhack(&self->vdbs, RefSeqMgr_VDbRelease, NULL);
        VDBManagerRelease(self->vmgr);
        KConfigRelease(self->kfg);
        free(self);
    }
    return 0;
}

LIB_EXPORT rc_t RefSeqMgr_Exists(const RefSeqMgr* cself, const char* accession, uint32_t accession_sz, char** path)
{
    rc_t rc = 0;

    if( cself == NULL || accession == NULL || accession_sz == 0 ) {
        rc = RC(rcAlign, rcIndex, rcAccessing, rcParam, rcNull);
    } else {
        rc = RefSeqMgr_FindTable(cself, accession, accession_sz, NULL, path);
    }
    return rc;
}

static
int CC RefSeq_Cmp(const void* item, const BSTNode* node)
{
    const String* i = (const String*)item;
    const RefSeq* n = (const RefSeq*)node;

    if( i->size == n->accession_sz ) {
        return strncasecmp(i->addr, n->accession, i->size);
    }
    return i->size - n->accession_sz;
}

static
int CC RefSeq_Sort(const BSTNode* item, const BSTNode* node)
{
    const RefSeq* i = (const RefSeq*)item;
    const RefSeq* n = (const RefSeq*)node;

    if( i->accession_sz == n->accession_sz ) {
        return strncasecmp(i->accession, n->accession, i->accession_sz);
    }
    return i->accession_sz - n->accession_sz;
}

static
void CC RefSeq_Unused( BSTNode *node, void *data )
{
    RefSeq* n = (RefSeq*)node;
    RefSeq** d = (RefSeq**)data;

    if( n->reader != NULL ) {
        if( *d == NULL || (*d)->usage > n->usage ) {
            *d = n;
        }
    }
}

LIB_EXPORT rc_t CC RefSeqMgr_Read(const RefSeqMgr* cself, const char* seq_id, uint32_t seq_id_sz,
                                  INSDC_coord_zero offset, INSDC_coord_len len,
                                  uint8_t* buffer, INSDC_coord_len* written)
{
    rc_t rc = 0;
    const RefSeq* obj;

    if( (rc = RefSeqMgr_GetSeq(cself, &obj, seq_id, seq_id_sz)) == 0 ) {
        rc = RefSeq_Read(obj, offset, len, buffer, written);
        RefSeq_Release(obj);
    }
    return rc;
}

LIB_EXPORT rc_t CC RefSeqMgr_GetSeq(const RefSeqMgr* cmgr, const RefSeq** cself, const char* seq_id, uint32_t seq_id_sz)
{
    rc_t rc = 0;
    RefSeq* obj = NULL;
    RefSeqMgr* mgr = (RefSeqMgr*)cmgr;

    if( cmgr == NULL || cself == NULL || seq_id == NULL ) {
        rc = RC(rcAlign, rcIndex, rcConstructing, rcParam, rcNull);
    } else {
        String s;
        StringInit(&s, seq_id, seq_id_sz, seq_id_sz);
        if( cmgr->last_used != NULL && RefSeq_Cmp(&s, &cmgr->last_used->dad) == 0 ) {
            obj = cmgr->last_used;
        } else {
            *cself = NULL;
            obj = (RefSeq*)BSTreeFind(&cmgr->tree, &s, RefSeq_Cmp);
            if( obj == NULL ) {
                if( (obj = calloc(1,  sizeof(*obj) + seq_id_sz)) == NULL ) {
                    rc = RC(rcAlign, rcIndex, rcConstructing, rcMemory, rcExhausted);
                } else {
                    memcpy(obj->accession, seq_id, seq_id_sz);
                    obj->accession_sz = seq_id_sz;
                    obj->mgr = cmgr;
                    if( (rc = BSTreeInsertUnique(&mgr->tree, &obj->dad, NULL, RefSeq_Sort)) != 0 ) {
                        free(obj);
                    }
                }
            }
        }
    }
    if( rc == 0 ) {
        if( obj->reader == NULL ) {
            const VTable* tbl;
            bool retry = false;

            do {
                if( retry || (cmgr->num_open_max > 0 && cmgr->num_open >= cmgr->num_open_max) ) {
                    RefSeq* old = NULL;
                    BSTreeForEach(&cmgr->tree, false, RefSeq_Unused, &old);
                    if( old != NULL ) {
                        TableReaderRefSeq_Whack(old->reader);
                        old->reader = NULL;
                        mgr->num_open--;
                    }
                }
                if( (rc = RefSeqMgr_FindTable(cmgr, obj->accession, obj->accession_sz, &tbl, NULL)) == 0 ) {
                    if( (rc = TableReaderRefSeq_MakeTable(&obj->reader, cmgr->vmgr, tbl,
                                                          cmgr->reader_options, cmgr->cache)) == 0 ) {
                        mgr->num_open++;
                    } else if( GetRCObject(rc) == rcMemory &&
                              (GetRCState(rc) == rcInsufficient || GetRCState(rc) == rcExhausted) ) {
                       retry = mgr->num_open > 0;
                       rc = 0;
                    }
                    VTableRelease(tbl);
                }
            } while( rc == 0 && retry);
        }
        if( rc == 0 ) {
            mgr->last_used = obj;
        }
    }
    if( rc == 0 ) {
        *cself = obj;
    } else {
        ALIGN_DBGERRP("SEQ_ID: '%.*s'", rc, seq_id_sz, seq_id);
    }
    return rc;
}

LIB_EXPORT rc_t CC RefSeq_Read(const RefSeq* cself, INSDC_coord_zero offset, INSDC_coord_len len,
                               uint8_t* buffer, INSDC_coord_len* written)
{
    rc_t rc = 0;

    if( cself == NULL || buffer == NULL || written == NULL ) {
        rc = RC(rcAlign, rcFile, rcReading, rcParam, rcNull);
    } else if( (rc = TableReaderRefSeq_Read(cself->reader, offset, len, buffer, written)) == 0 ) {
        ((RefSeq*)cself)->usage = ++((RefSeqMgr*)cself->mgr)->usage;
    }
    ALIGN_DBGERR(rc);
    return rc;
}

LIB_EXPORT rc_t CC RefSeq_Circular(const RefSeq* cself, bool* circular)
{
    rc_t rc = 0;

    if( cself == NULL ) {
        rc = RC(rcAlign, rcFile, rcReading, rcParam, rcNull);
    } else {
        rc = TableReaderRefSeq_Circular(cself->reader, circular);
    }
    ALIGN_DBGERR(rc);
    return rc;
}

LIB_EXPORT rc_t CC RefSeq_SeqLength(const RefSeq* cself, INSDC_coord_len* len)
{
    rc_t rc = 0;

    if( cself == NULL ) {
        rc = RC(rcAlign, rcFile, rcReading, rcParam, rcNull);
    } else {
        rc = TableReaderRefSeq_SeqLength(cself->reader, len);
    }
    ALIGN_DBGERR(rc);
    return rc;
}

LIB_EXPORT rc_t CC RefSeq_MD5(const RefSeq* cself, const uint8_t** md5)
{
    rc_t rc = 0;

    if( cself == NULL ) {
        rc = RC(rcAlign, rcFile, rcReading, rcParam, rcNull);
    } else {
        rc = TableReaderRefSeq_MD5(cself->reader, md5);
    }
    ALIGN_DBGERR(rc);
    return rc;
}

LIB_EXPORT rc_t CC RefSeq_Release(const RefSeq* cself)
{
    return 0;
}

