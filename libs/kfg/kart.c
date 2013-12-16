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

#include <kfg/kart.h>

#include <kfs/directory.h> /* KDirectoryOpenFileRead */
#include <kfs/file.h> /* KFile */
#include <kfs/gzip.h> /* KFileMakeGzipForRead */
#include <kfs/subfile.h> /* KFileMakeSubRead */

#include <klib/data-buffer.h> /* KDataBuffer */
#include <klib/rc.h>
#include <klib/refcount.h> /* KRefcount */
#include <klib/out.h> /* OUTMSG */

#include <strtol.h> /* strtou64 */
#include <sysalloc.h>

#include <assert.h>
#include <stdlib.h> /* free */
#include <string.h> /* memcmp */

#define RELEASE(type, obj) do { rc_t rc2 = type##Release(obj); \
    if (rc2 != 0 && rc == 0) { rc = rc2; } obj = NULL; } while (false)

struct KartItem {
    KRefcount refcount;

    const Kart *dad;

/*  String typeId; */
    String projId;
    String itemId;
    String accession;
    String name;
    String itemDesc;
};

static void KartItemWhack(KartItem *self) {
    assert(self);

    KartRelease(self->dad);

    memset(self, 0, sizeof *self);

    free(self);
}

/* AddRef
 * Release
 *  all objects are reference counted
 *  NULL references are ignored
 */
LIB_EXPORT rc_t CC KartItemAddRef(const KartItem *self) {
    if (self != NULL) {
        switch (KRefcountAdd(&self->refcount, "KartItem")) {
            case krefLimit:
                return RC(rcKFG, rcFile, rcAttaching, rcRange, rcExcessive);
        }
    }

    return 0;
}

LIB_EXPORT rc_t CC KartItemRelease(const KartItem *self) {
    if (self != NULL) {
        switch (KRefcountDrop(&self -> refcount, "KartItem")) {
            case krefWhack:
                KartItemWhack((KartItem*)self);
                break;
            case krefLimit:
                return RC(rcKFG, rcFile, rcReleasing, rcRange, rcExcessive);
        }
    }

    return 0;
}

static rc_t StringAsUint64(const String *self, uint64_t *pid) {
    uint64_t id = 0;

    char buffer[21] = "";
    size_t bytes = 0;
    char *end = NULL;

    assert(self);

    if (pid == NULL) {
        return RC(rcKFG, rcFile, rcAccessing, rcParam, rcNull);
    }

    *pid = 0;

    if (sizeof buffer - 1 < self->len) {
        return RC(rcKFG, rcFile, rcAccessing, rcBuffer, rcInsufficient);
    }

    bytes = string_copy(buffer, sizeof buffer, self->addr, self->len);
    if (bytes != self->len) {
        return RC(rcKFG, rcFile, rcAccessing, rcBuffer, rcInsufficient);
    }

    id = strtou64(buffer, &end, 0);
    if (end[0] != 0) {
        return RC(rcKFG, rcFile, rcAccessing, rcParam, rcInvalid);
    }

    *pid = id;

    return 0;
}

LIB_EXPORT rc_t CC KartItemProjIdNumber(const KartItem *self, uint64_t *pid) {
    if (self == NULL) {
        return RC(rcKFG, rcFile, rcAccessing, rcSelf, rcNull);
    }
    return StringAsUint64(&self->projId, pid);
}

LIB_EXPORT rc_t CC KartItemItemIdNumber(const KartItem *self, uint64_t *pid) {
    if (self == NULL) {
        return RC(rcKFG, rcFile, rcAccessing, rcSelf, rcNull);
    }
    return StringAsUint64(&self->itemId, pid);
}

static rc_t KartItemCheck(const KartItem *self, const String **elem) {
    if (elem == NULL) {
        return RC(rcKFG, rcFile, rcAccessing, rcParam, rcNull);
    }

    *elem = NULL;

    if (self == NULL) {
        return RC(rcKFG, rcFile, rcAccessing, rcSelf, rcNull);
    }

    return 0;
}

LIB_EXPORT rc_t CC KartItemProjId(const KartItem *self, const String **elem)
{
    rc_t rc = KartItemCheck(self, elem);
    if (rc == 0) {
        *elem = &self->projId;
    }
    return rc;
}
LIB_EXPORT rc_t CC KartItemItemId(const KartItem *self, const String **elem)
{
    rc_t rc = KartItemCheck(self, elem);
    if (rc == 0) {
        *elem = &self->itemId;
    }
    return rc;
}
LIB_EXPORT rc_t CC KartItemAccession(const KartItem *self, const String **elem)
{
    rc_t rc = KartItemCheck(self, elem);
    if (rc == 0) {
        *elem = &self->accession;
    }
    return rc;
}
LIB_EXPORT rc_t CC KartItemName(const KartItem *self, const String **elem)
{
    rc_t rc = KartItemCheck(self, elem);
    if (rc == 0) {
        *elem = &self->name;
    }
    return rc;
}
LIB_EXPORT rc_t CC KartItemItemDesc(const KartItem *self, const String **elem)
{
    rc_t rc = KartItemCheck(self, elem);
    if (rc == 0) {
        *elem = &self->itemDesc;
    }
    return rc;
}
/*LIB_EXPORT rc_t CC KartItemTypeId(const KartItem *self, const String **elem)
{
    rc_t rc = KartItemCheck(self, elem);
    if (rc == 0) {
        *elem = &self->typeId;
    }
    return rc;
}*/

/** Print KartItem using OUTMSG; if (self == NULL) then print the header */
LIB_EXPORT rc_t CC KartItemPrint(const KartItem *self) { /* AA-833 */
    if (self != NULL) {
        return OUTMSG(("'%S'\t'%S'\t'%S'\t'%S'\t'%S'\n", &self->projId,
            &self->itemId, &self->accession, &self->name, &self->itemDesc));
    }
    return 0;
}

struct Kart {
    KRefcount refcount;

    KDataBuffer mem;

    const char *text;
    uint64_t len;
};

static void KartWhack(Kart *self) {
    assert(self);

    KDataBufferWhack(&self->mem);

    memset(self, 0, sizeof *self);

    free(self);
}

/* AddRef
 * Release
 *  all objects are reference counted
 *  NULL references are ignored
 */
LIB_EXPORT rc_t CC KartAddRef(const Kart *self) {
    if (self != NULL) {
        switch (KRefcountAdd(&self->refcount, "Kart")) {
            case krefLimit:
                return RC(rcKFG, rcFile, rcAttaching, rcRange, rcExcessive);
        }
    }

    return 0;
}

LIB_EXPORT rc_t CC KartRelease(const Kart *self) {
    if (self != NULL) {
        switch (KRefcountDrop(&self -> refcount, "Kart")) {
            case krefWhack:
                KartWhack((Kart*)self);
                break;
            case krefLimit:
                return RC(rcKFG, rcFile, rcReleasing, rcRange, rcExcessive);
        }
    }

    return 0;
}

static rc_t KartItemInitFromKartRow(const Kart *self, const KartItem **item,
    const char *line, size_t len)
{
    rc_t rc = 0;
    int i = 0;
    KartItem *obj = NULL;
    assert(self && item && line && len);
    obj = calloc(1, sizeof *obj);
    if (obj == NULL) {
        return RC(rcKFG, rcData, rcAllocating, rcMemory, rcExhausted);
    }
    for (i = 0; ; ++i) {
        size_t l = 0;
        String *next = NULL;
        const char *p = string_chr(line, len, '|');
        if (p == NULL) {
            if (i != 4) {
                return RC(rcKFG, rcFile, rcParsing, rcFile, rcInsufficient);
            }
            l = len;
        }
        else {
            l = p - line;
        }
        switch (i) { /* AA-833 */
            case 0:
                next = &obj->projId;
                break;
            case 1:
                next = &obj->itemId;
                break;
            case 2:
                next = &obj->accession;
                break;
            case 3:
                next = &obj->name;
                break;
            case 4:
                next = &obj->itemDesc;
                break;
            default:
                return RC(rcKFG, rcFile, rcParsing, rcFile, rcExcessive);
                break;
        }
        assert(next);
        StringInit(next, line, l, l);
        if (l > len) {
            return RC(rcKFG, rcFile, rcParsing, rcFile, rcInvalid);
        }
        if (len == l) {
            break;
        }
        ++l;
        line += l;
        len -= l;
    }
    rc = KartAddRef(self);
    if (rc == 0) {
        obj->dad = self;
        *item = obj;
    }
    return rc;
}

LIB_EXPORT rc_t CC KartPrint(const Kart *self) {
    uint32_t len = self->mem.elem_count;
    if (self == NULL) {
        return RC(rcKFG, rcFile, rcLoading, rcSelf, rcNull);
    }
    return OUTMSG(("%.*s", len, self->mem.base));
}

LIB_EXPORT rc_t CC KartMakeNextItem(Kart *self, const KartItem **item) {
    size_t len = 0;
    const char *line = NULL;
    const char *next = NULL;

    if (item == NULL) {
        return RC(rcKFG, rcFile, rcLoading, rcParam, rcNull);
    }
    *item = NULL;
    if (self == NULL) {
        return RC(rcKFG, rcFile, rcLoading, rcSelf, rcNull);
    }

    while (self->len > 0
        && (self->text[0] == '\r' || self->text[0] == '\n'))
    {
        ++self->text;
        --self->len;
    }

    line = self->text;
    next = string_chr(self->text, self->len, '\n');
    if (next == NULL) {
        return RC(rcKFG, rcFile, rcLoading, rcFile, rcInsufficient);
    }

    len = next - self->text;
    if (*(next - 1) == '\r') {
        --len;
    }

    if (self->len >= next - self->text + 1) {
        self->len -= next - self->text + 1;
    }
    else {
        OUTMSG(("WARNING: STRING OVERFLOW DURING KART ROW PARSING"));
        self->len = 0;
    }

    self->text = next + 1;

    {
        const char end[] = "$end";
        if (string_cmp(line, string_size(line), end,
            sizeof end - 1, sizeof end - 1) == 0)
        {
            return 0;
        }
    }

    return KartItemInitFromKartRow(self, item, line, len);
}

static rc_t decode_kart(KDataBuffer *mem, const KFile *orig, size_t hdr_sz) {
    rc_t rc = 0;
    size_t num_read;
    uint64_t eof;
    assert(mem && orig && hdr_sz);
    rc = KFileSize ( orig, & eof );
    if ( rc == 0 )
    {
        const KFile *sub;
        rc = KFileMakeSubRead(&sub, orig, hdr_sz, eof - hdr_sz);
        if ( rc == 0 )
        {
            const KFile *gzip;
            rc = KFileMakeGzipForRead ( & gzip, sub );
            if ( rc == 0 )
            {
                rc = KDataBufferMakeBytes ( mem, 0 );
                if ( rc == 0 )
                {
                    size_t total, to_read;

                    /* after all of that, we're ready to decompress */
                    for ( total = 0; ; )
                    {
                        char *buff;

                        rc = KDataBufferResize ( mem,
                            total + 32 * 1024 );
                        if ( rc != 0 )
                            break;

                        buff = mem -> base;
                        to_read = ( size_t ) mem -> elem_count - total;

                        rc = KFileReadAll ( gzip, total,
                            & buff [ total ], to_read, & num_read );
                        if ( rc != 0 )
                            break;

                        total += num_read;
                        
                        if ( num_read < to_read )
                        {
                            buff [ total ] = 0;
                            mem -> elem_count = total;
                            break;
                        }
                    }
                }

                KFileRelease ( gzip );
            }

            KFileRelease ( sub );
        }
    }

    return rc;
}

static rc_t KartProcessHeader(Kart *self) {
    assert(self);

    self->text = self->mem.base;
    self->len = self->mem.elem_count;

    {
        const char version[] = "version ";
        size_t l = sizeof version - 1;
        if (string_cmp(version, l, self->text, self->len, l) != 0) {
            return RC(rcKFG, rcMgr, rcUpdating, rcFormat, rcUnrecognized);
        }

        self->text += l;
        self->len -= l;
    }

    {
        const char version[] = "1.0";
        size_t l = sizeof version - 1;
        if (string_cmp(version, l, self->text, l, l) != 0) {
            return RC(rcKFG, rcMgr, rcUpdating, rcFormat, rcUnsupported);
        }

        self->text += l;
        self->len -= l;
    }

    while (self->len > 0 && (self->text[0] == '\r' || self->text[0] == '\n')) {
        ++self->text;
        --self->len;
    }

    return 0;
}

#ifdef _DEBUGGING
static rc_t read_textkart(KDataBuffer *mem, const KFile *orig) {
    rc_t rc = 0;
    size_t num_read;
    uint64_t eof;
    assert(mem && orig);
    rc = KFileSize ( orig, & eof );
    if ( rc == 0 )
    {
        rc = KDataBufferMakeBytes ( mem, 0 );
        if ( rc == 0 ) {
            /* after all of that, we're ready to read */
            rc = KDataBufferResize(mem, eof);
            if ( rc != 0 )
                return rc;
            rc = KFileReadAll ( orig, 0, mem -> base, eof, & num_read );
            if ( rc != 0 )
                return rc;
        }
    }
    return rc;
}
KFG_EXTERN rc_t CC KartMakeText(const struct KDirectory *dir, const char *path,
    Kart **kart, bool *isKart)
{
    rc_t rc = 0;
    const KFile *f = NULL;

    if (dir == NULL || path == NULL || kart == NULL || isKart == NULL) {
        return RC(rcKFG, rcFile, rcReading, rcParam, rcNull);
    }

    *isKart = false;
    *kart = NULL;

    rc = KDirectoryOpenFileRead(dir, &f, path);
    if (rc != 0) {
        return rc;
    }

    {
        Kart *obj = NULL;

        *isKart = true;

        obj = calloc(1, sizeof *obj);
        if (obj == NULL) {
            return RC(rcKFG, rcData, rcAllocating, rcMemory, rcExhausted);
        }

        rc = read_textkart(&obj->mem, f);
        if (rc == 0) {
            rc = KartProcessHeader(obj);
        }
        if (rc == 0) {
            KRefcountInit(&obj->refcount, 0, "Kart", "MakeText", "kart");
            *kart = obj;
        }
        else {
            KartWhack(obj);
        }
    }

    RELEASE(KFile, f);
    return rc;
}
#endif

LIB_EXPORT rc_t KartMake(const KDirectory *dir, const char *path,
    Kart **kart, bool *isKart)
{
    rc_t rc = 0;
    const KFile *f = NULL;
    char hdr[8] = "";
    size_t num_read = 0;

    if (dir == NULL || path == NULL || kart == NULL || isKart == NULL) {
        return RC(rcKFG, rcFile, rcReading, rcParam, rcNull);
    }

    *isKart = false;
    *kart = NULL;

    rc = KDirectoryOpenFileRead(dir, &f, path);
    if (rc != 0) {
        return rc;
    }

    rc = KFileReadAll(f, 0, hdr, sizeof hdr, &num_read);
    if (rc == 0 && num_read == sizeof hdr &&
        memcmp(hdr, "ncbikart", sizeof hdr) == 0)
    {
        Kart *obj = NULL;

        *isKart = true;

        obj = calloc(1, sizeof *obj);
        if (obj == NULL) {
            return RC(rcKFG, rcData, rcAllocating, rcMemory, rcExhausted);
        }

        rc = decode_kart(&obj->mem, f, sizeof hdr);
        if (rc == 0) {
            rc = KartProcessHeader(obj);
        }
        if (rc == 0) {
            KRefcountInit(&obj->refcount, 0, "Kart", "Make", "kart");
            *kart = obj;
        }
        else {
            KartWhack(obj);
        }
    }

    RELEASE(KFile, f);
    return rc;
}
