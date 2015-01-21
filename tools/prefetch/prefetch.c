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
*
*/

/********** includes **********/

#include "prefetch.vers.h"

#include <kapp/main.h> /* KAppVersion */

#include <kfg/config.h> /* KConfig */
#include <kfg/kart.h> /* Kart */
#include <kfg/repository.h> /* KRepositoryMgr */

#include <vdb/database.h> /* VDatabase */
#include <vdb/dependencies.h> /* VDBDependencies */
#include <vdb/manager.h> /* VDBManager */

#include <kdb/manager.h> /* kptDatabase */

#include <vfs/manager.h> /* VFSManager */
#include <vfs/path.h> /* VPath */
#include <vfs/resolver.h> /* VResolver */

#include <kns/ascp.h> /* ascp_locate */
#include <kns/manager.h>
#include <kns/kns-mgr-priv.h>
#include <kns/http.h>

#include <kfs/file.h> /* KFile */
#include <kfs/gzip.h> /* KFileMakeGzipForRead */
#include <kfs/subfile.h> /* KFileMakeSubRead */

#include <klib/container.h> /* BSTree */
#include <klib/data-buffer.h> /* KDataBuffer */
#include <klib/log.h> /* PLOGERR */
#include <klib/out.h> /* KOutMsg */
#include <klib/printf.h> /* string_printf */
#include <klib/rc.h>
#include <klib/status.h> /* STSMSG */
#include <klib/text.h> /* String */

#include <strtol.h> /* strtou64 */
#include <sysalloc.h>

#include <assert.h>
#include <stdlib.h> /* free */
#include <string.h> /* memset */
#include <time.h> /* time */

#include <stdio.h> /* printf */

#define DISP_RC(rc, err) (void)((rc == 0) ? 0 : LOGERR(klogInt, rc, err))

#define DISP_RC2(rc, name, msg) (void)((rc == 0) ? 0 : \
    PLOGERR(klogInt, (klogInt,rc, "$(msg): $(name)","msg=%s,name=%s",msg,name)))

#define RELEASE(type, obj) do { rc_t rc2 = type##Release(obj); \
    if (rc2 != 0 && rc == 0) { rc = rc2; } obj = NULL; } while (false)

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

#define STS_TOP 0
#define STS_INFO 1
#define STS_DBG 2
#define STS_FIN 3

#define USE_CURL 0

#define rcResolver   rcTree
static bool NotFoundByResolver(rc_t rc) {
    if (GetRCModule(rc) == rcVFS) {
        if (GetRCTarget(rc) == rcResolver) {
            if (GetRCContext(rc) == rcResolving) {
                if (GetRCState(rc) == rcNotFound) {
                    return true;
                }
            }
        }
    }
    return false;
}

typedef enum {
    eOrderSize,
    eOrderOrig
} EOrder;
typedef enum {
    eForceNo, /* do not download found and complete objects */
    eForceYes,/* force download of found and complete objects */
    eForceYES /* force download; ignore lockes */
} EForce;
typedef enum {
    eRunTypeUnknown,
    eRunTypeList,     /* list sizes */
    eRunTypeDownload,
    eRunTypeGetSize    /* download ordered by sizes */
} ERunType;
typedef struct {
    const VPath *path;
    const String *str;
} VPathStr;
typedef struct {
    BSTNode n;
    char *path;
} TreeNode;
typedef struct {
    ERunType type;
    char *name;

    VPathStr local;
    const String *cache;
    const String *remote;

    const KFile *file;
    uint64_t remoteSz;

    bool undersized; /* remoteSz < min allowed size */
    bool oversized; /* remoteSz >= max allowed size */

    bool existing;

    /* path to the resolved object : either local or cache:
    should not be released */
    const char *path;

    VPath *accession;
    uint64_t project;

    const KartItem *kartItem;

    VResolver *resolver;
} Resolved;
typedef struct {
    Args *args;
    bool check_all;

    bool list_kart;
    bool list_kart_numbered;
    bool list_kart_sized;
    EOrder order;

    const char *rows;

    EForce force;
    KConfig *cfg;
    KDirectory *dir;

    const KRepositoryMgr *repoMgr;
    const VDBManager *mgr;
    VFSManager *vfsMgr;
    KNSManager *kns;

    VResolver *resolver;

    void *buffer;
    size_t bsize;

    bool undersized; /* remoteSz < min allowed size */
    bool oversized; /* remoteSz >= max allowed size */

    BSTree downloaded;

    size_t minSize;
    size_t maxSize;
    uint64_t heartbeat;

    bool noAscp;
    bool noHttp;

    bool forceAscpFail;

    bool ascpChecked;
    const char *ascp;
    const char *asperaKey;
    String *ascpMaxRate;

#ifdef _DEBUGGING
    const char *textkart;
#endif
} Main;
typedef struct {
    /* "plain" command line argument */
    const char *desc;

    const KartItem *item;

#ifdef _DEBUGGING
    const char *textkart;
#endif

    Resolved resolved;
    int number;

    Main *main; /* just a pointer, no refcount here, don't release it */
} Item;
typedef struct {
    const char *obj;
    bool done;
    Kart *kart;
    bool isKart;
} Iterator;
typedef struct {
    BSTNode n;
    Item *i;
} KartTreeNode;
/********** String extension **********/
static rc_t StringRelease(const String *self) {
    free((String*)self);
    return 0;
}

static char* StringCheck(const String *self, rc_t rc) {
    if (rc == 0 && self != NULL
        && self->addr != NULL && self->len != 0 && self->addr[0] != '\0')
    {
        return string_dup(self->addr, self->size);
    }
    return NULL;
}

static
bool _StringIsFasp(const String *self, const char **withoutScheme)
{
    const char fasp[] = "fasp://";
    const char *dummy = NULL;

    assert(self && self->addr);

    if (withoutScheme == NULL) {
        withoutScheme = &dummy;
    }

    *withoutScheme = NULL;

    if (memcmp(self->addr, fasp, sizeof fasp - 1) == 0) {
        *withoutScheme = self->addr + sizeof fasp - 1;
        return true;
    }
    return false;
}

/********** KFile extension **********/
static
rc_t _KFileOpenRemote(const KFile **self, KNSManager *kns, const char *path)
{
    rc_t rc = 0;
    assert(self);
    if (*self != NULL) {
        return 0;
    }
    rc = KNSManagerMakeReliableHttpFile(kns, self, NULL, 0x01010000, path);
    return rc;
}

/********** KDirectory extension **********/
static rc_t _KDirectoryMkTmpPrefix(const KDirectory *self,
    const String *prefix, char *out, size_t sz)
{
    size_t num_writ = 0;
    rc_t rc = string_printf(out, sz, &num_writ, "%S.tmp", prefix);
    if (rc != 0) {
        DISP_RC2(rc, "string_printf(tmp)", prefix->addr);
        return rc;
    }

    if (num_writ > sz) {
        rc = RC(rcExe, rcFile, rcCopying, rcBuffer, rcInsufficient);
        PLOGERR(klogInt, (klogInt, rc,
            "bad string_printf($(s).tmp) result", "s=%s", prefix->addr));
        return rc;
    }

    return rc;
}

static rc_t _KDirectoryMkTmpName(const KDirectory *self,
    const String *prefix, char *out, size_t sz)
{
    rc_t rc = 0;
    int i = 0;

    assert(prefix);

    while (rc == 0) {
        size_t num_writ = 0;
        rc = string_printf(out, sz, &num_writ,
            "%S.tmp.%d.tmp", prefix, rand() % 100000);
        if (rc != 0) {
            DISP_RC2(rc, "string_printf(tmp.rand)", prefix->addr);
            break;
        }

        if (num_writ > sz) {
            rc = RC(rcExe, rcFile, rcCopying, rcBuffer, rcInsufficient);
            PLOGERR(klogInt, (klogInt, rc,
                "bad string_printf($(s).tmp.rand) result",
                "s=%s", prefix->addr));
            return rc;
        }
        if (KDirectoryPathType(self, "%s", out) == kptNotFound) {
            break;
        }
        if (++i > 999) {
            rc = RC(rcExe, rcFile, rcCopying, rcName, rcUndefined);
            PLOGERR(klogInt, (klogInt, rc,
                "cannot generate unique tmp file name for $(name)", "name=%s",
                prefix->addr));
            return rc;
        }
    }

    return rc;
}

static rc_t _KDirectoryMkLockName(const KDirectory *self,
    const String *prefix, char *out, size_t sz)
{
    rc_t rc = 0;
    size_t num_writ = 0;

    assert(prefix);

    rc = string_printf(out, sz, &num_writ, "%S.lock", prefix);
    DISP_RC2(rc, "string_printf(lock)", prefix->addr);

    if (rc == 0 && num_writ > sz) {
        rc = RC(rcExe, rcFile, rcCopying, rcBuffer, rcInsufficient);
        PLOGERR(klogInt, (klogInt, rc,
            "bad string_printf($(s).lock) result", "s=%s", prefix->addr));
        return rc;
    }

    return rc;
}

static
rc_t _KDirectoryCleanCache(KDirectory *self, const String *local)
{
    rc_t rc = 0;

    char cache[PATH_MAX] = "";
    size_t num_writ = 0;

    assert(self && local && local->addr);

    if (rc == 0) {
        rc = string_printf(cache, sizeof cache, &num_writ, "%s.cache",
            local->addr);
        DISP_RC2(rc, "string_printf(.cache)", local->addr);
    }

    if (rc == 0 && KDirectoryPathType(self, "%s", cache) != kptNotFound) {
        STSMSG(STS_DBG, ("removing %s", cache));
        rc = KDirectoryRemove(self, false, "%s", cache);
    }

    return rc;
}

static rc_t _KDirectoryClean(KDirectory *self, const String *cache,
    const char *lock, const char *tmp, bool rmSelf)
{
    rc_t rc = 0;
    rc_t rc2 = 0;

    char tmpName[PATH_MAX] = "";
    const char *dir = tmpName;
    const char *tmpPfx = NULL;
    size_t tmpPfxLen = 0;

    assert(self && cache);

    rc = _KDirectoryMkTmpPrefix(self, cache, tmpName, sizeof tmpName);
    if (rc == 0) {
        char *slash = strrchr(tmpName, '/');
        if (slash != NULL) {
            if (strlen(tmpName) == slash + 1 - tmpName) {
                rc = RC(rcExe,
                    rcDirectory, rcSearching, rcDirectory, rcIncorrect);
                PLOGERR(klogInt, (klogInt, rc,
                    "bad file name $(path)", "path=%s", tmpName));
            }
            else {
                *slash = '\0';
                tmpPfx = slash + 1;
            }
        }
        else {
            rc = RC(rcExe, rcDirectory, rcSearching, rcChar, rcNotFound);
            PLOGERR(klogInt, (klogInt, rc,
                    "cannot extract directory from $(path)", "path=%s", dir));
        }
        tmpPfxLen = strlen(tmpPfx);
    }

    if (tmp != NULL && KDirectoryPathType(self, "%s", tmp) != kptNotFound) {
        rc_t rc3 = 0;
        STSMSG(STS_DBG, ("removing %s", tmp));
        rc3 = KDirectoryRemove(self, false, "%s", tmp);
        if (rc2 == 0 && rc3 != 0) {
            rc2 = rc3;
        }
    }

    if (rmSelf && KDirectoryPathType(self, "%s", cache->addr) != kptNotFound) {
        rc_t rc3 = 0;
        STSMSG(STS_DBG, ("removing %s", cache->addr));
        rc3 = KDirectoryRemove(self, false, "%s", cache->addr);
        if (rc2 == 0 && rc3 != 0) {
            rc2 = rc3;
        }
    }

    if (rc == 0) {
        uint32_t count = 0;
        uint32_t i = 0;
        KNamelist *list = NULL;
        STSMSG(STS_DBG, ("listing %s for old temporary files", dir));
        rc = KDirectoryList(self, &list, NULL, NULL, "%s", dir);
        DISP_RC2(rc, "KDirectoryList", dir);

        if (rc == 0) {
            rc = KNamelistCount(list, &count);
            DISP_RC2(rc, "KNamelistCount(KDirectoryList)", dir);
        }

        for (i = 0; i < count && rc == 0; ++i) {
            const char *name = NULL;
            rc = KNamelistGet(list, i, &name);
            if (rc != 0) {
                DISP_RC2(rc, "KNamelistGet(KDirectoryList)", dir);
            }
            else {
                if (strncmp(name, tmpPfx, tmpPfxLen) == 0) {
                    rc_t rc3 = 0;
                    STSMSG(STS_DBG, ("removing %s", name));
                    rc3 = KDirectoryRemove(self, false,
                        "%s%c%s", dir, '/', name);
                    if (rc2 == 0 && rc3 != 0) {
                        rc2 = rc3;
                    }
                }
            }
        }

        RELEASE(KNamelist, list);
    }

    if (lock != NULL && KDirectoryPathType(self, "%s", lock) != kptNotFound) {
        rc_t rc3 = 0;
        STSMSG(STS_DBG, ("removing %s", lock));
        rc3 = KDirectoryRemove(self, false, "%s", lock);
        if (rc2 == 0 && rc3 != 0) {
            rc2 = rc3;
        }
    }

    if (rc == 0 && rc2 != 0) {
        rc = rc2;
    }

    return rc;
}

/********** VResolver extension **********/
static rc_t V_ResolverRemote(const VResolver *self,
    VRemoteProtocols protocols, struct VPath const * accession,
    struct VPath const ** remote, struct VPath const ** cache)
{
    return VResolverQuery(self, protocols, accession, NULL, remote, cache);
}

static rc_t V_ResolverLocal(const VResolver *self,
    struct VPath const * accession, struct VPath const ** path )
{
    return VResolverQuery(self, eProtocolHttp, accession, path, NULL, NULL);
}

static rc_t _VResolverRemote(VResolver *self, VRemoteProtocols protocols,
    const char *name, const VPath *vaccession,
    const VPath **vremote, const String **remote,
    const String **cache)
{
    rc_t rc = 0;
    const VPath *vcache = NULL;
    assert(vaccession && vremote && !*vremote);
    rc = V_ResolverRemote(self, protocols, vaccession, vremote, &vcache);
    if (rc == 0) {
        char path[PATH_MAX] = "";
        size_t len = 0;
        rc = VPathReadUri(*vremote, path, sizeof path, &len);
        DISP_RC2(rc, "VPathReadUri(VResolverRemote)", name);
        if (rc == 0) {
            String local_str;
            char *query = string_chr(path, len, '?');
            if (query != NULL) {
                *query = '\0';
            }
            StringInit(&local_str, path, len, (uint32_t)len);
            RELEASE(String, *remote);
            rc = StringCopy(remote, &local_str);
            DISP_RC2(rc, "StringCopy(VResolverRemote)", name);
        }
    }
    else if (NotFoundByResolver(rc)) {
        PLOGERR(klogErr, (klogErr, rc, "'$(acc)' cannot be found.",
            "acc=%s", name));
    }
    else {
        DISP_RC2(rc, "Cannot resolve remote", name);
    }
    if (rc == 0 && cache != NULL) {
        String path_str;
        if (vcache == NULL) {
            rc = RC(rcExe, rcResolver, rcResolving, rcPath, rcNotFound);
            PLOGERR(klogInt, (klogInt, rc, "cannot get cache location "
                "for $(acc).", /* Try to cd out of protected repository.", */
                "acc=%s" , name));
        }
        if (rc == 0) {
            rc = VPathGetPath(vcache, &path_str);
            DISP_RC2(rc, "VPathGetPath(VResolverCache)", name);
        }
        if (rc == 0) {
            rc = StringCopy(cache, &path_str);
            DISP_RC2(rc, "StringCopy(VResolverCache)", name);
        }
    }
    RELEASE(VPath, vcache);
    return rc;
}

/********** VPathStr **********/
static rc_t VPathStrFini(VPathStr *self) {
    rc_t rc = 0;

    assert(self);

    VPathRelease(self->path);

    RELEASE(String, self->str);

    memset(self, 0, sizeof *self);

    return rc;
}

/********** TreeNode **********/
static int CC bstCmp(const void *item, const BSTNode *n) {
    const char* path = item;
    const TreeNode* sn = (const TreeNode*) n;

    assert(path && sn && sn->path);

    return strcmp(path, sn->path);
}

static int CC bstSort(const BSTNode* item, const BSTNode* n) {
    const TreeNode* sn = (const TreeNode*) item;

    return bstCmp(sn->path, n);
}

static void CC bstWhack(BSTNode* n, void* ignore) {
    TreeNode* sn = (TreeNode*) n;

    assert(sn);

    free(sn->path);

    memset(sn, 0, sizeof *sn);

    free(sn);
}

/********** NumIterator **********/

typedef enum {
    eNIBegin,
    eNINumber,
    eNIInterval,
    eNIDash,
    eNIComma,
    eNIBad,
    eNIForever,
    eNIEnd
} ENumIteratorState;
typedef struct {
    ENumIteratorState state;
    bool skip;
    const char *s;
    int32_t crnt;
    int32_t intEnd;
} NumIterator;
static void NumIteratorInit(NumIterator *self, const char *row) {
    assert(self);
    memset(self, 0, sizeof *self);
    self->crnt = self->intEnd = -1;
    self->s = row;
}

static int32_t NumIteratorGetNum(NumIterator *self) {
    int32_t n = 0;
    for (n = 0; *(self->s) >= '0' && *(self->s) <= '9'; ++(self->s)) {
        n = n * 10 + *(self->s) - '0';
    }
    return n;
}

static bool NumIteratorNext(NumIterator *self, int64_t crnt) {
    char c = '\0';
    assert(self);
    while (true) {
        self->skip = false;
        switch (self->state) {
            case eNIBegin:
            case eNIComma:
                if (self->s == NULL || *(self->s) == '\0') {
                    if (self->state == eNIBegin) {
                        self->state = eNIForever;
                        continue;
                    }
                    else {
                        self->state = eNIEnd;
                        continue;
                    }
                }
                c = *(self->s);
                ++(self->s);
                if (c == ',') {
                    self->state = eNIComma;
                    continue;
                }
                else if (c == '-') {
                    self->state = eNIDash;
                    continue;
                }
                else if (c >= '0' && c <= '9') {
                    --(self->s);
                    self->crnt = NumIteratorGetNum(self);
                    self->state = eNINumber;
                    if (self->crnt < crnt) {
                        continue;
                    }
                    else {
                        if (self->crnt > crnt) {
                            self->skip = true;
                        }
                        return true;
                    }
                    continue;
                }
                else {
                    self->state = eNIBad;
                    continue;
                }
            case eNIInterval:
                if (crnt <= self->intEnd) {
                    return true;
                }
          /* no break here */
            case eNINumber:
                if (self->crnt >= crnt) {
                    if (self->crnt > crnt) {
                        self->skip = true;
                    }
                    return true;
                }
                if (self->s == NULL || *(self->s) == '\0') {
                    self->state = eNIEnd;
                    continue;
                }
                c = *(self->s);
                ++(self->s);
                if (c == ',') {
                    self->state = eNIComma;
                    continue;
                }
                else if (c == '-') {
                    self->state = eNIDash;
                    continue;
                }
                else {
                    self->state = eNIBad;
                    continue;
                }
            case eNIDash:
                if (self->s == NULL || *(self->s) == '\0') {
                    self->state = eNIForever;
                    continue;
                }
                c = *(self->s);
                ++(self->s);
                if (c == ',' || c == '-') {
                    self->state = eNIForever;
                    continue;
                }
                else if (c >= '0' && c <= '9') {
                    --(self->s);
                    self->intEnd = NumIteratorGetNum(self);
                    self->state = eNIInterval;
                    if (crnt <= self->intEnd) {
                        return true;
                    }
                    else {
                        continue;
                    }
                }
                else {
                    self->state = eNIBad;
                    continue;
                }
            case eNIForever:
                return true;
            case eNIBad:
            case eNIEnd:
                return false;
        }
    }
}

/********** Resolved **********/
static rc_t ResolvedFini(Resolved *self) {
    rc_t rc = 0;

    assert(self);

    rc = VPathStrFini(&self->local);

    RELEASE(KFile, self->file);
    RELEASE(VPath, self->accession);
    RELEASE(VResolver, self->resolver);

    RELEASE(KartItem, self->kartItem);

    RELEASE(String, self->remote);
    RELEASE(String, self->cache);

    free(self->name);

    memset(self, 0, sizeof *self);

    return rc;
}

static void ResolvedReset(Resolved *self, ERunType type) {
    assert(self);

    memset(self, 0, sizeof *self);

    self->type = type;
}

/** isLocal is set to true when the object is found locally.
    i.e. does not need need not be [re]downloaded */
static rc_t ResolvedLocal(const Resolved *self,
    const KDirectory *dir, bool *isLocal, EForce force)
{
    rc_t rc = 0;
    uint64_t sRemote = 0;
    uint64_t sLocal = 0;
    const KFile *local = NULL;
    char path[PATH_MAX] = "";

    assert(isLocal && self);

    *isLocal = false;

    if (self->local.str == NULL) {
        return 0;
    }

    rc = VPathReadPath(self->local.path, path, sizeof path, NULL);
    DISP_RC(rc, "VPathReadPath");

    if (rc == 0 && KDirectoryPathType(dir, "%s", path) != kptFile) {
        if (force == eForceNo) {
            STSMSG(STS_TOP,
                ("%s (not a file) is found locally: consider it complete",
                 path));
            *isLocal = true;
        }
        else {
            STSMSG(STS_TOP,
                ("%s (not a file) is found locally and will be redownloaded",
                 path));
        }
        return 0;
    }

    if (rc == 0) {
        if (! _StringIsFasp(self->remote, NULL) && self->file != NULL) {
            rc = KFileSize(self->file, &sRemote);
            DISP_RC2(rc, "KFileSize(remote)", self->name);
        }
        else {
            sRemote = self->remoteSz;
        }
    }

    if (rc == 0) {
        rc = KDirectoryOpenFileRead(dir, &local, "%s", path);
        DISP_RC2(rc, "KDirectoryOpenFileRead", path);
    }

    if (rc == 0) {
        rc = KFileSize(local, &sLocal);
        DISP_RC2(rc, "KFileSize", path);
    }

    if (rc == 0) {
        if (sRemote == 0) {
            if (sLocal != 0) {
                if (force == eForceNo) {
                    *isLocal = true;
                    STSMSG(STS_INFO, ("%s (%,lu) is found", path, sLocal));
                }
                else {
                    STSMSG(STS_INFO,
                        ("%s (%,lu) is found and will be redownloaded",
                        path, sLocal));
                }
            }
            else if (sLocal == 0) {
                STSMSG(STS_INFO,
                    ("an empty %s (%,lu) is found and will be redownloaded",
                    path, sLocal));
            }
        }
        else if (sRemote == sLocal) {
            if (force == eForceNo) {
                *isLocal = true;
                STSMSG(STS_INFO, ("%s (%,lu) is found and is complete",
                    path, sLocal));
            }
            else {
                STSMSG(STS_INFO, ("%s (%,lu) is found and will be redownloaded",
                    path, sLocal));
            }
        }
        else {
            STSMSG(STS_TOP, ("%s (%,lu) is incomplete. Expected size is %,lu. "
                "It will be re-downloaded", path, sLocal, sRemote));
        }
    }

    RELEASE(KFile, local);

    return rc;
}

/********** Main **********/
static bool MainUseAscp(Main *self) {
    rc_t rc = 0;

    assert(self);

    if (self->ascpChecked) {
        return self->ascp != NULL;
    }

    self->ascpChecked = true;

    if (self->noAscp) {
        return false;
    }

    rc = ascp_locate(&self->ascp, &self->asperaKey, true, true);
    return rc == 0 && self->ascp && self->asperaKey;
}

static bool MainHasDownloaded(const Main *self, const char *local) {
    TreeNode *sn = NULL;

    assert(self);

    sn = (TreeNode*) BSTreeFind(&self->downloaded, local, bstCmp);

    return sn != NULL;
}

static rc_t MainDownloaded(Main *self, const char *path) {
    TreeNode *sn = NULL;

    assert(self);

    if (MainHasDownloaded(self, path)) {
        return 0;
    }

    sn = calloc(1, sizeof *sn);
    if (sn == NULL) {
        return RC(rcExe, rcStorage, rcAllocating, rcMemory, rcExhausted);
    }

    sn->path = string_dup_measure(path, NULL);
    if (sn->path == NULL) {
        bstWhack((BSTNode*) sn, NULL);
        sn = NULL;
        return RC(rcExe, rcStorage, rcAllocating, rcMemory, rcExhausted);
    }

    BSTreeInsert(&self->downloaded, (BSTNode*)sn, bstSort);

    return 0;
}

static rc_t MainDownloadFile(Resolved *self,
    Main *main, const char *to)
{
    rc_t rc = 0;
    KFile *out = NULL;
    size_t num_read = 0;
    uint64_t opos = 0;
    size_t num_writ = 0;
    uint64_t pos = 0;
    uint64_t prevPos = 0;

    assert(self && main);

    if (rc == 0) {
        STSMSG(STS_DBG, ("creating %s", to));
        rc = KDirectoryCreateFile(main->dir, &out,
                                  false, 0664, kcmInit | kcmParents, "%s", to);
        DISP_RC2(rc, "Cannot OpenFileWrite", to);
    }

    assert(self->remote);

    if (self->file == NULL) {
        rc = _KFileOpenRemote(&self->file, main->kns, self->remote->addr);
        if (rc != 0) {
            PLOGERR(klogInt, (klogInt, rc, "failed to open file for $(path)",
                "path=%s", self->remote->addr));
        }
    }

    STSMSG(STS_INFO, ("%s -> %s", self->remote->addr, to));
    do {
        bool print = pos - prevPos > 200000000;
        rc = Quitting();

        if (rc == 0) {
            if (print) {
                STSMSG(STS_FIN,
                    ("Reading %lu bytes from pos. %lu", main->bsize, pos));
            }
            rc = KFileRead(self->file,
                pos, main->buffer, main->bsize, &num_read);
            if (rc != 0) {
                DISP_RC2(rc, "Cannot KFileRead", self->remote->addr);
            }
            else {
                pos += num_read;
            }

            if (print) {
                prevPos = pos;
            }
        }

        if (rc == 0 && num_read > 0) {
            rc = KFileWrite(out, opos, main->buffer, num_read, &num_writ);
            DISP_RC2(rc, "Cannot KFileWrite", to);
            opos += num_writ;
        }
    } while (rc == 0 && num_read > 0);

    RELEASE(KFile, out);

    if (rc == 0) {
        STSMSG(STS_INFO, ("%s (%ld)", to, pos));
    }

    return rc;
}

/*  http://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByR.../SRR125365.sra
anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByR.../SRR125365.sra
*/
static rc_t MainDownloadAscp(const Resolved *self, Main *main,
    const char *to)
{
    const char *src = NULL;
    AscpOptions opt;
    assert(self && self->remote && self->remote->addr
        && main && main->ascp && main->asperaKey);
    memset(&opt, 0, sizeof opt);
    if (!_StringIsFasp(self->remote, &src)) {
        return RC(rcExe, rcFile, rcCopying, rcSchema, rcInvalid);
    }
    opt.target_rate
        = main->ascpMaxRate == NULL ? NULL : main->ascpMaxRate->addr;
    opt.name = self->name;
    opt.src_size = self->remoteSz;
    opt.heartbeat = main->heartbeat;
    opt.quitting = Quitting;
    return aspera_get(main->ascp, main->asperaKey, src, to, &opt);
}

static rc_t MainDownload(Resolved *self, Main *main) {
    bool canceled = false;
    rc_t rc = 0;
    KFile *flock = NULL;

    char tmp[PATH_MAX] = "";
    char lock[PATH_MAX] = "";

    assert(self
        && self->cache && self->cache->size && self->cache->addr && main);

    if (main->force != eForceYES &&
        MainHasDownloaded(main, self->cache->addr))
    {
        STSMSG(STS_INFO, ("%s has just been downloaded", self->cache->addr));
        return 0;
    }

    if (rc == 0) {
        rc = _KDirectoryMkLockName(main->dir, self->cache, lock, sizeof lock);
    }

    if (rc == 0) {
        rc = _KDirectoryMkTmpName(main->dir, self->cache, tmp, sizeof tmp);
    }

    if (KDirectoryPathType(main->dir, "%s", lock) != kptNotFound) {
        if (main->force != eForceYES) {
            KTime_t date = 0;
            rc = KDirectoryDate(main->dir, &date, "%s", lock);
            if (rc == 0) {
                time_t t = time(NULL) - date;
                if (t < 60 * 60 * 24) { /* 24 hours */
                    STSMSG(STS_DBG, ("%s found: canceling download", lock));
                    rc = RC(rcExe, rcFile, rcCopying, rcLock, rcExists);
                    PLOGERR(klogWarn, (klogWarn, rc,
                        "Lock file $(file) exists: download canceled",
                        "file=%s", lock));
                    return rc;
                }
                else {
                    STSMSG(STS_DBG, ("%s found and ignored as too old", lock));
                    rc = _KDirectoryClean(main->dir,
                        self->cache, NULL, NULL, true);
                }
            }
            else {
                STSMSG(STS_DBG, ("%s found", lock));
                DISP_RC2(rc, "KDirectoryDate", lock);
            }
        }
        else {
            STSMSG(STS_DBG, ("%s found and forced to be ignored", lock));
            rc = _KDirectoryClean(main->dir, self->cache, NULL, NULL, true);
        }
    }
    else {
        STSMSG(STS_DBG, ("%s not found", lock));
    }

    if (rc == 0) {
        STSMSG(STS_DBG, ("creating %s", lock));
        rc = KDirectoryCreateFile(main->dir, &flock,
            false, 0664, kcmInit | kcmParents, "%s", lock);
        DISP_RC2(rc, "Cannot OpenFileWrite", lock);
    }

    assert(!main->noAscp || !main->noHttp);

    if (rc == 0) {
        bool ascp = _StringIsFasp(self->remote, NULL);
        if (ascp) {
            STSMSG(STS_TOP, (" Downloading via fasp..."));
            if (main->forceAscpFail) {
                rc = 1;
            }
            else {
                rc = MainDownloadAscp(self, main, tmp);
            }
            if (rc == 0) {
                STSMSG(STS_TOP, (" fasp download succeed"));
            }
            else {
                rc_t rc = Quitting();
                if (rc != 0) {
                    canceled = true;
                }
                else {
                    STSMSG(STS_TOP, (" fasp download failed"));
                }
            }
        }
        if (!ascp || (rc != 0 && GetRCObject(rc) != rcMemory
                              && !canceled && !main->noHttp))
        {
            const VPath *vremote = NULL;
            STSMSG(STS_TOP, (" Downloading via http..."));
            if (ascp) {
                assert(self->resolver);
                RELEASE(String, self->remote);
                RELEASE(KFile, self->file);
                rc = _VResolverRemote(self->resolver, eProtocolHttp, self->name,
                    self->accession, &vremote, &self->remote, &self->cache);
            }
            if (rc == 0) {
                rc = MainDownloadFile(self, main, tmp);
            }
            RELEASE(VPath, vremote);
        }
    }

    RELEASE(KFile, flock);

    if (rc == 0) {
        STSMSG(STS_DBG, ("renaming %s -> %s", tmp, self->cache->addr));
        rc = KDirectoryRename(main->dir, true, tmp, self->cache->addr);
        if (rc != 0) {
            PLOGERR(klogInt, (klogInt, rc, "cannot rename $(from) to $(to)",
                "from=%s,to=%s", tmp, self->cache->addr));
        }
    }

    if (rc == 0) {
        rc = MainDownloaded(main, self->cache->addr);
    }

    if (rc == 0) {
        rc_t rc2 = _KDirectoryCleanCache(main->dir, self->cache);
        if (rc == 0 && rc2 != 0) {
            rc = rc2;
        }
    }

    {
        rc_t rc2 = _KDirectoryClean(main->dir, self->cache, lock, tmp, rc != 0);
        if (rc == 0 && rc2 != 0) {
            rc = rc2;
        }
    }

    return rc;
}

static rc_t MainDependenciesList(const Main *self,
    const char *path, const VDBDependencies **deps)
{
    rc_t rc = 0;
    bool isDb = true;
    const VDatabase *db = NULL;

    assert(self && path && deps);

    if ((VDBManagerPathType(self->mgr, "%s", path) & ~kptAlias) != kptDatabase) {
        return 0;
    }

    rc = VDBManagerOpenDBRead(self->mgr, &db, NULL, "%s", path);
    if (rc != 0) {
        if (rc == SILENT_RC(rcDB, rcMgr, rcOpening, rcDatabase, rcIncorrect)) {
            isDb = false;
            rc = 0;
        }
        else if (rc ==
            SILENT_RC(rcKFG, rcEncryptionKey, rcRetrieving, rcItem, rcNotFound))
        {
            STSMSG(STS_TOP, ("Cannot open encrypted file '%s'", path));
            isDb = false;
            rc = 0;
        }
        DISP_RC2(rc, "Cannot open database", path);
    }

    if (rc == 0 && isDb) {
        bool all = self->check_all || self->force != eForceNo;
        rc = VDatabaseListDependencies(db, deps, !all);
        DISP_RC2(rc, "VDatabaseListDependencies", path);
    }

    RELEASE(VDatabase, db);

    return rc;
}

/********** Item **********/
static rc_t ItemRelease(Item *self) {
    rc_t rc = 0;

    assert(self);

    rc = ResolvedFini(&self->resolved);
    RELEASE(KartItem, self->item);

    memset(self, 0, sizeof *self);

    return rc;
}

static rc_t ItemInit(Item *self, const char *obj) {
    assert(self);
    self->desc = obj;
    return 0;
}

static char* ItemName(const Item *self) {
    char *c = NULL;
    assert(self);
    if (self->desc != NULL) {
        return string_dup_measure(self->desc, NULL);
    }
    else {
        rc_t rc = 0;
        const String *elem = NULL;
        assert(self->item);

        rc = KartItemItemDesc(self->item, &elem);
        c = StringCheck(elem, rc);
        if (c != NULL) {
            return c;
        }

        rc = KartItemAccession(self->item, &elem);
        c = StringCheck(elem, rc);
        if (c != NULL) {
            return c;
        }

        rc = KartItemItemId(self->item, &elem);
        return StringCheck(elem, rc);
    }
}

static rc_t _ItemSetResolverAndAssessionInResolved(Item *item,
    VResolver *resolver, const KConfig *cfg, const KRepositoryMgr *repoMgr,
    const VFSManager *vfs)
{
    Resolved *resolved = NULL;
    rc_t rc = 0;

    assert(item && resolver && cfg && repoMgr && vfs);

    resolved = &item->resolved;

    if (item->desc != NULL) {
        rc = VFSManagerMakePath(vfs, &resolved->accession, "%s", item->desc);
        DISP_RC2(rc, "VFSManagerMakePath", item->desc);
        if (rc == 0) {
            rc = VResolverAddRef(resolver);
        }
        if (rc == 0) {
            resolved->resolver = resolver;
        }
    }
    else {
        uint64_t oid = 0;
        rc = KartItemProjIdNumber(item->item, &resolved->project);
        if (rc != 0) {
            DISP_RC(rc, "KartItemProjIdNumber");
            return rc;
        }
        rc = KartItemItemIdNumber(item->item, &oid);
        if (rc != 0) {
            DISP_RC(rc, "KartItemItemIdNumber");
            return rc;
        }
        else {
            const KRepository *p_protected = NULL;
            rc = KRepositoryMgrGetProtectedRepository(repoMgr, 
                (uint32_t)resolved->project, &p_protected);
            if (rc == 0) {
                rc = KRepositoryMakeResolver(p_protected,
                    &resolved->resolver, cfg);
                if (rc != 0) {
                    DISP_RC(rc, "KRepositoryMakeResolver");
                    return rc;
                }
            }
            else {
                PLOGERR(klogErr, (klogErr, rc,
                    "project '$(P)': cannot find protected repository", "P=%d",
                    resolved->project));
            }
            RELEASE(KRepository, p_protected);
        }
        if (rc == 0) {
            rc =
                VFSManagerMakeOidPath(vfs, &resolved->accession, (uint32_t)oid);
            if (rc != 0) {
                DISP_RC(rc, "VFSManagerMakeOidPath");
                return rc;
            }
        }
    }

    return rc;
}

/* resolve locations */
static rc_t _ItemResolveResolved(VResolver *resolver,
    VRemoteProtocols protocols, Item *item, const KRepositoryMgr *repoMgr,
    const KConfig *cfg, const VFSManager *vfs, KNSManager *kns, size_t minSize, size_t maxSize)
{
    Resolved *resolved = NULL;
    rc_t rc = 0;
    rc_t rc2 = 0;

    const VPath *vremote = NULL;

    assert(resolver && item);

    resolved = &item->resolved;

    memset(&resolved->local, 0, sizeof resolved->local);

    assert(resolved->accession == NULL);

    rc = _ItemSetResolverAndAssessionInResolved(item,
        resolver, cfg, repoMgr, vfs);

    if (rc == 0) {
        rc = V_ResolverLocal(resolved->resolver,
            resolved->accession, &resolved->local.path);
        if (rc == 0) {
            rc = VPathMakeString(resolved->local.path, &resolved->local.str);
            DISP_RC2(rc, "VPathMakeString(VResolverLocal)", resolved->name);
        }
        else if (NotFoundByResolver(rc)) {
            rc = 0;
        }
        else {
            DISP_RC2(rc, "VResolverLocal", resolved->name);
        }
    }

    if (rc == 0) {
        rc2 = 0;
        resolved->remoteSz = 0;
        assert(item->main);
        if ((minSize > 0 || maxSize > 0 || item->main->order == eOrderSize)
            && (protocols == eProtocolFasp ||
                protocols == eProtocolFaspHttp))
        {
            rc2 = _VResolverRemote(resolved->resolver, eProtocolHttp,
                resolved->name, resolved->accession,
                &vremote, &resolved->remote, NULL);
            if (rc2 != 0 && rc == 0) {
                rc = rc2;
            }
            else {
                rc_t rc3 = 0;
                if (resolved->file == NULL) {
                    rc3 = _KFileOpenRemote(&resolved->file,
                                           kns, resolved->remote->addr);
                    DISP_RC2(rc3,
                        "cannot open remote file", resolved->remote->addr);
                }

                RELEASE(VPath, vremote);

                if (rc3 == 0 && resolved->file != NULL) {
                    rc3 = KFileSize(resolved->file, &resolved->remoteSz);
                    if (rc3 != 0) {
                        DISP_RC2(rc3, "cannot get remote file size",
                            resolved->remote->addr);
                    }
                    else if (resolved->remoteSz >= maxSize) {
                        return rc;
                    }
                    else if (resolved->remoteSz < minSize) {
                        return rc;
                    }
                }
            }
        }

        if (rc2 == 0) {
            rc2 = _VResolverRemote(resolved->resolver, protocols,
                resolved->name, resolved->accession,
                &vremote, &resolved->remote, &resolved->cache);
            if (rc2 != 0 && rc == 0) {
                rc = rc2;
            }
        }

        if (rc == 0) {
            rc2 = 0;
            if (resolved->file == NULL) {
                assert(resolved->remote);
                if (!_StringIsFasp(resolved->remote, NULL)) {
                    rc2 = _KFileOpenRemote(
                        &resolved->file, kns, resolved->remote->addr);
                }
            }
            if (rc2 == 0 && resolved->file != NULL && resolved->remoteSz == 0) {
                rc2 = KFileSize(resolved->file, &resolved->remoteSz);
                DISP_RC2(rc2, "KFileSize(remote)", resolved->name);
            }
        }
    }

    RELEASE(VPath, vremote);
    return rc;
}

/* Resolved: resolve locations */
static rc_t ItemInitResolved(Item *self, VResolver *resolver,
    KDirectory *dir, bool ascp, const KRepositoryMgr *repoMgr,
    const KConfig *cfg, const VFSManager *vfs, KNSManager *kns, size_t minSize, size_t maxSize)
{
    Resolved *resolved = NULL;
    rc_t rc = 0;
    KPathType type = kptNotFound;
    VRemoteProtocols protocols = ascp ? eProtocolFaspHttp : eProtocolHttp;

    assert(self);

    resolved = &self->resolved;
    resolved->name = ItemName(self);
    
    assert(resolved->type != eRunTypeUnknown);

    type = KDirectoryPathType(dir, "%s", self->desc) & ~kptAlias;
    if (type == kptFile || type == kptDir) {
        resolved->path = self->desc;
        resolved->existing = true;
        if (resolved->type != eRunTypeDownload) {
            uint64_t s = -1;
            const KFile *f = NULL;
            rc = KDirectoryOpenFileRead(dir, &f, "%s", self->desc);
            if (rc == 0) {
                rc = KFileSize(f, &s);
            }
            if (s != -1) {
/*              OUTMSG(("%s\t%,zuB\n", self->desc, s)); */
                resolved->remoteSz = s;
            }
            else {
                OUTMSG(("%s\tunknown\n", self->desc));
            }
            RELEASE(KFile, f);
        }
        else {
            STSMSG(STS_TOP, ("'%s' is a local non-kart file", self->desc));
        }
        return 0;
    }

    rc = _ItemResolveResolved(resolver, protocols, self,
        repoMgr, cfg, vfs, kns, minSize, maxSize);

    STSMSG(STS_DBG, ("Resolve(%s) = %R:", resolved->name, rc));
    STSMSG(STS_DBG, ("local(%s)",
        resolved->local.str ? resolved->local.str->addr : "NULL"));
    STSMSG(STS_DBG, ("cache(%s)",
        resolved->cache ? resolved->cache->addr : "NULL"));
    STSMSG(STS_DBG, ("remote(%s:%,ld)",
        resolved->remote ? resolved->remote->addr : "NULL",
        resolved->remoteSz));

    if (rc == 0) {
        if (resolved->remoteSz >= maxSize) {
            resolved->oversized = true;
            return rc;
        }
        if (resolved->remoteSz < minSize) {
            resolved->undersized = true;
            return rc;
        }

        if (resolved->local.str == NULL
            && (resolved->cache == NULL || resolved->remote == NULL))
        {
            rc = RC(rcExe, rcPath, rcValidating, rcParam, rcNull);
            PLOGERR(klogInt, (klogInt, rc,
                "bad VResolverResolve($(acc)) result",
                "acc=%s", resolved->name));
        }
    }

    return rc;
}

/* resolve: locate */
static rc_t ItemResolve(Item *item, int32_t row) {
    Resolved *self = NULL;
    static int n = 0;
    rc_t rc = 0;
    bool ascp = false;

    assert(item && item->main);

    self = &item->resolved;
    assert(self->type);

    ++n;
    if (row > 0 && item->desc == NULL) {
        n = row;
    }

    item->number = n;

    ascp = MainUseAscp(item->main);
    if (self->type == eRunTypeList) {
        ascp = false;
    }

    rc = ItemInitResolved(item, item->main->resolver, item->main->dir, ascp,
        item->main->repoMgr, item->main->cfg, item->main->vfsMgr,
        item->main->kns, item->main->minSize, item->main->maxSize);

    return rc;
}

/* download if not found; obey size restriction */
static rc_t ItemDownload(Item *item) {
    bool isLocal = false;
    int n = 0;
    rc_t rc = 0;
    Resolved *self = NULL;
    assert(item && item->main);
    n = item->number;
    self = &item->resolved;
    assert(self->type);

    if (rc == 0) {
        bool skip = false;

        if (self->existing) {
            self->path = item->desc;
            return rc;
        }

        if (self->undersized) {
            STSMSG(STS_TOP,
               ("%d) '%s' (%,zu KB) is smaller than minimum allowed: skipped\n",
                n, self->name, self->remoteSz / 1024));
            skip = true;
        }
        else if (self->oversized) {
            STSMSG(STS_TOP,
                ("%d) '%s' (%,zu KB) is larger than maximum allowed: skipped\n",
                n, self->name, self->remoteSz / 1024));
            skip = true;
        }

        rc = ResolvedLocal(self, item->main->dir, &isLocal,
            skip ? eForceNo : item->main->force);

        if (rc == 0) {
            if (skip && !isLocal) {
                return rc;
            }
        }
    }

    if (rc == 0) {
        if (isLocal) {
            STSMSG(STS_TOP, ("%d) '%s' is found locally", n, self->name));
            if (self->local.str != NULL) {
                self->path = self->local.str->addr;
            }
        }
        else if (!_StringIsFasp(self->remote, NULL)
            && item->main->noHttp)
        {
            rc = RC(rcExe, rcFile, rcCopying, rcFile, rcNotFound);
            PLOGERR(klogErr, (klogErr, rc,
                "cannot download '$(name)' using requested transport",
                "name=%s", self->name));
        }
        else {
            STSMSG(STS_TOP, ("%d) Downloading '%s'...", n, self->name));
            rc = MainDownload(self, item->main);
            if (rc == 0) {
                STSMSG(STS_TOP,
                    ("%d) '%s' was downloaded successfully", n, self->name));
                if (self->cache != NULL) {
                    self->path = self->cache->addr;
                }
            }
            else if (rc != SILENT_RC(rcExe,
                rcProcess, rcExecuting, rcProcess, rcCanceled))
            {
                STSMSG(STS_TOP, ("%d) failed to download %s", n, self->name));
            }
        }
    }
    else {
        STSMSG(STS_TOP, ("%d) cannot locate '%s'", n, self->name));
    }

    return rc;
}

static rc_t ItemPrintSized(const Item *self, int32_t row, size_t size) {
    rc_t rc = 0;
    const String *projId = NULL;
    const String *itemId = NULL;
    const String *accession = NULL;
    const String *name = NULL;
    const String *itemDesc = NULL;
    assert(self && row);
    if (self->desc != NULL) {
        OUTMSG(("%d\t%s\t%,zuB\n", row, self->desc, size));
    }
    else {
        if (rc == 0) {
            rc = KartItemProjId(self->item, &projId);
        }
        if (rc == 0) {
            rc = KartItemItemId(self->item, &itemId);
        }
        if (rc == 0) {
            rc = KartItemAccession(self->item, &accession);
        }
        if (rc == 0) {
            rc = KartItemName(self->item, &name);
        }
        if (rc == 0) {
            rc = KartItemItemDesc(self->item, &itemDesc);
        }
        if (rc == 0) {
           rc = OUTMSG(("%d\t%S|%S|%S|%S|%S\t%,zuB\n",
               row, projId, itemId, accession, name, itemDesc, size));
        }
    }

    return rc;
}

/* resolve: locate; download if not found */
static rc_t ItemResolveResolvedAndDownloadOrProcess(Item *self, int32_t row) {
    rc_t rc = ItemResolve(self, row);
    if (rc != 0) {
        return rc;
    }

    assert(self);

    if (self->resolved.type == eRunTypeList) {
        return ItemPrintSized(self, row, self->resolved.remoteSz);
    }
    else if (self->resolved.type != eRunTypeDownload) {
        return rc;
    }

    return ItemDownload(self);
}

static rc_t ItemPostDownload(Item *item, int32_t row) {
    Resolved *resolved = NULL;
    rc_t rc = 0;
    rc_t rc2 = 0;
    const VDBDependencies *deps = NULL;
    uint32_t count = 0;
    uint32_t i = 0;

    assert(item && item->main);

    resolved = &item->resolved;

    if (rc == 0) {
        if (resolved->type == eRunTypeList) {
            return rc;
        }
        else if (resolved->oversized) {
            item->main->oversized = true;
        }
        else if (resolved->undersized) {
            item->main->undersized = true;
        }
        if (resolved->path != NULL) {
            rc = MainDependenciesList(item->main, resolved->path, &deps);
        }
    }

    /* resolve dependencies (refseqs) */
    if (rc == 0 && deps != NULL) {
        rc = VDBDependenciesCount(deps, &count);
        if (rc == 0) {
            STSMSG(STS_TOP, ("'%s' has %d%s dependenc%s",
                item->desc, count, item->main->check_all ? "" : " unresolved",
                count == 1 ? "y" : "ies"));
        }
        else {
            DISP_RC2(rc, "Failed to check %s's dependencies", item->desc);
        }
    }

    for (i = 0; i < count && rc == 0; ++i) {
        bool local = true;
        const char *seq_id = NULL;

        if (rc == 0) {
            rc = VDBDependenciesLocal(deps, &local, i);
            DISP_RC2(rc, "VDBDependenciesLocal", item->desc);
            if (local) {
                continue;
            }
        }

        if (rc == 0) {
            rc = VDBDependenciesSeqId(deps, &seq_id, i);
            DISP_RC2(rc, "VDBDependenciesSeqId", item->desc);
        }

        if (rc == 0) {
            size_t num_writ = 0;
            char ncbiAcc[512] = "";

            assert(seq_id);

            rc = string_printf(ncbiAcc, sizeof ncbiAcc, &num_writ,
                "ncbi-acc:%s?vdb-ctx=refseq", seq_id);
            DISP_RC2(rc, "string_printf(?vdb-ctx=refseq)", seq_id);
            if (rc == 0 && num_writ > sizeof ncbiAcc) {
                rc = RC(rcExe, rcFile, rcCopying, rcBuffer, rcInsufficient);
                PLOGERR(klogInt, (klogInt, rc,
                    "bad string_printf($(s)?vdb-ctx=refseq) result",
                    "s=%s", seq_id));
            }
    
            if (rc == 0) {
                Item *ditem = calloc(1, sizeof *ditem);
                if (ditem == NULL) {
                    return RC(rcExe,
                        rcStorage, rcAllocating, rcMemory, rcExhausted);
                }

                ditem->desc = ncbiAcc;
                ditem->main = item->main;

                ResolvedReset(&ditem->resolved, eRunTypeDownload);

                rc = ItemResolveResolvedAndDownloadOrProcess(ditem, 0);

                RELEASE(Item, ditem);
            }
        }
    }

    if (rc == 0 && rc2 != 0) {
        rc = rc2;
    }

    RELEASE(VDBDependencies, deps);

    return rc;
}

static rc_t ItemProcess(Item *item, int32_t row) {
    rc_t rc = 0;

    assert(item);

    /* resolve: locate; download if not found */
    rc = ItemResolveResolvedAndDownloadOrProcess(item, row);

    if (item->resolved.type != eRunTypeDownload) {
        return rc;
    }

    if (rc == 0) {
        rc = ItemPostDownload(item, row);
    }

    return rc;
}

/*********** Iterator **********/
static
rc_t IteratorInit(Iterator *self, const char *obj, const Main *main)
{
    rc_t rc = 0;

    KPathType type = kptNotFound;

    assert(self && main);
    memset(self, 0, sizeof *self);

#ifdef _DEBUGGING
    if (obj == NULL && main->textkart) {
        type = KDirectoryPathType(main->dir, "%s", main->textkart);
        if ((type & ~kptAlias) != kptFile) {
            rc = RC(rcExe, rcFile, rcOpening, rcFile, rcNotFound);
            DISP_RC(rc, main->textkart);
            return rc;
        }
        rc = KartMakeText(main->dir, main->textkart, &self->kart,
            &self->isKart);
        if (rc != 0) {
            if (!self->isKart) {
                rc = 0;
            }
            else {
                PLOGERR(klogErr, (klogErr, rc, "'$(F)' is not a text kart file",
                    "F=%s", main->textkart));
            }
        }
        return rc;
    }
#endif

    assert(obj);
    type = KDirectoryPathType(main->dir, "%s", obj);
    if ((type & ~kptAlias) == kptFile) {
        type = VDBManagerPathType(main->mgr, "%s", obj);
        if ((type & ~kptAlias) == kptFile) {
            rc = KartMake(main->dir, obj, &self->kart, &self->isKart);
            if (!self->isKart) {
                rc = 0;
            }
        }
    }

    if (rc == 0 && !self->isKart) {
        self->obj = obj;
    }

    return rc;
}

static rc_t IteratorNext(Iterator *self, Item **next, bool *done) {
    rc_t rc = 0;

    assert(self && next && done);

    *next = NULL;

    if (self->done) {
        *done = true;
        return 0;
    }

    *done = false;

    *next = calloc(1, sizeof **next);
    if (*next == NULL) {
        return RC(rcExe, rcStorage, rcAllocating, rcMemory, rcExhausted);
    }

    if (self->isKart) {
        rc = KartMakeNextItem(self->kart, &(*next)->item);
        if (rc != 0) {
            LOGERR(klogErr, rc, "Invalid kart file: cannot read next row");
        }
        else if ((*next)->item == NULL) {
            RELEASE(Item, *next);
            *next = NULL;
            *done = true;
        }

        if (rc == 0 && *done) {
            self->done = true;
        }
    }
    else {
        rc = ItemInit(*next, self->obj);

        self->done = true;
    }

    return rc;
}

static void IteratorFini(Iterator *self) {
    rc_t rc = 0;

    assert(self);

    RELEASE(Kart, self->kart);
}

/*********** Command line arguments **********/

static size_t _sizeFromString(const char *val) {
    size_t s = 0;

    for (s = 0; *val != '\0'; ++val) {
        if (*val < '0' || *val > '9') {
            break;
        }
        s = s * 10 + *val - '0';
    }

    if (*val == '\0' || *val == 'k' || *val == 'K') {
        s *= 1024L;
    }
    else if (*val == 'b' || *val == 'B') {
    }
    else if (*val == 'm' || *val == 'M') {
        s *= 1024L * 1024;
    }
    else if (*val == 'g' || *val == 'G') {
        s *= 1024L * 1024 * 1024;
    }
    else if (*val == 'u' || *val == 'U') {  /* unlimited */
        s = 0;
    }

    return s;
}

#define ASCP_OPTION "ascp-path"
#define ASCP_ALIAS  "a"
static const char* ASCP_USAGE[] =
{ "path to ascp program and private key file (asperaweb_id_dsa.putty)", NULL };

#define CHECK_ALL_OPTION "check-all"
#define CHECK_ALL_ALIAS  "c"
static const char* CHECK_ALL_USAGE[] = { "double-check all refseqs", NULL };

#define FORCE_OPTION "force"
#define FORCE_ALIAS  "f"
static const char* FORCE_USAGE[] = {
    "force object download - one of: no, yes, all.",
    "no [default]: skip download if the object if found and complete;",
    "yes: download it even if it is found and is complete;", "all: ignore lock "
    "files (stale locks or it is beeing downloaded by another process: "
    "use at your own risk!)", NULL };

#define FAIL_ASCP_OPTION "FAIL-ASCP"
#define FAIL_ASCP_ALIAS  "F"
static const char* FAIL_ASCP_USAGE[] = {
    "force ascp download fail to test ascp->http download combination" };

#define LIST_OPTION "list"
#define LIST_ALIAS  "l"
static const char* LIST_USAGE[] = { "list the content of a kart file", NULL };

#define NM_L_OPTION "numbered-list"
#define NM_L_ALIAS  "n"
static const char* NM_L_USAGE[] =
{ "list the content of a kart file with kart row numbers", NULL };

#define MINSZ_OPTION "min-size"
#define MINSZ_ALIAS  "N"
static const char* MINSZ_USAGE[] =
{ "minimum file size to download in KB (inclusive).", NULL };

#define ORDR_OPTION "order"
#define ORDR_ALIAS  "o"
static const char* ORDR_USAGE[] = { "kart prefetch order: one of: kart, size.",
    "(in kart order, by file size: smallest first), default: size", NULL };

#define HBEAT_OPTION "progress"
#define HBEAT_ALIAS  "p"
static const char* HBEAT_USAGE[] = {
    "time period in minutes to display download progress",
    "(0: no progress), default: 1", NULL };

#define ROWS_OPTION "rows"
#define ROWS_ALIAS  "R"
static const char* ROWS_USAGE[] =
{ "kart rows (default all).", "row list should be ordered", NULL };

#define SZ_L_OPTION "list-sizes"
#define SZ_L_ALIAS  "s"
static const char* SZ_L_USAGE[] =
{ "list the content of a kart file with target file sizes", NULL };

#define TRANS_OPTION "transport"
#define TRASN_ALIAS  "t"
static const char* TRANS_USAGE[] = { "transport: one of: ascp; http; both.",
    "(ascp only; http only; first try ascp, "
    "use http if cannot download by ascp).",
    "Default: both", NULL };

#define DEFAULT_MAX_FILE_SIZE "20G"
#define SIZE_OPTION "max-size"
#define SIZE_ALIAS  "X"
static const char* SIZE_USAGE[] = {
    "maximum file size to download in KB (exclusive).",
    "Default: " DEFAULT_MAX_FILE_SIZE, NULL };

#ifdef _DEBUGGING
#define TEXTKART_OPTION "text-kart"
static const char* TEXTKART_USAGE[] =
{ "To read a textual format kart file (DEBUG ONLY)", NULL };
#endif

static OptDef Options[] = {
    /*                                                    needs_value required*/
    { FORCE_OPTION    , FORCE_ALIAS    , NULL, FORCE_USAGE , 1, true, false }
   ,{ TRANS_OPTION    , TRASN_ALIAS    , NULL, TRANS_USAGE , 1, true, false }
   ,{ LIST_OPTION     , LIST_ALIAS     , NULL, LIST_USAGE  , 1, false,false }
   ,{ NM_L_OPTION     , NM_L_ALIAS     , NULL, NM_L_USAGE  , 1, false,false }
   ,{ SZ_L_OPTION     , SZ_L_ALIAS     , NULL, SZ_L_USAGE  , 1, false,false }
   ,{ ROWS_OPTION     , ROWS_ALIAS     , NULL, ROWS_USAGE  , 1, true, false }
   ,{ MINSZ_OPTION    , MINSZ_ALIAS    , NULL, MINSZ_USAGE , 1, true ,false }
   ,{ SIZE_OPTION     , SIZE_ALIAS     , NULL, SIZE_USAGE  , 1, true ,false }
   ,{ ORDR_OPTION     , ORDR_ALIAS     , NULL, ORDR_USAGE  , 1, true ,false }
   ,{ ASCP_OPTION     , ASCP_ALIAS     , NULL, ASCP_USAGE  , 1, true ,false }
   ,{ HBEAT_OPTION    , HBEAT_ALIAS    , NULL, HBEAT_USAGE , 1, true, false }
   ,{ FAIL_ASCP_OPTION, FAIL_ASCP_ALIAS, NULL, FAIL_ASCP_USAGE, 1, false, false}
#ifdef _DEBUGGING
   ,{ TEXTKART_OPTION , NULL           , NULL, TEXTKART_USAGE , 1, true , false}
#endif
   ,{ CHECK_ALL_OPTION, CHECK_ALL_ALIAS, NULL, CHECK_ALL_USAGE, 1, false, false}
};

static rc_t MainProcessArgs(Main *self, int argc, char *argv[]) {
    rc_t rc = 0;

    uint32_t pcount = 0;

    assert(self);

    rc = ArgsMakeAndHandle(&self->args, argc, argv, 1,
        Options, sizeof Options / sizeof (OptDef));
    if (rc != 0) {
        DISP_RC(rc, "ArgsMakeAndHandle");
        return rc;
    }

    do {
/* FORCE_OPTION goes first */
        rc = ArgsOptionCount (self->args, FORCE_OPTION, &pcount);
        if (rc != 0) {
            LOGERR(klogErr, rc, "Failure to get '" FORCE_OPTION "' argument");
            break;
        }

        if (pcount > 0) {
            const char *val = NULL;
            rc = ArgsOptionValue(self->args, FORCE_OPTION, 0, &val);
            if (rc != 0) {
                LOGERR(klogErr, rc,
                    "Failure to get '" FORCE_OPTION "' argument value");
                break;
            }
            if (val == NULL || val[0] == '\0') {
                rc = RC(rcExe, rcArgv, rcParsing, rcParam, rcInvalid);
                LOGERR(klogErr, rc,
                    "Unrecognized '" FORCE_OPTION "' argument value");
                break;
            }
            switch (val[0]) {
                case 'n':
                case 'N':
                    self->force = eForceNo;
                    break;
                case 'y':
                case 'Y':
                    self->force = eForceYes;
                    break;
                case 'a':
                case 'A':
                    self->force = eForceYES;
                    break;
                default:
                    rc = RC(rcExe, rcArgv, rcParsing, rcParam, rcInvalid);
                    LOGERR(klogErr, rc,
                        "Unrecognized '" FORCE_OPTION "' argument value");
                    break;
            }
            if (rc != 0) {
                break;
            }
        }

/* CHECK_ALL_OPTION goes after FORCE_OPTION */
        rc = ArgsOptionCount(self->args, CHECK_ALL_OPTION, &pcount);
        if (rc != 0) {
            LOGERR(klogErr,
                rc, "Failure to get '" CHECK_ALL_OPTION "' argument");
            break;
        }
        if (pcount > 0 || self->force != eForceNo) {
            self->check_all = true;
        }

/******* LIST OPTIONS BEGIN ********/
/* LIST_OPTION */
        rc = ArgsOptionCount(self->args, LIST_OPTION, &pcount);
        if (rc != 0) {
            LOGERR(klogErr,
                rc, "Failure to get '" LIST_OPTION "' argument");
            break;
        }
        if (pcount > 0) {
            self->list_kart = true;
        }

/* NM_L_OPTION */
        rc = ArgsOptionCount(self->args, NM_L_OPTION, &pcount);
        if (rc != 0) {
            LOGERR(klogErr,
                rc, "Failure to get '" NM_L_OPTION "' argument");
            break;
        }
        if (pcount > 0) {
            self->list_kart = self->list_kart_numbered = true;
        }

/* SZ_L_OPTION */
        rc = ArgsOptionCount(self->args, SZ_L_OPTION, &pcount);
        if (rc != 0) {
            LOGERR(klogErr,
                rc, "Failure to get '" SZ_L_OPTION "' argument");
            break;
        }
        if (pcount > 0) { /* self->list_kart is not set here! */
            self->list_kart_sized = true;
        }
/******* LIST OPTIONS END ********/

/* ASCP_OPTION */
        rc = ArgsOptionCount(self->args, ASCP_OPTION, &pcount);
        if (rc != 0) {
            LOGERR(klogErr,
                rc, "Failure to get '" ASCP_OPTION "' argument");
            break;
        }
        if (pcount > 0) {
            const char *val = NULL;
            rc = ArgsOptionValue(self->args, ASCP_OPTION, 0, &val);
            if (rc != 0) {
                LOGERR(klogErr, rc,
                    "Failure to get '" ASCP_OPTION "' argument value");
                break;
            }
            if (val != NULL) {
                char *sep = strchr(val, '|');
                if (sep == NULL) {
                    rc = RC(rcExe, rcArgv, rcParsing, rcParam, rcInvalid);
                    LOGERR(klogErr, rc,
             "ascp-path expected in the following format:\n"
             "--" ASCP_OPTION " \"<ascp-binary|private-key-file>\"\n"
             "Examples:\n"
             "--" ASCP_OPTION " \"/usr/bin/ascp|/etc/asperaweb_id_dsa.putty\"\n"
             "--" ASCP_OPTION " \"C:\\Program Files\\Aspera\\ascp.exe|C:\\Program Files\\Aspera\\etc\\asperaweb_id_dsa.putty\"\n");
                    break;
                }
                else {
                    self->ascp = string_dup(val, sep - val);
                    self->asperaKey = string_dup_measure(sep + 1, NULL);
                    self->ascpChecked = true;
                }
            }
        }

/* FAIL_ASCP_OPTION */
        rc = ArgsOptionCount(self->args, FAIL_ASCP_OPTION, &pcount);
        if (rc != 0) {
            LOGERR(klogErr,
                rc, "Failure to get '" FAIL_ASCP_OPTION "' argument");
            break;
        }
        if (pcount > 0) {
            self->forceAscpFail = true;
        }

/* HBEAT_OPTION */
        rc = ArgsOptionCount(self->args, HBEAT_OPTION, &pcount);
        if (rc != 0) {
            LOGERR(klogErr, rc, "Failure to get '" HBEAT_OPTION "' argument");
            break;
        }

        if (pcount > 0) {
            double f;
            const char *val = NULL;
            rc = ArgsOptionValue(self->args, HBEAT_OPTION, 0, &val);
            if (rc != 0) {
                LOGERR(klogErr, rc,
                    "Failure to get '" HBEAT_OPTION "' argument value");
                break;
            }
            f = atof(val) * 60000;
            self->heartbeat = (uint64_t)f;
        }

/* ORDR_OPTION */
        rc = ArgsOptionCount(self->args, ORDR_OPTION, &pcount);
        if (rc != 0) {
            LOGERR(klogErr, rc, "Failure to get '" ORDR_OPTION "' argument");
            break;
        }

        if (pcount > 0) {
            const char *val = NULL;
            rc = ArgsOptionValue(self->args, ORDR_OPTION, 0, &val);
            if (rc != 0) {
                LOGERR(klogErr, rc,
                    "Failure to get '" ORDR_OPTION "' argument value");
                break;
            }
            if (val != NULL && val[0] == 's') {
                self->order = eOrderSize;
            }
            else {
                self->order = eOrderOrig;
            }
        }

/* ROWS_OPTION */
        rc = ArgsOptionCount(self->args, ROWS_OPTION, &pcount);
        if (rc != 0) {
            LOGERR(klogErr,
                rc, "Failure to get '" ROWS_OPTION "' argument");
            break;
        }
        if (pcount > 0) {
            rc = ArgsOptionValue(self->args, ROWS_OPTION, 0, &self->rows);
            if (rc != 0) {
                LOGERR(klogErr, rc,
                    "Failure to get '" ROWS_OPTION "' argument value");
                break;
            }
        }

/* MINSZ_OPTION */
        {
            const char *val = "0";
            rc = ArgsOptionCount(self->args, MINSZ_OPTION, &pcount);
            if (rc != 0) {
                LOGERR(klogErr,
                    rc, "Failure to get '" MINSZ_OPTION "' argument");
                break;
            }
            if (pcount > 0) {
                rc = ArgsOptionValue(self->args, MINSZ_OPTION, 0, &val);
                if (rc != 0) {
                    LOGERR(klogErr, rc,
                        "Failure to get '" MINSZ_OPTION "' argument value");
                    break;
                }
            }
            self->minSize = _sizeFromString(val);
        }

/* SIZE_OPTION */
        {
            const char *val = DEFAULT_MAX_FILE_SIZE;
            rc = ArgsOptionCount(self->args, SIZE_OPTION, &pcount);
            if (rc != 0) {
                LOGERR(klogErr,
                    rc, "Failure to get '" SIZE_OPTION "' argument");
                break;
            }
            if (pcount > 0) {
                rc = ArgsOptionValue(self->args, SIZE_OPTION, 0, &val);
                if (rc != 0) {
                    LOGERR(klogErr, rc,
                        "Failure to get '" SIZE_OPTION "' argument value");
                    break;
                }
            }
            self->maxSize = _sizeFromString(val);
            if (self->maxSize == 0) {
                rc = RC(rcExe, rcArgv, rcParsing, rcParam, rcInvalid);
                LOGERR(klogErr, rc, "Maximum requested file size is zero");
                break;
            }
        }

        if (self->maxSize > 0 && self->minSize > self->maxSize) {
            rc = RC(rcExe, rcArgv, rcParsing, rcParam, rcInvalid);
            LOGERR(klogErr, rc, "Minimum file size is larger than maximum");
            break;
        }

/* TRANS_OPTION */
        rc = ArgsOptionCount(self->args, TRANS_OPTION, &pcount);
        if (rc != 0) {
            LOGERR(klogErr, rc, "Failure to get '" TRANS_OPTION "' argument");
            break;
        }

        if (pcount > 0) {
            const char *val = NULL;
            rc = ArgsOptionValue(self->args, TRANS_OPTION, 0, &val);
            if (rc != 0) {
                LOGERR(klogErr, rc,
                    "Failure to get '" TRANS_OPTION "' argument value");
                break;
            }
            assert(val);
            if (val[0] == 'a') {
                self->noHttp = true;
            }
            else if (val[0] == 'h') {
                self->noAscp = true;
            }
        }

#ifdef _DEBUGGING
/* TEXTKART_OPTION */
        rc = ArgsOptionCount(self->args, TEXTKART_OPTION, &pcount);
        if (rc != 0) {
            LOGERR(klogErr, rc,
                "Failure to get '" TEXTKART_OPTION "' argument");
            break;
        }

        if (pcount > 0) {
            const char *val = NULL;
            rc = ArgsOptionValue(self->args, TEXTKART_OPTION, 0, &val);
            if (rc != 0) {
                LOGERR(klogErr, rc,
                    "Failure to get '" TEXTKART_OPTION "' argument value");
                break;
            }
            self->textkart = val;
        }
#endif
    } while (false);

    STSMSG(STS_FIN, ("heartbeat = %ld Milliseconds", self->heartbeat));

    return rc;
}

const char UsageDefaultName[] = "prefetch";
rc_t CC UsageSummary(const char *progname) {
    return OUTMSG((
        "Usage:\n"
        "  %s [options] <SRA accession | kart file> [...]\n"
        "  Download SRA or dbGaP files and their dependencies\n"
        "\n"
        "  %s [options] <SRA file> [...]\n"
        "  Check SRA file for missed dependencies "
                                           "and download them\n"
        "\n"
        "  %s --list <kart file> [...]\n"
        "  List the content of a kart file\n"
        , progname, progname, progname));
}

rc_t CC Usage(const Args *args) {
    rc_t rc = 0;
    int i = 0;

    const char *progname = UsageDefaultName;
    const char *fullpath = UsageDefaultName;

    if (args == NULL) {
        rc = RC(rcExe, rcArgv, rcAccessing, rcSelf, rcNull);
    }
    else {
        rc = ArgsProgram(args, &fullpath, &progname);
    }
    if (rc != 0) {
        progname = fullpath = UsageDefaultName;
    }

    UsageSummary(progname);
    OUTMSG(("\n"));

    OUTMSG(("Options:\n"));
    for (i = 0; i < sizeof(Options) / sizeof(Options[0]); i++ ) {
        const char *param = NULL;

        if (Options[i].aliases != NULL) {
            if (strcmp(Options[i].aliases, FAIL_ASCP_ALIAS) == 0) {
                continue; /* debug option */
            }
            if (strcmp(Options[i].aliases, ASCP_ALIAS) == 0) {
                param = "ascp-binary|private-key-file";
            }
            else if (strcmp(Options[i].aliases, FORCE_ALIAS) == 0 ||
                strcmp(Options[i].aliases, HBEAT_ALIAS) == 0 ||
                strcmp(Options[i].aliases, HBEAT_ALIAS) == 0 ||
                strcmp(Options[i].aliases, ORDR_ALIAS) == 0 ||
                strcmp(Options[i].aliases, TRASN_ALIAS) == 0)
            {
                param = "value";
            }
            else if (strcmp(Options[i].aliases, ROWS_ALIAS) == 0) {
                param = "rows";
            }
            else if (strcmp(Options[i].aliases, SIZE_ALIAS) == 0
                  || strcmp(Options[i].aliases, MINSZ_ALIAS) == 0)
            {
                param = "size";
            }
        }
#ifdef _DEBUGGING
        else if (strcmp(Options[i].name, TEXTKART_OPTION) == 0) {
            param = "value";
        }
#endif

        if (Options[i].aliases != NULL &&
            (strcmp(Options[i].aliases, TRASN_ALIAS) == 0 ||
             strcmp(Options[i].aliases, CHECK_ALL_ALIAS) == 0))
        {
            OUTMSG(("\n"));
        }

        HelpOptionLine(Options[i].aliases, Options[i].name,
            param, Options[i].help);
    }

    OUTMSG(("\n"));
    HelpOptionsStandard();
    HelpVersion(fullpath, KAppVersion());

    return rc;
}

ver_t CC KAppVersion(void) { return PREFETCH_VERS; }

/******************************************************************************/

/********** KartTreeNode **********/
static void CC bstKrtWhack(BSTNode *n, void *ignore) {
    KartTreeNode *sn = (KartTreeNode*) n;

    assert(sn);

    ItemRelease(sn->i);

    memset(sn, 0, sizeof *sn);

    free(sn);
}

static int CC bstKrtSort(const BSTNode *item, const BSTNode *n) {
    const KartTreeNode *sn1 = (const KartTreeNode*)item;
    const KartTreeNode *sn2 = (const KartTreeNode*)n;

    assert(sn1 && sn2 && sn1->i && sn2->i);

    return sn1->i->resolved.remoteSz > sn2->i->resolved.remoteSz;
}

static void CC bstKrtDownload(BSTNode *n, void *data) {
    rc_t rc = 0;

    const KartTreeNode *sn = (const KartTreeNode*) n;
    assert(sn && sn->i);

    rc = ItemDownload(sn->i);

    if (rc == 0) {
        rc = ItemPostDownload(sn->i, sn->i->number);
    }
}

/*********** Finalize Main object **********/
static rc_t MainFini(Main *self) {
    rc_t rc = 0;

    assert(self);

    RELEASE(KConfig, self->cfg);
    RELEASE(VResolver, self->resolver);
    RELEASE(VDBManager, self->mgr);
    RELEASE(KDirectory, self->dir);
    RELEASE(KRepositoryMgr, self->repoMgr);
    RELEASE(KNSManager, self->kns);
    RELEASE(VFSManager, self->vfsMgr);
    RELEASE(Args, self->args);

    BSTreeWhack(&self->downloaded, bstWhack, NULL);

    free(self->buffer);

    free((void*)self->ascp);
    free((void*)self->asperaKey);
    free(self->ascpMaxRate);

    memset(self, 0, sizeof *self);

    return rc;
}

/*********** Initialize Main object **********/
static rc_t MainInit(int argc, char *argv[], Main *self) {
    rc_t rc = 0;

    assert(self);
    memset(self, 0, sizeof *self);

    self->heartbeat = 60000;
/*  self->heartbeat = 69; */

    BSTreeInit(&self->downloaded);

    if (rc == 0) {
        rc = MainProcessArgs(self, argc, argv);
    }

    if (rc == 0) {
        self->bsize = 1024 * 1024;
        self->buffer = malloc(self->bsize);
        if (self->buffer == NULL) {
            rc = RC(rcExe, rcData, rcAllocating, rcMemory, rcExhausted);
        }
    }

    if (rc == 0) {
        rc = VFSManagerMake(&self->vfsMgr);
        DISP_RC(rc, "VFSManagerMake");
    }

    if ( rc == 0 )
    {
        rc = VFSManagerGetKNSMgr (self->vfsMgr, & self->kns);
        DISP_RC(rc, "VFSManagerGetKNSMgr");
    }

    if (rc == 0) {
        VResolver *resolver = NULL;
        rc = VFSManagerGetResolver(self->vfsMgr, &resolver);
        DISP_RC(rc, "VFSManagerGetResolver");
        VResolverRemoteEnable(resolver, vrAlwaysEnable);
        RELEASE(VResolver, resolver);
    }

    if (rc == 0) {
        rc = KConfigMake(&self->cfg, NULL);
        DISP_RC(rc, "KConfigMake");
    }

    if (rc == 0) {
        rc = KConfigMakeRepositoryMgrRead(self->cfg, &self->repoMgr);
        DISP_RC(rc, "KConfigMakeRepositoryMgrRead");
    }

    if (rc == 0) {
        rc = VFSManagerMakeResolver(self->vfsMgr, &self->resolver, self->cfg);
        DISP_RC(rc, "VFSManagerMakeResolver");
    }

    if (rc == 0) {
        VResolverRemoteEnable(self->resolver, vrAlwaysEnable);
        VResolverCacheEnable(self->resolver, vrAlwaysEnable);
    }

    if (rc == 0) {
        rc = VDBManagerMakeRead(&self->mgr, NULL);
        DISP_RC(rc, "VDBManagerMakeRead");
    }

    if (rc == 0) {
        rc = KDirectoryNativeDir(&self->dir);
        DISP_RC(rc, "KDirectoryNativeDir");
    }

    if (rc == 0) {
        srand((unsigned)time(NULL));
    }

    return rc;
}

/*********** Process one command line argument **********/
static rc_t MainRun(Main *self, const char *arg, const char *realArg) {
    ERunType type = eRunTypeDownload;
    static bool maxSzPrntd = false;
    rc_t rc = 0;
    Iterator it;
    assert(self && realArg);
    memset(&it, 0, sizeof it);

    if (rc == 0) {
        rc = IteratorInit(&it, arg, self);
    }

    if (self->list_kart_sized) {
        type = eRunTypeList;
    }
    else if (self->order == eOrderSize) {
        if (rc == 0 && it.kart == NULL) {
            type = eRunTypeDownload;
        }
        else {
            type = eRunTypeGetSize;
        }
    }
    else {
        type = eRunTypeDownload;
    }

    if (rc == 0) {
        BSTree trKrt;
        BSTreeInit(&trKrt);
        if (self->list_kart) {
            if (it.kart != NULL) {
                if (self->list_kart_numbered) {
                    rc = KartPrintNumbered(it.kart);
                }
                else {
                    rc = KartPrint(it.kart);
                }
            }
            else {
                PLOGMSG(klogWarn, (klogWarn,
                    "'$(F)' is invalid or not a kart file",
                    "F=%s", realArg));
            }
        }
        else {
            size_t total = 0;
            const char *row = self->rows;
            int64_t n = 1;
            NumIterator nit;
            NumIteratorInit(&nit, row);
            if (type == eRunTypeList) {
                self->maxSize = ~0;
                if (it.kart != NULL) {
                    OUTMSG((
                      "row\tproj-id|item-id|accession|name|item-desc\tsize\n"));
                }
            }
            else {
                if (!maxSzPrntd) {
                    maxSzPrntd = true;
                    if (self->maxSize == 0) {
                        OUTMSG((
                            "Maximum file size download limit is unlimited\n"));
                    }
                    else {
                        OUTMSG(("Maximum file size download limit is %,zuKB\n",
                             self->maxSize / 1024));
                    }
                }
                if (it.kart != NULL) {
                    OUTMSG(("Downloading kart file '%s'\n", realArg));
                    if (type == eRunTypeGetSize) {
                        OUTMSG(("Checking sizes of kart files...\n"));
                    }
                }
                OUTMSG(("\n"));
            }
                
            for (n = 1; ; ++n) {
                rc_t rc2 = 0;
                rc_t rc3 = 0;
                bool done = false;
                Item *item = NULL;
                rc_t rcq = Quitting();
                if (rcq != 0) {
                    if (rc == 0) {
                        rc = rcq;
                    }
                    break;
                }
                rc2 = IteratorNext(&it, &item, &done);
                if (rc2 != 0 || done) {
                    if (rc == 0 && rc2 != 0) {
                        rc = rc2;
                    }
                    break;
                }
                done = ! NumIteratorNext(&nit, n);
                if (done) {
                    break;
                }
                if (!nit.skip) {
                    item->main = self;
                    ResolvedReset(&item->resolved, type);

                    rc3 = ItemProcess(item, (int32_t)n);
                    if (rc3 != 0) {
                        if (rc == 0) {
                            rc = rc3;
                        }
                    }
                    else {
                        if (item->resolved.undersized &&
                            type == eRunTypeGetSize)
                        {
                            STSMSG(STS_TOP,
               ("%d) '%s' (%,zu KB) is smaller than minimum allowed: skipped\n",
                n, item->resolved.name, item->resolved.remoteSz / 1024));
                        }
                        else if (item->resolved.oversized &&
                             type == eRunTypeGetSize)
                        {
                            STSMSG(STS_TOP,
                ("%d) '%s' (%,zu KB) is larger than maximum allowed: skipped\n",
                n, item->resolved.name, item->resolved.remoteSz / 1024));
                        }
                        else {
                            total += item->resolved.remoteSz;

                            if (item != NULL) {
                                if (type == eRunTypeGetSize) {
                                    KartTreeNode *sn = calloc(1, sizeof *sn);
                                    if (sn == NULL) {
                                        return RC(rcExe, rcStorage,
                                           rcAllocating, rcMemory, rcExhausted);
                                    }
                                    if (item->resolved.remoteSz == 0) {
                                        /* remoteSz is unknown:
                     add it to the end of download list preserving kart order */
                                        item->resolved.remoteSz
                                            = (~0ul >> 1) + n + 1;
                                    }
                                    sn->i = item;
                                    item = NULL;
                                    BSTreeInsert(&trKrt, (BSTNode*)sn,
                                        bstKrtSort);
                                }
                            }
                        }
                    }
                }
                else {
                    RELEASE(Item, item);
                }
            }

            if (type == eRunTypeList) {
                if (it.kart != NULL && total > 0) {
                    OUTMSG(("--------------------\ntotal\t%,zuB\n\n", total));
                }
            }
            else if (type == eRunTypeGetSize) {
                OUTMSG(("\nDownloading the files\n\n", realArg));
                BSTreeForEach(&trKrt, false, bstKrtDownload, NULL);
            }
        }
        BSTreeWhack(&trKrt, bstKrtWhack, NULL);
    }
    if (it.isKart) {
        if (self->list_kart) {
            rc_t rc2 = OUTMSG(("\n"));
            if (rc2 != 0 && rc == 0) {
                rc = rc2;
            }
        }
        else if (rc == 0) {
            uint16_t number = 0;
            rc = KartItemsProcessed(it.kart, &number);
            if (rc == 0 && number == 0) {
                PLOGMSG(klogWarn, (klogWarn,
                    "kart file '$(F)' is empty", "F=%s", realArg));
            }
        }
    }
    IteratorFini(&it);

    return rc;
}

/*********** Main **********/
rc_t CC KMain(int argc, char *argv[]) {
    rc_t rc = 0;
    uint32_t pcount = 0;

    Main pars;
    rc = MainInit(argc, argv, &pars);

    if (rc == 0) {
        rc = ArgsParamCount(pars.args, &pcount);
    }
    if (rc == 0 && pcount == 0) {
#ifdef _DEBUGGING
        if (!pars.textkart)
#endif
          rc = UsageSummary(UsageDefaultName);
    }

    if (rc == 0) {
        uint32_t i = ~0;

#ifdef _DEBUGGING
        if (pars.textkart) {
            rc = MainRun(&pars, NULL, pars.textkart);
        }
#endif

        for (i = 0; i < pcount; ++i) {
            const char *obj = NULL;
            rc_t rc2 = ArgsParamValue(pars.args, i, &obj);
            DISP_RC(rc2, "ArgsParamValue");
            if (rc2 == 0) {
                rc2 = MainRun(&pars, obj, obj);
                if (rc2 != 0 && rc == 0) {
                    rc = rc2;
                }
            }
        }

        if (pars.undersized || pars.oversized) {
            OUTMSG(("\n"));
            if (pars.undersized) {
                OUTMSG((
               "Download of some files was skipped because they are too small\n"
                ));
            }
            if (pars.oversized) {
                OUTMSG((
               "Download of some files was skipped because they are too large\n"
                ));
            }
            OUTMSG((
               "You can change size download limit by setting\n"
             "--" MINSZ_OPTION " and --" SIZE_OPTION " command line arguments\n"
            ));
        }
    }

    {
        rc_t rc2 = MainFini(&pars);
        if (rc2 != 0 && rc == 0) {
            rc = rc2;
        }
    }

    return rc;
}
