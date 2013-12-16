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

#include "prefetch.vers.h"

#include <kapp/main.h> /* KAppVersion */

#include <kfg/config.h> /* KConfig */
#include <kfg/kart.h> /* Kart */
#include <kfg/repository.h> /* KRepositoryMgr */

#include <vdb/database.h> /* VDatabase */
#include <vdb/dependencies.h> /* VDBDependencies */
#include <vdb/manager.h> /* VDBManager */

#include <vfs/manager.h> /* VFSManager */
#include <vfs/path.h> /* VPath */
#include <vfs/resolver.h> /* VResolver */

#include <kns/ascp.h> /* ascp_locate */
#include <kns/curl-file.h> /* KCurlFileMake */

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

typedef struct {
    const VPath *path;
    const String *str;
} VPathStr;
static rc_t VPathStrFini(VPathStr *self) {
    rc_t rc = 0;

    assert(self);

    VPathRelease(self->path);

    free((void*)self->str);

    memset(self, 0, sizeof *self);

    return rc;
}

typedef struct {
    /* "plain" command line argument */
    const char *desc;

    const KartItem *item;

#ifdef _DEBUGGING
    const char *textkart;
#endif
} Item;
static rc_t ItemInit(Item *self, const char *obj) {
    assert(self);
    self->desc = obj;
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

static rc_t V_ResolverRemote ( const VResolver * self,
    VRemoteProtocols protocols, struct VPath const * accession,
    struct VPath const ** remote, struct KFile const ** opt_file_rtn,
    struct VPath const ** cache)
{
    return VResolverQuery(self, protocols, accession, NULL, remote, cache);
}

static rc_t _VResolverRemote(VResolver *self, VRemoteProtocols protocols,
    const char *name, const VPath *vaccession,
    const VPath **vremote, const String **remote,
    const String **cache)
{
    rc_t rc = 0;
    const VPath *vcache = NULL;
    assert(vaccession && vremote && !*vremote);
    rc = V_ResolverRemote(self, protocols, vaccession, vremote, NULL, &vcache);
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
            StringInit(&local_str, path, len, len);
            rc = StringCopy(remote, &local_str);
            DISP_RC2(rc, "StringCopy(VResolverRemote)", name);
        }
    }
    else if (NotFoundByResolver(rc)) {
        PLOGERR(klogErr, (klogErr, rc, "$(acc) cannot be found.",
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
                "for $(acc). Try to cd out of protected repository.",
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

static rc_t V_ResolverLocal ( const VResolver * self,
    struct VPath const * accession, struct VPath const ** path )
{
    return VResolverQuery(self, eProtocolHttp, accession, path, NULL, NULL);
}

static rc_t _VResolverResolve(VResolver *self, VRemoteProtocols protocols,
    const Item *item, VPathStr *local, const String **remote,
    const String **cache, const KFile **file, VPath **vaccession,
    uint64_t *project, const KRepositoryMgr *repoMgr, const KConfig *cfg,
    const char *name)
{
    VResolver *resolver = NULL;
    const char *accession = NULL;
    rc_t rc = 0;
    rc_t rc2 = 0;
    VFSManager *vfs = NULL;

    const VPath **vlocal = NULL;
    const VPath *vremote = NULL;

    if (item == NULL) {
        return RC(rcExe, rcResolver, rcResolving, rcParam, rcNull);
    }
    if (self == NULL) {
        return RC(rcExe, rcResolver, rcResolving, rcSelf, rcNull);
    }
    if (local != NULL) {
        memset(local, 0, sizeof *local);
        vlocal = &local->path;
    }

    assert(vaccession && !*vaccession && project && !*project);

    if (rc == 0) {
        rc = VFSManagerMake(&vfs);
    }

    accession = item->desc;
    if (accession == NULL) {
        uint64_t oid = 0;
        rc = KartItemProjIdNumber(item->item, project);
        if (rc != 0) {
            DISP_RC(rc, "KartItemProjIdNumber");
            return rc;
        }
        rc = KartItemItemIdNumber(item->item, &oid);
        if (rc != 0) {
            DISP_RC(rc, "KartItemItemIdNumber");
            return rc;
        }
        {
            const KRepository *p_protected = NULL;
            rc = KRepositoryMgrGetProtectedRepository(repoMgr, 
                *project, &p_protected);
            if (rc == 0) {
                rc = KRepositoryMakeResolver(p_protected, &resolver, cfg);
                if (rc == 0) {
                    self = resolver;
                }
            }
            RELEASE(KRepository, p_protected);
        }
        rc = VFSManagerMakeOidPath(vfs, vaccession, oid);
    }
    else {
        rc = VFSManagerMakePath(vfs, vaccession, accession);
        DISP_RC2(rc, "VFSManagerMakePath", accession);
    }

    RELEASE(VFSManager, vfs);

    if (rc == 0 && local != NULL) {
        rc = V_ResolverLocal(self, *vaccession, vlocal);
        if (rc == 0) {
            rc = VPathMakeString(*vlocal, &local->str);
            DISP_RC2(rc, "VPathMakeString(VResolverLocal)", name);
        }
        else if (NotFoundByResolver(rc)) {
            rc = 0;
        }
        else {
            DISP_RC2(rc, "VResolverLocal", name);
        }
    }

    if (remote != NULL) {
        rc2 = _VResolverRemote(self, protocols, name,
            *vaccession, &vremote, remote, cache);
        if (rc2 != 0 && rc == 0) {
            rc = rc2;
        }
    }

    RELEASE(VPath, vremote);
    RELEASE(VResolver, resolver);
    return rc;
}

typedef enum {
    eForceNo, /* do not download found and complete objects */
    eForceYes,/* force download of found and complete objects */
    eForceYES /* force download; ignore lockes */
} Force;
typedef struct {
    BSTNode n;
    char *path;
} TreeNode;
typedef struct {
    Args *args;
    bool check_all;
    bool list_kart;
    Force force;
    KConfig *cfg;
    const VDBManager *mgr;
    KDirectory *dir;
    const KRepositoryMgr *repoMgr;
    VResolver *resolver;

    void *buffer;
    size_t bsize;

    BSTree downloaded;

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
static rc_t StringRelease(const String *self) {
    free((String*)self);
    return 0;
}

typedef struct {
    VPathStr local;
    const String *cache;
    const String *remote;

    const KFile *file;
    uint64_t remoteSz;

    bool existing;

    /* path to the resolved object : either local or cache:
    should not be released */
    const char *path;

    char *name;
    VPath *accession;
    uint64_t project;
} Resolved;
static rc_t ResolvedFini(Resolved *self) {
    rc_t rc = 0;

    assert(self);

    rc = VPathStrFini(&self->local);

    RELEASE(KFile, self->file);
    RELEASE(VPath, self->accession);

    RELEASE(String, self->remote);
    RELEASE(String, self->cache);

    free(self->name);

    memset(self, 0, sizeof *self);

    return rc;
}

static rc_t ResolvedInit(Resolved *self,
    const Item *item, VResolver *resolver, KDirectory *dir, bool ascp,
    const KRepositoryMgr *repoMgr, const KConfig *cfg)
{
    rc_t rc = 0;
    KPathType type = kptNotFound;
    VRemoteProtocols protocols = ascp ? eProtocolFaspHttp : eProtocolHttp;

    assert(self && item);

    memset(self, 0, sizeof *self);

    self->name = ItemName(item);

    type = KDirectoryPathType(dir, item->desc) & ~kptAlias;
    if (type == kptFile || type == kptDir) {
        self->path = item->desc;
        self->existing = true;
        return 0;
    }

    rc = _VResolverResolve(resolver, protocols, item, &self->local,
        &self->remote, &self->cache, &self->file, &self->accession,
        &self->project, repoMgr, cfg, self->name);

    if (rc == 0 && self->file != NULL) {
        rc_t rc2 = KFileSize(self->file, &self->remoteSz);
        DISP_RC2(rc2, "KFileSize(remote)", self->name);
    }


    STSMSG(STS_DBG, ("Resolve(%s) = %R:", self->name, rc));
    STSMSG(STS_DBG, ("local(%s)",
        self->local.str ? self->local.str->addr : "NULL"));
    STSMSG(STS_DBG, ("cache(%s)", self->cache ? self->cache->addr : "NULL"));
    STSMSG(STS_DBG, ("remote(%s:%,ld)",
        self->remote ? self->remote->addr : "NULL", self->remoteSz));

    if (rc == 0) {
        if (self->local.str == NULL
            && (self->cache == NULL || self->remote == NULL))
        {
            rc = RC(rcExe, rcPath, rcValidating, rcParam, rcNull);
            PLOGERR(klogInt, (klogInt, rc,
                "bad VResolverResolve($(acc)) result", "acc=%s", self->name));
        }
    }

    return rc;
}

static bool _StringIsFasp(const String *self, const char **withoutScheme) {
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

/** isLocal is set to true when the object is found locally.
    i.e. does not need need not be [re]downloaded */
static rc_t ResolvedLocal(const Resolved *self,
    const KDirectory *dir, bool *isLocal, Force force)
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

    if (rc == 0 && KDirectoryPathType(dir, path) != kptFile) {
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

    if (rc == 0 && ! _StringIsFasp(self->remote, NULL) && self->file != NULL) {
        rc = KFileSize(self->file, &sRemote);
        DISP_RC2(rc, "KFileSize(remote)", self->name);
    }

    if (rc == 0) {
        rc = KDirectoryOpenFileRead(dir, &local, path);
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

static rc_t _KDirectoryMkTmpPrefix(const KDirectory *self,
    const String *prefix, char *out, size_t sz)
{
    size_t num_writ = 0;
    rc_t rc = string_printf(out, sz, &num_writ,
        "%.*s.tmp", prefix->size, prefix->addr);
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
            "%.*s.tmp.%d.tmp", prefix->size, prefix->addr, rand() % 100000);
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
        if (KDirectoryPathType(self, out) == kptNotFound) {
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

    rc = string_printf(out, sz, &num_writ,
        "%.*s.lock", prefix->size, prefix->addr);
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

    if (rc == 0 && KDirectoryPathType(self, cache) != kptNotFound) {
        STSMSG(STS_DBG, ("removing %s", cache));
        rc = KDirectoryRemove(self, false, cache);
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

    if (tmp != NULL && KDirectoryPathType(self, tmp) != kptNotFound) {
        rc_t rc3 = 0;
        STSMSG(STS_DBG, ("removing %s", tmp));
        rc3 = KDirectoryRemove(self, false, tmp);
        if (rc2 == 0 && rc3 != 0) {
            rc2 = rc3;
        }
    }

    if (rmSelf && KDirectoryPathType(self, cache->addr) != kptNotFound) {
        rc_t rc3 = 0;
        STSMSG(STS_DBG, ("removing %s", cache->addr));
        rc3 = KDirectoryRemove(self, false, cache->addr);
        if (rc2 == 0 && rc3 != 0) {
            rc2 = rc3;
        }
    }

    if (rc == 0) {
        uint32_t count = 0;
        uint32_t i = 0;
        KNamelist *list = NULL;
        STSMSG(STS_DBG, ("listing %s for old temporary files", dir));
        rc = KDirectoryList(self, &list, NULL, NULL, dir);
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

    if (lock != NULL && KDirectoryPathType(self, lock) != kptNotFound) {
        rc_t rc3 = 0;
        STSMSG(STS_DBG, ("removing %s", lock));
        rc3 = KDirectoryRemove(self, false, lock);
        if (rc2 == 0 && rc3 != 0) {
            rc2 = rc3;
        }
    }

    if (rc == 0 && rc2 != 0) {
        rc = rc2;
    }

    return rc;
}

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

/* #if ! WINDOWS #endif */
    free(sn->path);

    memset(sn, 0, sizeof *sn);

    free(sn);
}
/*
static bool _KConfigAscpDisabled(const KConfig *self) {
    bool disabled = false;
    const char path[] = "tools/ascp/disabled";
    rc_t rc = KConfigReadBool(self, path, &disabled);
    if (rc != 0) {
        if (rc != SILENT_RC(rcKFG, rcNode, rcOpening, rcPath, rcNotFound)) {
            DISP_RC(rc, path);
        }
        else {
            STSMSG(STS_DBG, ("'%s': not found in configuration", path));
        }
        disabled = false;
    }
    else {
        STSMSG(STS_DBG, ("'%s' = '%s'", path, disabled ? "true" : "false"));
    }
    return disabled;
}

static String* _KConfigAscpString(const KConfig *self,
    const char *path, const char *name)
{
    String *ascp = NULL;
    rc_t rc = KConfigReadString(self, path, &ascp);
    if (rc == 0) {
        assert(ascp);
        STSMSG(STS_INFO, ("Using %s from configuration: '%s'",
            name, ascp->addr));
        return ascp;
    }
    else {
        if (rc != SILENT_RC(rcKFG, rcNode, rcOpening, rcPath, rcNotFound)) {
            DISP_RC(rc, path);
        }
        else {
            STSMSG(STS_DBG, ("'%s': not found in configuration", path));
        }
        free(ascp);
        return NULL;
    }
}

static String* _StringMake(const char *s) {
    size_t len = string_size(s);
    size_t tmp = 0;
    String *str = malloc(sizeof *str + len + 1);
    if (str == NULL) {
        return NULL;
    }
    StringInit(str,
        (char*)str + sizeof *str, len, len + 1);
    tmp = string_copy((char*)str + sizeof *str, len + 1, s, len + 1);
    assert(tmp == len + 1);
    return str;
}

static String* _SystemHelp(String* prev, const char *command) {
    int value = 0;
    if (prev != NULL) {
        return prev;
    }
    STSMSG(STS_DBG, ("Checking '%s'", command));
    value = _SilentSystem("%s -h", command);
    if (value == 0) {
        STSMSG(STS_INFO, ("Using '%s'", command));
        return _StringMake(command);
    }
    else {
        STSMSG(STS_DBG, ("'%s': not found", command));
        return NULL;
    }
}

static String* _KDirectoryFileFound(const KDirectory *self,
    String* prev, const char *path)
{
    KPathType type = kptNotFound;
    if (prev != NULL) {
        return prev;
    }
    STSMSG(STS_DBG, ("Checking '%s'", path));
    type = KDirectoryPathType(self, path);
    if ((type & ~kptAlias) == kptFile) {
        STSMSG(STS_DBG, ("'%s': found", path));
        return _StringMake(path);
    }
    else {
        STSMSG(STS_DBG, ("'%s': not found", path));
        return NULL;
    }
}
*/
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

typedef struct {
    const char *obj;
    bool done;
    Kart *kart;
    bool isKart;
} Iterator;
static
rc_t IteratorInit(Iterator *self, const char *obj, const Main *main)
{
    rc_t rc = 0;

    KPathType type = kptNotFound;

    assert(self && main);
    memset(self, 0, sizeof *self);

#ifdef _DEBUGGING
    if (obj == NULL && main->textkart) {
        type = KDirectoryPathType(main->dir, main->textkart);
        if ((type & ~kptAlias) != kptFile) {
            rc = RC(rcExe, rcFile, rcOpening, rcFile, rcNotFound);
            DISP_RC(rc, main->textkart);
            return rc;
        }
        rc = KartMakeText(main->dir, main->textkart, &self->kart,
            &self->isKart);
        if (!self->isKart) {
            rc = 0;
        }
        return rc;
    }
#endif

    assert(obj);
    type = KDirectoryPathType(main->dir, obj);
    if ((type & ~kptAlias) == kptFile) {
        type = VDBManagerPathType(main->mgr, obj);
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

static rc_t IteratorNext(Iterator *self, Item *next, bool *done) {
    rc_t rc = 0;

    assert(self && next && done);

    memset(next, 0, sizeof *next);

    if (self->done) {
        *done = true;
        return 0;
    }

    *done = false;

    if (self->isKart) {
        RELEASE(KartItem, next->item);
        rc = KartMakeNextItem(self->kart, &next->item);
        if (rc != 0) {
            LOGERR(klogErr, rc, "Invalid kart file row");
        }
        else if (next->item == NULL) {
            *done = true;
        }

        if (rc == 0 && *done) {
            self->done = true;
        }
    }
    else {
        rc = ItemInit(next, self->obj);

        self->done = true;
    }

    return rc;
}

static void IteratorFini(Iterator *self) {
    rc_t rc = 0;

    assert(self);

    RELEASE(Kart, self->kart);
}

static rc_t ResolvedDownloadFile(Resolved *self,
    Main *main, const char *to)
{
    rc_t rc = 0;
    KFile *out = NULL;
    uint64_t pos = 0;
    size_t num_read = 0;
    uint64_t opos = 0;
    size_t num_writ = 0;

    assert(self && main);

    if (rc == 0) {
        STSMSG(STS_DBG, ("creating %s", to));
        rc = KDirectoryCreateFile(main->dir, &out,
            false, 0664, kcmInit | kcmParents, to);
        DISP_RC2(rc, "Cannot OpenFileWrite", to);
    }

    assert(self->remote);

    if (self->file == NULL) {
#define USE_CURL 1
#if USE_CURL
        rc = KCurlFileMake(&self->file, self->remote->addr, false);
#else
        rc =
            KNSManagerMakeHttpFile(kns, &self->file, NULL, 0x01010000, "%S", s);
#endif
        if (rc != 0) {
            PLOGERR(klogInt, (klogInt, rc, "failed to open file for $(path)",
                "path=%s", self->remote->addr));
        }
    }

    STSMSG(STS_INFO, ("%s -> %s", self->remote->addr, to));
    do {
        rc = Quitting();

        if (rc == 0) {
            STSMSG(STS_FIN,
                ("> Reading %lu bytes from pos. %lu", main->bsize, pos));
            rc = KFileRead(self->file,
                pos, main->buffer, main->bsize, &num_read);
            if (rc != 0) {
                DISP_RC2(rc, "Cannot KFileRead", self->remote->addr);
            }
            else {
                STSMSG(STS_FIN,
                    ("< Read %lu bytes from pos. %lu", num_read, pos));
                pos += num_read;
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
static rc_t ResolvedDownloadAscp(const Resolved *self, Main *main,
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

static rc_t ResolvedDownload(Resolved *self, Main *main) {
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

    if (KDirectoryPathType(main->dir, lock) != kptNotFound) {
        if (main->force != eForceYES) {
            KTime_t date = 0;
            rc = KDirectoryDate(main->dir, &date, lock);
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
            false, 0664, kcmInit | kcmParents, lock);
        DISP_RC2(rc, "Cannot OpenFileWrite", lock);
    }

    assert(!main->noAscp || !main->noHttp);

    if (rc == 0) {
        bool ascp = _StringIsFasp(self->remote, NULL);
        if (ascp) {
            STSMSG(STS_TOP, (" Downloading via fasp..."));
            rc = ResolvedDownloadAscp(self, main, tmp);
            if (rc == 0 && main->forceAscpFail) {
                rc = 1;
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
            RELEASE(String, self->remote);
            RELEASE(KFile, self->file);
            rc = _VResolverRemote(main->resolver, eProtocolHttp, self->name,
                self->accession, &vremote, &self->remote, &self->cache);
            if (rc == 0) {
                rc = ResolvedDownloadFile(self, main, tmp);
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

static
rc_t ResolvedResolve(Resolved *self, Main *main, const Item *item)
{
    static int n = 0;
    rc_t rc = 0;
    bool isLocal = false;
    bool ascp = false;

    assert(self && main);

    ++n;

    ascp = MainUseAscp(main);

    rc = ResolvedInit(self, item,
        main->resolver, main->dir, ascp, main->repoMgr, main->cfg);

    if (rc == 0) {
        if (self->existing) {
            self->path = item->desc;
            return rc;
        }

        rc = ResolvedLocal(self, main->dir, &isLocal, main->force);
    }

    if (rc == 0) {
        if (isLocal) {
            STSMSG(STS_TOP, ("%d) %s is found locally", n, self->name));
            if (self->local.str != NULL) {
                self->path = self->local.str->addr;
            }
        }
        else if (!_StringIsFasp(self->remote, NULL) && main->noHttp) {
            rc = RC(rcExe, rcFile, rcCopying, rcFile, rcNotFound);
            PLOGERR(klogErr, (klogErr, rc,
                "cannot download '$(name)' using requested transport",
                "name=%s", self->name));
        }
        else {
            STSMSG(STS_TOP, ("%d) Downloading %s...", n, self->name));
            rc = ResolvedDownload(self, main);
            if (rc == 0) {
                STSMSG(STS_TOP,
                    ("%d) %s was downloaded successfully", n, self->name));
                if (self->cache != NULL) {
                    self->path = self->cache->addr;
                }
            }
            else if (rc !=
                SILENT_RC(rcExe, rcProcess, rcExecuting, rcProcess, rcCanceled))
            {
                STSMSG(STS_TOP, ("%d) failed to download %s", n, self->name));
            }
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

    rc = VDBManagerOpenDBRead(self->mgr, &db, NULL, path);
    if (rc != 0) {
        if (rc == SILENT_RC(rcDB, rcMgr, rcOpening, rcDatabase, rcIncorrect)) {
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

static rc_t MainExecute(Main *self, const Item *item) {
    rc_t rc = 0;
    rc_t rc2 = 0;
    const VDBDependencies *deps = NULL;
    uint32_t count = 0;
    uint32_t i = 0;

    Resolved resolved;

    assert(self && item);

    if (item->desc == NULL) {
        OUTMSG(("Requested kart row download\n"));
    }

    /* resolve: locate; download if not found */
    rc = ResolvedResolve(&resolved, self, item);

    if (rc == 0) {
        rc = MainDependenciesList(self, resolved.path, &deps);
    }

    /* resolve dependencies (refseqs) */
    if (rc == 0 && deps != NULL) {
        rc = VDBDependenciesCount(deps, &count);
        if (rc == 0) {
            STSMSG(STS_TOP, ("'%s' has %d%s dependenc%s",
                item->desc, count, self->check_all ? "" : " unresolved",
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
                Resolved resolved;
                Item item;
                item.desc = ncbiAcc;
                rc = ResolvedResolve(&resolved, self, &item);
                rc2 = ResolvedFini(&resolved);
            }
        }
    }

    if (rc == 0 && rc2 != 0) {
        rc = rc2;
    }

    RELEASE(VDBDependencies, deps);

    {
        rc_t rc2 = ResolvedFini(&resolved);
        if (rc == 0 && rc2 != 0) {
            rc = rc2;
        }
    }

    return rc;
}

static rc_t MainFini(Main *self) {
    rc_t rc = 0;

    assert(self);

    RELEASE(KConfig, self->cfg);
    RELEASE(VResolver, self->resolver);
    RELEASE(VDBManager, self->mgr);
    RELEASE(KDirectory, self->dir);
    RELEASE(KRepositoryMgr, self->repoMgr);
    RELEASE(Args, self->args);

    BSTreeWhack(&self->downloaded, bstWhack, NULL);

    free(self->buffer);

    free((void*)self->ascp);
    free((void*)self->asperaKey);
    free(self->ascpMaxRate);

    memset(self, 0, sizeof *self);

    return rc;
}

#define CHECK_ALL_OPTION "check-all"
#define CHECK_ALL_ALIAS  "c"
static const char* CHECK_ALL_USAGE[] = { "double-check all refseqs", NULL };

#define FORCE_OPTION "force"
#define FORCE_ALIAS  "f"
static const char* FORCE_USAGE[] = {
    "force object download - one of: no, yes, all.",
    "no [default]: skip download if the object if found and complete;",
    "yes: download it even if it is found and is complete;", "all: ignore lock "
    "files (stale locks or it is beeing downloaded by another process)", NULL };

#define FAIL_ASCP_OPTION "FAIL-ASCP"
#define FAIL_ASCP_ALIAS  "F"
static const char* FAIL_ASCP_USAGE[] = {
    "force ascp download fail to test ascp->http download combination" };


#define LIST_OPTION "list"
#define LIST_ALIAS  "l"
static const char* LIST_USAGE[] = { "list the content of a kart file", NULL };

#define HBEAT_ALIAS  "p"
#define HBEAT_OPTION "progress"
static const char* HBEAT_USAGE[] = {
    "time period in minutes to display download progress",
    "(0: no progress), default: 1", NULL };

#define TRANS_OPTION "transport"
#define TRASN_ALIAS  "t"
static const char* TRANS_USAGE[] = { "transport: one of: ascp; http; both.",
    "(ascp only; http only; first try ascp, "
    "use http if cannot download by ascp).",
    "Default: both", NULL };

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
   ,{ HBEAT_OPTION    , HBEAT_ALIAS    , NULL, HBEAT_USAGE , 1, true, false }
   ,{ CHECK_ALL_OPTION, CHECK_ALL_ALIAS, NULL, CHECK_ALL_USAGE, 1, false, false}
   ,{ FAIL_ASCP_OPTION, FAIL_ASCP_ALIAS, NULL, FAIL_ASCP_USAGE, 1, false, false}
#ifdef _DEBUGGING
   ,{ TEXTKART_OPTION , NULL           , NULL, TEXTKART_USAGE , 1, true , false}
#endif
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

        rc = ArgsOptionCount(self->args, CHECK_ALL_OPTION, &pcount);
        if (rc != 0) {
            LOGERR(klogErr,
                rc, "Failure to get '" CHECK_ALL_OPTION "' argument");
            break;
        }
        if (pcount > 0 || self->force != eForceNo) {
            self->check_all = true;
        }

        rc = ArgsOptionCount(self->args, FAIL_ASCP_OPTION, &pcount);
        if (rc != 0) {
            LOGERR(klogErr,
                rc, "Failure to get '" FAIL_ASCP_OPTION "' argument");
            break;
        }
        if (pcount > 0) {
            self->forceAscpFail = true;
        }

        rc = ArgsOptionCount(self->args, LIST_OPTION, &pcount);
        if (rc != 0) {
            LOGERR(klogErr,
                rc, "Failure to get '" LIST_OPTION "' argument");
            break;
        }
        if (pcount > 0) {
            self->list_kart = true;
        }

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
            self->heartbeat = f;
        }

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
            if (strcmp(Options[i].aliases, FORCE_ALIAS) == 0) {
                param = "value";
            }
        }
        HelpOptionLine(Options[i].aliases, Options[i].name,
            param, Options[i].help);
    }

    OUTMSG(("\n"));
    HelpOptionsStandard();
    HelpVersion(fullpath, KAppVersion());

    return rc;
}

static rc_t MainInit(int argc, char *argv[], Main *self) {
    rc_t rc = 0;
    VFSManager *mgr = NULL;

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
        rc = VFSManagerMake(&mgr);
        DISP_RC(rc, "VFSManagerMake");
    }

    if (rc == 0) {
        VResolver *resolver = NULL;
        rc = VFSManagerGetResolver(mgr, &resolver);
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
        rc = VFSManagerMakeResolver(mgr, &self->resolver, self->cfg);
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
        srand(time(NULL));
    }

    RELEASE(VFSManager, mgr);

    return rc;
}

ver_t CC KAppVersion(void) { return PREFETCH_VERS; }

/******************************************************************************/

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
            rc_t rc2 = 0;
            Iterator it;
            memset(&it, 0, sizeof it);
            if (rc2 == 0) {
                rc2 = IteratorInit(&it, NULL, &pars);
            }
            if (rc2 == 0) {
                if (pars.list_kart) {
                    rc2 = KartPrint(it.kart);
                }
                else {
                    while (true) {
                        rc_t rc3 = 0;
                        bool done = false;
                        Item item;
                        rc2 = IteratorNext(&it, &item, &done);
                        if (rc2 != 0 || done == true) {
                            break;
                        }
                        rc3 = MainExecute(&pars, &item);
                        if (rc == 0) {
                            if (rc2 != 0) {
                                rc = rc2;
                            }
                            else if (rc3 != 0) {
                                rc = rc3;
                            }
                        }
                    }
                }
            }
        }
#endif

        for (i = 0; i < pcount; ++i) {
            rc_t rc2 = 0;
            Iterator it;
            const char *obj = NULL;

            memset(&it, 0, sizeof it);

            rc2 = ArgsParamValue(pars.args, i, &obj);
            DISP_RC(rc2, "ArgsParamValue");

            if (rc2 == 0) {
                rc2 = IteratorInit(&it, obj, &pars);
            }

            if (rc2 == 0) {
                if (pars.list_kart) {
                    rc2 = KartPrint(it.kart);
                }
                else {
                    while (true) {
                        rc_t rc3 = 0;
                        bool done = false;
                        Item item;
                        rc2 = IteratorNext(&it, &item, &done);
                        if (rc2 != 0 || done == true) {
                            break;
                        }
                        rc3 = MainExecute(&pars, &item);
                        if (rc == 0) {
                            if (rc2 != 0) {
                                rc = rc2;
                            }
                            else if (rc3 != 0) {
                                rc = rc3;
                            }
                        }
                    }
                }
            }

            if (pars.list_kart && it.isKart) {
                rc_t rc2 = OUTMSG(("\n"));
                if (rc2 != 0 && rc == 0) {
                    rc = rc2;
                }
            }

            IteratorFini(&it);
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
