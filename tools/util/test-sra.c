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

#include <kapp/main.h> /* KMain */

#include <kfg/config.h> /* KConfig */
#include <kfg/kart.h> /* Kart */
#include <kfg/repository.h> /* KRepositoryMgr */

#include <vfs/manager.h> /* VFSManager */
#include <vfs/path.h> /* VPath */
#include <vfs/resolver.h> /* VResolver */

#include <vdb/database.h> /* VDatabase */
#include <vdb/dependencies.h> /* VDBDependencies */
#include <vdb/manager.h> /* VDBManager */

#include <kns/curl-file.h> /* KCurlFileMake */
#include <kns/manager.h> /* KNSManager */

#include <kdb/manager.h> /* kptDatabase */

#include <kfs/directory.h> /* KDirectory */
#include <kfs/file.h> /* KFile */

#include <klib/log.h> /* KLogHandlerSet */
#include <klib/out.h> /* KOutMsg */
#include <klib/printf.h> /* string_vprintf */
#include <klib/rc.h>
#include <klib/text.h> /* String */

#include <sysalloc.h>

#include <assert.h>
#include <ctype.h> /* isprint */
#include <stdlib.h> /* calloc */
#include <string.h> /* memset */

VFS_EXTERN rc_t CC VResolverProtocols ( VResolver * self,
    VRemoteProtocols protocols );

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

#define RELEASE(type, obj) do { rc_t rc2 = type##Release(obj); \
    if (rc2 != 0 && rc == 0) { rc = rc2; } obj = NULL; } while (false)

typedef enum {
    eCfg = 1,
/*  eType = 2, */
    eResolve = 2,
    eDependMissing = 4,
    eDependAll = 8,
    eCurl = 16,
    eAll = 32
} Type;
typedef struct {
    KConfig *cfg;
    KDirectory *dir;
    const VDBManager *mgr;
    const KRepositoryMgr *repoMgr;
    VResolver *resolver;
    uint8_t tests;
    bool recursive;
    bool noVDBManagerPathType;
    bool noRfs;

    bool allowCaching;
    VResolverEnableState cacheState;
} Main;

uint32_t CC KAppVersion(void) { return 0; }

const char UsageDefaultName[] = "test-sra";

rc_t CC UsageSummary(const char *prog_name) {
    return KOutMsg(
        "Usage:\n"
        "    %s [+crdDa] [-crdDa] [-R] [-N] [-C] [options] name [ name... ]\n",
        prog_name);
}

rc_t CC Usage(const Args *args) {
    rc_t rc2 = 0;

    const char *progname, *fullpath;
    rc_t rc = ArgsProgram(args, &fullpath, &progname);
    if (rc != 0) {
        progname = fullpath = UsageDefaultName;
    }

    rc2 = UsageSummary(progname);
    if (rc == 0 && rc2 != 0) {
        rc = rc2;
    }

/*      "  t - test object's VDB type\n" */
    rc2 = KOutMsg("\n"
        "Test [SRA] object, resolve it, print dependencies, configuration\n\n"
        "[+tests] - add tests\n"
        "[-tests] - remove tests\n\n"
        "Tests:\n"
        "  c - print configuration\n"
        "  k - print curl info\n"
        "  r - call VResolver\n"
        "  d - call ListDependencies(missing)\n"
        "  D - call ListDependencies(all)\n"
        "  a - all tests\n\n"
        "If no tests were specified then all tests will be run\n\n"
        "-R - check object type recursively\n"
        "-N - do not call VDBManagerPathType\n"
        "-C - do not disable caching (default: from configuration)\n\n"
        "More options:\n");
    if (rc == 0 && rc2 != 0) {
        rc = rc2;
    }

    HelpOptionsStandard();

    return rc;
}

static bool testArg(const char *arg, uint8_t *testOn, uint8_t *testOff) {
    int j = 1;
    uint8_t *res = NULL;

/*  const char tests[] = "ctrdDa"; */
    const char tests[] = "crdDka";

    assert(arg && testOn && testOff);
    if (arg[0] != '+' && arg[0] != '-') {
        return false;
    }

    if (arg[0] == '-' &&
        arg[1] != '\0' && strchr(tests, arg[1]) == NULL)
    {
        return false;
    }

    res = arg[0] == '-' ? testOff : testOn;

    for (j = 1; arg[j] != '\0'; ++j) {
        char *c = strchr(tests, arg[j]);
        if (c != NULL) {
            int offset = c - tests;
            *res |= 1 << offset;
        }
    }

    return true;
}

static uint8_t Turn(uint8_t in, uint8_t tests, bool on) {
    uint8_t c = 1;
    for (c = 1; c < eAll; c <<= 1) {
        if (tests & c) {
            if (on) {
                in |= c;
            }
            else {
                in &= ~c;
            }
        }
    }
    return in;
}

static uint8_t processTests(uint8_t testsOn, uint8_t testsOff) {
    uint8_t tests = 0;

    bool allOn = false;
    bool allOff = false;

    if (testsOn & eAll && testsOff & eAll) {
        testsOn &= ~eAll;
        testsOff &= ~eAll;
    }
    else if (testsOn & eAll) {
        allOn = true;
    }
    else if (testsOff & eAll) {
        allOff = true;
    }

    if (allOn) {
        tests = ~0;
        tests = Turn(tests, testsOff, false);
    }
    else if (allOff) {
        tests = Turn(tests, testsOn, true);
    }
    else if (testsOn != 0 || testsOff != 0) {
        tests = Turn(tests, testsOff, false);
        tests = Turn(tests, testsOn, true);
    }
    else {
        tests = ~0;
    }

    return tests;
} 

static bool MainHasTest(const Main *self, Type type) {
    assert(self);
    return self->tests & type;
}

static void MainPrint(const Main *self) {
    return;

    assert(self);

    if (MainHasTest(self, eCfg)) {
        KOutMsg("eCfg\n");
    }

/*  if (MainHasTest(self, eType)) {
        KOutMsg("eType\n");
    }*/

    if (MainHasTest(self, eResolve)) {
        KOutMsg("eResolve\n");
    }

    if (MainHasTest(self, eDependMissing)) {
        KOutMsg("eDependMissing\n");
    }

    if (MainHasTest(self, eDependAll)) {
        KOutMsg("eDependAll\n");
    }
}

static rc_t MainInitObjects(Main *self) {
    rc_t rc = 0;

    VFSManager* mgr = NULL;
    VResolver *resolver = NULL;

    if (rc == 0) {
        rc = KDirectoryNativeDir(&self->dir);
    }

    if (rc == 0) {
        rc = VDBManagerMakeRead(&self->mgr, NULL);
    }

    if (rc == 0) {
        rc = KConfigMake(&self->cfg, NULL);
    }

    if (rc == 0) {
        rc = KConfigMakeRepositoryMgrRead(self->cfg, &self->repoMgr);
    }

    if (rc == 0) {
        rc = VFSManagerMake(&mgr);
    }

    if (rc == 0) {
        rc = VFSManagerGetResolver(mgr, &resolver);
    }

    if (rc == 0) {
        if (!self->allowCaching) {
            self->cacheState = VResolverCacheEnable(resolver, vrAlwaysDisable);
        }
    }

    if (rc == 0) {
        rc = VFSManagerMakeResolver(mgr, &self->resolver, self->cfg);
    }


    RELEASE(VResolver, resolver);
    RELEASE(VFSManager, mgr);

    return rc;
}

static
void _MainInit(Main *self, int argc, char *argv[], int *argi, char **argv2)
{
    int i = 0;

    uint8_t testsOn = 0;
    uint8_t testsOff = 0;

    assert(self && argv && argi && argv2);

    *argi = 0;
    argv2[(*argi)++] = argv[0];

    for (i = 1; i < argc; ++i) {
        if (!testArg(argv[i], &testsOn, &testsOff)) {
            argv2[(*argi)++] = argv[i];
        }
    }

    self->tests = processTests(testsOn, testsOff);

    MainPrint(self);
}

static char** MainInit(Main *self, int argc, char *argv[], int *argi) {
    char **argv2 = calloc(argc, sizeof *argv2);

    assert(self);

    memset(self, 0, sizeof *self);

    if (argv2 != NULL) {
        _MainInit(self, argc, argv, argi, argv2);
    }

    return argv2;
}

static rc_t MainPrintConfig(const Main *self) {
    rc_t rc = 0;

    assert(self);

    if (rc == 0) {
        rc = KConfigPrint(self->cfg, 0);
        if (rc != 0) {
            OUTMSG(("KConfigPrint() = %R", rc));
        }
        OUTMSG(("\n"));
    }

    return rc;
}

static
rc_t _KDBPathTypePrint(const char *head, KPathType type, const char *tail)
{
    rc_t rc = 0;
    rc_t rc2 = 0;
    assert(head && tail);
    {
        rc_t rc2 = OUTMSG(("%s", head));
        if (rc == 0 && rc2 != 0) {
            rc = rc2;
        }
    }
    switch (type) {
        case kptNotFound:
            rc2 = OUTMSG(("NotFound"));
            break;
        case kptBadPath:
            rc2 = OUTMSG(("BadPath"));
            break;
        case kptFile:
            rc2 = OUTMSG(("File"));
            break;
        case kptDir:
            rc2 = OUTMSG(("Dir"));
            break;
        case kptCharDev:
            rc2 = OUTMSG(("CharDev"));
            break;
        case kptBlockDev:
            rc2 = OUTMSG(("BlockDev"));
            break;
        case kptFIFO:
            rc2 = OUTMSG(("FIFO"));
            break;
        case kptZombieFile:
            rc2 = OUTMSG(("ZombieFile"));
            break;
        case kptDataset:
            rc2 = OUTMSG(("Dataset"));
            break;
        case kptDatatype:
            rc2 = OUTMSG(("Datatype"));
            break;
        case kptDatabase:
            rc2 = OUTMSG(("Database"));
            break;
        case kptTable:
            rc2 = OUTMSG(("Table"));
            break;
        case kptIndex:
            rc2 = OUTMSG(("Index"));
            break;
        case kptColumn:
            rc2 = OUTMSG(("Column"));
            break;
        case kptMetadata:
            rc2 = OUTMSG(("Metadata"));
            break;
        case kptPrereleaseTbl:
            rc2 = OUTMSG(("PrereleaseTbl"));
            break;
        default:
            rc2 = OUTMSG(("unexpectedFileType(%d)", type));
            assert(0);
            break;
    }
    if (rc == 0 && rc2 != 0) {
        rc = rc2;
    }
    {
        rc_t rc2 = OUTMSG(("%s", tail));
        if (rc == 0 && rc2 != 0) {
            rc = rc2;
        }
    }
    return rc;
}

static bool isprintString(const unsigned char *s) {
    assert(s);

    while (*s) {
        int c = *(s++);
        if (!isprint(c)) {
            return false;
        }
    }

    return true;
}

static rc_t printString(const char *s) {
    rc_t rc = 0;

    const unsigned char *u = (unsigned char*)s;

    assert(u);

    if (isprintString(u)) {
        return OUTMSG(("%s", u));
    }

    while (*u) {
        rc_t rc2 = 0;
        int c = *(u++);
        if (isprint(c)) {
            rc2 = OUTMSG(("%c", c));
        }
        else {
            rc2 = OUTMSG(("\\%03o", c));
        }
        if (rc == 0 && rc2 != 0) {
            rc = rc2;
        }
    }

    return rc;
}

static rc_t _KDirectoryReport(const KDirectory *self,
    const char *name, int64_t *size, KPathType *type, bool *alias)
{
    rc_t rc = 0;
    const KFile *f = NULL;

    bool dummyB = false;
    int64_t dummy = 0;

    KPathType dummyT = kptNotFound;;
    if (type == NULL) {
        type = &dummyT;
    }

    if (alias == NULL) {
        alias = &dummyB;
    }
    if (size == NULL) {
        size = &dummy;
    }

    *type = KDirectoryPathType(self, name);

    if (*type & kptAlias) {
        OUTMSG(("alias|"));
        *type &= ~kptAlias;
        *alias = true;
    }

    rc = _KDBPathTypePrint("", *type, " ");

    if (*type == kptFile) {
        rc = KDirectoryOpenFileRead(self, &f, name);
        if (rc != 0) {
            OUTMSG(("KDirectoryOpenFileRead("));
            printString(name);
            OUTMSG((")=%R ", rc));
        }
        else {
            uint64_t sz = 0;
            rc = KFileSize(f, &sz);
            if (rc != 0) {
                OUTMSG(("KFileSize(%s)=%R ", name, rc));
            }
            else {
                OUTMSG(("%,lu ", sz));
                *size = sz;
            }
        }
    }
    else {
        OUTMSG(("- "));
    }

    RELEASE(KFile, f);

    return rc;
}

static rc_t _VDBManagerReport(const VDBManager *self,
    const char *name, KPathType *type)
{
    KPathType dummy = kptNotFound;;

    if (type == NULL) {
        type = &dummy;
    }

    *type = VDBManagerPathType(self, name);

    *type &= ~kptAlias;

    return _KDBPathTypePrint("", *type, " ");
}

static
rc_t _KDirectoryFileHeaderReport(const KDirectory *self, const char *path)
{
    rc_t rc = 0;
    char hdr[8] = "";
    const KFile *f = NULL;
    size_t num_read = 0;
    size_t i = 0;

    assert(self && path);

    rc = KDirectoryOpenFileRead(self, &f, path);
    if (rc != 0) {
        OUTMSG(("KDirectoryOpenFileRead(%s) = %R\n", path, rc));
        return rc;
    }

    rc = KFileReadAll(f, 0, hdr, sizeof hdr, &num_read);
    if (rc != 0) {
        OUTMSG(("KFileReadAll(%s, 8) = %R\n", path, rc));
    }

    for (i = 0; i < num_read && rc == 0; ++i) {
        if (isprint(hdr[i])) {
            OUTMSG(("%c", hdr[i]));
        }
        else {
            OUTMSG(("\\X%02X", hdr[i]));
        }
    }
    OUTMSG((" "));

    RELEASE(KFile, f);
    return rc;
}

static rc_t MainReport(const Main *self,
    const char *name, int64_t *size, KPathType *type, bool *alias)
{
    rc_t rc = 0;

    assert(self);

    rc = _KDirectoryReport(self->dir, name, size, type, alias);

    if (!self->noVDBManagerPathType) { /* && MainHasTest(self, eType)) { */
        _VDBManagerReport(self->mgr, name, type);
    }

    if (type != NULL && *type == kptFile) {
        _KDirectoryFileHeaderReport(self->dir, name);
    }

    return rc;
}

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
      ePathLocal
    , ePathRemote
    , ePathCache
} EPathType;
static rc_t MainPathReport(const Main *self, rc_t rc, const VPath *path,
    EPathType type, const char *name, const VPath* remote, int64_t *size,
    bool fasp, const KFile *fRemote)
{
    switch (type) {
        case ePathLocal:
            OUTMSG(("Local: "));
            break;
        case ePathRemote:
            OUTMSG(("Remote %s: ", fasp ? "fasp" : "http"));
            break;
        case ePathCache:
            OUTMSG(("Cache %s: ", fasp ? "fasp" : "http"));
            if (remote == NULL) {
                OUTMSG(("skipped\n"));
                return rc;
            }
            break;
    }
    if (rc != 0) {
        if (NotFoundByResolver(rc)) {
            OUTMSG(("not found\n"));
            rc = 0;
        }
        else {
            switch (type) {
                case ePathLocal:
                    OUTMSG(("VResolverLocal(%s) = %R\n", name, rc));
                    break;
                case ePathRemote:
                    OUTMSG(("VResolverRemote(%s) = %R\n", name, rc));
                    break;
                case ePathCache:
                    OUTMSG(("VResolverCache(%s) = %R\n", name, rc));
                    break;
            }
        }
    }
    else {
        const String *s = NULL;
        rc_t rc = VPathMakeString(path, &s);
        if (rc == 0) {
            OUTMSG(("%S ", s));
            switch (type) {
                case ePathLocal:
                case ePathCache:
                    rc = MainReport(self, s->addr, size, NULL, NULL);
                    break;
                case ePathRemote: {
                    uint64_t sz = 0;
                    if (!fasp && fRemote != NULL) {
                        rc = KFileSize(fRemote, &sz);
                        if (rc != 0) {
                            OUTMSG(("KFileSize(%s)=%R ", name, rc));
                        }
                        else {
                            OUTMSG(("%,lu ", sz));
                            *size = sz;
                        }
                    }
                    break;
                }
            }
            OUTMSG(("\n"));
        }
        else {
            const char *s = "";
            switch (type) {
                case ePathLocal:
                    s = "Local";
                    break;
                case ePathCache:
                    s = "Cache";
                    break;
                case ePathRemote:
                    s = "Remote";
                    break;
                default:
                    assert(0);
                    break;
            }
            OUTMSG(("VPathMakeUri(VResolver%s(%s)) = %R\n", s, name, rc));
        }
        if (type == ePathCache) {
            char cachecache[PATH_MAX] = "";
            if (rc == 0) {
                rc = string_printf(cachecache,
                    sizeof cachecache, NULL, "%s.cache", s->addr);
                if (rc != 0) {
                    OUTMSG(("string_printf(%s) = %R\n", s->addr, rc));
                }
            }

            if (rc == 0) {
                OUTMSG(("Cache.cache %s: ", fasp ? "fasp" : "http"));
                OUTMSG(("%s ", cachecache));
                rc = MainReport(self, cachecache, NULL, NULL, NULL);
                OUTMSG(("\n"));
            }
        }
        free((void*)s);
    }
    return rc;
}

static rc_t MainResolveLocal(const Main *self, const VResolver *resolver,
    const char *name, const VPath *acc, int64_t *size)
{
    rc_t rc = 0;

    const VPath* local = NULL;

    assert(self);

    if (resolver == NULL) {
        resolver = self->resolver;
    }
/*
    OUTMSG(("Local: "));*/

    rc = VResolverLocal(resolver, acc, &local);
    rc = MainPathReport(self,
        rc, local, ePathLocal, name, NULL, NULL, false, NULL);

    RELEASE(VPath, local);

    return rc;
}

static rc_t MainResolveRemote(const Main *self, VResolver *resolver,
    const char *name, const VPath* acc, const VPath **remote, int64_t *size,
    bool fasp)
{
    rc_t rc = 0;

    const KFile* f = NULL;

    assert(self && size);

    if (resolver == NULL) {
        resolver = self->resolver;
    }

    rc = VResolverRemote(resolver,
        fasp ? eProtocolFaspHttp : eProtocolHttp, acc, remote, &f);
    rc = MainPathReport(self,
        rc, *remote, ePathRemote, name, NULL, size, fasp, f);
    RELEASE(KFile, f);

    return rc;
}

static rc_t MainResolveCache(const Main *self, const VResolver *resolver,
    const char *name, const VPath* remote, bool fasp)
{
    rc_t rc = 0;
    const VPath* cache = NULL;

    assert(self);
    
    if (remote == NULL) {
        rc = MainPathReport(self,
            rc, cache, ePathCache, name, remote, NULL, fasp, NULL);
    }
    else {
        uint64_t file_size = 0;

        if (resolver == NULL) {
            resolver = self->resolver;
        }

        if (!self->allowCaching) {
            VResolverCacheEnable(resolver, self->cacheState);
        }
        rc = VResolverCache(resolver, remote, &cache, file_size);
        rc = MainPathReport(self,
            rc, cache, ePathCache, name, remote, NULL, fasp, NULL);

        RELEASE(VPath, cache);
        if (!self->allowCaching) {
            VResolverCacheEnable(resolver, vrAlwaysDisable);
        }
    }

    return rc;
}

typedef enum {
      eQueryAll
    , eQueryLocal
    , eQueryRemote
} EQueryType ;
static rc_t VResolverQueryByType(const Main *self, const VResolver *resolver,
    const char *name, const VPath *query,
    bool fasp, VRemoteProtocols protocols, EQueryType type)
{
    rc_t rc = 0;
    rc_t rc2 = 0;
    const VPath *local = NULL;
    const VPath *remote = NULL;
    const VPath *cache = NULL;
    const VPath **pLocal = NULL;
    const VPath **pRemote = NULL;
    const VPath **pCache = NULL;

    switch (type) {
        case eQueryLocal:
        case eQueryAll:
            pLocal = &local;
        default:
            break;
    }

    switch (type) {
        case eQueryRemote:
        case eQueryAll:
            pRemote = &remote;
            pCache = &cache;
        default:
            break;
    }

    if (pCache != NULL && !self->allowCaching) {
        VResolverCacheEnable(resolver, vrAlwaysEnable);
    }

    rc = VResolverQuery(resolver, protocols, query, pLocal, pRemote, pCache);
    OUTMSG(("\nVResolverQuery(%s, %s, local%s, remote%s, cache%s)= %R\n",
        name, protocols == eProtocolHttp ? "Http" : "FaspHttp", 
        pLocal == NULL ? "=NULL" : "", pRemote == NULL ? "=NULL" : "",
        pCache == NULL ? "=NULL" : "", rc));
    if (rc == 0) {
        if (local != NULL) {
            rc2 = MainPathReport(self,
                rc, local, ePathLocal, name, NULL, NULL, false, NULL);
        }
        if (remote != NULL) {
            rc2 = MainPathReport(self,
                rc, remote, ePathRemote, name, NULL, NULL, fasp, NULL);
        }
        if (cache != NULL) {
            rc2 = MainPathReport(self,
                rc, cache, ePathCache, name, remote, NULL, fasp, NULL);
        }
    }

    RELEASE(VPath, local);
    RELEASE(VPath, remote);
    RELEASE(VPath, cache);

    if (pCache != NULL && !self->allowCaching) {
        VResolverCacheEnable(resolver, self->cacheState);
    }

    return rc;
}

static rc_t MainResolveQuery(const Main *self, const VResolver *resolver,
    const char *name, const VPath *query, bool fasp)
{
    rc_t rc = 0;
    rc_t rc2 = 0;

    VRemoteProtocols protocols = eProtocolHttp;
    if (fasp) {
        protocols = eProtocolFaspHttp;
    }
    else {
        protocols = eProtocolHttp;
    }

    if (resolver == NULL) {
        resolver = self->resolver;
    }

    rc2 = VResolverQueryByType(self, resolver, name, query, fasp, protocols,
        eQueryAll);
    if (rc2 != 0 && rc == 0) {
        rc = rc2;
    }

    rc2 = VResolverQueryByType(self, resolver, name, query, fasp, protocols,
        eQueryLocal);
    if (rc2 != 0 && rc == 0) {
        rc = rc2;
    }

    rc2 = VResolverQueryByType(self, resolver, name, query, fasp, protocols,
        eQueryRemote);
    if (rc2 != 0 && rc == 0) {
        rc = rc2;
    }

    return rc;
}

static rc_t MainResolve(const Main *self, const KartItem *item,
    const char *name, int64_t *localSz, int64_t *remoteSz)
{
    rc_t rc = 0;

    VPath* acc = NULL;
    VResolver* resolver = NULL;

    assert(self);

    if (rc == 0) {
        VFSManager *mgr = NULL;
        rc = VFSManagerMake(& mgr);
        if (rc != 0) {
            OUTMSG(("VFSManagerMake = %R\n", rc));
        }
        else {
            if (item == NULL) {
                rc = VFSManagerMakePath(mgr, &acc, name);
                if (rc != 0) {
                    OUTMSG(("VFSManagerMakePath(%s) = %R\n", name, rc));
                }
            }
            else {
                const KRepository *p_protected = NULL;
                uint64_t oid = 0;
                uint64_t project = 0;
                if (rc == 0) {
                    rc = KartItemProjIdNumber(item, &project);
                    if (rc != 0) {
                        OUTMSG(("KartItemProjectIdNumber = %R\n", rc));
                    }
                }
                if (rc == 0) {
                    rc = KartItemItemIdNumber(item, &oid);
                    if (rc != 0) {
                        OUTMSG(("KartItemItemIdNumber = %R\n", rc));
                    }
                }
                if (rc == 0) {
                    rc = VFSManagerMakeOidPath(mgr, &acc, oid);
                    if (rc != 0) {
                        OUTMSG(("VFSManagerMakePath(%d) = %R\n", oid, rc));
                    }
                }
                if (rc == 0) {
                    rc = KRepositoryMgrGetProtectedRepository(self->repoMgr, 
                        project, &p_protected);
                    if (rc != 0) {
                        OUTMSG((
                            "KRepositoryMgrGetProtectedRepository(%d) = %R\n",
                            project, rc));
                    }
                }
                if (rc == 0) {
                    rc = KRepositoryMakeResolver(p_protected, &resolver,
                        self->cfg);
                    if (rc != 0) {
                        OUTMSG((
                            "KRepositoryMakeResolver(%d) = %R\n", project, rc));
                    }
                }
                RELEASE(KRepository, p_protected);
            }
            RELEASE(VFSManager, mgr);
        }
    }

    if (rc == 0) {
        const VPath* remote = NULL;

        rc_t rc2 = MainResolveLocal(self, resolver, name, acc, localSz);
        if (rc2 != 0 && rc == 0) {
            rc = rc2;
        }

        rc2 = MainResolveRemote(self, resolver, name, acc, &remote, remoteSz,
            false);
        if (rc2 != 0 && rc == 0) {
            rc = rc2;
        }

        rc2 = MainResolveCache(self, resolver, name, remote, false);
        if (rc2 != 0 && rc == 0) {
            rc = rc2;
        }

        rc2 = MainResolveRemote(self, resolver, name, acc, &remote, remoteSz,
            true);
        if (rc2 != 0 && rc == 0) {
            rc = rc2;
        }

        rc2 = MainResolveCache(self, resolver, name, remote, true);
        if (rc2 != 0 && rc == 0) {
            rc = rc2;
        }

        rc2 = MainResolveQuery(self, resolver, name, acc, false);
        if (rc2 != 0 && rc == 0) {
            rc = rc2;
        }

        rc2 = MainResolveQuery(self, resolver, name, acc, true);
        if (rc2 != 0 && rc == 0) {
            rc = rc2;
        }

        RELEASE(VPath, remote);
    }

    RELEASE(VPath, acc); 
    RELEASE(VResolver, resolver);

    return rc;
}

static
rc_t MainDepend(const Main *self, const char *name, bool missing)
{
    rc_t rc = 0;

    const VDatabase *db = NULL;
    const VDBDependencies* dep = NULL;
    uint32_t count = 0;

    if (rc == 0) {
        rc = VDBManagerOpenDBRead(self->mgr, &db, NULL, name);
        if (rc != 0) {
            if (rc == SILENT_RC(rcVFS,rcMgr,rcOpening,rcDirectory,rcNotFound)) {
                return 0;
            }
            OUTMSG(("VDBManagerOpenDBRead(%s) = %R\n", name, rc));
        }
    }

    if (rc == 0) {
        if (!self->allowCaching) {
            VResolverCacheEnable(self->resolver, self->cacheState);
        }
        rc = VDatabaseListDependenciesWithCaching(db,
            &dep, missing, !self->allowCaching);
        if (rc != 0) {
            OUTMSG(("VDatabaseListDependencies(%s, %s) = %R\n",
                name, missing ? "missing" : "all", rc));
        }
        if (!self->allowCaching) {
            VResolverCacheEnable(self->resolver, vrAlwaysDisable);
        }
    }

    if (rc == 0) {
        rc = VDBDependenciesCount(dep, &count);
        if (rc != 0) {
            OUTMSG(("VDBDependenciesCount(%s, %s) = %R\n",
                name, missing ? "missing" : "all", rc));
        }
        else {
            OUTMSG(("VDBDependenciesCount(%s)=%d\n",
                missing ? "missing" : "all", count));
        }
    }

    if (rc == 0) {
        uint32_t i = 0;
        rc_t rc2 = 0;
        for (i = 0; i < count; ++i) {
            bool b = true;
            const char *s = NULL;
            KPathType type = kptNotFound;

            OUTMSG((" %6d\t", i + 1));

            rc2 = VDBDependenciesSeqId(dep, &s, i);
            if (rc2 == 0) {
                OUTMSG(("seqId=%s,", s));
            }
            else {
                OUTMSG(("VDBDependenciesSeqId(%s, %s, %i)=%R ",
                    name, missing ? "missing" : "all", i, rc2));
                if (rc == 0) {
                    rc = rc2;
                }
            }

            rc2 = VDBDependenciesName(dep, &s, i);
            if (rc2 == 0) {
                OUTMSG(("name=%s,", s));
            }
            else {
                OUTMSG(("VDBDependenciesName(%s, %s, %i)=%R ",
                    name, missing ? "missing" : "all", i, rc2));
                if (rc == 0) {
                    rc = rc2;
                }
            }

            rc2 = VDBDependenciesCircular(dep, &b, i);
            if (rc2 == 0) {
                OUTMSG(("circular=%s,", b ? "true" : "false"));
            }
            else {
                OUTMSG(("VDBDependenciesCircular(%s, %s, %i)=%R ",
                    name, missing ? "missing" : "all", i, rc2));
                if (rc == 0) {
                    rc = rc2;
                }
            }

            rc2 = VDBDependenciesType(dep, &type, i);
            if (rc2 == 0) {
                rc2 = _KDBPathTypePrint("type=", type, ",");
                if (rc2 != 0 && rc == 0) {
                    rc = rc2;
                }
            }
            else {
                OUTMSG(("VDBDependenciesType(%s, %s, %i)=%R ",
                    name, missing ? "missing" : "all", i, rc2));
                if (rc == 0) {
                    rc = rc2;
                }
            }

            rc2 = VDBDependenciesLocal(dep, &b, i);
            if (rc2 == 0) {
                OUTMSG(("local=%s,", b ? "local" : "remote"));
            }
            else {
                OUTMSG(("VDBDependenciesLocal(%s, %s, %i)=%R ",
                    name, missing ? "missing" : "all", i, rc2));
                if (rc == 0) {
                    rc = rc2;
                }
            }

            rc2 = VDBDependenciesPath(dep, &s, i);
            if (rc2 == 0) {
                OUTMSG(("\n\tpathLocal=%s,", s == NULL ? "notFound" : s));
            }
            else {
                OUTMSG(("VDBDependenciesPath(%s, %s, %i)=%R ",
                    name, missing ? "missing" : "all", i, rc2));
                if (rc == 0) {
                    rc = rc2;
                }
            }

            rc2 = VDBDependenciesPathRemote(dep, &s, i);
            if (rc2 == 0) {
                if (s == NULL) {
                    OUTMSG(("\n\tpathRemote: notFound "));
                }
                else {
                    OUTMSG(("\n\tpathRemote: %s ", s));
                    if (!self->noRfs) {
                        const KFile *f = NULL;
                        rc2 = KCurlFileMake(&f, s, false);
                        if (rc2 != 0) {
                            OUTMSG(("KCurlFileMake=%R", rc2));
                            if (rc == 0) {
                                rc = rc2;
                            }
                        }
                        if (rc2 == 0) {
                            uint64_t sz = 0;
                            rc2 = KFileSize(f, &sz);
                            if (rc2 != 0) {
                                OUTMSG(("KFileSize=%R", rc2));
                                if (rc == 0) {
                                    rc = rc2;
                                }
                            }
                            else {
                                OUTMSG(("%,lu", sz));
                            }
                        }
                        RELEASE(KFile, f);
                    }
                }
            }
            else {
                OUTMSG(("VDBDependenciesPathRemote(%s, %s, %i)=%R ",
                    name, missing ? "missing" : "all", i, rc2));
                if (rc == 0) {
                    rc = rc2;
                }
            }

            rc2 = VDBDependenciesPathCache(dep, &s, i);
            if (rc2 == 0) {
                OUTMSG(("\n\tpathCache: %s ", s == NULL ? "notFound" : s));
                if (s != NULL) {
                    char cachecache[PATH_MAX] = "";
                    rc2 = MainReport(self, s, NULL, NULL, NULL);
                    OUTMSG(("\n"));
                    if (rc == 0) {
                        rc = rc2;
                    }
                    if (rc2 == 0) {
                        rc2 = string_printf(cachecache,
                            sizeof cachecache, NULL, "%s.cache", s);
                        if (rc2 != 0) {
                            if (rc == 0) {
                                rc = rc2;
                            }
                            OUTMSG(("string_printf(%s) = %R\n", s, rc2));
                        }
                    }
                    if (rc == 0) {
                        OUTMSG(("\tpathCache.cache: "));
                        OUTMSG(("%s ", cachecache));
                        rc = MainReport(self, cachecache, NULL, NULL, NULL);
                    }
                }
            }
            else {
                OUTMSG(("VDBDependenciesPathCache(%s, %s, %i)=%R ",
                    name, missing ? "missing" : "all", i, rc2));
                if (rc == 0) {
                    rc = rc2;
                }
            }

            OUTMSG(("\n"));
        }
    }

    RELEASE(VDBDependencies, dep);
    RELEASE(VDatabase, db);

    return rc;
}

static rc_t PrintCurl() {
    KNSManager *mgr = NULL;

    rc_t rc = KNSManagerMake(&mgr);
    if (rc != 0) {
        OUTMSG(("KNSManagerMake = %R\n", rc));
    }

    if (rc == 0) {
        rc_t rc = KNSManagerAvail(mgr);
        OUTMSG(("KNSManagerAvail = %R", rc));

        if (rc == 0) {
            const char *version_string = NULL;
            rc = KNSManagerCurlVersion(mgr, &version_string);
            if (rc == 0) {
                OUTMSG((". Curl Version = %s\n", version_string));
            }
            else {
                OUTMSG((". KNSManagerCurlVersion = %R\n", rc));
            }
        }
    }

    RELEASE(KNSManager, mgr);

    return rc;
}

#define kptKartITEM (kptAlias - 1)

static
rc_t MainExec(const Main *self, const KartItem *item, const char *aArg, ...)
{
    rc_t rc = 0;
    rc_t rce = 0;

    KPathType type = kptNotFound;
    bool alias = false;
    int64_t directSz = -1;
    int64_t localSz = -1;
    int64_t remoteSz = -1;
    size_t num_writ = 0;
    char arg[PATH_MAX] = "";

    va_list args;
    va_start(args, aArg);

    assert(self);

    if (item != NULL) {
        type = kptKartITEM;
    }

    else {
        rc = string_vprintf(arg, sizeof arg, &num_writ, aArg, args);
        if (rc != 0) {
            OUTMSG(("string_vprintf(%s)=%R\n", aArg, rc));
            return rc;
        }
        assert(num_writ < sizeof arg);

        OUTMSG(("\n"));
        rc = printString(arg);
        if (rc != 0) {
            OUTMSG(("printString=%R\n", rc));
            return rc;
        }
        OUTMSG((" "));
        rc = MainReport(self, arg, &directSz, &type, &alias);
        OUTMSG(("\n"));
    }

    if (self->recursive && type == kptDir && !alias) {
        uint32_t i = 0;
        uint32_t count = 0;
        KNamelist *list = NULL;
        rc = KDirectoryList(self->dir, &list, NULL, NULL, arg);
        if (rc != 0) {
            OUTMSG(("KDirectoryList(%s)=%R ", arg, rc));
        }
        else {
            rc = KNamelistCount(list, &count);
            if (rc != 0) {
                OUTMSG(("KNamelistCount(KDirectoryList(%s))=%R ", arg, rc));
            }
        }
        for (i = 0; i < count && rc == 0; ++i) {
            const char *name = NULL;
            rc = KNamelistGet(list, i, &name);
            if (rc != 0) {
                OUTMSG(("KNamelistGet(KDirectoryList(%s), %d)=%R ",
                    arg, i, rc));
            }
            else {
                rc_t rc2 = MainExec(self, NULL, "%s/%s", arg, name);
                if (rc2 != 0 && rce == 0) {
                    rce = rc2;
                }
            }
        }
        RELEASE(KNamelist, list);
    }
    else {
        bool isKart = false;
        Kart *kart = NULL;
        if (type == kptFile) {
            rc_t rc = KartMake(self->dir, arg, &kart, &isKart);
            if (rc != 0) {
                OUTMSG(("KartMake = %R\n", rc));
            }
        }

        if (isKart) {
            const KartItem *item = NULL;
            while (true) {
                rc_t rc2 = 0;
                RELEASE(KartItem, item);
                rc2 = KartMakeNextItem(kart, &item);
                if (rc2 != 0) {
                    OUTMSG(("KartMakeNextItem = %R\n", rc2));
                    if (rce == 0) {
                        rce = rc2;
                    }
                    break;
                }
                if (item == NULL) {
                    break;
                }
                rc2 = MainExec(self, item, NULL);
                if (rc2 != 0 && rce == 0) {
                    rce = rc2;
                }
            }
            KartRelease(kart);
            kart = NULL;
        }
        else {
            if (MainHasTest(self, eResolve)) {
                rc_t rc2 = MainResolve(self, item, arg, &localSz, &remoteSz);
                if (rc == 0 && rc2 != 0) {
                    rc = rc2;
                }
            }

            if (item == NULL) {
                if (type == kptDatabase || type == kptNotFound) {
                    if (MainHasTest(self, eDependMissing)) {
                        rc_t rc2 = MainDepend(self, arg, true);
                        if (rc == 0 && rc2 != 0) {
                            rc = rc2;
                        }
                    }

                    if (MainHasTest(self, eDependAll)) {
                        rc_t rc2 = MainDepend(self, arg, false);
                        if (rc == 0 && rc2 != 0) {
                            rc = rc2;
                        }
                    }
                }
            }

            if (MainHasTest(self, eResolve) && (
                (directSz != -1 && localSz != -1 && directSz != localSz) ||
                (remoteSz != -1 && localSz != -1 && localSz != remoteSz))
               )
            {
                OUTMSG(("FILE SIZES DO NOT MATCH: "));
                if (directSz != -1 && localSz != -1 && directSz != remoteSz)
                {
                    OUTMSG(("direct=%ld != remote=%ld. ", directSz, remoteSz));
                }
                if (remoteSz != -1 && localSz != -1 && localSz != remoteSz) {
                    OUTMSG(("local=%ld != remote=%ld. ", localSz, remoteSz));
                }
                OUTMSG(("\n"));
            }

            OUTMSG(("\n"));
        }
    }

    if (rce != 0 && rc == 0) {
        rc = rce;
    }
    return rc;
}

static rc_t MainFini(Main *self) {
    rc_t rc = 0;

    assert(self);

    RELEASE(VResolver, self->resolver);
    RELEASE(KConfig, self->cfg);
    RELEASE(KRepositoryMgr, self->repoMgr);
    RELEASE(VDBManager, self->mgr);
    RELEASE(KDirectory, self->dir);

    return rc;
}

#define OPTION_CACHE "allow-caching"
#define ALIAS_CACHE  "C"
static const char* USAGE_CACHE[] = { "do not disable caching", NULL };

#define OPTION_NO_RFS "no-rfs"
static const char* USAGE_NO_RFS[]
    = { "do not check remote file size for dependencies", NULL };

#define OPTION_NO_VDB "no-vdb"
#define ALIAS_NO_VDB  "N"
static const char* USAGE_NO_VDB[] = { "do not call VDBManagerPathType", NULL };

#define OPTION_REC "recursive"
#define ALIAS_REC  "R"
static const char* USAGE_REC[] = { "check object type recursively", NULL };

OptDef Options[] = {                             /* needs_value, required */
    { OPTION_CACHE , ALIAS_CACHE , NULL, USAGE_CACHE , 1, false, false },
    { OPTION_NO_RFS, NULL        , NULL, USAGE_NO_RFS, 1, false, false },
    { OPTION_NO_VDB, ALIAS_NO_VDB, NULL, USAGE_NO_VDB, 1, false, false },
    { OPTION_REC   , ALIAS_REC   , NULL, USAGE_REC   , 1, false, false }
};

rc_t CC KMain(int argc, char *argv[]) {
    rc_t rc = 0;
    uint32_t pcount = 0;
    uint32_t i = 0;
    Args *args = NULL;
    rc_t rc3 = 0;
    int argi = 0;

    Main prms;
    char **argv2 = MainInit(&prms, argc, argv, &argi);

    if (rc == 0) {
        rc = ArgsMakeAndHandle(&args, argi, argv2, 1,
            Options, sizeof Options / sizeof Options[0]);
    }

    if (rc == 0) {
        rc = ArgsOptionCount(args, OPTION_CACHE, &pcount);
        if (rc) {
            LOGERR(klogErr, rc, "Failure to get '" OPTION_CACHE "' argument");
        }
        else {
            if (pcount > 0) {
                prms.allowCaching = true;
            }
        }
    }

    if (rc == 0) {
        rc = MainInitObjects(&prms);
    }

    if (MainHasTest(&prms, eCfg)) {
        rc_t rc2 = MainPrintConfig(&prms);
        if (rc == 0 && rc2 != 0) {
            rc = rc2;
        }
    }

    if (MainHasTest(&prms, eCurl)) {
        PrintCurl();
    }

    if (rc == 0) {
        rc = ArgsOptionCount(args, OPTION_NO_RFS, &pcount);
        if (rc) {
            LOGERR(klogErr, rc, "Failure to get '" OPTION_NO_RFS "' argument");
        }
        else {
            if (pcount > 0) {
                prms.noRfs = true;
            }
        }
    }

    if (rc == 0) {
        rc = ArgsOptionCount(args, OPTION_NO_VDB, &pcount);
        if (rc) {
            LOGERR(klogErr, rc, "Failure to get '" OPTION_NO_VDB "' argument");
        }
        else {
            if (pcount > 0) {
                prms.noVDBManagerPathType = true;
            }
        }
    }

    if (rc == 0) {
        rc = ArgsOptionCount(args, OPTION_REC, &pcount);
        if (rc) {
            LOGERR(klogErr, rc, "Failure to get '" OPTION_REC "' argument");
        }
        else {
            if (pcount > 0) {
                prms.recursive = true;
            }
        }
    }

    if (rc == 0) {
        rc = ArgsParamCount(args, &pcount);
    }

    for (i = 0; i < pcount; ++i) {
        const char *name = NULL;
        rc3 = ArgsParamValue(args, i, &name);
        if (rc3 == 0) {
            rc_t rc2 = MainExec(&prms, NULL, name);
            if (rc == 0 && rc2 != 0) {
                rc = rc2;
            }
        }
    }
    if (rc == 0 && rc3 != 0) {
        rc = rc3;
    }

    RELEASE(Args, args);

    {
        rc_t rc2 = MainFini(&prms);
        if (rc == 0 && rc2 != 0) {
            rc = rc2;
        }
    }
    free(argv2);

    OUTMSG(("\n"));

    return rc;
}
