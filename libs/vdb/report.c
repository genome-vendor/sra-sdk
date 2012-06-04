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

#include <klib/report.h> /* ReportInit */
#include <vdb/report.h> /* ... */
#include <klib/namelist.h> /* KNamelistRelease */
#include <klib/out.h> /* OUTMSG */
#include <klib/klib-priv.h> /* ReportFuncs */

#include <vdb/manager.h> /* VDBManagerVersion */
#include <vdb/vdb-priv.h> /* VTableOpenKTableRead */
#include <vdb/database.h> /* VDatabaseRelease */
#include <vdb/table.h> /* VTableOpenParentRead */
#include <vdb/dependencies.h> /* VDBDependenciesRelease */

#include <kdb/kdb-priv.h> /* KDatabaseGetPath */
#include <kdb/manager.h> /* kptDatabase */
#include <kdb/database.h> /* KDatabaseRelease */
#include <kdb/table.h> /* KTableRelease */

#include <kfs/dyload.h> /* KDyld */
#include <kfs/file.h> /* KFileRead */
#include <kfs/nullfile.h> /* KFileMakeNullUpdate */
#include <kfs/md5.h> /* KMD5SumFmt */

#include <atomic.h>

#include <stdarg.h> /* va_start */
#include <stdio.h> /* sprintf */
#include <stdlib.h> /* malloc */
#include <string.h> /* memset */
#include <limits.h> /* PATH_MAX */
#include <assert.h>

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

#define RELEASE(type, obj) do { rc_t rc2 = type##Release(obj); \
    if (rc2 && !rc) { rc = rc2; } obj = NULL; } while (false)

/*
 * An unrecoverable error happened.
 * We can help to solve it
 * by reporting information about known application execution environment.
 */

typedef struct Report {
    const VDBManager* mgr;
    const VDatabase* db;
    const VTable* table;
} Report;

static Report *report_singleton;

static rc_t CC ReportRelease ( void )
{
    rc_t rc = 0;
    Report *prev_report, *cur_report;

    cur_report = report_singleton;
    do
    {
        prev_report = cur_report;
        cur_report = atomic_test_and_set_ptr ( ( void* volatile* ) & report_singleton, NULL, prev_report );
    }
    while ( cur_report != prev_report );

    if ( cur_report != NULL )
    {
        /* cleanup */
        VTableRelease ( cur_report -> table );
        VDatabaseRelease ( cur_report -> db );
        VDBManagerRelease ( cur_report -> mgr );
        memset ( cur_report, 0, sizeof * cur_report );
    }

    return rc;
}

static rc_t CC ReportObj ( const ReportFuncs *f, uint32_t indent, const char *path );
static rc_t CC ReportSOFTWARE ( const ReportFuncs *f, uint32_t indent, const char *argv_0, const char *date, ver_t tool_ver );

static rc_t ReportGet(Report** self)
{
    rc_t rc = 0;

    static bool latch;
    if ( ! latch )
    {
        static Report self;

        rc = ReportInitVDB ( ReportObj, ReportSOFTWARE, ReportRelease );
        if ( rc == 0 )
        {
            report_singleton = & self;
            latch = true;
        }
    }

    * self = ( Report* ) report_singleton;

    return rc;
}

#define report ( * f -> report )
#define reportData ( * f -> reportData )
#define reportData1 ( * f -> reportData1 )
#define reportOpen ( * f -> reportOpen )
#define reportOpen1 ( * f -> reportOpen1 )
#define reportClose ( * f -> reportClose )
#define reportClose1 ( * f -> reportClose1 )
#define reportError ( * f -> reportError )
#define reportErrorStr ( * f -> reportErrorStr )
#define reportErrorStrImpl ( * f -> reportErrorStrImpl )
#define reportErrorStrInt ( * f -> reportErrorStrInt )
#define reportError3Str ( * f -> reportError3Str )

static
rc_t ReportBuild(const ReportFuncs *f, uint32_t indent, bool dynamic, const VDBManager* mgr)
{
    rc_t rc = 0;
    if (!dynamic)
    {   report(indent, "Build", 1, "static", 's', "true"); }
    else {
        KNamelist* list = NULL;
        reportOpen(indent, "Build", 1, "static", 's', "false");
        if (mgr) {
            rc_t rc2 = VDBManagerListExternalSchemaModules(mgr, &list);
            if (rc2 != 0) {
                reportError
                    (indent + 1, rc2, "VDBManagerListExternalSchemaModules");
                if (rc == 0 && rc2 != 0)
                {   rc = rc2; }
            }
            else {
                uint32_t count = 0;
                rc2 = KNamelistCount(list, &count);
                if (rc2 != 0) {
                    reportErrorStr(indent + 1, rc2, "KNamelistCount", "origin",
                        "VDBManagerListExternalSchemaModules");
                    if (rc == 0 && rc2 != 0)
                    {   rc = rc2; }
                }
                else {
                    int64_t i = 0;
                    for (i = 0; i < count && rc2 == 0; ++i) {
                        const char* name = NULL;
                        rc2 = KNamelistGet(list, i, &name);
                        if (rc2 != 0) {
                            reportErrorStr(
                                indent + 1, rc2, "KNamelistGet", "origin",
                                "VDBManagerListExternalSchemaModules");
                            if (rc == 0 && rc2 != 0)
                            {   rc = rc2; }
                        }
                        else {
                            report(indent + 1, "Module", 1, "name", 's', name);
                        }
                    }
                }
            }
        }
        RELEASE(KNamelist, list);
        reportClose(indent, "Build");
    }
    return rc;
}

static rc_t ReportDepend( const ReportFuncs *f, uint32_t indent, const VDatabase* db) {
    rc_t rc = 0;

    const VDBDependencies* dep = NULL;

    const char tag[] = "Dependencies";

    assert(db);

    reportOpen(indent, tag, 0);

    rc = VDatabaseListDependencies(db, &dep, true);
    if (rc != 0) {
        reportError(indent + 1, rc, "VDatabaseListDependencies");
    }
    else {
        uint32_t count = 0;
        rc = VDBDependenciesCount(dep, &count);
        if (rc != 0) { 
            reportError(indent + 1, rc, "VDBDependenciesCount");
        }
        else {
            uint32_t i = 0;
            const char tag[] = "Missing";
            reportOpen(indent + 1, tag, 1, "count", 'd', count);
            for (i = 0; i < count; ++i) {
                const char* seq_id = "";
                rc = VDBDependenciesSeqId(dep, &seq_id, i);
                if (rc != 0) {
                    reportErrorStrInt(indent + 2, rc, "VDBDependenciesSeqId",
                        "origin", "VDatabaseListDependencies", "idx", i);
                }
                else {
                    report(indent + 2, "Dependency", 2,
                        "index", 'd', i, "seq_id", 's', seq_id);
                }
            }
            reportClose(indent + 1, tag);
        }
    }

    reportClose(indent, tag);

    RELEASE(VDBDependencies, dep);

    return rc;
}

typedef struct Total {
    int64_t sz;
    int64_t files;
} Total;

static rc_t CC visitor(const KDirectory* dir,
    uint32_t type, const char* name, void* data)
{
    rc_t rc = 0;
    Total* total = (Total*) data;
    if (type & kptAlias)
    {   return rc; }
    assert(total);
    switch (type) {
        case kptFile: {
            uint64_t size = 0;
            rc = KDirectoryFileSize(dir, &size, name);
            if (rc == 0) {
                total->sz += size;
            }
            ++total->files;
            break;
        }
        case kptDir: 
            rc = KDirectoryVisit(dir, false, visitor, total, name);
            break;
        default:
            rc = RC(rcExe, rcDirectory, rcVisiting, rcType, rcUnexpected);
            break;
    }
    return rc;
}

static rc_t ReportDir(const ReportFuncs *f, uint32_t indent, const KTable* tbl) {
    rc_t rc = 0;
    const KDirectory* dir = NULL;
    if (tbl == NULL) {
        report(indent, "Error", 1, "KTable" , 's', "NULL");
        return rc;
    }
    rc = KTableOpenDirectoryRead(tbl, &dir);
    if (rc != 0) {
        reportError(indent, rc, "KTableOpenDirectoryRead");
    }
    else {
        Total total;
        memset(&total, 0, sizeof total);
        rc = KDirectoryVisit(dir, false, visitor, &total, NULL);
        report(indent, "Directory", 2,
            "size", 'l', total.sz, "files", 'l', total.files);
    }
    RELEASE(KDirectory, dir);
    return rc;
}

#define OBJ_OPEN(indent,count,path,type) reportOpen(indent, "Object", count, \
    "path", 's', path, "type", 's', type
#define OBJ_P_OPEN(indent,count,path,type,file_type) \
    OBJ_OPEN(indent, count, path, type), "fs_type", 's', file_type

#define OBJ(indent,path,type) \
    OBJ_OPEN(indent, 2, path, type) )
#define OBJ_P(indent,path,type,file_type) \
    OBJ_P_OPEN(indent, 3, path, type, file_type) )
#define OBJ_P_A(indent,path,type,file_type) \
    OBJ_P_OPEN(indent, 4, path, type, file_type), "alias", 's', "true")
#define OBJ_P_S(indent,path,type,file_type,size) \
    OBJ_P_OPEN(indent, 4, path, type, file_type), "size", 'l', size)
#define OBJ_P_S_A(indent,path,type,file_type,size) \
    OBJ_P_OPEN(indent, 5, path, type, file_type), "size", 'l', size, \
                                                  "alias", 's', "true")

static rc_t CC ReportObj(const ReportFuncs *f, uint32_t indent, const char *object) {
    Report* self = NULL;
    const char* fullpath = NULL;
    const KDatabase* kdb = NULL;
    const KTable* ktbl = NULL;
    const VDatabase* db = NULL;
    KPathType type = kptNotFound;
    KPathType file_type = kptNotFound;
    bool alias = false;
    uint64_t size = 0;
    bool size_unknown = true;

    rc_t rc = ReportGet(&self);
    assert(self);

    if (self->db) {
        type = kptDatabase;
        db = self->db;
    }
    else if (self->table)
    {
        rc_t rc2 = VTableOpenParentRead(self->table, &db);
        if (rc2)
        {
            if (rc == 0)
            {
                rc = rc2;
            }
        }
        else if (!db)
        {
            type = kptTable;
            rc2 = VTableGetKTableRead(self->table, &ktbl);
            if (rc2)
            {
                if (rc == 0)
                {
                    rc = rc2;
                }
            }
            else
            {
                rc2 = KTableGetPath(ktbl, &fullpath);
            }
        }
    }

    if (db) {
        rc_t rc2 = VDatabaseOpenKDatabaseRead(db, &kdb);
        type = kptDatabase;
        if (rc2) {
            if (rc == 0)
            {   rc = rc2; }
        }
        else {
            rc2 = KDatabaseGetPath(kdb, &fullpath);
            if (rc2) {
                if (rc == 0)
                {   rc = rc2; }
            }
        }
    }

    if (fullpath) {
        KDirectory* dir = NULL;
        rc_t rc2 = KDirectoryNativeDir(&dir);
        if (rc2) {
            if (rc == 0)
            {   rc = rc2; }
        }
        else {
            file_type = KDirectoryPathType(dir, fullpath);
            alias = file_type & kptAlias;
            file_type &= ~kptAlias;
            if (file_type == kptFile) {
                rc2 = KDirectoryFileSize(dir, &size, fullpath);
                if (rc2) {
                    if (rc == 0)
                    {   rc = rc2; }
                }
                else {  size_unknown = false; }
            }
        }
        RELEASE(KDirectory, dir);
    }

    if (object || type != kptNotFound) {
        const char* path
            = fullpath ? fullpath : object ? object : "not set";
        const char* stype = type == kptTable ? "table" : 
            type == kptDatabase ? "database" : "unknown";
        const char* sfile_type = file_type == kptFile ? "archive" : 
            file_type == kptDir ? "dir" : "unexpected";

        if (fullpath && !size_unknown) {
            if (alias)
            { OBJ_P_S_A(indent, path, stype, sfile_type, size); }
            else
            { OBJ_P_S  (indent, path, stype, sfile_type, size); }
        }
        else if (fullpath && size_unknown) {
            if (alias)
            { OBJ_P_A  (indent, path, stype, sfile_type); }
            else
            { OBJ_P    (indent, path, stype, sfile_type); }
        }
        else
        {     OBJ      (indent, path, stype); }

        if (!db)
        {   db = self->db; }

        if (db) {
            rc_t rc2 = ReportDepend(f, indent + 1, db);
            if (rc == 0)
            {   rc = rc2; }
        }
        if (file_type == kptDir) {
            rc_t rc2 = ReportDir(f, indent + 1, ktbl);
            if (rc == 0)
            {   rc = rc2; }
        }

        reportClose(indent, "Object");
    }

    if (db != self->db)
    {   RELEASE(VDatabase, db); }
    RELEASE(KTable, ktbl);
    RELEASE(KDatabase, kdb);

    return rc;
}

#ifdef _STATIC

static
rc_t md5(const char* path, uint8_t digest[16], const KDirectory* dir)
{
    const KFile* kf = NULL;
    rc_t rc = KDirectoryOpenFileRead(dir, &kf, path);
    if (rc == 0) {
        KFile* fnull = NULL;
        rc = KFileMakeNullUpdate(&fnull);
        if (rc == 0) {
            KMD5SumFmt* fmt = NULL;
            rc = KMD5SumFmtMakeUpdate(&fmt, fnull);
            if (rc == 0) {
                const KFile* md5 = NULL;
                rc = KFileMakeNewMD5Read(&md5, kf, fmt, path);
                if (rc == 0) {
                    uint64_t ps = 0;
                    char buffer[512];
                    size_t read = 0;
                    do {
                        rc = KFileRead(md5, ps, buffer, sizeof buffer, &read);
                        if (rc == 0)
                        {   ps += read; }
                    } while (rc == 0 && read > 0);
                    if (rc == 0) {
                        bool bin;
                        rc = KMD5SumFmtFind(fmt, path, digest, &bin);
                    }
                }
                RELEASE(KFile, md5);
            }
            RELEASE(KMD5SumFmt, fmt);
        }
/*      RELEASE(KFile, fnull); fnull is released by KMD5SumFmt* fmt */
    }
/*  RELEASE(KFile, kf); kf is released by KFile* md5 */
    return rc;
}

static
rc_t ReportAlias(const ReportFuncs *f, uint32_t indent, const char* alias, const KDirectory* dir)
{
    char resolved[PATH_MAX + 1];
    rc_t rc
        = KDirectoryResolveAlias(dir, false, resolved, sizeof resolved, alias);
    if (rc == 0) {
        const char tag[] = "Alias";
        uint32_t type = KDirectoryPathType(dir, resolved);
        if (type & kptAlias) {
            reportOpen(indent, tag, 1, "resolved", 's', resolved);
            rc = ReportAlias(f, indent + 1, resolved, dir);
            reportClose(indent, tag);
        }
        else
        {   report(indent, tag, 1, "resolved", 's', resolved); }
    }
    return rc;
}

static rc_t ReportBinary(const ReportFuncs *f, uint32_t indent, const char* argv0) {
    rc_t rc = 0;
#ifdef _STATIC
    KDyld *dyld = NULL;
    assert(argv0);
    rc = KDyldMake(&dyld);
    if (rc != 0) {
        reportError(indent + 1, rc, "KDyldMake");
    }
    else {
        const KDirectory* dir = NULL;
        rc = KDyldHomeDirectory(dyld, &dir, (fptr_t) ReportFinalize);
        if (rc != 0) {
            reportError(indent + 1, rc, "KDyldHomeDirectory");
        }
        else {
            char binary[PATH_MAX + 1];
            const char* name = strpbrk(argv0, "/\\");
            const char* last_name = name;
            if (last_name)
            {   ++last_name; }
            while (name) {
                name = strpbrk(last_name, "/\\");
                if (name) {
                    last_name = name;
                    if (last_name)
                    {   ++last_name; }
                }
            }
            name = last_name ? last_name : argv0;
            rc = KDirectoryResolvePath(dir, true, binary, sizeof binary, name);
            if (rc != 0) {
                reportErrorStr(indent + 1, rc, "KDirectoryResolvePath",
                    "origin", "KDyldHomeDirectory");
            }
            else {
                bool found = false;
                const char tag[] = "Binary";
                const char* sType = NULL;
                uint8_t digest[16];
                uint32_t type = KDirectoryPathType(dir, binary);
                switch (type & ~kptAlias) {
                    case kptFile:
                        sType = type & kptAlias ? "alias" : "file";
                        found = true;
                        break;
                    case kptNotFound:
                        sType = "not found";
                        break;
                    default:
                        sType = "unknown";
                        break;
                }
                if (found)
                {   rc = md5(name, digest, dir); }
                if (type & kptAlias) {
                    if (found && rc == 0)  {
                        reportOpen(indent, tag, 3, "path", 's', binary,
                            "type", 's', sType, "md5", 'M', digest);
                    }
                    else {
                        reportOpen(indent, tag, 2, "path", 's', binary,
                            "type", 's', sType);
                    }
                    if (rc == 0 && type & kptAlias)
                    {   rc = ReportAlias(f, indent + 1, name, dir); }
                    reportClose(indent, tag);
                }
                else {
                    if (found && rc == 0)  {
                        report(indent, tag, 3, "path", 's', binary,
                            "type", 's', sType, "md5", 'M', digest);
                    }
                    else {
                        report(indent, tag, 2, "path", 's', binary,
                            "type", 's', sType);
                    }
                }
            }
        }
        RELEASE(KDirectory, dir);
    }
    RELEASE(KDyld, dyld);
#endif
    return rc;
}
#endif

static rc_t CC ReportSOFTWARE(const ReportFuncs *f, uint32_t indent, const char *argv_0, const char *date, ver_t tool_ver ) {
    rc_t rc = 0;

    Report* self = NULL;
    ReportGet(&self);
    assert(self);

    reportOpen(indent, "SOFTWARE", 0);

    if (self->mgr) {
        uint32_t version = 0;
        rc = VDBManagerVersion(self->mgr, &version);
        if (rc != 0) {
            reportOpen(indent + 1, "Library", 0);
            reportError(indent + 2, rc, "VDBManagerVersion");
            reportClose(indent + 1, "Library");
        }
        else { report(indent + 1, "VDBLibrary", 1, "vers", 'V', version); }
    }

    {
#ifdef _STATIC
        bool dynamic = false;
#else
        bool dynamic = true;
#endif

        rc_t rc2 = ReportBuild(f, indent + 1, dynamic, self->mgr);
        if (rc == 0 && rc2 != 0)
        {   rc = rc2; }
    }

    if (argv_0) {
        const char tag[] = "Tool";
#ifdef _STATIC
        reportOpen(indent + 1, tag, 3, "date", 's', date,
                   "name", 's', argv_0, "vers", 'V', tool_ver);
        {
            rc_t rc2 = ReportBinary(f, indent + 2, argv_0);
            if (rc == 0 && rc2 != 0)
            {   rc = rc2; }
        }
        reportClose(indent + 1, tag);
#else
        report(indent + 1, tag, 3, "date", 's', date,
                   "name", 's', argv_0, "vers", 'V', tool_ver);
#endif
    }

    reportClose(indent, "SOFTWARE");

    return rc;
}



/* SetVDBManager
 *  remember the manager in use
 */
LIB_EXPORT rc_t CC ReportSetVDBManager ( const VDBManager *mgr )
{
    rc_t rc = 0;

    Report* self = NULL;
    ReportGet(&self);
    if ( self != NULL )
    {
        rc = VDBManagerAddRef ( mgr );
        if ( rc == 0 )
        {
            rc = VDBManagerRelease ( self -> mgr );
            if ( rc != 0 )
                VDBManagerRelease ( mgr );
            else
                self -> mgr = mgr;
        }
    }
    return rc;
}


/* SetDatabase
 *  call it if you work with Database
 *
 *  "path" [ IN ] - path to the database that is used to access it
 */
LIB_EXPORT rc_t CC ReportResetDatabase ( const char *path, const VDatabase *db )
{
    rc_t rc = 0;

    Report* self = NULL;
    ReportGet(&self);
    if ( self != NULL )
    {
        VDatabaseRelease ( self -> db ), self -> db = NULL;
        VTableRelease ( self -> table ), self -> table = NULL;

        rc = VDatabaseAddRef ( db );
        if ( rc == 0 )
            self -> db = db;
    }

    return rc;
}


/* SetTable
 *  call it if you work with Table
 *
 *  "path" [ IN ] - path to the table that is used to access it
 */
LIB_EXPORT rc_t CC ReportResetTable ( const char *path, const VTable *tbl )
{
    rc_t rc = 0;

    Report* self = NULL;
    ReportGet(&self);
    if ( self != NULL )
    {
        VDatabaseRelease ( self -> db ), self -> db = NULL;
        VTableRelease ( self -> table ), self -> table = NULL;

        rc = VTableAddRef ( tbl );
        if ( rc == 0 )
            self -> table = tbl;
    }

    return rc;
}
