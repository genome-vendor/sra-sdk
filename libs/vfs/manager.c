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

#include <vfs/extern.h>
#include <vfs/manager.h>
#include <vfs/path.h>
#include "path-priv.h"

#include <krypto/key.h>
#include <krypto/encfile.h>
#include <krypto/wgaencrypt.h>

#include <kfg/config.h>

#include <kfs/directory.h>
#include <kfs/file.h>
#include <kfs/sra.h>
#include <kfs/tar.h>
#include <kfs/dyload.h>
#include <kfs/kfs-priv.h>
#include <kfs/nullfile.h>

#include <klib/refcount.h>
#include <klib/out.h>
#include <klib/log.h>
#include <klib/rc.h>

#include <sysalloc.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>


#ifdef _DEBUGGING
#define MGR_DEBUG(msg) DBGMSG(DBG_VFS,DBG_FLAG(DBG_VFS_MGR), msg)
#else
#define MGR_DEBUG(msg)
#endif


/*--------------------------------------------------------------------------
 * VFSManager
 */

/* currently expected to be a singleton and not use a vtable but
 * be fully fleashed out here */
struct VFSManager
{
    KRefcount refcount;

    /* the current directory in the eyes of the O/S when created */
    KDirectory * cwd;

    /* the underlying operating systems view of the path of the 
     * current working directory */
    VPath * cpath;


    /* configuration manager */
    KConfig * cfg;
};
static const char kfsmanager_classname [] = "VFSManager";

static 
VFSManager * singleton = NULL;


/* Destroy
 *  destroy file
 */
LIB_EXPORT rc_t CC VFSManagerDestroy ( VFSManager *self )
{
    if ( self == NULL )
        return RC ( rcFS, rcFile, rcDestroying, rcSelf, rcNull );

    KRefcountWhack (&self->refcount, kfsmanager_classname);

    KDirectoryRelease (self->cwd);

    VPathRelease (self->cpath);

    KConfigRelease (self->cfg);

    free (self);
    singleton = NULL;
    return 0;
}

/* AddRef
 *  creates a new reference
 *  ignores NULL references
 */
LIB_EXPORT rc_t CC VFSManagerAddRef ( const VFSManager *self )
{
    if (self != NULL)
    {
        switch (KRefcountAdd (&self->refcount, kfsmanager_classname))
        {
        case krefOkay:
            break;
        case krefZero:
            return RC (rcFS, rcMgr, rcAttaching, rcRefcount, rcIncorrect);
        case krefLimit:
            return RC (rcFS, rcMgr, rcAttaching, rcRefcount, rcExhausted);
        case krefNegative:
            return RC (rcFS, rcMgr, rcAttaching, rcRefcount, rcInvalid);
        default:
            return RC (rcFS, rcMgr, rcAttaching, rcRefcount, rcUnknown);
        }
    }
    return 0;
}

/* Release
 *  discard reference to file
 *  ignores NULL references
 */
LIB_EXPORT rc_t CC VFSManagerRelease ( const VFSManager *self )
{
    rc_t rc = 0;
    if (self != NULL)
    {
        switch (KRefcountDrop (&self->refcount, kfsmanager_classname))
        {
        case krefOkay:
        case krefZero:
            break;
        case krefWhack:
            rc = VFSManagerDestroy ((VFSManager*)self);
            break;
        case krefNegative:
            return RC (rcFS, rcMgr, rcAttaching, rcRefcount, rcInvalid);
        default:
            rc = RC (rcFS, rcMgr, rcAttaching, rcRefcount, rcUnknown);
            break;            
        }
    }
    return rc;
}



LIB_EXPORT rc_t CC VFSManagerGetConfigPWFd (const VFSManager * self, 
                                            char * b, size_t bz, size_t * pz)
{
    const KConfigNode * node;
    size_t oopsy;
    size_t z;
    rc_t rc;

    *pz = 0;

    rc = KConfigOpenNodeRead (self->cfg, &node, "krypto/pwfd");
    if (rc == 0)
    {
        rc = KConfigNodeRead (node, 0, b, bz-1, &z, &oopsy);
        if (rc == 0)
        {
            if (oopsy != 0)
                rc = RC (rcKrypto, rcMgr, rcReading, rcBuffer, rcInsufficient);
            else
            {
                b[z] = '\0';
                *pz = z;
            }
        }
        KConfigNodeRelease (node);
    }
    return rc;
}


LIB_EXPORT rc_t CC VFSManagerGetConfigPWFile (const VFSManager * self, 
                                              char * b, size_t bz, size_t * pz)
{
    const char * env;
    const KConfigNode * node;
    size_t oopsy;
    size_t z;
    rc_t rc;

    if (pz)
        *pz = 0;

    env = getenv ("VDB_PWFILE");

    if (env)
    {
        z = string_copy (b, bz, env, string_size (env));

        /* force a NUL that string_copy might have omitted 
         * even if this truncates the path */
        b[bz-1] = '\0';

        if (pz)
            *pz = z;
       
        return 0;
    }

    rc = KConfigOpenNodeRead (self->cfg, &node, "krypto/pwfile");
    if (rc == 0)
    {
        rc = KConfigNodeRead (node, 0, b, bz-1, &z, &oopsy);
        if (rc == 0)
        {
            if (oopsy != 0)
                rc = RC (rcKrypto, rcMgr, rcReading, rcBuffer, rcInsufficient);
            else
            {
                b[z] = '\0';
                *pz = z;
            }
        }
        KConfigNodeRelease (node);
    }
    return rc;
}


/* OpenFileRead
 *  opens an existing file with read-only access
 *
 *  "f" [ OUT ] - return parameter for newly opened file
 *
 *  "path" [ IN ] - NUL terminated string in directory-native
 *  character set denoting target file
 */
LIB_EXPORT rc_t CC VFSManagerOpenFileRead (const VFSManager *self,
                                           KFile const **f,
                                           const VPath * path)
{
    /* -----
     * this is a first pass that only opens files directory referenced from 
     * the ced or have a sysdir root; that is it uses KSysDir and KSysFile
     * only.
     */
    const KFile * file = NULL;
    size_t num_read;
    char pbuff [4096];
    rc_t rc;

    if ((f == NULL) || (path == NULL))
        return RC (rcFS, rcMgr, rcOpening, rcParam, rcNull);

    *f = NULL;

    if (self == NULL)
        return RC (rcFS, rcMgr, rcOpening, rcSelf, rcNull);

    rc = VPathReadPath (path, pbuff, sizeof pbuff, &num_read);
    if (rc == 0)
    {
        /* handle a few special case path names
         * This probably needs to be system specifica eventually
         */
        if (strncmp ("/dev/", pbuff, sizeof "/dev/" - 1) == 0)
        {

            if (strcmp ("/dev/stdin", pbuff) == 0)
                rc = KFileMakeStdIn (&file);
            else if (strcmp ("/dev/null", pbuff) == 0)
                rc = KFileMakeNullRead (&file);
            else if (strncmp ("/dev/fd/", pbuff, sizeof "/dev/fd/" - 1) == 0)
            {
                char * pc;
                size_t ix;

                pc = pbuff + sizeof "/dev/fd/" - 1;

                for (ix = 0; isdigit (pc[ix]); ++ix)
                    ;

                if ((ix > 0)&&(pc[ix] == '\0'))
                {
                    int fd = atoi (pc);

                    rc = KFileMakeFDFileRead (&file, fd);
                }
            }
        }
        if ((rc == 0)&&(file == NULL))
        {
            char rbuff [4096];

            rc = KDirectoryResolvePath (self->cwd, true, rbuff, sizeof rbuff, pbuff);
            if (rc == 0)
            {
                uint32_t type;

                type = KDirectoryPathType (self->cwd, rbuff);
                switch (type & ~kptAlias)
                {
                case kptNotFound:
                    rc = RC (rcFS, rcMgr, rcOpening, rcFile, rcNotFound);
                    break;

                case kptFile:
                    rc = KDirectoryOpenFileRead (self->cwd, &file, rbuff);
                    break;

                case kptBadPath:
                    rc = RC (rcFS, rcMgr, rcOpening, rcFile, rcInvalid);
                    break;

                case kptDir:
                case kptCharDev:
                    case kptBlockDev:
                case kptFIFO:
                case kptZombieFile:
                    rc = RC (rcFS, rcMgr, rcOpening, rcFile, rcIncorrect);
                    break;

                default:
                    rc = RC (rcFS, rcMgr, rcOpening, rcFile, rcUnknown);
                    break;
                }
            }
        }
    }
    if (rc == 0)
    {
        size_t z;
        char obuff [4096 + 2];

        if (VPathOption (path, vpopt_encrypted, obuff, sizeof obuff, &z) == 0)
        {
            const KFile * pwfile;

            /* -----
             * get the path to and open the pwfile
             *
             * first check the option for pwfile in the VPath
             * then check the option for pwfd
             * then check the environment
             * then check the configuration
             */

            if (VPathOption (path, vpopt_pwpath, obuff, sizeof obuff, &z) == 0)
                rc = KDirectoryOpenFileRead (self->cwd, &pwfile, obuff);

            else if (VPathOption (path, vpopt_pwfd, obuff, sizeof obuff, &z) == 0)
                rc = KFileMakeFDFileRead (&pwfile, atoi (obuff));

            else if ((rc = VFSManagerGetConfigPWFile(self, obuff, sizeof obuff, &z)) == 0)
                rc = KDirectoryOpenFileRead (self->cwd, &pwfile, obuff);

            else
                rc = RC (rcFS, rcPath, rcConstructing, rcParam, rcUnsupported);

            /* if have opened a password file */
            if (rc == 0)
            {
                /* read the password into obuff */
                rc = KFileRead (pwfile, 0, obuff, sizeof obuff - 1, &z);

                /* we're done with the password now */
                KFileRelease (pwfile);

                if (rc == 0)
                {
                    char * pc;
                    size_t tz;
                    char tbuff [VFS_KRYPTO_PASSWORD_MAX_SIZE];

                    /* -----
                     * trim back the contents of the file to
                     * a single ASCII/UTF-8 text line
                     * We actually only check for the two normal
                     * end of line characters so it could have other
                     * control characters...
                     */

                    obuff[z] = '\0';

                    pc = string_chr (obuff, z, '\r');
                    if (pc)
                    {
                        *pc = '\0';
                        z = pc - obuff;
                    }
                    pc = string_chr (obuff, z, '\n');
                    if (pc)
                    {
                        *pc = '\0';
                        z = pc - obuff;
                    }

                    if (z == 0)
                        rc = RC (rcVFS, rcMgr, rcOpening, rcEncryptionKey, rcTooShort);

                    else if  (VFS_KRYPTO_PASSWORD_MAX_SIZE < z) /* pwz came in as greater than 4096 */
                        rc = RC (rcVFS, rcMgr, rcOpening, rcEncryptionKey, rcTooLong);

                    /* -----
                     * pre-read 4kb from the 'encrypted file'
                     */
                    else
                        rc = KFileRead (file, 0, tbuff, sizeof tbuff, &tz);

                    if (rc == 0)
                    {
                        const KFile * encfile;

                        /* is this the header of an ecnrypted file? */
                        if (KFileIsEnc (tbuff, tz) == 0)
                        {
                            KKey key;

                            /* create the AES Key */
                            rc = KKeyInitRead (&key, kkeyAES128, obuff, z);
/*                             obuff[z] = '\0'; */
                            if (rc == 0)
                            {
                                rc = KEncFileMakeRead (&encfile, file, &key);
                                if (rc == 0)
                                {
                                    KFileRelease (file);
                                    *f = encfile;
                                    return 0;
                                }
                            }
                        }
                        else if (KFileIsWGAEnc (tbuff, tz) == 0)
                        {
                            rc = KFileMakeWGAEncRead (&encfile, file, obuff, z);
                            if (rc == 0)
                            {
                                KFileRelease (file);
                                *f = encfile;
                                return 0;
                            }
                        }
                        else
                            rc = RC (rcFS, rcPath, rcConstructing, rcFile, rcWrongType);
                    }
                }
            }
            if (rc)
                KFileRelease (file);
        }
        else
        {
            *f = file;
            return 0;
        }
    }
    return rc;
}
static
rc_t get_config_password (const VFSManager * self, char * b, size_t s, size_t * p)
{
    const KFile * pwfile;
    size_t z;
    rc_t rc;
    char obuff [4096 + 16];

    assert (b && s && p);

    rc = VFSManagerGetConfigPWFile(self, obuff, sizeof obuff, &z);
    if (rc == 0)
    {
        VPath * vpath;

        rc = VPathMake (&vpath, obuff);
        if (rc)
            ;

        else
        {
            rc =  VFSManagerOpenFileRead (self, &pwfile, vpath);
            if (rc)
                ;

            else
            {
                VPathRelease (vpath);

                /* read the password into obuff */
                rc = KFileReadAll (pwfile, 0, obuff, sizeof obuff - 1, &z);
            
                KFileRelease (pwfile);

                if (rc == 0)
                {
                    char * pc;

                    /* -----
                     * trim back the contents of the file to
                     * a single ASCII/UTF-8 text line
                     * We actually only check for the two normal
                     * end of line characters so it could have other
                     * control characters...
                     */
                    obuff[z] = '\0';
                    pc = string_chr (obuff, z, '\r');
                    if (pc)
                    {
                        *pc = '\0';
                        z = pc - obuff;
                    }
                    pc = string_chr (obuff, z, '\n');
                    if (pc)
                    {
                        *pc = '\0';
                        z = pc - obuff;
                    }

                    if (z > 4096) /* arbitrary maximum size */
                        rc = RC (rcKrypto, rcEncryptionKey, rcRetrieving, rcEncryptionKey, rcExcessive);

                    else if (s < z)
                        rc = RC (rcKrypto, rcEncryptionKey, rcRetrieving, rcBuffer, rcInsufficient);

                    else
                    {
                        memcpy (b, obuff, z);
                        if (s>z)
                            b[z] = '\0';
                        *p = z;
                    }
                }
            }
        }
    }
    return rc;
}
LIB_EXPORT rc_t CC VFSManagerOpenDirectoryRead (const VFSManager *self,
                                                KDirectory const **d,
                                                const VPath * path)
{
    const KFile * file = NULL;
    size_t num_read;
    char pbuff [4096]; /* path buffer */
    rc_t rc;

    if ((d == NULL) || (path == NULL))
        return RC (rcFS, rcMgr, rcOpening, rcParam, rcNull);

    *d = NULL;

    if (self == NULL)
        return RC (rcFS, rcMgr, rcOpening, rcSelf, rcNull);

    rc = VPathReadPath (path, pbuff, sizeof pbuff, &num_read);
    if (rc == 0)
    {
        char rbuff [4096]; /* resolved path buffer */

        rc = KDirectoryResolvePath (self->cwd, true, rbuff, sizeof rbuff, pbuff);
        if (rc == 0)
        {
            uint32_t type;

            type = KDirectoryPathType (self->cwd, rbuff);
            switch (type & ~kptAlias)
            {
            case kptNotFound:
                rc = RC (rcFS, rcMgr, rcOpening, rcFile, rcNotFound);
                break;

            case kptFile:
                rc = VFSManagerOpenFileRead (self, &file, path);
                if (rc == 0)
                {
                    rc = KFileRandomAccess(file);
                    if (rc == 0)
                    {
                        size_t tz;
                        char tbuff [4096];

                        rc = KFileRead (file, 0, tbuff, sizeof tbuff, &tz);
                        if (rc == 0)
                        {
                            const KFile * f;
                            size_t kz;
                            char kbuff [4096];

                            /* handle possible encryption */
                            if (KFileIsEnc (tbuff, tz) == 0)
                            {
                                rc = get_config_password (self, kbuff, sizeof kbuff, &kz);
                                if (rc == 0)
                                {
                                    KKey key;

                                    rc = KKeyInitRead (&key, kkeyAES128, kbuff, kz);
                                    if (rc == 0)
                                    {
                                        rc = KEncFileMakeRead (&f, file, &key);
                                        if (rc == 0)
                                        {
                                            KFileRelease (file);
                                            file = f;
                                        }
                                    }
                                }
                            }
                            else if (KFileIsWGAEnc (tbuff, tz) == 0)
                            {
                                rc = get_config_password (self, kbuff, sizeof kbuff, &kz);
                                if (rc == 0)
                                {
                                    rc = KFileMakeWGAEncRead (&f, file, kbuff, tz);
                                    if (rc == 0)
                                    {
                                        KFileRelease (file);
                                        file = f;
                                    }
                                }
                            }

                            /* now make an archive directory */
                            if (KFileIsSRA (tbuff, tz) == 0)
                            {
                                rc = KDirectoryOpenSraArchiveReadUnbounded_silent_preopened
                                    (self->cwd, d, false, file, rbuff);
                            }
                            else
                            {
                                rc = KDirectoryOpenTarArchiveRead_silent_preopened
                                    (self->cwd, d, false, file, rbuff);
                            }
                        }
                    }
                }
                break;

            case kptBadPath:
                rc = RC (rcFS, rcMgr, rcOpening, rcFile, rcInvalid);
                break;

            case kptDir:
                rc = KDirectoryOpenDirRead (self->cwd, d, false, rbuff);
                return rc;

            case kptCharDev:
            case kptBlockDev:
            case kptFIFO:
            case kptZombieFile:
                rc = RC (rcFS, rcMgr, rcOpening, rcFile, rcIncorrect);
                break;

            default:
                rc = RC (rcFS, rcMgr, rcOpening, rcFile, rcUnknown);
                break;
            }
        }
    }
    return rc;
}


/* OpenFileWrite
 *  opens an existing file with write access
 *
 *  "f" [ OUT ] - return parameter for newly opened file
 *
 *  "update" [ IN ] - if true, open in read/write mode
 *  otherwise, open in write-only mode
 *
 *  "path" [ IN ] - NUL terminated string in directory-native
 *  character set denoting target file
 */
LIB_EXPORT rc_t CC VFSManagerOpenFileWrite (const VFSManager *self,
                                            KFile **f, bool update,
                                            const VPath * path )
{
    /* -----
     * this is a first pass that only opens files directory referenced from 
     * the ced or have a sysdir root; that is it uses KSysDir and KSysFile
     * only.
     */
    KFile * file = NULL;
    size_t num_read;
    char pbuff [4096];
    rc_t rc;

    if ((f == NULL) || (path == NULL))
        return RC (rcFS, rcMgr, rcOpening, rcParam, rcNull);

    *f = NULL;

    if (self == NULL)
        return RC (rcFS, rcMgr, rcOpening, rcSelf, rcNull);

    rc = VPathReadPath (path, pbuff, sizeof pbuff, &num_read);
    if (rc == 0)
    {
        /* handle a few special case path names
         * This probably needs to be system specifica eventually
         */
        if (strncmp ("/dev/", pbuff, sizeof "/dev/" - 1) == 0)
        {

            if (strcmp ("/dev/stdout", pbuff) == 0)
                rc = KFileMakeStdOut (&file);
            else if (strcmp ("/dev/stderr", pbuff) == 0)
                rc = KFileMakeStdErr (&file);
            else if (strcmp ("/dev/null", pbuff) == 0)
                rc = KFileMakeNullUpdate (&file);
            else if (strncmp ("/dev/fd/", pbuff, sizeof "/dev/fd/" - 1) == 0)
            {
                char * pc;
                size_t ix;

                pc = pbuff + sizeof "/dev/fd/" - 1;

                for (ix = 0; isdigit (pc[ix]); ++ix)
                    ;

                if ((ix > 0)&&(pc[ix] == '\0'))
                {
                    int fd = atoi (pc);

                    rc = KFileMakeFDFileWrite (&file, update, fd);
                }
            }
        }
        if ((rc == 0)&&(file == NULL))
        {
            char rbuff [4096];

            rc = KDirectoryResolvePath (self->cwd, true, rbuff, sizeof rbuff, pbuff);
            if (rc == 0)
            {
                uint32_t type;

                type = KDirectoryPathType (self->cwd, rbuff);
                switch (type & ~kptAlias)
                {
                case kptNotFound:
                    rc = RC (rcFS, rcMgr, rcOpening, rcFile, rcNotFound);
                    break;

                case kptFile:
                    rc = KDirectoryOpenFileWrite (self->cwd, f, update, rbuff);
                    break;

                case kptBadPath:
                    rc = RC (rcFS, rcMgr, rcOpening, rcFile, rcInvalid);
                    break;
                case kptDir:
                case kptCharDev:
                case kptBlockDev:
                case kptFIFO:
                case kptZombieFile:
                    rc = RC (rcFS, rcMgr, rcOpening, rcFile, rcIncorrect);
                    break;

                default:
                    rc = RC (rcFS, rcMgr, rcOpening, rcFile, rcUnknown);
                    break;
                }
            }
        }
    }
    if (rc == 0)
    {
        size_t z;
        char obuff [4096+1];

        if (VPathOption (path, vpopt_encrypted, obuff, sizeof obuff, &z) == 0)
        {
            const KFile * pwfile;

            if (VPathOption (path, vpopt_pwpath, obuff, sizeof obuff, &z) == 0)
                rc = KDirectoryOpenFileRead (self->cwd, &pwfile, obuff);

            else if (VPathOption (path, vpopt_pwfd, obuff, sizeof obuff, &z) == 0)
                rc = KFileMakeFDFileRead (&pwfile, atoi (obuff));

            else if ((rc = VFSManagerGetConfigPWFile(self, obuff, sizeof obuff, &z)) == 0)
                rc = KDirectoryOpenFileRead (self->cwd, &pwfile, obuff);

            else
                rc = RC (rcFS, rcPath, rcConstructing, rcParam, rcUnsupported);

            if (rc == 0)
            {
                rc = KFileRead (pwfile, 0, obuff, sizeof obuff - 1, &z);
                
                KFileRelease (pwfile);

                if (rc == 0)
                {
                    KKey key;
                    KFile * encfile;
                    char * pc;

                    obuff[z] = '\0';
                    pc = string_chr (obuff, z, '\r');
                    if (pc)
                    {
                        *pc = '\0';
                        z = pc - obuff;
                    }
                    pc = string_chr (obuff, z, '\n');
                    if (pc)
                    {
                        *pc = '\0';
                        z = pc - obuff;
                    }

                    rc = KKeyInitUpdate (&key, kkeyAES128, obuff, z);
                    if (rc == 0)
                    {
                        rc = KEncFileMakeWrite (&encfile, file, &key);
                        if (rc == 0)
                        {
                            *f = encfile;
                            return 0;
                        }
                    }
                }
            }
            if (rc)
                KFileRelease (file);
        }
        else
        {
            *f = file;
            return 0;
        }
    }
    return rc;
}


/* CreateFile
 *  opens a file with write access
 *
 *  "f" [ OUT ] - return parameter for newly opened file
 *
 *  "update" [ IN ] - if true, open in read/write mode
 *  otherwise, open in write-only mode
 *
 *  "access" [ IN ] - standard Unix access mode, e.g. 0664
 *
 *  "mode" [ IN ] - a creation mode ( see explanation above ).
 *
 *  "path" [ IN ] VPath representing the path, URL or URN of the desired file
 */
LIB_EXPORT rc_t CC VFSManagerCreateFile ( const VFSManager *self, KFile **f,
    bool update, uint32_t access, KCreateMode mode, const VPath * path )
{
    /* -----
     * this is a first pass that only opens files directory referenced from 
     * the ced or have a sysdir root; that is it uses KSysDir and KSysFile
     * only.
     */
    KFile * file = NULL;
    size_t num_read;
    rc_t rc;
    bool file_created = false;
    char pbuff [4096];
    char rbuff [4096];

    if ((f == NULL) || (path == NULL))
        return RC (rcFS, rcMgr, rcOpening, rcParam, rcNull);

    *f = NULL;

    if (self == NULL)
        return RC (rcFS, rcMgr, rcOpening, rcSelf, rcNull);

    rc = VPathReadPath (path, pbuff, sizeof pbuff, &num_read);
    if (rc == 0)
    {

        /* handle a few special case path names
         * This probably needs to be system specifica eventually
         */
        if (strncmp ("/dev/", pbuff, sizeof "/dev/" - 1) == 0)
        {

            if (strcmp ("/dev/stdout", pbuff) == 0)
                rc = KFileMakeStdOut (&file);
            else if (strcmp ("/dev/stderr", pbuff) == 0)
                rc = KFileMakeStdErr (&file);
            else if (strcmp ("/dev/null", pbuff) == 0)
                rc = KFileMakeNullUpdate (&file);
            else if (strncmp ("/dev/fd/", pbuff, sizeof "/dev/fd/" - 1) == 0)
            {
                char * pc;
                size_t ix;

                pc = pbuff + sizeof "/dev/fd/" - 1;

                for (ix = 0; isdigit (pc[ix]); ++ix)
                    ;

                if ((ix > 0)&&(pc[ix] == '\0'))
                {
                    int fd = atoi (pc);

                    rc = KFileMakeFDFileWrite (&file, update, fd);
                }
            }
        }
        if ((rc == 0)&&(file == NULL))
        {
            rc = KDirectoryResolvePath (self->cwd, true, rbuff, sizeof rbuff, pbuff);
            if (rc == 0)
            {
                uint32_t type;

                type = KDirectoryPathType (self->cwd, rbuff);
                switch (type & ~kptAlias)
                {
                case kptNotFound:
                case kptFile:
                    rc = KDirectoryCreateFile (self->cwd, &file, update, access, mode,
                                               rbuff);
                    if (rc == 0)
                        file_created = true;
                    break;

                case kptBadPath:
                    rc = RC (rcFS, rcMgr, rcOpening, rcFile, rcInvalid);
                    break;
                case kptDir:
                case kptCharDev:
                case kptBlockDev:
                case kptFIFO:
                case kptZombieFile:
                    rc = RC (rcFS, rcMgr, rcOpening, rcFile, rcIncorrect);
                    break;

                default:
                    rc = RC (rcFS, rcMgr, rcOpening, rcFile, rcUnknown);
                    break;
                }
            }
        }
    }
    if (rc == 0)
    {
        size_t z;
        char obuff [4096];

        if (VPathOption (path, vpopt_encrypted, obuff, sizeof obuff, &z) == 0)
        {
            const KFile * pwfile;

            if (VPathOption (path, vpopt_pwpath, obuff, sizeof obuff, &z) == 0)
                rc = KDirectoryOpenFileRead (self->cwd, &pwfile, obuff);

            else if (VPathOption (path, vpopt_pwfd, obuff, sizeof obuff, &z) == 0)
                rc = KFileMakeFDFileRead (&pwfile, atoi (obuff));
            else
                rc = RC (rcFS, rcPath, rcConstructing, rcParam, rcUnsupported);

            if (rc == 0)
            {
                rc = KFileRead (pwfile, 0, obuff, sizeof obuff, &z);
                
                KFileRelease (pwfile);

                if (rc == 0)
                {
                    KKey key;
                    KFile * encfile;
                    char * pc;

                    pc = string_chr (obuff, z, '\r');
                    if (pc)
                    {
                        *pc = '\0';
                        z = pc - obuff;
                    }
                    pc = string_chr (obuff, z, '\n');
                    if (pc)
                    {
                        *pc = '\0';
                        z = pc - obuff;
                    }

                    rc = KKeyInitUpdate (&key, kkeyAES128, obuff, z);

                    obuff[z] = '\0';

                    rc = KEncFileMakeWrite (&encfile, file, &key);
                    if (rc == 0)
                    {
                        *f = encfile;
                        return 0;
                    }
                }
            }
            if (rc)
                KFileRelease (file);
        }
        else
        {
            *f = file;
            return 0;
        }
    }
    if (rc && file_created)
        KDirectoryRemove (self->cwd, true, rbuff);
    return rc;
}


/* Remove
 *  remove an accessible object from its directory
 *
 *  "force" [ IN ] - if true and target is a directory,
 *  remove recursively
 *
 *  "path" [ IN ] - NUL terminated string in directory-native
 *  character set denoting target object
 */
LIB_EXPORT rc_t CC VFSManagerRemove ( const VFSManager *self, bool force,
                                      const VPath * path )
{
    /* -----
     * this is a first pass that only opens files directory referenced from 
     * the ced or have a sysdir root; that is it uses KSysDir and KSysFile
     * only.
     */
    size_t num_read;
    char pbuff [4096];
    rc_t rc;

    if (path == NULL)
        return RC (rcFS, rcMgr, rcOpening, rcParam, rcNull);

    if (self == NULL)
        return RC (rcFS, rcMgr, rcOpening, rcSelf, rcNull);

    rc = VPathReadPath (path, pbuff, sizeof pbuff, &num_read);
    if (rc == 0)
    {
        char rbuff [4096];
    
        rc = KDirectoryResolvePath (self->cwd, true, rbuff, sizeof rbuff, pbuff);
        if (rc == 0)
        {
            uint32_t type;

            type = KDirectoryPathType (self->cwd, rbuff);
            switch (type & ~kptAlias)
            {
            case kptNotFound:
                break;

            case kptFile:
            case kptDir:
            case kptCharDev:
            case kptBlockDev:
            case kptFIFO:
            case kptZombieFile:
                rc = KDirectoryRemove (self->cwd, force, rbuff);
                break;

            case kptBadPath:
                rc = RC (rcFS, rcMgr, rcOpening, rcFile, rcInvalid);
                break;
/*                 rc = RC (rcFS, rcMgr, rcOpening, rcFile, rcIncorrect); */
/*                 break; */

            default:
                rc = RC (rcFS, rcMgr, rcOpening, rcFile, rcUnknown);
                break;
            }
        }
    }
    return rc;
}


/* Make
 */
LIB_EXPORT rc_t CC VFSManagerMake ( VFSManager ** pmanager )
{
    if (pmanager == NULL)
        return RC (rcFS, rcMgr, rcConstructing, rcParam, rcNull);

    if (singleton)
    {
        *pmanager = singleton;
    }
    else
    {
        VFSManager * obj;
        rc_t rc;

        obj = calloc (1, sizeof (*obj));
        if (obj == NULL)
            return RC (rcFS, rcMgr, rcConstructing, rcMemory, rcExhausted);

        KRefcountInit (&obj->refcount, 1, kfsmanager_classname, "init", 
                       kfsmanager_classname);

        rc = KDirectoryNativeDir (&obj->cwd);
        if (rc)
        {
            obj->cwd = NULL;
            VFSManagerDestroy (obj);
            return rc;
        }

        rc = VPathMakeCurrentPath (&obj->cpath);
        if (rc)
        {
            obj->cpath = NULL;
            VFSManagerDestroy (obj);
            return rc;
        }

        rc = KConfigMake (&obj->cfg, NULL);
        if (rc)
        {
            obj->cfg = NULL;
            VFSManagerDestroy (obj);
            return rc;
        }

        *pmanager = singleton = obj;

    }
    return 0;
}


LIB_EXPORT rc_t CC VFSManagerGetCWD (const VFSManager * self, KDirectory ** cwd)
{
    rc_t rc;

    if (cwd == NULL)
        return RC (rcFS, rcMgr, rcAccessing, rcParam, rcNull);

    *cwd = NULL;

    if (self == NULL)
        return RC (rcFS, rcMgr, rcAccessing, rcSelf, rcNull);

    rc = KDirectoryAddRef (self->cwd);
    if (rc)
        return rc;

    *cwd = self->cwd;

    return 0;
}

LIB_EXPORT rc_t CC VFSManagerGetKryptoPassword (const VFSManager * self,
                                                char * password,
                                                size_t max_size,
                                                size_t * size)
{
    rc_t rc;

    if (self == NULL)
        rc = RC (rcVFS, rcMgr, rcAccessing, rcSelf, rcNull);

    else if ((password == NULL) || (max_size == 0) || (size == NULL))
        rc = RC (rcVFS, rcMgr, rcAccessing, rcParam, rcNull);

    else
        rc = get_config_password (self, password, max_size, size);

    return rc;
}


#if 0 /* not using this now */
LIB_EXPORT rc_t CC VFSManagerSetKryptoPassword (const VFSManager * self,
                                                char * password,
                                                size_t size)
{
    rc_t rc;

    if (self == NULL)
        rc = RC (rcVFS, rcMgr, rcAccessing, rcSelf, rcNull);

    else if ((password == NULL) || (size == NULL))
        rc = RC (rcVFS, rcMgr, rcAccessing, rcParam, rcNull);

    else
    {
        KFile * pwfile;
        size_t size;
        char buffer [8192];
        

        rc = VFSManagerGetConfigPWFile (self, buffer, sizeof buffer, &size);
        if (rc)
            PLOGERR (klogErr,
                     (klogErr,
                      "failed to obtain configured path for password file "
                      "'$(P)'", "P=%s", buffer));

        else
        {
            VPath * vpath;

            buffer[size] = '\0';

            rc = VPathMake (&vpath, buffer);
            if (rc)
                PLOGERR (klogErr,
                         (klogErr,
                          "failed to create vpath for password file '$(P)'",
                          "P=%s", buffer));

            else
            {
                KPath * file;

                rc = VFSManagerCreateFile (self, &file, false, 0600,
                                           kcmInit|kcmParents, vpath);
                if (rc)
                    PLOGERR (klogErr,
                             (klogErr,
                              "failed to create password file '$(P)'",
                              "P=%s", buffer));

                else
                {
                    size_t writ;

                    rc = KFileWriteAll (file, 0, password, size, &writ);
                    if (rc)
                        PLOGERR (klogErr,
                                 (klogErr,
                                  "failed to write to password file '$(P)'",
                                  "P=%s", buffer));
                    else
                    {
                        if (writ != size)
                        {
                            rc = KFileSetSize (file. 0);
                            if (rc)
                                PLOGERR (klogErr,
                                         (klogErr,
                                          "failed to emplty corrupt "
                                          "password file '$(P)'",
                                          "P=%s", buffer));

                            rc = RC (rcVfs, rcMgr, rcWriting,
                                     rcInsufficient);
                        }
                    }
                    KFileRelease (file);

                    if (rc)
                    {
                        rc_t orc;

                        orc = KDirectoryRemove (self->cwd, force, buffer);
                        if (orc)
                            PLOGERR (klogErr,
                                     (klogErr,
                                      "failed to delete bad password file "
                                      "'$(P)'", "P=%s", buffer));
                    }
                }
                VPathRelease (vpath);
            }
        }
    }

    return rc;
}
#endif

LIB_EXPORT rc_t CC VFSManagerUpdateKryptoPassword (const VFSManager * self, 
                                                   const char * password,
                                                   size_t size)
{
    static const char temp_extension [] = ".tmp";
    rc_t rc;

    if (self == NULL)
        rc = RC (rcVFS, rcEncryptionKey, rcUpdating, rcSelf, rcNull);

    else if ((password == NULL) || (size == 0))
        rc = RC (rcVFS, rcEncryptionKey, rcUpdating, rcParam, rcNull);

    else if (size > VFS_KRYPTO_PASSWORD_MAX_SIZE)
        rc = RC (rcVFS, rcEncryptionKey, rcUpdating, rcSize, rcExcessive);

    else if ((string_chr (password, size, '\n') != NULL) ||
             (string_chr (password, size, '\r') != NULL))
        rc = RC (rcVFS, rcEncryptionKey, rcUpdating, rcEncryptionKey, rcInvalid);

    else
    {
        size_t old_password_file_size;
        char old_password_file [8193];
        
        rc = VFSManagerGetConfigPWFile (self, old_password_file,
                                        sizeof old_password_file - 1,
                                        &old_password_file_size);
        if (rc)
            LOGERR (klogErr, rc, "failed to obtain configured path for password file");

        else if (old_password_file_size >= (sizeof old_password_file - 1))
        {
            rc = RC (rcVFS, rcEncryptionKey, rcUpdating, rcPath, rcExcessive);
            PLOGERR (klogErr,
                     (klogErr, rc, "configured path too long for function "
                      "'$(P)' '${F}'", "P=%s,F=%s",
                      old_password_file, __func__));
        }
        else
        {
            KPathType ftype;
            bool old_exists;

            old_password_file[old_password_file_size] = '\0';
            ftype = KDirectoryPathType (self->cwd, old_password_file);

            switch (ftype)
            {
            case kptNotFound:
                old_exists = false;
                break;

            case kptBadPath:
                rc = RC (rcVFS, rcEncryptionKey, rcUpdating, rcPath, rcInvalid);
                break;

            case kptFile:
                old_exists = true;
                break;

            case kptDir:
            case kptCharDev:
            case kptBlockDev:
            case kptFIFO:
            case kptZombieFile:
            case kptDataset:
            case kptDatatype:
                rc = RC (rcVFS, rcEncryptionKey, rcUpdating, rcPath, rcIncorrect);
                break;

            default:
                rc = RC (rcVFS, rcEncryptionKey, rcUpdating, rcPath, rcCorrupt);
                break;
            }

            if (rc)
                PLOGERR (klogErr,
                         (klogErr, rc, "can not use configured path for "
                          "password file '$(P)'", "P=%s", old_password_file));

            else
            {
                VPath * vold;
                size_t new_password_file_size;
                char new_password_file [sizeof old_password_file + sizeof temp_extension];
                size_t password_dir_size;
                char password_dir [sizeof old_password_file];
                bool save_old_password;
                char * pc;

                memcpy (password_dir, old_password_file, old_password_file_size);
                memcpy (new_password_file, old_password_file, old_password_file_size);
                memcpy (new_password_file + old_password_file_size, temp_extension, sizeof temp_extension);
                new_password_file_size = old_password_file_size + sizeof temp_extension - 1;

                pc = string_rchr (password_dir, old_password_file_size, '/');
                if (pc == NULL)
                {
                    password_dir[0] = '.';
                    pc = password_dir+1;
                }
                *pc = '\0';
                password_dir_size = pc - password_dir;

                rc = VPathMake (&vold, old_password_file);
                if (rc)
                    PLOGERR (klogErr,
                             (klogErr, rc, "could not create vpath for "
                              "password file '$(P)'", "P=%s",
                              old_password_file));

                else
                {
                    VPath * vnew;

                    rc = VPathMake (&vnew, new_password_file);
                    if (rc)
                        PLOGERR (klogErr,
                                 (klogErr, rc, "counld not create vpath for "
                                  "password file '$(P)'", "P=%s",
                                  new_password_file));

                    else
                    {
                        const KFile * fold = NULL;
                        KFile * fnew = NULL;

                        if (old_exists)
                        {
                            rc = VFSManagerOpenFileRead (self, &fold, vold);

                            if (rc)
                                PLOGERR (klogErr,
                                         (klogErr, rc, "unable to open existing "
                                          "password file '$(P)'", "P=%s",
                                          old_password_file));
                        }
                        

                        if (rc == 0)
                        {
                            rc = VFSManagerCreateFile (self, &fnew, false, 0600,
                                                       kcmInit|kcmParents,
                                                       vnew);
                            if (rc)
                                PLOGERR (klogErr,
                                         (klogErr, rc, "unable to open temporary "
                                          "password file '$(P)'", "P=%s",
                                          new_password_file));

                            else
                            {
                                uint64_t writ;
                                size_t this_writ;

                                rc = KFileWriteAll (fnew, 0, password, size, &this_writ);
                                if (rc)
                                    PLOGERR (klogErr,
                                             (klogErr, rc, "unable to write "
                                              "password to temporary password "
                                              "file '$(P)'", "P=%s",
                                              new_password_file));

                                else if (this_writ != size)
                                {
                                    rc = RC (rcVFS, rcEncryptionKey, rcWriting,
                                             rcFile, rcInsufficient);
                                    PLOGERR (klogErr,
                                             (klogErr, rc, "unable to write complete "
                                              "password to temporary password "
                                              "file '$(P)'", "P=%s",
                                              new_password_file));
                                }

                                else
                                {
                                    writ = this_writ;

                                    rc = KFileWriteAll (fnew, this_writ, "\n", 1, &this_writ);
                                    if (rc)
                                        PLOGERR (klogErr,
                                                 (klogErr, rc, "unable to write "
                                                  "password to temporary password "
                                                  "file '$(P)'", "P=%s",
                                                  new_password_file));

                                    else if (this_writ != 1)
                                    {
                                        rc = RC (rcVFS, rcEncryptionKey, rcWriting,
                                                 rcFile, rcInsufficient);
                                        PLOGERR (klogErr,
                                                 (klogErr, rc, "unable to write complete "
                                                  "password to temporary password "
                                                  "file '$(P)'", "P=%s",
                                                  new_password_file));
                                    }

                                    else
                                    {
                                        bool do_rename;

                                        do_rename = true;
                                        ++writ;

                                        if (old_exists)
                                        {
                                            uint64_t read;
                                            size_t this_read;
                                            char buffer [VFS_KRYPTO_PASSWORD_MAX_SIZE+4];

                                            rc = KFileReadAll (fold, 0, buffer,
                                                               sizeof buffer, &this_read);
                                            if (rc)
                                                ;

                                            else
                                            {
                                                read = this_read;
                                                /* look for duplicated password */
                                                if (read > size)
                                                {
                                                    char cc;

                                                    cc = buffer[size];
                                                    if (((cc == '\n') || (cc == '\r')) &&
                                                        (memcmp (buffer, password, size) == 0))
                                                    {
                                                        do_rename = false;
                                                    }
                                                }
                                                if (read)
                                                    rc = KFileWriteAll (fnew, writ, buffer, read, &this_writ);

                                                if (rc)
                                                    ;
                                                else if (do_rename)
                                                {
                                                    writ += this_writ;

                                                    do
                                                    {
                                                        rc = KFileReadAll (fold, read, buffer,
                                                                   sizeof buffer, &this_read);
                                                        if (rc)
                                                            ;

                                                        else if (this_read == 0)
                                                            break;

                                                        else
                                                        {
                                                            rc = KFileWriteAll (fnew, writ, buffer,
                                                                                this_read, &this_writ);
                                                            if (rc)
                                                                ;

                                                            else if (this_read != this_writ)
                                                            {
                                                                rc = RC (rcVFS, rcEncryptionKey, rcWriting,
                                                                         rcFile, rcInsufficient);
                                                                PLOGERR (klogErr,
                                                                         (klogErr, rc,
                                                                          "unable to write complete "
                                                                          "password to temporary password "
                                                                          "file '$(P)'", "P=%s",
                                                                          new_password_file));
                                                            }

                                                            else
                                                            {
                                                                read += this_read;
                                                                writ += this_writ;
                                                            }
                                                        }
                                                    } while (rc == 0);
                                                }
                                            }
                                            KFileRelease (fold);
                                            fold = NULL;
                                        }

                                        KFileRelease (fnew);
                                        fnew = NULL;

                                        if (rc == 0)
                                        {
                                            if (do_rename)
                                            {
                                                rc = KDirectoryRename (self->cwd, true, 
                                                                       new_password_file,
                                                                       old_password_file);
                                            }
                                            else
                                            {
                                                KDirectoryRemove (self->cwd, true, new_password_file);
                                            }

#if !WINDOWS
                                            if (rc == 0)
                                            {
                                                uint32_t access;

                                                rc = KDirectoryAccess (self->cwd,
                                                                       &access, password_dir);
                                                if (rc)
                                                    ;

                                                else
                                                {
                                                    if (access & 0027)
                                                        rc = RC (rcVFS, rcEncryptionKey, rcUpdating, rcDirectory, rcExcessive);
                                                }
                                            }
#endif
                                        }
                                    }
                                }
                                KFileRelease (fnew);
                            }
                            KFileRelease (fold);
                        }
                        VPathRelease (vold);
                    }
                    VPathRelease (vnew);
                }
            }
        }
    }
    return rc;
}


LIB_EXPORT rc_t CC VFSManagerGetCPath (const VFSManager * self, VPath ** cpath)
{
    rc_t rc;

    if (cpath == NULL)
        return RC (rcFS, rcMgr, rcAccessing, rcParam, rcNull);

    *cpath = NULL;

    if (self == NULL)
        return RC (rcFS, rcMgr, rcAccessing, rcSelf, rcNull);

    rc = VPathAddRef (self->cpath);
    if (rc)
        return rc;

    *cpath = self->cpath;

    return 0;
}
