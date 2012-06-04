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
#ifndef _h_vfs_manager_
#define _h_vfs_manager_

#ifndef _h_vfs_extern_
#include <vfs/extern.h>
#endif

#ifndef _h_klib_defs_
#include <klib/defs.h>
#endif


#ifndef _h_kfs_defs_
#include <kfs/defs.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct KDirectory;
struct KFile;
struct VPath;


/*--------------------------------------------------------------------------
 * VFSManager
 */
typedef struct VFSManager VFSManager;

#define VFS_KRYPTO_PASSWORD_MAX_SIZE 4096


/* AddRef
 * Release
 *  ignores NULL references
 */
VFS_EXTERN rc_t CC VFSManagerAddRef ( const VFSManager *self );
VFS_EXTERN rc_t CC VFSManagerRelease ( const VFSManager *self );

/* OpenFileRead
 *  opens an existing file with read-only access
 *
 *  "f" [ OUT ] - return parameter for newly opened file
 *
 *  "path" [ IN ] VPath representing the path, URL or URN of the desired file
 */
VFS_EXTERN rc_t CC VFSManagerOpenFileRead (const VFSManager *self, 
                                           struct KFile const **f,
                                           const struct VPath * path);

VFS_EXTERN rc_t CC VFSManagerOpenDirectoryRead ( const VFSManager *self,
    struct KDirectory const **d, const struct VPath * path );



/* OpenFileWrite
 *  opens an existing file with write access
 *
 *  "f" [ OUT ] - return parameter for newly opened file
 *
 *  "update" [ IN ] - if true, open in read/write mode
 *  otherwise, open in write-only mode
 *
 *  "path" [ IN ] VPath representing the path, URL or URN of the desired file
 */
VFS_EXTERN rc_t CC VFSManagerOpenFileWrite (const VFSManager *self,
                                            struct KFile **f,
                                            bool update,
                                            const struct VPath * path);

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
VFS_EXTERN rc_t CC VFSManagerCreateFile (const VFSManager *self, 
                                         struct KFile **f,
                                         bool update, uint32_t access,
                                         KCreateMode mode,
                                         const struct VPath * path );


/* Remove
 *  remove an accessible object from its directory
 *
 *  "force" [ IN ] - if true and target is a directory,
 *  remove recursively
 *
 *  "path" [ IN ] - NUL terminated string in directory-native
 *  character set denoting target object
 */
VFS_EXTERN rc_t CC VFSManagerRemove (const VFSManager *self, bool force,
                                     const struct VPath * path );


/* Make
 */
VFS_EXTERN rc_t CC VFSManagerMake ( VFSManager ** pmanager );


/* GetConfigPWFile
 */
VFS_EXTERN rc_t CC VFSManagerGetConfigPWFile (const VFSManager * self, 
                                              char * b, size_t bz, size_t * pz);


/* GetCWD
 */
VFS_EXTERN rc_t CC VFSManagerGetCWD (const VFSManager * self, struct KDirectory ** cwd);


VFS_EXTERN rc_t CC VFSManagerGetKryptoPassword (const VFSManager * self, char * new_password, size_t max_size, size_t * size);

#if 0 /* obsoleted? */
VFS_EXTERN rc_t CC VFSManagerSetKryptoPassword (const VFSManager * self, char * new_password, size_t size);
VFS_EXTERN rc_t CC VFSManagerResetKryptoPassword (const VFSManager * self, 
                                                  char * new_password, size_t size,
                                                  char * new_password, size_t size);
#else

/*
  NULL value for self
  RC (rcVFS, rcEncryptionKey, rcUpdating, rcSelf, rcNull);

  NULL value for password or 0 value for size
  RC (rcVFS, rcEncryptionKey, rcUpdating, rcParam, rcNull);

  size greater than VFS_KRYPTO_PASSWORD_MAX_SIZE
  RC (rcVFS, rcEncryptionKey, rcUpdating, rcSize, rcExcessive);

  illegal CR or LF (NL) in the password
  RC (rcVFS, rcEncryptionKey, rcUpdating, rcEncryptionKey, rcInvalid);

  path/file name for password too long for function as written
  RC (rcVFS, rcEncryptionKey, rcUpdating, rcPath, rcExcessive);

  existing password path/file name is not a file
  RC (rcVFS, rcEncryptionKey, rcUpdating, rcPath, rcIncorrect);

  unknown file type for configured path/file name
  RC (rcVFS, rcEncryptionKey, rcUpdating, rcPath, rcCorrupt);

  incomplete writes to temporary password file
  RC (rcVFS, rcEncryptionKey, rcWriting, rcFile, rcInsufficient);




  other errors from KFS and KLIB
*/



VFS_EXTERN rc_t CC VFSManagerUpdateKryptoPassword (const VFSManager * self, 
                                                   const char * password,
                                                   size_t size);

#endif
#ifdef __cplusplus
}
#endif

#endif /* _h_kfs_manager_ */
