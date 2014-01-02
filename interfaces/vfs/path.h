/*===========================================================================
*
*                            Public Domain Notice
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

#ifndef _h_vfs_path_
#define _h_vfs_path_

#ifndef _h_vfs_extern_
#include <vfs/extern.h>
#endif

#ifndef _h_klib_defs_
#include <klib/defs.h>
#endif

#include <stdarg.h>

#ifdef __cplusplus
extern "C" {
#endif


/*--------------------------------------------------------------------------
 * forwards
 */
struct String;
struct VFSManager;


/*--------------------------------------------------------------------------
 * VPath
 *  represents an unbound object query key
 *  may be created from a simple file-system path,
 *  a more formal URN or URL,
 *  or other modes of creation
 *
 *  a path will have these parts:
 *    scheme       : a scheme for retrieval
 *    auth         : login name for authentication
 *    host         : authoritative source
 *    port         : port for connecting with host
 *    path         : host-relative path
 *    query        : parameters for interpretation
 *    fragment     : internal component of object
 *    proj         : project id
 *    name         : alternate or primary name
 *
 *  file-system paths with no modifying parameters
 *  will be given standard "file" scheme. paths having
 *  parameters will be given the scheme "ncbi-file".
 *
 *  standard networking schemes ( "http", "ftp", etc. )
 *  are supported.
 *
 *  NCBI accessions are given the scheme "ncbi-acc".
 *
 *  NCBI remote object id paths receive scheme "ncbi-obj".
 */
typedef struct VPath VPath;


/* MakePath
 *  make a path object from a string conforming to
 *  either a standard POSIX path or a URI
 *
 *  "new_path" [ OUT ] - return parameter for new path object
 *
 *  "path_str" [ IN ] - a UTF-8 NUL-terminated string
 *  representing a POSIX path or URI, or
 *  a string_printf compatible format string
 *
 *  "path_fmt" [ IN ] and "args" [ IN ] - a UTF-8 NUL-terminated fmt string
 *  compatible with string_vprintf, plus argument list
 *
 *  Examples:
 *      "ncbi-file:/home/my-name/data-files"
 *      "ncbi-file://win-server/archive/secure/read12345?encrypted"
 *      "ncbi-file:///c/scanned-data/0001/file.sra?enc&pwfile=/c/Users/JamesMcCoy/ncbi.pwd"
 */
VFS_EXTERN rc_t CC VFSManagerMakePath ( struct VFSManager const * self,
    VPath ** new_path, const char * path_str, ... );
VFS_EXTERN rc_t CC VFSManagerVMakePath ( struct VFSManager const * self,
    VPath ** new_path, const char * path_fmt, va_list args );


/* MakeSysPath
 *  make a path object from an OS native filesystem path string
 *
 *  "new_path" [ OUT ] - return parameter for new path object
 *
 *  "sys_path" [ IN ] - a UTF-8 NUL-terminated string
 *  representing a native filesystem path
 *
 *  "wide_sys_path" [ IN ] - a wide NUL-terminated string
 *  representing a native filesystem path, where
 *  wchar_t is either USC-2 or UTF-32 depending upon libraries
 */
VFS_EXTERN rc_t CC VFSManagerMakeSysPath ( struct VFSManager const * self,
    VPath ** new_path, const char * sys_path );
VFS_EXTERN rc_t CC VFSManagerWMakeSysPath ( struct VFSManager const * self,
    VPath ** new_path, const wchar_t * wide_sys_path );


/* MakeAccPath - TEMPORARY
 *  takes a textual accession representation
 *  creates a VPath representing an accession
 *
 *  "new_path" [ OUT ] - return parameter for new path object
 *
 *  "acc" [ IN ] - a NUL-terminated ASCII fmt string
 */
VFS_EXTERN rc_t CC VFSManagerMakeAccPath ( struct VFSManager const * self,
    VPath ** new_path, const char * acc, ... );
VFS_EXTERN rc_t CC VFSManagerVMakeAccPath ( struct VFSManager const * self,
    VPath ** new_path, const char * fmt, va_list args );


/* MakeOidPath - TEMPORARY
 *  takes an integer oid
 *  creates a VPath representing an obj-id
 *
 *  "new_path" [ OUT ] - return parameter for new path object
 *
 *  "oid" [ IN ] - a non-zero object id
 */
VFS_EXTERN rc_t CC VFSManagerMakeOidPath ( struct VFSManager const * self,
    VPath ** new_path, uint32_t oid );


/* AddRef
 * Release
 *  ignores NULL references
 */
VFS_EXTERN rc_t CC VPathAddRef ( const VPath *self );
VFS_EXTERN rc_t CC VPathRelease ( const VPath *self );


/* IsFSCompatible
 *  asks if the path can be used with the OS' filesystems
 */
VFS_EXTERN bool CC VPathIsFSCompatible ( const VPath * self );


/* FromUri
 *  asks if the path was created from a formal URI
 */
VFS_EXTERN bool CC VPathFromUri ( const VPath * self );


/* Read*
 *  read the various parts
 *  copies out data into user-supplied buffer
 *
 *  "buffer" [ OUT ] and "buffer_size" [ IN ] - output buffer
 *  for data read. if sufficient space is available, the copy
 *  will be NUL-terminated.
 *
 *  "num_read" [ OUT, NULL OKAY ] - optional return parameter
 *  for the number of valid bytes in "buffer" after a successful
 *  read. on failure due to insufficient buffer, contains the
 *  number of bytes required for transfer.
 */
VFS_EXTERN rc_t CC VPathReadUri ( const VPath * self,
    char * buffer, size_t buffer_size, size_t * num_read );
VFS_EXTERN rc_t CC VPathReadScheme ( const VPath * self,
    char * buffer, size_t buffer_size, size_t * num_read );
VFS_EXTERN rc_t CC VPathReadAuth ( const VPath * self,
    char * buffer, size_t buffer_size, size_t * num_read );
VFS_EXTERN rc_t CC VPathReadHost ( const VPath * self,
    char * buffer, size_t buffer_size, size_t * num_read );
VFS_EXTERN rc_t CC VPathReadPortName ( const VPath * self,
    char * buffer, size_t buffer_size, size_t * num_read );
VFS_EXTERN rc_t CC VPathReadPath ( const VPath * self,
    char * buffer, size_t buffer_size, size_t * num_read );
VFS_EXTERN rc_t CC VPathReadQuery ( const VPath * self,
    char * buffer, size_t buffer_size, size_t * num_read );
VFS_EXTERN rc_t CC VPathReadParam ( const VPath * self, const char * param,
    char * buffer, size_t buffer_size, size_t * num_read );
VFS_EXTERN rc_t CC VPathReadFragment ( const VPath * self,
    char * buffer, size_t buffer_size, size_t * num_read );


/* MakeUri
 *  convert a VPath into a URI
 */
VFS_EXTERN rc_t CC VPathMakeUri ( const VPath * self,
    struct String const ** uri );


/* MakeString
 *  convert a VPath into a String
 *  respects original source of path,
 *  i.e. does not add scheme unnecessarily
 */
VFS_EXTERN rc_t CC VPathMakeString ( const VPath * self,
    struct String const ** str );


/* Get*
 *  retrieves internal parts
 *  returns pointers to internal String data
 *  Strings remain valid while "self" is valid
 */
VFS_EXTERN rc_t CC VPathGetScheme ( const VPath * self, struct String * str );
VFS_EXTERN rc_t CC VPathGetAuth ( const VPath * self, struct String * str );
VFS_EXTERN rc_t CC VPathGetHost ( const VPath * self, struct String * str );
VFS_EXTERN rc_t CC VPathGetPortName ( const VPath * self, struct String * str );
VFS_EXTERN uint16_t CC VPathGetPortNum ( const VPath * self );
VFS_EXTERN rc_t CC VPathGetPath ( const VPath * self, struct String * str );
VFS_EXTERN rc_t CC VPathGetQuery ( const VPath * self, struct String * str );
VFS_EXTERN rc_t CC VPathGetParam ( const VPath * self, const char * param, struct String * str );
VFS_EXTERN rc_t CC VPathGetFragment ( const VPath * self, struct String * str );
/* TEMPORARY */
VFS_EXTERN uint32_t CC VPathGetOid ( const VPath * self );

/* legacy support */
#define VPathMake LegacyVPathMake
VFS_EXTERN rc_t VPathMake ( VPath ** new_path, const char * posix_path );
#define VPathMakeFmt LegacyVPathMakeFmt
rc_t VPathMakeFmt ( VPath ** new_path, const char * fmt, ... );
#define VPathMakeVFmt LegacyVPathMakeVFmt
rc_t VPathMakeVFmt ( VPath ** new_path, const char * fmt, va_list args );
#define VPathMakeSysPath LegacyVPathMakeSysPath
VFS_EXTERN rc_t VPathMakeSysPath ( VPath ** new_path, const char * sys_path );


#ifdef __cplusplus
}
#endif

#endif /* _h_vfs_path_ */
