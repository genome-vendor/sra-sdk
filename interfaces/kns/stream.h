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

#ifndef _h_kns_stream_
#define _h_kns_stream_

#ifndef _h_kns_extern_
#include <kns/extern.h>
#endif

#ifndef _h_klib_defs_
#include <klib/defs.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif


/*--------------------------------------------------------------------------
 * KStream
 *  the stream is defined to have no concept of size,
 *  and to not support any form of random access
 */
typedef struct KStream KStream;


/* AddRef
 * Release
 *  ignores NULL references
 */
KNS_EXTERN rc_t CC KStreamAddRef ( const KStream *self );
KNS_EXTERN rc_t CC KStreamRelease ( const KStream *self );


/* Read
 *  read data from stream
 *
 *  "buffer" [ OUT ] and "bsize" [ IN ] - return buffer for read
 *
 *  "num_read" [ OUT ] - return parameter giving number of bytes
 *  actually read. when returned value is zero and return code is
 *  also zero, interpreted as end of stream.
 */
KNS_EXTERN rc_t CC KStreamRead ( const KStream *self,
    void *buffer, size_t bsize, size_t *num_read );

/* ReadAll
 *  read from stream until "bsize" bytes have been retrieved
 *  or until end-of-input
 *
 *  "buffer" [ OUT ] and "bsize" [ IN ] - return buffer for read
 *
 *  "num_read" [ OUT ] - return parameter giving number of bytes
 *  actually read. when returned value is zero and return code is
 *  also zero, interpreted as end of stream.
 */
KNS_EXTERN rc_t CC KStreamReadAll ( const KStream *self,
    void *buffer, size_t bsize, size_t *num_read );


/* Write
 *  send data to stream
 *
 *  "buffer" [ IN ] and "size" [ IN ] - data to be written
 *
 *  "num_writ" [ OUT, NULL OKAY ] - optional return parameter
 *  giving number of bytes actually written
 */
KNS_EXTERN rc_t CC KStreamWrite ( KStream *self,
    const void *buffer, size_t size, size_t *num_writ );

/* WriteAll
 *  write to stream until "size" bytes have been transferred
 *  or until no further progress can be made
 *
 *  "buffer" [ IN ] and "size" [ IN ] - data to be written
 *
 *  "num_writ" [ OUT, NULL OKAY ] - optional return parameter
 *  giving number of bytes actually written
 */
KNS_EXTERN rc_t CC KStreamWriteAll ( KStream *self,
    const void *buffer, size_t size, size_t *num_writ );


/* MakeStdIn
 *  creates a read-only stream on stdin
 */
KNS_EXTERN rc_t CC KStreamMakeStdIn ( const KStream **std_in );

/* MakeStdOut
 * MakeStdErr
 *  creates a write-only stream on stdout or stderr
 */
KNS_EXTERN rc_t CC KStreamMakeStdOut ( KStream **std_out );
KNS_EXTERN rc_t CC KStreamMakeStdErr ( KStream **std_err );


#ifdef __cplusplus
}
#endif

#endif /* _h_kns_stream_ */
