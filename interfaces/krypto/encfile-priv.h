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

#ifndef _h_krypto_encfile_priv_
#define _h_krypto_encfile_priv_

#ifndef _h_krypto_extern_
#include <krypto/extern.h>
#endif

#ifndef _h_klib_defs_
#include <klib/defs.h>
#endif

#include <krypto/encfile.h>

#ifdef __cplusplus
extern "C" {
#endif

#define DEBUG_STS(msg)     DBGMSG(DBG_KRYPTO,DBG_FLAG(DBG_KRYPTO_STS),msg)
#define DEBUG_CFG(msg)     DBGMSG(DBG_KRYPTO,DBG_FLAG(DBG_KRYPTO_CFG_,msg)
#define DEBUG_ENCRYPT(msg) DBGMSG(DBG_KRYPTO,DBG_FLAG(DBG_KRYPTO_ENCRYPT),msg)
#define DEBUG_DECRYPT(msg) DBGMSG(DBG_KRYPTO,DBG_FLAG(DBG_KRYPTO_DECRYPT),msg)

/* -----
 * Encrypted file structure:
 *   - File Header
 *   - one or more data blocks
 *   - file footer
 *
 * File Header:
 *   - file signature "NCBInenc"
 *   - byte order flag 
 *   - version
 *
 * Data Block:
 *   - rkey - randomly generated 32 byte key for this block (encrypted using user key)
 *   - encrypted block of 32768 bytes size above says how many are really used (encrypted using rkey and salt above)
 *   - block-offset + (valid bytes in block % 32768)
 *   - crc-32 (includes phantom block start offset as initial seed)
 *
 * File Footer:
 *   - footer signature "foot"
 *   - checksum of crcs
 */


/* ----------------------------------------------------------------------
 * Header - the file header
 * all constant values for the first version
 */
typedef char KEncFileSig [8];
typedef uint32_t Endian_t;
typedef uint32_t KEncFileVersion;


typedef struct KEncFileHeader KEncFileHeader;
struct KEncFileHeader
{
    KEncFileSig     file_sig;   /* "NCBInenc" */
    Endian_t        byte_order; /* do we byte swap on read? */
    KEncFileVersion version;    /* simple incrementation starting at 1 */
};


/* ----------------------------------------------------------------------
 * KEncFileBlock
 *    The body of the file is blocks containing a portion of the decrypted
 *    file.  These are an ordered sequence with the last block being the
 *    same size as the rest but with only some of the data portion being
 *    a part of the file.
 *
 *    An crypted file is longer than an unencrypted file by 
 *       a constant: the lengths of the header and the footer
 *       proportinally by the length of the block key and crc
 */

/* -----
 * Key  the header for an encrypted block
 *
 * when initialized the first 38 bytes should be set to random data.
 * valid is a count of how many bytes in the block are valid data
 * offset is the offset of this block with in the decrypted file
 */
typedef uint8_t KEncFileKey [32];


/* -----
 * We sized the data portion of a block to match the KPageFile
 * structure allowing a KBufFile in front of a KEncFile to
 * operate in a fairly efficient manner
 */
#define ENC_DATA_BLOCK_SIZE     (32*1024)
typedef uint8_t KEncFileData [ENC_DATA_BLOCK_SIZE];

typedef uint16_t KEncFileOffValid;

typedef uint64_t KEncFileBlockId;
typedef uint16_t KEncFileBlockValid;

/* -----
 * we use the same 32 bit CRC as the rest of the project
 */
typedef uint32_t KEncFileCRC;


/*
 * NOTE:
 * The size of data + u + id + crc + crc_copy must remain divisible
 * by the size of key
 */
typedef struct KEncFileBlock KEncFileBlock;
struct KEncFileBlock
{
    KEncFileKey         key;  /* encrypted with the user key */
    KEncFileData        data; /* encrypted with block key */
    union
    {
        KEncFileBlockValid  valid; /* obscured and encrypted */
        uint8_t bytes [16];        /* mostly fill */
    } u;
    KEncFileBlockId     id;        /* plain text */
    KEncFileCRC         crc;       /* plain text */
    KEncFileCRC         crc_copy;  /* plain text */
};


/* ----------------------------------------------------------------------
 * Foot - the ending of an encrypted file: 
 *   these are in plan text for non-decryption validation of the whole file
 */
typedef uint64_t KEncFileFooter_t;
typedef struct KEncFileFooter KEncFileFooter;
struct KEncFileFooter
{
    KEncFileFooter_t block_count;  /* how many blocks do we have? */
    KEncFileFooter_t crc_checksum; /* sum of crcs of all blocks */
};

/* ----------
 * Update mode is read/write mode where seeking within the file is allowed.
 *
 * NOTE this is in the private interface because it is not actually working
 * yet.
 */
KRYPTO_EXTERN rc_t CC KEncFileMakeUpdate (struct KFile ** pself, 
                                          struct KFile * encrypted,
                                          const struct KKey * key);


#ifdef __cplusplus
}
#endif

#endif /* _h_krypto_encfile_priv_ */
