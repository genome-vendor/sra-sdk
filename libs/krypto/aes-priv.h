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
 */

#ifndef _h_krypto_aes_priv_
#define _h_krypto_aes_priv_

/*
 * This header file was written to integrate the public domain AES code
 * with the SRA project
 */

#define AES_ENCRYPT	1
#define AES_DECRYPT	0

#define AES_MAXNR (14)
#define AES_BLOCK_SIZE (16)


#ifdef  __cplusplus
extern "C" {
#endif

typedef uint8_t AES_BYTE;
typedef uint32_t AES_WORD;


typedef struct AES_KEY AES_KEY;

struct AES_KEY
{
    AES_WORD rd_key [sizeof (AES_WORD) * (AES_MAXNR + 1)];
    uint32_t rounds;
};


int AES_set_encrypt_key(const uint8_t *userKey, const uint32_t bits,
                        AES_KEY *key);

int AES_set_decrypt_key(const uint8_t *userKey, const uint32_t bits,
                        AES_KEY *key);

void AES_encrypt(const uint8_t *in, uint8_t *out,
                 const AES_KEY *key);
void AES_decrypt(const uint8_t *in, uint8_t *out,
                 const AES_KEY *key);


int AESx86_set_encrypt_key(const uint8_t *userKey, const uint32_t bits,
                        AES_KEY *key);

int AESx86_set_decrypt_key(const uint8_t *userKey, const uint32_t bits,
                        AES_KEY *key);

void AESx86_encrypt(const uint8_t *in, uint8_t *out,
                 const AES_KEY *key);
void AES_decrypt(const uint8_t *in, uint8_t *out,
                 const AES_KEY *key);


/* Most useful for testnig against Appendix A of the FIPS-197 document */
#define DEBUG_KEYEXP(msg) DBGMSG(DBG_AES,DBG_FLAG(DBG_AES_KEYEXP), msg)
#define DEBUG_CIPHER(msg) DBGMSG(DBG_AES,DBG_FLAG(DBG_AES_CIPHER), msg)
#define DEBUG_INVKEYEXP(msg) DBGMSG(DBG_AES,DBG_FLAG(DBG_AES_INVKEYEXP), msg)
#define DEBUG_INVCIPHER(msg) DBGMSG(DBG_AES,DBG_FLAG(DBG_AES_INVCIPHER), msg)

#define DEBUG_CIPHER_MVECTOR(M,v) DBGMSG(DBG_AES,DBG_FLAG(DBG_AES_CIPHER), \
                                        ("%s:\t%0.8x %0.8x %0.8x %0.8x\n",M, \
                                         v.columns[0],v.columns[1],v.columns[2],v.columns[3]))
#define DEBUG_INVCIPHER_MVECTOR(M,v) DBGMSG(DBG_AES,DBG_FLAG(DBG_AES_INVCIPHER), \
                                           ("%s:\t%0.8x %0.8x %0.8x %0.8x\n",M, \
                                            v.columns[0],v.columns[1],v.columns[2],v.columns[3]))

#define DEBUG_CIPHER_VECTOR(M,v) DBGMSG(DBG_AES,DBG_FLAG(DBG_AES_CIPHER), \
                                        ("%s:\t%0.8x %0.8x %0.8x %0.8x\n",M, \
                                         bswap_32(v.columns[0]),bswap_32(v.columns[1]),bswap_32(v.columns[2]),bswap_32(v.columns[3])))
#define DEBUG_INVCIPHER_VECTOR(M,v) DBGMSG(DBG_AES,DBG_FLAG(DBG_AES_INVCIPHER), \
                                        ("%s:\t%0.8x %0.8x %0.8x %0.8x\n",M, \
                                         bswap_32(v.columns[0]),bswap_32(v.columns[1]),bswap_32(v.columns[2]),bswap_32(v.columns[3])))

#ifdef  __cplusplus
}
#endif

#endif /* #ifndef _h_krypto_aes_priv_ */


