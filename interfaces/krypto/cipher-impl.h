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

#ifndef _h_krypto_cipher_impl_
#define _h_krypto_cipher_impl_

#include <krypto/extern.h>
#include <klib/refcount.h>
#include <klib/text.h>

#ifdef __cplusplus
extern "C" {
#endif


/*--------------------------------------------------------------------------
 * forwards
 */
typedef union KBlockCipher_vt KBlockCipher_vt;


/*--------------------------------------------------------------------------
 * CIPHER
 *  
 */

#ifndef CIPHER_IMPL
#define CIPHER_IMPL struct KBlockCipher
#endif

typedef struct KBlockCipher_vt_v1 KBlockCipher_vt_v1;
struct KBlockCipher_vt_v1
{
    /* version == 1.x */
    uint32_t maj;
    uint32_t min;

    /* start minor version == 0 */
    rc_t (CC * destroy)(CIPHER_IMPL * self);

    rc_t (CC * block_size)(const CIPHER_IMPL * self, size_t * bytes);

    rc_t (CC * key_size)(const CIPHER_IMPL * self, size_t * bytes);

    rc_t (CC * set_encrypt_key)(const CIPHER_IMPL * self,
                                void * key,
                                const char * user_key,
                                uint32_t user_key_bits);

    rc_t (CC * set_decrypt_key)(const CIPHER_IMPL * self,
                                void * key,
                                const char * user_key,
                                uint32_t user_key_bits);

    rc_t (CC * encrypt)(const CIPHER_IMPL * self, const void * in, void * out,
                        const void * key);
    rc_t (CC * decrypt)(const CIPHER_IMPL * self, const void * in, void * out,
                        const void * key);

    /* end minor version == 0 */

    /* ANY NEW ENTRIES MUST BE REFLECTED IN libs/krypto/cipher.c
       BY BOTH THE CORRESPONDING MESSAGE DISPATCH FUNCTION(s) AND
       VTABLE VALIDITY CHECKS IN CIPHERInit */
};

union KBlockCipher_vt
{
    KBlockCipher_vt_v1 v1;
};


#ifdef __cplusplus
}
#endif

#endif /* _h_krypto_cipher_impl_ */
