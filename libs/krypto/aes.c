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
#include <krypto/extern.h>
#include <krypto/cipher.h>
#include <krypto/cipher-priv.h>
#include "aes-priv.h"

#include <klib/rc.h>

#define rcCipher rcNoTarg
#define rcBlockCipher rcNoTarg

#include <assert.h>

typedef struct KBlockCipher KAESBlockCipher;

/* instantiate for class functions that use the VT */
#define CIPHER_IMPL struct KBlockCipher
#include <krypto/cipher-impl.h>

static const char KAESCipherClassName[] = "KAESCipher";

/* trying to force the CC stuff to all work out */


static
rc_t CC KAESCipherBlockSize (const struct KBlockCipher * self,
                             size_t * bytes)
{
    assert (self);
    assert (bytes);

    *bytes = AES_BLOCK_SIZE;

    return 0;
}


/* ----------------------------------------------------------------------
 * KeySize
 *   How large is the stored key for this cipher?  Not the user key used
 *   to create this key (key schedule)
 *
 *   This is needed by KCipher to know how large the KCipher objecr is at
 *   allocation and to know how much of a buffer each decryption/encryption is
 */
KRYPTO_EXTERN rc_t CC KAESCipherKeySize (const struct KBlockCipher * self,
                                         size_t * bytes)
{
    assert (self);
    assert (bytes);

    *bytes = sizeof (AES_KEY);

    return 0;
}


/* ----------------------------------------------------------------------
 * SetEncryptKey
 *   The KCipher calls this to have the block cipher build an encryption
 *   key in the KCipher object
 *
 */
KRYPTO_EXTERN rc_t CC KAESCipherSetEncryptKey (const KBlockCipher * self, 
                                               void * encrypt_key,
                                               const char * user_key,
                                               uint32_t user_key_size)
{
    int iii;
    rc_t rc;

    assert (self);
    assert (encrypt_key);
    assert (user_key);
    assert (user_key_size != 0);

    iii = AES_set_encrypt_key((const unsigned char *)user_key, (int)user_key_size * 8,
                              encrypt_key);
    switch (iii)
    {
    default: /* not in the code when this was written */
        rc = RC (rcKrypto, rcCipher, rcUpdating, rcEncryptionKey, rcUnknown);
        break;

    case -1: /* bad parameters */
        rc = RC (rcKrypto, rcCipher, rcUpdating, rcParam, rcInvalid);
        break;

    case -2: /* bad bit count */
        rc = RC (rcKrypto, rcCipher, rcUpdating, rcParam, rcIncorrect);
        break;

    case 0:
        /* all is copasectic */
        rc = 0;
        break;
    }
    return rc;
}


/* ----------------------------------------------------------------------
 * SetDecryptKey
 *   The KCipher calls this to have the block cipher build an decryption
 *   key in the KCipher object
 *
 */
KRYPTO_EXTERN rc_t CC KAESCipherSetDecryptKey (const KBlockCipher * self, 
                                               void * decrypt_key,
                                               const char * user_key,
                                               uint32_t user_key_size)
{
    rc_t rc;
    int iii;

    assert (self);
    assert (decrypt_key);
    assert (user_key);
    assert (user_key_size != 0);

    iii = AES_set_decrypt_key((const unsigned char *)user_key, (int)user_key_size * 8,
                              decrypt_key);
    switch (iii)
    {
    default: /* not in the code when this was written */
        rc = RC (rcKrypto, rcCipher, rcUpdating, rcEncryptionKey, rcUnknown);
        break;

    case -1: /* bad parameters */
        rc = RC (rcKrypto, rcCipher, rcUpdating, rcParam, rcInvalid);
        break;

    case -2: /* bad bit count */
        rc = RC (rcKrypto, rcCipher, rcUpdating, rcParam, rcIncorrect);
        break;

    case 0:
        /* all is copasectic */
        rc = 0;
        break;
    }
    return rc;

}


/* ----------------------------------------------------------------------
 * Encrypt
 *
 *   Perform an encryption of a single block.  Chained modes and stream
 *   cipher modes will call this multiple times.
 *
 */
KRYPTO_EXTERN rc_t CC KAESCipherEncrypt (const KBlockCipher * self,
                                         const void * in, void * out,
                                         const void * key)
{
    assert (self);
    assert (in);
    assert (out);
    assert (key);

    AES_encrypt (in, out, key);
    return 0;
}


/* ----------------------------------------------------------------------
 * Decrypt
 *
 *   Perform a decryption of a single block.  Chained modes and stream
 *   cipher modes will call this multiple times.
 *
 */
KRYPTO_EXTERN rc_t CC KAESCipherDecrypt (const KBlockCipher * self,
                                           const void * in, void * out,
                                           const void * key)
{
    assert (self);
    assert (in);
    assert (out);
    assert (key);

    AES_decrypt (in, out, key);
    return 0;
}




/* ----------------------------------------------------------------------
 * Alloc
 *
 *   Allocate the space for an object.  The derived class has to pass in it's
 *   size plus the size of the derived classes instane name so this method can
 *   know how much to allocate
 */
rc_t KAESCipherAlloc (KBlockCipher ** obj)
{
    return KBlockCipherAlloc ((KBlockCipher**)obj,
                              sizeof ** obj +
                              sizeof KAESCipherClassName);
}


static const
KBlockCipher_vt_v1 aes_vt = 
{
    1, 0,

    KBlockCipherDestroy,
    KAESCipherBlockSize,
    KAESCipherKeySize,
    KAESCipherSetEncryptKey,
    KAESCipherSetDecryptKey,
    KAESCipherEncrypt,
    KAESCipherDecrypt
};


/* ----------------------------------------------------------------------
 * Init
 *
 *   Initialize the fields of this object.  The derived class will call this
 *   during it's initialization.
 *
 * self      object to initialze
 * vt        the virtual table of the derived class
 * mgr       the cipher manager that is the construction factory block cipher
 *           objects hold references to the manager while the manager merely
 *           points at the block ciphers when all block ciphers are destroyed
 *           the manager loses its references and it too can be destroyed if not
 *           held elsewhere
 * name      ASCIZ c-string the name of this class
 */
rc_t KAESCipherInit (KBlockCipher * self,
                     const struct KCipherManager * mgr)
{
    assert (self);
    assert (mgr);
    return KBlockCipherInit ((KBlockCipher*)self, 
                             (const union KBlockCipher_vt *)&aes_vt,
                             KAESCipherClassName);
}

