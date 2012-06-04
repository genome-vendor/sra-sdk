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

#ifndef _h_krypto_cipher_priv_
#define _h_krypto_cipher_priv_

#include <krypto/extern.h>
#include <krypto/cipher.h>

#include <klib/refcount.h>
#include <klib/text.h>

struct KCipherManager;
union KBlockCipher_vt;

/* ======================================================================
 */
typedef struct KBlockCipher KBlockCipher;

struct KBlockCipher
{
    KRefcount refcount;
    const union KBlockCipher_vt * vt;
    char * name;
    bool uses_openssl;
};

/* ----------------------------------------------------------------------
 * BlockSize
 *   How large is a block for this cipher?
 *
 */
KRYPTO_EXTERN rc_t CC KBlockCipherBlockSize (const KBlockCipher * self,
                                             size_t * bytes);


/* ----------------------------------------------------------------------
 * KeySize
 *   How large is the stored key for this cipher?  Not the user key used
 *   to create this key (key schedule)
 *
 *   This is needed by KCipher to know how large the KCipher objecr is at
 *   allocation and to know how much of a buffer each decryption/encryption is
 */
KRYPTO_EXTERN rc_t CC KBlockCipherKeySize (const KBlockCipher * self,
                                           size_t * bytes);


/* ----------------------------------------------------------------------
 * SetEncryptKey
 *   The KCipher calls this to have the block cipher build an encryption
 *   key in the KCipher object
 *
 */
KRYPTO_EXTERN rc_t CC KBlockCipherSetEncryptKey (const KBlockCipher * self, 
                                                 void * encrypt_key,
                                                 const char * user_key,
                                                 uint32_t user_key_size);


/* ----------------------------------------------------------------------
 * SetDecryptKey
 *   The KCipher calls this to have the block cipher build an decryption
 *   key in the KCipher object
 *
 */
KRYPTO_EXTERN rc_t CC KBlockCipherSetDecryptKey (const KBlockCipher * self, 
                                                 void * decrypt_key,
                                                 const char * user_key,
                                                 uint32_t user_key_size);


/* ----------------------------------------------------------------------
 * Encrypt
 *
 *   Perform an encryption of a single block.  Chained modes and stream
 *   cipher modes will call this multiple times.
 *
 */
KRYPTO_EXTERN rc_t CC KBlockCipherEncrypt (const KBlockCipher * self,
                                           const void * in, void * out,
                                           const void * key);


/* ----------------------------------------------------------------------
 * Decrypt
 *
 *   Perform a decryption of a single block.  Chained modes and stream
 *   cipher modes will call this multiple times.
 *
 */
KRYPTO_EXTERN rc_t CC KBlockCipherDecrypt (const KBlockCipher * self,
                                           const void * in, void * out,
                                           const void * key);


/* ----------------------------------------------------------------------
 * Alloc
 *
 *   Allocate the space for an object.  The derived class has to pass in it's
 *   size plus the size of the derived classes instane name so this method can
 *   know how much to allocate
 */
KRYPTO_EXTERN rc_t CC KBlockCipherAlloc (KBlockCipher ** obj, size_t z);


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
KRYPTO_EXTERN rc_t CC KBlockCipherInit (KBlockCipher * self,
                                        const union KBlockCipher_vt * vt,
                                        const char * name);


/* ----------------------------------------------------------------------
 * Destroy
 *   base class destruction called during the derived class destruction
 */
KRYPTO_EXTERN rc_t CC KBlockCipherDestroy (KBlockCipher * self);


/* ----------------------------------------------------------------------
 * AddRef
 *   add a new owner to this class.  THis will mean another instance of 
 *   KCipher used this Block Cipher
 */
KRYPTO_EXTERN rc_t CC KBlockCipherAddRef (const KBlockCipher * self);


/* ----------------------------------------------------------------------
 * Release
 *   
 */
KRYPTO_EXTERN rc_t CC KBlockCipherRelease (const KBlockCipher * self);



/* ======================================================================
 */

struct KCipher
{
    KRefcount refcount;
    const KBlockCipher * block_cipher;
    struct String name;

    void * encrypt_key;
    void * decrypt_key;
    void * encrypt_ivec;
    void * decrypt_ivec;
};


rc_t CC KCipherAlloc (KCipher ** obj,
                      const KBlockCipher * block_cipher);

rc_t CC KCipherInit (KCipher * self,
                     const KBlockCipher * block_cipher);


rc_t KAESCipherInit (KBlockCipher * self,
                     const struct KCipherManager * mgr);

rc_t KAESCipherAlloc (KBlockCipher ** obj);

rc_t KAESx86CipherInit (KBlockCipher * self,
                     const struct KCipherManager * mgr);

rc_t KAESx86CipherAlloc (KBlockCipher ** obj);

rc_t KAESNCBICipherInit (KBlockCipher * self,
                        const struct KCipherManager * mgr);

rc_t KAESNCBICipherAlloc (KBlockCipher ** obj);

rc_t KAESNCBIniCipherInit (KBlockCipher * self,
                        const struct KCipherManager * mgr);

rc_t KAESNCBIniCipherAlloc (KBlockCipher ** obj);



KRYPTO_EXTERN rc_t  KCipherManagerMakeAES (struct KCipherManager * self, KCipher ** obj);




#ifdef __cplusplus
}
#endif


#endif /* #ifndef _h_krypto_cipher_priv_ */
