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
#include <krypto/cipher-impl.h>
#include <krypto/ciphermgr.h>


#include <klib/defs.h>
#include <klib/refcount.h>
#include <klib/rc.h>

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define rcBlockCipher rcNoTarg
#define rcCipher rcNoTarg


const char KBlockCipherClassName[] = "KBlockCipher";

/* ======================================================================
 */
/* methods from the virtual table */

/* ----------------------------------------------------------------------
 * Destroy
 *   base class destruction called during the derived class destruction
 */
KRYPTO_EXTERN rc_t CC KBlockCipherDestroy (KBlockCipher * self)
{
    assert (self);
    KRefcountWhack (&self->refcount, KBlockCipherClassName);
    return 0;
}


/* ----------------------------------------------------------------------
 * Whack
 *   start the destruction of a KBlockCIpher form the most most derived class 
 *   in the VT
 *
 *   This is called not directly like other vt methods but from class method
 *   Release which is not in the vt
 */
static rc_t KBlockCipherWhack (KBlockCipher * self)
{
    switch (self->vt->v1.maj)
    {
    case 1:
        return self->vt->v1.destroy (self);
    }
    return RC (rcKrypto, rcCipher, rcDestroying, rcInterface, rcBadVersion);
}


/* ----------------------------------------------------------------------
 * BlockSize
 *   How large is a block for this cipher?
 *
 */
KRYPTO_EXTERN rc_t CC KBlockCipherBlockSize (const KBlockCipher * self,
                                          size_t * bytes)
{
    if (self == NULL)
        return RC (rcKrypto, rcCipher, rcAccessing, rcSelf, rcNull);

    if (bytes == NULL)
        return RC (rcKrypto, rcCipher, rcAccessing, rcParam, rcNull);

    switch (self->vt->v1.maj)
    {
    case 1:
        return self->vt->v1.block_size (self, bytes);
    }
    return RC (rcKrypto, rcCipher, rcAccessing, rcInterface, rcBadVersion);
}


/* ----------------------------------------------------------------------
 * BlockSize
 *   How large is the stored key for this cipher?  Not the user key used
 *   to create this key (key schedule)
 *
 *   This is needed by KCipher to know how large the KCipher objecr is at
 *   allocation and to know how much of a buffer each decryption/encryption is
 */
KRYPTO_EXTERN rc_t CC KBlockCipherKeySize (const KBlockCipher * self,
                                           size_t * bytes)
{
    if (self == NULL)
        return RC (rcKrypto, rcCipher, rcAccessing, rcSelf, rcNull);

    if (bytes == NULL)
        return RC (rcKrypto, rcCipher, rcAccessing, rcParam, rcNull);

    switch (self->vt->v1.maj)
    {
    case 1:
        return self->vt->v1.key_size (self, bytes);
    }
    return RC (rcKrypto, rcCipher, rcAccessing, rcInterface, rcBadVersion);
}


/* ----------------------------------------------------------------------
 * SetEncryptKey
 *   The KCipher calls this to have the block cipher build an encryption
 *   key in the KCipher object
 *
 */
KRYPTO_EXTERN rc_t CC KBlockCipherSetEncryptKey (const KBlockCipher * self, 
                                                 void * encrypt_key,
                                                 const char * user_key,
                                                 uint32_t user_key_size)
{
    if (self == NULL)
        return RC (rcKrypto, rcCipher, rcAccessing, rcSelf, rcNull);

    if ((user_key == NULL)||(user_key_size == 0))
        return RC (rcKrypto, rcCipher, rcAccessing, rcParam, rcNull);

    switch (self->vt->v1.maj)
    {
    case 1:
        return self->vt->v1.set_encrypt_key (self, encrypt_key, user_key,
                                             user_key_size);
    }
    return RC (rcKrypto, rcCipher, rcAccessing, rcInterface, rcBadVersion);
}


/* ----------------------------------------------------------------------
 * SetDecryptKey
 *   The KCipher calls this to have the block cipher build an decryption
 *   key in the KCipher object
 *
 */
KRYPTO_EXTERN rc_t CC KBlockCipherSetDecryptKey (const KBlockCipher * self, 
                                                 void * decrypt_key,
                                                 const char * user_key,
                                                 uint32_t user_key_size)
{
    if (self == NULL)
        return RC (rcKrypto, rcCipher, rcAccessing, rcSelf, rcNull);

    if ((user_key == NULL)||(user_key_size == 0))
        return RC (rcKrypto, rcCipher, rcAccessing, rcParam, rcNull);

    switch (self->vt->v1.maj)
    {
    case 1:
        return self->vt->v1.set_decrypt_key (self, decrypt_key, user_key, user_key_size);
    }
    return RC (rcKrypto, rcCipher, rcAccessing, rcInterface, rcBadVersion);
}


/* ----------------------------------------------------------------------
 * Encrypt
 *
 *   Perform an encryption of a single block.  Chained modes and stream
 *   cipher modes will call this multiple times.
 *
 */
KRYPTO_EXTERN rc_t CC KBlockCipherEncrypt (const KBlockCipher * self,
                                           const void * in, void * out,
                                           const void * key)
{
    if (self == NULL)
        return RC (rcKrypto, rcCipher, rcAccessing, rcSelf, rcNull);

    if ((key == NULL)||(in == NULL)||(out == NULL))
        return RC (rcKrypto, rcCipher, rcAccessing, rcParam, rcNull);

    switch (self->vt->v1.maj)
    {
    case 1:
        return self->vt->v1.encrypt (self, in, out, key);
    }
    return RC (rcKrypto, rcCipher, rcAccessing, rcInterface, rcBadVersion);
}


/* ----------------------------------------------------------------------
 * Decrypt
 *
 *   Perform a decryption of a single block.  Chained modes and stream
 *   cipher modes will call this multiple times.
 *
 */
KRYPTO_EXTERN rc_t CC KBlockCipherDecrypt (const KBlockCipher * self,
                                           const void * in, void * out,
                                           const void * key)
{
    if (self == NULL)
        return RC (rcKrypto, rcCipher, rcAccessing, rcSelf, rcNull);

    if ((key == NULL)||(in == NULL)||(out == NULL))
        return RC (rcKrypto, rcCipher, rcAccessing, rcParam, rcNull);

    switch (self->vt->v1.maj)
    {
    case 1:
        return self->vt->v1.decrypt (self, in, out, key);
    }
    return RC (rcKrypto, rcCipher, rcAccessing, rcInterface, rcBadVersion);
}


/* ----------------------------------------------------------------------
 * Alloc
 *
 *   Allocate the space for an object.  The derived class has to pass in it's
 *   size plus the size of the derived classes instane name so this method can
 *   know how much to allocate
 */
KRYPTO_EXTERN rc_t CC KBlockCipherAlloc (KBlockCipher ** obj, size_t z)
{
    void * mem;

    mem = malloc (z);
    if (mem)
    {
        *obj = mem;
        return 0;
    }
    return RC (rcKrypto, rcBlockCipher, rcConstructing, rcMemory, rcExhausted);
}


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
                                        const char * name)
{
    KRefcountInit (&self->refcount, 1, KBlockCipherClassName, "init", name);
    self->vt = vt;
    self->name = (char*)(self+1);
    strcpy (self->name, name);
    return 0;
}


/* ----------------------------------------------------------------------
 * AddRef
 *   add a new owner to this class.  THis will mean another instance of 
 *   KCipher used this Block Cipher
 */
KRYPTO_EXTERN rc_t CC
KBlockCipherAddRef (const KBlockCipher * self)
{
    if (self)
    {
        switch (KRefcountAdd (&self->refcount, KBlockCipherClassName))
        {
        case krefLimit:
            return RC (rcKrypto, rcCipher, rcAttaching, rcRange, rcExcessive);
        }
    }
    return 0;
}


/* ----------------------------------------------------------------------
 * Release
 *   
 */
KRYPTO_EXTERN rc_t CC
KBlockCipherRelease (const KBlockCipher * self)
{
    if ( self != NULL )
    {
        switch (KRefcountDrop ( & self -> refcount, KBlockCipherClassName))
        {
        case krefWhack:
            return KBlockCipherWhack ((KBlockCipher *)self);
        case krefLimit:
            return RC ( rcDB, rcColumn, rcReleasing, rcRange, rcExcessive );
        }
    }
    return 0;
}


