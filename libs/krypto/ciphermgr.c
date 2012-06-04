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

#define USE_AES_NI false

#include <krypto/extern.h>
#include <krypto/ciphermgr.h>
#include <krypto/cipher.h>
#include <krypto/cipher-priv.h>

#include <klib/rc.h>
#include <klib/refcount.h>
#include <klib/debug.h>

#include <stdlib.h>
#include <assert.h>

#ifdef _DEBUGGING
#define MGR_DEBUG(msg) DBGMSG(DBG_KFS,DBG_FLAG(DBG_KFS_MGR), msg)
#else
#define MGR_DEBUG(msg)
#endif



static const char kciphermanager_classname [] = "KCipherManager";


/*--------------------------------------------------------------------------
 * KCipherManager
 */
/* currently expected to be a singleton and not use a vtable but
 * be fully fleashed out here */
static 
KCipherManager * singleton = NULL;

struct KCipherManager
{
    KRefcount refcount;

    /* reference counting is from the cipher to the manager
     * not form the manager to the cipher
     *
     * when we decide to support more than AES and especially when we
     * decide to allow foreign block ciphers we'll need a better
     * structure here...
     */
    const KBlockCipher * block_ciphers [kcipher_count];
};


static
rc_t KCipherManagerAlloc (KCipherManager ** ppobj)
{
    KCipherManager * pobj;

    assert (ppobj);

    pobj = calloc (sizeof *pobj, 1);
    if (pobj)
    {
        *ppobj = pobj;
        return 0;
    }

    *ppobj = NULL;
    return RC (rcKrypto, rcMgr, rcConstructing, rcMemory, rcExhausted);
}


static
rc_t KCipherManagerInit (KCipherManager * self)
{
    uint32_t index;

    KRefcountInit (&self->refcount, 1, kciphermanager_classname, "init",
                   "singleton");

    for (index = 0; index < kcipher_count; ++index)
        self->block_ciphers[index] = NULL;

    return 0;
}


/* Destroy
 *  destroy
 */
LIB_EXPORT rc_t CC KCipherManagerDestroy ( KCipherManager *self )
{
    rc_t rc = 0;

    if ( self == NULL )
        rc = RC ( rcKrypto, rcMgr, rcDestroying, rcSelf, rcNull );
    else
    {
        assert (self == singleton);

        /* no return value */
        KRefcountWhack (&self->refcount, kciphermanager_classname);

        {
            unsigned ix;
            
            for (ix = 0; ix < kcipher_count; ++ix)
                if (self->block_ciphers[ix])
                {
                    KBlockCipherRelease (self->block_ciphers[ix]);
                    self->block_ciphers[ix] = NULL;
                }
        }
        free (self);
        singleton = NULL;
    }
    return rc;
}


/* AddRef
 *  creates a new reference
 *  ignores NULL references
 */
LIB_EXPORT rc_t CC KCipherManagerAddRef ( const KCipherManager *self )
{
    if (self != NULL)
    {
        switch (KRefcountAdd (&self->refcount, kciphermanager_classname))
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
 *  discard reference to manager
 *  ignores NULL references
 */
LIB_EXPORT rc_t CC KCipherManagerRelease ( const KCipherManager *self )
{
    rc_t rc = 0;
    if (self != NULL)
    {
        switch (KRefcountDrop (&self->refcount, kciphermanager_classname))
        {
        case krefOkay:
        case krefZero:
            break;
        case krefWhack:
            rc = KCipherManagerDestroy ((KCipherManager*)self);
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


/* Make
 *  we have a shared singleton for the cipher manager
 *  first call actually makes the managerblo
 *  subsequent calls get added references
 */
LIB_EXPORT rc_t CC KCipherManagerMake (KCipherManager ** mgr)
{
    rc_t rc = 0;

    if (mgr == NULL)
        return RC (rcKrypto, rcMgr, rcConstructing, rcSelf, rcNull);
    *mgr = NULL;

    if (singleton)
    {
        rc = KCipherManagerAddRef (singleton);
        if (rc == 0)
        {
            *mgr = singleton;
            return 0;
        }
    }
    else
    {
        rc = KCipherManagerAlloc (&singleton);
        if (rc == 0)
        {
            rc = KCipherManagerInit (singleton);
            if (rc == 0)
            {
                *mgr = singleton;
                return 0;
            }

            KCipherManagerDestroy (singleton);
            singleton = NULL;
        }
    }
    return rc;
}

static
rc_t KCipherManagerMakeBlockCipher (const KCipherManager * self,
                                    const KBlockCipher ** bcipher,
                                    kcipher_type type)
{
    rc_t rc = 0;

    assert (self);
    assert (bcipher);

    *bcipher = NULL;

    if (type >= kcipher_count)
        rc = RC (rcKrypto, rcMgr, rcConstructing, rcParam, rcInvalid);

    else
    {
        const KBlockCipher * c_block_cipher;

        c_block_cipher = self->block_ciphers[type];

        if (c_block_cipher == NULL)
        {
            KBlockCipher * block_cipher;
            switch (type)
            {
            default:
                return RC (rcKrypto, rcMgr, rcConstructing, rcParam, rcInvalid);

            case kcipher_AES:
                rc = KAESCipherAlloc (&block_cipher);
                if (rc == 0)
                {
                    rc = KAESCipherInit (block_cipher, self);
                    if (rc == 0)
                    {
                        ((KCipherManager*)self)->block_ciphers[type] = block_cipher;
                        c_block_cipher = block_cipher;
                    }
                    else
                        free (block_cipher);
                }
                break;
#if USE_NCBI_AES
            case kcipher_AES_ncbi_ni:
#if USE_AES_NI
                rc = KAESNCBIniCipherAlloc (&block_cipher);
                if (rc == 0)
                {
                    rc = KAESNCBIniCipherInit (block_cipher, self);
                    if (rc == 0)
                    {
                        ((KCipherManager*)self)->block_ciphers[type] = block_cipher;
                        c_block_cipher = block_cipher;
                    }
                    else
                        free (block_cipher);
                }
                break;
#endif
            case kcipher_AES_ncbi:
                rc = KAESNCBICipherAlloc (&block_cipher);
                if (rc == 0)
                {
                    rc = KAESNCBICipherInit (block_cipher, self);
                    if (rc == 0)
                    {
                        ((KCipherManager*)self)->block_ciphers[type] = block_cipher;
                        c_block_cipher = block_cipher;
                    }
                    else
                        free (block_cipher);
                }
                break;
#endif
            }
        }

        if (rc == 0)
        {
            rc = KBlockCipherAddRef (c_block_cipher);
            if (rc == 0)
            {
                *bcipher = c_block_cipher;
                return 0;
            }
        }
    }
    return rc;
}


static
rc_t KCipherManagerMakeCipherInt (const KCipherManager *self,
                                  const KBlockCipher * pbc,
                                  KCipher ** pcipher)
{
    KCipher * pc;
    rc_t rc;

    assert (self);
    assert (pbc);
    assert (pcipher);

    *pcipher = NULL;

    rc = KCipherAlloc (&pc, pbc);
    if (rc == 0)
    {
        rc = KCipherInit (pc, pbc);
    
        if (rc == 0)
        {
            *pcipher = pc;
            return 0;
        }
        KBlockCipherRelease (pbc);
    }

    return rc;
}


LIB_EXPORT
rc_t CC KCipherManagerMakeCipher (const KCipherManager * self,
                                  struct KCipher ** pcipher,
                                  kcipher_type type)
{
    const KBlockCipher * pbc;
    rc_t rc;
    
    if (self == NULL)
        return RC (rcKrypto, rcMgr, rcConstructing, rcSelf, rcNull);

    if (pcipher == NULL)
        return RC (rcKrypto, rcMgr, rcConstructing, rcParam, rcNull);

    rc = KCipherManagerMakeBlockCipher (self, &pbc, type);
    if (rc == 0)
    {
        rc = KCipherManagerMakeCipherInt (self, pbc, pcipher);
        KBlockCipherRelease (pbc);
    }
    return rc;
}



