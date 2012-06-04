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

#include <klib/refcount.h>
#include <klib/rc.h>
#include <klib/out.h>

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define rcBlockCipher rcNoTarg
#define rcCipher rcNoTarg


#define MAX_BLOCK_SIZE (1024)

/* ======================================================================
 */
const char KCipherClassName[] = "KCipher";



KRYPTO_EXTERN rc_t CC KCipherAddref (const KCipher * self)
{
    if (self)
    {
        switch (KRefcountAdd (&self->refcount, KCipherClassName))
        {
        case krefLimit:
            return RC (rcKrypto, rcCipher, rcAttaching, rcRange, rcExcessive);
        }
    }
    return 0;
}

static
rc_t KCipherWhack (KCipher * self)
{
    rc_t rc;

    rc = KBlockCipherRelease (self->block_cipher);

    free (self);
    return rc;
}


KRYPTO_EXTERN rc_t CC KCipherRelease (const KCipher * self)
{
    if ( self != NULL )
    {
        switch ( KRefcountDrop ( & self -> refcount, KCipherClassName ) )
        {
        case krefWhack:
            return KCipherWhack ((KCipher *)self);
        case krefLimit:
            return RC ( rcDB, rcColumn, rcReleasing, rcRange, rcExcessive );
        }
    }
    return 0;
}


KRYPTO_EXTERN rc_t CC KCipherBlockSize (const KCipher * self, size_t * bytes)
{
    if (self == NULL)
        return RC (rcKrypto, rcCipher, rcAccessing, rcSelf, rcNull);

    if (bytes == NULL)
        return RC (rcKrypto, rcCipher, rcAccessing, rcParam, rcNull);

    return (KBlockCipherBlockSize (self->block_cipher, bytes));
}



KRYPTO_EXTERN rc_t CC KCipherSetEncryptKey (KCipher * self, const void * user_key, size_t user_key_size)
{
    if (self == NULL)
        return RC (rcKrypto, rcCipher, rcAccessing, rcSelf, rcNull);

    if ((user_key == NULL)||(user_key_size == 0))
        return RC (rcKrypto, rcCipher, rcAccessing, rcParam, rcNull);

    return KBlockCipherSetEncryptKey (self->block_cipher, self->encrypt_key,
                                      user_key, user_key_size);
}


KRYPTO_EXTERN rc_t CC KCipherSetDecryptKey (KCipher * self, const void * user_key, size_t user_key_size)
{
    if (self == NULL)
        return RC (rcKrypto, rcCipher, rcAccessing, rcSelf, rcNull);

    if ((user_key == NULL)||(user_key_size == 0))
        return RC (rcKrypto, rcCipher, rcAccessing, rcParam, rcNull);

    return KBlockCipherSetDecryptKey (self->block_cipher, self->decrypt_key,
                                      user_key, user_key_size);
}


/*
 * Set the ivec (Initialization vector or feedback) for the cipher
 * this is done automatically for the longer runs defined below.
 *
 * the size of ivec  must match KCipherBlockSize
 *
 * the ivec is copied into the cipher not used in place
 */
KRYPTO_EXTERN rc_t CC KCipherSetEncryptIVec (KCipher * self, const void * ivec)
{
    if (self == NULL)
        return RC (rcKrypto, rcCipher, rcAccessing, rcSelf, rcNull);

    else if (ivec == NULL)
        return RC (rcKrypto, rcCipher, rcAccessing, rcParam, rcNull);

    else
    {
        size_t bytes;
        rc_t rc;

        rc = KBlockCipherBlockSize (self->block_cipher, &bytes);
        if (rc == 0)
            memcpy (self->encrypt_ivec, ivec, bytes);
        return rc;
    }
}


KRYPTO_EXTERN rc_t CC KCipherSetDecryptIVec (KCipher * self, const void * ivec)
{
    if (self == NULL)
        return RC (rcKrypto, rcCipher, rcAccessing, rcSelf, rcNull);

    else if (ivec == NULL)
        return RC (rcKrypto, rcCipher, rcAccessing, rcParam, rcNull);

    else
    {
        size_t bytes;
        rc_t rc;

        rc = KBlockCipherBlockSize (self->block_cipher, &bytes);
        if (rc == 0)
            memcpy (self->decrypt_ivec, ivec, bytes);
        return rc;
    }
}



/*
 * 'in' can equal 'out'
 */
KRYPTO_EXTERN rc_t CC KCipherEncrypt (const KCipher * self, const void * in,
                                      void * out)
{
    rc_t rc;

    if (self == NULL)
        return RC (rcKrypto, rcCipher, rcAccessing, rcSelf, rcNull);

    else if ((in == NULL) || (out == NULL))
        return RC (rcKrypto, rcCipher, rcAccessing, rcParam, rcNull);

    else
    {
        rc = KBlockCipherEncrypt (self->block_cipher, in, out,
                                  self->encrypt_key);
    }
    return rc;
}


KRYPTO_EXTERN rc_t CC KCipherDecrypt (const KCipher * self, const void * in,
                                      void * out)
{
    rc_t rc;

    if (self == NULL)
        return RC (rcKrypto, rcCipher, rcAccessing, rcSelf, rcNull);

    else if ((in == NULL) || (out == NULL))
        return RC (rcKrypto, rcCipher, rcAccessing, rcParam, rcNull);

    else
    {
        rc = KBlockCipherDecrypt (self->block_cipher, in, out,
                                  self->decrypt_key);
    }
    return rc;
}


/* ====================
 * longer runs of multiple blocks.
 *
 * The algorithms are well defined and standard in most cases
 *
 * These aremore or elss equivalent to class functions as they do not depend upon
 * the operation of the cipher and the algorithms are independent of anything about
 * the cipher other than its block size.
 *
 * PT: plain text block
 * CT: cipher text block
 * EK: encryption key
 * DK: decryption key (might be sthe same as EK)
 * ENC: encrypt cipher function on a block using a key
 * DEC: decrypt cipher function on a block using a key
 * IV: initialization vector - used as feedback for chaining
 * N:  number used once (nonce)
 * FB: feedback is the next IV in a chained/feedback mode
 */

typedef rc_t (CC * BlockFunc)(const KBlockCipher * self, const void * in, void * out,
                              const void * key);


/* -----
 * NOTE:
 * 'in' can be the same as 'out' but other overlaps are dangers as a block at a
 * time is written. The code does not look for overlaps at this point.
 */

/* ----------
 * Electronic Code Book - simple cipher with no chaining feedback  just iterate
 * simple encrypt/decrypt with the plain, text, cipher text and key/
 *
 * CT = ENC (PT,EK)
 * PT = DEC (CT,DK)
 */

/* -----
 * NOTE: currently an implmentation detail limits us to 8192 bit cipher block
 * size.  Changing MAX_BLOCK_SIZE in cipher.c can up that limit without 
 * causing any other compatibility issues. 
 *
 * Two local byte arrays are defined on the stack of 1024 bytes or 8192 bits.
 */
/*
 * NOTE: if in and out overlap incorrectly this will fail
 */

static rc_t KBlockCipherECB (const KBlockCipher * self, const void * _in, void * _out,
                             size_t len, const void * key, bool encrypt)
{
    size_t block_size;
    rc_t rc;

    assert (self);
    assert (_in);
    assert (_out);
    assert (key);
    assert (len > 0);
    assert (encrypt == true || encrypt == false);

    rc = KBlockCipherBlockSize (self, &block_size);
    if (rc)
        ;

    else if (block_size == 0)
        rc = RC (rcKrypto, rcCipher, rcEncoding, rcSize, rcInvalid);

    else if (block_size > MAX_BLOCK_SIZE)
        rc = RC (rcKrypto, rcCipher, rcEncoding, rcSize, rcIncorrect);

    else
    {
        BlockFunc func;
        const char * in;
        char * out;
        size_t full_blocks;
        size_t partial_block_count;
        char temp [MAX_BLOCK_SIZE];

        if (encrypt)
            func = KBlockCipherEncrypt;
        else
            func = KBlockCipherDecrypt;

        in = _in;
        out = _out;

        full_blocks = len / block_size;
        partial_block_count = len % block_size;

        while (full_blocks --)
        {
            rc = func (self, in, out, key);
            in += block_size;
            out += block_size;
        }
        if (partial_block_count)
        {
            /* use the temp to prevent read/write past 
             * claimed size of in and out */
            memmove (temp, in, partial_block_count);
            rc = func (self, temp, temp, key);
            if (rc == 0)
                memmove (out, temp, partial_block_count);
        }
    }
    return rc;
}


KRYPTO_EXTERN
rc_t CC KCipherEncryptECB (KCipher * self, const void * in, void * out, size_t len)
{
    rc_t rc;

    if (self == NULL)
        rc = RC (rcKrypto, rcCipher, rcEncoding, rcSelf, rcNull);

    else if ((in == NULL) || (out == NULL))
        rc = RC (rcKrypto, rcCipher, rcEncoding, rcParam, rcNull);

    if (len == 0)
        return 0;

    return KBlockCipherECB (self->block_cipher, in, out, len, self->encrypt_key, true);
}

KRYPTO_EXTERN rc_t CC KCipherDecryptECB (KCipher * self, const void * in, void * out, size_t len)
{
    rc_t rc;

    if (self == NULL)
        rc = RC (rcKrypto, rcCipher, rcEncoding, rcSelf, rcNull);

    else if ((in == NULL) || (out == NULL))
        rc = RC (rcKrypto, rcCipher, rcEncoding, rcParam, rcNull);

    if (len == 0)
        return 0;

    return KBlockCipherECB (self->block_cipher, in, out, len, self->decrypt_key, false);
}

/* ----------
 * Cipher-Block Chaining
 * CT = (FB = ENC (PT^IV, EK))
 * PT = DEC ((FB = CT), DK)
 */
static void out_block (uint8_t * b)
{
    int ii;
    for (ii = 0; ii < 16 ; ++ii)
        KOutMsg ("%x ", b[ii]);
}


KRYPTO_EXTERN
rc_t CC KCipherEncryptCBC (KCipher * self, const void * _in, void * _out, size_t len)
{
    rc_t rc;

    if (self == NULL)
        rc = RC (rcKrypto, rcCipher, rcEncoding, rcSelf, rcNull);

    else if ((_in == NULL) || (_out == NULL))
        rc = RC (rcKrypto, rcCipher, rcEncoding, rcParam, rcNull);

    else
    {
        size_t block_size;

        rc = KBlockCipherBlockSize (self->block_cipher, &block_size);
        if (rc)
            ;

        else if (block_size == 0)
            rc = RC (rcKrypto, rcCipher, rcEncoding, rcSize, rcInvalid);

        else if (block_size > MAX_BLOCK_SIZE)
            rc = RC (rcKrypto, rcCipher, rcEncoding, rcSize, rcIncorrect);

        else
        {
            const char * in;
            char * out;
            char * ivec;
            size_t full_blocks;
            size_t partial_block_count;
            char temp [MAX_BLOCK_SIZE];

            in = _in;
            out = _out;
            ivec = self->encrypt_ivec;

            full_blocks = len / block_size;
            partial_block_count = len % block_size;

            while (full_blocks --)
            {
                uint32_t idx;

                for (idx = 0; idx < block_size; ++idx)
                    temp[idx] = in[idx] ^ ivec[idx];

                rc = KBlockCipherEncrypt (self->block_cipher, temp, ivec, self->encrypt_key);
                if (rc)
                    return rc;

                memmove (out, ivec, block_size);

                in += block_size;
                out += block_size;
            }
            if (partial_block_count)
            {
                uint32_t idx;

                memset (temp, 0, block_size);

                for (idx = 0; idx < partial_block_count; ++idx)
                    temp[idx] = in[idx] ^ ivec[idx];

                rc = KBlockCipherEncrypt (self->block_cipher, temp, ivec, self->encrypt_key);
                if (rc)
                    return rc;

                memmove (out, ivec, block_size);
            }

        }
    }
    return rc;
}



KRYPTO_EXTERN
rc_t CC KCipherDecryptCBC (KCipher * self, const void * _in, void * _out, size_t len)
{
    rc_t rc;

    if (self == NULL)
        rc = RC (rcKrypto, rcCipher, rcEncoding, rcSelf, rcNull);

    else if ((_in == NULL) || (_out == NULL))
        rc = RC (rcKrypto, rcCipher, rcEncoding, rcParam, rcNull);

    else
    {
        size_t block_size;

        rc = KBlockCipherBlockSize (self->block_cipher, &block_size);
        if (rc)
            ;

        else if (block_size == 0)
            rc = RC (rcKrypto, rcCipher, rcEncoding, rcSize, rcInvalid);

        else if (block_size > MAX_BLOCK_SIZE)
            rc = RC (rcKrypto, rcCipher, rcEncoding, rcSize, rcIncorrect);

        else
        {
            const char * in;
            char * out;
            char * ivec;
            size_t full_blocks;
            size_t partial_block_count;
            char temp [MAX_BLOCK_SIZE];

            in = _in;
            out = _out;
            ivec = self->decrypt_ivec;

            full_blocks = len / block_size;
            partial_block_count = len % block_size;

            while (full_blocks --)
            {
                uint32_t idx;

                rc = KBlockCipherDecrypt (self->block_cipher, in, temp, self->decrypt_key);
                if (rc)
                    return rc;

                for (idx = 0; idx < block_size; ++idx)
                    out[idx] = temp[idx] ^ ivec[idx];

                memmove (ivec, in, block_size);


                in += block_size;
                out += block_size;
            }
            if (partial_block_count)
            {
                uint32_t idx;

                memset (temp, 0, block_size);
                memmove (temp, in, partial_block_count);

                rc = KBlockCipherDecrypt (self->block_cipher, temp, temp, self->decrypt_key);
                if (rc)
                    return rc;

                for (idx = 0; idx < block_size; ++idx)
                    out[idx] = temp[idx] ^ ivec[idx];

                memmove (ivec, in, block_size);
            }

        }
    }
    return rc;
}



/* ----------
 * Cipher Feedback
 * CT = (FB = PT) ^ ENC (IV, EK))
 * PT = (FB = CT) ^ DEC (IV, DK)
 */

#if 0
KRYPTO_EXTERN rc_t CC KCipherEncryptCFB (KCipher * self, const void * in, void * out, size_t len,
                                      uint32_t * bits)
{
    if (self->uses_openssl)
    {
        return KCryptoCipherEncryptCFB (self->block_cipher, in, out, len, bits);
    }
    return RC (rcKrypto, rcCipher, rcEncoding, rcFunction, rcUnsupported);
}


KRYPTO_EXTERN rc_t CC KCipherDecryptCFB (KCipher * self, const void * in, void * out, size_t len, uint32_t * bits)
{
    if (self->uses_openssl)
    {
        return KCryptoCipherDecryptCFB (self->block_cipher, in, out, len, bits);
    }
    return RC (rcKrypto, rcCipher, rcEncoding, rcFunction, rcUnsupported);
}
#endif

/* ----------
 * Propagating cipher-block chaining
 * FB = PT ^ (CT = ENC ((PT^IV), EK))
 * FB = CT ^ (PT = DEC (CT,DK) ^ IV)
 */

/* not implmented yet */

/* ----------
 * Output Feedback
 * CT = PT ^ (FB = ENC (IV, EK))
 * PT = CT ^ (FB = DEC (IV, DK))
 */

#if 0
KRYPTO_EXTERN
rc_t CC KCipherEncryptOFB (KCipher * self, const void * in, void * out, size_t len,
                        uint32_t * bits)
{
    if (self->uses_openssl)
    {
        return KCryptoCipherEncryptOFB (self->block_cipher, in, out, len, bits);
    }
    return RC (rcKrypto, rcCipher, rcEncoding, rcFunction, rcUnsupported);
}

KRYPTO_EXTERN
rc_t CC KCipherDecryptOFB (KCipher * self, const void * in, void * out, size_t len,
                        uint32_t * bits)
{
    if (self->uses_openssl)
    {
        return KCryptoCipherDecryptOFB (self->block_cipher, in, out, len, bits);
    }
    return RC (rcKrypto, rcCipher, rcEncoding, rcFunction, rcUnsupported);
}
#endif

/* Counter
 * IV is a nonce and not re-used as FB
 * CT = PT ^ ENC (N, EK)
 * PT = CT ^ ENC (N, DK) 
 * Note decrypt is encrypt.
 * nonce is a function that given an iv generates the next iv
 */
#if 0
KRYPTO_EXTERN rc_t CC KCipherEncryptCTR (KCipher * self, const void * in, void * out, 
                                   void(*nonce)(void*iv), void * iv)
{
    return RC (rcKrypto, rcCipher, rcEncoding, rcFunction, rcUnsupported);
}
#endif

rc_t CC KCipherAlloc (KCipher ** obj, const KBlockCipher * block_cipher)
{
    void * mem;
    size_t bsize;
    size_t ksize;
    rc_t rc;

    if ((obj == NULL) || (block_cipher == NULL))
        rc = RC (rcKrypto, rcCipher, rcConstructing, rcParam, rcNull);
    else
    {
        rc = KBlockCipherBlockSize (block_cipher, &bsize);
        if (rc == 0)
        {
            rc = KBlockCipherKeySize (block_cipher, &ksize);
            if (rc == 0)
            {
                mem = malloc (sizeof (**obj) + 2 * bsize + 2 * ksize);
                if (mem)
                {
                    *obj = mem;
                    return 0;
                }
                rc =  RC (rcKrypto, rcCipher, rcConstructing, rcMemory,
                          rcExhausted);
            }
        }
    }
    return rc;
}


rc_t CC KCipherInit (KCipher * self, const KBlockCipher * block_cipher)
{
    rc_t rc;

    if (self == NULL)
        rc = RC (rcKrypto, rcCipher, rcConstructing, rcSelf, rcNull);

    else if (block_cipher == NULL)
        rc = RC (rcKrypto, rcCipher, rcConstructing, rcParam, rcNull);

    else
    {
        size_t bsize;

        rc = KBlockCipherBlockSize (block_cipher, &bsize);
        if (rc == 0)
        {
            size_t ksize;

            rc = KBlockCipherKeySize (block_cipher, &ksize);
            if (rc == 0)
            {
                rc = KBlockCipherAddRef (block_cipher);
                if (rc == 0)
                {
                    KRefcountInit (&self->refcount, 1, KCipherClassName,
                                   "init", block_cipher->name);

                    self->block_cipher = block_cipher;

                    StringInitCString (&self->name, block_cipher->name);

                    self->encrypt_key = (char*)(self + 1);
                    self->decrypt_key = ((char*)self->encrypt_key) + ksize;
                    self->encrypt_ivec = ((char*)self->decrypt_key) + ksize;
                    self->decrypt_ivec = ((char*)self->encrypt_ivec) + bsize;
                }
            }
        }
    }
    return rc;
}



