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

#include <vdb/extern.h>
 
#define TRACK_REFERENCES 0

#include "page-map.h"
#include "blob-headers.h"
#include "blob.h"
#include "blob-priv.h"
#include <klib/rc.h>
#include <klib/defs.h>
#include <byteswap.h>
#include <klib/data-buffer.h>
#include <klib/vlen-encode.h>
#include <vdb/schema.h>
#include <vdb/xform.h>
#include <klib/log.h>
#include <sysalloc.h>
#include <bitstr.h>

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <stdio.h> /* temp for debugging */

#if _DEBUGGING
void VBlobCheckIntegrity ( const VBlob *self )
{
    if ( self != NULL )
    {
        rc_t rc = KDataBufferCheckIntegrity ( & self -> data );
        if ( rc != 0 )
        {
            fprintf ( stderr, "AAAAAH!\n" );
        }
    }
}
#endif

rc_t VBlobNew ( VBlob **lhs, int64_t start_id, int64_t stop_id, const char *name ) {
    VBlob *y;
    
    if ( name == NULL )
        name = "";
#if VBLOG_HAS_NAME
    *lhs = y = malloc(sizeof(*y) + strlen(name));
#else
    *lhs = y = calloc(1, sizeof(*y));
#endif
    if (y) {
        KRefcountInit(&y->refcount, 1, "VBlob", "new", name);
        y->start_id = start_id;
        y->stop_id = stop_id;
        y->data.elem_bits = 1;
        y->byte_order = vboNative;
#if VBLOG_HAS_NAME
        y->pm = NULL;
        y->headers = NULL;
        y->spmc = NULL;
        memset(&y->data, 0, sizeof(y->data));
        y->no_cache = 0;
        strcpy(&(((char *)y->name)[0]), name);
#endif
        
        return 0;
    }
    return RC(rcVDB, rcBlob, rcConstructing, rcMemory, rcExhausted);
}

static rc_t VBlobDestroy( VBlob *that ) {
    if (that->spmc) {
        int i;
        
        for (i = 0; i != that->spmc->n; ++i)
            PageMapRelease(that->spmc->pm[i]);
        free(that->spmc);
    }
    KDataBufferWhack(&that->data);
    BlobHeadersRelease(that->headers);
    PageMapRelease(that->pm);
    free(that);
    return 0;
}

#ifdef VBlobAddRef
rc_t VBlobReleaseLast ( VBlob *self )
#else
rc_t VBlobAddRef ( const VBlob *self )
{
    if ( self != NULL )
        KRefcountAdd(&self->refcount, "VBlob");
    return 0;
}

rc_t VBlobRelease ( const VBlob *self )
#endif
{
    rc_t rc = 0;

    if ( self != 0 )
    {
        switch ( KRefcountDrop(&self->refcount, "VBlob") )
        {
        case krefWhack:
            return VBlobDestroy ( self );
        case krefNegative:
            rc = RC (rcVDB, rcBlob, rcDestroying, rcBlob, rcExcessive);
            PLOGERR ( klogInt, (klogInt, rc, "Released a blob $(B) with no more references",
                      PLOG_P(self)));
            break;
        }
    }
    return rc;
}

static rc_t decode_header_byte_v2(
                                  uint8_t header_byte,
                                  uint8_t *variant,
                                  uint8_t *adjust,
                                  VByteOrder *byte_order
) {
    *adjust = (8 - (header_byte & 7)) & 7; header_byte >>= 3;
    *byte_order = (header_byte & 1) ? vboBigEndian : vboLittleEndian;
    header_byte >>= 1;
    *variant = header_byte & 3;
    header_byte >>= 2;
    return header_byte == 2 ? 0 : RC(rcVDB, rcBlob, rcReading, rcData, rcBadVersion);
}

static rc_t decode_header_v2_0(
                               const uint8_t *src,
                               uint64_t ssize,
                               uint32_t *hdr_size,
                               uint32_t *map_size,
                               uint32_t *offset
) {
    *offset = 3;
    if (ssize < *offset)
        return RC(rcVDB, rcBlob, rcConstructing, rcData, rcInsufficient);
    
    *hdr_size = src[1];
    *map_size = src[2];
    return 0;
}

static rc_t decode_header_v2_1(
                               const uint8_t *src,
                               uint64_t ssize,
                               uint32_t *hdr_size,
                               uint32_t *map_size,
                               uint32_t *offset
) {
    *offset = 4;
    if (ssize < *offset)
        return RC(rcVDB, rcBlob, rcConstructing, rcData, rcInsufficient);
    
    *hdr_size = src[1];
    *map_size = (uint32_t)src[2] | ((uint32_t)src[3] << 8);
    return 0;
}

static rc_t decode_header_v2_2(
                               const uint8_t *src,
                               uint64_t ssize,
                               uint32_t *hdr_size,
                               uint32_t *map_size,
                               uint32_t *offset
) {
    *offset = 6;
    if (ssize < *offset)
        return RC(rcVDB, rcBlob, rcConstructing, rcData, rcInsufficient);
    
    *hdr_size = src[1];
    *map_size = (uint32_t)src[2] | ((uint32_t)src[3] << 8) | ((uint32_t)src[4] << 16) | ((uint32_t)src[5] << 24);
    return 0;
}

static rc_t decode_header_v2_3(
                               const uint8_t *src,
                               uint64_t ssize,
                               uint32_t *hdr_size,
                               uint32_t *map_size,
                               uint32_t *offset
) {
    *offset = 9;
    if (ssize < *offset)
        return RC(rcVDB, rcBlob, rcConstructing, rcData, rcInsufficient);
    
    *hdr_size = (uint32_t)src[1] | ((uint32_t)src[2] << 8) | ((uint32_t)src[3] << 16) | ((uint32_t)src[4] << 24);
    *map_size = (uint32_t)src[5] | ((uint32_t)src[6] << 8) | ((uint32_t)src[7] << 16) | ((uint32_t)src[8] << 24);
    return 0;
}

static rc_t decode_header_v2(
                             const uint8_t *src,
                             uint64_t ssize,
                             uint32_t *hdr_size,
                             uint32_t *map_size,
                             uint32_t *offset,
                             uint8_t *adjust,
                             VByteOrder *byte_order
) {
    rc_t rc;
    uint8_t variant;
    
    if (ssize == 0)
        return RC(rcVDB, rcBlob, rcConstructing, rcData, rcInsufficient);

    rc = decode_header_byte_v2(src[0], &variant, adjust, byte_order);
    if (rc)
        return rc;

    switch (variant) {
    case 0:
        return decode_header_v2_0(src, ssize, hdr_size, map_size, offset);
    case 1:
        return decode_header_v2_1(src, ssize, hdr_size, map_size, offset);
    case 2:
        return decode_header_v2_2(src, ssize, hdr_size, map_size, offset);
    case 3:
        return decode_header_v2_3(src, ssize, hdr_size, map_size, offset);
    default:
        return RC(rcVDB, rcBlob, rcConstructing, rcData, rcBadVersion);
    }
}

static rc_t encode_header_v1(
                             uint8_t *dst,
                             uint64_t dsize,
                             uint64_t *used,
                             uint32_t row_length,
                             bitsz_t data_size,
                             VByteOrder byte_order
) {
    /* byte-order goes in bits 0..1 */
    uint8_t header_byte = byte_order & 3;
    if ( header_byte == vboNative )
    {
#if __BYTE_ORDER == __LITTLE_ENDIAN
        header_byte = ( uint8_t) vboLittleEndian;
#else
        header_byte = ( uint8_t) vboBigEndian;
#endif
    }

    /* blob size adjust goes in bits 2..4 */
    header_byte |= ( ( 8 - ( data_size & 7 ) ) & 7 ) << 2;
    
    /* row-length code goes in bits 5..6 */
    if ( row_length == 1 ) {
        header_byte |= 3 << 5;
        * used = 1;
        if ( dsize < * used )
            return RC(rcVDB, rcBlob, rcConstructing, rcBuffer, rcInsufficient);
        dst[0] = header_byte;
    }
    else if (row_length < 0x100) {
        *used = 2;
        if (dsize < *used)
            return RC(rcVDB, rcBlob, rcConstructing, rcBuffer, rcInsufficient);
        dst[0] = header_byte;
        dst[1] = ( uint8_t ) row_length;
    }
    else if (row_length < 0x10000) {
        header_byte |= 1 << 5;
        *used = 3;
        if (dsize < *used)
            return RC(rcVDB, rcBlob, rcConstructing, rcBuffer, rcInsufficient);
        dst[0] = header_byte;
        dst[1] = ( uint8_t ) row_length;
        dst[2] = ( uint8_t ) ( row_length >> 8 );
    }
    else {
        header_byte |= 2 << 5;
        *used = 5;
        if (dsize < *used)
            return RC(rcVDB, rcBlob, rcConstructing, rcBuffer, rcInsufficient);
        dst[0] = header_byte;
        dst[1] = ( uint8_t ) row_length;
        dst[2] = ( uint8_t ) ( row_length >> 8 );
        dst[3] = ( uint8_t ) ( row_length >> 16 );
        dst[4] = ( uint8_t ) ( row_length >> 24 );
    }
    return 0;
}

static rc_t encode_header_v2(
                             uint8_t *dst,
                             uint64_t dsize,
                             uint64_t *used,
                             uint64_t hdr_size,
                             uint64_t map_size,
                             bitsz_t data_size
) {
#if __BYTE_ORDER == __LITTLE_ENDIAN
    uint8_t header_byte = 0x80 | ( (uint8_t)data_size & 7 );
#else
    uint8_t header_byte = 0x88 | ( (uint8_t)data_size & 7 );
#endif
    
    assert(hdr_size >> 32 == 0);
    assert(map_size >> 32 == 0);
    
    if ((hdr_size >> 8) == 0) {
        if ((map_size >> 8) == 0) {
            *used = 3;
            if (dsize < *used)
                return RC(rcVDB, rcBlob, rcConstructing, rcBuffer, rcInsufficient);
            
            dst[0] = header_byte;
            dst[1] = hdr_size;
            dst[2] = map_size;
        }
        else if ((map_size >> 16) == 0) {
            *used = 4;
            if (dsize < *used)
                return RC(rcVDB, rcBlob, rcConstructing, rcBuffer, rcInsufficient);
            
            dst[0] = header_byte | 0x10;
            dst[1] = hdr_size;
            dst[2] = map_size;
            dst[3] = map_size >> 8;
        }
        else {
            *used = 6;
            if (dsize < *used)
                return RC(rcVDB, rcBlob, rcConstructing, rcBuffer, rcInsufficient);

            dst[0] = header_byte | 0x20;
            dst[1] = hdr_size;
            dst[2] = map_size;
            dst[3] = map_size >> 8;
            dst[4] = map_size >> 16;
            dst[5] = map_size >> 24;
        }
    }
    else {
        *used = 9;
        if (dsize < *used)
            return RC(rcVDB, rcBlob, rcConstructing, rcBuffer, rcInsufficient);

        dst[0] = header_byte | 0x30;

        dst[1] = hdr_size;
        dst[2] = hdr_size >> 8;
        dst[3] = hdr_size >> 16;
        dst[4] = hdr_size >> 24;

        dst[5] = map_size;
        dst[6] = map_size >> 8;
        dst[7] = map_size >> 16;
        dst[8] = map_size >> 24;
    }
    return 0;
}

static
rc_t VBlobCreateFromData_v2(
                            VBlob **lhs,
                            const KDataBuffer *data,
                            int64_t start_id, int64_t stop_id,
                            uint32_t elem_bits
) {
    const uint8_t *src = data->base;
    uint64_t ssize = data->elem_count;
    uint32_t hsize;
    uint32_t msize;
    uint32_t offset;
    VByteOrder byte_order;
    uint8_t adjust;
    VBlob *y;
    uint32_t data_offset;
    bitsz_t databits;
    uint32_t elem_count;
    rc_t rc;
    
    rc = decode_header_v2(src, ssize, &hsize, &msize, &offset, &adjust, &byte_order);
    if (rc)
        return rc;

    if (ssize < offset + hsize + msize)
        return RC(rcVDB, rcBlob, rcConstructing, rcData, rcInsufficient);
    
    src += offset;
    data_offset = offset + hsize + msize;
    assert(data_offset <= ssize);
    ssize -= data_offset;
    databits = (ssize << 3) - adjust;
    assert(databits % elem_bits == 0);
    elem_count = (uint32_t)( databits / elem_bits );

    rc = VBlobNew(&y, start_id, stop_id, NULL);
    TRACK_BLOB (VBlobNew, y);
    if (rc == 0) {
        if (hsize)
            rc = BlobHeadersCreateFromData(&y->headers, src, hsize);
        if (rc == 0) {
            rc = PageMapDeserialize(&y->pm, src + hsize, msize, BlobRowCount(y));
            if (rc == 0) {
                KDataBufferSub(data, &y->data, data_offset, ssize);
                y->data.elem_bits = elem_bits;
                y->data.elem_count = elem_count;
                y->byte_order = byte_order;
                *lhs = y;
                return 0;
            }
        }
        (void)VBlobRelease(y);
        TRACK_BLOB (VBlobRelease, y);
    }
    return rc;
}

static
rc_t VBlobCreateFromData_v1(
                            VBlob **lhs,
                            const KDataBuffer *data,
                            int64_t start_id, int64_t stop_id,
                            uint32_t elem_bits
) {
    const uint8_t *src = data->base;
    uint64_t ssize = data->elem_count;
    uint8_t header;
    rc_t rc;
    VBlob *y;
    VByteOrder byte_order;
    uint32_t offset;
    int adjust;
    int rls; /* row length size */
    uint64_t row_len;
    bitsz_t databits;
    
    if (ssize == 0)
        return RC(rcVDB, rcBlob, rcConstructing, rcData, rcInsufficient);
    
    header = *src;
    byte_order = header & 3; header >>= 2;
    adjust = header & 7;
    header >>= 3;
    rls = header & 3;
    
    /* convert rls from a code to an actual length */
    rls = "\x01\x02\x04\x00" [ rls ];

    /* adjust offset */
    offset = rls + 1;

    /* handle special code where row length is implicitly 1 */
    if ( rls == 0 )
        row_len = 1;

    /* ensure sufficient header bytes */
    else if ( ssize < offset )
        return RC(rcVDB, rcBlob, rcConstructing, rcData, rcInsufficient);
    else
    {
        /* produce little-endian 64-bit row-length */
        row_len = 0;
        memcpy ( & row_len, & src [ 1 ], rls );

#if __BYTE_ORDER != __LITTLE_ENDIAN
        /* correct for big-endian */
        row_len = bswap_64 ( row_len );
#endif
    }

    ssize -= offset;
    databits = (ssize << 3) - adjust;
    assert(databits % elem_bits == 0);

    rc = VBlobNew(&y, start_id, stop_id, NULL);
    TRACK_BLOB (VBLobNew, y);
    if (rc == 0) {

        uint64_t row_count = BlobRowCount ( y );

        /* test for badly formed row-length */
        if ( rls == 4 )
        {
            assert ( row_len != 0 );
            if ( row_len * row_count != databits / elem_bits )
            {
                /* we can fix a length if we know the count */
                if ( row_count != 0 )
                    row_len = ( databits / elem_bits ) / row_count;
                else
                {
                    /* rely on code to handle legacy blobs in prod-cmn.c:VFunctionProdCallByteswap */
                    row_len = 0;
                }
            }
        }

        if ( row_len != 0 )
            rc = PageMapNewFixedRowLength( &y->pm, row_count, row_len );
        if (rc == 0) {
            KDataBufferSub(data, &y->data, offset, ssize);
            y->data.elem_bits = elem_bits;
            y->data.elem_count = (uint32_t)( databits / elem_bits );
            y->byte_order = byte_order;
			
            *lhs = y;
            return 0;
        }
        /* like a call to VBlobRelease (y); */
        TRACK_BLOB (VBlobRelease-free, y);
        free(y);
    }
    return rc;
}

rc_t VBlobSerialize ( const VBlob *self, KDataBuffer *result ) {
    uint64_t sz;
    rc_t rc;
    bitsz_t data_bits = KDataBufferBits(&self->data);
    uint64_t data_bytes = KDataBufferBytes(&self->data);
    uint32_t row_length;
    
    if (self->headers == NULL && (row_length = PageMapHasSimpleStructure(self->pm)) != 0) {
        rc = KDataBufferResize(result, 5 + data_bytes);
        if (rc == 0) {

#if _DEBUGGING && 1
            /* temporary assert that we are setting byte_order properly
               in the future, we may allow some functions to issue other
               byte orders, although there is no conceivable reason to do so */
#if __BYTE_ORDER == __LITTLE_ENDIAN
            assert ( self -> byte_order == vboNative || self -> byte_order == vboLittleEndian );
#else
            assert ( self -> byte_order == vboNative || self -> byte_order == vboBigEndian );
#endif
#endif
            rc = encode_header_v1(result->base, result->elem_count, &sz, row_length, data_bits, self->byte_order);
            if (rc == 0) {
                memcpy(&((uint8_t *)result->base)[sz], self->data.base, data_bytes);
                result->elem_count = sz + data_bytes;
            }
        }
    }
    else {
        KDataBuffer headers;
        KDataBuffer pagemap;
        
        rc = KDataBufferMakeBytes(&headers, 0);
        if (rc == 0) {
            if (self->headers)
                rc = BlobHeadersSerialize(self->headers, &headers, 0, &sz);
            else
                sz = 0;
            if (rc == 0) {
                headers.elem_count = sz;
                rc = KDataBufferMakeBytes(&pagemap, 0);
                if (rc == 0) {
                    if (self->pm)
                        rc = PageMapSerialize(self->pm, &pagemap, 0, &sz);
                    else
                        sz = 0;
                    if (rc == 0) {
                        pagemap.elem_count = sz;
                        rc = KDataBufferResize(result, 9 + data_bytes + headers.elem_count + pagemap.elem_count);
                        if (rc == 0) {
                            rc = encode_header_v2(result->base, result->elem_count, &sz, headers.elem_count, pagemap.elem_count, data_bits);
                            if (rc == 0) {
                                memcpy(&((uint8_t *)result->base)[sz], headers.base, headers.elem_count);
                                sz += headers.elem_count;
                                memcpy(&((uint8_t *)result->base)[sz], pagemap.base, pagemap.elem_count);
                                sz += pagemap.elem_count;
                                memcpy(&((uint8_t *)result->base)[sz], self->data.base, data_bytes);
                                result->elem_count = sz + data_bytes;
                            }
                        }
                    }
                    KDataBufferWhack(&pagemap);
                }
            }
        }
        KDataBufferWhack(&headers);
    }
    
    return rc;
}

rc_t VBlobCreateFromData ( struct VBlob **lhs,
                         int64_t start_id, int64_t stop_id,
                         const KDataBuffer *src,
                         uint32_t elem_bits )
{
    VBlob *y = NULL;
    rc_t rc;
    
    assert(lhs);
    assert(src);
    assert(src->elem_bits == 8);
    assert(src->bit_offset == 0);

    *lhs = 0;

    if ((((const uint8_t *)src->base)[0] & 0x80) == 0)
        rc = VBlobCreateFromData_v1(&y, src, start_id, stop_id, elem_bits);
    else
        rc = VBlobCreateFromData_v2(&y, src, start_id, stop_id, elem_bits);

    if (rc == 0)
        *lhs = y;

    return rc;
}

rc_t VBlobCreateFromSingleRow (
			      struct VBlob **lhs,
			      int64_t start_id, int64_t stop_id,
			      const KDataBuffer *src,
			      VByteOrder byte_order )
{
    VBlob *y;
    rc_t rc;
    
    rc = VBlobNew(&y, start_id, stop_id, NULL);
    TRACK_BLOB (VBlobNew, y);
    if (rc == 0) {
        assert(src->elem_count >> 32 == 0);
        rc = PageMapNewSingle(&y->pm, BlobRowCount(y), (uint32_t)src->elem_count);
        if (rc == 0) {
            rc = KDataBufferSub(src, &y->data, 0, UINT64_MAX);
            if (rc == 0) {
                y->byte_order = byte_order;
                *lhs = y;
                return 0;
            }
        }
        /* should add a release/free? */
    }
    return rc;
}

bool VBlobIsSingleRow( const struct VBlob *self ) {
    return self->pm && PageMapFastRowCount(self->pm) == BlobRowCount(self) ? true : false;
}

uint32_t VBlobFixedRowLength( const struct VBlob *self ) {
    return self->pm ? PageMapFixedRowLength(self->pm) : 0;
}

#define COMPARE(FORCE, BITS, DBASE, DOFF, SBASE, SOFF, LENGTH) \
    (((FORCE == 0) || (BITS & 7) != 0) ? ( \
        bitcmp(DBASE, DOFF * BITS + FORCE, SBASE, SOFF * BITS, LENGTH * BITS)) : ( \
        memcmp(((const char *)DBASE)+((DOFF * BITS) >> 3), \
               ((const char *)SBASE)+((SOFF * BITS) >> 3), \
               ((LENGTH * BITS) >> 3))))

#define COPY(FORCE, BITS, DBASE, DOFF, SBASE, SOFF, LENGTH) \
    (((FORCE == 0) || (BITS & 7) != 0) ? ( \
        bitcpy(DBASE, DOFF * BITS + FORCE, SBASE, SOFF * BITS, LENGTH * BITS)) : ( \
        (void)memcpy(((      char *)DBASE)+((DOFF * BITS) >> 3), \
               ((const char *)SBASE)+((SOFF * BITS) >> 3), \
               ((LENGTH * BITS) >> 3))))

rc_t VBlobAppendRow(VBlob *self,
                    elem_count_t *last_offset,
                    elem_count_t *last_length,
                    const KDataBuffer *src,
                    elem_count_t offset,
                    elem_count_t length,
                    row_count_t repeat_count
                    )
{
    rc_t rc;
    
    if (!PageMapHasRows(self->pm) || length != *last_length ||
        COMPARE(self->data.bit_offset, self->data.elem_bits,
                self->data.base, *last_offset,
                src->base, offset,
                length) != 0
        )
    {
        *last_offset = self->data.elem_count;
        rc = KDataBufferResize(&self->data, *last_offset + length);
        if (rc == 0) {
            COPY(self->data.bit_offset, self->data.elem_bits,
                 self->data.base, *last_offset,
                 src->base, offset, length);
            rc = PageMapAppendRows(self->pm, length, repeat_count, false);
        }
        *last_length = length;
    }
    else
        rc = PageMapAppendRows(self->pm, length, repeat_count, true);
    
    return rc;
}

static rc_t VBlobGetLastRow(VBlob *self, elem_count_t *offset, elem_count_t *length) {
    
    *length = PageMapLastLength(self->pm);
    *offset = self->data.elem_count - *length;

    return 0;
}

rc_t VBlobAppend(VBlob *self, const VBlob *other) {
    rc_t rc;
    row_count_t offset;
    row_count_t length;
    
    if (self->headers)
        return RC(rcVDB, rcBlob, rcConcatenating, rcSelf, rcInconsistent);
    if (other->headers)
        return RC(rcVDB, rcBlob, rcConcatenating, rcParam, rcInvalid);

    if (self->stop_id + 1 != other->start_id)
        return RC(rcVDB, rcBlob, rcConcatenating, rcId, rcOutofrange);

    if (other->data.elem_bits != self->data.elem_bits)
        return RC(rcVDB, rcBlob, rcConcatenating, rcData, rcInvalid);

#if 0
    fprintf(stderr, "appending %u(%u) (length: %u) to %u(%u) (length: %u) %s\n",
            (unsigned)other->start_id, (unsigned)BlobRowCount(other),
            (unsigned)other->data.elem_count,
            (unsigned)self->start_id, (unsigned)BlobRowCount(self),
            (unsigned)self->data.elem_count,
            self->name);
#endif
    
    rc = VBlobGetLastRow(self, &offset, &length);
    if (rc == 0) {
        PageMapIterator iter;
        
        rc = PageMapNewIterator(other->pm, &iter, 0, -1);
        if (rc == 0) {
            KDataBuffer orig;
            
            rc = KDataBufferMakeWritable(&self->data , &orig);
            if (rc == 0) {
                row_count_t row_count;
                KDataBufferWhack(&self->data);
                self->data = orig;
                do {
                    row_count = PageMapIteratorRepeatCount(&iter);
                    rc = VBlobAppendRow(self, &offset, &length, &other->data,
                                        PageMapIteratorDataOffset(&iter),
                                        PageMapIteratorDataLength(&iter),
                                        row_count);
                } while (rc == 0 && PageMapIteratorAdvance(&iter, row_count));
                if (rc == 0) {
                    self->stop_id = other->stop_id;
                    self->no_cache |= other->no_cache;
                }
            }
        }
    }
    return rc;
}

rc_t VBlobSubblob( const struct VBlob *self,struct VBlob **sub, int64_t start_id )
{
    rc_t rc;
    KDataBuffer  kd;
    PageMapIterator pmi;
    
    if (start_id < self->start_id || start_id > self->stop_id)
        return RC(rcVDB, rcBlob, rcConverting, rcId, rcOutofrange);
    
    rc=PageMapNewIterator(self->pm,&pmi, 0, -1);
    if(rc == 0){
        if(PageMapIteratorAdvance(&pmi,start_id-self->start_id)){
            row_count_t numrep = PageMapIteratorRepeatCount(&pmi);
            elem_count_t offset = PageMapIteratorDataOffset(&pmi);
            elem_count_t length = PageMapIteratorDataLength(&pmi);
            
#if 0
            fprintf(stderr, "splitting %u(%u) (offset: %u, length: %u) from %s\n",
                    (unsigned)start_id, (unsigned)numrep,
                    (unsigned)offset, (unsigned)length,
                    self->name);
#endif
            
            rc = KDataBufferSub(&self->data, &kd, offset, length);
            if(rc == 0){
		int64_t	stop_id;

		if(length > 0) stop_id = start_id + numrep - 1;
		else           stop_id = start_id; /*** HACK - 0 sized data may be a sign that real data is somewhere else ***/

                rc = VBlobCreateFromSingleRow(sub, start_id, stop_id, &kd, self->byte_order);
                KDataBufferWhack(&kd);
            }
        } else {
            rc = RC(rcVDB, rcBlob, rcConverting, rcId, rcOutofrange);
        }
    }
    return rc;
}

