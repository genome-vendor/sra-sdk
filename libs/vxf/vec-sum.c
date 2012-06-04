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
#include <klib/defs.h>
#include <klib/rc.h>
#include <vdb/xform.h>
#include <vdb/schema.h>
#include <sysalloc.h>

#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

typedef void (*math_f)(void *dst, const void *src, uint32_t vec_length, uint32_t vec_count);
typedef struct self_t {
    uint32_t vec_length;
    math_f f;
} self_t;

#define FUNC(VALTYPE) F_ ## VALTYPE

#define FUNC_DEF(VALTYPE) \
static void FUNC(VALTYPE)(void *Dst, const void *Src, uint32_t vec_length, uint32_t vec_count) { \
    VALTYPE sum; \
    VALTYPE *dst = (VALTYPE *)Dst; \
    const VALTYPE *src = (const VALTYPE *)Src; \
    uint32_t i; \
    uint32_t j; \
    uint32_t k; \
    \
    for (i = k = 0; i != vec_count; ++i) { \
        for (sum = 0, j = 0; j != vec_length; ++j, ++k) \
            sum += src[k]; \
        dst[i] = sum; \
    } \
}

FUNC_DEF(float)
FUNC_DEF(double)
FUNC_DEF(uint8_t)
FUNC_DEF(uint16_t)
FUNC_DEF(uint32_t)
FUNC_DEF(uint64_t)
FUNC_DEF(int8_t)
FUNC_DEF(int16_t)
FUNC_DEF(int32_t)
FUNC_DEF(int64_t)

static
rc_t CC array_func(
                void *Self,
                const VXformInfo *info,
                void *dst,
                const void *src,
                uint64_t elem_count
) {
    const self_t *self = Self;
    
    assert(elem_count % self->vec_length == 0);
    assert((elem_count / self->vec_length) >> 32 == 0);
    self->f(dst, src, self->vec_length, (uint32_t)(elem_count / self->vec_length));
    return 0;
}

static
void CC vxf_vec_sum_wrapper( void *ptr )
{
	free( ptr );
}

/*
 function < type T, U32 dim >
 T vec_sum #1.0 ( T [ dim ] in )
 */
VTRANSFACT_IMPL(vdb_vec_sum, 1, 0, 0) (
                                       const void *Self,
                                       const VXfactInfo *info,
                                       VFuncDesc *rslt,
                                       const VFactoryParams *cp,
                                       const VFunctionParams *dp
) {
    self_t *self;
    rc_t rc = 0;
    
    self = malloc(sizeof(*self));
    if (self == NULL)
        return RC(rcVDB, rcFunction, rcConstructing, rcMemory, rcExhausted);
    
    rslt->self = self;
    rslt->whack = vxf_vec_sum_wrapper;
    rslt->variant = vftArray;
    rslt->u.af = array_func;
    
    self->vec_length = dp->argv[0].fd.td.dim;
    
    switch (info->fdesc.desc.intrinsic_bits) {
    case 8:
        switch (info->fdesc.desc.domain) {
        case vtdInt:
            self->f = FUNC(int8_t);
            break;
        case vtdUint:
            self->f = FUNC(uint8_t);
            break;
        default:
            rc = RC(rcVDB, rcFunction, rcConstructing, rcParam, rcInvalid);
        }
        break;
    case 16:
        switch (info->fdesc.desc.domain) {
        case vtdInt:
            self->f = FUNC(int16_t);
            break;
        case vtdUint:
            self->f = FUNC(uint16_t);
            break;
        default:
            rc = RC(rcVDB, rcFunction, rcConstructing, rcParam, rcInvalid);
        }
        break;
    case 32:
        switch (info->fdesc.desc.domain) {
        case vtdInt:
            self->f = FUNC(int32_t);
            break;
        case vtdUint:
            self->f = FUNC(uint32_t);
            break;
        case vtdFloat:
            self->f = FUNC(float);
            break;
        default:
            rc = RC(rcVDB, rcFunction, rcConstructing, rcParam, rcInvalid);
        }
        break;
    case 64:
        switch (info->fdesc.desc.domain) {
        case vtdInt:
            self->f = FUNC(int64_t);
            break;
        case vtdUint:
            self->f = FUNC(uint64_t);
            break;
        case vtdFloat:
            self->f = FUNC(double);
            break;
        default:
            rc = RC(rcVDB, rcFunction, rcConstructing, rcParam, rcInvalid);
        }
        break;
    default:
        rc = RC(rcVDB, rcFunction, rcConstructing, rcParam, rcInvalid);
    }
    if (rc)
        free(self);
    return 0;
}
