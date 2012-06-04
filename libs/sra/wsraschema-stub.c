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

#include <sra/sraschema.h>
#include <sra/sradb-priv.h>
#include <kfs/dyload.h>
#include <klib/rc.h>

#include <sysalloc.h>
#include <os-native.h>

#include "sra-priv.h"

#include <string.h>
#include <assert.h>


/*--------------------------------------------------------------------------
 * SRASchema
 *  a schema object pre-loaded with default SRA schema
 */

#if ! _STATIC
/* LoadLibrary
 *  loads libsraschema and fills out function pointers
 */
static KDylib *libsraschema = NULL;

static struct
{
    rc_t ( CC * sraSchemaMake ) ( struct VSchema **schema, struct VDBManager const *mgr );

} imports;

static
rc_t SRASchemaLoadLibrary ( void )
{
    KDyld *dl;
    rc_t rc = KDyldMake ( & dl );
    if ( rc == 0 )
    {
        rc = KDyldLoadLib ( dl, & libsraschema, LPFX "wsra-schema" SHLX );
        if ( rc == 0 )
        {
            KSymAddr *sym;

            /* resolve symbols */
            rc = KDylibSymbol ( libsraschema, & sym, "SRASchemaMake" );
            if ( rc == 0 )
            {
                KSymAddrAsFunc ( sym, ( fptr_t* ) & imports . sraSchemaMake );
                KSymAddrRelease ( sym );
            }

            /* bail on error */
            if ( rc != 0 )
            {
                KDylibRelease ( libsraschema );
                libsraschema = NULL;
                memset ( & imports, 0, sizeof imports );
            }
        }
        KDyldRelease ( dl );
    }

    return rc;
}
#endif

/* Make
 *  create an instance of the default SRA schema
 */
rc_t CC VDBManagerMakeSRASchema ( struct VDBManager const *self, struct VSchema **schema )
{
#if _STATIC
    return SRASchemaMake ( schema, self );
#else
    if ( libsraschema == NULL )
    {
        rc_t rc = SRASchemaLoadLibrary ();
        if ( rc != 0 )
            return rc;
    }

    assert ( imports . sraSchemaMake != NULL );
    return ( * imports . sraSchemaMake ) ( schema, self );
#endif
}

rc_t CC SRAMgrMakeSRASchema ( const SRAMgr *self, struct VSchema **schema )
{
    if ( self != NULL )
        return VDBManagerMakeSRASchema ( self -> vmgr, schema );
    return RC ( rcSRA, rcMgr, rcAccessing, rcSelf, rcNull );
}
