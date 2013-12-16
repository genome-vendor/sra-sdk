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

#include "keyring.h"

#include <klib/rc.h>
#include <klib/refcount.h>

#include <kns/stream.h>

#include "keyring-priv.h"

#include <sysalloc.h>
#include <stdlib.h>
#include <string.h>
 
#define rcTarget rcNoTarg

struct KKeyRing
{
    KRefcount refcount;

    bool read_only;
    
    KStream* ipc;
};

static
rc_t CC KKeyRingInit(KKeyRing* self)
{
    memset(self, 0, sizeof(KKeyRing));

    KRefcountInit ( & self -> refcount, 1, "KKeyRing", "init", "" );
    
    return StartKeyRing(&self->ipc);
}

static
rc_t CC KKeyRingWhack(KKeyRing* self)
{
    rc_t rc = KStreamRelease(self->ipc);
    free(self);
    return rc;
}

extern
rc_t CC KKeyRingAddRef ( const KKeyRing *self )
{
    if ( self != NULL )
    {
        switch ( KRefcountAdd ( & self -> refcount, "KKeyRing" ) )
        {
        case krefLimit:
            return RC ( rcKFG, rcTarget, rcAttaching, rcRange, rcExcessive );
        }
    }
    return 0;
}

extern 
rc_t CC KKeyRingRelease ( KKeyRing *self, bool shutdown_server )
{
    if ( self != NULL )
    {
        switch ( KRefcountDrop ( & self -> refcount, "KKeyRing" ) )
        {
        case krefWhack:
            return KKeyRingWhack ( ( KKeyRing* ) self );
        break;
        case krefLimit:
            return RC ( rcKFG, rcTarget, rcReleasing, rcRange, rcExcessive );
        }
    }
    return 0;
}

LIB_EXPORT rc_t CC KKeyRingMakeRead( const KKeyRing** cself )
{
    KKeyRing** self = (KKeyRing**)cself;
    rc_t rc = KKeyRingMakeUpdate(self);
    if (rc == 0)
        (*self)->read_only = true;
    return rc;
}

LIB_EXPORT 
rc_t CC KKeyRingMakeUpdate(KKeyRing** self)
{
    KKeyRing* obj;
    rc_t rc;
    
    if ( self == NULL )
        rc = RC ( rcKFG, rcTarget, rcCreating, rcParam, rcNull );
    else
    {
        obj = malloc(sizeof(KKeyRing));
        if (obj == NULL)
            rc = RC ( rcKFG, rcTarget, rcCreating, rcMemory, rcExhausted );
        {
            rc = KKeyRingInit(obj);
            if (rc == 0)
               *self = obj;
            else
                free(obj);
        }
    }
    
    return rc;
}

