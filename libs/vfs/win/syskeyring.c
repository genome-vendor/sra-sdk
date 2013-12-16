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

#include "keyring-priv.h"

#include <kfg/config.h>
#include <klib/text.h>

rc_t StartKeyRing(struct KStream** ipc)
{
    KConfig* kfg;
    rc_t rc = KConfigMake(&kfg, NULL);
    if (rc == 0)
    {
        const KConfigNode *node;
        char path[] = "$(APPPATH)";
        char buf[4096];
        size_t num_read;
    
        rc_t rc=KConfigOpenNodeRead(kfg, &node, path, string_measure(path, NULL), buf);
        if (rc == 0) 
        {
            rc = KConfigNodeRead(node, 0, buf, sizeof(buf), &num_read, NULL);
            if (rc == 0)
            {
printf("apppath='%s'\n", buf);        
            }
            KConfigNodeRelease(node);
        }
        rc = KConfigRelease(kfg);
    }
    return rc;
}
