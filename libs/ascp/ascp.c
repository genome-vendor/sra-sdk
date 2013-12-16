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
* ==============================================================================
*
*/

#include "ascp-priv.h"

#include <kfg/config.h> /* KConfig */

#include <kfs/directory.h> /* KDirectory */

#include <klib/log.h> /* LOGERR */
#include <klib/out.h> /* OUTMSG */
#include <klib/printf.h> /* string_vprintf */
#include <klib/rc.h>
#include <klib/status.h> /* STSMSG */

#include <sysalloc.h> /* free */

#include <assert.h>
#include <stdlib.h> /* free */
#include <string.h> /* memset */

#define DISP_RC(rc, err) (void)((rc == 0) ? 0 : LOGERR(klogInt, rc, err))

#define RELEASE(type, obj) do { rc_t rc2 = type##Release(obj); \
    if (rc2 != 0 && rc == 0) { rc = rc2; } obj = NULL; } while (false)

#define STS_INFO 1
#define STS_FIN 3

static bool _StringStartsWith(const String *self, const char *buf) {
    size_t len = 0;
    assert(self && buf);
    len = string_size(buf);
    assert(len);
/*printf("_StringStartsWith(%.*s, %s)\n", self->len, self->addr, buf);*/
    if (self->len < len) {
        return false;
    }
    return string_cmp(self->addr, self->len, buf, len, len) == 0;
}

static
bool _StringHas(const String *self, const char *buf, String *res)
{
    String dummy;
    size_t len = 0;
    assert(self && buf);
    if (res == NULL) {
        res = &dummy;
    }
    len = string_size(buf);
    assert(len);
    StringInit(res, self->addr, self->size, self->len);
    while (true) {
        if (res->len < len) {
            StringInit(res, NULL, 0, 0);
            return false;
        }
        if (_StringStartsWith(res, buf)) {
            res->len = res->size = len;
            return true;
        }
        res->len = res->size = res->len - 1;
        ++res->addr;
    }
}

typedef struct {
    EAscpState s;
    char *msg;
    bool failure;
} SAscpState;
static void SAscpStateFini(SAscpState *self) {
    assert(self);

/* #if ! WINDOWS #endif */
    free(self->msg);

    memset(self, 0, sizeof *self);
}

static bool sStatus = false;

static rc_t parseAscpLine(const String *s, SAscpState *state, const char *name)
{
    bool debug = true;
    if (!sStatus) {
        debug = false;
    }

    assert(s && state);

    memset(state, 0, sizeof *state);

    if (_StringHas(s, "CHILD: ", NULL)) {
        if (sStatus) {
            OUTMSG(("%.*s\n",  s->len, s->addr));
        }
        state->s = eChild;
    }
    else if (_StringStartsWith(s, "Cannot open log file: ")) {
        if (debug) OUTMSG(("matched: LOG: '%.*s'\n", s->len, s->addr));
        state->s = eLog;
    }
    else if (_StringStartsWith(s, "The server's host key is not")) {
        if (debug) OUTMSG(("matched: KeySTART: '%.*s'\n", s->len, s->addr));
        state->s = eKeyStart;
    }
    else if (_StringHas(s, "no guarantee that the server is th", NULL)) {
        if (debug) OUTMSG(("matched: KeyIN: '%.*s'\n", s->len, s->addr));
        state->s = eKeyIn;
    }
    else if (_StringHas(s, "think it is.", NULL)) {
        if (debug) OUTMSG(("matched: KeyIN: '%.*s'\n", s->len, s->addr));
        state->s = eKeyIn;
    }
    else if (_StringHas(s, "The server's rsa2 key fingerprint ", NULL)) {
        if (debug) OUTMSG(("matched: KeyIN: '%.*s'\n", s->len, s->addr));
        state->s = eKeyIn;
    }
    else if (_StringHas(s, "ssh-rsa 1024 ", NULL)) {
        if (debug) OUTMSG(("matched: KeyIN: '%.*s'\n", s->len, s->addr));
        state->s = eKeyIn;
    }
    else if (_StringHas(s, "If you trust this host, enter ", NULL)) {
        if (debug) OUTMSG(("matched: KeyIN: '%.*s'\n", s->len, s->addr));
        state->s = eKeyIn;
    }
    else if (_StringHas(s, "PuTTY's cache and carry on connect", NULL)) {
        if (debug) OUTMSG(("matched: KeyIN: '%.*s'\n", s->len, s->addr));
        state->s = eKeyIn;
    }
    else if (_StringHas(s, "If you want to carry on connecting", NULL)) {
        if (debug) OUTMSG(("matched: KeyIN: '%.*s'\n", s->len, s->addr));
        state->s = eKeyIn;
    }
    else if (_StringHas(s, "adding the key to the cache, enter", NULL)) {
        if (debug) OUTMSG(("matched: KeyIN: '%.*s'\n", s->len, s->addr));
        state->s = eKeyIn;
    }
    else if (_StringHas(s, " you do not trust this host, press", NULL)) {
        if (debug) OUTMSG(("matched: KeyIN: '%.*s'\n", s->len, s->addr));
        state->s = eKeyIn;
    }
    else if (_StringHas(s, "connection.", NULL)) {
        if (debug) OUTMSG(("matched: KeyIN: '%.*s'\n", s->len, s->addr));
        state->s = eKeyIn;
    }
    else if (_StringHas(s, "Store key in cache? (y/n) ", NULL)) {
        if (debug) OUTMSG(("matched: KeyEND: '%.*s'\n", s->len, s->addr));
        state->s = eKeyEnd;
    }
    else if (string_chr(s->addr, s->len, '%') != NULL) {
        if (debug) OUTMSG(("matched: PROGRESS: '%.*s'\n", s->len, s->addr));
        state->s = eProgress;
        if ((state->msg = string_dup(s->addr, s->len)) == NULL) {
            return RC(rcExe, rcStorage, rcAllocating, rcMemory, rcExhausted);
        }
    }
    else if (_StringStartsWith(s, "Completed: ")) {
        if (debug) OUTMSG(("matched: COMPLETED: '%.*s'\n", s->len, s->addr));
        state->s = eCompleted;
        if ((state->msg = string_dup(s->addr, s->len)) == NULL) {
            return RC(rcExe, rcStorage, rcAllocating, rcMemory, rcExhausted);
        }
    }
    else if (_StringStartsWith(s, "Partial Completion: ")) {
        if (debug) OUTMSG(("matched: END: '%.*s'\n", s->len, s->addr));
        state->s = eEnd;
        if ((state->msg = string_dup(s->addr, s->len)) == NULL) {
            return RC(rcExe, rcStorage, rcAllocating, rcMemory, rcExhausted);
        }
    }    
    else if (_StringStartsWith(s, "Connection abandoned.")) {
/* printed in caller     printf("matched: END: '%.*s'\n", s->len, s->addr); */
        state->s = eFailed;
        if ((state->msg = string_dup(s->addr, s->len)) == NULL) {
            return RC(rcExe, rcStorage, rcAllocating, rcMemory, rcExhausted);
        }
    }
    else if (_StringHas(s, "failed to open connection to remot", NULL)) {
/* printed in caller       printf("matched: END: '%.*s'\n", s->len, s->addr); */
        state->s = eFailed;
        if ((state->msg = string_dup(s->addr, s->len)) == NULL) {
            return RC(rcExe, rcStorage, rcAllocating, rcMemory, rcExhausted);
        }
    }
    else if (_StringHas(s, "exiting", NULL)) {
        state->s = eFailed;
        if ((state->msg = string_dup(s->addr, s->len)) == NULL) {
            return RC(rcExe, rcStorage, rcAllocating, rcMemory, rcExhausted);
        }
    }
    else if (_StringStartsWith(s, "Session Stop  (Error: Disk write ")) {
   if (debug) OUTMSG(("matched: Disk write failed: '%.*s'\n", s->len, s->addr));
        state->s = eWriteFailed;
        if ((state->msg = string_dup(s->addr, s->len)) == NULL) {
            return RC(rcExe, rcStorage, rcAllocating, rcMemory, rcExhausted);
        }
    }
    else if (_StringStartsWith(s, "Session Stop ")) {
        if (debug) OUTMSG(("matched: COMPLETED: '%.*s'\n", s->len, s->addr));
        state->s = eFailed;
        if ((state->msg = string_dup(s->addr, s->len)) == NULL) {
            return RC(rcExe, rcStorage, rcAllocating, rcMemory, rcExhausted);
        }
    }
    else if (_StringHas(s, " bits/sec), in 1 file", NULL)) {
        if (debug) OUTMSG(("matched: END: '%.*s'\n", s->len, s->addr));
        state->s = eEnd;
        if ((state->msg = string_dup(s->addr, s->len)) == NULL) {
            return RC(rcExe, rcStorage, rcAllocating, rcMemory, rcExhausted);
        }
    }
    else if (_StringStartsWith(s, name)) {
        /* in the beginning:
           line starting by dest file name, then some white characters */
        if (debug) OUTMSG(("matched: PROGRESS: '%.*s'\n", s->len, s->addr));
        state->s = eProgress;
        if ((state->msg = string_dup(s->addr, s->len)) == NULL) {
            return RC(rcExe, rcStorage, rcAllocating, rcMemory, rcExhausted);
        }
    }
    else {
        OUTMSG(("LINE = (%d) '%.*s'\n", s->len, s->len, s->addr));
/*      assert(0); */
    }
    return 0;
}

rc_t ascpParse(const char *buf, size_t len, const char *filename,
    EAscpState *state, String *line)
{
    bool failure = false;
    const char *p = buf;
    int64_t l = len;
    assert(buf && len && filename && state && line);
    StringInit(line, NULL, 0, 0);
    while (true) {
        const char *n = string_chr(p, l, '\n');
        const char *r = string_chr(p, l, '\r');
        if (n == NULL) {
            if (r != NULL) {
                n = r;
            }
        }
        else {
            if (r != NULL) {
                if (r < n) {
                    n = r;
                }
            }
        }
        if (n != NULL) {
            StringInit(line, p, n - p, n - p);
            l -= n - p + 1;
        }
        else {
            StringInit(line, p, l, l);
        }
        if (line->addr && line->len && line->addr[line->len - 1] == '\r') {
            line->len = line->size = line->len - 1;
        }
        if (line->addr && line->len && line->addr[0] == '\r') {
            ++line->addr;
            line->len = line->size = line->len - 1;
        }
        if (line->len != 0) {
            SAscpState full;
            rc_t rc = parseAscpLine(line, &full, filename);
            if (rc != 0) {
                return rc;
            }
            switch (full.s) {
                case eChild:
                    break;
                case eUnknown:
                    switch (*state) {
                        case eKeyStart:
                        case eKeyMayBeIn:
                        case eKeyIn:
                            *state = eKeyMayBeIn;
                            break;
                        case eCompleted:
                        case eFailed:
                        case eWriteFailed:
                            *state = eEnd;
                            /* report to user */
                            break;
                        case eProgress:
                            if (sStatus) {
                                OUTMSG(("\n"));
                            }
/*                          no break; */
                        default:
                            *state = eUnknown;
                            /* report to user */
                            break;
                    }
                    break;
                case eFailed:
                    if (*state == eProgress) {
                        if (sStatus) {
                            OUTMSG(("\n"));
                        }
                    }
                    failure = true;
                    *state = full.s;
                    if (sStatus) {
                        OUTMSG(("%s\n", full.msg));
                    }
/*                  no break; */
                    break;
                case eWriteFailed:
                    if (*state == eProgress) {
                        if (sStatus) {
                            OUTMSG(("\n"));
                        }
                    }
                    failure = true;
                    *state = full.s;
                    if (sStatus) {
                        OUTMSG(("%s\n", full.msg));
                    }
/*                  no break; */
                    break;
                case eCompleted:
                    if (*state == eProgress) {
                        if (sStatus) {
                            OUTMSG(("\n"));
                        }
                    }
                    failure = false;
                    *state = full.s;
                    if (sStatus) {
                        OUTMSG(("%s\n", full.msg));
                    }
/*                  no break; */
                    break;
                case eProgress:
                    if (*state == eProgress) {
                        if (sStatus) {
                            OUTMSG(("\r"));
                        }
                    }
                    *state = full.s;
                    if (sStatus) {
                        OUTMSG(("%s", full.msg));
                    }
                    break;
                case eEnd:
                    if (*state == eProgress) {
                        if (sStatus) {
                            OUTMSG(("\n"));
                        }
                    }
                    *state = full.s;
                    if (sStatus) {
                        OUTMSG(("%s\n", full.msg));
                    }
                    /* report to user */
                    break;
                default:
                    *state = full.s;
                    break;
            }
            SAscpStateFini(&full);
        }
        if (n == NULL || l <= 0) {
            break;
        }
        if (*state == eKeyEnd) {
            String end;
            if (_StringHas(line, "Store key in cache? (y/n) ", &end)) {
                if (n > end.addr + end.len) {
                    l += n - end.addr + end.len;
                    n = end.addr + end.len - 1;
                }
            }
        }
        p = n + 1;
        if (p >= buf + len) {
            break;
        }
    }
    if (sStatus) {
        STSMSG(STS_FIN, ("%.*s", len, buf));
    }
    return 0;
}

static bool _KConfigAscpDisabled(const KConfig *self, bool status) {
    bool disabled = false;
    const char path[] = "tools/ascp/disabled";
    rc_t rc = KConfigReadBool(self, path, &disabled);
    if (rc != 0) {
        if (rc != SILENT_RC(rcKFG, rcNode, rcOpening, rcPath, rcNotFound)) {
            DISP_RC(rc, path);
        }
        else {
            if (status) {
                STSMSG(STS_DBG, ("'%s': not found in configuration", path));
            }
        }
        disabled = false;
    }
    else {
        if (status) {
            STSMSG(STS_DBG, ("'%s' = '%s'", path, disabled ? "true" : "false"));
        }
    }
    return disabled;
}

static String* _KConfigAscpString(const KConfig *self,
    const char *path, const char *name)
{
    String *ascp = NULL;
    rc_t rc = KConfigReadString(self, path, &ascp);
    if (rc == 0) {
        assert(ascp);
/*      STSMSG(STS_INFO, ("Using %s from configuration: '%s'",
            name, ascp->addr)); */
        return ascp;
    }
    else {
        if (rc != SILENT_RC(rcKFG, rcNode, rcOpening, rcPath, rcNotFound)) {
            DISP_RC(rc, path);
        }
        else {
/*          STSMSG(STS_DBG, ("'%s': not found in configuration", path)); */
        }
        free(ascp);
        return NULL;
    }
}

/******************************************************************************/

static int _SilentSystem(const char *fmt, ...) {
    rc_t rc = 0;
    char buffer[4096];
    size_t num_writ = 0;
    va_list args;
    va_start(args, fmt);
    rc = string_vprintf(buffer, sizeof buffer, &num_writ, fmt, args);
    va_end(args);
    if (rc != 0) {
        LOGERR(klogInt, rc, "while making ascp command line");
        return 1;
    }
    return silent_system(buffer);
}

static bool _SystemHelp(const char *command, bool status) {
    int value = 0;
    if (status) {
        STSMSG(STS_DBG, ("Checking '%s'", command));
    }
    value = _SilentSystem("\"%s\" -h", command);
    if (value == 0) {
        if (status) {
            STSMSG(STS_INFO, ("Using '%s'", command));
        }
        return true;
    }
    else {
        if (status) {
            STSMSG(STS_DBG, ("'%s': not found", command));
        }
        return false;
    }
}

static rc_t _KConfigGetAscp(const KConfig *self,
    const char **ascp_bin, const char **private_file)
{
    String *bin = NULL;
    String *key = NULL;
    assert(self && ascp_bin && private_file);
    *ascp_bin = *private_file = NULL;
    bin = _KConfigAscpString(self, "tools/ascp/path", "ascp");
    key = _KConfigAscpString(self, "tools/ascp/key", "Aspera key");
    if (bin != NULL && key != NULL) {
        *ascp_bin = string_dup_measure(bin->addr, NULL);
        *private_file = string_dup_measure(key->addr, NULL);
        free(bin);
        free(key);
        if (*ascp_bin == NULL || *private_file == NULL) {
            free((void*)*ascp_bin);
            free((void*)*private_file);
            *ascp_bin = *private_file = NULL;
            return RC(rcNS, rcStorage, rcAllocating, rcMemory, rcExhausted);
        }
        return 0;
    }
    free(bin);
    free(key);
    return 0;
}
static bool _KDirectoryFileFound(const KDirectory *self,
    const char *path, bool status)
{
    KPathType type = kptNotFound;
    if (status) {
        STSMSG(STS_DBG, ("Checking '%s'", path));
    }
    type = KDirectoryPathType(self, path);
    if ((type & ~kptAlias) == kptFile) {
        if (status) {
            STSMSG(STS_DBG, ("'%s': found", path));
        }
        return true;
    }
    else {
        if (status) {
            STSMSG(STS_DBG, ("'%s': not found", path));
        }
        return false;
    }
}

LIB_EXPORT rc_t CC ascp_locate(const char **ascp_bin, const char **private_file,
    bool use_config, bool status)
{
    rc_t rc = 0;
    KConfig *cfg = NULL;
    if (ascp_bin == NULL || private_file == NULL) {
        return RC(rcNS, rcFile, rcCopying, rcParam, rcNull);
    }
    *ascp_bin = *private_file = NULL;
    rc = KConfigMake(&cfg, NULL);
    if (rc != 0) {
        return rc;
    }
    if (_KConfigAscpDisabled(cfg, status)) {
        if (status) {
            STSMSG(STS_INFO, ("Use of Aspera transfer is disabled "
                "by the configuration, using HTTP transfer"));
        }
    }
    else {
        KDirectory *dir = NULL;
        const char *bin = NULL;
        const char *key = NULL;
        rc = _KConfigGetAscp(cfg, ascp_bin, private_file);
        if (*ascp_bin != NULL) {
            assert(*private_file && !rc);
            RELEASE(KConfig, cfg);
            return 0;
        }
        rc = KDirectoryNativeDir(&dir);
        if (rc != 0) {
            return rc;
        }
        while (ascp_path(&bin, &key)) {
            if (_SystemHelp(bin, status)) {
                if (_KDirectoryFileFound(dir, key, status)) {
                    *ascp_bin = string_dup_measure(bin, NULL);
                    *private_file = string_dup_measure(key, NULL);
                    if (*ascp_bin == NULL || *private_file == NULL) {
                        free((void*)*ascp_bin);
                        free((void*)*private_file);
                        *ascp_bin = *private_file = NULL;
                        return RC(rcNS,
                            rcStorage, rcAllocating, rcMemory, rcExhausted);
                    }
                    break;
                }
            }
        }
        RELEASE(KDirectory, dir);
    }
    RELEASE(KConfig, cfg);
    return rc;
}

LIB_EXPORT rc_t CC aspera_get(
    const char *ascp_bin, const char *private_file, const char *src,
    const char *dest, AscpOptions *opt)
{
    AscpOptions dummy;
    bool status = false;
    uint64_t prev = 0;
    int attempt = 0;
    KDirectory *dir = NULL;
    TQuitting *quitting = NULL;
    rc_t rc = KDirectoryNativeDir(&dir);
    if (rc != 0) {
        return rc;
    }
    if (ascp_bin == NULL || private_file == NULL ||
        src == NULL || dest == NULL)
    {
        return RC(rcNS, rcFile, rcCopying, rcParam, rcNull);
    }
    if (opt == NULL) {
        memset(&dummy, 0, sizeof dummy);
        opt = &dummy;
    }

    sStatus = status = opt->status;
    quitting = opt->quitting;

    while (true) {
        rc = run_ascp(ascp_bin, private_file, src, dest, opt);
        if (rc == 0) {
            if (status) {
                STSMSG(STS_DBG, ("ascp finished with success"));
            }
            break;
        }
        else if (rc == SILENT_RC(rcExe,
            rcProcess, rcExecuting, rcMemory, rcExhausted))
        {
            if (status) {
                STSMSG(STS_DBG, ("ascp failed: %R", rc));
            }
            break;
        }
        else {
            rc_t rc = 0;
            uint64_t size = 0;
            if (quitting != NULL) {
                rc = quitting();
                if (rc != 0) {
                    break;
                }
            }
            if (status) {
                STSMSG(STS_DBG, ("ascp failed: %R", rc));
            }
            rc = KDirectoryFileSize(dir, &size, dest);
            if (rc != 0 || size < prev) {
            if (status) {
STSMSG(0, ("KDirectoryFileSize after ascp run: "
"rc = %ld, size = %ld", rc, size));
            }
                break;
            }
            else if (size > prev) {
                if (status) {
                    STSMSG(STS_INFO, ("  fasp download failed. %ld bytes "
                        "received so far. Retrying...", size));
                }
                attempt = 0;
                prev = size;
            } else {
                if (attempt++ > 3) {
                    break;
                }
                if (status) {
                    STSMSG(STS_INFO, ("  fasp download failed. %ld bytes "
                        "received so far. Retrying %d...", size, attempt));
                }
            }
        }
    }
    RELEASE(KDirectory, dir);
    return rc;
}

rc_t mkAscpCmd(const char *ascp_bin, const char *private_file,
    const char *src, const char *dest, const AscpOptions *opt,
    char *const argv[], size_t argvSz)
{
    rc_t rc = 0;
    return rc;
}
