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

#include <klib/extern.h>
#include <sysalloc.h>
#include <stdint.h>

LIB_EXPORT const uint32_t cp1252 [ 128 ] =
{
    /*0x80   0x81   0x82   0x83   0x84   0x85   0x86   0x87*/
    0x20AC,0x0000,0x201A,0x0192,0x201E,0x201E,0x201E,0x201E,
    /*0x88   0x89   0x8A   0x8B   0x8C   0x8D   0x8E   0x8F*/
    0x201E,0x201E,0x201E,0x201E,0x201E,0x0000,0x201E,0x0000,
    /*0x90   0x91   0x92   0x93   0x94   0x95   0x96   0x97*/
    0x0000,0x2018,0x2019,0x201C,0x201D,0x2022,0x2013,0x2014,
    /*0x98   0x99   0x9A   0x9B   0x9C   0x9D   0x9E   0x9F*/
    0x02DC,0x2122,0x0161,0x203A,0x0153,0x0000,0x017E,0x0178,
    /*0xA0   0xA1   0xA2   0xA3   0xA4   0xA5   0xA6   0xA7*/
    0x00A0,0x00A1,0x00A2,0x00A3,0x00A4,0x00A5,0x00A6,0x00A7,
    /*0xA8   0xA9   0xAA   0xAB   0xAC   0xAD   0xAE   0xAF*/
    0x00A8,0x00A9,0x00AA,0x00AB,0x00AC,0x00AD,0x00AE,0x00AF,
    /*0xB0   0xB1   0xB2   0xB3   0xB4   0xB5   0xB6   0xB7*/
    0x00B0,0x00B1,0x00B2,0x00B3,0x00B4,0x00B5,0x00B6,0x00B7,
    /*0xB8   0xB9   0xBA   0xBB   0xBC   0xBD   0xBE   0xBF*/
    0x00B8,0x00B9,0x00BA,0x00BB,0x00BC,0x00BD,0x00BE,0x00BF,
    /*0xC0   0xC1   0xC2   0xC3   0xC4   0xC5   0xC6   0xC7*/
    0x00C0,0x00C1,0x00C2,0x00C3,0x00C4,0x00C5,0x00C6,0x00C7,
    /*0xC8   0xC9   0xCA   0xCB   0xCC   0xCD   0xCE   0xCF*/
    0x00C8,0x00C9,0x00CA,0x00CB,0x00CC,0x00CD,0x00CE,0x00CF,
    /*0xD0   0xD1   0xD2   0xD3   0xD4   0xD5   0xD6   0xD7*/
    0x00D0,0x00D1,0x00D2,0x00D3,0x00D4,0x00D5,0x00D6,0x00D7,
    /*0xD8   0xD9   0xDA   0xDB   0xDC   0xDD   0xDE   0xDF*/
    0x00D8,0x00D9,0x00DA,0x00DB,0x00DC,0x00DD,0x00DE,0x00DF,
    /*0xE0   0xE1   0xE2   0xE3   0xE4   0xE5   0xE6   0xE7*/
    0x00E0,0x00E1,0x00E2,0x00E3,0x00E4,0x00E5,0x00E6,0x00E7,
    /*0xE8   0xE9   0xEA   0xEB   0xEC   0xED   0xEE   0xEF*/
    0x00E8,0x00E9,0x00EA,0x00EB,0x00EC,0x00ED,0x00EE,0x00EF,
    /*0xF0   0xF1   0xF2   0xF3   0xF4   0xF5   0xF6   0xF7*/
    0x00F0,0x00F1,0x00F2,0x00F3,0x00F4,0x00F5,0x00F6,0x00F7,
    /*0xF8   0xF9   0xFA   0xFB   0xFC   0xFD   0xFE   0xFF*/
    0x00F8,0x00F9,0x00FA,0x00FB,0x00FC,0x00FD,0x00FE,0x00FF
};