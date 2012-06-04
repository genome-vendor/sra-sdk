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

#ifndef _h_ktst_unit_test_
#define _h_ktst_unit_test_

#include <ktst/unit_test_suite.hpp>

////////////////////////////////////////////////////////////////////////////////

#define FIXTURE_TEST_CASE(test_name, F) \
namespace ncbi { namespace NK { \
class C##test_name : public TestCase, F { \
public: \
    C##test_name(void* globalFixture, F* fixture) : TestCase(#test_name) \
        , _globalFixture \
            (static_cast<AUTO_TEST_CASE_FIXTURE*>(globalFixture)) \
       {} \
    void test_method(); \
private: \
    AUTO_TEST_CASE_FIXTURE* GET_GLOBAL_FIXTURE(void) const \
    { return _globalFixture; } \
    const F* GET_FIXTURE(void) const { return this; } \
    F* GET_FIXTURE(void) { return this; } \
    AUTO_TEST_CASE_FIXTURE* _globalFixture; \
}; \
class C##test_name##Invoker : TestInvoker { \
public: C##test_name##Invoker(void) : TestInvoker(#test_name) \
{ GetTestSuite()->Add(this); } \
private: virtual void Run(void* globalFixture) { \
        F fixture; \
        C##test_name t(globalFixture, &fixture); \
        try { \
            t.test_method(); \
            SetErrorCounter(t.GetErrorCounter()); \
        } catch (...) { \
            SetErrorCounter(t.GetErrorCounter()); \
            if (GetErrorCounter() == 0)\
            {\
                SetErrorCounter(1);\
            }\
            throw; \
        } \
    } \
}; \
static C##test_name##Invoker _##test_name##instance; \
} } /* ncbi::NK */ \
void ncbi::NK::C##test_name::test_method()

#define TEST_CASE(test_name) FIXTURE_TEST_CASE(test_name, ncbi::NK::Empty)

#define PROCESS_TEST_CASE(test_name, rc, timeout) \
    namespace ncbi { namespace NK { static void test_name##_impl(); }}\
    TEST_CASE(test_name) \
    {\
        REQUIRE_EQ(TestEnv::RunProcessTestCase(test_name##_impl, timeout), rc);\
    }\
    void ::ncbi::NK::test_name##_impl()

//TODO: PROCESS_FIXTURE_TEST_CASE

#define FIXTURE_TEST_SUITE( suite_name, F ) \
typedef F AUTO_TEST_CASE_FIXTURE; \
int suite_name(int argc, char* argv[]) { \
    ncbi::NK::TestEnv args(argc, argv); \
    if (args.catch_system_errors) { \
        args.set_handlers(); \
    } \
    ncbi::NK::counter_t ec = ncbi::NK::Main<AUTO_TEST_CASE_FIXTURE>(argc, argv, #suite_name); \
    return ec == 0 ? 0 : -ec; /* positive rc represents the signal that killed the process */ \
} 

#define TEST_SUITE( suite_name ) \
FIXTURE_TEST_SUITE(suite_name, ncbi::NK::Empty)

#endif// _h_ktst_unit_test_
