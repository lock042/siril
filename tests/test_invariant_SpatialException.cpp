#include <gtest/gtest.h>
#include <string>
#include <cstring>
#include <vector>

// Forward declare the SpatialException class from the production code
namespace htmesh {
class SpatialException;
}

// Include the actual production header
#include "subprojects/htmesh/SpatialException.h"

class BufferOverflowSecurityTest : public ::testing::TestWithParam<std::string> {};

TEST_P(BufferOverflowSecurityTest, NoBufferOverflowOnOversizedInput) {
    // Invariant: Buffer reads never exceed declared length.
    // SpatialException must safely handle strings larger than internal buffer.
    
    std::string payload = GetParam();
    
    // Attempt to construct exception with potentially oversized input.
    // The test passes if no crash/overflow occurs and the object remains valid.
    try {
        htmesh::SpatialException exc(payload.c_str());
        
        // Verify the exception object is in a valid state by calling what()
        // If buffer overflow occurred, this would likely crash or return garbage.
        const char* msg = exc.what();
        ASSERT_NE(msg, nullptr);
        
        // Verify the returned message is null-terminated and reasonable length
        size_t msg_len = std::strlen(msg);
        EXPECT_LT(msg_len, 10000) << "Message length suspiciously large";
        
    } catch (const std::exception& e) {
        // Exception is acceptable; overflow would cause segfault, not throw
        SUCCEED();
    }
}

INSTANTIATE_TEST_SUITE_P(
    AdversarialInputs,
    BufferOverflowSecurityTest,
    ::testing::Values(
        std::string("Normal error message"),                           // Valid input
        std::string(256, 'A'),                                         // Boundary: 256 bytes
        std::string(1024, 'B'),                                        // 1KB payload
        std::string(10000, 'C'),                                       // 10KB payload (10x typical)
        std::string(100000, 'D')                                       // 100KB payload (extreme)
    )
);

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}