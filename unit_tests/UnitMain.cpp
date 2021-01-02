#include <gtest/gtest.h>
#include "stk_util/parallel/Parallel.hpp"
 
int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    stk::parallel_machine_init(&argc, &argv);

    auto status = RUN_ALL_TESTS();

    stk::parallel_machine_finalize();
    return status;
}
