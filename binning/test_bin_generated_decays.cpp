/*
 * Quick tests for bin_generated_decays
 *
 * Just uses asserts and print statements
 */

#include "bin_generated_decays.cpp"

#include <cassert>
#include <functional>

void test(void)
{
    std::cout << "test running " << std::endl;
    assert(true);
}

int main(void)
{
    std::vector<std::function<void()>> testCases{&test};

    for (auto fcn = testCases.begin(); fcn != testCases.end(); ++fcn) {
        (*fcn)();
    }

    return 0;
}
