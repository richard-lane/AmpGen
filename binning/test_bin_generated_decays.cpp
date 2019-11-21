/*
 * Quick tests for bin_generated_decays
 *
 * Just uses exceptions and print statements
 *
 */

#include "bin_generated_decays.cpp"

#include <cassert>
#include <exception>
#include <functional>

/*
 * Exception thrown if a test fails
 *
 */
class TestException : public std::exception
{
  private:
    std::string message{""};

  public:
    TestException(std::string testName) { message = "Test " + testName + " failed"; }
    const char* what() const noexcept override { return message.c_str(); }
};

/*
 * Test that the splitVectors function correctly splits a vector of 3 10-element vectors into a vector of 3 (4, 4, 2)
 * element vectors
 *
 */
void testSplitVectors(void)
{
    std::string testName = "testSplitVectors";

    std::vector<std::vector<double>> v{{1, 2, 3, 4, 5, 6, 7, 8, 9, 10},
                                       {2, 4, 6, 8, 10, 12, 14, 16, 18, 20},
                                       {10, 20, 30, 40, 50, 60, 70, 80, 90, 100}};

    std::vector<std::vector<std::vector<double>>> splitV = splitVectors(v, 4);
    std::vector<std::vector<std::vector<double>>> expectedSplitV{{{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10}},
                                                                 {{2, 4, 6, 8}, {10, 12, 14, 16}, {18, 20}},
                                                                 {{10, 20, 30, 40}, {50, 60, 70, 80}, {90, 100}}};
    if (splitV != expectedSplitV) {
        throw TestException(testName);
    }

    std::cout << testName + " Passed" << std::endl;
}

int test_bin_generated_decays(void)
{
    std::vector<std::function<void()>> testCases{&testSplitVectors};

    for (auto fcn = testCases.begin(); fcn != testCases.end(); ++fcn) {
        (*fcn)();
    }

    return 0;
}
