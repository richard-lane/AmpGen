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
 * Gets caught by the test runner if raised during a testcase
 *
 */
class TestException : public std::exception
{
  public:
    std::string message{""};
    TestException(std::string testName) { message = "Test " + testName + " failed"; }
    const char* what() const noexcept override { return message.c_str(); }
};

inline void check(bool condition, std::string testName)
{
    if (!condition) {
        throw TestException(testName);
    }

    std::cout << testName + " Passed" << std::endl;
}

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
    check(splitV == expectedSplitV, testName);
}

/*
 * Test vector averaging works
 */
void testVectorAvg(void)
{
    std::string         testName = "testVectorAvg";
    std::vector<double> myVector{1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21};
    double              expectedAvg{11};

    double avg = vectorAvg(myVector);

    check(avg == expectedAvg, testName);
}

/*
 * Test vector std dev works
 */
void testVectorStdDev(void)
{
    std::string         testName = "testVectorStdDev";
    std::vector<double> myVector{1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21};
    double              expectedStdDev{6.32455532};

    double stdDev = vectorStdDev(myVector);

    check(std::abs(stdDev - expectedStdDev) < 1e-8, testName);
}

/*
 * Test sorting a vector of vectors of pairs
 */
void testSortVectorofVectorofPairs(void)
{
    std::string testName = "testSortVectorofVectorofPairs";

    std::vector<std::vector<std::pair<double, double>>> vectors{
        {std::make_pair(1, 5), std::make_pair(3, 3), std::make_pair(2, 4), std::make_pair(5, 1), std::make_pair(4, 2)},
        {std::make_pair(10, 20), std::make_pair(10, 30), std::make_pair(5, 25)}};

    std::vector<std::vector<std::pair<double, double>>> expectedSortedVectors{
        {std::make_pair(5, 1), std::make_pair(4, 2), std::make_pair(3, 3), std::make_pair(2, 4), std::make_pair(1, 5)},
        {std::make_pair(10, 20), std::make_pair(5, 25), std::make_pair(10, 30)}};

    std::vector<std::vector<std::pair<double, double>>> sortedVectors(vectors);
    sortVectorOfPairs(sortedVectors);

    check(sortedVectors == expectedSortedVectors, testName);
}

void test_bin_generated_decays(void)
{
    // List your functions to test here
    std::vector<std::function<void()>> testCases{
        &testSplitVectors, &testVectorAvg, &testVectorStdDev, &testSortVectorofVectorofPairs};

    for (auto fcn = testCases.begin(); fcn != testCases.end(); ++fcn) {
        try {
            (*fcn)();
        } catch (TestException e) {
            std::cout << e.message << std::endl;
        }
    }
}
