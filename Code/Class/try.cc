#include <iostream>
#include <vector>
#include <algorithm>

int main() {
    std::vector<double> numbers = {1, 2, 4, 5};
    double x = 2; // Set the threshold value

    // Use binary search to find the position of the first element greater than x
    auto it = std::upper_bound(numbers.begin(), numbers.end(), x);

    // Erase elements before the position found by binary search
    numbers.erase(numbers.begin(), it);

    // Print the remaining elements
    for (auto num : numbers) {
        std::cout << num << " ";
    }
    std::cout << std::endl;

    return 0;
}
