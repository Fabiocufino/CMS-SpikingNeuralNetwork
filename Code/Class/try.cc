#include <iostream>
#include <vector>
#include <algorithm>

std::vector<std::pair<int, int>> uniquePairsInRange(const std::vector<std::pair<int, int>>& my_pairs, int start_idx, int end_idx) {
    std::vector<std::pair<int, int>> range_subset(my_pairs.begin() + start_idx, my_pairs.begin() + end_idx + 1);

    // Sort the subrange
    std::sort(range_subset.begin(), range_subset.end());

    // Remove consecutive duplicates
    auto it = std::unique(range_subset.begin(), range_subset.end(), [](const auto& a, const auto& b) {
        return a.first == b.first && a.second == b.second;
    });

    // Resize the vector to remove the duplicates
    range_subset.resize(std::distance(range_subset.begin(), it));

    return range_subset;
}

int main() {
    std::vector<std::pair<int, int>> my_pairs = {{0, 1}, {0, 1}, {0, 2}, {1, 2}, {2, 3}, {0, 1}, {1, 2}, {1, 2}};
    int start_idx = 2;
    int end_idx = 6;
    std::vector<std::pair<int, int>> unique_pairs = uniquePairsInRange(my_pairs, start_idx, end_idx);

    // Print the unique pairs
    for (const auto& pair : unique_pairs) {
        std::cout << "(" << pair.first << ", " << pair.second << ") ";
    }
    std::cout << std::endl;

    return 0;
}