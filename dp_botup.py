def select_sections(occupancy_probability):
    """
    Use bottom up dynamic programming to find the minimum total occupancy for the selected n sections
    Solve from the bottom which is the second last row to the top iteratively with the base case which
    is the last row of the memoization table.

    Input:
        occupancy_probability: input list of n rows and m columns of numbers
    Output:
        total occupancy for the selected n sections, with the sections_location in list of n tuples in the form of (i, j)

    let n = number of rows in occupancy_probability
    let m = number of columns in occupancy_probability

    Time Complexity: O(nm)  <-  O(nm+n+m+n+n),dominated by iterating through the n+1 columns and m rows of the
                                            occupancy_probability for the bottom up approach.
    Aux Space complexity: O(nm)  <-  O(nm+n) where create of memoization table which is in size of n*m
                                            and n is from output sections_location list which is in size n
    Total Space Complexity: O(nm)  <-  0(nm)+O(nm) =  Aux + input

    """
    n, m = len(occupancy_probability), len(occupancy_probability[0])

    # create a memoization table with n+1 rows and m columns
    memo = [[-1] * m for _ in range(n+1)]    # O(nm) aux space complexity and O(n) time complexity

    # base case
    for j in range(m):  # O(m) time complexity
        memo[n][j] = 0

    # bottom up starting from the second last row and sum the below row to the current row
    # and going towards to the first row with last row as base case.
    for i in range(n-1, -1, -1):  # O(nm) time complexity
        for j in range(m):
            min_sum = float('inf')
            if min_sum > memo[i+1][j]:  # adjacent vertical
                min_sum = memo[i+1][j]
            # to check if it is not before the first column
            if j - 1 != -1:
                if min_sum > memo[i+1][j-1]:  # adjacent left diagonal
                    min_sum = memo[i+1][j-1]
            # to check if it is not over the last column
            if j + 1 != m:
                if min_sum > memo[i+1][j+1]:  # adjacent right diagonal
                    min_sum = memo[i+1][j+1]
            memo[i][j] = occupancy_probability[i][j] + min_sum

    # find the minimum sum position in the first row
    minimum_total_occupancy = memo[0][0]
    start = 0
    for j in range(m):    # O(m) time complexity
        if minimum_total_occupancy > memo[0][j]:
            minimum_total_occupancy = memo[0][j]
            start = j

    # backtrack from the minimum starting point from first row to the second last row of the memoization table
    i, j = 0, start
    sections_location = []
    sections_location.append((i, j))
    for i in range(n-1):    # O(n) time complexity
        if j+1 == m:    # to check if it is not over the last column
            if memo[i+1][j] < memo[i+1][j-1]:
                sections_location.append((i+1, j))
            elif memo[i+1][j-1] < memo[i+1][j]:
                sections_location.append((i+1, j-1))
                j = j-1
        elif j-1 == -1:   # to check if it is not before the first column
            if memo[i+1][j] < memo[i+1][j+1]:
                sections_location.append((i+1, j))
            elif memo[i+1][j+1] < memo[i+1][j]:
                sections_location.append((i+1, j+1))
                j = j+1
        else:
            if memo[i+1][j-1] < memo[i+1][j] and memo[i+1][j-1] < memo[i+1][j+1]:
                sections_location.append((i+1, j-1))
                j = j - 1
            elif memo[i+1][j+1] < memo[i+1][j] and memo[i+1][j+1] < memo[i+1][j-1]:
                sections_location.append((i+1, j+1))
                j = j + 1
            else:
                sections_location.append((i+1, j))
    return [minimum_total_occupancy, sections_location]
