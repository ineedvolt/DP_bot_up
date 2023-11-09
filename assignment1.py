"""
FIT2004 2023 S1 Assignment 1
Author: Ngeoh Khai Vin
Question 1
"""
class Vertex:
    def __init__(self, label) -> None:
        """
        Constructor for Vertex
        Input:
            label: the label of the vertex(u)
        """
        self.label = label  # to check its label(u)
        self.discovered = False
        self.visited = False
        self.previous = None
        self.edges = []
        self.distance = 0
        self.passenger = False  # to check this location include passenger or not

    def get_edges(self):
        return self.edges

    def get_distance(self):
        return self.distance

    def __int__(self):
        return self.get_distance()

    def __str__(self) -> str:
        return str(self.label) + "| " + ", ".join(str(e) for e in self.edges)


class Edge:
    def __init__(self, u, v, w, x):
        """
        Constructor for Edge
        Input:
            u: the starting vertex
            v: the ending vertex
            w: the weight of the edge for solo
            x: the weight of the edge for carpool
        """
        self.u = u
        self.v = v
        self.w = w  # solo weight
        self.x = x  # carpool weight

    def get_u(self):
        return self.u

    def get_v(self):
        return self.v

    def get_w(self):
        return self.w

    def get_x(self):
        return self.x

    def __str__(self):
        return f"({self.v.label}, {self.w}, {self.x})"

class MinHeap:
    """
    referenced the implementation from FIT1008 and FIT2004 PASS
    """

    def __init__(self, size) -> None:
        """
        Builds an empty min heap ADT with a fixed size

        Input:
            size: number of vertex

        let L = key locations(vertex)
        Time complexity: O(|L|)
        Aux space complexity: O(|L|)
        """
        self.array = [0 for i in range(size+1)]
        self.length = 0

    def swap(self, x, y):
        """
        Swap two vertex's position in the MinHeap array
        Time Complexity: O(1)
        """
        self.array[x], self.array[y] = self.array[y], self.array[x]

    def rise(self, element):
        """
        Adjusts the position of the element accordingly upwards
        let L = key locations(vertex)
        Time Complexity: O(log |L|)
        """
        parent = element // 2
        while parent >= 1:
            if int(self.array[parent]) > int(self.array[element]):
                self.swap(parent, element)
                element = parent
                parent = element // 2
            else:
                break

    def sink(self, element):
        """
        Adjusts the position of the element accordingly downwards

        let L = key locations(vertex)
        Time Complexity: O(log |L|)
        """
        child = 2 * element
        while child <= self.length:
            if child < self.length and int(self.array[child + 1]) < int(self.array[child]):
                child += 1
            if int(self.array[element]) > int(self.array[child]):
                self.swap(element, child)
                element = child
                child = 2 * element
            else:
                break

    def serve(self):
        """
        Returns the smallest number in the MinHeap's array

        let L = key locations(vertex)
        Time Complexity: O(log |L|)
        Aux Space complexity: O(1)
        """
        item = self.array[1]
        self.swap(1, self.length)
        self.length -= 1
        self.sink(1)
        return item

    def queue(self, item):
        """
        Adds an item into the MinHeap

        let L = key locations(vertex)
        Time Complexity: O(log |L|)
        Aux space complexity: O(1), in-place insertion of item and reordering of heap
        """
        self.length += 1
        self.array[len(self)] = item
        self.rise(self.length)


    def update(self, element):
        """
        Updates the distance of the vertex in the MinHeap

        let L = key locations(vertex)
        Time Complexity: O(log |L|)
        Aux Space complexity: O(1)
        """
        for i in range(1, self.length + 1):
            if self.array[i] == element:
                self.array[i] = element
                self.rise(i)
                self.sink(i)
                break


    def __len__(self):
        return self.length

    def __str__(self):
        return str(self.array)

class Graph:
    def __init__(self, roads):
        """
        Builds a graph ADT based on the roads.

        let L = vertices = key locations
        let R = edges = roads
        Time complexity: O(|L|+|R|) because of iterating through all routes to get the range of vertices, either the range of vertices is
        bigger than number or edges or the other way around.
        Aux space complexity: O(|L|+|R|) storing all vertices and edges as form of adjacency list
        """

        self.vertices = []
        self.add_vertices(roads)
        self.add_edges(roads)

    def add_vertices(self, roads):
        """
        adds the key locations in roads as vertices into an adjacency list

        let L = key locations
        let R = roads
        Time complexity: O(|L| + |R|) iterating through the roads and adding all vertices from 0 to max vertex
        Aux space complexity: O(|L|) storing all vertices
        """
        # to get the max vertex and later easy to add all vertices by range of the max vertex
        max_location = 0
        for r in roads:
            if r[0] > max_location:
                max_location = r[0]
            if r[1] > max_location:
                max_location = r[1]
        for i in range(0, max_location + 1):
            self.vertices.append(Vertex(i))

    def add_edges(self, roads):
        """
        adds all edges in roads as edges into vertex u instances.

        let L = key locations
        let R = roads
        Time complexity:  O(|L| + |R|) to iterate through all routes in roads
        Aux space complexity: O(|R|) to store all the edges in corresponding source vertex u.
        """
        for r in roads:
            u = self.vertices[r[0]]
            v = self.vertices[r[1]]
            w = r[2]
            x = r[3]
            route = Edge(u, v, w, x)
            u.edges.append(route)

    def Dijkstra(self, start, end):
        """
        Pathfinding algorithm to find the shortest route from start to end in a graph with no negative edges.
        It will return when it reaches the end vertex or the visited vertex but no passenger vertex.

        Input:
            start: the starting vertex
            end: the ending vertex
        Return:
            the shortest path from start to end with distance

        let L = vertices = key locations
        let R = edges = roads

        Time complexity: O(|R| log |L|)
        Aux Space complexity: O(|L|) because of queuing |L| vertices into MinHeap
        Total Space complexity: O(|L|)

        """
        discovered = MinHeap(len(self.vertices))
        discovered.queue(start)
        self.vertices[start.label].discovered = True
        while len(discovered) > 0:
            u = discovered.serve()
            u.visited = True
            # a variable to remind to use carpool weight if m == 1
            m = 0
            # check if there is a passenger in the previous vertices
            check_passenger = u
            while check_passenger is not None:
                if check_passenger.passenger == True:
                    m = 1
                    break
                else:
                    check_passenger = check_passenger.previous

            # return if current vertex == end
            if u == end:
                distance = u.distance
                path = []
                while u.previous is not None:
                    id = u.label
                    path.insert(0, id)
                    u = u.previous
                path.insert(0, u.label)
                return path, distance

            for e in u.get_edges():
                v = e.get_v()
                # use carpool weight
                if m == 1:
                    if v.discovered == False:
                        v.discovered = True
                        v.distance = u.distance + e.get_x()
                        v.previous = u
                        discovered.queue(v)
                    elif v.visited == False:
                        if v.distance >= u.distance + e.get_x():
                            v.distance = u.distance + e.get_x()
                            v.previous = u
                            discovered.update(v)
                    else:
                        # if v is visited and not a passenger location as visited passenger location has already used carpool weight,
                        # then return the current path from start to current vertex
                        if v.passenger == False:
                            distance = u.distance
                            path = []
                            path.append(u)
                            return path, distance

                # use solo weight
                elif m == 0:
                    if v.discovered == False:
                        v.discovered = True
                        v.distance = u.distance + e.get_w()
                        v.previous = u
                        discovered.queue(v)
                    elif v.visited == False:
                        if v.distance > u.distance + e.get_w():
                            v.distance = u.distance + e.get_w()
                            v.previous = u
                            discovered.update(v)
        return None, None

    def Dijkstra_Carpool(self, start, end):
        """
        Pathfinding algorithm to find the shortest route using carpool weight
        from start to end in a graph with no negative edges but will be called only when
        Dijkstra() returns a path which the last vertex is not the end vertex because it
        reached a visited vertex but no passenger vertex meant that its distance already finalized
        but it is actually the solo weight but not visited with carpool weight yet.

        Input:
            start: the starting vertex
            end: the ending vertex
        Return:
            the shortest path from start to end with distance

        let L = vertices = key locations
        let R = edges = roads

        Time complexity: O(|R| log |L|)
        Aux Space complexity: O(|L|) because of queuing |L| vertices into MinHeap
        Total Space complexity: O(|L|)
        """
        discovered = MinHeap(len(self.vertices))
        discovered.queue(start)
        self.vertices[start.label].discovered = True
        while len(discovered) > 0:
            u = discovered.serve()
            u.visited = True
            if u == end:
                d2 = u.distance
                path = []
                while u.previous is not None:
                    id = u.label
                    path.insert(0, id)
                    u = u.previous
                return path, d2
            for e in u.get_edges():
                v = e.get_v()
                if v.discovered == False:
                    v.discovered = True
                    v.distance = u.distance + e.get_x()
                    v.previous = u
                    discovered.queue(v)
                elif v.visited == False:
                    if v.distance > u.distance + e.get_x():
                        v.distance = u.distance + e.get_x()
                        v.previous = u
                        discovered.update(v)
        return None, None


    def backtracking_for_Dijkstra(self, from_vertex, start, end):
        """
        backtrack through vertex.previous until it is None and append into the path and also return shortest path
        from start to end if there is.

        Input:
            from_vertex: the vertex that Dijkstra() returns if it reached a visited vertex but no passenger vertex
            start: the starting vertex
            end: the ending vertex
        Return:
            the shortest path from start to visited vertex with distance and the shortest path from start to end with distance

        let L = vertices = key locations

        Time complexity: O(|L|) because of iterating the vertex.previous until it is None
        Aux Space complexity: O(|L|) because of appending |L| vertices into path
        """
        distance = from_vertex.distance
        path = []
        while from_vertex.previous is not None:
            id = from_vertex.label
            path.insert(0, id)
            from_vertex = from_vertex.previous
        path.insert(0, from_vertex.label)
        # also return the current shortest path from start to end(if there is) to later compare with the path.
        endv = self.vertices[end.label]
        distance2 = endv.distance
        path2 = []
        while endv.previous is not None:
            id2 = endv.label
            path2.insert(0, id2)
            endv = endv.previous
        if endv.label != self.vertices[start.label].label:
            return path, distance, None, None
        else:
            path2.insert(0, endv.label)
            return path, distance, path2, distance2

    def __str__(self):
        return "\n".join(str(vertex) for vertex in self.vertices)

    def __len__(self):
        return len(self.vertices)

def optimalRoute(start, end, passengers, roads):
    """
    Provides the shortest path from start to end and if the path passed through a passenger locations, the later path
    could use carpool lane else only could use non-carpool lane.

    Input:
        start: the starting vertex
        end: the ending vertex
        passengers: a list of passenger locations
        roads: a list of tuple(a,b,c,d) where a is start vertex, b is end vertex, c is the solo weight and d is the carpool weight
    Return:
        the shortest path from start to end

    let L = vertices = key locations
    let R = edges = roads
    Time complexity: O(|R| log |L|)  dominated by Dijkstra pathfinding algorithm
    Aux Space complexity: O(|L| + |R|) by adjacency list representation of graph
    """
    my = Graph(roads)  # O(|L| + |R|) aux space complexity

    # O(|P|) time complexity, loop through input passengers and change the passenger attribute of corresponding vertex to True
    for i in passengers:
        my.vertices[i].passenger = True

    # check if the end or start exceeds the number of vertices
    if len(my) < end or len(my) < start:
        return None

    path, d1 = my.Dijkstra(my.vertices[start], my.vertices[end])  # O(|R| log |L|) time complexity

    # if the path returned by Dijkstra() is the shortest path from start to end, return the path,
    # else Dijkstra_Carpool() to find the shortest path from last vertex in res to end.
    if path[len(path) - 1] == end or path is None:
        return path
    else:
        path, d1, path3, d3 = my.backtracking_for_Dijkstra(path[len(path)-1], my.vertices[start], my.vertices[end])
        my2 = Graph(roads)  # O(|L| + |R|) aux space complexity
        path2, d2 = my2.Dijkstra_Carpool(my2.vertices[path[len(path) - 1]], my2.vertices[end])  # O(|R| log |L|) time complexity

    # if the another path(path3) returned by Dijkstra() is the shortest path from start to end and not None, return the path,
    # else return the path from Dijkstra() with the path(path2) from Dijkstra_Carpool().
    if d3 is not None:
        if d1 + d2 < d3:
            return path + path2
        else:
            return path3
    else:
        return path + path2


"""
Question 2
"""

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
