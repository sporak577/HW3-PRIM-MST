import numpy as np
import heapq
from typing import Union

"""
I use chatgpt to help me identify code errors and give me inspiration
on how to construct my code. 
"""

class Graph:

    def __init__(self, adjacency_mat: Union[np.ndarray, str]):
        """
    
        Unlike the BFS assignment, this Graph class takes an adjacency matrix as input. `adjacency_mat` 
        can either be a 2D numpy array of floats or a path to a CSV file containing a 2D numpy array of floats.

        In this project, we will assume `adjacency_mat` corresponds to the adjacency matrix of an undirected graph.
    
        """
        
        if type(adjacency_mat) == str:
            """
            if the adj matrix is a string, it assumes the string represents a file path.
            it then loads the adjacency matrix from this file using the method _load etc. 
            """
            self.adj_mat = self._load_adjacency_matrix_from_csv(adjacency_mat)
        #
        elif type(adjacency_mat) == np.ndarray:
            """
            if the matrix is a numpy array, it assumes it's already an adjacency matrix 
            and assigns it direclty to self.adj_mat.
            """
            self.adj_mat = adjacency_mat
        else: 
            """
            if the matrix is neither a string nor a NumPy array, it raises an error. 
            """
            raise TypeError('Input must be a valid path or an adjacency matrix')
        """
        clas is working with a minimum spanning tree, initializing self.mst to None
        """
        self.mst = None

    def _load_adjacency_matrix_from_csv(self, path: str) -> np.ndarray:
        with open(path) as f:
            return np.loadtxt(f, delimiter=',')

    def construct_mst(self):
        """
    
        TODO: Given `self.adj_mat`, the adjacency matrix of a connected undirected graph, implement Prim's 
        algorithm to construct an adjacency matrix encoding the minimum spanning tree of `self.adj_mat`. 
            
        `self.adj_mat` is a 2D numpy array of floats. Note that because we assume our input graph is
        undirected, `self.adj_mat` is symmetric. Row i and column j represents the edge weight between
        vertex i and vertex j. An edge weight of zero indicates that no edge exists. 
        
        This function does not return anything. Instead, store the adjacency matrix representation
        of the minimum spanning tree of `self.adj_mat` in `self.mst`. We highly encourage the
        use of priority queues in your implementation. Refer to the heapq module, particularly the 
        `heapify`, `heappop`, and `heappush` functions.

        """
        self.mst = None
        adj_matrix = self.adj_mat

        """
        computes the minimum spanning tree of a graph using Prim's algorithm.

        parameters:
        - adjacency matrix: 2D list where adj_matrix[i][j] is the weight of the edge between i and j. 
            if there is no edge -> 0. 
        
        returns:
        - mst_matrix: 2D numpy array representing the adjacency matrix of the MST.
        """
        #number of nodes
        n = len(adj_matrix)
        #track visited nodes
        visited = set()
        #a min heap ensures that the smallest element is always at the top (index 0)
        min_heap = []
        #initialize mst matrix with zeros
        mst_matrix = np.zeros((n,n))

        #start from node 0 (arbitrary choice)
        visited.add(0)

        #adding all edges from node 0 to the heap
        for j in range(n):
            if adj_matrix[0, j] > 0:
                #when using heapq.heappush, heapq rearranged the heap to maintain smallest element property
                heapq.heappush(min_heap, (adj_matrix[0,j], 0, j)) #(weight, from, to)

        #ensures that not all nodes are in the mst and checks if the heap is not empty
        while len(visited) < n and min_heap:
                #removes the smallest edge and assigns it to weight, u, v
            weight, u, v = heapq.heappop(min_heap)
            
            if v in visited: 
                #continue skips the current interation of the loop and moves on to the next iteration
                #jumps immediately back to the condition of the while loop
                continue 

            #add the edge to the mst adjacency matrix 
            mst_matrix[u, v] = weight
            mst_matrix[v, u] = weight
                
            visited.add(v)

            #add new edges to the heap
            for j in range(n):
                if j not in visited and adj_matrix[v, j] > 0:
                    heapq.heappush(min_heap, (adj_matrix[v, j], v, j))
        
        self.mst = mst_matrix

        







