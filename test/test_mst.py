import pytest
import numpy as np
from mst import Graph
from sklearn.metrics import pairwise_distances

"""
as always using chatgpt to help with some issues
"""

def check_mst(adj_mat: np.ndarray, 
              mst: np.ndarray, 
              expected_weight: int, 
              allowed_error: float = 0.0001):
    """
    
    Helper function to check the correctness of the adjacency matrix encoding an MST.
    Note that because the MST of a graph is not guaranteed to be unique, we cannot 
    simply check for equality against a known MST of a graph. 

    Arguments:
        adj_mat: adjacency matrix of full graph
        mst: adjacency matrix of proposed minimum spanning tree
        expected_weight: weight of the minimum spanning tree of the full graph
        allowed_error: allowed difference between proposed MST weight and `expected_weight`

    TODO: Add additional assertions to ensure the correctness of your MST implementation. For
    example, how many edges should a minimum spanning tree have? Are minimum spanning trees
    always connected? What else can you think of?

    """

    def approx_equal(a, b):
        return abs(a - b) < allowed_error

    total = 0
    for i in range(mst.shape[0]):
        for j in range(i+1):
            total += mst[i, j]
    assert approx_equal(total, expected_weight), 'Proposed MST has incorrect expected weight'

    'assert all nodes are connected'
    visited = set()
    def dfs(node):
        visited.add(node)
        #the row is node, we extract the column index and the matrix entry
        #enumerating basically iterates over all elements in that row, returning index and value
        for neighbor, weight in enumerate(mst[node]):
            if weight > 0 and neighbor not in visited:
                dfs(neighbor)
    
    dfs(0) #starting from node 0 
    #asserting that length of visited nodes is number of total nodes
    assert len(visited) == len(mst)

    'assert mst is not cyclic by ensuring it has exactly n-1 edges'
    #np.sum counts all nonzero edges in the adjacency matrix, if it were cyclic it would have the
    #exact same amount of edges as nodes. 
    assert np.sum(mst != 0) // 2 == len(mst) - 1

    

def test_mst_small():
    """
    
    Unit test for the construction of a minimum spanning tree on a small graph.
    
    """
    file_path = './data/small.csv'
    g = Graph(file_path)
    g.construct_mst()
    check_mst(g.adj_mat, g.mst, 8)


def test_mst_single_cell_data():
    """
    
    Unit test for the construction of a minimum spanning tree using single cell
    data, taken from the Slingshot R package.

    https://bioconductor.org/packages/release/bioc/html/slingshot.html

    """
    file_path = './data/slingshot_example.txt'
    coords = np.loadtxt(file_path) # load coordinates of single cells in low-dimensional subspace
    dist_mat = pairwise_distances(coords) # compute pairwise distances to form graph
    g = Graph(dist_mat)
    g.construct_mst()
    check_mst(g.adj_mat, g.mst, 57.263561605571695)


def test_mst_student():
    """
    
    TODO: Write at least one unit test for MST construction.

    Ensuring MST handles negative edge weights correctly
    
    """
    adj_mat = np.array([
        [0, 3, 0, 7, 0],
        [3, 0, 1, 5, 0],
        [0, 1, 0, 2, 8],
        [7, 5, 2, 0, 6],
        [0, 0, 8, 6, 0]
    ])

    g = Graph(adj_mat)
    g.construct_mst()
    #manually checked the expected weight by drawing the spanning tree on paper -> 12
    check_mst(adj_mat, g.mst, 12)

    pass
