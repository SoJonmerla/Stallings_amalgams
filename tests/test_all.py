# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 15:57:20 2026

@author: jfm1v22
"""

from stallings_amalgams import Graph, get_red_precover
from stallings_amalgams.groups import create_D2n,create_Cn
from stallings_amalgams.visualization import visualize
import numpy as np

def test_empty_graph():
    G = Graph(["a", "b"])

    assert G.n_verts == 1
    assert G.basepoint == 0
    assert G.labels == {"a", "b","a^-1","b^-1"}
    
def test_add_edge():
    G = Graph(["a"])

    G.add_edge(0, 0, "a")

    assert G.mat["a"][0, 0] == 1
    assert G.mat["a^-1"][0, 0] == 1
    
def test_remove_vertex():
    G = Graph(["a"])

    G.add_edge(0, -1, "a")

    G.remove_vertex(1)

    assert G.n_verts == 1
    
def test_read_word():
    G = Graph(["a"])

    G.add_edge(0, -1, "a")

    assert G.read_word(0, ["a"]) == 1
    
def test_read_word_failure():
    G = Graph(["a"])

    assert G.read_word(0, ["a"]) is None
    
def test_read_empty_word():
    G = Graph(["a"])

    assert G.read_word(0, []) == 0
    
def test_D12_relations():
    D = create_D2n(6)

    assert D.read_word(0,["a" for i in range(6)]) == 0
    assert D.read_word(0,["b","b"]) == 0
    assert D.read_word(0,["a^-1"]) == D.read_word(0,["b","a","b"])
    
def test_precover_of_factor():
    """
    Normal core of relative cayley graph of the amalgam by one of the factor 1 must be
    rel cayley graph of factor 1 by factor 1 wedged by relative cayley graph of factor 2 
    by amalgamated subgroup.

    Returns
    -------
    None
    """
    A = {0:0,6: 2}

    H1 = [
        ["a"],
        ["b"]
    ]
    H2 = [
        ["c"]
    ]
    
    
    
    D12 = create_D2n(6)

    C4 = create_Cn(4)

    # visualize(Graph.rel_cayley(C4, A.values()))
    # visualize(get_red_precover(D12, C4, A, H1))
    goal1=Graph.rel_cayley(C4, A.values())
    goal1.add_edge(0, 0, "a")
    goal1.add_edge(0, 0, "b")
    
    goal2 = Graph.rel_cayley(D12, A.keys())
    goal2.add_edge(0, 0, "c")
    assert get_red_precover(D12, C4, A, H1) == goal1
    assert get_red_precover(D12, C4, A, H2) == goal2
    
def test_precover_of_all():
    D24 = create_D2n(12)

    C8 = create_Cn(8)
    
    A = {0:0,12:4}
    
    H = ["a","b","c"]
    
    goal = Graph(labels = ["a","b","c"],G = [[["a","b","c"]]])
    
    assert get_red_precover(D24, C8, A, H) == goal
    
def test_precover_of_triv():
    D24 = create_D2n(12)

    C8 = create_Cn(8)
    
    A = {0:0,12:4}
    
    H = []
    
    goal = Graph(labels = ["a","b","c"])
    
    assert get_red_precover(D24, C8, A, H) == goal
    

    