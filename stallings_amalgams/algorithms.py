import numpy as np
from .graph import Graph

def Step1(G1,G2,A,H):
    """
    Produces the wedge of loops labelled with the 
    generators of H, i.e., the graph Gamma1 of the article.
    
    Param:
        G1 Graph: Graph corresponding to the Cayley graph of the group G1
        G2 Graph: Same but for G2
        A list of tuples of list: tuples with one element in G1 and the other in G2
        H list of str: Each str corresponds to a word in the generators of G1 and G2

    Output:
        Graph object corresponding to the mentioned wedge.
    """
    Gam1=Graph(G1.labels|G2.labels)
    for gen in H:
        Gtemp=Graph(Gam1.labels)
        for i,label in enumerate(gen[:-1]):
            Gtemp.add_edge(i,-1,label)
        Gtemp.add_edge(Gtemp.n_verts-1,0,gen[-1]) 
        Gam1.Wedge(Gtemp, 0,0)
    return(Gam1)

def Step2(G,G1,G2,A):
    """
    Folds and cuts hairs of the output of Step1 given in Gam1
    
    Param:
        Gam1 Graph: Output of Step1
        
    Output:
        None. G is changed

    """
    G.fold()
    G.cut_hairs()
    return(G)

def Step3(G,G1,G2,A):
    """
    

    Parameters
    ----------
    G1 : TYPE
        DESCRIPTION.
    G2 : TYPE
        DESCRIPTION.
    A : Dict
        A should be a dict with keys elements of G1 and values the element of G2 that A glues.

    Returns
    -------
    None.

    """    
    comps1 = G.monochromatic_components(G1)
    comps2 = G.monochromatic_components(G2)
    
    # Wedge a copy of G1 on each G1 comp of G
    for C in comps1:
        for u in C:
            break
        G.Wedge(G1,u,G1.basepoint)
    
    # same for G2
    for C in comps2:
        for u in C:
            break
        G.Wedge(G2,u,G2.basepoint)
    G.fold()
    return(G)

    
def Step4(G,G1,G2,A):
    elements1 = G1.get_cosets()
    elements2 = G2.get_cosets()
    bichromatic = G.bichromatic_vertices(G1,G2)
    vs_to_glue=[]
    while bichromatic:
        v = bichromatic.pop()
        for a in A:
            u1 = G.read_word(v, elements1[a])
            u2 = G.read_word(v, elements2[A[a]])
            if u1 is not None and u2 is not None and u1 != u2:
                vs_to_glue.append([u1,u2])
                
    G.glue_pairs(vs_to_glue)
    G.fold()
    return(G)

    
def Step5(G, G1, G2, A):

    while True:

        components1 = G.monochromatic_components(G1)
        components2 = G.monochromatic_components(G2)

        components = [
            (C, G1) for C in components1
        ] + [
            (C, G2) for C in components2
        ]

        removed = False

        for C, group in components:

            redundant, _, K = G.is_redundant_component(
                C, G1, G2, A
            )

            if redundant:
                G.remove_redundant_component(C, group, G1, G2)
                removed = True
                break

        if not removed:
            break
        
    if not G.bichromatic_vertices(G1, G2) and (G.stabilizer(G.basepoint, G1) == {G1.basepoint} or G.stabilizer(G.basepoint, G2) == {G2.basepoint}):
        G.n_verts = 1
        G.basepoint = 0
        
        for label in G.mat:
            G.mat[label] = np.zeros((1, 1), dtype=int)
    return(G)

        
def Step6(G,G1,G2,A):
    # ----------------------------------------------------------
    # Case 1: basepoint is G1-monochromatic
    # ----------------------------------------------------------
    if G.is_monochromatic_vertex(G.basepoint, G1, G2):
        K = G.stabilizer(G.basepoint, G1)
        L = K & set(A.keys())
        G1elms=G1.get_cosets()
        G2elms=G2.get_cosets()
        if len(L)>1:
            rel = Graph.rel_cayley(G2,L)
            G.Wedge(rel, G.basepoint, rel.basepoint)
            for a in A:
                if a not in L:
                    G.glue([G.read_word(G.basepoint, G1elms[a]),
                               G.read_word(G.basepoint, G2elms[A[a]])]
                              )
    
    # ----------------------------------------------------------
    # Case 2: basepoint is G2-monochromatic
    # ----------------------------------------------------------

    elif G.is_monochromatic_vertex(G.basepoint, G2, G1):
        K = G.stabilizer(G.basepoint, G2)
        L = K & set(A.values())
        G1elms=G1.get_cosets()
        G2elms=G2.get_cosets()
        if len(L)>1:
            rel = Graph.rel_cayley(G1,L)
            G.Wedge(rel, G.basepoint, rel.basepoint)
            for a in A:
                if A[a] not in L:
                    G.glue([G.read_word(G.basepoint, G1elms[a]),
                               G.read_word(G.basepoint, G2elms[A[a]])]
                    )
    return(G)

    
    
def get_red_precover(G1,G2,A,H):
    
    Output=Step1(G1,G2,A,H)

    Output=Step2(Output,G1, G2, A)

    Output=Step3(Output,G1, G2, A)

    Output=Step4(Output,G1, G2, A)

    Output=Step5(Output,G1, G2, A)

    Output=Step6(Output,G1, G2, A)

    return(Output)