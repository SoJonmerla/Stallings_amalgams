import numpy as np
from .graph import Graph

def Step1(
    G1: Graph,
    G2: Graph,
    A: dict[int, int],
    H: list[list[str]]) -> Graph:
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
    G=Graph(G1.labels|G2.labels)
    for gen in H:
        Gtemp=Graph(G.labels)
        for i,label in enumerate(gen[:-1]):
            Gtemp.add_edge(i,-1,label)
        Gtemp.add_edge(Gtemp.n_verts-1,0,gen[-1]) 
        G.Wedge(Gtemp, 0,0)
    return(G)

def Step2(
    G: Graph,
    G1: Graph,
    G2: Graph,
    A: dict[int, int]) -> Graph:
    
    """
    Folds and cuts hairs of the output of Step1 given in G
    
    Param:
        G Graph: Output of Step1
        
    Output:
        G Graph

    """
    G.fold()
    G.cut_hairs()
    return(G)

def Step3(
    G: Graph,
    G1: Graph,
    G2: Graph,
    A: dict[int, int]) -> Graph:
    """
    Complete each monochromatic component with a copy of the
    corresponding factor's Cayley graph.

    For every G1-monochromatic component of G, attach a copy of G1 by
    identifying an arbitrary vertex of the component with the basepoint
    of G1. The analogous operation is then performed for every
    G2-monochromatic component.

    Finally, fold the resulting graph until it is well-labelled.

    Parameters
    ----------
    G : Graph
    The graph modified by this step of the reduced-precover
    construction.

    G1 : Graph
    Cayley graph of the first finite factor.

    G2 : Graph
    Cayley graph of the second finite factor.

    A : dict[int, int]
    Dictionary representing the identification of the amalgamated
    subgroup. Its keys represent elements of G1, and its values
    represent the corresponding elements of G2.

    This parameter is included for consistency with the complete
    algorithm, although it is not used directly in Step 3.

    Returns
    -------
    The graph G that is also modified in place.
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

    
def Step4(
    G: Graph,
    G1: Graph,
    G2: Graph,
    A: dict[int, int]) -> Graph:
    """
    Perform Step 4 of the reduced-precover construction.

    At every bichromatic vertex, compare the endpoints obtained by
    reading words representing corresponding elements of the
    amalgamated subgroup in the two factors.

    If the two words are both readable but terminate at different
    vertices, those terminal vertices are marked for identification.
    After all bichromatic vertices have been examined, the marked
    vertex pairs are glued and the resulting graph is folded.

    Parameters
    ----------
    G : Graph
    The graph on which Step 4 is performed. It is modified in place.

    G1 : Graph
    Cayley graph of the first finite factor.

    G2 : Graph
    Cayley graph of the second finite factor.

    A : dict[int, int]
    Dictionary representing the identification of the amalgamated
    subgroup. Each key represents an element of G1, and the
    corresponding value represents the identified element of G2.

    Returns
    -------
    Graph
    The modified and folded graph ``G``.
    """
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

    
def Step5(
    G: Graph,
    G1: Graph,
    G2: Graph,
    A: dict[int, int]) -> Graph:
    """
    Perform Step 5 of the reduced-precover construction.

    Redundant monochromatic components are removed one at a time.
    After each removal, the monochromatic components are recomputed
    because deleting a component may change the structure of the graph.

    If the resulting graph has no bichromatic vertices and the
    stabilizer of the basepoint is trivial in at least one factor,
    the graph is replaced by the trivial one-vertex graph.

    Parameters
    ----------
    G : Graph
    The graph on which Step 5 is performed. It is modified in place.

    G1 : Graph
    Cayley graph of the first finite factor.

    G2 : Graph
    Cayley graph of the second finite factor.

    A : dict[int, int]
    Dictionary representing the identification of the amalgamated
    subgroup. Each key represents an element of G1, and its value
    represents the corresponding element of G2.

    Returns
    -------
    Graph
    The modified graph ``G``.
    """
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
        
    if not G.bichromatic_vertices(G1, G2) and (
        (G.stabilizer(G.basepoint, G1) == {G1.basepoint}
         and G.is_monochromatic_vertex(G.basepoint, G1, G2))
        or
        (G.stabilizer(G.basepoint, G2) == {G2.basepoint}
         and G.is_monochromatic_vertex(G.basepoint, G2, G1))
    ):        
        G.n_verts = 1
        G.basepoint = 0
        
        for label in G.mat:
            G.mat[label] = np.zeros((1, 1), dtype=int)
    return(G)

        
def Step6(
    G: Graph,
    G1: Graph,
    G2: Graph,
    A: dict[int, int]) -> Graph:
    """
    Perform Step 6 of the reduced-precover construction.
    
    If the basepoint is monochromatic with respect to one factor, compute
    its stabilizer in that factor and intersect the stabilizer with the
    amalgamated subgroup. When this intersection is nontrivial, attach the
    appropriate relative Cayley graph of the other factor at the basepoint.
    
    The remaining corresponding elements of the amalgamated subgroup are
    then identified by gluing the endpoints of their representative words.
    
    Parameters
    ----------
    G : Graph
        The graph on which Step 6 is performed. It is modified in place.
    
    G1 : Graph
        Cayley graph of the first finite factor.
    
    G2 : Graph
        Cayley graph of the second finite factor.
    
    A : dict[int, int]
        Dictionary representing the identification of the amalgamated
        subgroup. Each key represents an element of G1, and its corresponding
        value represents the identified element of G2.
    
    Returns
    -------
    G
        The graph ``G`` tha is also modified in place.
    """
    # ----------------------------------------------------------
    # Case 1: basepoint is G1-monochromatic
    # ----------------------------------------------------------
    if G.is_monochromatic_vertex(G.basepoint, G1, G2):
        K = G.stabilizer(G.basepoint, G1)
        L = K & set(A.keys())
        G1elms=G1.get_cosets()
        G2elms=G2.get_cosets()
        if len(L)>1:
            L = {A[i] for i in L}
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
            L = {i for i in A if A[i] in L}
            rel = Graph.rel_cayley(G1,L)
            G.Wedge(rel, G.basepoint, rel.basepoint)
            for a in A:
                if A[a] not in L:
                    G.glue([G.read_word(G.basepoint, G1elms[a]),
                               G.read_word(G.basepoint, G2elms[A[a]])]
                    )
    return(G)

    
    
def get_red_precover(G1: Graph,G2: Graph,A: dict[int,int],H: list[list[str]]) -> Graph:
    
    Output=Step1(G1,G2,A,H)

    Output=Step2(Output,G1, G2, A)

    Output=Step3(Output,G1, G2, A)

    Output=Step4(Output,G1, G2, A)

    Output=Step5(Output,G1, G2, A)

    Output=Step6(Output,G1, G2, A)

    return(Output)


def is_subgroup_member(
    G1: Graph,
    G2: Graph,
    A: dict[int, int],
    subgroup_generators: list[list[str]],
    word: list[str],
) -> bool:
    """
    Determine whether a word represents an element of a subgroup.

    The ambient group is the amalgamated free product G1 *_A G2. The
    subgroup is generated by ``subgroup_generators``. Membership is
    decided by constructing the reduced precover of the subgroup and
    checking whether ``word`` labels a closed path at its basepoint.

    Parameters
    ----------
    G1 : Graph
        Cayley graph of the first finite factor.

    G2 : Graph
        Cayley graph of the second finite factor.

    A : dict[int, int]
        Identification of the amalgamated subgroup in the two factors.
        Each key represents an element of G1, and the corresponding value
        represents the identified element of G2.

    subgroup_generators : list[list[str]]
        Words generating the subgroup whose membership problem is being
        considered. Each word is represented as a list of edge labels.

    word : list[str]
        Word to test for membership, represented as a list of edge labels.

    Returns
    -------
    bool
        ``True`` if ``word`` represents an element of the subgroup
        generated by ``subgroup_generators``; otherwise ``False``.
    """
    reduced_precover = get_red_precover(
        G1,
        G2,
        A,
        subgroup_generators,
    )

    endpoint = reduced_precover.read_word(
        reduced_precover.basepoint,
        word,
    )

    return endpoint == reduced_precover.basepoint

def subgroup_index(
    G1: Graph,
    G2: Graph,
    A: dict[int, int],
    subgroup_generators: list[list[str]],
) -> int | None:
    """
    Return the index of the generated subgroup, or ``None`` if the index
    is infinite.

    The ambient group is the amalgamated free product G1 *_A G2, and the
    subgroup is generated by ``subgroup_generators``.

    The function constructs the reduced precover representing the subgroup.
    The subgroup has finite index precisely when this reduced precover is a
    complete coset graph: from every vertex, every label of the ambient
    generating set must be readable.

    Parameters
    ----------
    G1 : Graph
    Cayley graph of the first finite factor.

    G2 : Graph
    Cayley graph of the second finite factor.

    A : dict[int, int]
    Identification of the amalgamated subgroup in the two factors.
    Each key represents an element of G1, and the corresponding value
    represents the identified element of G2.

    subgroup_generators : list[list[str]]
    Words generating the subgroup. Each word is represented as a list
    of edge labels.

    Returns
    -------
    int or None
    The index of the generated subgroup if the index is finite;
    otherwise ``None``.

    """
    reduced_precover = get_red_precover(G1,G2,A,subgroup_generators)
    if all(
    reduced_precover.mat[label].any(axis=1).all()
    for label in reduced_precover.labels
    ): # I.e. graph is saturated/complete
        return reduced_precover.n_verts

    return None

def is_free(
    G1: Graph,
    G2: Graph,
    A: dict[int, int],
    subgroup_generators: list[list[str]],
) -> bool:
    """
    Determine whether a finitely generated subgroup of G1 *_A G2 is free.

    Parameters
    ----------
    G1 : Graph
    Cayley graph of the first finite factor.

    G2 : Graph
    Cayley graph of the second finite factor.

    A : dict[int, int]
    Identification of the amalgamated subgroup in the two factors.

    subgroup_generators : list[list[str]]
    Words generating the subgroup.

    Returns
    -------
    bool
    ``True`` if the generated subgroup is free; otherwise ``False``.
    """
    reduced_precover = get_red_precover(G1,G2,A,subgroup_generators)
    G1comps = reduced_precover.monochromatic_components(G1)
    G2comps = reduced_precover.monochromatic_components(G2)
    components = [
        (C, G1) for C in G1comps
    ] + [
        (C, G2) for C in G2comps
    ]
    return all(
        reduced_precover.stabilizer(next(iter(component)),factor)=={factor.basepoint}
        for component,factor in components
        )


def product_graph(G1: Graph, G2: Graph) -> Graph:
    """
    Construct the synchronous labelled product of two graphs.

    The vertices of the product are ordered pairs ``(v1, v2)``, where
    ``v1`` is a vertex of ``G1`` and ``v2`` is a vertex of ``G2``. The
    pair ``(v1, v2)`` is represented internally by the integer

    v1 * G2.n_verts + v2.

    There is an edge

    (v1, v2) --label--> (u1, u2)

    precisely when both coordinate graphs contain the corresponding
    labelled edges

    v1 --label--> u1

    and

    v2 --label--> u2.

    Parameters
    ----------
    G1 : Graph
    First labelled graph.

    G2 : Graph
    Second labelled graph.

    Returns
    -------
    Graph
    The synchronous labelled product of ``G1`` and ``G2``.

    Notes
    -----
    This implementation assumes that ``add_edge`` automatically adds the
    corresponding inverse-labelled edge. Therefore, only one label from
    each inverse pair is processed explicitly.
    """
    n1 = G1.n_verts
    n2 = G2.n_verts
    N = n1*n2
    common_labels = G1.labels & G2.labels
    prod = Graph(common_labels, np.array([[ [] for i in range(N) ] for j in range(N)]))
    for v1 in range(n1):
        for u1 in range(n2): # add edges (v1,u1) -> (v2, u2) when appropiate
            v = v1*n2 + u1
            for label in common_labels:
                if label.endswith("^-1"):
                    continue
                
                link1 = np.nonzero(G1.mat[label][v1])[0]
                link2 = np.nonzero(G2.mat[label][u1])[0]
                for v2 in link1:
                    for u2 in link2:
                        prod.add_edge(v, v2*n2 + u2,label)
    prod.basepoint = (
    G1.basepoint * n2 + G2.basepoint
    )

    return prod

    

                                   
            


## Implement normality, malnormality, presentation??
    
    
    
    
    
    