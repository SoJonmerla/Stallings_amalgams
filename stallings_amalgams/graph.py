from __future__ import annotations

import numpy as np
import copy
import networkx as nx

# ======================================================================================================================

class Graph:
    def __init__(
        self,
        labels: list[str] | set[str],
        G: np.ndarray | None = None,
        basepoint: int = 0,
    ) -> None:
        """
        Initialize a finite directed labelled graph.
         
        The graph is represented by a dictionary of adjacency matrices,
        with one matrix for each edge label. For every supplied positive
        label, the corresponding inverse label is added automatically.
         
        Parameters
        ----------
        labels : list[str] | set[str]
        Edge labels used in the graph. A corresponding inverse label
        of the form ``"label^-1"`` is added automatically when it is
        not already present.
         
        G : numpy.ndarray, optional
        A predefined labelled incidence matrix. The entry ``G[i][j]``
        should contain the labels of the directed edges from vertex
        ``i`` to vertex ``j``.
         
        If omitted, the graph is initialized with one vertex and no
        edges.
         
        basepoint : int, default=0
        Index of the distinguished basepoint vertex.
         
        Raises
        ------
        ValueError
        If the graph has no vertices, if the basepoint is outside the
        valid vertex range, or if no edge labels are supplied.
        TypeError
        If a label is not a string.
         
        Notes
        -----
        For each edge labelled ``x`` from vertex ``i`` to vertex ``j``,
        the graph also stores an inverse edge labelled ``x^-1`` from
        vertex ``j`` to vertex ``i``.
        """
        if not labels:
            raise ValueError("At least one edge label must be supplied.")
            
        if not all(isinstance(label, str) for label in labels):
            raise TypeError("Every edge label must be a string.")
            
        self.n_verts = 1 if G is None else len(G)  # Number of vertices
        
        if self.n_verts == 0:
            raise ValueError("The graph must have at least one vertex.")
        
        if not 0 <= basepoint < self.n_verts:
            raise ValueError(
                f"basepoint must be between 0 and {self.n_verts - 1}"
            )
            
        self.labels=set(labels)
        
        self.labels |= {f"{label}^-1" for label in self.labels if label[-3:]!="^-1"}
        
        self.mat = self.initialize(self.n_verts, self.labels)  # Incidence matrices for each label
        
        self.basepoint = basepoint
        
        if G is not None:  # Matrices already given
            self.predefined(G)
            
            
    def initialize(self, n_verts, labels):
        """
        Initialize the dictionary representing edge relations between vertices.
        
        Parameters
        ----------
        labels : set[str]
            Labels used for the edges.
        n_verts : int
            Number of vertices in the graph.
        
        Returns
        -------
        dict
            A dictionary whose keys are edge labels and whose associated values
            are NumPy matrices. A nonzero entry in row i and column j indicates
            that there is a directed edge from vertex i to vertex j. The label
            of that edge is the corresponding dictionary key.
        """
        d = {label: np.zeros([n_verts, n_verts], dtype=int) for label in labels}
        return d

    def predefined(self, g):
        """
        Register a predefined graph.
        
        Parameters
        ----------
        g : list
            An incidence matrix in which each entry is a list containing the
            labels of the directed edges from vertex i to vertex j.
        """
        index = np.argwhere(np.array(np.array(g, dtype=object), dtype=bool))
        for i, j,_ in index:
            for label in g[i][j]:
                self.add_edge(i, j, label)
        
        
    def __eq__(self, other):
        if not isinstance(other, Graph):
            return NotImplemented
    
        if self.labels != other.labels:
            return False
    
        if self.n_verts != other.n_verts:
            return False
    
        if self.basepoint != other.basepoint:
            return False
    
        for label in self.labels:
            if not np.array_equal(self.mat[label], other.mat[label]):
                return False
    
        return True
                
    def get_cosets(self) -> dict[int, list[str]]:
        """
        Returns elements of the graph self via dict of elements with value a word on self.labels to get there.
        -------
        None.

        """

    
        elements = {}
        seen = set()
    
        start = (self.basepoint,[]) # vertex and how to get there from basepoint
        stack = [start]
        seen.add(start[0])

        while stack:
            v = stack.pop()
            elements[v[0]]=v[1]

            for label in self.labels:
                w = self.next_vertex(v[0],label)
                if w is None:
                    continue
                if w not in seen:
                    seen.add(w)
                    stack.append((w,v[1]+[label]))

        return elements
    
    def add_edge(self,  vert_ini: int,  vert_end: int, label: str) -> None:
        """
        Add a directed edge, together with its inverse, between the specified
        initial and terminal vertices.
        
        Parameters
        ----------
        vert_ini : int
            Index of the initial vertex. If ``vert_ini == -1``, a new initial
            vertex is created.
        
        vert_end : int
            Index of the terminal vertex. If ``vert_end == -1``, a new terminal
            vertex is created.
        
        label : str
            Label assigned to the directed edge. If label not in self.labels,
            adds it to the set of labels.
        
        Notes
        -----
        If an edge labelled ``x`` is added from ``vert_ini`` to ``vert_end``,
        an inverse edge labelled ``x^-1`` is also added from ``vert_end`` to
        ``vert_ini``.
        """
        inverse=label[:-3] if label[-3:]=="^-1" else label + "^-1"
        n_new=0
        if vert_ini ==-1:
            n_new+=1
            self.n_verts+=1
            vert_ini=self.n_verts-1
        if vert_end==-1:
            n_new+=1
            self.n_verts+=1
            vert_end=self.n_verts-1        
        if n_new>0:
            for letter in self.mat:
                a=self.mat[letter]
                self.mat[letter]=np.full((self.n_verts,self.n_verts),0)
                self.mat[letter][:self.n_verts-n_new,:self.n_verts-n_new]=a
        if label not in self.labels:
            self.mat[label] = np.full((self.n_verts,self.n_verts),0)
            self.mat[inverse] = np.full((self.n_verts,self.n_verts),0)
            self.labels |= {label,inverse}
        
        self.mat[label][vert_ini,  vert_end]|= 1
        self.mat[inverse][vert_end,vert_ini]|=1
    def remove_vertex(self, vertex: int) -> None:
        """
        Remove a vertex and all edges incident to it.
    
        The row and column corresponding to the vertex are removed from
        every label-specific adjacency matrix. Vertices with larger indices
        are subsequently reindexed.
    
        Parameters
        ----------
        vertex : int
            Index of the vertex to remove.
    
        Raises
        ------
        IndexError
            If ``vertex`` is not a valid vertex index.
    
        ValueError
            If ``vertex`` is the graph's basepoint.
    
        Notes
        -----
        Removing a vertex decreases by one the indices of all vertices that
        originally had indices greater than ``vertex``. The basepoint index
        is adjusted accordingly.
        """
        if not 0 <= vertex < self.n_verts:
            raise IndexError(
                f"vertex must be between 0 and {self.n_verts - 1}"
            )
    
        if vertex == self.basepoint:
            raise ValueError("The basepoint vertex cannot be removed.")
    
        for label in self.mat:
            self.mat[label] = np.delete(
                self.mat[label],
                vertex,
                axis=0,
            )
            self.mat[label] = np.delete(
                self.mat[label],
                vertex,
                axis=1,
            )
    
        self.n_verts -= 1
    
        if self.basepoint > vertex:
            self.basepoint -= 1
    
    def remove_redundant_component(self, C, group, G1, G2):
        """
        Remove a redundant monochromatic component as in
        Definition 6.13 / the construction preceding Lemma 6.17.
        """
    
        C = set(C)
    
        # Special case: the whole graph is the component.
        if C == set(range(self.n_verts)):
    
            self.n_verts = 1
            self.basepoint = 0
    
            for label in self.mat:
                self.mat[label] = np.zeros((1, 1), dtype=int)
    
            return
    
        # Bichromatic vertices survive.
        VB = self.bichromatic_vertices(G1, G2, C)
    
        # These are the vertices that must actually be removed.
        VM = C - VB
    
        # Remove all edges belonging to C.
        for label in group.labels:
            matrix = self.mat[label]
    
            for u in C:
                for v in C:
                    matrix[u, v] = 0
    
        # Remove monochromatic vertices.
        for v in sorted(VM, reverse=True):
            self.remove_vertex(v)

    def plot_matrix(self):
        """
        Print incidence matrices in readable way
        """
        result = [[[] for _ in range(self.n_verts)] for _ in range(self.n_verts)]
        for label in self.mat:
            matrix = self.mat[label]
            if np.any(matrix):
                indexs = np.argwhere(matrix)
                for i, j in indexs:
                    result[i][j].append(label)
        if result==[[[]]]:
            print("The graph only has one vertex: \n")
        output=""
        max_width = max(len(str(cell))  for row in result for cell in row)
        for row in result:
            row_str = ' '.join(str(cell).center(max_width) for cell in row)
            output += f"{row_str}\n\n"
        
        print(output)


    def glue(self, vertices: list[int]) -> None:
        """
        Glue the specified vertices into a single vertex.
    
        The vertex with the smallest index survives. All other vertices in
        ``vertices`` are deleted, and their incident edges are transferred
        to the surviving vertex.
    
        Parameters
        ----------
        vertices : list[int]
            Indices of the vertices to glue.
    
        Raises
        ------
    
        IndexError
            If any supplied vertex index is invalid.
    
        Notes
        -----
        If the basepoint belongs to ``vertices``, the surviving vertex
        becomes the new basepoint. Otherwise, the basepoint index is adjusted
        to account for deleted vertices with smaller indices.
        """
        vertices = sorted(set(vertices))
    
        if len(vertices) <= 1:
            return
    
        if any(vertex < 0 or vertex >= self.n_verts for vertex in vertices):
            raise IndexError(
                f"Vertex indices must be between 0 and {self.n_verts - 1}."
            )
    
        survivor = vertices[0]
        removed_vertices = vertices[1:]
    
        # Work out the new basepoint before deleting and reindexing vertices.
        if self.basepoint in vertices:
            new_basepoint = survivor
        else:
            deleted_before_basepoint = sum(
                vertex < self.basepoint
                for vertex in removed_vertices
            )
            new_basepoint = self.basepoint - deleted_before_basepoint
    
        # Delete vertices in descending order to preserve the remaining
        # vertex indices during deletion.
        for label in self.mat:
            for vertex in reversed(removed_vertices):
                self.mat[label][survivor] |= self.mat[label][vertex]
                self.mat[label][:, survivor] |= self.mat[label][:, vertex]
    
                self.mat[label] = np.delete(
                    self.mat[label],
                    vertex,
                    axis=0,
                )
                self.mat[label] = np.delete(
                    self.mat[label],
                    vertex,
                    axis=1,
                )
    
        self.n_verts -= len(removed_vertices)
        self.basepoint = new_basepoint
        
        
    def glue_pairs(self,pairs: list[list[int,int]]) -> None:
        """
        Glue all pairs ov vertices in pairs. Ensures reindicing. Pairs must be lists of list
        """
        while pairs:
            pair = pairs.pop()
            self.glue(pair)
            survivor = min(pair[0],pair[1])
            removed = max(pair[0],pair[1])
            for i in range(len(pairs)):
                if pairs[i][0] == removed:
                    pairs[i][0]=survivor
                    
                if pairs[i][0] > removed:
                    pairs[i][0] -=1
                if pairs[i][1] == removed:
                    pairs[0]=survivor
                    
                if pairs[i][1] > removed:
                    pairs[i][1] -=1
     
    
    def fold(self) -> None:
        """
        Folds the graph of self until no further folding are possible.
        
        """
        while True:
            fold_found=False
            for label in self.mat:
                if not label.endswith("^-1"):
                    for vertex in range(self.n_verts):
                        R=np.nonzero(self.mat[label][vertex])[0]
                        if len(R)>1:
                            self.glue(R)
                            fold_found=True
                            break
                            
                            
                        R=np.nonzero(self.mat[label][:,vertex])[0]
                        if len(R)>1:
                            self.glue(R)
                            fold_found=True
                            break
                            
            if not fold_found:
                return
                    
        
    def cut_hairs(self) -> None:
        """
        Remove all hairs from the graph.
    
        A hair is a non-basepoint vertex of valence at most one. When a hair
        is removed, another vertex may become a hair, so the procedure is
        repeated until no hairs remain.
    
        The basepoint is never removed.
    
        Notes
        -----
        The graph stores inverse edges explicitly. Therefore, examining the
        outgoing edges for every label also detects edges that are incoming
        with respect to the corresponding positive label.
    
        This method modifies the graph in place.
        """
        while True:
            hair_found = False
    
            for vertex in range(self.n_verts):
                if vertex == self.basepoint:
                    continue
    
                valence = 0
    
                for label in self.mat:
                    valence += np.count_nonzero(
                        self.mat[label][vertex]
                    )
    
                    if valence > 1:
                        break
    
                if valence <= 1:
                    self.remove_vertex(vertex)
                    hair_found = True
                    break
    
            if not hair_found:
                return
            
            
    def is_monochromatic_vertex(self, vertex, G1, G2):
        """
        True exactly when every edge incident to `vertex` is labelled
        by a generator from `G1` (or its inverse).
    
        A vertex with no incident edges is treated as non-monochromatic.
        """
        generators = set(G1.labels)
        other_generators = set(G2.labels)
    
        has_edge = False
    
        for label in self.labels:
            # Check both incoming and outgoing edges.
            incident = (
                np.any(self.mat[label][vertex, :]) or
                np.any(self.mat[label][:, vertex])
            )
    
            if not incident:
                continue
    
            has_edge = True
    
            # label belongs to the other factor
            if label in other_generators:
                return False

        return has_edge

    def monochromatic_vertices(
        self,
        G1: Graph,
        G2: Graph,
    ) -> set[int]:
        """
        Return the vertices that are monochromatic with respect to ``G1``.
    
        A vertex is considered ``G1``-monochromatic if every edge incident
        to it is labelled by a generator of ``G1`` or by the inverse of such
        a generator. An isolated vertex is not considered monochromatic.
    
        Parameters
        ----------
        G1 : Graph
            Cayley graph of the factor whose labels determine whether a
            vertex is monochromatic.
    
        G2 : Graph
            Cayley graph of the other factor.
    
        Returns
        -------
        set[int]
            Indices of the vertices that are monochromatic with respect to
            ``G1``.
        """
        return {
            v for v in range(self.n_verts)
            if self.is_monochromatic_vertex(v,  G1, G2)
        }
    def bichromatic_vertices(self,G1,G2,component = None):
        """
        Parameters
        ----------
        G1 : graph
        G2 : graph

        Returns
        -------
        set
            Bichromatic vertices of self contained in component. 
            If no component given, assumes component is all graph

        """
        if component is None:
            component = {i for i in range(self.n_verts)}
        return {i for i in component if i not in self.monochromatic_vertices(G1, G2)|self.monochromatic_vertices(G2, G1)}
            
    def monochromatic_components(self, G: Graph) -> list[dict]:
        """
        Return the connected X-components of the graph. Where X is the gens of G
    
        generators:
            set/list of positive generators belonging to one factor.
    
        Returns:
            list of sets of vertices. Components consisting of a single
            isolated vertex are not returned.
        """
        generators = set(G.labels)
    
        # We only need the positive labels. The inverse edges are already
        # stored in self.mat.
        adjacency = {v: set() for v in range(self.n_verts)}
    
        for label in generators:
            if label not in self.mat:
                continue
    
            matrix = self.mat[label]
    
            for u, v in np.argwhere(matrix):
                adjacency[u].add(v)
                adjacency[v].add(u)
    
        components = []
        seen = set()
    
        for start in range(self.n_verts):
            if start in seen or not adjacency[start]:
                continue
    
            component = set()
            stack = [start]
            seen.add(start)
    
            while stack:
                v = stack.pop()
                component.add(v)
    
                for w in adjacency[v]:
                    if w not in seen:
                        seen.add(w)
                        stack.append(w)
    
            components.append(component)

        return components
    
    def is_redundant_component(self, C, G1, G2, A):
        """
        Determine whether the monochromatic component C is redundant
        according to Definition 6.13.
        
        Returns:
            (True/False, factor, K)
        """
    
        components1 = self.monochromatic_components(G1)
        components2 = self.monochromatic_components(G2)
    
    
        # Determine which factor C belongs to.
        if C in components1:
            group = G1
            A_elements = set(A.keys())
        elif C in components2:
            group = G2
            A_elements = set(A.values())
        else:
            raise ValueError("C is not a monochromatic component.")
    
    
        # ------------------------------------------------------------
        # Case 2: at least two monochromatic components.
        # ------------------------------------------------------------
    
        VB = self.bichromatic_vertices(G1, G2, C)
    
        # Definition 6.13(2) starts with theta in VB(C).
        if not VB:
            return False, group, None
    
        theta = next(iter(VB))
    
        K = self.stabilizer(theta, group)
    
        # (i) K <= A
        if not K <= A_elements:
            return False, group, K
    
        # (ii) |VB(C)| = [A : K]
        index = len(A_elements) // len(K)
    
        if len(VB) != index:
            return False, group, K
    
        # (iii) basepoint condition
        if self.basepoint not in C:
            return True, group, K
    
        if self.basepoint in VB and K == {group.basepoint}:
            return True, group, K
    
        return False, group, K
        
    def next_vertex(self, vertex: int, label: str) -> int:
        """
        Follow one directed edge labelled `label`.
    
        Returns None if no such edge exists.
        Raises ValueError if the graph is not well-labelled.
        """
        if label not in self.mat:
            raise ValueError(f"Unknown label: {label}")
    
        targets = np.nonzero(self.mat[label][vertex])[0]
    
        if len(targets) == 0:
            return None
    
        if len(targets) > 1:
            raise ValueError(
                f"Graph is not well-labelled: vertex {vertex} "
                f"has multiple outgoing {label}-edges."
                "Consider folding."
            )
    
        return int(targets[0])
    
    def read_word(self, start: int, word: list[str]) -> int | None:
        """
        Read a word from `start`. Word must be a list of gen and inverses.
    
        Returns:
            terminal vertex, or None if the path is not defined.
        """
        vertex = start
    
        for label in word:
            vertex = self.next_vertex(vertex, label)
    
            if vertex is None:
                return None
    
        return vertex
    
    def stabilizer(self, vertex: int,group: Graph) -> set:
        """
        Return the subgroup of `group` represented by loops at `vertex` in self. 
    
        Assumes group is a graph object representing its cayley graph' 
        Furhthermore, assumes self is already a cover of group!!
        """
        K = set()
        elements = group.get_cosets()

        for g,word in elements.items():    
            endpoint = self.read_word(vertex, word)
    
            if endpoint == vertex:
                K.add(g)
    
        return K
    
    def rel_cayley(G,H: dict[int]):
        """
        Get rel cayley graph of H in G. H given by subset of elements of G

        Returns
        Graph object corresponding to Cay(G,H)
        
        This works because cayley graph with identifications is still G-based. 
        So G-based + (stabilizer = H) + saturated  implies isomorphic to rel cayley graph.
        G-based comes from the fact that if a path p in the cayley graph after identification
        labels 1, then we can find a path in 
        original cayley graph with same label that has to be closed. So folding 
        closes the path p.

        """
        D = copy.deepcopy(G)
        H = set(H) | {G.basepoint}
        D.glue(list(H))
        D.fold()
        return D
    
    def Wedge(self,G2,u,v):
        """
        Does the disjoint union of the two graphs self, G2 and then 
        identifies the vertex v of self with the first vertex of G2.
        So basically, the wedge along the vertices v and 0 of self and 
        G2 respectively.  Assumes gens of self and G2 are different
        Param:
        self,G2 Graph: Graph objects
            v int: Index of vertex of self 
        Output:
            The Graph object corresponding to the mentioned wedge
        """
        n1=self.n_verts
        self.n_verts+= G2.n_verts
        N = self.n_verts
        self.labels|=G2.labels
        for label in self.labels:
            zero_matrix = np.full((N,N),0)
            if label in self.mat:
                zero_matrix[:n1,:n1] = self.mat[label]
            if label in G2.mat:
                zero_matrix[n1:,n1:] = G2.mat[label]
            self.mat[label]=zero_matrix
        self.glue([u,n1+v])

    def isomorphic_cayley_graphs(G1, G2):
        """
        Check whether two Cayley graphs are identical up to a renumbering
        of vertices.
    
        Parameters
        ----------
        G1, G2 : Graph
    
        Returns
        -------
        bool
        """
        if G1.labels != G2.labels:
            return False
    
        H1 = nx.MultiDiGraph()
        H2 = nx.MultiDiGraph()
    
        H1.add_nodes_from(range(G1.n_verts))
        H2.add_nodes_from(range(G2.n_verts))
    
        for label in G1.labels:
            M = G1.mat[label]
    
            for i, j in zip(*M.nonzero()):
                H1.add_edge(i, j, label=label)
    
        for label in G2.labels:
            M = G2.mat[label]
    
            for i, j in zip(*M.nonzero()):
                H2.add_edge(i, j, label=label)
    
        edge_match = nx.algorithms.isomorphism.categorical_multiedge_match(
            "label", None
        )
    
        return nx.is_isomorphic(H1, H2, edge_match=edge_match)