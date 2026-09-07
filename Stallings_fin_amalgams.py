import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import my_networkx as my_nx
import copy
import time

# ======================================================================================================================

class Graph:

    def __init__(self,  labels, G = np.array([[0]]),basepoint=0):
        self.n_verts = len(G)  # Número de nodos
        self.labels=set(labels)
        self.labels |= {f"{label}^-1" for label in self.labels if label[-3:]!="^-1"}
        self.mat = self.initialize(self.n_verts, self.labels)  # Representación del grafo
        self.basepoint = basepoint
        if G.any():  # Grafi predefinido en matriz de incidencia
            self.predefined(G)
    def initialize(self, n_verts, labels):
        """
        Esta función inicializa el diccionario que representa relaciones de aristas entre nodos
        :param:
            labels str: letras del mismo valor
            n_ verts int: número de nodos que hay en el grafo
        :return:
            d dict: diccionario donde las claves son letras y los valores asociados son matrices numpy. Si el valor de
                la fila i y columna j es distinto de cero, indica que existe una arista dirigida desde el vertice i al
                vertice j y la letra de dicha arista correspondrá con la clave en el diccionario.
                OJO: el valor 1 indica a y el valor 1/2 indica a^-1???? BORRAR?
        """
        d = {label: np.zeros([n_verts, n_verts], dtype=int) for label in labels}
        return d

    def predefined(self, g):
        """
        Registrar el grafo predefinido
        :param:
            g list: matriz incidencia donde cada elemento es una lista que indica los labels de vertice i a vertice j
        """
        index = np.argwhere(np.array(np.array(g, dtype=object), dtype=bool))
        for i, j in index:
            for label in g[i][j]:
                self.add_edge(i, j, label)
                
    def get_cosets(self):
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
    
    def add_edge(self,  vert_ini,  vert_end, label):
        """
        Añadir la arista direccional (y su inversa) según el vertice inicial y final indicado con su label
        :param:
             vert_ini int: Indice del vertice saliente. 
             vert_ini=-1: Arista con origen un nuevo vertice.
             vert_end int: Indice del vertice entrante.
             vert_end=-1: End vertex of edge is a new vertex to be defined.
       

            label str: la letra que tiene dicha arista
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
                
        self.mat[label][vert_ini,  vert_end]|= 1
        self.mat[inverse][vert_end,vert_ini]|=1
    def remove_vertex(self, vertex):
        for label in self.mat:
            self.mat[label] = np.delete(self.mat[label], vertex, axis=0)
            self.mat[label] = np.delete(self.mat[label], vertex, axis=1)
    
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
        A partir del grafo actual, imprimir la matriz de incidencia
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


    def glue(self,  L):
        """
        Pegar los vertices de self indexados por los int de L, vertice identificado tendr'a indice el menor de L'
        :param:
             L list of int: Vertices of self to glue.            
        """
        L.sort()
        n=L[0]
        L1=list(L[1:])
        L1.reverse()
        for label in self.mat:
            for i in L1:
                self.mat[label][n]|= self.mat[label][ i]
                self.mat[label][:,  n]|= self.mat[label][:,  i]
                self.mat[label] = np.delete(self.mat[label],  i, axis=0)
                self.mat[label] = np.delete(self.mat[label],  i, axis=1)
        self.n_verts=len(self.mat[label])
        
        
    def glue_pairs(self,pairs):
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
     
    
    def fold(self):
        """
        Folds the graph of self until no further folding are possible.
        
        """
        for label in self.mat:
            if not label.endswith("^-1"):
                for vertex in range(len(self.mat[label])):
                    R=np.nonzero(self.mat[label][vertex])[0]
                    if len(R)>1:
                        self.glue(R)
                        self.fold()
                        return
                    R=np.nonzero(self.mat[label][:,vertex])[0]
                    if len(R)>1:
                        self.glue(R)
                        self.fold()
                        return
                
        
    def cut_hairs(self):
        for vertex in range(self.n_verts):
            if vertex == self.basepoint:
                continue
            valence=0
            for label in self.mat:
                valence+=len(np.nonzero(self.mat[label][vertex])[0])
                if valence>1:
                    break
            else: # This will be executed if the if statement in the previous loop is never satisfied.
                self.remove_vertex(vertex)
                return self.cut_hairs()
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

    def monochromatic_vertices(self, G1, G2):
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
            
    def monochromatic_components(self, G):
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
        
    def next_vertex(self, vertex, label):
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
    
    def read_word(self, start, word):
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
    
    def stabilizer(self, vertex,group):
        """
        Return the subgroup of `group` represented by loops at `vertex`. 
    
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
    
    def rel_cayley(G,H):
        """
        Get rel cayley graph of H in G. H given by subset of elements of G

        Returns
        Graph object corresponding to Cay(G,H)
        
        This works because cayley graph with identifications is still G-based. 
        So G-based + stablizer =H + saturated  implies isomorphic to rel cayley graph.
        G-based comes from the fact that if a path p in the cayley graph after identification
        labels 1, then we can find a path in 
        original cayley graph with same label that has to be closed. So folding 
        closes the the path p.

        """
        D = copy.deepcopy(G)
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
        
        # Output=Graph(list(self.labels|G2.labels))
        # n1, n2= self.n_verts, G2.n_verts
        # N = n1 + n2
        # Output.n_verts = N
        # for label in Output.labels:
        #     zero_matrix = np.full((N,N),0)
        #     if label in self.mat:
        #         zero_matrix[:n1,:n1] = self.mat[label]
        #     if label in G2.mat:
        #         zero_matrix[n1:,n1:] = G2.mat[label]
        #     Output.mat[label]=zero_matrix

        # Output.glue([u,n1+v])
        # # Output.fold()
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
            Gam1=Wedge(Gam1, Gtemp, 0,0)
        return(Gam1)

    def Step2(self,G1,G2,A):
        """
        Folds and cuts hairs of the output of Step1 given in Gam1
        
        Param:
            Gam1 Graph: Output of Step1
            
        Output:
            None. Self is changed

        """
        self.fold()
        self.cut_hairs()
    
    
    def Step3(self,G1,G2,A):
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
        comps1 = self.monochromatic_components(G1)
        comps2 = self.monochromatic_components(G2)
        
        # Wedge a copy of G1 on each G1 comp of self
        for C in comps1:
            for u in C:
                break
            self.Wedge(G1,u,G1.basepoint)
        
        # same for G2
        for C in comps2:
            for u in C:
                break
            self.Wedge(G2,u,G2.basepoint)
        self.fold()
        
    def Step4(self,G1,G2,A):
        elements1 = G1.get_cosets()
        elements2 = G2.get_cosets()
        bichromatic = self.bichromatic_vertices(G1,G2)
        vs_to_glue=[]
        while bichromatic:
            v = bichromatic.pop()
            for a in A:
                u1 = self.read_word(v, elements1[a])
                u2 = self.read_word(v, elements2[A[a]])
                if u1 is not None and u2 is not None and u1 != u2:
                    vs_to_glue.append([u1,u2])
                    
        self.glue_pairs(vs_to_glue)
        self.fold()
        
    def Step5(self, G1, G2, A):

        while True:
    
            components1 = self.monochromatic_components(G1)
            components2 = self.monochromatic_components(G2)
    
            components = [
                (C, G1) for C in components1
            ] + [
                (C, G2) for C in components2
            ]
    
            removed = False
    
            for C, group in components:
    
                redundant, _, K = self.is_redundant_component(
                    C, G1, G2, A
                )
    
                if redundant:
                    self.remove_redundant_component(C, group, G1, G2)
                    removed = True
                    break
    
            if not removed:
                break
            
        if not self.bichromatic_vertices(G1, G2) and (self.stabilizer(self.basepoint, G1) == {G1.basepoint} or self.stabilizer(self.basepoint, G2) == {G2.basepoint}):
            self.n_verts = 1
            self.basepoint = 0
            
            for label in self.mat:
                self.mat[label] = np.zeros((1, 1), dtype=int)
            
            
    def Step6(self,G1,G2,A):
        # ----------------------------------------------------------
        # Case 1: basepoint is G1-monochromatic
        # ----------------------------------------------------------
        if self.is_monochromatic_vertex(self.basepoint, G1, G2):
            K = self.stabilizer(self.basepoint, G1)
            L = K & set(A.keys())
            G1elms=G1.get_cosets()
            G2elms=G2.get_cosets()
            if len(L)>1:
                rel = Graph.rel_cayley(G2,L)
                self.Wedge(rel, self.basepoint, rel.basepoint)
                for a in A:
                    if a not in L:
                        self.glue([self.read_word(self.basepoint, G1elms[a]),
                                   self.read_word(self.basepoint, G2elms[A[a]])]
                                  )
        
        # ----------------------------------------------------------
        # Case 2: basepoint is G2-monochromatic
        # ----------------------------------------------------------

        elif self.is_monochromatic_vertex(self.basepoint, G2, G1):
            K = self.stabilizer(self.basepoint, G2)
            L = K & set(A.values())
            G1elms=G1.get_cosets()
            G2elms=G2.get_cosets()
            if len(L)>1:
                rel = Graph.rel_cayley(G1,L)
                self.Wedge(rel, self.basepoint, rel.basepoint)
                for a in A:
                    if A[a] not in L:
                        self.glue([self.read_word(self.basepoint, G1elms[a]),
                                   self.read_word(self.basepoint, G2elms[A[a]])]
                        )
        
    def get_red_precover(G1,G2,A,H):
        
        Output=Graph.Step1(G1,G2,A,H)

        Output.Step2(G1, G2, A)

        Output.Step3(G1, G2, A)

        Output.Step4(G1, G2, A)

        Output.Step5(G1, G2, A)

        Output.Step6(G1, G2, A)

        return(Output)
        
        
    def visualizee(self):
        G= nx.MultiDiGraph()
        G.add_nodes_from(range(self.n_verts))
        pos_edges=[]
        for label in self.mat:
            if not label.endswith("^-1"):
                for origin,row in enumerate(self.mat[label]):
                    pos_edges.extend([(origin,end,{"label":label}) for end in np.nonzero(row)[0]])
        G.add_edges_from(pos_edges)
        # print(pos_edges)
        pos = nx.kamada_kawai_layout(G)
        # print(pos)
        fig, ax = plt.subplots(figsize=(15, 10))

        nx.draw_networkx_nodes(G, pos, ax=ax,node_size=25,node_color='white')
        nx.draw_networkx_labels(G, pos, ax=ax,font_size=20)
    
        
        straight_edges = [edge for edge in G.edges if (edge[1],edge[0]) not in G.edges()]
        loop_edges = [edge for edge in G.edges if edge[0]==edge[1]]
        curved_edges = [edge for edge in G.edges if (edge[1],edge[0]) in G.edges() and edge[0]!=edge[1]]
        Multiple_edges={}
            
        R=""
        curves=[edge for edge in G.edges(data=True) if (edge[1],edge[0]) in G.edges() and edge[0]!=edge[1]]
        loops=[edge for edge in G.edges(data=True) if edge[0]==edge[1]]
        for edge in curves:
            if (edge[0],edge[1]) in Multiple_edges:
                Multiple_edges[(edge[0],edge[1])]+=edge[2]["label"]
            else:
                Multiple_edges[(edge[0],edge[1])]=[edge[2]["label"]]
        for edge in curves + loops:
            R += str(edge) + "\n"
        nx.draw_networkx_edges(G, pos,ax=ax, edgelist=straight_edges, width = 2, arrowsize = 20)
        arc_rad = 0.25
        nx.draw_networkx_edges(G, pos,ax=ax, edgelist=curved_edges, connectionstyle=f'arc3, rad = {arc_rad}', width = 2,arrowsize=20)
        nx.draw(G, pos,ax=ax, edgelist=loop_edges,nodelist=[], node_size=600 , node_color='skyblue', font_size=25, font_weight='bold', arrowstyle='->',connectionstyle='arc3, rad = 4', width = 2, arrowsize = 20)

        edge_labels = nx.get_edge_attributes(G,'label')
        curved_edge_labels = {edge[:-1]: edge_labels[edge] for edge in curved_edges}
        straight_edge_labels = {edge[:-1]: edge_labels[edge] for edge in straight_edges}
        loop_edge_labels= {edge[:-1]: edge_labels[edge] for edge in loop_edges}
        my_nx.my_draw_networkx_edge_labels(G, pos, ax=ax, edge_labels=curved_edge_labels,rad = arc_rad, font_size=25)
        nx.draw_networkx_edge_labels(G, pos, ax=ax, edge_labels=straight_edge_labels, font_size=25)
        # nx.draw_networkx_edge_labels(G, pos, ax=ax, edge_labels=loop_edge_labels,label_pos=0.8)
        for loop in loop_edges:
            x,y=pos[loop[0]]
            plt.text(x-0.02, y+0.33, edge_labels[loop], fontsize=21, zorder=10,color="black")
        nx.draw_networkx_edges(G, pos,ax=ax, edgelist=[(0,1,0)], width = 2, arrowsize = 20)#REMOVE
        nx.draw_networkx_edge_labels(G, pos, ax=ax, edge_labels={(0,1):"y"}, font_size=25)#REMOVE
        plt.title("Directed Graph Visualization")
        plt.show()    
        # print(self.datatotex())
            
    def visualize(self):
            G= nx.MultiGraph()
            G.add_nodes_from(range(self.n_verts))
            pos_edges=[]
            for label in self.mat:
                if not label.endswith("^-1"):
                    for origin,row in enumerate(self.mat[label]):
                        pos_edges.extend([(origin,end,label) for end in np.nonzero(row)[0]])
            for edge in pos_edges:
                G.add_edge(edge[0],edge[1],label=edge[2])
            pos = nx.kamada_kawai_layout(G)
            nx.draw(G, pos, with_labels=True, node_size=100, node_color='skyblue', font_size=10, font_weight='bold', arrowstyle='->', arrowsize=10)
            edge_labels = {(source, target): data['label'] for source, target, data in G.edges(data=True)}
            nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels,label_pos=0.5, font_size=8)
    
            plt.title("Directed Graph Visualization")
            plt.show()            
    # def visualize(self):
    #     G= nx.DiGraph()
    #     G.add_nodes_from(range(self.n_verts))
    #     pos_edges=[]
    #     for label in self.mat:
    #         if not label.endswith("^-1"):
    #             for origin,row in enumerate(self.mat[label]):
    #                 pos_edges.extend([(origin,end,label) for end in np.nonzero(row)[0]])
    #     for edge in pos_edges:
    #         G.add_edge(edge[0],edge[1],label=edge[2])
    #     pos = nx.kamada_kawai_layout(G)
    #     nx.draw(G, pos, with_labels=True, node_size=1000, node_color='skyblue', font_size=10, font_weight='bold', arrowstyle='->', arrowsize=10,connectionstyle='arc3, rad = 0.1')
    #     edge_labels = {(source, target): data['label'] for source, target, data in G.edges(data=True)}
    #     nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels,label_pos=0.8, font_size=8)

    #     plt.title("Directed Graph Visualization")
    #     plt.show()
    def datatotex(self):
        f = ""
        G= nx.MultiDiGraph()
        G.add_nodes_from(range(self.n_verts))
        pos_edges=[]
        for label in self.mat:
            if not label.endswith("^-1"):
                for origin,row in enumerate(self.mat[label]):
                    pos_edges.extend([(origin,end,{"label":label}) for end in np.nonzero(row)[0]])
        G.add_edges_from(pos_edges)
        pos = nx.kamada_kawai_layout(G)
        for key in pos:
            f+="\draw (5*{},5*{}) node[circle,scale = 1.3,thick,blue!60,inner sep = 0.4pt,outer sep = 1.7pt]({}){{{}}};\n".format(pos[key][0],pos[key][1],key,key)
        edges={}
        for edge in pos_edges:
            if edge[:-1] not in edges and edge[:-1][::-1] not in edges:
                edges[edge[:-1]]=[edge[-1]["label"]]
            elif edge[:-1] in edges:
                edges[edge[:-1]]+=[edge[-1]["label"]]
            elif edge[:-1][::-1] in edges:
                edges[edge[:-1][::-1]]+=[edge[-1]["label"]+"-1"]
        for pair in edges:
            n=len(edges[pair])
            if pair[0]==pair[1]:
                x,y = pos[pair[0]][0],pos[pair[0]][1]
                zeta = 0
                if x == 0:
                    zeta = 90*y/abs(y)
                elif x > 0:
                    zeta = np.arctan(y/x)*360/(2*np.pi)
                if x < 0:
                    zeta = (np.arctan(y/x) + np.pi)*360/(2*np.pi)
                if n==1:                        
                    f+="\path [thick,draw=black,-{{Stealth}}]\n ({}) to[loop,min distance=15mm,in={}+30,out={}-30,looseness=5] node[scale=1.2,thick,fill=white,circle, anchor=center, pos=0.5,inner sep=1pt,minimum size=4pt]{{${}$}} ({});\n".format(pair[0],zeta,zeta,edges[pair][0],pair[1])
                else:
                    for i, label in enumerate(edges[pair]):
                        f+= "\path [thick,draw=black,-{{Stealth}}]\n ({}) to[loop,min distance=20mm,in={},out={},looseness=5] node[scale=1.2,thick,fill=white,circle, anchor=center, pos=0.5,inner sep=1pt,minimum size=4pt]{{${}$}} ({});\n".format(pair[0],zeta+90-(i+1)*180/(n+2),zeta+90-(i+2)*180/(n+2),label,pair[1])

            else:
                if n==1:
                    f+="\path [thick,draw=black,-{{Stealth}}]\n ({}) --node[fill=white, anchor=center, pos=0.5,inner sep=0.5pt,minimum size=4pt]{{${}$}} ({});\n".format(pair[0],edges[pair][0],pair[1])
                else:
                    for i,label in enumerate(edges[pair]):
                        if label.endswith("-1"):
                            f+= "\path [thick,draw=black,-{{Stealth}}]\n ({}) to[bend right = {}] node[scale=1.2,thick,fill=white, anchor=center, pos=0.5,inner sep=0.7pt,minimum size=4pt]{{${}$}} ({});\n".format(pair[1],90-(i+1)*180/(n+1),label[:-2],pair[0])
                        else:
                            f+= "\path [thick,draw=black,-{{Stealth}}]\n ({}) to[bend left = {}] node[scale=1.2,thick,fill=white, anchor=center, pos=0.5,inner sep=0.7pt,minimum size=4pt]{{${}$}} ({});\n".format(pair[0],90-(i+1)*180/(n+1),label,pair[1])
        print(f)
def Wedge(G1,G2,u,v):
    """
    Does the disjoint union of the two graphs G1, G2 and then 
    identifies the vertex v of G1 with the first vertex of G2.
    So basically, the wedge along the vertices v and 0 of G1 and 
    G2 respectively. 
    Param:
        G1,G2 Graph: Graph objects
        v int: Index of vertex of G1 
    Output:
        The Graph object corresponding to the mentioned wedge
    """
    Output=Graph(list(G1.labels|G2.labels))
    n1, n2= G1.n_verts, G2.n_verts
    N = n1 + n2
    Output.n_verts = N
    for label in Output.mat:
        zero_matrix = np.full((N,N),0)
        if label in G1.mat:
            zero_matrix[:n1,:n1] = G1.mat[label]
        if label in G2.mat:
            zero_matrix[n1:,n1:] = G2.mat[label]
        Output.mat[label]=zero_matrix

    Output.glue([u,n1+v])
    # Output.fold()
    return(Output)




D12=Graph(["a","b"],np.array([[[] for _ in range(12)]for _ in range(12)]))
for i in range(5):
    D12.add_edge(i,i+1,label="a")
D12.add_edge(5, 0, "a")
for i in range(5):
    D12.add_edge(i+6+1,i+6,label="a")
D12.add_edge(6, 11, "a")

for i in range(6):
    D12.add_edge(i,i+6,label="b")
# D12.visualize()

C4=Graph(["c"],np.array([[[] for _ in range(4)]for _ in range(4)]))
for i in range(4):
    C4.add_edge(i,(i+1)%4,label="c")

A = {6:2}
H=[["b"],["a","b","a^-1"],["a","c","a","c","a","c"]]
Output = Graph.get_red_precover(D12,C4,A,H)
Output.visualizee()

#Possible completion to saturated...

# Output.add_edge(11, 6, "c")
# Output.add_edge(6, 13, "c")
# Output.add_edge(13, 8, "c")
# Output.add_edge(8, 11, "c")
# Output.add_edge(12, 7, "c")
# Output.add_edge(7, 12, "c")
# Output.Step3(D12,C4,A)
# Output.Step4(D12,C4,A)
# Output.Step5(D12,C4,A)
# Output.Step6(D12,C4,A)
# Output.visualizee()