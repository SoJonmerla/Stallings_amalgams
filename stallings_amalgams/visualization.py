import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from .graph import Graph

def my_draw_networkx_edge_labels(
    G,
    pos,
    edge_labels=None,
    label_pos=0.5,
    font_size=10,
    font_color="k",
    font_family="sans-serif",
    font_weight="normal",
    alpha=None,
    bbox=None,
    horizontalalignment="center",
    verticalalignment="center",
    ax=None,
    rotate=True,
    clip_on=True,
    rad=0
):
    """Draw edge labels.

    Parameters
    ----------
    G : graph
        A networkx graph

    pos : dictionary
        A dictionary with nodes as keys and positions as values.
        Positions should be sequences of length 2.

    edge_labels : dictionary (default={})
        Edge labels in a dictionary of labels keyed by edge two-tuple.
        Only labels for the keys in the dictionary are drawn.

    label_pos : float (default=0.5)
        Position of edge label along edge (0=head, 0.5=center, 1=tail)

    font_size : int (default=10)
        Font size for text labels

    font_color : string (default='k' black)
        Font color string

    font_weight : string (default='normal')
        Font weight

    font_family : string (default='sans-serif')
        Font family

    alpha : float or None (default=None)
        The text transparency

    bbox : Matplotlib bbox, optional
        Specify text box properties (e.g. shape, color etc.) for edge labels.
        Default is {boxstyle='round', ec=(1.0, 1.0, 1.0), fc=(1.0, 1.0, 1.0)}.

    horizontalalignment : string (default='center')
        Horizontal alignment {'center', 'right', 'left'}

    verticalalignment : string (default='center')
        Vertical alignment {'center', 'top', 'bottom', 'baseline', 'center_baseline'}

    ax : Matplotlib Axes object, optional
        Draw the graph in the specified Matplotlib axes.

    rotate : bool (deafult=True)
        Rotate edge labels to lie parallel to edges

    clip_on : bool (default=True)
        Turn on clipping of edge labels at axis boundaries

    Returns
    -------
    dict
        `dict` of labels keyed by edge

    Examples
    --------
    >>> G = nx.dodecahedral_graph()
    >>> edge_labels = nx.draw_networkx_edge_labels(G, pos=nx.spring_layout(G))

    Also see the NetworkX drawing examples at
    https://networkx.org/documentation/latest/auto_examples/index.html

    See Also
    --------
    draw
    draw_networkx
    draw_networkx_nodes
    draw_networkx_edges
    draw_networkx_labels
    """
    import matplotlib.pyplot as plt
    import numpy as np

    if ax is None:
        ax = plt.gca()
    if edge_labels is None:
        labels = {(u, v): d for u, v, d in G.edges(data=True)}
    else:
        labels = edge_labels
    text_items = {}
    for (n1, n2), label in labels.items():
        (x1, y1) = pos[n1]
        (x2, y2) = pos[n2]
        (x, y) = (
            x1 * label_pos + x2 * (1.0 - label_pos),
            y1 * label_pos + y2 * (1.0 - label_pos),
        )
        pos_1 = ax.transData.transform(np.array(pos[n1]))
        pos_2 = ax.transData.transform(np.array(pos[n2]))
        linear_mid = 0.5*pos_1 + 0.5*pos_2
        d_pos = pos_2 - pos_1
        rotation_matrix = np.array([(0,1), (-1,0)])
        ctrl_1 = linear_mid + rad*rotation_matrix@d_pos
        ctrl_mid_1 = 0.5*pos_1 + 0.5*ctrl_1
        ctrl_mid_2 = 0.5*pos_2 + 0.5*ctrl_1
        bezier_mid = 0.5*ctrl_mid_1 + 0.5*ctrl_mid_2
        (x, y) = ax.transData.inverted().transform(bezier_mid)

        if rotate:
            # in degrees
            angle = np.arctan2(y2 - y1, x2 - x1) / (2.0 * np.pi) * 360
            # make label orientation "right-side-up"
            if angle > 90:
                angle -= 180
            if angle < -90:
                angle += 180
            # transform data coordinate angle to screen coordinate angle
            xy = np.array((x, y))
            trans_angle = ax.transData.transform_angles(
                np.array((angle,)), xy.reshape((1, 2))
            )[0]
        else:
            trans_angle = 0.0
        # use default box of white with white border
        if bbox is None:
            bbox = dict(boxstyle="round", ec=(1.0, 1.0, 1.0), fc=(1.0, 1.0, 1.0))
        if not isinstance(label, str):
            label = str(label)  # this makes "1" and 1 labeled the same

        t = ax.text(
            x,
            y,
            label,
            size=font_size,
            color=font_color,
            family=font_family,
            weight=font_weight,
            alpha=alpha,
            horizontalalignment=horizontalalignment,
            verticalalignment=verticalalignment,
            rotation=trans_angle,
            transform=ax.transData,
            bbox=bbox,
            zorder=1,
            clip_on=clip_on,
        )
        text_items[(n1, n2)] = t

    ax.tick_params(
        axis="both",
        which="both",
        bottom=False,
        left=False,
        labelbottom=False,
        labelleft=False,
    )

    return text_items

def visualize(H: Graph,
    loop_rad: float = 0.5,
    loop_label_offset: float = 0.13,
    arc_rad: float = 0.25
    ) -> None:
    """
    Visualize a finite directed labelled graph.
     
    Only edges with positive labels are drawn. Inverse-labelled edges
    are omitted because they represent the reverse orientations of the
    corresponding positive-labelled edges.
     
    Edges between distinct vertices are drawn as straight edges unless
    an edge also exists in the opposite direction. In that case, the
    edges are curved to make both orientations visible. Loops are drawn
    separately, and their labels are positioned manually.
     
    The distinguished basepoint is displayed as a larger vertex with a
    black border.
     
    Parameters
    ----------
    H : Graph
    The labelled graph to visualize.
     
    loop_rad : float, default=4.0
    Curvature parameter used when drawing loops.
     
    loop_label_offset : float, default=0.30
    Vertical distance between a vertex and the first loop label at
    that vertex.
     
    arc_rad : float, default=0.25
    Curvature parameter used for edges whose reverse orientation
    is also present.
     
    Raises
    ------
    ValueError
    If the graph has no vertices, if the basepoint is invalid, or
    if one of the layout parameters is negative.
     
    Notes
    -----
    This function modifies neither the graph nor its adjacency matrices.
    It displays the resulting figure using Matplotlib and returns no
    value.
    """
    if H.n_verts < 1:
        raise ValueError("The graph must contain at least one vertex.")
        
    if not 0 <= H.basepoint < H.n_verts:
        raise ValueError(
        f"The basepoint must be between 0 and {H.n_verts - 1}."
        )
    if loop_rad < 0:
        raise ValueError("loop_rad must be non-negative.")
        
    if loop_label_offset < 0:
        raise ValueError("loop_label_offset must be non-negative.")
        
    if arc_rad < 0:
        raise ValueError("arc_rad must be non-negative.")
        
    G= nx.MultiDiGraph()
    G.add_nodes_from(range(H.n_verts))
    # Add one directed edge for each positive label.
    pos_edges = []
    for label in H.mat:
        if label.endswith("^-1"):
            continue
       
        for origin, row in enumerate(H.mat[label]):
            for end in np.nonzero(row)[0]:
                pos_edges.append(
                    (origin, end, {"label": label})
                )
       
    G.add_edges_from(pos_edges)
    
    # Layout
    pos = nx.kamada_kawai_layout(G)

    fig, ax = plt.subplots(figsize=(15, 10))

    # ------------------------------------------------------------------
    # Nodes
    # ------------------------------------------------------------------

    nx.draw_networkx_nodes(
        G,
        pos,
        ax=ax,
        node_size=25,
        node_color="white"
    )

    # Highlight the basepoint
    nx.draw_networkx_nodes(
        G,
        pos,
        ax=ax,
        nodelist=[H.basepoint],
        node_size=100,
        node_color="white",
        edgecolors="black",
        linewidths=2
    )

    # Vertex labels
    nx.draw_networkx_labels(
        G,
        pos,
        ax=ax,
        font_size=20
    )
    # ------------------------------------------------------------------
    # Straight, curved and loop edges
    # ------------------------------------------------------------------
    all_edges = list(G.edges(keys=True))
    
    loop_edges = [
        edge for edge in all_edges
        if edge[0] == edge[1]
    ]
    
    non_loop_edges = [
        edge for edge in all_edges
        if edge[0] != edge[1]
    ]
    
    straight_edges = [
        edge for edge in non_loop_edges
        if (edge[1], edge[0]) not in G.edges()
    ]
    
    curved_edges = [
        edge for edge in non_loop_edges
        if (edge[1], edge[0]) in G.edges()
    ]
    
    
    nx.draw_networkx_edges(
        G,
        pos,
        ax=ax,
        edgelist=straight_edges,
        width=2,
        arrowsize=20
    )
    nx.draw_networkx_edges(
        G,
        pos,
        ax=ax,
        edgelist=curved_edges,
        connectionstyle=f"arc3, rad={arc_rad}",
        width=2,
        arrowsize=20
    )
    
    nx.draw_networkx_edges(
        G,
        pos,
        ax=ax,
        edgelist=loop_edges,
        connectionstyle=f"arc3, rad={loop_rad}",
        width=2,
        arrowsize=20,
        node_size=600
    )

    ## Edge labels ##
    edge_labels = nx.get_edge_attributes(G,'label')
    curved_edge_labels = {edge[:-1]: edge_labels[edge] for edge in curved_edges}
    straight_edge_labels = {edge[:-1]: edge_labels[edge] for edge in straight_edges}
    my_draw_networkx_edge_labels(G, pos, ax=ax, edge_labels=curved_edge_labels,rad = arc_rad, font_size=20)
    nx.draw_networkx_edge_labels(G, pos, ax=ax, edge_labels=straight_edge_labels, font_size=20)
    # ------------------------------------------------------------------
    # Labels on loops
    #
    # NetworkX's automatic loop-label placement is not very reliable,
    # so the labels are positioned manually.
    # ------------------------------------------------------------------
    
    for i, (u, v, k) in enumerate(loop_edges):
        x, y = pos[u]
    
        ax.text(
            x,
            y + loop_label_offset + 0.08 * i,
            edge_labels[(u, v, k)],
            fontsize=21,
            ha="center",
            va="center",
            zorder=10,
            color="black"
        )
    
    plt.title("Directed Graph Visualization")
    plt.show()    
    # print(H.datatotex())
        

def datatotex(H):
    f = ""
    G= nx.MultiDiGraph()
    G.add_nodes_from(range(H.n_verts))
    pos_edges=[]
    for label in H.mat:
        if not label.endswith("^-1"):
            for origin,row in enumerate(H.mat[label]):
                pos_edges.extend([(origin,end,{"label":label}) for end in np.nonzero(row)[0]])
    G.add_edges_from(pos_edges)
    pos = nx.kamada_kawai_layout(G)
    for key in pos:
        f+=f"\\draw (5*{pos[key][0]},5*{pos[key][1]}) node[circle,scale = 1.3,thick,blue!60,inner sep = 0.4pt,outer sep = 1.7pt]({key}){{{key}}};\n"
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
                f+=f"\\path [thick,draw=black,-{{Stealth}}]\n ({pair[0]}) to[loop,min distance=15mm,in={zeta}+30,out={zeta}-30,looseness=5] node[scale=1.2,thick,fill=white,circle, anchor=center, pos=0.5,inner sep=1pt,minimum size=4pt]{{${edges[pair][0]}$}} ({pair[1]});\n"
            else:
                for i, label in enumerate(edges[pair]):
                    f+= f"\\path [thick,draw=black,-{{Stealth}}]\n ({pair[0]}) to[loop,min distance=20mm,in={zeta+90-(i+1)*180/(n+2)},out={zeta+90-(i+2)*180/(n+2)},looseness=5] node[scale=1.2,thick,fill=white,circle, anchor=center, pos=0.5,inner sep=1pt,minimum size=4pt]{{${label}$}} ({pair[1]});\n"

        else:
            if n==1:
                f+=f"\\path [thick,draw=black,-{{Stealth}}]\n ({pair[0]}) --node[fill=white, anchor=center, pos=0.5,inner sep=0.5pt,minimum size=4pt]{{${edges[pair][0]}$}} ({pair[1]});\n"
            else:
                for i,label in enumerate(edges[pair]):
                    if label.endswith("-1"):
                        f+= f"\\path [thick,draw=black,-{{Stealth}}]\n ({pair[1]}) to[bend right = {90-(i+1)*180/(n+1)}] node[scale=1.2,thick,fill=white, anchor=center, pos=0.5,inner sep=0.7pt,minimum size=4pt]{{${label[:-2]}$}} ({pair[0]});\n"
                    else:
                        f+= f"\\path [thick,draw=black,-{{Stealth}}]\n ({pair[0]}) to[bend left = {90-(i+1)*180/(n+1)}] node[scale=1.2,thick,fill=white, anchor=center, pos=0.5,inner sep=0.7pt,minimum size=4pt]{{${label}$}} ({pair[1]});\n"
    print(f)