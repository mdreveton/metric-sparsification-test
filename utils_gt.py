import graph_tool.all as gt
import numpy as np

def reweighting( G, weighting_type ='Adamic-Adard', weighted = False):
    """
    Reweight the edges of the graph G according to the specified weighting type.
    Options are 'Adamard' or 'Jaccard'.
    """

    if weighting_type.lower() == 'adamic-adard':
        # Adamard reweighting
        weights = []
        for e in G.edges():
            u, v = e
            weights.append( AdamicAdard(G, u, v, weighted = weighted) )
    elif weighting_type.lower() == 'jaccard':
        # Jaccard reweighting
        weights = [ ]
        for e in G.edges():
            u, v = e
            weights.append( Jaccard(G, u, v, weighted = weighted) )
    else:
        raise ValueError("Unknown weighting type: {}".format(weighting_type))
    #G.edge_properties['weight'].set_values( weights )

    weights = G.new_edge_property("double", vals = weights)
    G.edge_properties['weight'] = weights
    return G


def AdamicAdard( G, u, v, weighted = True):
    """
    Compute the Adard similarity between two vertices u and v in graph G
    """
    
    neighbors_u = set(G.get_all_neighbors(u))
    neighbors_v = set(G.get_all_neighbors(v))

    intersection = len(neighbors_u.intersection(neighbors_v))
    if intersection == 0:
        return 0.0
    
    if weighted:
        # If weighted, we use the inverse of the logarithm of the degree
        sum_inv_log_degrees = sum(1 / np.log( G.get_out_degrees([common_neighbor], eweight=G.edge_properties['weight']) for common_neighbor in neighbors_u.intersection(neighbors_v) if G.get_out_degrees([common_neighbor], eweight=G.edge_properties['weight'] ) > 0) )
    else:
        sum_inv_log_degrees = sum(1 / np.log( G.get_out_degrees([common_neighbor], eweight=None) ) for common_neighbor in neighbors_u.intersection(neighbors_v) if G.get_out_degrees([common_neighbor], eweight=None) > 0)
    
    return sum_inv_log_degrees[0] if sum_inv_log_degrees[0] > 0 else 0.0


def Jaccard( G, u, v, weighted = True ):
    """
    Compute the weighted Jaccard similarity between two vertices u and v in graph G
    """
    
    neighbors_u = set(G.get_all_neighbors(u))
    neighbors_v = set(G.get_all_neighbors(v))

    if weighted:
        intersection = sum( np.min( [ G.edge_properties['weight'][u, common_neighbor], G.edge_properties['weight'][v, common_neighbor] ] ) for common_neighbor in neighbors_u.intersection(neighbors_v))
        union = sum( np.max( [ G.edge_properties['weight'][u, common_neighbor], G.edge_properties['weight'][v, common_neighbor] ] ) for common_neighbor in neighbors_u.intersection(neighbors_v))
    else:
        # If not weighted, we just count the edges
        intersection = len(neighbors_u.intersection(neighbors_v))
        union = len(neighbors_u.union(neighbors_v))
    
    if union == 0:
        return 0.0
    return intersection / union