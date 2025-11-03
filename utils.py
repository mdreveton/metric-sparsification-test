#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  5 18:19:38 2025

@author: dreveton
"""



import numpy as np
import scipy as sp
import igraph as ig

""" 
Utility functions for graph processing and analysis.
"""

def get_largest_component( G ):
    """
    Get the largest connected component of the graph G.
    
    G: igraph.Graph
        A graph.
    
    Returns:
        G_largest: igraph.Graph
            The largest connected component of G.
    """
    components = G.components()
    largest_component = max(components, key=len)
    G_largest = G.subgraph(largest_component)
    
    return G_largest 



def get_neighbors( G, u ):
    return G.neighbors(u)

def get_degree( G, u, weighted = True ):
    
    if weighted:
        return G.strength( u, weights = 'weight' )
    else:
        return G.degree( u )





def reweighting( G, similarity = 'Jaccard', weighted = True ):
    """
    Reweight the edges of the graph G according to the specified similarity measure.
    
    G: igraph.Graph
        A graph with edges weighted as distances.
    similarity: str
        The type of similarity to use for reweighting ('Jaccard' or 'Adamic-Adard').
    weighted: bool
        Whether the graph is weighted or not.
    
    Returns:
        W: scipy sparse matrix
            The reweighted adjacency matrix.
    """

    #TODO REWRITE THIS FUNCTION TO USE IGRAH AND NOT SCIPY SPARSE MATRIX
    
    if not isinstance(G, ig.Graph):
        raise ValueError("Input must be an igraph.Graph object.")
    
    G_similarity= G.copy()
    
    for e in G_similarity.es:
        u = e.source
        v = e.target 
        G_similarity[u,v] = compute_edgeSimilarity(G, u, v, weighted=weighted, similarity=similarity)
    #W = G.get_adjacency_sparse( attribute = 'weight' )
    
    
    return G_similarity #ig.Graph().Weighted_Adjacency( _reweighting(W, similarity=similarity, weighted=weighted), mode = 'undirected' )



def compute_edgeSimilarity( G, u, v, weighted = True, similarity = 'Jaccard' ):
    
    if similarity.lower() == 'jaccard':
        """
        Compute the weighted Jaccard similarity between two vertices u and v in graph G
        """
        neighbors_u = set( G.neighbors(u) ).union( {u,v} )
        neighbors_v = set( G.neighbors(u) ).union( {v,v} )
        
        if weighted == True:
            intersection = sum( np.min( [ G[u, common_neighbor], G[v, common_neighbor] ] ) for common_neighbor in neighbors_u.intersection(neighbors_v))
            union = sum( np.max( [ G[u, common_neighbor], G[v, common_neighbor] ] ) for common_neighbor in neighbors_u.union(neighbors_v) )
            
            if G[u,u] == 0 and G[v,v] == 0:
                intersection += 1
        
        else:
            intersection = len(neighbors_u.intersection(neighbors_v))
            union = len(neighbors_u.union(neighbors_v))
            
        return intersection / union 
    
    
    elif similarity.lower() == 'adamic-adard':
        """
        Compute the Adard similarity between two vertices u and v in graph G
        """
        
        neighbors_u = set( G.neighbors(u) )
        neighbors_v = set( G.neighbors(v) )
        if weighted:
             # If weighted, we use the inverse of the logarithm of the degree
             sum_inv_log_degrees = sum( 1 / np.log( get_degree ( G, common_neighbor, weighted=True) ) for common_neighbor in neighbors_u.intersection(neighbors_v) if get_degree( G, common_neighbor, weighted=True) > 0)
        else:
             sum_inv_log_degrees = sum(1 / np.log( get_degree ( G, common_neighbor, weighted=False) ) for common_neighbor in neighbors_u.intersection(neighbors_v) if get_degree( G, common_neighbor, weighted=False) > 0)
     
        return sum_inv_log_degrees if sum_inv_log_degrees > 0 else 0.0




def distance_to_similarity_graph( G_distance ):
    G_similarity = G_distance.copy()
    
    for e in G_similarity.es:
        G_similarity[ e.source, e.target ] = 1 / ( G_distance[e.source, e.target] + 1 )
    
    return G_similarity

def similarity_to_distance_graph( G_similarity ):
    G_distance = G_similarity.copy()
    
    for e in G_distance.es:
        G_distance[ e.source, e.target ] = 1 / G_similarity[e.source, e.target] - 1.0 
    
    return G_distance


def estimate_weight_pdf( G ):
    return sp.stats.gaussian_kde( G.es['weight'] )
    




# =============================================================================
# OLD FUNCTIONS USING SCIPY SPARSE MATRICES
# =============================================================================

def _get_neighbors( W, u ):
    
    return list( W[u].nonzero( )[ 1 ] )

def _get_degree( W, u, weighted = True ):
    
    if weighted:
        return np.sum( W[u] )
    else:
        return np.sum( _get_neighbors( W, u ) )



def _reweighting( W, similarity ='Jaccard', weighted = True ):
    """
    Reweight the edges of the weighted matrix W according to the specified weighting type.
    Options for similarity are 'Adamard' or 'Jaccard'.
    """
    weights = []
    
    for (u,v) in zip( *W.nonzero() ):
        weights.append( _compute_edgeSimilarity( W, u, v, weighted = weighted, similarity = similarity)  )
    
    W[ W.nonzero()] = weights
    return W



def _compute_edgeSimilarity( W, u, v, weighted = True, similarity = 'Jaccard'):
    
    if similarity.lower() == 'jaccard':
        """
        Compute the weighted Jaccard similarity between two vertices u and v in graph G
        """
        neighbors_u = set(_get_neighbors(W, u)).union({u,v})
        neighbors_v = set(_get_neighbors(W, v)).union({u,v})
        W[u,u] = 1.0
        W[v,v] = 1.0
        if weighted == True:
            intersection = sum( np.min( [ W[u, common_neighbor], W[v, common_neighbor] ] ) for common_neighbor in neighbors_u.intersection(neighbors_v))
            union = sum( np.max( [ W[u, common_neighbor], W[v, common_neighbor] ] ) for common_neighbor in neighbors_u.union(neighbors_v) )
        else:
            intersection = len(neighbors_u.intersection(neighbors_v))
            union = len(neighbors_u.union(neighbors_v))
            
        return intersection / union 
    
    
    elif similarity.lower() == 'adamic-adard':
        """
        Compute the Adard similarity between two vertices u and v in graph G
        """
        
        neighbors_u = set(_get_neighbors( W, u ))
        neighbors_v = set(_get_neighbors( W, v ))
        if weighted:
             # If weighted, we use the inverse of the logarithm of the degree
             sum_inv_log_degrees = sum( 1 / np.log( _get_degree ( W, common_neighbor, weighted=True) ) for common_neighbor in neighbors_u.intersection(neighbors_v) if _get_degree( W, common_neighbor, weighted=True) > 0)
        else:
             sum_inv_log_degrees = sum(1 / np.log( _get_degree ( W, common_neighbor, weighted=False) ) for common_neighbor in neighbors_u.intersection(neighbors_v) if _get_degree( W, common_neighbor, weighted=False) > 0)
     
        return sum_inv_log_degrees if sum_inv_log_degrees > 0 else 0.0


        
