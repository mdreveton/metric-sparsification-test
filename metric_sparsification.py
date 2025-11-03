#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 22 11:31:54 2025

Backbone Subgraph - Fast Implementation
==================

Computes the shortest path distance backbone on a weighted graph.
These algorithms work with edges weighted as distances.


@author: dreveton
"""


"""
import time

start_time = time.time()
distances = G.distances( source=current_sources, target = target_left, weights='weight' ) #Compute the geodesic distances between sources vertices
end_time = time.time()

print("--- %s seconds ---" % (end_time - start_time))


"""


import numpy as np
#import rustworkx as rx
import igraph as ig
from tqdm import tqdm


__all__ = [
    "metric_backbone",
]


__kinds__ = ['metric']
__algorithms__ = ['rustworkx', 'dijkstra']


def removeSemiMetric( G, max_order = 2, verbose = False ):
    
    G_star = G.copy()
    
    for order in tqdm( range(1, max_order+1) ):
        G_star = removeSemiMetricOfGivenOrder( G_star, order = order )
        
    return G_star



def removeSemiMetricOfGivenOrder( G, order = 1 ):
    
    G_star = G.copy( )
    G_star.vs["name"] = [ i for i in range( G.vcount() ) ]
    
    edges_to_remove = [ ]
    
    neighborhoods = G.neighborhood( G.vs, order=order )
    
    for edge in tqdm( G.es ):
        u, v = edge.source, edge.target
        common_neighbors = set( neighborhoods[ u ] ).intersection( set( neighborhoods[ v ] ) )
        
        if order == 1:
            for neigh in common_neighbors:
                if G[u,v] > G[u, neigh ] + G[neigh, v]:
                    edges_to_remove.append( edge )
                    
        else:
            neighborhood_graph = G_star.subgraph( common_neighbors )
            new_index_of_u = neighborhood_graph.vs.select( name_eq = u )[0].index
            new_index_of_v = neighborhood_graph.vs.select( name_eq = v )[0].index
            if neighborhood_graph[new_index_of_u, new_index_of_v] > neighborhood_graph.distances(new_index_of_u, new_index_of_v, weights='weight' )[0][0]:
                edges_to_remove.append ( edge )

    G_star.delete_edges( edges_to_remove )
    
    return G_star
    


def metric_backbone( G, method = 'direct', start_by_removing_semiMetricEdges = True, verbose = False ):
    """
    Compute the metric backbone of a distance graph G.

    G: igraph.Graph
        A graph with edges weighted as distances.
    
    Returns:
        G_mb: igraph.Graph
            The metric backbone of the graph G.
    """
    if not isinstance(G, ig.Graph):
        raise TypeError("Input must be an igraph.Graph object.")
    
    if G.vcount() > 30000 and method == 'direct':
        print( 'The direct method will likely lead to a memory crash. I will use the batch one instead' )
        method = 'batch'
    
    if method == 'auto':
        method = 'batch'
        start_by_removing_semiMetricEdges = True
    
    if method == 'direct':
        distances = G.distances( weights = 'weight' ) #Compute all pairs shortest path distances
        G_mb = G.copy()
        edges_to_remove = [ e for e in G.es if G[e.source, e.target] > distances[e.source][e.target] ]
        #for e in G.es:
        #    if G[e.source, e.target] != distances[e.source][e.target]:
        #        edges_to_remove.append(e)
        G_mb.delete_edges( edges_to_remove )
    
    else:
        if start_by_removing_semiMetricEdges:
            G_mb = removeSemiMetric( G, max_order = 2, verbose = verbose )
        else:
            G_mb = G.copy()
        
        if method == 'batch':
            G_mb = DistancesByBatch( G_mb, maximal_space_limitation = 100000000 )
    
    return G_mb



def allPairShortestPathGraph( G, verbose = False ):
    
    apsp = ig.Graph( n = G.vcount() )
    unique_edges_id = list( )
    
    if verbose:
        iteration = tqdm( range( G.vcount() ) )
    else:
        iteration = range( G.vcount() )

    for source in iteration:
        path = G.get_shortest_paths(source, to = [i for i in range(source, G.vcount())], weights='weight', output='epath')
    
        for u in range( len(path) ):
            unique_edges_id += path[ u ]
        unique_edges_id = list( set(unique_edges_id) )
    
    edges = [ ]
    weights = [ ]
    for elt in unique_edges_id:
        e = G.es[elt]
        edges.append( (e.source, e.target) )
        weights.append( e['weight'] )
    apsp.add_edges( edges, attributes = { 'weight' : weights } )
    
    return apsp


def DistancesByBatch( G, maximal_space_limitation = 100000000, verbose = True ):
    
    maximal_space_limitation = int( maximal_space_limitation )
    n = G.vcount()
    batch_size = maximal_space_limitation // n    
    number_of_batchs = int( np.ceil(n / batch_size ) )
    
    target_left = [ i for i in range(n) ]
    G_mb = G.copy()
    
    if verbose:
        iteration = tqdm( range( number_of_batchs ) )
        print( 'We will start the distance computations batch by batch' )
    
    else:
        iteration = range( number_of_batchs )

    for dummy in iteration:
    #while len( target_left ) > 0:
        edges_to_remove = [ ]

        current_sources = [ target_left[i] for i in range( batch_size ) ]
        distances = G_mb.distances( source=current_sources, target = target_left, weights='weight' ) #Compute the geodesic distances between sources vertices
        
        for u in current_sources:
            index_shift = dummy * batch_size
            index_u = u - index_shift
            temp = set( target_left )
            edges_to_remove += [  (min(u,v), max(u,v) ) for v in set(G_mb.neighbors(u)).intersection( temp ) if G[u,v] > distances[ index_u ][ v - index_shift ] ]
        #Note to myself: v - dummy * batch is actually equal to target_left.index( v )
        #and similarly u - dummy * batch is equal to current_sources.index( u )
        
        target_left = [x for x in target_left if x not in current_sources]

        G_mb.delete_edges( edges_to_remove )

    return G_mb
    


def union_shortest_path_trees( G, source_vertices = 'all', target_vertices = 'all', verbose = False ):
    """ 
    Compute the union of shortest path trees from source vertices to target vertices.
    
    G: igraph.Graph
        A graph with edges weighted as distances.
    source_vertices: list
        List of source vertices.
    target_vertices: list
        List of target vertices. Default is all
    
    Returns:
        G_union: igraph.Graph
            The union of shortest path trees.
    """
    
    G_union_spt = ig.Graph( n = G.vcount() )
    
    if source_vertices=='all':
        source_vertices = G.vs.indices
        apsp = True
    else:
        apsp = False
    
    if target_vertices == 'all':
        target_vertices = G.vs.indices


    if verbose:
        iteration = tqdm( source_vertices )
    else:
        iteration = source_vertices
    
    for source in iteration:
        if apsp:
            SPT_source = makeShortestPathTree( G, [source], target_vertices = range(source , G.vcount() ) )
        else:
            SPT_source = makeShortestPathTree( G, [source], target_vertices = target_vertices )

        G_union_spt = G_union_spt.union( SPT_source )
    
    return G_union_spt


def makeShortestPathTree( G, source_vertices, target_vertices = 'all' , verbose = False ):
    
    if target_vertices == 'all':
        target_vertices = G.vs.indices

    SPT = ig.Graph( n = G.vcount() )
    unique_edges_id = list( )
    
    if verbose:
        iteration = tqdm( source_vertices )
    else:
        iteration = source_vertices


    for source in iteration:
        path = G.get_shortest_paths(source, weights='weight', output='epath')

        for u in range( len(path) ):
            unique_edges_id += path[ u ]
        
        unique_edges_id = list( set(unique_edges_id) )
    
    edges = [ ]
    weights = [ ]
    for elt in unique_edges_id:
        e = G.es[elt]
        edges.append( (e.source, e.target) )
        weights.append( e['weight'] )
    SPT.add_edges( edges, attributes = { 'weight' : weights } )
    
    return SPT
    


"""
def metric_backbone_fast_rx( G ):
    g = rx.networkx_converter( G )
    APSP_distances = rx.all_pairs_dijkstra_path_lengths( g , edge_cost_fn = lambda edge: edge['weight'] )

    G_mb = G.copy()
    for edge in G.edges:
        u = edge[ 0 ]
        v = edge[ 1 ]
        if G[ u ][ v ]['weight'] > APSP_distances[ u ][v ]:
            G_mb.remove_edge( u, v )
    return G_mb 
""" 

def thresholdGraph( G, threshold, graph_type = 'distance' ):
    """
    Take a DISTANCE graph and return the threshold subgraph
    """
    G_threshold = G.copy()
    
    if graph_type == 'distance':
        for edge in G.edges:
            u = edge[ 0 ]
            v = edge[ 1 ]
            if G[ u ][ v ][ 'weight' ] > threshold:
                G_threshold.remove_edge( u, v )
                
    elif graph_type == 'similarity':
        for edge in G.edges:
            u = edge[ 0 ]
            v = edge[ 1 ]
            if G[ u ][ v ][ 'weight' ] > threshold:
                G_threshold.remove_edge( u, v )
    else:
        raise TypeError( 'Graph type not implemented' )
        
    return G_threshold




def _check_for_kind(kind):
    """
    Check for available metric functions.
    """
    if kind not in __kinds__:
        raise TypeError("Metric not found for this algorithm. Try 'metric' or 'ultrametric',")