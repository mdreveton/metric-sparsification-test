#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  5 22:14:00 2025

@author: dreveton
"""

import graph_tool.all as gt
import utils as utils
import scipy as sp 
import numpy as np
import igraph as ig 
import time

import random as random 

from tqdm import tqdm 

import metric_sparsification as metric_sparsification



transportation_graphs = ['chicago_road'  ]
informational_graphs = [ 'word_adjacency/darwin', 'word_adjacency/french', 'word_adjacency/spanish', 'word_adjacency/japanese', 'dblp_cite', 'foldoc', 'google', 'word_assoc', 'cora', 'movielens_100k', 'webkb/webkb_cornell_cocite', 'webkb/webkb_washington_cocite', 'webkb/webkb_texas_cocite', 'webkb/webkb_wisconsin_cocite' ]
collaboration_graphs = ['cofe', 'physics_collab/arXiv' ]
social_graphs = ['sp_infectious', 'wiki_rfa', 'anybeat', 'escorts', 'marvel_universe', 'digg_reply' ]
biological_graphs = ['tree-of-life', 'genetic_multiplex', 'fly_hemibrain']
technological_graphs = ['internet_as', 'caida_as'] 



#Note to myself: largest graph (in terms of number of vertices ) is is digg_reply



"""


graph_statistic = dict( )
results_clever = dict( )
results_naive = dict( )

for graph_name in tqdm( informational_graphs ):
    print( 'Current graph : ', graph_name )

    G_original = retrieve_from_graphtool( graph_name )
    weighted = G_original.is_weighted()
    
    
    G_similarity = utils.reweighting( G_original, similarity = 'Jaccard', weighted = weighted )
    G_distance = utils.similarity_to_distance_graph( G_similarity )

    start_time = time.time()
    
    G_mb = metric_sparsification.metric_backbone( G_distance )
    computation_time = time.time() - start_time
    print("Computing the metric backbone took --- %s seconds ---" % (computation_time) )
    
    graph_statistic[ graph_name ] = getGraphsStatistics( G_original, computation_time = computation_time )

    results_clever[ graph_name ] = getResults(G_distance, G_mb, method = 'clever' )
    results_naive[ graph_name ] = getResults(G_distance, G_mb, method = 'naive' )
    

G_mb_uspt = metric_sparsification.union_shortest_path_trees( G_distance, source_vertices = 'all', target_vertices = 'all', verbose = True )

estimate_mb_density( G_distance )

relative_error = 

"""



def metric_backbone_density_estimation( G, method = 'naive', n_samples = 100 ):
    
    if method == 'naive':
        edges_id = random.choices( G.es.indices, k = n_samples )
        
        success = 0
        for e in edges_id:
            e = G.es[e]
            source = e.source
            target = e.target
            distance = G.distances( source, target, weights = 'weight' )[ 0 ]
            if distance == G[ source, target ]:
                success += 1
        #error = sp.stats.binomtest( k = success, n = n_samples, p = success / n_samples ).proportion_ci()
        #return success / n_samples, error
        return success / n_samples * G.density( )
    
    else:
        kernel = sp.stats.gaussian_kde( G.es['weight'] ) 
        n = G.vcount( )
        p0 = G.density( ) * n / np.log( n )

        return G.density() * kernel.integrate_box_1d( 0, 1/(p0 * kernel.pdf(0) ) )

            
        
    
    

def getGraphsStatistics( G, computation_time = None ):
    graph_statistic = dict( )
    
    graph_statistic[ 'n' ] = G.vcount()
    graph_statistic[ 'E' ] = G.ecount()
    graph_statistic[ 'mean degree' ] = np.mean( G.degree( G.vs.indices ) )
    graph_statistic[ 'std degree' ] = np.std( G.degree( G.vs.indices ) )
    graph_statistic[ 'clustering coefficient' ] = G.transitivity_undirected()
    
    if isinstance( computation_time , (int, float) ):
        graph_statistic[ 'metric backbone computation time' ] = computation_time

    return graph_statistic


def getResults( G, G_mb, method = 'clever' ):
    
    
    p = G.density()
    p_mb = G_mb.density( )
    
    predicted_p_mb = metric_backbone_density_estimation( G, method = method )
    
    relative_error = np.abs( predicted_p_mb - p_mb ) / p_mb * 100

    ratio_deleted_edges = (1 - p_mb / p ) * 100
    predicted_ratio_deleted_edges = ( 1 - predicted_p_mb / p ) * 100

    result = dict( )
    result[ 'original density' ] = G.density( )
    result[ 'metric backbone density' ] = G_mb.density( )
    result[ 'predicted p_mb' ] = predicted_p_mb
    result[ 'relative error density prediction' ] = relative_error
    result[ 'empirical ratio deleted edges' ] = ratio_deleted_edges
    result[ 'predicted ratio deleted edges' ] = predicted_ratio_deleted_edges

    return result 


def retrieve_from_graphtool( graph_name ):
    """
    Retrieve a graph from the graph-tool collection by its name.
    
    graph_name: str
        The name of the graph in the graph-tool collection.
    
    Returns:
        g: graph_tool.Graph
            The graph retrieved from the collection.
    """
    
    g = gt.collection.ns[graph_name]
    
    try:
        try:
            A = gt.adjacency(g, weight=g.edge_properties['weight'])
            weighted = True
            
        except KeyError:
            A = gt.adjacency(g)
            weighted = False
    
        A = (A+A.T)/2
        G_original = ig.Graph().Weighted_Adjacency( A, mode = 'undirected' )
        G_original = utils.get_largest_component( G_original )
        return G_original.simplify( loops=True, combine_edges=sum)
        
    except KeyError:
        raise ValueError(f"Graph '{graph_name}' not found in the collection.")
    
    

def estimate_mb_density( G ):
    """ 
    Estimate the edge density of the metric backbone of a graph distance G.
        using the kernel density estimation of the empirical weight distribution 
        and the original edge density of G

    G: igraph.Graph
        A graph with edges weighted as distances.
    """

    kernel = utils.estimate_weight_pdf( G )
    original_edge_density = G.density()
    n = G.vcount()
    p0 = original_edge_density * n / np.log(n)
    return original_edge_density * kernel.integrate_box_1d(0,1/(p0 * kernel.pdf(0)) )



"""

results = dict( )
graph_statistics = dict( )

for graph_name in transportation_graphs:
    print( graph_name )
    
    g = gt.collection.ns[ graph_name ]
    try:
        A = gt.adjacency(g, weight=g.edge_properties['weight'])
        weighted = True
    except KeyError:
        A = gt.adjacency(g)
        weighted = False
    
    A = (A+A.T)/2
    G_original = ig.Graph().Weighted_Adjacency( A, mode = 'undirected' )
    G_original = utils.get_largest_component( G_original )

    G_similarity = utils.reweighting( G_original, similarity = 'Jaccard', weighted = weighted )
    G_mb = metric_sparsification.metric_backbone( G_similarity )

    W = utils.reweighting( A, similarity = 'Jaccard', weighted = weighted )
    D = utils.similarityToDistanceGraph( W )

    
    G = ig.Graph().Weighted_Adjacency( D, mode = 'undirected' )
    distances = G.distances( weights = 'weight' )

    D_mb = D.copy()
    dummy =0
    for (u,v) in zip( *W.nonzero() ):
        if D[u,v] != distances[u][v]:
            D_mb[u,v] = 0.0
            dummy+=1
    D_mb.eliminate_zeros()
    
    final_edge_density = utils.estimate_edgeDensity( D_mb )
    original_edge_density = utils.estimate_edgeDensity( D )
    
    n = D_mb.shape[ 0 ]
    p0 = original_edge_density * n / np.log(n)

    kernel = utils.estimate_weight_pdf( D )
    
    predicted_edge_density = original_edge_density * kernel.integrate_box_1d(0,1/(p0 * kernel.pdf(0)) )
    
    
    relative_error = np.abs( predicted_edge_density - final_edge_density ) / final_edge_density * 100
    
    ratio_deleted_edges = (1 - final_edge_density / original_edge_density ) * 100
    predicted_ratio_deleted_edges = ( 1 - predicted_edge_density / original_edge_density ) * 100

    results_given_graph = dict( )
    results_given_graph[ 'p' ] = original_edge_density
    results_given_graph[ 'empirical p_mb' ] = final_edge_density
    results_given_graph[ 'predicted p_mb' ] = predicted_edge_density
    results_given_graph[ 'relative error prediction p_mb' ] = relative_error
    results_given_graph[ 'empirical ratio deleted edges' ] = ratio_deleted_edges
    results_given_graph[ 'predicted ratio deleted edges' ] = predicted_ratio_deleted_edges
    
    results[ graph_name ] = results_given_graph
    
    
    
"""



