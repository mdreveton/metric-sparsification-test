#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  9 20:30:50 2025

@author: dreveton
"""

import multiprocessing as mp
#import rustworkx as rx
import igraph as ig
from tqdm import tqdm


"""
import experiments_ErdosRenyi as er
"""



def makeUnionShortestPathTree_mp( G, source_vertices, target_vertices = 'all' , verbose = True, n_threads = 'auto' ):
    
    if n_threads == 'auto':
        n_threads = max( mp.cpu_count() - 1, 1 )
    
    args = [ [i for i in range(thread * len(source_vertices) // n_threads, (thread+1) * len(source_vertices) // n_threads )] for thread in range(n_threads)  ]
    
    pool = mp.Pool()
    for i in range( n_threads ):
        pool.apply_async( makeShortestPathTree, args = ( G, args[i] ), callback = log_result )
    pool.close()
    pool.join()
    
    unique_edges_id = list( )
    for u in range( len( result_list ) ):
        unique_edges_id += result_list[ u ]
        unique_edges_id = list( set(unique_edges_id) )
    unique_edges_id = set( unique_edges_id )

    edges = [ ]
    weights = [ ]
    for elt in unique_edges_id:
        e = G.es[elt]
        edges.append( (e.source, e.target) )
        weights.append( e['weight'] )
        
    SPT = ig.Graph( n = G.vcount() )
    SPT.add_edges( edges, attributes = { 'weight' : weights } )

    return SPT

    
    

def makeShortestPathTree( G, source_vertices, target_vertices = 'all' , verbose = True ):
    
    #G = arg[0]
    #source_vertices = arg[1]
    
    if target_vertices == 'all':
        target_vertices = G.vs.indices

    #SPT = ig.Graph( n = G.vcount() )
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
        
    return unique_edges_id

result_list = [ ]
def log_result(result):
    # This is called whenever foo_pool(i) returns a result.
    # result_list is modified only by the main process, not the pool workers.
    result_list.append(result)


"""
if __name__ == '__main__':
    n = 50000
    p = 2*np.log(n) / n
    G = ig.Graph.Erdos_Renyi( n, p )
    G = er.add_weights( G )
    
    n_threads = 8
    source_vertices = [i for i in range( n_threads * 200 ) ]
    args = [ [i for i in range(thread * len(source_vertices) // n_threads, (thread+1) * len(source_vertices) // n_threads )] for thread in range(n_threads)  ]
    
    pool = mp.Pool()
    for i in range( n_threads ):
        pool.apply_async(makeShortestPathTree, args = (G,args[i]), callback = log_result)
    pool.close()
    pool.join()
"""


"""
    args = [ (G, [i for i in range(thread * len(source_vertices) // n_threads, (thread+1) * len(source_vertices) // n_threads )]) for thread in range(n_threads)  ]
    with mp.Pool(n_threads) as p:
        results = p.map(makeShortestPathTree, args) 

    unique_edges_id = list( )
    for u in range( len(results) ):
        unique_edges_id += results[ u ]
        unique_edges_id = list( set(unique_edges_id) )
"""
