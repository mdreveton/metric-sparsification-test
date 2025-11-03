#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  4 17:26:34 2025

@author: dreveton
"""

import igraph as ig
import numpy as np
import scipy as sp
import random as random 
from tqdm import tqdm 



"""

n = 200000
n_clusters = 2
alpha = [ 1 / n_clusters for i in range( n_clusters ) ]

block_sizes = [ int( alpha[i] * n ) for i in range( n_clusters ) ]

pin = 5
pout = 1 
B = (pin-pout) * np.eye( n_clusters ) + pout * np.ones( n_clusters )
asymptotic_rate = np.log(n) / n 
B = np.array( [ [1,2],[2,8] ])


G = ig.Graph.SBM( np.sum(block_sizes), B * asymptotic_rate, block_sizes, directed=False, loops=False)


lambdain = 1
lambdaout = 1
weights_parameters = (lambdain - lambdaout) * np.eye( n_clusters ) + lambdaout * np.ones( n_clusters )

G = add_weights( G, weights_parameters, block_sizes, weight_distribution = 'exponential' )


source_vertices = [ i for i in range( block_sizes[0]) ]
target_vertices = [ i for i in range( block_sizes[0]+1, block_sizes[0]+block_sizes[1] ) ]

costs_mean = np.zeros((2,2))
costs_std = np.zeros((2,2))

costs_mean[0,0], costs_std[0,0] = estimateCostShortestPath( G, 100, source_vertices, source_vertices )
costs_mean[0,1], costs_std[0,1] = estimateCostShortestPath( G, 100, source_vertices, target_vertices )
costs_mean[1,1], costs_std[1,1] = estimateCostShortestPath( G, 100, target_vertices, target_vertices )
costs_mean[1,0], costs_std[1,0] = estimateCostShortestPath( G, 100, target_vertices, source_vertices )

T = np.multiply( weights_parameters, B) @ np.diag( alpha )
vals, vecs = np.linalg.eigh(T)
rho = np.max( vals )


gumbels = sp.stats.gumbel_r.rvs( size = 1000 ) + sp.stats.gumbel_r.rvs( size = 1000 ) - sp.stats.gumbel_r.rvs( size = 1000 )

rescaled_costs = distributionalLimitRescaledCosts( G, 100, source_vertices, target_vertices, rho, verbose = True ) 
sp.stats.ks_2samp(rescaled_costs, gumbels )

"""





def estimateCostShortestPath( G, n_samples, source_vertices, target_vertices, verbose = True ):
    
    if verbose:
        iteration = tqdm( range( n_samples ) )
    else:
        iteration = range( n_samples )
        
    cost = [ ]
    for i in iteration:
        #(u,v) = edges_to_use[ i ]
        u = random.choice( source_vertices )
        v = random.choice( target_vertices )
        distances = G.distances( source = u, target = v, weights = 'weight' )
        cost.append( distances[0][0] )
        
    cost = np.asarray( cost, dtype='float64' )
    return np.mean( cost ), np.std( cost ) / np.sqrt( n_samples )
    

def distributionalLimitRescaledCosts( G, n_samples, source_vertices, target_vertices, rho, verbose = True ):
    
    if verbose:
        iteration = tqdm( range( n_samples ) )
    else:
        iteration = range( n_samples )

    rescaled_costs = []
    n = G.vcount( )
    
    for i in iteration:
        #(u,v) = edges_to_use[ i ]
        u = random.choice( source_vertices )
        v = random.choice( target_vertices )
        distances = G.distances( source = u, target = v, weights = 'weight' )
        rescaled_costs.append( np.log(n) * ( rho * distances[0][0] - 1 ) )

    return rescaled_costs


def add_weights( G, weights_parameters, block_sizes, weight_distribution = 'exponential' ):
    
    
    #A = G.get_adjacency_sparse().asfptype().tolil() #.asfptype() is important to obtain the matrix with float (and not int) values
    
    A = sp.sparse.triu( G.get_adjacency_sparse() ).tolil().asfptype()
    
    # triu retains only the upper diagonal and .asfptype() is important to obtain the matrix with float (and not int) values
    
    n_clusters = len( block_sizes )
    
    for a in range( n_clusters ):
        row_index = int( np.sum( [ block_sizes[a] for dummy in range(a) ] ) )
            
        for b in range( n_clusters ):
            colum_index = int( np.sum( [ block_sizes[b] for dummy in range(b) ] ) )
            temp = A[ row_index:row_index+block_sizes[a], colum_index:colum_index+block_sizes[b] ].tolil( )
            #print( temp.shape )
            edges_between_community_a_b = list( temp.nonzero() )
            edges_between_community_a_b = edges_between_community_a_b
            edges_between_community_a_b[ 0 ] = edges_between_community_a_b[ 0 ] + row_index 
            edges_between_community_a_b[ 1 ] += colum_index
            A[ tuple( edges_between_community_a_b ) ] = generateWeights( weight_distribution, weights_parameters[ a, b ],  len( edges_between_community_a_b[0] ) )
            #sp.stats.expon.rvs( scale = 1 / weights_parameters[ a, b ], size = len( edges_between_community_a_b[0] ) )

    return ig.Graph.Weighted_Adjacency( A + A.T, mode='undirected')

def generateWeights( distribution, parameter, size ):

    if distribution == 'exponential':
        return sp.stats.expon.rvs( scale = 1 / parameter, size = size )
    elif distribution == 'uniform':
        return sp.stats.uniform.rvs( loc = 0, scale = 1 / parameter, size = size )
    else:
        raise ValueError( "Unknown distribution: ", distribution )