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

import matplotlib.pyplot as plt


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


"""
# =============================================================================
# TESTING DISTRIBUTIONAL LIMIT FOR AN ER GRAPH
# =============================================================================

n = 200000
n_clusters = 2
alpha = [ 1 / n_clusters for i in range( n_clusters ) ]
block_sizes = [ int( alpha[i] * n ) for i in range( n_clusters ) ]

asymptotic_rate = np.log(n)
B = np.array( [ [4,4],[4,4] ] )

lambdain = 1
lambdaout = 1
weights_parameters = (lambdain - lambdaout) * np.eye( n_clusters ) + lambdaout * np.ones( n_clusters )

T = np.multiply( weights_parameters, B) @ np.diag( alpha ) * asymptotic_rate - np.ones( n_clusters )
vals, vecs = np.linalg.eigh( T.T )
rho_theoretic = np.max( vals ) / np.log(n)
rho_range = np.linspace( start = 2, stop = 6, num=41 )

n_trials = 10
ks_distance_er = np.zeros( (n_trials, len(rho_range )) )

c = 0 

for trial in tqdm( range( n_trials ) ):
    G = ig.Graph.SBM( np.sum(block_sizes), B * asymptotic_rate / n, block_sizes, directed=False, loops=False)
    G = add_weights( G, weights_parameters, block_sizes, weight_distribution = 'exponential' )
    ks_distance_er[trial,:] = testingDistributionalLimit( G, rho_range, n_samples = 100, source_vertices = [i for i in range(n)], target_vertices = [i for i in range(n)], c = 0 ) 

ks_distance_er_mean = np.mean(ks_distance_er, axis=0)
ks_distance_er_std = np.std(ks_distance_er, axis=0)

fileName = 'KS_rho_er_n_' + str(n) + '_B_' + str(B) + '_alpha_' + str(alpha) + '.pdf'
plotFigure( rho_range, ks_distance_er_mean, accuracy_err = ks_distance_er_std, methods = None, 
           vertical_line = rho_theoretic,
              xticks = [2,3,4,5,6], yticks = None,
              xlabel = "rho", ylabel = "KS distance",
              savefig = False, fileName = fileName,
              logscale = False )



# =============================================================================
# TESTING DISTRIBUTIONAL LIMIT FOR AN SBM GRAPH
# =============================================================================


n = 200000
n_clusters = 2
alpha = [ 1 / n_clusters for i in range( n_clusters ) ]
block_sizes = [ int( alpha[i] * n ) for i in range( n_clusters ) ]

labels = []
for k in range( n_clusters ):
    labels += [ k for i in range( block_sizes[k] ) ] 

asymptotic_rate = np.log(n)**2
B = np.array( [ [1,2],[2,4] ] )

lambdain = 1
lambdaout = 1
weights_parameters = (lambdain - lambdaout) * np.eye( n_clusters ) + lambdaout * np.ones( n_clusters )


T = np.multiply( weights_parameters, B) @ np.diag( alpha ) * asymptotic_rate - np.ones( n_clusters )
vals, vecs = np.linalg.eigh(T.T)
rho_theoretic = np.max( vals ) / asymptotic_rate 

start = np.min( T @ np.ones( n_clusters ) ) / asymptotic_rate 
stop = np.max( T @ np.ones( n_clusters ) ) / asymptotic_rate 
rho_range = np.linspace( start = start, stop = stop, num=30 )

pi = vecs[:,1]
c = np.log( np.sum( [ pi[k]**2 / alpha[k] for k in range(n_clusters) ] ) )

n_trials = 10

------  
ks_distance_sbm = np.zeros( (n_trials, len(rho_range )) )

for trial in tqdm(range( n_trials )):
    G = ig.Graph.SBM( np.sum(block_sizes), B * asymptotic_rate / n, block_sizes, directed=False, loops=False)
    G = add_weights( G, weights_parameters, block_sizes, weight_distribution = 'exponential' )
    ks_distance_sbm[trial,:] = testingDistributionalLimit( G, rho_range * asymptotic_rate, n_samples = 100, source_vertices = [i for i in range(n)], target_vertices = [i for i in range(n)] , c =c ) 

ks_distance_sbm_mean = np.mean(ks_distance_sbm, axis=0)
ks_distance_sbm_std = np.std(ks_distance_sbm, axis=0)

fileName = 'KS_rho_SBM_n_' + str(n) + '_B_' + str(B) + '_alpha_' + str(alpha) + '.pdf'
plotFigure( rho_range, ks_distance_sbm_mean, accuracy_err = ks_distance_sbm_std, methods = None, 
           vertical_line = rho_theoretic,
              xticks = [4,4.5, 5.0, 5.5, 6], yticks = None,
              xlabel = "rho", ylabel = "KS distance",
              savefig = False, fileName = fileName )
--------
mean_limit_sbm = np.zeros( (n_trials, len(rho_range )) )

for trial in tqdm(range( n_trials )):
    G = ig.Graph.SBM( np.sum(block_sizes), B * asymptotic_rate / n, block_sizes, directed=False, loops=False)
    G = add_weights( G, weights_parameters, block_sizes, weight_distribution = 'exponential' )
    mean_limit_sbm[trial,:] = testingLimit( G, rho_range, n_samples = 100, source_vertices = [i for i in range(n)], target_vertices = [i for i in range(n)] ) 

mean_limit_sbm_mean = np.mean(mean_limit_sbm, axis=0)
mean_limit_sbm_std = np.std(mean_limit_sbm, axis=0)

fileName = 'finding_rho_SBM_n_' + str(n) + '_B_' + str(B) + '_alpha_' + str(alpha) + '.pdf'
plotFigure( rho_range, np.abs( mean_limit_sbm_mean ), accuracy_err = mean_limit_sbm_std, methods = None, 
           vertical_line = rho_theoretic,
              xticks = None, yticks = None,
              xlabel = "rho", ylabel = "Empirical average",
              savefig = False, fileName = fileName )
-----------


ks_distance_sbm = dict()
ks_distance_sbm['00'] = np.zeros( (n_trials, len(rho_range )) )
ks_distance_sbm['01'] = np.zeros( (n_trials, len(rho_range )) )
ks_distance_sbm['11'] = np.zeros( (n_trials, len(rho_range )) )

for trial in tqdm(range( n_trials )):
    G = ig.Graph.SBM( np.sum(block_sizes), B * asymptotic_rate / n, block_sizes, directed=False, loops=False)
    G = add_weights( G, weights_parameters, block_sizes, weight_distribution = 'exponential' )
    ks_distance_sbm['00'][trial,:] = testingDistributionalLimit( G, rho_range, n_samples = 50, source_vertices = [i for i in range(n) if labels[i]==0], target_vertices = [i for i in range(n) if labels[i]==0 ], c=c )
    ks_distance_sbm['01'][trial,:] = testingDistributionalLimit( G, rho_range, n_samples = 50, source_vertices = [i for i in range(n) if labels[i]==0], target_vertices = [i for i in range(n) if labels[i]==1], c=c )
    ks_distance_sbm['11'][trial,:] = testingDistributionalLimit( G, rho_range, n_samples = 50, source_vertices = [i for i in range(n) if labels[i]==1], target_vertices = [i for i in range(n) if labels[i]==1], c=c )

ks_distance_sbm_mean = dict( )
ks_distance_sbm_std = dict( )
for key in ks_distance_sbm.keys( ):
    ks_distance_sbm_mean[key] = np.mean(ks_distance_sbm[key], axis=0)
    ks_distance_sbm_std[key] = np.std(ks_distance_sbm[key], axis=0)


fileName = 'KS_rho_SBM_per_community_n_' + str(n) + '_B_' + str(B) + '_alpha_' + str(alpha) + '.pdf'
plotFigure( rho_range, ks_distance_sbm_mean, accuracy_err = ks_distance_sbm_std, methods = ks_distance_sbm_mean.keys( ), 
           vertical_line = rho_theoretic,
              xticks = None, yticks = None,
              xlabel = r"$ \rho $", ylabel = "KS distance",
              savefig = False, fileName = fileName )

-----------

mean_limit_sbm = dict()
mean_limit_sbm['00'] = np.zeros( (n_trials, len(rho_range )) )
mean_limit_sbm['01'] = np.zeros( (n_trials, len(rho_range )) )
mean_limit_sbm['11'] = np.zeros( (n_trials, len(rho_range )) )

for trial in tqdm(range( n_trials )):
    G = ig.Graph.SBM( np.sum(block_sizes), B * asymptotic_rate / n, block_sizes, directed=False, loops=False)
    G = add_weights( G, weights_parameters, block_sizes, weight_distribution = 'exponential' )
    mean_limit_sbm['00'][trial,:] = testingLimit( G, rho_range, n_samples = 50, source_vertices = [i for i in range(n) if labels[i]==0], target_vertices = [i for i in range(n) if labels[i]==0 ] )
    mean_limit_sbm['01'][trial,:] = testingLimit( G, rho_range, n_samples = 50, source_vertices = [i for i in range(n) if labels[i]==0], target_vertices = [i for i in range(n) if labels[i]==1] )
    mean_limit_sbm['11'][trial,:] = testingLimit( G, rho_range, n_samples = 50, source_vertices = [i for i in range(n) if labels[i]==1], target_vertices = [i for i in range(n) if labels[i]==1] )

mean_limit_sbm_mean = dict( )
mean_limit_sbm_std = dict( )
for key in ks_distance_sbm.keys( ):
    mean_limit_sbm_mean[key] = np.abs( np.mean(mean_limit_sbm[key], axis=0) )
    mean_limit_sbm_std[key] = np.std(mean_limit_sbm[key], axis=0)


fileName = 'finding_rho_SBM_per_community_n_' + str(n) + '_B_' + str(B) + '_alpha_' + str(alpha) + '.pdf'
plotFigure( rho_range, mean_limit_sbm_mean, accuracy_err = mean_limit_sbm_std, methods = ks_distance_sbm_mean.keys( ), 
           vertical_line = rho_theoretic,
              xticks = None, yticks = None,
              xlabel = r"$ \rho $", ylabel = r"$ | \rho \, C(u,v) - 1 | $",
              savefig = False, fileName = fileName )


# =============================================================================
# TESTING SEVA COMPUTATIONS
# =============================================================================


n_samples = 100
tau0 = [ ]
tau1 = [ ]

communities = [ [ i for i in range(n) if labels[i]==k] for k in range(n_clusters) ]
for i in tqdm( range( n_samples ) ):
    u = random.choice( communities[ 0 ] )
    tau0.append( q_nearest_neighbors_cost(G, u, int( np.sqrt(n) ) ) )
    v = random.choice( communities[ 1 ] )
    tau1.append( q_nearest_neighbors_cost(G, v, int( np.sqrt(n) ) ) )

tau0_mean = np.mean( tau0 )
tau1_mean = np.mean( tau1 )


-----------------
OLD STUFFS
n = 250000
n_clusters = 2
alpha = [ 1 / n_clusters for i in range( n_clusters ) ]
block_sizes = [ int( alpha[i] * n ) for i in range( n_clusters ) ]

labels = []
for k in range( n_clusters ):
    labels += [ k for i in range( block_sizes[k] ) ] 

asymptotic_rate = np.log(n) / n 
B = np.array( [ [4,4],[4,8] ] )

G = ig.Graph.SBM( np.sum(block_sizes), B * asymptotic_rate, block_sizes, directed=False, loops=False)
lambdain = 1
lambdaout = 1
weights_parameters = (lambdain - lambdaout) * np.eye( n_clusters ) + lambdaout * np.ones( n_clusters )
G = add_weights( G, weights_parameters, block_sizes, weight_distribution = 'exponential' )

T = np.multiply( weights_parameters, B) @ np.diag( alpha )
vals, vecs = np.linalg.eigh(T)
rho_theoretic = np.max( vals )

rho_range = np.linspace( start = np.min( T @ np.ones( n_clusters ) ), stop = np.max( T @ np.ones( n_clusters ) ), num=30 )

p_values_SBM = dict( )
n_samples = 200
p_values_SBM['all'] = testingDistributionalLimit( G, rho_range, n_samples, source_vertices = [i for i in range(n)], target_vertices = [i for i in range(n)], verbose = True ) 
p_values_SBM['0-0'] = testingDistributionalLimit( G, rho_range, n_samples, source_vertices = [i for i in range(n) if labels[i]==0], target_vertices = [i for i in range(n) if labels[i]==0], verbose = True ) 
p_values_SBM['0-1'] = testingDistributionalLimit( G, rho_range, n_samples, source_vertices = [i for i in range(n) if labels[i]==0], target_vertices = [i for i in range(n) if labels[i]==1], verbose = True ) 
p_values_SBM['1-1'] = testingDistributionalLimit( G, rho_range, n_samples, source_vertices = [i for i in range(n) if labels[i]==1], target_vertices = [i for i in range(n) if labels[i]==1], verbose = True ) 


fileName = 'pvalues_finding_rho_sbm_n_' + str(n) + '_B_' + str(B) + '_alpha_' + str(alpha) + '.pdf'
plotFigure( rho_range, p_values_SBM, accuracy_err = None, methods = p_values_SBM.keys( ), 
           vertical_line = rho_theoretic,
              xticks = [4,4.5,5,5.5,6], yticks = None,
              xlabel = "rho", ylabel = "p value",
              savefig = False, fileName = fileName,
              logscale = False )
"""


def q_nearest_neighbors_cost(G, u, q):
    """
    Returns the costs of the shortest paths from vertex u to its q nearest neighbors.
    
    Parameters:
        G (igraph.Graph): A weighted graph.
        u (int or str): The vertex index or name.
        q (int): Number of nearest neighbors to consider.
    
    Returns:
        list of tuples: Each tuple is (vertex, cost), sorted by increasing cost.
    """
    # Compute shortest path distances from u to all vertices
    distances = G.distances(source=u, weights="weight")[0]
    
    # Create list of (vertex, distance) excluding self
    vertex_dist = [(v, d) for v, d in enumerate(distances) if v != G.vs.find(u).index]
    
    # Filter out unreachable vertices (distance == inf)
    vertex_dist = [(v, d) for v, d in vertex_dist if d < float('inf')]
    
    # Sort by distance (ascending)
    vertex_dist.sort(key=lambda x: x[1])
    
    # Take the q nearest neighbors
    q_nearest = vertex_dist[:q]
    
    return q_nearest[-1][1]


import math
from mpmath import digamma

def sum_inverse_fraction(n):
    """
    Compute S = 3/25 * sum_{k=0}^{n-1} 1 / ((3/5 + k)*(6/5 + k))
    using both direct summation and closed form.
    """
    # Direct numerical sum
    direct_sum = sum(1 / ((3/5 + k)*(6/5 + k)) for k in range(n))
    
    # Closed form using digamma function
    a, b = 3/5, 6/5
    closed_form = (5/3) * ((digamma(n + a) - digamma(a)) - (digamma(n + b) - digamma(b)))
    
    return 3/25 * float(direct_sum), 3/25 * float(closed_form)



def testingNearestNeighborCost( G, asymptotic_rate, rho_range, n_samples, source_vertices, gamma_shape, verbose = False ):

    result = [ ]
    costs = [ ]
    n = G.vcount( )

    if verbose:
        iteration = tqdm( range( n_samples ))
    else:
        iteration = range( n_samples )
    for _ in iteration:
        u = random.choice( source_vertices )
        costs.append( q_nearest_neighbors_cost( G, u, q = int(np.sqrt(n)) ) )

    for rho in rho_range:
        logGamma = - np.log( sp.stats.gamma.rvs( size = 5000, a = gamma_shape ) )
        #gumbel = sp.stats.gumbel_r.rvs( size = 5000 )
        rescaled_costs = [ rho * asymptotic_rate * cost - np.log(n) / 2 for cost in costs ]
        result.append( sp.stats.ks_2samp(rescaled_costs, logGamma ).statistic )
    
    return np.asarray( result )



def samplingCosts( G, n_samples, source_vertices, target_vertices, verbose = True ):
    
    if verbose:
        iteration = tqdm( range( n_samples ) )
    else:
        iteration = range( n_samples )
    
    costs = [ ]
    for i in iteration:
        u = random.choice( source_vertices )
        v = random.choice( target_vertices )
        distances = G.distances( source = u, target = v, weights = 'weight' )
        costs.append( distances[0][0] )
    costs = np.asarray( costs, dtype='float64' )

    return costs


def testingDistributionalLimit( G, rho_range, n_samples, source_vertices, target_vertices, c, verbose = False, shift_source = 1, shift_target =1 ):

    result = [ ]
    costs = samplingCosts( G, n_samples, source_vertices, target_vertices, verbose = verbose )

    for rho in rho_range:
        #gumbels = sp.stats.gumbel_r.rvs( size = 5000 ) + sp.stats.gumbel_r.rvs( size = 5000 ) - sp.stats.gumbel_r.rvs( size = 5000 ) - c * np.ones( 5000 )
        limitDistribution = - np.log( sp.stats.gamma.rvs( size = 5000, a = shift_source ) ) - np.log( sp.stats.gamma.rvs( size = 5000, a = shift_target ) ) - sp.stats.gumbel_r.rvs( size = 5000 ) - c * np.ones( 5000 )
        # if shift_source = shift_target = 1, then limitDistribution is gumbel1 + gumbel2 - gumbel3 - c
        rescaled_costs = [  ( rho * cost - np.log( G.vcount( ) ) ) for cost in costs ]
        result.append( sp.stats.ks_2samp(rescaled_costs, limitDistribution ).statistic )
    
    return np.asarray( result )



def testingLimit( G, rho_range, n_samples, source_vertices, target_vertices, verbose = False ):

    result = [ ]
    costs = samplingCosts( G, n_samples, source_vertices, target_vertices, verbose = verbose )

    for rho in rho_range:
        rescaled_costs = [ ( rho * cost / np.log(n) - 1 ) for cost in costs ]
        result.append( np.mean( rescaled_costs ) )
    
    return np.asarray( result )



def estimateCostShortestPath( G, n_samples, source_vertices, target_vertices, 
                             verbose = True ):
    
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
        
        
        
        
# =============================================================================
# CODE TO PLOT FIGURES
# =============================================================================
SIZE_TITLE = 24
SIZE_LABELS = 24
SIZE_TICKS = 20
SIZE_LEGEND = 18

def plotFigure( x, accuracy_mean, accuracy_err = None, methods = None, 
               vertical_line = None, 
               xticks = None, yticks = None,
               xlabel = "x", ylabel = "Accuracy",
               savefig = False, fileName = "fig.pdf",
               logscale = False ):

    if logscale:
        fig, ax = plt.subplots( 1, 1  )
        ax.set_xscale("log")

    if methods is None and accuracy_err is None:
        plt.plot( x, accuracy_mean, linestyle = '-.', marker = '.' )
        
    elif methods is None and accuracy_err is not None:
        plt.errorbar( x, accuracy_mean, yerr = accuracy_err, linestyle = '-.' )
        
    elif methods is not None and accuracy_err is None:
        for method in methods:
            plt.plot( x, accuracy_mean[ method ], linestyle = '-.', marker = '.', label = method )
            legend = plt.legend( loc=0,  fancybox = True, fontsize = SIZE_LEGEND )
            plt.setp( legend.get_title(),fontsize = SIZE_LEGEND )
    
    elif methods is not None and accuracy_err is not None:
        for method in methods:
            plt.errorbar( x, accuracy_mean[ method ], yerr = accuracy_err[ method ], linestyle = '-.', label = method )
            legend = plt.legend( loc=0,  fancybox = True, fontsize = SIZE_LEGEND )
            plt.setp( legend.get_title(),fontsize = SIZE_LEGEND )
        
    if vertical_line is not None:
        plt.axvline( x = vertical_line )
    plt.xlabel( xlabel, fontsize = SIZE_LABELS )
    plt.ylabel( ylabel, fontsize = SIZE_LABELS )
    
    if xticks != None:
        plt.xticks( xticks, fontsize = SIZE_TICKS )
    else:
        plt.xticks( fontsize = SIZE_TICKS )
    
    if yticks != None:
        plt.yticks( yticks, fontsize = SIZE_TICKS )
    else:
        plt.yticks( fontsize = SIZE_TICKS )
    if(savefig):
        plt.savefig( fileName, bbox_inches='tight' )
    plt.show( )
