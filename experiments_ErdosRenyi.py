#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 22 11:41:16 2025

@author: dreveton
"""


import igraph as ig
import numpy as np
import random as random
import scipy as sp
from scipy.stats import binomtest

from tqdm import tqdm

import matplotlib.pyplot as plt



SIZE_TITLE = 24
SIZE_LABELS = 24
SIZE_TICKS = 20
SIZE_LEGEND = 18


#import math as m
#EULER_MASCHERONI_CONSTANT = round( (1.-m.gamma(1+1.e-8))*1.e14 )*1.e-6
EULER_MASCHERONI_CONSTANT = 0.5772156649


"""

# =============================================================================
# CRITICAL ERDOS-RENYI
# =============================================================================

#EXPERIMENT 1: DIFFERENT WEIGHT DISTRIBUTIONS WITH NON-ZERO PDF AT 0


n = 250000
n_graphs = 10
n_samples = 100

#epsilon_range = [ 0.01, 0.1, 1, 10 ]

epsilon_range = np.logspace( -2, np.log10(2), num=10 )


weight_distributions = [ 'exponential', 'uniform', 'gaussian', 'cauchy' ]

B_mb = dict( )
confidence_intervals = dict( )

for dummy in tqdm( range( len( weight_distributions ) ) ):
    weight_distribution = weight_distributions[ dummy ]
    B_mb[ weight_distribution ], confidence_intervals[ weight_distribution ] = sparsification_criticalErdosRenyi_varyingEspilon( n, n_samples, n_graphs, epsilon_range, weight_distribution = weight_distribution, verbose = False )
    

plotFigure( epsilon_range, B_mb, curve_err = confidence_intervals, 
           methods = weight_distributions, 
               xticks = [0,0.5,1,1.5,2], yticks = None,
               xlabel = "$epsilon$", ylabel = "$B_{mb}$",
               savefig = False, 
               fileName = "fig.pdf" )


filename = 'n_' + str(n) + '_n_graphs_' + str(n_graphs) + '_n_samples_' + str(n_samples) + '.pdf'
fig, ax = plt.subplots(1, 1  )
ax.set_xscale("log")
for weight_distribution in weight_distributions:
    ax.errorbar(epsilon_range, B_mb[ weight_distribution ], yerr = confidence_intervals[ weight_distribution ], linestyle = '-.', label = weight_distribution )
legend = plt.legend( loc=0,  fancybox = True, fontsize = SIZE_LEGEND )
plt.setp( legend.get_title(),fontsize = SIZE_LEGEND )
plt.xlabel( '$epsilon$', fontsize = SIZE_LABELS )
plt.ylabel( '$B^{mb}$', fontsize = SIZE_LABELS )
plt.xticks( fontsize = SIZE_TICKS )
plt.yticks( fontsize = SIZE_TICKS )
plt.savefig( filename, bbox_inches='tight' )
plt.show()





#EXPERIMENT 2: WEIGHT DISTRIBUTION SUCH THAT THE EDGE-DENSITY OF THE METRIC BACKBONE IS ARBITRARY (BETWEEN 0 AND 1 TIMES LOG N / N)


n = 250000
n_graphs = 10
n_samples = 100

epsilon_range = np.logspace( -2, np.log10(2), num=10 )

weight_distribution = 'made_up_mixture'

weightDistribution_parameters = [0.2, 0.4, 0.6 ]

B_mb = dict( )
confidence_intervals = dict( )

for dummy in tqdm( range( len( weightDistribution_parameters ) ) ):
    parameter = weightDistribution_parameters[ dummy ]
    B_mb[ parameter ], confidence_intervals[ parameter ] = sparsification_criticalErdosRenyi_varyingEspilon( n, n_samples, n_graphs, epsilon_range, weight_distribution = weight_distribution, weightDistribution_parameter = [parameter], verbose = False )
    

filename = 'made_up_n_' + str(n) + '_n_graphs_' + str(n_graphs) + '_n_samples_' + str(n_samples) + '.pdf'

plotFigure( epsilon_range, B_mb, curve_err = confidence_intervals, 
           methods = weightDistribution_parameters, 
               xticks = [0,0.5,1,1.5,2], yticks = None,
               xlabel = "$epsilon$", ylabel = "$B_{mb}$",
               savefig = False, 
               fileName = filename )

filename = 'made_up_n_' + str(n) + '_n_graphs_' + str(n_graphs) + '_n_samples_' + str(n_samples) + 'logscale.pdf'
fig, ax = plt.subplots(1, 1  )
ax.set_xscale("log")
for parameter in weightDistribution_parameters:
    ax.errorbar(epsilon_range, B_mb[ parameter ], yerr = confidence_intervals[ parameter ], linestyle = '-.', label = parameter )
legend = plt.legend( loc=0,  fancybox = True, fontsize = SIZE_LEGEND )
plt.setp( legend.get_title(),fontsize = SIZE_LEGEND )
plt.xlabel( '$epsilon$', fontsize = SIZE_LABELS )
plt.ylabel( '$B^{mb}$', fontsize = SIZE_LABELS )
plt.xticks( fontsize = SIZE_TICKS )
plt.yticks( fontsize = SIZE_TICKS )
plt.savefig( filename, bbox_inches='tight' )
plt.show()



# =============================================================================
# LIMIT OF THE SHORTEST PATH COSTS
# =============================================================================


n = 200000
n_graphs = 1
n_samples = 100

epsilon_range = np.logspace( -2, np.log10(2), num=10 )

weight_distributions = [ 'exponential' ]
weightDistribution_parameters = [ 1 ]

rescaled_cost_mb = dict( )
confidence_intervals = dict( )

for dummy in tqdm( range( len( weight_distributions ) ) ):
    weight_distribution = weight_distributions[ dummy ]
    rescaled_cost_mb[ weight_distribution ], confidence_intervals[ weight_distribution ] = cost_criticalErdosRenyi_varyingEspilon( n, n_samples, n_graphs, epsilon_range, 
                                               weight_distribution = weight_distribution, weightDistribution_parameter = weightDistribution_parameters, 
                                               verbose = False, rescaling = 'exact' )

filename = 'made_up_n_' + str(n) + '_n_graphs_' + str(n_graphs) + '_n_samples_' + str(n_samples) + '.pdf'

plotFigure( epsilon_range, rescaled_cost_mb, curve_err = confidence_intervals, 
           methods = weight_distributions, 
               xticks = [0,0.5,1,1.5,2], yticks = None,
               xlabel = "$epsilon$", ylabel = "Rescaled cost",
               savefig = False, 
               fileName = filename )


----------------------

n_range = [ 100000, 150000, 200000 ]
n_graphs = 1
n_samples = 100

epsilon = 2 

weight_distributions = [ 'uniform' ]
weightDistribution_parameters = [ 1 ]

rescaled_cost_mb = dict( )
confidence_intervals = dict( )

for dummy in tqdm( range( len( weight_distributions ) ) ):
    weight_distribution = weight_distributions[ dummy ]
    rescaled_cost_mb[ weight_distribution ], confidence_intervals[ weight_distribution ] = cost_criticalErdosRenyi_varyingNumberVertices( n_range, n_samples, n_graphs, epsilon, 
                                               weight_distribution = weight_distribution, weightDistribution_parameter = weightDistribution_parameters, 
                                               verbose = False, rescaling = 'exact' )
    

filename = 'ER_cost_epsilon_' + str(epsilon) + '_n_graphs_' + str(n_graphs) + '_n_samples_' + str(n_samples) + '.pdf'

plotFigure( n_range, rescaled_cost_mb, curve_err = confidence_intervals, 
           methods = weight_distributions, 
               xticks = n_range, yticks = None,
               xlabel = "$n$", ylabel = "Rescaled cost",
               savefig = False, 
               fileName = filename )

----------------------

n = 200000
epsilon = 2
n_samples = 100
n_pvalues = 10

weight_distributions = [ 'exponential', 'uniform', 'gaussian', 'cauchy' ]
weightDistribution_parameter = [1]

p = (1+epsilon) * np.log(n) / n

pvalues = dict( )

for weight_distribution in weight_distributions:
    
    G = ig.Graph.Erdos_Renyi( n, p )
    G = add_weights( G, weight_distribution = weight_distribution, parameter = weightDistribution_parameter )
    result = [ ]
    
    for _ in tqdm( range(n_pvalues) ):
        vertex_pairs = []
        for _ in range( n_samples ):
            u = random.randrange( 0, G.vcount() )
            v = random.randrange( 0, G.vcount() )
            vertex_pairs.append( (u,v) )
        
        costs = shortestPathCost( G, vertex_pairs, rescaling = None, verbose = True )
        rescaled_costs = np.log(n) * ( (1+epsilon) * _get_pdf_at_0( weight_distribution, weightDistribution_parameter ) * np.asarray( costs ) - 1.0 )
        
        gumbels = sp.stats.gumbel_r.rvs( size = 1000 ) + sp.stats.gumbel_r.rvs( size = 1000 ) - sp.stats.gumbel_r.rvs( size = 1000 )
        result.append( sp.stats.ks_2samp(rescaled_costs, gumbels ).pvalue )
    
    pvalues[ weight_distribution ] = result



to_plot = [ pvalues[ weight_distribution ] for weight_distribution in weight_distributions ]

fig, ax = plt.subplots()
ax.set_ylabel('p values')

bplot = ax.boxplot(to_plot,
                   tick_labels=weight_distributions)  # will be used to label x-ticks

plt.show( )


"""


"""
OLD STUFFS
-----------------------------------
n = 500000


probabilities_edge_kept = dict( )
confidence_interval = dict( )

n_samples = 500
n_graphs = 1
epsilon = 0.1
p = (1+epsilon) * np.log(n) / n

G = ig.Graph.Erdos_Renyi( n, p )
G = exponential_weights( G, exponential_rate = 1 ) 
#G = lognormal_weights( G ) 
#G = normal_weights( G ) 

edge_list = G.get_edgelist( )



method = 0
success = 0

for i in tqdm( range( n_samples ) ):
    (u,v) = edge_list[ i ]
    distances = G.distances( source = u, target = v, weights = 'weight' )
    if distances[0][0] == G[u,v]:
        success += 1
        
ratio = success / n_samples
error = binomtest( k = success, n = n_samples, p = ratio ).proportion_ci()







success_0 = 0
success_1 = 0
mismatch_01 = [ ]
mismatch_10 = [ ]


graph = sp.sparse.csr_array( G.get_adjacency_sparse( attribute = 'weight' ) )

for i in tqdm( range( n_samples ) ):
    (u,v) = edge_list[ i ]
    distances = G.distances( source = u, target = v, weights='weight' )
    if distances[0][0] == G[u,v]:
        success_0 += 1
        
    dist_matrix = sp.sparse.csgraph.dijkstra( csgraph = graph, directed = False, indices = u, return_predecessors = False, limit = G[u,v] )
    
    if dist_matrix[v] == G[u,v]:
        success_1 += 1
        
    if distances[0][0] == G[u,v] and not dist_matrix[v] == G[u,v]:
           mismatch_01.append( (u,v) )

    if dist_matrix[v] == G[u,v] and not distances[0][0] == G[u,v] :
           mismatch_10.append( (u,v) )

        
ratio = success / n_samples



-----------------

success = 0
n_samples = 500

for iter in tqdm( range( n_samples ) ):
    (u,v) = edge_list[ np.random.randint(len(edge_list)) ]
    distances = G.distances( source = u, target = v, weights='weight' )
    if distances[0][0] == G[u,v]:
        success += 1
        
ratio = success / n_samples




# =============================================================================
# CORE-PERIPHERY: EXTREMAL CASE WHERE CORE = CLIQUE
# =============================================================================

n = 250000
beta = 2/3 
p = 2 * np.log(n) / n
clique_size = int( n**(beta) )

G = plantedClique( n, p, clique_size )
G = add_weights( G, weight_distribution = 'exponential', parameter = [1] )

costs = [ ]

n_samples = 100 
for i in tqdm( range(n_samples) ):
    u = random.randrange( 0, clique_size )
    #v = random.randrange( 0, clique_size )
    v = random.randrange( clique_size, G.vcount() )

    costs.append( G.distances( source = u, target = v, weights = 'weight' )[0][0] )

print( np.mean( costs) )

(1-beta) * np.log(n) / (n*p) + beta * np.log(n) / (n**(beta))

"""


def shortestPathCost( G, vertex_pairs, rescaling = None, verbose = False ):
    
    if verbose:
        iteration = tqdm( range( len( vertex_pairs ) ) )
    else:
        iteration = range( len( vertex_pairs ) )

    costs = []
    
    for i in iteration:
        (u,v) = vertex_pairs[ i ]
        distances = G.distances( source = u, target = v, weights = 'weight' )
        
        if isinstance(rescaling, (int,float)):
            costs.append( rescaling * distances[0][0] )
        
        else:
            costs.append( distances[0][0] )

    return costs


def cost_criticalErdosRenyi_varyingNumberVertices( n_range, n_samples, n_graphs, epsilon, 
                                           weight_distribution = 'exponential', weightDistribution_parameter = [1], 
                                           verbose = False, rescaled = False, rescaling = 'estimated' ):
    
    cost_mean = np.zeros( len(n_range) )
    cost_std = np.zeros( len(n_range) )

    if verbose:
        iterations_epsilon = tqdm( range( len(n_range) ) )
    else:
        iterations_epsilon = range( len(n_range) )
    
    
    for dummy in iterations_epsilon:
        n = n_range[ dummy ]
        p = (1+epsilon) * np.log(n) / n
        cost = [ ] 
        for _ in range(n_graphs):
            G = ig.Graph.Erdos_Renyi( n, p )
            G = add_weights( G, weight_distribution = weight_distribution, parameter = weightDistribution_parameter )
            
            #edges_to_use = random.sample( G.get_edgelist( ), 100 )
            
            for i in range( n_samples ):
                #(u,v) = edges_to_use[ i ]
                u = random.randrange( 0, G.vcount() )
                v = random.randrange( 0, G.vcount() )
                distances = G.distances( source = u, target = v, weights = 'weight' )
                cost.append( distances[0][0] )            

        cost = np.asarray( cost, dtype='float64' )
        if isinstance(rescaling, (int, float) ):
            cost = rescaling * cost
            
        elif rescaling == True:
            pdf_at_0 = _get_pdf_at_0( weight_distribution, weightDistribution_parameter )
            estimated_rescaling = G.density() * pdf_at_0 / ( np.log(n) / n )
            #print( 'The estimated rescaling factor is equal to : ', estimated_rescaling )
            cost = estimated_rescaling * cost
            
        elif rescaling == 'exact':
            pdf_at_0 = _get_pdf_at_0( weight_distribution, weightDistribution_parameter )
            exact_rescaling = (1+epsilon) * pdf_at_0
            #print( 'The estimated rescaling factor is equal to : ', estimated_rescaling )
            cost = exact_rescaling * cost
        
        cost_mean[ dummy ] = np.mean( cost )
        cost_std[ dummy ] = np.std( cost ) / np.sqrt( n_samples * n_graphs )
        #print( epsilon )
        
    return cost_mean, cost_std




def cost_criticalErdosRenyi_varyingEspilon( n, n_samples, n_graphs, epsilon_range, 
                                           weight_distribution = 'exponential', weightDistribution_parameter = [1], 
                                           verbose = False, rescaled = False, rescaling = 'estimated' ):
    
    cost_mean = np.zeros( len(epsilon_range) )
    cost_std = np.zeros( len(epsilon_range) )

    if verbose:
        iterations_epsilon = tqdm( range( len(epsilon_range ) ) )
    else:
        iterations_epsilon = range( len(epsilon_range ) )
    
    
    for dummy in iterations_epsilon:
        epsilon = epsilon_range[ dummy ]
        p = (1+epsilon) * np.log(n) / n
        cost = [ ] 
        for _ in range(n_graphs):
            G = ig.Graph.Erdos_Renyi( n, p )
            G = add_weights( G, weight_distribution = weight_distribution, parameter = weightDistribution_parameter )
            
            #edges_to_use = random.sample( G.get_edgelist( ), 100 )
            
            for i in range( n_samples ):
                #(u,v) = edges_to_use[ i ]
                u = random.randrange( 0, G.vcount() )
                v = random.randrange( 0, G.vcount() )
                distances = G.distances( source = u, target = v, weights = 'weight' )
                cost.append( distances[0][0] )            

        cost = np.asarray( cost, dtype='float64' )
        if isinstance(rescaling, (int, float) ):
            cost = rescaling * cost
            
        elif rescaling == True:
            pdf_at_0 = _get_pdf_at_0( weight_distribution, weightDistribution_parameter )
            estimated_rescaling = G.density() * pdf_at_0 / ( np.log(n) / n )
            if verbose:
                print( 'The estimated rescaling factor is equal to : ', estimated_rescaling )
            cost = estimated_rescaling * cost
            
        elif rescaling == 'exact':
            pdf_at_0 = _get_pdf_at_0( weight_distribution, weightDistribution_parameter )
            exact_rescaling = (1+epsilon) * pdf_at_0
            if verbose:
                print( 'The estimated rescaling factor is equal to : ', estimated_rescaling )
            cost = exact_rescaling * cost -  EULER_MASCHERONI_CONSTANT / np.log(n)
        
        cost_mean[ dummy ] = np.mean( cost )
        cost_std[ dummy ] = np.std( cost ) / np.sqrt( n_samples * n_graphs )
        if verbose:
            print( epsilon )
        
    return cost_mean, cost_std
    



def sparsification_criticalErdosRenyi_varyingEspilon( n, n_samples, n_graphs, epsilon_range, weight_distribution = 'exponential', weightDistribution_parameter = [], verbose = False ):
    
    B_mb = np.zeros( len(epsilon_range) )
    confidence_interval = np.zeros( (2,len(epsilon_range) ) )

    if verbose:
        iterations_epsilon = tqdm( range( len(epsilon_range ) ) )
    else:
        iterations_epsilon = range( len(epsilon_range ) )
        
    for dummy in iterations_epsilon:
        epsilon = epsilon_range[ dummy ]
        p = (1+epsilon) * np.log(n) / n
        
        success = 0
        
        for _ in range(n_graphs):
            G = ig.Graph.Erdos_Renyi( n, p )
            G = add_weights( G, weight_distribution = weight_distribution, parameter = weightDistribution_parameter )
            success += estimate_probability_edge_kept( G, n_samples, return_number_success = True, verbose = verbose )
        
        ratio = success / ( n_samples * n_graphs )
        error = binomtest( k = success, n = n_samples * n_graphs, p = ratio ).proportion_ci()
        
        B_mb[ dummy ] = ratio * (1+epsilon)
        confidence_interval[ 0, dummy ] = np.abs( error[0] - ratio ) * (1+epsilon)
        confidence_interval[ 1, dummy ] = ( error[1] - ratio ) * (1+epsilon)
        print( epsilon )
        if error[1] - ratio < 0:
            print('This quantity should not be negative') #This error message should never happens, I am just making sure. 
        
    return B_mb, confidence_interval

        
def estimate_probability_edge_kept( G, n_samples, return_number_success = False, verbose = False ):
    
    edges_to_use = random.sample( G.get_edgelist( ), 100 )
    success = 0
    
    if verbose:
        iterations_samples = tqdm( range( n_samples ) )
    else:
        iterations_samples = range( n_samples )
        
    for i in iterations_samples:
        (u,v) = edges_to_use[ i ]
        distances = G.distances( source = u, target = v, weights = 'weight' )
        if distances[0][0] == G[u,v]:
            success += 1
    
    if return_number_success:
        return success
    else:
        return success / n_samples
    

def plantedClique( n, p, clique_size ):
    G2 = ig.Graph.Erdos_Renyi( n, p )
    G1 = ig.Graph.Full( n = clique_size )
    
    return ig.Graph.union(G1,G2)



def add_weights( G, weight_distribution = 'exponential', parameter = [1] ):
    rng = np.random.default_rng()
    
    if weight_distribution == 'exponential':
        weights = rng.exponential( scale = parameter[ 0 ], size = G.ecount() )
        
    elif weight_distribution == 'lognormal':
        weights = rng.lognormal( mean = 0, sigma = 1, size = G.ecount() )
        
    elif weight_distribution == 'normal' or weight_distribution.lower() == 'gaussian':
        weights = np.abs( rng.standard_normal( size = G.ecount() ) )
    
    elif weight_distribution.lower() == 'cauchy':
        weights = np.abs( rng.standard_cauchy( size = G.ecount() ) )
        
    elif weight_distribution.lower() == 'gamma':
        weights = np.abs( rng.standard_gamma( shape = 1, size = G.ecount() ) )

    elif weight_distribution.lower() == 'uniform':
        weights = np.abs( rng.uniform( low=0, high = parameter[0], size = G.ecount() ) )
        
    elif weight_distribution == 'made_up_mixture':
        nu = parameter[ 0 ]
        sizes = rng.multinomial( G.ecount(), [nu, 1-nu] )
        unif_weights = rng.random( size = sizes[0] ) * nu
        other_distrib_weights = np.ones( sizes[1]) + np.random.default_rng().pareto( 1.5, size = sizes[1] )
        weights = np.concat( (unif_weights, other_distrib_weights ) )
        rng.shuffle(weights)
    else:
        raise TypeError('This weight distribution is not implemented')

    G.es['weight'] = weights
    return G


def _get_pdf_at_0( weight_distribution, weightDistribution_parameter ):
    
    if weight_distribution == 'exponential':
        return weightDistribution_parameter[ 0 ]
    
    elif weight_distribution.lower() == 'cauchy':
        return 2 / ( np.pi )
    
    elif weight_distribution.lower() == 'normal' or weight_distribution.lower() == 'gaussian':
        return 2 / ( np.sqrt( 2 * np.pi ) )

    elif weight_distribution.lower() == 'uniform':
        return 1 / weightDistribution_parameter[ 0 ]

    else:
        return TypeError( 'Not implemented yet' )


def _exponential_weights( G, exponential_rate = 1 ):
    rng = np.random.default_rng()
    weights = rng.exponential(scale= 1 / exponential_rate, size = G.ecount() )
    G.es['weight'] = weights
    return G


def _lognormal_weights( G, mean = 0, sigma = 1 ):
    rng = np.random.default_rng()
    weights = rng.lognormal( mean = mean, sigma = sigma, size = G.ecount() )
    G.es['weight'] = weights
    return G


def _normal_weights( G, mean = 0, sigma = 1 ):
    rng = np.random.default_rng()
    weights = np.abs( rng.normal( loc = mean, scale = sigma, size = G.ecount() ) )
    G.es['weight'] = weights
    return G


def plotFigure( x, curve, curve_err = None, methods = None, 
               xticks = None, yticks = None,
               xlabel = "x", ylabel = "Accuracy",
               savefig = False, fileName = "fig.pdf" ):
    
    if methods is None and curve_err is None:
        plt.plot( x, curve, linestyle = '-.', marker = '.' )
        
    elif methods is None and curve_err is not None:
        plt.errorbar( x, curve, yerr = curve_err, linestyle = '-.' )
        
    elif methods is not None and curve_err is None:
        for method in methods:
            plt.plot( x, curve[ method ], linestyle = '-.', marker = '.', label = method )
            legend = plt.legend( loc=0,  fancybox = True, fontsize = SIZE_LEGEND )
            plt.setp( legend.get_title(),fontsize = SIZE_LEGEND )
    
    elif methods is not None and curve_err is not None:
        for method in methods:
            plt.errorbar( x, curve[ method ], yerr = curve_err[ method ], linestyle = '-.', label = method )
            legend = plt.legend( loc=0,  fancybox = True, fontsize = SIZE_LEGEND )
            plt.setp( legend.get_title(),fontsize = SIZE_LEGEND )

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

