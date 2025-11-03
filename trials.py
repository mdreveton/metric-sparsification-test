#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  9 22:25:32 2025

@author: dreveton
"""

from trial_multiprocessing import makeUnionShortestPathTree_mp

import experiments_ErdosRenyi as er

import numpy as np
import igraph as ig


if __name__ == '__main__':
    
    n = 50000
    p = 2*np.log(n) / n
    G = ig.Graph.Erdos_Renyi( n, p )
    G = er.add_weights( G )
    source_vertices = [i for i in range(n)]

    SPT = makeUnionShortestPathTree_mp( G, source_vertices) # this would run multiprocessing code
