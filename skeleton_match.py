#! /usr/bin/env python

class SkeletonMatch(object):
"""
implement Oscar's skeleton matching alogrithm in this class
"""
def __init__(self, skel1, skel2):
    if skel1 != None and skel2 != None :
        self.skel1 = skel1
        self.skel2 = skel2
    else:
        print 'need input two skeleton to match'
