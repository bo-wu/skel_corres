#! /usr/bin/env python
import os, sys
import numpy as np
from graph_tool import Graph, search

class DepthVisitor(search.DFSVisitor):
    """
    inherented from graph_tool for depth first visitor
    """
    def __init__(self, name, pred, time):
        self.name = name
        self.pred = pred
        self.time = time


class SkeletonData(object):
    """
    class to store and process skeleton data, like generated from starlab mean curvature skeleton
    """
    def __init__(self, fname):
        if fname != None:
            self.fname = fname
            self.read_file(fname)
            self._parse_data()

    def read_file(self, fname, dim=3):
        if fname == None:
            print 'please input skeleton file name'
            sys.exit(0)
        elif os.path.isfile(fname):
            self.verts = []
            self.edges = []
            with open(fname) as sf:
                for line in sf:
                    line = line.strip('\n')
                    line = line.split(' ')
                    if line[0] == '#':
                        continue
                    elif line[0] == 'v':
                        self.verts.append([x for x in line[1:(dim+1)]])
                    elif line[0] == 'e':
                        self.edges.append([x for x in line[1:3]])
                    else:
                        print 'not support this format'
                        sys.exit(0)
        else:
            print 'no such flie', fname
            sys.exit(0)

    def _parse_data(self):
        """
        extract interal points(degree>2) and endpoints(degree=1)
        extract segments
        """
        if self.verts == None or self.edges == None:
            print 'please first call read_file function'
        else:
            self.verts = np.array(self.verts, dtype=np.float)
            self.edges = np.array(self.edges, dtype=np.int)
            terminal_index = []
            junction_index = []
            ######### attention verts of edge start from 1 #########
            self.edges -= 1
            #######################################################
            self.skel_graph = Graph(directed=False)
            self.skel_graph.add_vertex(len(self.verts))
            for edge in self.edges :
                self.skel_graph.add_edge(self.skel_graph.vertex(edge[0]), self.skel_graph.vertex(edge[1]))

            for v in self.skel_graph.vertices():
                if v.out_degree() == 2 :
                    continue
                elif v.out_degree() == 1 :
                    terminal_index.append(int(v))
                elif v.out_degree() > 2 :
                    junction_index.append(int(v))

            self.terminal = self.verts[terminal_index]
            self.junction = self.verts[junction_index]

            """
            edge_vert_index = self.edges.flatten()
            print 'edge vertex index dtype', edge_vert_index.dtype
            if 0 in edge_vert_index:
                print 'vertex start from 0'
            else:
                print 'vertex start from 1'
            print 'skeleton vertex num', self.skel_graph.num_vertices()
            print 'skeleton edge num', self.skel_graph.num_edges()
            """
    
    def filter_short_branch(self):
        """
        filter out very short branches
        need to delete nodes, switch with the last one
        """
        pass

    def write_file(self, file_path='./'):
        """
        maybe need to save file after filter
        """
        full_name = file_path + self.fname
        with open(full_name, 'w') as f:
            pass
        pass


if __name__ == '__main__':
    skeleton = SkeletonData('data/chair_skeleton/115_ckel.cg')
    #print 'out degree==2', skeleton.skel_graph.vertex(5).out_degree() == 2
