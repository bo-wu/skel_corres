#! /usr/bin/env python
import os, sys
import numpy as np
import openmesh as om
from graph_tool import Graph, search, topology
from sklearn.neighbors import NearestNeighbors

class DepthVisitor(search.DFSVisitor):
    """
    inherented from graph_tool for depth first visitor
    """
    def __init__(self):
        self.node_path = []

    def reset(self):
        self.node_path = []

    def get_short_branch(self, min_length=30):
        if len(self.node_path) < min_length:
            print 'short path', self.node_path
            return self.node_path
        else:
            return []

    def discover_vertex(self, v):
        if v.out_degree() > 2:
            raise search.StopSearch
        else:
            self.node_path.append(int(v))

    def examine_vertex(self, v):
        pass
    def non_tree_edge(self, e):
        pass
    def black_target(self, e):
        pass


class SkeletonData(object):
    """
    class to store and process skeleton data, like generated from starlab mean curvature skeleton
    """
    def __init__(self, fname=None, mesh_name=None, filter_sb=False):
        """
        @param filter_sb: if filter out Short Branch
        """
        if fname != None:
            self.skel_name = fname
            self.read_skel_file(fname)
            self._filter_short_branch(filter=filter_sb)
            self._parse_data()
            self.mesh_name = mesh_name

    def read_skel_file(self, fname, dim=3):
        if fname == None:
            print 'please input skeleton file name'
            sys.exit(0)
        elif os.path.isfile(fname):
            self.verts_init = []
            self.edges_init = []
            with open(fname) as sf:
                for line in sf:
                    line = line.strip('\n')
                    line = line.split(' ')
                    if line[0] == '#':
                        continue
                    elif line[0] == 'v':
                        self.verts_init.append([x for x in line[1:(dim+1)]])
                    #### attention!! verts of edge start from 1 in files ####
                    elif line[0] == 'e':
                        self.edges_init.append([int(x)-1 for x in line[1:3]])
                    else:
                        print 'not support this format'
                        sys.exit(0)
        else:
            print 'no such flie', fname
            sys.exit(0)

    def _filter_short_branch(self, filter=False):
        """
        filter out very short branches: do this maybe not right for some models, for models with flat part, it is right
        I will test how this effect the final matching results
        need to delete nodes, switch with the last one then delete last
        """
        if filter == False:
            self.verts = self.verts_init
            self.edges = self.edges_init
        else:
            init_graph = Graph(directed=False)
            init_graph.add_vertex(len(self.verts_init))
            for edge in self.edges_init:
                init_graph.add_edge(init_graph.vertex(edge[0]), init_graph.vertex(edge[1]))

            terminal_node = []
            for v in init_graph.vertices():
                if v.out_degree() == 1:
                    terminal_node.append(v)

            visitor = DepthVisitor()
            short_nodes = []
            for tn in terminal_node:
                search.dfs_search(init_graph, tn, visitor)
                tmp_node = visitor.get_short_branch()
                visitor.reset()
                for n in tmp_node:
                    short_nodes.append(n)

            ## get edges on the short paths
            short_nodes = list(set(short_nodes))
            short_edges = []
            for v in reversed(sorted(short_nodes)):
                for ve in init_graph.vertex(v).out_edges():
                    short_edges.append(ve)

            ## delete edges first, then vertex
            short_edges = list(set(short_edges))
            for e in short_edges:
                init_graph.remove_edge(e)

            temp_verts = self.verts_init[:]
            v_num = len(self.verts_init)
            print 'deleting vertex',
            for v in reversed(sorted(short_nodes)):
                print v,
                temp_verts[int(v)] = temp_verts[v_num-1]
                init_graph.remove_vertex(v, fast=True)
                v_num -= 1
            print 'deleting related edges' # already done above, just info user
            ######## new vertices and edges ########
            self.verts = temp_verts[:v_num]
            self.edges = []
            for e in init_graph.edges():
                self.edges.append([int(e.source()), int(e.target())])


    def _parse_data(self):
        """
        extract interal points(degree>2) and endpoints(degree=1)
        extract segments
        """
        if self.verts == None or self.edges == None:
            print 'please first call read_skel_file function'
        else:
            self.verts = np.array(self.verts, dtype=np.float)
            self.edges = np.array(self.edges, dtype=np.int)
            terminal_index = []
            junction_index = []
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
            self.terminal_index = terminal_index
            self.junction_index = junction_index
            self.feature_node_index = junction_index + terminal_index 

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
    
    def calc_edge_length(self):
        """
        calc edge length and make it edge property map in graph-tool
        """
        vec = self.verts[self.edges[:,0]] - self.verts[self.edges[:,1]]
        edge_length = np.sqrt(np.sum(vec**2, axis=-1))
        self.edge_length_map = self.skel_graph.new_edge_property("double")
        self.edge_length_map.a = edge_length
    
    def calc_node_centricity(self):
        """
        calc node centricity of feature nodes(terminal and junction nodes)
        T1 in Oscar's EG 2010 paper
        """
        node_centricity = []
        max_length = []
        for n_idx in self.feature_node_index:
            dist = topology.shortest_distance(self.skel_graph, self.skel_graph.vertex(n_idx))
            node_centricity.append(dist.a.mean())
            max_length.append(float(dist.a.max()))
        node_centricity = np.array(node_centricity)
        max_length = np.array(max_length)
        self.node_centricity = node_centricity / np.max(max_length)

    def calc_skel_radius(self, dim=3):
        if self.mesh_name == None:
            print 'please set mesh_name before calc_skel_radius'
        elif os.path.isfile(self.mesh_name):
            mesh = om.TriMesh()
            assert om.read_mesh(mesh, self.mesh_name)
            mesh_vertices = np.array((mesh.n_vertices, dim), dtype=float)
            for vh in mesh.vertices():
                

    def write_file(self, file_path='./'):
        """
        maybe need to save file after filter
        same as starlab mean curvature skeleotn
        """
        file_name = os.path.basename(self.skel_name)
        full_name = file_path + file_name
        v_num = len(self.verts)
        e_num = len(self.edges)
        first_line = '# D:3 ' + 'NV:' + str(v_num) + ' NE:' + str(e_num) + '\n'
        with open(full_name, 'w') as f:
            f.write(first_line)
            for v in self.verts:
                line = 'v ' + str(v[0]) + ' ' + str(v[1]) + ' ' + str(v[2]) + '\n'
                f.write(line)

            for e in self.edges:
                line = 'e ' + str(e[0]+1) + ' ' + str(e[1]+1) + '\n'
                f.write(line)


if __name__ == '__main__':
    skeleton = SkeletonData('data/chair_skeleton/115_ckel.cg')
    skeleton.write_file()
    #print 'out degree==2', skeleton.skel_graph.vertex(5).out_degree() == 2

