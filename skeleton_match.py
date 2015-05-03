#! /usr/bin/env python
from graph_tool import Graph, util, search, topology
import numpy as np
import itertools


class SkeletonMatch(object):
    """
    implement Oscar's skeleton matching alogrithm in this class
    """
    def __init__(self, skel1, skel2):
        if skel1 != None and skel2 != None :
            self.skel1 = skel1
            self.skel2 = skel2
            self.skel1.calc_node_centricity()
            self.skel2.calc_node_centricity()
            self.skel1.calc_skel_radius()
            self.skel2.calc_skel_radius()
            self.skel1.calc_path_length_ratio()
            self.skel2.calc_path_length_ratio()
            self.skel1.calc_path_radius_ratio()
            self.skel2.calc_path_radius_ratio()
            # use index instead of real value
            skel1_index = np.arange(len(self.skel1.feature_node_index))
            skel2_index = np.arange(len(self.skel2.feature_node_index))
            junc1_num = len(skel1.junction_index)
            junc2_num = len(skel2.junction_index)
            
            #candidate matched pairs
            junction_pairs = []
            junc_term_pairs = []
            terminal_pairs = []
            for i, j in itertools.product(skel1_index, skel2_index):
                if self.match_node_centricity(c1=i, c2=j, threhold=.5):
                    if i < junc1_num and j < junc2_num: # only junction nodes
                        junction_pairs.append([i,j])
                    elif i >= junc1_num and j >= junc2_num: # with junction nodes
                        terminal_pairs.append([i,j])
                    else:
                        junc_term_pairs.append([i,j])

            self.junction_pairs = np.array(junction_pairs)
            self.terminal_pairs = np.array(terminal_pairs)
            self.junc_term_pairs = np.array(junc_term_pairs)
            self.all_junc_pairs = np.vstack((self.junction_pairs, self.junc_term_pairs))

            self.vote_tree = Graph(directed=False)
            self.node_pair = self.vote_tree.new_vertex_property("vector<short>")

            self._construct_search_tree()
        else:
            print 'need input two skeleton to match'


    def _construct_search_tree(self, prev_pairs=np.array([]), junc_pairs=np.array([]), term_pairs=np.array([])):
        """
        recursively consturct search tree
        @param prev_pairs record that already on the path
        @param junc_pairs record that left part junction pairs (on current tree level)
        @param term_pairs record that left terminal pairs on current tree level
        """
        # root of the tree
        if len(prev_pairs) == 0:
            v1 = self.vote_tree.add_vertex()
            #first level, only junction pairs (both are junction)
            for n, pair in enumerate(self.junction_pairs):
                new_prev = pair.reshape(-1,2)  # to use len(for level one), need to change shape
                v2 = self._construct_search_tree(prev_pairs=new_prev)
                self.vote_tree.add_edge(v1, v2)
            return v1

        elif len(prev_pairs) == 1:  # first level
            v1 = self.vote_tree.add_vertex()
            self.node_pair[v1] = prev_pairs[-1,:]

            check_junc = True
           # curr_junc = self.all_junc_pairs.copy()
           # counter = 0
            for n, pair in enumerate(self.all_junc_pairs):
                if pair[0] not in prev_pairs[:,0] and pair[1] not in prev_pairs[:,1]:
                    new_prev = np.vstack((prev_pairs, pair))
                    check_junc = False
                    v2 = self._construct_search_tree(prev_pairs=new_prev)
                    if v2 != None:
                        self.vote_tree.add_edge(v1, v2)
           #     else:
           #         curr_junc = np.delete(curr_junc, n-counter, 0)
           #         counter += 1

            if check_junc:
                for n, pair in enumerate(self.terminal_pairs):
                    new_prev = np.vstack((prev_pairs, pair))
                    v2 = self._construct_search_tree(prev_pairs=new_prev)
                    if v2 != None:
                        self.vote_tree.add_edge(v1, v2)

            return v1

        elif len(prev_pairs) > 1:  # above level two
            if self.match_length_radius(n1=prev_pairs[-1,0], n2=prev_pairs[-1,1], matched_pairs=prev_pairs[:-1]):
                v1 = self.vote_tree.add_vertex()
                self.node_pair[v1] = prev_pairs[-1,:]

                check_junc = True
                for pair in self.all_junc_pairs:
                    if pair[0] not in prev_pairs[:,0] and pair[1] not in prev_pairs[:,1]:
                        new_prev = np.vstack((prev_pairs, pair))
                        check_junc = False
                        v2 = self._construct_search_tree(prev_pairs=new_prev)
                        if v2 != None:
                            self.vote_tree.add_edge(v1, v2)

                if check_junc:
                    for pair in self.terminal_pairs:
                        if pair[0] not in prev_pairs[:,0] and pair[1] not in prev_pairs[:,1]:
                            new_prev = np.vstack((prev_pairs, pair))
                            v2 = self._construct_search_tree(prev_pairs=new_prev)
                            if v2 != None:
                                self.vote_tree.add_edge(v1, v2)

                return v1
            else:
                return None



    def match_node_centricity(self, c1, c2, threhold=.5):
        """
        match node centricity of the graph
        """
        node_cent1 = self.skel1.node_centricity[c1]
        node_cent2 = self.skel2.node_centricity[c2]
        match_result = abs(node_cent1 - node_cent2) / (node_cent1 + node_cent2)
        threhold *= 0.5
        return match_result < threhold


    def match_length_radius(self, n1, n2, matched_pairs, threhold=.5):
        """
        match path length and radius from n1/n2 to nodes that already in matched_pairs
        """
        path_len1 = self.skel1.path_length_ratio[n1, matched_pairs[:,0]]
        path_len2 = self.skel2.path_length_ratio[n2, matched_pairs[:,1]]
        length_match = abs(path_len1 - path_len2) / (path_len1 + path_len2)
        threhold *= 0.5
        #if all satisfied
        if np.all(length_match < threhold):
            path_rad1 = self.skel1.path_radius_ratio[n1, matched_pairs[:,0]]
            path_rad2 = self.skel2.path_radius_ratio[n2, matched_pairs[:,1]]
            radius_match = 2 * abs(path_rad1 - path_rad2) / (path_rad1 + path_rad2)
            if np.all(radius_match < threhold):
                return True
            else:
                return False
        else:
            return False


    def match_topology_consistency(self, n1, n2, matched_pairs):
        """
        match skeleton topology consistency
        """
        if len(matched_pairs) > 0:
            junct1 = matched_pairs[matched_pairs[:,0] < len(self.skel1.junction_index), 0]
            junct2 = matched_pairs[matched_pairs[:,1] < len(self.skel2.junction_index), 1]
            if len(junct1) < 1 or len(junct2) < 1:
                return False
            else:
                idx1 = np.argmin(self.skel1.path_to_junction[n1, junct1])
                idx2 = np.argmin(self.skel2.path_to_junction[n2, junct2])
                return [junct1[idx1], junct2[idx2]] in matched_pairs.tolist()
        else:
            print 'none in matched_pairs'
            return False


    def match_spatial_configuration(self):
        """
        match spatial configuration
        """
        pass


    def elector_vote(self):
        """
        use elector vote to find better correspondence
        """
        vote_matrix = np.zeros((len(self.skel1.feature_node_index), len(self.skel2.feature_node_index)))
        for v in self.vote_tree.vertices():
            if v.out_degree() < 2:
                v_list, e_list = topology.shortest_path(self.vote_tree, self.vote_tree.vertex(0), v)
                if len(v_list) > 4:
                    for v in v_list[1:]:
                        pair = self.node_pair[v]
                        vote_matrix[pair[0], pair[1]] += 1
        
        self.vote_matrix = vote_matrix


if __name__ == '__main__':
    from skeleton_data import SkeletonData
    from display_skeleton import DrawSkeleton
    from mayavi import mlab
    skel_name1 = './data/chair_skeleton/1_ckel.cg'
    mesh_name1 = './data/chair/1.off'
    skel_name2 = './data/chair_skeleton/2_ckel.cg'
    mesh_name2 = './data/chair/2.off'
    sskel1 = SkeletonData(fname=skel_name1, mesh_name=mesh_name1, filter_sb=True)
    sskel2 = SkeletonData(fname=skel_name2, mesh_name=mesh_name2, filter_sb=True)
    skel_match = SkeletonMatch(skel1=sskel1, skel2=sskel2)
    print 'tree vertex num', skel_match.vote_tree.num_vertices()
    skel_match.elector_vote()

    mlab.figure(1)
    draw_skel1 = DrawSkeleton(sskel1)
    draw_skel2 = DrawSkeleton(sskel2)
    draw_skel1.draw_all(point_visible=True)
    draw_skel1.draw_feature_node()
    draw_skel2.draw_all(point_visible=True)
    draw_skel2.draw_feature_node()

    mlab.figure(2)
    mlab.imshow(skel_match.vote_matrix)
    print skel_match.vote_matrix
    
    mlab.show()
                


               # for n, pair in enumerate(curr_pairs):
               #     new_prev = np.vstack((prev_pairs, pair))
               #     new_curr = np.delete(curr_pairs, n, 0)
               #     v2 = self._construct_search_tree(prev_pairs=new_prev, curr_pairs=new_curr)
               #     if v2 != None:
               #         self.vote_tree.add_edge(v1, v2)
