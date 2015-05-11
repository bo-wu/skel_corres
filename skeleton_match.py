#! /usr/bin/env python
from graph_tool import Graph, util, topology
import numpy as np
import itertools


class SkeletonMatch(object):
    """
    implement Oscar's skeleton matching alogrithm in this class
    """
    def __init__(self, skel1, skel2):
        if skel1 is not None and skel2 is not None :
            self.skel1 = skel1
            self.skel2 = skel2
            self.skel1.calc_skel_properties()
            self.skel2.calc_skel_properties()
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
            #self.all_junc_pairs = np.vstack((self.junction_pairs, self.junc_term_pairs))

            self.vote_tree = Graph(directed=False)
            self.node_pair = self.vote_tree.new_vertex_property("vector<short>")

            self._construct_voting_tree()
        else:
            print 'need input two skeleton to match'


    def _construct_voting_tree(self, prev_pairs=np.array([])):
        """
        recursively consturct voting tree
        @param prev_pairs record that already on the path
        @param junc_pairs record that left part junction pairs (on current tree level)
        @param term_pairs record that left terminal pairs on current tree level

        now limits: at least one junction pair
        """
        # root of the tree
        if len(prev_pairs) == 0:
            v1 = self.vote_tree.add_vertex()
            #first level, only junction pairs (both are junction)
            for n, pair in enumerate(self.junction_pairs):
                new_prev = pair.reshape(-1,2)  # to use len(for level one), need to change shape
                print 'adding subtree', n+1, '/', len(self.junction_pairs)
                v2 = self._construct_voting_tree(prev_pairs=new_prev)
                self.vote_tree.add_edge(v1, v2)
            return v1                    # return the root

        elif len(prev_pairs) == 1:  # first level
            v1 = self.vote_tree.add_vertex()
            self.node_pair[v1] = prev_pairs.flatten()
            """
            priority order: junction pairs, termianl pairs, junc-term pairs
            """

            check_junc = True
            #curr_junc = self.all_junc_pairs.copy()
            #counter = 0
            #prepare for next(second) level
            for n, pair in enumerate(self.junction_pairs):
                if pair[0] not in prev_pairs[:,0] and pair[1] not in prev_pairs[:,1]:
                    new_prev = np.vstack((prev_pairs, pair))
                    check_junc = False
                    v2 = self._construct_voting_tree(prev_pairs=new_prev)
                    if v2 is not None: 
                        self.vote_tree.add_edge(v1, v2)
           #     else:
           #         curr_junc = np.delete(curr_junc, n-counter, 0)
           #         counter += 1

            # it is sure that that should be some terminal_pairs
            # but in case
            check_term = False     # if allow mix junc and term 
            if check_junc:
                for n, pair in enumerate(self.terminal_pairs):
                    new_prev = np.vstack((prev_pairs, pair))
                    check_term = False
                    v2 = self._construct_voting_tree(prev_pairs=new_prev)
                    if v2 is not None:
                        self.vote_tree.add_edge(v1, v2)

            if check_junc and check_term:
                for n, pair in enumerate(self.junc_term_pairs):
                    new_prev = np.vstack((prev_pairs, pair))
                    v2 = self._construct_voting_tree(prev_pairs=new_prev)
                    if v2 is not None:
                        self.vote_tree.add_edge(v1, v2)

            return v1                   # return the first level of the tree

        elif 4 > len(prev_pairs) > 1:  # above level two
            #if satisfy T2 (length and radius prune)
            if self.match_length_radius(n1=prev_pairs[-1,0], n2=prev_pairs[-1,1], matched_pairs=prev_pairs[:-1]) and self.match_topology_consistency(n1=prev_pairs[-1,0], n2=prev_pairs[-1,1], matched_pairs=prev_pairs[:-1]):
                v1 = self.vote_tree.add_vertex()
                self.node_pair[v1] = prev_pairs.flatten()

                check_junc = True
                for pair in self.junction_pairs:
                    if pair[0] not in prev_pairs[:,0] and pair[1] not in prev_pairs[:,1]:
                        new_prev = np.vstack((prev_pairs, pair))
                        check_junc = False
                        v2 = self._construct_voting_tree(prev_pairs=new_prev)
                        if v2 is not None:
                            self.vote_tree.add_edge(v1, v2)
                
                check_term = False    # if allow mix junc and term
                if check_junc:
                    for pair in self.terminal_pairs:
                        if pair[0] not in prev_pairs[:,0] and pair[1] not in prev_pairs[:,1]:
                            new_prev = np.vstack((prev_pairs, pair))
                            check_term = False
                            v2 = self._construct_voting_tree(prev_pairs=new_prev)
                            if v2 is not None:
                                self.vote_tree.add_edge(v1, v2)

                if check_junc and check_term:
                    for pair in self.junc_term_pairs:
                        if pair[0] not in prev_pairs[:,0] and pair[1] not in prev_pairs[:,1]:
                            new_prev = np.vstack((prev_pairs, pair))
                            v2 = self._construct_voting_tree(prev_pairs=new_prev)
                            if v2 is not None:
                                self.vote_tree.add_edge(v1, v2)

                return v1             # return second and above level tree
            else:
                return None             # fail to match

        elif len(prev_pairs) >= 4:
            if self.match_length_radius(n1=prev_pairs[-1,0], n2=prev_pairs[-1,1], matched_pairs=prev_pairs[:-1]) and self.match_topology_consistency(n1=prev_pairs[-1,0], n2=prev_pairs[-1,1], matched_pairs=prev_pairs[:-1]):
                print 'len(prev_pairs) >= 4',
                if self.match_spatial_configuration(n1=prev_pairs[-1,0], n2=prev_pairs[-1,1], matched_pairs=prev_pairs[:-1]):
                    v1 = self.vote_tree.add_vertex()
                    self.node_pair[v1] = prev_pairs.flatten()

                    check_junc = True
                    for pair in self.junction_pairs:
                        if pair[0] not in prev_pairs[:,0] and pair[1] not in prev_pairs[:,1]:
                            new_prev = np.vstack((prev_pairs, pair))
                            check_junc = False
                            v2 = self._construct_voting_tree(prev_pairs=new_prev)
                            if v2 is not None:
                                self.vote_tree.add_edge(v1, v2)

                    check_term = False   # if allow mix junction and terminal
                    if check_junc:
                        for pair in self.terminal_pairs:
                            if pair[0] not in prev_pairs[:,0] and pair[1] not in prev_pairs[:,1]:
                                new_prev = np.vstack((prev_pairs, pair))
                                check_term = False
                                v2 = self._construct_voting_tree(prev_pairs=new_prev)
                                if v2 is not None:
                                    self.vote_tree.add_edge(v1, v2)

                    if check_junc and check_term:
                        for pair in self.junc_term_pairs:
                            if pair[0] not in prev_pairs[:,0] and pair[1] not in prev_pairs[:,1]:
                                new_prev = np.vstack((prev_pairs, pair))
                                v2 = self._construct_voting_tree(prev_pairs=new_prev)
                                if v2 is not None:
                                    self.vote_tree.add_edge(v1, v2)
                    print 'succeed'
                    return v1
                else:
                    print 'failed'
                    return None
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
            radius_match = abs(path_rad1 - path_rad2) / (path_rad1 + path_rad2)
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
                print 'no junction node in already matched pairs'
                return False
            else:
                idx1 = np.argmin(self.skel1.path_to_junction[n1, junct1])
                idx2 = np.argmin(self.skel2.path_to_junction[n2, junct2])
                return [junct1[idx1], junct2[idx2]] in matched_pairs.tolist()
        else:
            print 'none in matched_pairs'
            return False


    def match_spatial_configuration(self, n1, n2, matched_pairs, threhold=3.0):
        """
        match spatial configuration
        """
        #need test if can be inverse
        skel1_vectors = self.skel1.normalized_verts[matched_pairs[-3:,0]] - self.skel1.normalized_verts[n1]
        skel2_vectors = self.skel2.normalized_verts[matched_pairs[-3:,1]] - self.skel2.normalized_verts[n2]
        a = np.dot(skel2_vectors, np.linalg.inv(skel1_vectors))
        u, s, v = np.linalg.svd(a)
        r = np.dot(u, v)
        if np.linalg.det(r) < 0:
            r *= -1.0
        res1 = np.linalg.norm(a-r)
        print 'res1', res1,
        if res1 > threhold:
            return False
        else:
            a = np.dot(skel1_vectors, np.linalg.inv(skel2_vectors))
            u, s, v = np.linalg.svd(a)
            r = np.dot(u, v)
            if np.linalg.det(r) < 0:
                r *= -1.0
            res2 = np.linalg.norm(a-r)
            print 'res2', res2,
            if res2 > threhold:
                return False

        return max(res1, res2) <= threhold


    def elector_vote(self):
        """
        use elector vote to find better correspondence
        """
        vote_matrix = np.zeros((len(self.skel1.feature_node_index), len(self.skel2.feature_node_index)))
        for v in self.vote_tree.vertices():
            if v.out_degree() < 2:
                pairs = self.node_pair[v]
                if len(pairs) >= 8:
                    temp_pairs = pairs.a.reshape(-1,2)
                    for pair in temp_pairs:
                        vote_matrix[pair[0], pair[1]] += 1

        self.vote_matrix = vote_matrix


if __name__ == '__main__':
    from skeleton_data import SkeletonData
    from display_skeleton import DrawSkeleton
    from mayavi import mlab
    #skel_pair = [11, 14]
    skel_pair = [1, 2]
    skel_name1 = './data/chair_skeleton/'+str(skel_pair[0])+'_ckel.cg'
    mesh_name1 = './data/chair/'+str(skel_pair[0])+'.off'
    skel_name2 = './data/chair_skeleton/'+str(skel_pair[1])+'_ckel.cg'
    mesh_name2 = './data/chair/'+str(skel_pair[1])+'.off'
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

    #print 'junction_pairs', skel_match.junction_pairs
    #print 'terminal_pairs', skel_match.terminal_pairs
    #print 'junc-term pairs', skel_match.junc_term_pairs

    
    mlab.show()
                


               # for n, pair in enumerate(curr_pairs):
               #     new_prev = np.vstack((prev_pairs, pair))
               #     new_curr = np.delete(curr_pairs, n, 0)
               #     v2 = self._construct_voting_tree(prev_pairs=new_prev, curr_pairs=new_curr)
               #     if v2 is not None:
               #         self.vote_tree.add_edge(v1, v2)
