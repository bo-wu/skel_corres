#! /usr/bin/env python
from graph_tool import Graph, util
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
            centricity_matched_pairs = []
            junc_pairs = []
            for i, j in itertools.product(skel1_index, skel2_index):
                if self.match_node_centricity(c1=i, c2=j, threhold=.5):
                    centricity_matched_pairs.append([i,j])
            self.centricity_matched_pairs = np.array(centricity_matched_pairs)

            self.vote_tree = Graph(directed=False)
            self.node_pair = self.vote_tree.new_vertex_property("vector<short>")

            self._construct_search_tree(curr_pairs=self.centricity_matched_pairs)
        else:
            print 'need input two skeleton to match'


    def _construct_search_tree(self, prev_pairs=np.array([]), curr_pairs=np.array([])):
        """
        recursively consturct search tree
        @param prev_pairs record that already on the path
        @param curr_pairs record that left pairs (which pass first test)
        """
        # root of the tree
        if len(prev_pairs) == 0:
            v1 = self.vote_tree.add_vertex()
            for n, pair in enumerate(curr_pairs):
                new_prev = pair.reshape(-1,2)  # to use len(for level one), need to change shape
                new_curr = np.delete(curr_pairs, n, 0)
                v2 = self._construct_search_tree(prev_pairs=new_prev, curr_pairs=new_curr)
                if v2 != None:
                    self.vote_tree.add_edge(v1, v2)
            return v1

        elif len(prev_pairs) == 1:  # first level
            v1 = self.vote_tree.add_vertex()
            self.node_pair[v1] = prev_pairs[-1,:]
            for n, pair in enumerate(curr_pairs):
                new_prev = np.vstack((prev_pairs, pair))
                new_curr = np.delete(curr_pairs, n, 0)
                v2 = self._construct_search_tree(prev_pairs=new_prev, curr_pairs=new_curr)
                if v2 != None:
                    self.vote_tree.add_edge(v1, v2)
            return v1

        elif len(prev_pairs) > 1:  # above level two
            if self.match_length_radius(n1=prev_pairs[-1,0], n2=prev_pairs[-1,1], matched_pairs=prev_pairs[:-1]):
                v1 = self.vote_tree.add_vertex()
                self.node_pair[v1] = prev_pairs[-1,:]
                for n, pair in enumerate(curr_pairs):
                    new_prev = np.vstack((prev_pairs, pair))
                    new_curr = np.delete(curr_pairs, n, 0)
                    v2 = self._construct_search_tree(prev_pairs=new_prev, curr_pairs=new_curr)
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
        pass


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
    print 'centricity matched pairs', skel_match.centricity_matched_pairs
    print 'tree vertex num', skel_match.vote_tree.num_vertices()
   # skel_match.match_skeleton()
   # mlab.figure(1)
   # draw_skel1 = DrawSkeleton(sskel1)
   # draw_skel2 = DrawSkeleton(sskel2)
   # draw_skel1.draw_all(point_visible=True)
   # draw_skel1.draw_feature_node()
   # draw_skel2.draw_all(point_visible=True)
   # draw_skel2.draw_feature_node()
   # mlab.show()
