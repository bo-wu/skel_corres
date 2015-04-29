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
            centricity_matched_pairs = []
            for i, j in itertools.product(skel1_index, skel2_index):
                if self.match_node_centricity(c1=i, c2=j, threhold=.5):
                    centricity_matched_pairs.append([i,j])
            self.centricity_matched_pairs = np.array(centricity_matched_pairs)
        else:
            print 'need input two skeleton to match'


    def match_skeleton(self):
        """
        match all related metrics
        """
        vote_tree = Graph(directed=False)
        #plan to store all pairs from root util here
        node_pair = vote_tree.new_vertex_property("vector<short>")
        ## root for all the other nodes ##
        v0 = vote_tree.add_vertex()
        vote_tree.add_vertex(len(self.centricity_matched_pairs))
        v_num = vote_tree.num_vertices()
        #from the second node (first is root)
        for i in xrange(1, v_num):     
            graph_v = vote_tree.vertex(i)
            node_pair[graph_v] = self.centricity_matched_pairs[i-1]
            vote_tree.add_edge(v0, graph_v)

        for 


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
    skel_match.match_skeleton()
    mlab.figure(1)
    draw_skel1 = DrawSkeleton(sskel1)
    draw_skel2 = DrawSkeleton(sskel2)
    draw_skel1.draw_all(point_visible=True)
    draw_skel1.draw_feature_node()
    draw_skel2.draw_all(point_visible=True)
    draw_skel2.draw_feature_node()
    mlab.show()
