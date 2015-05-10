#! /usr/bin/env python
from mayavi import mlab
import numpy as np
from skeleton_data import SkeletonData

class DrawSkeleton(object):
    def __init__(self, skel_data):
        self.skel_data = skel_data
        self.program_id = np.random.randint(100, size=1)

    def draw_all(self, point_visible=True, color=(.8,.8,.0)):
        #mlab.figure(self.program_id[0])
        #mlab.clf()
        pts = mlab.points3d(self.skel_data.verts[:,0], self.skel_data.verts[:,1], self.skel_data.verts[:,2], color=color, scale_factor=.020, resolution=20)
        pts.mlab_source.dataset.lines = self.skel_data.edges
        pts.visible = point_visible

       ## lines = mlab.pipeline.stripper(pts)
       ## mlab.pipeline.surface(lines, color=(.8,.8,.8), line_width=8, )

        tube = mlab.pipeline.tube(pts, tube_radius=.01, tube_sides=20)
        mlab.pipeline.surface(tube, color=color)

    def draw_terminal(self):
        #mlab.figure(self.program_id[0])
        pts = mlab.points3d(self.skel_data.terminal[:,0], self.skel_data.terminal[:,1], self.skel_data.terminal[:,2], color=(.0,.0,.7), scale_factor=.03, resolution=20)

    def draw_junction(self):
        #mlab.figure(self.program_id[0])
        pts = mlab.points3d(self.skel_data.junction[:,0], self.skel_data.junction[:,1], self.skel_data.junction[:,2], color=(.7,0.,0.), scale_factor=.035, resolution=20)

    def draw_feature_node(self):
        label_text = []
        junc_len = len(self.skel_data.junction_index)
        for i in xrange(len(self.skel_data.feature_node_index)):
            label_text.append(str(i))
        #mlab.figure(self.program_id[0])
        pts = mlab.points3d(self.skel_data.feature_node[:junc_len,0], self.skel_data.feature_node[:junc_len,1], self.skel_data.feature_node[:junc_len,2], color=(.7,.0,0.0), scale_factor=.135, resolution=20)
        mlab.points3d(self.skel_data.feature_node[junc_len:,0], self.skel_data.feature_node[junc_len:,1], self.skel_data.feature_node[junc_len:,2], color=(.0,.0,0.7), scale_factor=.1, resolution=20)
        for i in xrange(len(label_text)):
            mlab.text3d(self.skel_data.feature_node[i,0], self.skel_data.feature_node[i,1], self.skel_data.feature_node[i,2], label_text[i], scale=.15)

if __name__ == '__main__':
    skeleton = SkeletonData('data/psb_skeleton/115_ckel.cg', True)
    #skeleton = SkeletonData('data/chair_skeleton/115_ckel.cg', True)
    draw_skel = DrawSkeleton(skeleton)
    draw_skel.draw_all(point_visible=True)
    draw_skel.draw_feature_node()
    #draw_skel.draw_terminal()
    #draw_skel.draw_junction()
    mlab.show()
