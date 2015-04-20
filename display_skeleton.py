#! /usr/bin/env python
from mayavi import mlab
from skeleton_data import SkeletonData

class DrawSkeleton(object):
    def __init__(self, skel_data):
        self.skel_data = skel_data

    def draw_all(self, point_visible=False):
        mlab.figure(1)
        mlab.clf()
        pts = mlab.points3d(self.skel_data.verts[:,0], self.skel_data.verts[:,1], self.skel_data.verts[:,2], scale_factor=.025, resolution=20)
        pts.mlab_source.dataset.lines = self.skel_data.edges
        pts.visible = point_visible

       ## lines = mlab.pipeline.stripper(pts)
       ## mlab.pipeline.surface(lines, color=(.8,.8,.8), line_width=8, )

        tube = mlab.pipeline.tube(pts, tube_radius=.01, tube_sides=20)
        mlab.pipeline.surface(tube, color=(.8, .8, .0))

    def draw_terminal(self):
        mlab.figure(1)
        pts = mlab.points3d(self.skel_data.terminal[:,0], self.skel_data.terminal[:,1], self.skel_data.terminal[:,2], color=(.0,.0,.7), scale_factor=.03, resolution=20)

    def draw_junction(self):
        mlab.figure(1)
        pts = mlab.points3d(self.skel_data.junction[:,0], self.skel_data.junction[:,1], self.skel_data.junction[:,2], color=(.7,0.,0.), scale_factor=.035, resolution=20)

if __name__ == '__main__':
    skeleton = SkeletonData('data/psb_skeleton/115_ckel.cg')
    draw_skel = DrawSkeleton(skeleton)
    draw_skel.draw_all(point_visible=False)
    draw_skel.draw_terminal()
    draw_skel.draw_junction()
    mlab.show()
