from tree import *
from random import *

#======================================================================
#            CREATE VTK OUTPUT OF TREE FOR USE WITH PARAVIEW
#                        by C.P. Dullemond
#                              2014
#======================================================================


#
# Count the nr of nodes within a node
# Note that countnodes(0) must (will) be equal to len(q.nodelist)
#
def countnodes(q,nodeid):
    nr = 1
    if q.nodelist[nodeid].leaf == 0:
        for iz in range(0,2):
            for iy in range(0,2):
                for ix in range(0,2):
                    id = q.nodelist[nodeid].childrenids[ix,iy,iz]
                    nr = nr + countnodes(q,id)
    return nr


#
# Count the nr of leafs within a node
#
def countleafs(q,nodeid):
    nr = q.nodelist[nodeid].leaf
    if q.nodelist[nodeid].leaf == 0:
        for iz in range(0,2):
            for iy in range(0,2):
                for ix in range(0,2):
                    id = q.nodelist[nodeid].childrenids[ix,iy,iz]
                    nr = nr + countleafs(q,id)
    return nr


#
# Write node corners
#
def vtkwritecorners(q,f,nodeid,alsobranches=0):
    if q.nodelist[nodeid].leaf == 1 or alsobranches == 1:
        xc = [0.,0.,0.]
        xc[:] = q.nodelist[nodeid].xc
        size  = q.nodelist[nodeid].size
        for iz in range(0,2):
            for iy in range(0,2):
                for ix in range(0,2):
                    f.write(str(xc[0]+0.5*size*(2*ix-1))+' '+str(xc[1]+0.5*size*(2*iy-1))+' '+str(xc[2]+0.5*size*(2*iz-1))+'\n')
    if q.nodelist[nodeid].leaf == 0:
        for iz in range(0,2):
            for iy in range(0,2):
                for ix in range(0,2):
                    id = q.nodelist[nodeid].childrenids[ix,iy,iz]
                    vtkwritecorners(q,f,id,alsobranches)

#
# Write node cells
#
def vtkwritecells(q,f,nodeid,ipt,alsobranches=0):
    if q.nodelist[nodeid].leaf == 1 or alsobranches == 1:
        f.write('8 ')
        for i in range(8):
            f.write(str(ipt)+' ')
            ipt = ipt + 1
        f.write('\n')
    if q.nodelist[nodeid].leaf == 0:
        for iz in range(0,2):
            for iy in range(0,2):
                for ix in range(0,2):
                    id = q.nodelist[nodeid].childrenids[ix,iy,iz]
                    ipt = vtkwritecells(q,f,id,ipt,alsobranches)
    return ipt

#
# Write node masses
#
def vtkwritemass(q,f,nodeid,alsobranches=0):
    if q.nodelist[nodeid].leaf == 1 or alsobranches == 1:
        f.write(str(q.nodelist[nodeid].mass)+'\n')
    if q.nodelist[nodeid].leaf == 0:
        for iz in range(0,2):
            for iy in range(0,2):
                for ix in range(0,2):
                    id = q.nodelist[nodeid].childrenids[ix,iy,iz]
                    vtkwritemass(q,f,id,alsobranches)

#
# Write a VTK dataset for visualization of (a part of) the tree
# You can use such a file with e.g. PARAVIEW to visualize the tree
# in 3D.
#
def vtk(q,nodeid,alsobranches=0):
    #
    # Count the number of nodes
    #
    if alsobranches == 1:
        nr = countnodes(q,nodeid)
    else:
        nr = countleafs(q,nodeid)
    #
    # Open file and write header
    #
    f = open('tree.vtk','w')
    f.write('# vtk DataFile Version 1.0\n')
    f.write('Tree (Python)\n')
    f.write('ASCII\n\n')
    f.write('DATASET UNSTRUCTURED_GRID\n')
    f.write('POINTS ')
    f.write(str(8*nr))
    f.write('  float\n')
    vtkwritecorners(q,f,nodeid,alsobranches)
    f.write('\n')
    f.write('CELLS ')
    f.write(str(nr)+' '+str(9*nr)+'\n')
    vtkwritecells(q,f,nodeid,0,alsobranches)
    f.write('\n')
    f.write('CELL_TYPES '+str(nr)+'\n')
    for i in range(nr):
        f.write('11\n')
    f.write('\n')
    f.write('CELL_DATA '+str(nr)+'\n')
    f.write('SCALARS mass float\n')
    f.write('LOOKUP_TABLE default\n')
    vtkwritemass(q,f,nodeid,alsobranches)
    f.write('\n')
    f.close()


nparticles = 1000
particles = []
for i in range(nparticles):
    x = random()
    y = random()
    z = random()
    particles.append([x,y,z,1.0])

q=TreeClass(particles)
q.insertallparticles()
q.computemultipoles(0)
vtk(q,0,0)


# Now use PARAVIEW to load the file tree.vtk and view in 3D
