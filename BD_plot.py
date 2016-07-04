#       # Save particle trajectory #0 for visualization
#       if args.plot == True and i == 0:
#           traj = numpy.append(traj, [[trajectory.x, trajectory.y, trajectory.z]], axis = 0)   




# # PLOT
# if args.plot:
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     ax.set_aspect('equal')
#     u = numpy.linspace(0, 2 * numpy.pi, 16)
#     v = numpy.linspace(0, numpy.pi, 16)
#     x = starting_radius * numpy.outer(numpy.cos(u), numpy.sin(v))
#     y = starting_radius * numpy.outer(numpy.sin(u), numpy.sin(v))
#     z = starting_radius * numpy.outer(numpy.ones(numpy.size(u)), numpy.cos(v))
#     Axes3D.plot_wireframe(ax, x,y,z)
#     x = r_max * numpy.outer(numpy.cos(u), numpy.sin(v))
#     y = r_max * numpy.outer(numpy.sin(u), numpy.sin(v))
#     z = r_max * numpy.outer(numpy.ones(numpy.size(u)), numpy.cos(v))
#     Axes3D.plot_wireframe(ax, x,y,z, colors = ['r'])
#     Axes3D.scatter(ax, start_pos[:,0], start_pos[:,1], start_pos[:,2], zdir='z', c='b')
#     Axes3D.scatter(ax, end_pos[:,0], end_pos[:,1], end_pos[:,2], zdir='z',  c='r')
#     Axes3D.plot(ax, traj[:,0], traj[:,1], traj[:,2], zdir='z', c='k')
#     plt.show()