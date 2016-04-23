from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np

map = Basemap(projection='robin', resolution = 'l', area_thresh = 1000.0,
              lat_0=0, lon_0=-130)
map.drawcoastlines()
map.drawcountries()
map.fillcontinents(color = 'gray')
map.drawmapboundary()
map.drawmeridians(np.arange(0, 360, 30))
map.drawparallels(np.arange(-90, 90, 30))

def get_marker_color(magnitude):
    if magnitude < 3.0:
        return ('go')
    elif magnitude < 5.0:
        return ('yo')
    else:
        return ('ro')

# Variable size dots:
#  go through each lat, lon, plot it individually, calculating size dynamically
#  magnitude 1.0 = min dot size; larger quakes scaled by magnitude
min_marker_size = 4
# for lon, lat, mag in zip(lons, lats, magnitudes):
#
#     msize = mag * min_marker_size
#     marker_string = get_marker_color(mag)
#     map.plot(x, y, marker_string, markersize=msize)
x,y = map(-122.74, 38.78)
map.plot(x, y, 'ro', markersize=100)
title_string = "Earthquakes of Magnitude 1.0 or Greater\n"
#title_string += "%s through %s" % (timestrings[-1], timestrings[0])
plt.title(title_string)

plt.show()





# from mpl_toolkits.basemap import Basemap
# import matplotlib.pyplot as plt
# import numpy as np
# import matplotlib.animation as animation
# map = Basemap(projection='merc',
#               ellps='WGS84',
#               llcrnrlon=-160,urcrnrlon=-60,llcrnrlat=10,urcrnrlat=80,
#               resolution="l")
# map.drawcoastlines()
# map.drawcountries()
# #map.fillcontinents(color = 'gray')
# map.drawmapboundary()
# map.drawlsmask(land_color='coral',ocean_color='aqua',lakes=True)
# map.drawmeridians(np.arange(0, 360, 30))
# map.drawparallels(np.arange(-90, 90, 30))
#
# x,y = map(0, 0)
# point = map.plot(x, y, 'ro', markersize=5)[0]
#
# def init():
#     point.set_data([], [])
#     return point,
#
# # animation function.  This is called sequentially
# def animate(i):
#     lons, lats =  np.random.random_integers(-130, 130, 2)
#     x, y = map(lons, lats)
#     point.set_data(x, y)
#     return point,
#
# # call the animator.  blit=True means only re-draw the parts that have changed.
# anim = animation.FuncAnimation(plt.gcf(), animate, init_func=init,
#                                frames=20, interval=500, blit=False)
#

"""
Draw a graph with matplotlib.
You must have matplotlib for this to work.
"""


# try:
#     import matplotlib.pyplot as plt
#     import matplotlib.colors as colors
#     import matplotlib.cm as cmx
#     import numpy as np
# except:
#     raise
#
# import networkx as nx
#
# G=nx.path_graph(8)
# #Number of edges is 7
# values = range(7)
# # These values could be seen as dummy edge weights
#
# jet = cm = plt.get_cmap('jet')
# cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
# scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
# colorList = []
#
# for i in range(7):
#   colorVal = scalarMap.to_rgba(values[i])
#   colorList.append(colorVal)
#
#
# nx.draw(G,edge_color=colorList)
# plt.savefig("simple_path.png") # save as png
# plt.show() # display
#

# cut = 1.05
# xmax = cut * max(xx for xx,yy in pos.values())
# xmin =  min(xx for xx,yy in pos.values())
# xmin = xmin - (cut * xmin)
#
#
# ymax = cut * max(yy for xx,yy in pos.values())
# ymin = (cut) * min(yy for xx,yy in pos.values())
# ymin = ymin - (cut * ymin)
#
# plt.xlim(-11563.71285631819,10402200.831035526)
# plt.ylim(-68403.0940226051,4684384.722747428)
#
# plt.show()


    # pos = dict()
    #
    # for pos_node in network.nodes():
    #     # Normalize the lat and lon values
    #     x,y = m(float(network.node[pos_node]['lon']),
    #             float(network.node[pos_node]['lat']))
    #
    #     pos[pos_node] = [x,y]
    #
    # #m.drawmapboundary("aqua")
    # #m.fillcontinents('#555555')
    # m.drawlsmask(land_color='coral',ocean_color='aqua',lakes=True)
    # #m.bluemarble()
    #
    # # First pass - Green lines
    # nx.draw_networkx_edges(network,pos,edgelist=network.edges(),
    #         width=1,
    #         edge_color="green",
    #         alpha=0.5,
    #         arrows=False)
    #
    # nx.draw_networkx_nodes(network,
    #         pos,
    #         linewidths=1,
    #         node_size=10,
    #         with_labels=False,
    #         node_color = "green")
    #
    # #m.bluemarble()
    # #plt.title=title
    #
    # # Adjust the plot limits
    # cut = 1.05
    # xmax = cut * max(xx for xx,yy in pos.values())
    # xmin =  min(xx for xx,yy in pos.values())
    # xmin = xmin - (cut * xmin)
    #
    #
    # ymax = cut * max(yy for xx,yy in pos.values())
    # ymin = (cut) * min(yy for xx,yy in pos.values())
    # ymin = ymin - (cut * ymin)
    #
    # plt.xlim(xmin,xmax)
    # plt.ylim(ymin,ymax)
    #
    # plt.axis('off')
    # plt.show()
    # plt.close()
    #
    #
