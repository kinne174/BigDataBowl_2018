import pandas as pd
import os
from collections import namedtuple
from string_to_LineString import s2LS
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from shapely.geometry import LineString
from load_data import Data
from route_tree import rotate_route

gameId = '2017091101'

header = r'C:\Users\Mitch\Documents\UofM\Fall 2018\NFL\Data'
save_header = r'C:\Users\Mitch\Documents\UofM\Fall 2018\NFL\Manuscript\figs'
routes_df = pd.read_csv(os.path.join(header, r'Routes_{}.csv'.format(gameId)))
play_df = routes_df[['LineString', 'position', 'route']]

#organize routes
ROUTEDOC = namedtuple('RouteDoc', 'LineString position route')
routes_info = []
for _, row in play_df.iterrows():
    routes_info.append(ROUTEDOC(LineString=s2LS(row['LineString']), position=row['position'], route=row['route']))

route_names = ['out',]# 'post', 'slant',]# 'flat', 'wheel', 'streak', 'corner', 'comeback', 'curl', 'bubble.block', 'out']
route_colors = ['g',]# 'r', 'b']

route_linestrings = [[doc.LineString for doc in routes_info if doc.route == r] for r in route_names]

fig = plt.figure()
# fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
ax = fig.add_subplot(111, xlim=(-25, 25), ylim=(-5, 35))
ax.set_aspect(aspect='equal', adjustable='box')
# plt.title(' and '.join(['{} in {}'.format(r, c) for r, c in zip(route_names, route_colors)]))

for ind, ls_list in enumerate(route_linestrings):
    for ls in ls_list:
        assert isinstance(ls, LineString)

        x_pre, _ = ls.xy
        running_right = x_pre[0] < x_pre[6]

        coords = rotate_route(ls, running_right)

        assert isinstance(coords, LineString)

        x, y = coords.xy

        ax.plot(x, y, route_colors[ind], alpha=0.3)

rect1 = patches.Rectangle((-25, 0), 50, 5, linewidth=1, edgecolor='k', facecolor='none')
rect2 = patches.Rectangle((-25, 0), 50, 10, linewidth=1, edgecolor='k', facecolor='none')
rect3 = patches.Rectangle((-25, 0), 50, 15, linewidth=1, edgecolor='k', facecolor='none')
rect4 = patches.Rectangle((-25, 0), 50, 20, linewidth=1, edgecolor='k', facecolor='none')
rect5 = patches.Rectangle((-25, 0), 50, 25, linewidth=1, edgecolor='k', facecolor='none')

ax.add_patch(rect1)
ax.add_patch(rect2)
ax.add_patch(rect3)
ax.add_patch(rect4)
ax.add_patch(rect5)

plt.axis('off')
fig.show()
plt.close(fig)
# fig.savefig(fname=os.path.join(save_header, r'game_outs.png'))
