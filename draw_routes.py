import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import argparse
import os
from collections import namedtuple
from string_to_LineString import s2LS
from matplotlib.patches import Patch

parser = argparse.ArgumentParser()
parser.add_argument('--gameId', type=int, help = 'gameId of interest')
parser.add_argument('--playId', type=int, help = 'playId of interest')
parser.add_argument('--play', type=int, help='play in dictionary')
args = parser.parse_args()

#get route information
play_list = [(2017091101, 565), (2017091710, 303)]

gameId = play_list[args.play][0]
playId = play_list[args.play][1]

# USE_CMDLINE = False
# gameId = str(args.gameId) if USE_CMDLINE else '2017091004'
# playId = (args.playId) if USE_CMDLINE else '2293'

header = r'C:\Users\Mitch\Documents\UofM\Fall 2018\NFL\Data'
routes_df = pd.read_csv(os.path.join(header, r'Routes_{}.csv'.format(gameId)))
play_df = routes_df.loc[routes_df['playId'] == int(playId), ['LineString', 'position', 'route']]

#organize routes
ROUTEDOC = namedtuple('RouteDoc', 'LineString position route')
routes_info = []
for _, row in play_df.iterrows():
    routes_info.append(ROUTEDOC(LineString=s2LS(row['LineString']), position=row['position'], route=row['route']))

def _animate_init():
    field.set_edgecolor('none')
    hash_lines.set_edgecolor('none')
    goal_lines.set_edgecolor('none')

    for ra, rd in zip(Route_annotations, routes_info):
        x, y = rd.LineString.coords.xy
        x = x[0]
        y = y[0]
        ra.set_position((x, y))

    for Rd, rd in zip(Route_drawings, routes_info):
        Rd[0].set_data([], [])

    return (field, hash_lines, goal_lines, *Route_drawings, *Route_annotations)

def _animate(i):
    field.set_edgecolor('k')
    hash_lines.set_edgecolor('k')
    goal_lines.set_edgecolor('k')

    for ra, rd in zip(Route_annotations, routes_info):
        x, y = rd.LineString.coords.xy
        x = x[0]
        y = y[0]
        ra.set_position((x, y))
        ra.xy = (x, y)

    for Rd, rd in zip(Route_drawings, routes_info):
        x, y = rd.LineString.coords.xy
        x = x[:i]
        y = y[:i]
        Rd[0].set_data(x, y)
        Rd[0].axes.axis([-10, 130, -5, 58.3])

    return (field, hash_lines, goal_lines, *Route_drawings, *Route_annotations)

fig = plt.figure()
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
ax = fig.add_subplot(111, xlim=(-10, 130), ylim=(-5, 58.3))
ax.set_aspect(aspect='equal', adjustable='box')
#ax.set_title('gameId: {}, playId: {}'.format(str(gameId), str(playId)))

Route_annotations = [ax.annotate(rd.route, xy=list(rd.LineString.coords)[0], xytext=list(rd.LineString.coords)[0], size='large') for rd in routes_info]

color_dic = {"WR": 'k', "TE": 'orange', "RB": 'g'}

Route_drawings = [ax.plot([], [], color_dic[rd.position], lw=3) for rd in routes_info]

legend_elements = [Patch(facecolor='k', edgecolor='k', label='WR'),
                   Patch(facecolor='orange', edgecolor='orange', label='TE'),
                   Patch(facecolor='g', edgecolor='g', label='RB')]
ax.legend(handles=legend_elements)

#draw field
field = plt.Rectangle(xy=(0, 0), width=120, height=53.3, fill=False, lw=2)
hash_lines = plt.Rectangle(xy=(10, 23.36667), width = 100, height=6.6, fill=False, lw=1, linestyle='--')
goal_lines = plt.Rectangle(xy=(10, 0), width=100, height=53.3, fill=False, lw=1)
ax.add_patch(field)
ax.add_patch(hash_lines)
ax.add_patch(goal_lines)

animation = animation.FuncAnimation(fig, _animate, frames=len(routes_info[0].LineString.coords), interval=50, init_func=_animate_init, blit=False, repeat=True, repeat_delay=500*2)

plt.show()






