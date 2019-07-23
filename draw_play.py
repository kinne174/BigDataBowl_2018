import matplotlib.pyplot as plt
import matplotlib.animation as animation
from load_data import Data
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--gameId', type=int, help = 'gameId of interest')
parser.add_argument('--playId', type=int, help = 'playId of interest')
parser.add_argument('--play', type=int, help='play in dictionary')
args = parser.parse_args()


#get the play from all tracking game data
def subset_by_playId(tracking_df, playId):
    return tracking_df.loc[tracking_df['playId'] == int(playId)]

#update objects by updating frame id
def _animate(i):
    global field, ax, fig, tracking_df

    i += 1
    Offense_positions.set_data(list(tracking_df['x'].loc[(tracking_df['side'] == 'offense') & (tracking_df['frame.id'] == i) & (tracking_df['displayName'] != 'football')]), list(tracking_df['y'].loc[(tracking_df['side'] == 'offense') & (tracking_df['frame.id'] == i) & (tracking_df['displayName'] != 'football')]))
    Defense_positions.set_data(list(tracking_df['x'].loc[(tracking_df['side'] == 'defense') & (tracking_df['frame.id'] == i) & (tracking_df['displayName'] != 'football')]), list(tracking_df['y'].loc[(tracking_df['side'] == 'defense') & (tracking_df['frame.id'] == i) & (tracking_df['displayName'] != 'football')]))
    Football_position.set_data(list(tracking_df['x'].loc[(tracking_df['frame.id'] == i) & (tracking_df['displayName'] == 'football')]), list(tracking_df['y'].loc[(tracking_df['frame.id'] == i) & (tracking_df['displayName'] == 'football')]))
    field.set_edgecolor('k')
    hash_lines.set_edgecolor('k')
    goal_lines.set_edgecolor('k')

    for pa, ind_row in zip(Position_annotations, tracking_df[(tracking_df['frame.id'] == i) & (tracking_df['displayName'] != 'football')].iterrows()):
        row = ind_row[1]
        pa.set_position((row['x'], row['y']))

    return (Offense_positions, Defense_positions, Football_position, field, hash_lines, goal_lines, *Position_annotations)

#initialize objects inside animation
def _animate_init():
    global field, tracking_df

    Offense_positions.set_data([], [])
    Defense_positions.set_data([], [])
    Football_position.set_data([], [])
    field.set_edgecolor('none')
    hash_lines.set_edgecolor('none')
    goal_lines.set_edgecolor('none')
    for pa, ind_row in zip(Position_annotations, tracking_df[(tracking_df['frame.id'] == 1) & (tracking_df['displayName'] != 'football')].iterrows()):
        row = ind_row[1]
        pa.set_position((row['x'], row['y']))

    return (Offense_positions, Defense_positions, Football_position, field, hash_lines, goal_lines, *Position_annotations)

D = Data()

#define play and game
play_list = [(2017091101, 565), (2017091710, 303), (2017092406, 1802), (2017091004, 2293)]

problem_tuples = [(2017100108, 4098, 4),]

gameId = str(2017092406)
playId = str(1802)

# USE_CMDLINE = False
# gameId = str(args.gameId) if USE_CMDLINE else '2017092406'
# playId = (args.playId) if USE_CMDLINE else '1802'

D.load_tracking_info(gameId)
tracking_df_dic = D.get_tracking_info()
players_df = D.get_all_players_info()
tracking_df = tracking_df_dic[gameId]

tracking_df = subset_by_playId(tracking_df, playId)

if tracking_df.empty:
    raise Exception('The playId does not exist in that game!')

fig = plt.figure()
fig.subplots_adjust(left=0, right=1,bottom=0,top=1)
ax = fig.add_subplot(111, xlim=(-10, 130), ylim=(-5, 58.3))
ax.set_aspect(aspect='equal', adjustable='box')
ax.set_title('gameId: {}, playId: {}'.format(str(gameId), str(playId)))

Offense_positions, = ax.plot([], [], 'ro', ms=7)
Defense_positions, = ax.plot([], [], 'go', ms=7)
Football_position, = ax.plot([], [], 'mo', ms=5)

Position_annotations = [ax.annotate(row['position'], xy=(row['x'], row['y'])) for _, row in tracking_df[(tracking_df['frame.id'] == 1) & (tracking_df['displayName'] != 'football')].iterrows()]
for pa in Position_annotations:
    pa.set_animated(True)

#draw field
field = plt.Rectangle(xy=(0, 0), width=120, height=53.3, fill=False, lw=2)
hash_lines = plt.Rectangle(xy=(10, 23.36667), width = 100, height=6.6, fill=False, lw=1, linestyle='--')
goal_lines = plt.Rectangle(xy=(10, 0), width=100, height=53.3, fill=False, lw=1)
ax.add_patch(field)
ax.add_patch(hash_lines)
ax.add_patch(goal_lines)

animation = animation.FuncAnimation(fig, _animate, frames=max(tracking_df['frame.id']), interval=50, init_func=_animate_init, blit=True, repeat=True, repeat_delay=50)

plt.show()