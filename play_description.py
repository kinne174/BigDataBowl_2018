#build a dictionary of parameters to define what might be useful in play

#maybe first try ignoring pre snap defense stuff and just define a catchable ball by how much seperation there is and direction of defenders when ball is thrown
#then say which combination of routes works best given offense alignment (shotgun and wr positioning), time between snap and throw, deepness of the routes,

#might need to add into the model somehow how long it takes for the quarterback to get sacked depending on how many people rush/block

#can label routes success or failures depending on seperation at common release times
#
# to do for each passing play:
#   -label success based on catching the ball or 'drop', ow no success
#   -measure time between ball passed and ball caught/incomplete/intercepted. If sack do nothing
#   -measure time between ball snapped and ball passed
#   -when ball is passed the distance of the closest two defenders from the ball caught/incomplete/interception point, their direction relative to ball point and speed
#   -when ball is passed the distance of the intended wide receiver (closest one when ball is caught/incomplete/intercepted), direction and speed too
#   -distance ball caught/incomplete/intercepted is from closest sideline
#   -measure distance (pos and neg) of wide receivers from ball.
#   -depth of route up the field (interaction possibly with some time stuff?)
#   -blitz?
#   -offense personnel
#   -closest defender at pass

import numpy as np
from scipy.spatial import distance_matrix
import pandas as pd
from shapely.geometry import LineString
from string_to_LineString import s2LS
import matplotlib.pyplot as plt

#the play descriptors I wanted to exract and the type of data they are
play_parameters = {
    'blitz': bool, #whether 5+ players rush the quarterback
    'wr_bunch': bool, #whether three wr are together in a triangle
    'success': bool, #was the play a success complete pass or a drop
    'shotgun': bool, #was the quarterback in shotgun
    'time_between_snap_and_pass': int, #num frames between when the ball was snapped and when the ball was thrown
    'time_between_pass_and_outcome': int, #num frames between when the ball was passed and the outcome
    'closest_defenders_at_pass': list, #four element tuples of speed, distance, and direction relative to outcome point when ball is thrown, and position
    'closest_defenders_at_outcome': list, #four element tuples of speed, distance, and direction relative to outcome point at outcome, and position
    'intended_receiver_at_pass': tuple, #three element tuple of speed, distance, direction relative to outcome point when ball is thrown
    'all_receivers': list, # four element tuples of position, route, distance from ball at snap, depth
    'closest_sideline': float, #distance to the closest sideline from outcome point
    'num_closest_defenders_to_qb': float, #number of defenders close to quarterback at pass to measure pressure
    'in_pocket': bool, #whether the quarterback passed outside of the pocket
    'thrown_away': bool, #whether the quarterback threw the ball away
    'yards_after_catch': int,  # how many yards gained after the catch
    'total_yards': int,  # total number of yards gained on the catch
    'in_redzone': bool,  # whether the play took place in the redzone
    'pass_success': bool,  # whether the play gained enough yards to be considered a success
    'num_DBs': int, #number of cornerbacks and safeties on defense
}

'''
This class contains all functions of play parameters that I am interested in extracting from the data
'''
class PlayDescription:
    def __init__(self, Data, gameId, playId):
        #anything that I use more than twice in the functions I put here
        #since only one game and one play I can subset the tracking and plays data
        self.Data = Data
        self.gameId = gameId
        self.playId = playId

        self.Data.load_tracking_info(str(self.gameId))
        tracking_dic = self.Data.get_tracking_info()
        self.tracking_df = tracking_dic[str(self.gameId)]
        self.tracking_df = self.tracking_df.loc[self.tracking_df['playId'] == int(self.playId)]

        all_plays_df = self.Data.get_all_plays_info()
        self.play_row = all_plays_df.loc[(all_plays_df['gameId'] == int(self.gameId)) & (all_plays_df['playId'] == int(self.playId))]

        self.players_df = self.Data.get_all_players_info()

        #filter plays to make sure they have proper frame notes, added any bad plays to issues and updated them accordingly in the files, sometimes by deleting
        try:
            self.pass_forward_frame = int(self.tracking_df.loc[(self.tracking_df['event'] == 'pass_forward') | (self.tracking_df['event'] == 'pass_shovel') , 'frame.id'].mean())
            self.ball_snap_frame = int(self.tracking_df.loc[self.tracking_df['event'] == 'ball_snap', 'frame.id'].mean())
            self.outcome_frame = int(self.tracking_df.loc[(self.tracking_df['event'] == 'pass_outcome_incomplete') | (self.tracking_df['event'] == 'pass_outcome_caught') | (self.tracking_df['event'] == 'pass_outcome_interception') | (self.tracking_df['event'] == 'pass_outcome_touchdown') | (self.tracking_df['event'] == 'qb_sack'), 'frame.id'].mean())
        except ValueError:
            raise Exception("Cannot find one of these frames in gameId {} and playId {}!".format(str(self.gameId), str(self.playId)))

        #read in Routes.csv for the game and play
        filename = r'C:\Users\Mitch\Documents\UofM\Fall 2018\NFL\Data\Routes_{}.csv'.format(str(self.gameId))
        self.route_df = pd.read_csv(filename)
        self.route_df = self.route_df.loc[self.route_df['playId'] == int(self.playId)]

        #ball info is used a lot so include it in init
        ball_df = self.tracking_df.loc[self.tracking_df['team'] == 'ball']
        self.xy_ball_outcome = ball_df.loc[ball_df['frame.id'] == self.outcome_frame, ['x', 'y']]
        self.xy_ball_snap = ball_df.loc[ball_df['frame.id'] == self.ball_snap_frame, ['x', 'y']]

        self.closest_nflId = -1


    def is_blitz(self):
        #define a blitz by if more than 5 defenders are across the line of scrimmage when the ball is thrown

        #get 'x' of ball from tracking_df at los and pass_forward
        ball_df = self.tracking_df.loc[self.tracking_df['team'] == 'ball']
        try:
            x_ball_los = float(ball_df['x'].iloc[0])
        except IndexError:
            raise Exception("Empty ball dataframe for gameId {} and playId {}!".format(str(self.gameId), str(self.playId)))
        x_ball_pass = float(ball_df.loc[ball_df['frame.id'] == self.pass_forward_frame, 'x'])
        defense_is_rushing_left = x_ball_pass < x_ball_los

        #get defenders df and reduce to pass_forward frame x positions
        defenders_series = self.tracking_df.loc[(self.tracking_df['side'] == 'defense') & (self.tracking_df['frame.id'] == self.pass_forward_frame), 'x']

        #count number of blitzers past los
        if defense_is_rushing_left:
            number_of_blitzers = int(np.sum([x_def < x_ball_los for x_def in defenders_series]))
        else:
            number_of_blitzers = int(np.sum([x_def > x_ball_los for x_def in defenders_series]))

        #define blitz as atleast 5 rushers
        return number_of_blitzers >= 5

    def is_wr_bunch(self):
        #define a wr_bunch if there are atleast 3wr on the field and three of them are within 2 yards of each other

        #get df with just wr at ball snap with 'x' and 'y' columns
        wr_df = self.tracking_df.loc[(self.tracking_df['frame.id'] == self.ball_snap_frame) & (self.tracking_df['position'] == 'WR'), ['x', 'y']]

        #if theres not enough wrs to form a bunch return false
        if len(wr_df.index) < 3:
            return False

        #get distance matrix for each pair of receivers x2
        coord_mat = [[row['x'], row['y']] for _,row in wr_df.iterrows()]
        dist_mat = distance_matrix(coord_mat, coord_mat, p=2)

        #if atleast six of the distances are close to each other that means that there are 3 wide receivers within two yards of each other
        #have to add the diagonal since those will be zeroes
        out = np.where(dist_mat <= 3)[0]
        return len(out) >= (6 + len(wr_df.index))

    def is_success(self):
        #was the pass completed or or touchdown

        #get all rows where pass completed or touchdown, if none exist then pass was not completed otherwise pass was completed
        pass_completed_df = self.tracking_df.loc[(self.tracking_df['event'] == 'pass_outcome_caught') | (self.tracking_df['event'] == 'pass_outcome_touchdown')]

        return not bool(pass_completed_df.empty)

    def is_shotgun(self):
        #from play dataset was the quarterback lined up in shotgun

        offense_formation = self.play_row['offenseFormation'].to_string(index=False)
        comments = self.play_row['playDescription'].to_string(index=False)
        return bool('SHOTGUN' in offense_formation or 'Shotgun' in comments)

    def time_between_snap_and_pass(self):
        return int(self.pass_forward_frame - self.ball_snap_frame)

    def time_between_pass_and_outcome(self):
        return int(self.outcome_frame - self.pass_forward_frame)

    def closest_defenders_at_pass(self):
        #get information such as speed, distance, positions and direction relative to outcome point of closest defenders at pass

        #get 'x' 'y' of defenders at forward pass
        xy_defenders_pass = self.tracking_df.loc[(self.tracking_df['frame.id'] == self.pass_forward_frame) & (self.tracking_df['side'] == 'defense'), ['x', 'y', 'nflId', 'position']]

        #get nflIds of two closest defenders to outcome point at pass
        ball_coords_mat = [[float(self.xy_ball_outcome['x']), float(self.xy_ball_outcome['y'])]]
        defenders_coords_mat = [[row['x'], row['y']] for _, row in xy_defenders_pass.iterrows()]
        dist_mat = distance_matrix(ball_coords_mat, defenders_coords_mat, p=2)
        argpart = np.argpartition(dist_mat, kth=2, axis=None)[:2]
        closest_defenders = xy_defenders_pass.iloc[argpart]


        #get speeds, distances, and directions when pass thrown
        speeds = [float(self.tracking_df.loc[(self.tracking_df['frame.id'] == self.pass_forward_frame) & (self.tracking_df['nflId'] == row['nflId']), 's']) for _, row in closest_defenders.iterrows()]
        distances = [dist_mat[0, a] for a in argpart]
        tracking_directions = [float(self.tracking_df.loc[(self.tracking_df['frame.id'] == self.pass_forward_frame) & (self.tracking_df['nflId'] == row['nflId']), 'dir']) for _, row in closest_defenders.iterrows()]
        player_direction_points = [(5*np.cos(np.deg2rad(tracking_directions[i] + 90)) + closest_defenders['x'].iloc[i], 5*np.sin(np.deg2rad(tracking_directions[i] + 90)) + closest_defenders['y'].iloc[i]) for i in range(len(tracking_directions))]

        angles_relative_to_outcome = self._angle_relative_to_([(row['x'], row['y']) for _, row in self.xy_ball_outcome.iterrows()], [(row['x'], row['y']) for _, row in closest_defenders.iterrows()], player_direction_points)

        positions = [(row['nflId'], row['position']) for _, row in closest_defenders.iterrows()]

        #get height, weight top speed
        nflIds_ = [row['nflId'] for _, row in closest_defenders.iterrows()]
        h_w_s = [(row['nflId'], row['Height_cm'], row['Weight'], row['top_speed']) for _, row in self.players_df.loc[self.players_df['nflId'].isin(nflIds_)].iterrows()]

        def _sort_helper(elem):
            out_helper = [i for i in range(len(positions)) if int(positions[i][0]) == int(elem[0])]
            return out_helper[0]

        h_w_s.sort(key=_sort_helper)

        #possibly add in angle compared to wide receiver the wide receiver??
        return [(float(speeds[i]), round(float(distances[i]), 3), round(float(angles_relative_to_outcome[i]), 3), str(positions[i][1]),
                 float(h_w_s[i][1]), float(h_w_s[i][2]), float(h_w_s[i][3])) for i in range(len(angles_relative_to_outcome))]

    def closest_defenders_at_outcome(self):
        #get information such as speed, distance, positions and direction relative to outcome point of closest defenders at outcome

        #get 'x' 'y' of defenders at outcome
        xy_defenders_outcome = self.tracking_df.loc[(self.tracking_df['frame.id'] == self.outcome_frame) & (self.tracking_df['side'] == 'defense'), ['x', 'y', 'nflId', 'position']]

        #get nflIds of two closest defenders
        ball_coords_mat = [[float(self.xy_ball_outcome['x']), float(self.xy_ball_outcome['y'])]]
        defenders_coords_mat = [[row['x'], row['y']] for _, row in xy_defenders_outcome.iterrows()]
        dist_mat = distance_matrix(ball_coords_mat, defenders_coords_mat, p=2)
        argpart = np.argpartition(dist_mat, kth=2, axis=None)[:2]
        closest_defenders = xy_defenders_outcome.iloc[argpart]


        #get speeds, distances, and directions when outcome occurs
        speeds = [float(self.tracking_df.loc[(self.tracking_df['frame.id'] == self.outcome_frame) & (self.tracking_df['nflId'] == row['nflId']), 's']) for _, row in closest_defenders.iterrows()]
        distances = [dist_mat[0, a] for a in argpart]
        tracking_directions = [float(self.tracking_df.loc[(self.tracking_df['frame.id'] == self.outcome_frame) & (self.tracking_df['nflId'] == row['nflId']), 'dir']) for _, row in closest_defenders.iterrows()]
        player_direction_points = [(5*np.cos(np.deg2rad(tracking_directions[i] + 90)) + closest_defenders['x'].iloc[i], 5*np.sin(np.deg2rad(tracking_directions[i] + 90)) + closest_defenders['y'].iloc[i]) for i in range(len(tracking_directions))]

        angles_relative_to_outcome = self._angle_relative_to_([(row['x'], row['y']) for _, row in self.xy_ball_outcome.iterrows()], [(row['x'], row['y']) for _, row in closest_defenders.iterrows()], player_direction_points)

        positions = [(row['nflId'], row['position']) for _, row in closest_defenders.iterrows()]

        #get height, weight top speed
        nflIds_ = [row['nflId'] for _, row in closest_defenders.iterrows()]
        h_w_s = [(row['nflId'], row['Height_cm'], row['Weight'], row['top_speed']) for _, row in self.players_df.loc[self.players_df['nflId'].isin(nflIds_)].iterrows()]

        def _sort_helper(elem):
            out_helper = [i for i in range(len(positions)) if int(positions[i][0]) == int(elem[0])]
            return out_helper[0]

        h_w_s.sort(key=_sort_helper)

        #possibly add in angle compared to wide receiver the wide receiver??
        return [(float(speeds[i]), round(float(distances[i]), 3), round(float(angles_relative_to_outcome[i]), 3), str(positions[i][1]),
                 float(h_w_s[i][1]), float(h_w_s[i][2]), float(h_w_s[i][3])) for i in range(len(angles_relative_to_outcome))]

    def _angle_relative_to_(self, xy_ball, xy_defenders, xy_player_directions):
        #this function uses the law of cosines to analyze the angle between where the outcome point is and the players direction
        #uses three points: the defenders xy, the outcome point xy, and the xy of five yards in front of where the defender is travelling
        A_list, B_list, C_list = [], [], []

        for i in range(len(xy_player_directions)):
            A_list.append(np.linalg.norm(np.subtract(xy_player_directions[i], xy_ball[0])))
            B_list.append(np.linalg.norm(np.subtract(xy_defenders[i], xy_player_directions[i])))
            C_list.append(np.linalg.norm(np.subtract(xy_defenders[i], xy_ball[0])))


        A_angles_list = [np.rad2deg(np.arccos((-A_list[ind]**2 + B_list[ind]**2 + C_list[ind]**2)/(2*B_list[ind]*C_list[ind]))) for ind in range(len(A_list))]

       # [self._draw_player_triangles(A_angles_list[i], xy_ball[0], xy_defenders[i], xy_player_directions[i]) for i in range(len(A_list))]

        return A_angles_list

    def _draw_player_triangles(self, angle, xy_ball, xy_defender, xy_player_direction):
        #function for error checking, make sure the angles match what the picture shows

        def midpoint(p1, p2):
            return ((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2)

        all_coords = [xy_ball, xy_defender, xy_player_direction]
        xs = [t[0] for t in all_coords]
        ys = [t[1] for t in all_coords]

        centers = [midpoint(xy_ball, xy_player_direction), midpoint(xy_defender, xy_player_direction), midpoint(xy_defender, xy_ball)]

        fig, ax = plt.subplots()
        ax.set_title('A angle is {}'.format(angle))

        ax.plot(xs + [xs[0]], ys + [ys[0]], 'b-')

        for i, txt in enumerate(['A', 'B', 'C']):
            ax.annotate(txt, centers[i])

        for i, txt in enumerate(['ball', 'defender', 'direction']):
            ax.annotate(txt, (xs[i], ys[i]))

        fig.show()
        plt.close()


    def intended_receiver_info(self):
        #get information such as speed, distance, and angle between direction and ball point of intended receiver at the pass

        #get xy of all possible receivers in play
        wrterb_df = self.tracking_df.loc[((self.tracking_df['position'] == 'WR') | (self.tracking_df['position'] == 'TE') | (self.tracking_df['position'] == 'RB') | (self.tracking_df['position'] == 'FB')) & (self.tracking_df['frame.id'] == self.outcome_frame)]
        wrterb_xy = [(row['nflId'], row['x'], row['y']) for _, row in wrterb_df.iterrows()]

        ball_xy = (float(self.xy_ball_outcome['x']), float(self.xy_ball_outcome['y']))

        #find closest wide receiver at ball outcome
        closest_nflId_dist = min([(wrterb[0], np.linalg.norm(np.subtract(wrterb[1:], ball_xy))) for wrterb in wrterb_xy], key= lambda t: t[1])
        self.closest_nflId = closest_nflId_dist[0]

        #get row of intented receiver at pass forward
        closest_receiver_at_pass_forward = self.tracking_df.loc[(self.tracking_df['nflId'] == self.closest_nflId) & (self.tracking_df['frame.id'] == self.pass_forward_frame)]

        #get distance to ball, speed and direction of receiver
        distance = closest_nflId_dist[1]
        speed = float(closest_receiver_at_pass_forward['s'])
        player_angle = float(closest_receiver_at_pass_forward['dir'])

        angle_relative_to_ball = self._angle_relative_to_([ball_xy], [(float(closest_receiver_at_pass_forward['x']), float(closest_receiver_at_pass_forward['y']))],
        [(5*np.cos(np.deg2rad(player_angle + 90)) + float(closest_receiver_at_pass_forward['x']), 5*np.sin(np.deg2rad(player_angle + 90)) + float(closest_receiver_at_pass_forward['y']))])[0]

        #height weight and top_speed
        h_w_s = (float(self.players_df.loc[self.players_df['nflId'] == self.closest_nflId, 'Height_cm']), float(self.players_df.loc[self.players_df['nflId'] == self.closest_nflId, 'Weight']), float(self.players_df.loc[self.players_df['nflId'] == self.closest_nflId, 'top_speed']))

        return (float(speed), round(float(distance), 3), round(float(angle_relative_to_ball), 3), h_w_s[0], h_w_s[1], h_w_s[2])

    def all_receivers_info(self):
        #get positions of receivers
        positions = [(row['nflId'], row['position']) for _, row in self.route_df.iterrows()]
        #get routes of receivers
        route_names = [(row['nflId'], row['route']) for _, row in self.route_df.iterrows()]

        #get difference between receiver 'y' and ball 'y' at snap
        ball_xy = (float(self.xy_ball_snap['x']), float(self.xy_ball_snap['y']))
        wrterb_snap = self.tracking_df.loc[((self.tracking_df['position'] == 'WR') | (self.tracking_df['position'] == 'TE')| (self.tracking_df['position'] == 'RB')| (self.tracking_df['position'] == 'FB'))
        & (self.tracking_df['frame.id'] == self.ball_snap_frame)]
        wrterb_xy = [(row['nflId'], row['x'], row['y']) for _, row in wrterb_snap.iterrows()]

        #need post snap ball 'x' to see which way receivers are running
        ball_post_snap_x = float(self.tracking_df.loc[(self.tracking_df['team'] == 'ball') & (self.tracking_df['frame.id'] == (self.ball_snap_frame + 2)), 'x'])
        running_right = ball_xy[0] > ball_post_snap_x
        #positive distances for receivers lined up on the left side of the ball, negative on the right
        differences = [(wrterb_tup[0], ((-1)**running_right)*np.subtract(ball_xy[1], wrterb_tup[2])) for wrterb_tup in wrterb_xy]

        #get linestrings of routes run and get the max 'x'
        route_linestrings = [(row['nflId'], s2LS(row['LineString'])) for _, row in self.route_df.iterrows()]
        route_depths = [(t[0], t[1].bounds[2]-t[1].bounds[0]) for t in route_linestrings]

        #nfl players height, weight and top speed
        nflIds_ = self.route_df['nflId'].tolist()
        h_w_s = [(row['nflId'], row['Height_cm'], row['Weight'], row['top_speed']) for _, row in self.players_df.loc[self.players_df['nflId'].isin(nflIds_)].iterrows()]

        def _sort_helper(elem):
            out_helper = [d[1] for d in differences if int(d[0]) == int(elem[0])]
            return out_helper[0]

        #sort the data based on where the receivers were lined up at the pass
        positions.sort(key=_sort_helper, reverse=True)
        route_names.sort(key=_sort_helper, reverse=True)
        differences.sort(key=lambda x: x[1], reverse=True)
        route_depths.sort(key=_sort_helper, reverse=True)
        h_w_s.sort(key=_sort_helper, reverse=True)

        out_list = [(positions[i][1], route_names[i][1], round(float(differences[i][1]), 3), round(float(route_depths[i][1]), 3), bool(route_depths[i][0] == self.closest_nflId),
                     float(h_w_s[i][1]), float(h_w_s[i][2]), float(h_w_s[i][3])) for i in range(len(positions))]

        while len(out_list) < 5:
            sign_list = [np.sign(ol[2]) for ol in out_list]
            try:
                first_1 = sign_list.index(1)
            except ValueError:
                first_1 = 0
            out_list.insert(first_1, ('TE', 'bubble.block', 0, 0, False, 182.8, 250.0, 8.0))

        out_list.sort(key= lambda t: t[2], reverse=True)

        return out_list

    def closest_sideline(self):
        #get outcome spot of ball
        ball_x, ball_y = float(self.xy_ball_outcome['x']), float(self.xy_ball_outcome['y'])
        right_sideline = 0 #compare to y
        left_sideline = 53.3 #compare to y
        left_endzone = 0 #compare to x
        right_endzone = 120 #compare to x
        return round(float(min(abs(ball_x - left_endzone), abs(ball_x - right_endzone), abs(ball_y - right_sideline), abs(ball_y - left_sideline))), 3)

    def count_num_close_defenders(self, radius=3):
        #find how many defenders were in close proximity to quarterback at pass thrown

        try:
        #get 'x' 'y' of qb at pass forward
            qb_pass_xy = [float(self.tracking_df.loc[(self.tracking_df['frame.id'] == self.pass_forward_frame) & (self.tracking_df['position'] == 'QB'), 'x']),
                      float(self.tracking_df.loc[(self.tracking_df['frame.id'] == self.pass_forward_frame) & (self.tracking_df['position'] == 'QB'), 'y'])]
        except TypeError:
            raise Exception("There is something wrong with QB x or y in gameId {} and playId {}!".format(str(self.gameId), str(self.playId)))

        #get dataframe of defenders 'x' 'y'
        defenders_xy_df = self.tracking_df.loc[(self.tracking_df['frame.id'] == self.pass_forward_frame) & (self.tracking_df['side'] == 'defense'), ['x', 'y']]
        defenders_xy = [(row['x'], row['y']) for _, row in defenders_xy_df.iterrows()]

        #get distance between each and return the distances less than radius
        distances = [np.linalg.norm(np.subtract(qb_pass_xy, d_xy)) for d_xy in defenders_xy]
        return int(np.sum([d < radius for d in distances]))

    def in_pocket(self):
        #define the pocket as about 3 yards to either side of where the ball was snapped from and see if qb was within this range when ball was thrown
        y_ball_snap = [row['y'] for _, row in self.xy_ball_snap.iterrows()][0]
        y_qb = float(self.tracking_df.loc[(self.tracking_df['frame.id'] == self.pass_forward_frame) & (self.tracking_df['position'] == 'QB'), 'y'])

        if np.isnan(y_ball_snap):
            y_qb_snap = float(self.tracking_df.loc[(self.tracking_df['frame.id'] == self.ball_snap_frame) & (self.tracking_df['position'] == 'QB'), 'y'])
            return bool(not ((y_qb < y_qb_snap - 3 and y_qb < y_qb_snap + 3) or (y_qb > y_qb_snap - 3 and y_qb > y_qb_snap + 3)))


        return bool(not ((y_qb < y_ball_snap - 3 and y_qb < y_ball_snap + 3) or (y_qb > y_ball_snap - 3 and y_qb > y_ball_snap + 3)))

    def thrown_away(self):
        #based on how far away the intended receiver was when the ball gets to the outcome point

        xy_closest_receiver = (float(self.tracking_df.loc[(self.tracking_df['nflId'] == self.closest_nflId) & (self.tracking_df['frame.id'] == self.outcome_frame), 'x']),
                               float(self.tracking_df.loc[(self.tracking_df['nflId'] == self.closest_nflId) & (self.tracking_df['frame.id'] == self.outcome_frame), 'y']))
        xy_ball = (float(self.xy_ball_outcome['x']), float(self.xy_ball_outcome['y']))

        dist = np.linalg.norm(np.subtract(xy_ball, xy_closest_receiver))

        return bool(dist >= 7) #7 pretty arbitrary

    def yards_after_catch(self):
        #this is directly from the plays.csv
        if self.is_success():
            try:
                return int(self.play_row['YardsAfterCatch'])
            except ValueError:
                return 0
        return 0

    def total_yards(self):
        return int(self.play_row['PlayResult'])

    def in_redzone(self):
        out = (self.play_row['possessionTeam'].to_string(index=False, header=False) != self.play_row['yardlineSide'].to_string(index=False, header=False)) and (float(self.play_row['yardlineNumber']) <= 20)
        return out

    def pass_success(self):
        #based on certain conditions determine if the play should be considered a success

        # if the ball was not completed no success
        if not self.is_success():
            return False

        if self.play_row['PassResult'].to_string(index=False, header=False) == 'IN':
            return False

        if int(self.play_row['down']) == 1:
            #if team gained at least half of what it takes to have a better than 50% chance of making first down or atleast 8 yards consider a success
            if float(self.play_row['PlayResult']) >= 8:
                return True
            if float(self.play_row['PlayResult']) >= (float(self.play_row['yardsToGo']) - 4)/2:
                return True
            return False

        elif int(self.play_row['down']) == 2:
            #if team gained enough yards to make it aleast a 3rd and 4 or atleast 8 yards consider a success
            if float(self.play_row['PlayResult']) >= 8:
                return True
            if float(self.play_row['yardsToGo']) - float(self.play_row['PlayResult']) <= 4:
                return True
            return False

        elif int(self.play_row['down']) == 3:
            #if the team converts consider a success
            return float(self.play_row['PlayResult']) >= float(self.play_row['yardsToGo'])

        elif int(self.play_row['down']) == 4:
            #if the team converts consider a sucess
            return (float(self.play_row['PlayResult']) >= float(self.play_row['yardsToGo'])) and (not self.play_row['isSTPlay'].to_string(index=False, header=False) == 'TRUE')

        else:
            #down equals zero - no sucess
            return False

    def num_DBs(self):
        #directly from the plays.csv, how many defensive backs were on the field
        d_personnel = self.play_row['personnel.defense'].to_string(index=False, header=False)
        if d_personnel == 'NA':
            return -1
        d_personnel_list = d_personnel.rstrip().split(' ')
        try:
            db_ind = d_personnel_list.index('DB')
        except ValueError:
            return 0
        return int(d_personnel_list[db_ind - 1])

from load_data import Data

def define_plays(pI_gI_dic):
    #function for parallelizing the code, updates the dictionary of play parameters

    pI = pI_gI_dic[0]
    gI = pI_gI_dic[1]
    dic_of_plays = pI_gI_dic[2]

    data = Data()

    # pass if play resulted in a sack
    plays = data.get_all_plays_info()
    if plays.loc[(plays['gameId'] == int(gI)) & (plays['playId'] == int(pI)), 'PassResult'].to_string(index=False).strip() == 'S':
        return {'playId': 'sack'}

    PD = PlayDescription(Data=data, gameId=gI, playId=pI)

    play_description_dic = {
        'gameId': int(gI),  # the gameId
        'playId': int(pI),  # the playId
        'blitz': PD.is_blitz() if 'blitz' not in dic_of_plays else dic_of_plays['blitz'],  # whether 5+ players rush the quarterback
        'wr_bunch': PD.is_wr_bunch() if 'wr_bunch' not in dic_of_plays else dic_of_plays['wr_bunch'],  # whether three wr are together in a triangle
        'success': PD.is_success() if 'success' not in dic_of_plays else dic_of_plays['success'],  # was the play a success complete pass or a drop
        'shotgun': PD.is_shotgun() if 'shotgun' not in dic_of_plays else dic_of_plays['shotgun'],  # was the quarterback in shotgun
        'time_between_snap_and_pass': PD.time_between_snap_and_pass() if 'time_between_snap_and_pass' not in dic_of_plays else dic_of_plays['time_between_snap_and_pass'],  # num frames between when the ball was snapped and when the ball was thrown
        'time_between_pass_and_outcome': PD.time_between_pass_and_outcome() if 'time_between_pass_and_outcome' not in dic_of_plays else dic_of_plays['time_between_pass_and_outcome'], # num frames between when the ball was passed and the outcome
        'closest_defenders_at_pass': PD.closest_defenders_at_pass(), # if 'closest_defenders_at_pass' not in dic_of_plays else dic_of_plays['closest_defenders_at_pass'],  # seven element tuples of speed, distance, direction relative to outcome point when ball is thrown, position, height, weight, speed
        'closest_defenders_at_outcome': PD.closest_defenders_at_outcome(), # if 'closest_defenders_at_outcome' not in dic_of_plays else dic_of_plays['closest_defenders_at_outcome'],  # seven element tuples of speed, distance and direction relative to outcome point at outcome, position, height, weight, speed
        'intended_receiver_at_pass': PD.intended_receiver_info(), # if 'intended_receiver_at_pass' not in dic_of_plays else dic_of_plays['intended_receiver_at_pass'],  # three element tuple of speed, distance, direction relative to outcome point when ball is thrown
        'all_receivers': PD.all_receivers_info(), # if 'all_receivers' not in dic_of_plays else dic_of_plays['all_receivers'],   # eight element tuples of position, route, distance from ball at snap, depth, intended receiver, height, weight, speed
        'closest_sideline': PD.closest_sideline() if 'closest_sideline' not in dic_of_plays else dic_of_plays['closest_sideline'],  # distance to the closest sideline from outcome point
        'num_closest_defenders_to_qb': PD.count_num_close_defenders() if 'num_closest_defenders_to_qb' not in dic_of_plays else dic_of_plays['num_closest_defenders_to_qb'],  #number of defenders close to quarterback at pass to measure pressure
        'in_pocket': PD.in_pocket() if 'in_pocket' not in dic_of_plays else dic_of_plays['in_pocket'], #whether the quarterback passed outside of the pocket
        'thrown_away': PD.thrown_away() if 'thrown_away' not in dic_of_plays else dic_of_plays['thrown_away'], #whether the quarterback threw the ball away
        'yards_after_catch': PD.yards_after_catch() if 'yards_after_catch' not in dic_of_plays else dic_of_plays['yards_after_catch'], #how many yards gained after the catch
        'total_yards': PD.total_yards() if 'total_yards' not in dic_of_plays else dic_of_plays['total_yards'], #total number of yards gained on the catch
        'in_redzone': PD.in_redzone() if 'in_redzone' not in dic_of_plays else dic_of_plays['in_redzone'], #whether the play took place in the redzone
        'pass_success': PD.pass_success() if 'pass_success' not in dic_of_plays else dic_of_plays['pass_success'], #whether the play gained enough yards to be considered a success
        'num_DBs': PD.num_DBs() if 'num_DBs' not in dic_of_plays else dic_of_plays['num_DBs'], #number of cornerbacks and safeties on defense
    }

    del PD
    del data

    return play_description_dic

if __name__ == '__main__':
    #if this file is run update the play parameters for each game and play
    #issues would come up from time to time and I would note those in issues.txt
    import os
    import json
    import multiprocessing as mp

    dir_path = r'C:\Users\Mitch\Documents\UofM\Fall 2018\NFL\Data'
    all_gameIds = [f[-14:-4] for f in os.listdir(dir_path) if f[0] == 'R']

    all_play_descriptions = []
    cpus = mp.cpu_count()

    for gI in all_gameIds:
        rf = r'C:\Users\Mitch\Documents\UofM\Fall 2018\NFL\Data\Routes_{}.csv'.format(gI)
        route_df = pd.read_csv(rf)
        unique_playIds = route_df['playId'].unique()

        output_filename = r'C:\Users\Mitch\Documents\UofM\Fall 2018\NFL\Data\play_descriptions_{}.json'.format(str(gI))

        if os.path.exists(output_filename):

        # uncomment this if have to stop half way through
        #     continue
        # current_list_of_dics = [{}] * len(unique_playIds)

        # uncomment this is starting over
            with open(output_filename) as json_file:
                current_list_of_dics = json.load(json_file)
        else:
            current_list_of_dics = [{}]*len(unique_playIds)

        pI_gI_dic_list = [(pI, gI, dic) for pI, dic in zip(unique_playIds, current_list_of_dics)]

        #parallelize to increase speed
        with mp.Pool(processes=cpus-1) as pool:
            parallel_out = pool.map(define_plays, pI_gI_dic_list)

        with open(output_filename, 'w') as of:
            json.dump(parallel_out, of)


