import os
import pandas as pd
import multiprocessing as mp
import numpy as np
'''
This class is used to load the initial datasets. It also has functions to edit the tracking information to include positions and whether the player is on offense or defense.
I also use this class to segment the wide receiver and running back tracking information from the data on all pass plays.
'''
class Data:
    def __init__(self):
        self.header = r'C:\Users\Mitch\Documents\UofM\Fall 2018\NFL\Data'
        self.game_info = pd.read_csv(os.path.join(self.header, r'games.csv'))
        self.players_info = pd.read_csv(os.path.join(self.header, r'players.csv'))
        self.plays_info = pd.read_csv(os.path.join(self.header, r'plays.csv'))
        self.tracking_info = {}
        self.wr_te_routes = {}
        self.rb_routes = {}

        self.top_speed_dic = {}
        for _, nflId in self.players_info['nflId'].iteritems():
            self.top_speed_dic[nflId] = -1

    def get_all_games_info(self):
        return self.game_info

    def get_all_players_info(self):
        return self.players_info

    def get_all_plays_info(self):
        return self.plays_info

    def get_tracking_info(self):
        return self.tracking_info

    def _assign_offense(self, team, playId, gameId):
        #Get the home team from games.csv and get who has possesion of the ball from plays.csv
        #Then based on if the player's team is home or away in the tracking info assign offense or defense

        #the ball rows in the tracking info needs its own team
        if team == 'ball':
            return 'ball'
        homeTeam = self.game_info.loc[self.game_info['gameId'] == int(gameId), 'homeTeamAbbr'].to_string(header=False, index=False)
        possTeam = self.plays_info.loc[(self.plays_info['playId'] == int(playId)) & (self.plays_info['gameId'] == int(gameId)), 'possessionTeam'].to_string(header=False, index=False)
        flag = homeTeam == possTeam
        if flag:
            out = 'offense' if team == 'home' else 'defense'
        else:
            out = 'offense' if team == 'away' else 'defense'
        return out

    def _assign_position(self, nflId):
        #get position of player from players.csv
        if np.isnan(nflId):
            return 'ball'
        return self.players_info.loc[self.players_info['nflId'] == nflId, 'PositionAbbr'].to_string(index=False)

    def assign_top_speed(self):
        tracking_files = [f for f in os.listdir(self.header) if f[0] == 't']

        for tf in tracking_files:
            temp_df = pd.read_csv(os.path.join(self.header, tf), header=0)
            temp_df = temp_df.loc[(temp_df['position'] == 'RB') | (temp_df['position'] == 'FB') | (temp_df['position'] == 'WR') | (temp_df['position'] == 'TE')
                                  | (temp_df['position'] == 'CB') | (temp_df['position'] == 'ILB')| (temp_df['position'] == 'OLB') | (temp_df['position'] == 'SS')
                                  | (temp_df['position'] == 'FS') | (temp_df['position'] == 'DE')]

            unique_nflIds = temp_df['nflId'].unique().astype(np.int32)
            for un in unique_nflIds:
                max_speed = temp_df.loc[temp_df['nflId'] == un, 's'].max()
                try:
                    if un in self.top_speed_dic:
                        self.top_speed_dic[un] = max(max_speed, self.top_speed_dic[un])
                except ValueError:
                    pass

        top_speed_df = pd.DataFrame({'nflId': list(self.top_speed_dic.keys()), 'top_speed': list(self.top_speed_dic.values())})
        self.players_info.merge(top_speed_df, left_on='nflId', right_on='nflId').to_csv(os.path.join(self.header, r'players.csv'), index=False)


    def _update_offense_and_position_wrapper(self, tf):
        #wrapper function to update tracking csv to include offense/defense and player positions
        temp_df = pd.read_csv(os.path.join(self.header, tf), header=0)
        gameId = tf[16:-4]
        #if 'side' not in temp_df.columns:
        temp_df['side'] = temp_df.apply(lambda row: self._assign_offense(row['team'], row['playId'], gameId), axis=1)
        #if 'position' not in temp_df.columns:
        temp_df['position'] = temp_df.apply(lambda row: self._assign_position(row['nflId']), axis=1)
        temp_df = self._drop_unamed_cols_in_tracking(temp_df)
        temp_df.to_csv(os.path.join(self.header, tf), index=False)

    def update_offense_and_position(self):
        #parrellize this

        #raise Exception("this has already been done.")
        cpus = mp.cpu_count()
        tracking_files = [f for f in os.listdir(self.header) if f[0] == 't']
        with mp.Pool(processes=cpus-1) as p:
            p.map(self._update_offense_and_position_wrapper, tracking_files)

    def _drop_unamed_cols_in_tracking(self, df):
        if isinstance(df, pd.DataFrame):
            bad_column_names = ['Unnamed: 0.1.1.1.1.1', 'Unnamed: 0.1.1.1.1', 'Unnamed: 0.1.1.1', 'Unnamed: 0.1.1', 'Unnamed: 0.1', 'Unnamed: 0']
            for bcn in bad_column_names:
                if bcn in df.columns.tolist():
                    df.drop(bcn, axis=1)
        return df

    def load_tracking_info(self, gameId):
        #get tracking.csv into python

        if not isinstance(gameId, str):
            raise Exception('gameId is not a string!')
        ending = r'tracking_gameId_{}.csv'.format(gameId)
        filename = os.path.join(self.header, ending)
        if os.path.exists(filename):
            self.tracking_info[gameId] = pd.read_csv(filename)
        else:
            raise Exception('file {} not found!'.format(filename))

    def load_wr_te_routes(self):
        #load wr and te routes under the restriction of a pass play and only up to 50 frames after the ball is snapped
        #many plays did not have frame notes like when the ball was snapped or when the ball was completed so I skipped those
        #create and save this new df to python object
        if not bool(self.tracking_info):
            raise Exception("No game has been loaded into tracking info!")
        for gameId, tracking_df in self.tracking_info.items():
            df = pd.DataFrame()
            red_plays = self.plays_info[self.plays_info['gameId'] == int(gameId)]
            playIds = red_plays.loc[(red_plays['PassResult'] == 'C') | (red_plays['PassResult'] == 'I') | (red_plays['PassResult'] == 'IN') | (red_plays['PassResult'] == 'S'), 'playId']
            for pI in playIds:
                red_tracking_df = tracking_df[(tracking_df['playId'] == int(pI)) & ((tracking_df['position'] == 'WR') | (tracking_df['position'] == 'TE'))]
                try:
                    ball_snap_frame = int(red_tracking_df.loc[red_tracking_df['event'] == 'ball_snap', 'frame.id'].mean())
                    outcome_frame = int(red_tracking_df.loc[(red_tracking_df['event'] == 'pass_outcome_incomplete') | (red_tracking_df['event'] == 'pass_outcome_caught') | (red_tracking_df['event'] == 'pass_outcome_interception') | (red_tracking_df['event'] == 'pass_outcome_touchdown') | (red_tracking_df['event'] == 'qb_sack'), 'frame.id'].mean())
                except ValueError:
                    continue
                end_frame = ball_snap_frame + 50 if outcome_frame - ball_snap_frame > 50 else outcome_frame
                ball_list_x = tracking_df.loc[(tracking_df['playId'] == int(pI)) & (tracking_df['position'] == 'ball') & (tracking_df['frame.id'] >= ball_snap_frame), 'x'].tolist()
                if not ball_list_x:
                    continue

                #update if the offense was running right
                running_right = ball_list_x[0] > ball_list_x[6]
                only_routes_df = red_tracking_df.loc[(red_tracking_df['frame.id'] >= ball_snap_frame) & (red_tracking_df['frame.id'] <= end_frame), ['nflId', 'gameId', 'playId', 'x', 'y', 'position']]
                only_routes_df.reset_index(drop=True, inplace=True)
                only_routes_df = only_routes_df.join(pd.DataFrame({'running_right': [running_right]*len(only_routes_df.index)}))
                df = df.append(only_routes_df)
            self.wr_te_routes[gameId] = df

    def load_rb_routes(self):
        if not bool(self.tracking_info):
            raise Exception("No game has been loaded into tracking info!")
        for gameId, tracking_df in self.tracking_info.items():
            df = pd.DataFrame()
            red_plays = self.plays_info[self.plays_info['gameId'] == int(gameId)]
            playIds = red_plays.loc[(red_plays['PassResult'] == 'C') | (red_plays['PassResult'] == 'I') | (red_plays['PassResult'] == 'IN') | (red_plays['PassResult'] == 'S'), 'playId']
            for pI in playIds:
                red_tracking_df = tracking_df[(tracking_df['playId'] == int(pI)) & ((tracking_df['position'] == 'RB') | (tracking_df['position'] == 'FB'))]
                if red_tracking_df.empty:
                    continue
                try:
                    ball_snap_frame = int(red_tracking_df.loc[red_tracking_df['event'] == 'ball_snap', 'frame.id'].mean())
                    outcome_frame = int(red_tracking_df.loc[(red_tracking_df['event'] == 'pass_outcome_incomplete') | (red_tracking_df['event'] == 'pass_outcome_caught') | (red_tracking_df['event'] == 'pass_outcome_interception') | (red_tracking_df['event'] == 'pass_outcome_touchdown') | (red_tracking_df['event'] == 'qb_sack'), 'frame.id'].mean())
                except ValueError:
                    continue
                end_frame = ball_snap_frame + 50 if outcome_frame - ball_snap_frame > 50 else outcome_frame
                ball_list_x = tracking_df.loc[(tracking_df['playId'] == int(pI)) & (tracking_df['position'] == 'ball') & (tracking_df['frame.id'] >= ball_snap_frame), 'x'].tolist()
                if not ball_list_x:
                    continue
                running_right = ball_list_x[0] > ball_list_x[6]
                only_routes_df = red_tracking_df.loc[(red_tracking_df['frame.id'] >= ball_snap_frame) & (red_tracking_df['frame.id'] <= end_frame), ['nflId', 'gameId', 'playId', 'x', 'y', 'position']]
                only_routes_df.reset_index(drop=True, inplace=True)
                only_routes_df = only_routes_df.join(pd.DataFrame({'running_right': [running_right]*len(only_routes_df.index)}))
                df = df.append(only_routes_df)
            self.rb_routes[gameId] = df

    def get_wr_te_routes(self):
        return self.wr_te_routes

    def get_rb_routes(self):
        return self.rb_routes


if __name__ == '__main__':
    #if this is run update the offense/defense and positions
    import datetime
    now = datetime.datetime.now()
    print(now)
    data = Data()
    # data.update_offense_and_position()
    data.assign_top_speed()
    now = datetime.datetime.now()
    print(now)
