import json
import os
import pandas as pd
import numpy as np

DATA_directory = r'C:\Users\Mitch\Documents\UofM\Fall 2018\NFL\Data'
all_json_filenames = [f for f in os.listdir(DATA_directory) if f[-5:] == '.json']
single_items = ['gameId', 'playId', 'blitz', 'wr_bunch', 'success', 'shotgun', 'time_between_snap_and_pass', 'time_between_pass_and_outcome',
                'closest_sideline', 'num_closest_defenders_to_qb', 'in_pocket', 'thrown_away', 'yards_after_catch', 'total_yards',
                'in_redzone', 'pass_success', 'num_DBs',]

def remove_sacks(play_list):
    #takes in the list of play dictionaries for the game, removes any play where there was a sack
    out_list = [PD for PD in play_list if PD['playId'] != 'sack']
    return out_list

#check for NaN or None in final dataframe

#need to think of a way to determine if the ball data is bad. some plays it looked like ball was a bit laggy or just in the wrong spots
#good website: https://newonlinecourses.science.psu.edu/stat504/node/150/

#stuff to look at that might be a little odd:
#   -large distance of defenders at outcome.
#   -number of rb's and te's being 3 or more
#   -might need to rethink how wide receivers data is organized, not clear if the predictors should
#    just have the wide receivers in order or if the predictors should be for a certain area and left blank if no wide receiver in that area

if __name__ == '__main__':

    for num_doc, j_file in enumerate(all_json_filenames):
        with open(os.path.join(DATA_directory, j_file), 'r') as jsonfile:
            current_play_list = json.load(jsonfile)

            current_play_list = remove_sacks(current_play_list)

            game_rows = []
            gameId = current_play_list[0]['gameId']

            for play_dic in current_play_list:
                new_row = {}

                for key, value in play_dic.items():
                    #if the data is not a list then it is fine to just add to the new dataframe
                    #otherwise need to do some organizing
                    if key in single_items:
                        new_row[key] = value
                        continue

                    #going to do WR 1, 2, 3 on the right and left and two middle positions. Will also include distance just in case I want to do something different
                    #in the future. Will make it fit as best I can.

                    if key == 'closest_defenders_at_pass':
                        for tup_ind, tup in enumerate(play_dic[key]):
                            new_row['cDef{}_pass_distance'.format(tup_ind + 1)] = tup[1]
                            new_row['cDef{}_pass_speed'.format(tup_ind + 1)] = tup[0]
                            new_row['cDef{}_pass_direction'.format(tup_ind + 1)] = tup[2]
                            new_row['cDef{}_pass_position'.format(tup_ind + 1)] = tup[3]
                            new_row['cDef{}_pass_height'.format(tup_ind + 1)] = tup[4]
                            new_row['cDef{}_pass_weight'.format(tup_ind + 1)] = tup[5]
                            new_row['cDef{}_pass_top_speed'.format(tup_ind + 1)] = tup[6]

                    if key == 'closest_defenders_at_outcome':
                        for tup_ind, tup in enumerate(play_dic[key]):
                            new_row['cDef{}_outcome_distance'.format(tup_ind + 1)] = tup[1]
                            new_row['cDef{}_outcome_speed'.format(tup_ind + 1)] = tup[0]
                            new_row['cDef{}_outcome_direction'.format(tup_ind + 1)] = tup[2]
                            new_row['cDef{}_outcome_position'.format(tup_ind + 1)] = tup[3]
                            new_row['cDef{}_outcome_height'.format(tup_ind + 1)] = tup[4]
                            new_row['cDef{}_outcome_weight'.format(tup_ind + 1)] = tup[5]
                            new_row['cDef{}_outcome_top_speed'.format(tup_ind + 1)] = tup[6]

                    if key == 'intended_receiver_at_pass':
                        new_row['intended_distance'] = play_dic[key][1]
                        new_row['intended_speed'] = play_dic[key][0]
                        new_row['intended_direction'] = play_dic[key][2]
                        new_row['intended_height'] = play_dic[key][3]
                        new_row['intended_weight'] = play_dic[key][4]
                        new_row['intended_top_speed'] = play_dic[key][5]

                    if key == 'all_receivers':
                        all_receivers_list = play_dic[key]
                        while len(all_receivers_list) > 5:
                            block_bubble_indices = [i for i in range(len(all_receivers_list)) if all_receivers_list[i][1] == 'bubble.block']
                            if block_bubble_indices:
                                remove_ind = np.random.choice(block_bubble_indices, 1)[0]
                            else:
                                remove_ind = np.random.choice(list(range(len(all_receivers_list))), 1)[0]
                            del all_receivers_list[remove_ind]

                        position_names = ['L2', 'L1', 'M', 'R1', 'R2']

                        distances = [tup[2] for tup in all_receivers_list]

                        for pn_ind, pn in enumerate(position_names):
                            new_row['player{}_position'.format(pn)] = all_receivers_list[pn_ind][0]
                            new_row['player{}_route'.format(pn)] = all_receivers_list[pn_ind][1]
                            new_row['player{}_depth'.format(pn)] = all_receivers_list[pn_ind][3]
                            new_row['player{}_distance'.format(pn)] = distances[pn_ind]
                            new_row['player{}_intended'.format(pn)] = all_receivers_list[pn_ind][4]
                            new_row['player{}_height'.format(pn)] = all_receivers_list[pn_ind][5]
                            new_row['player{}_weight'.format(pn)] = all_receivers_list[pn_ind][6]
                            new_row['player{}_top_speed'.format(pn)] = all_receivers_list[pn_ind][7]

                        #wanted to organize the receivers by where they were on the field
                        #8 positions L3, L2, L1, M1, M2, R1, R2, R3
                        #These positions would be defined by ranges across the line of scrimmage
                        #this code helps rearrange them according to a one rule if there were receivers in the same ranges
                        #   -if more than one receiver occupied a space move the one closer to the center inwards
                        # all_receivers_list = play_dic[key]
                        # all_receivers_list.sort(key= lambda t: t[2], reverse=True)
                        # distances = [tup[2] for tup in all_receivers_list]
                        #
                        # if any([np.isnan(d) for d in distances]):
                        #     new_row = {}
                        #     break
                        #
                        # position_names = ['L3', 'L2', 'L1', 'M1', 'M2', 'R1', 'R2', 'R3']
                        # position_fill = [-1] * len(position_names)
                        # position_distances = [(19, 100), (12, 19), (5, 12), (0, 5), (-5, 0), (-12, -5), (-19, -12),
                        #                       (-100, -19)]
                        #
                        # position_ind = 0
                        # pos_dist_ind = 0
                        # d_ind = 0
                        # flag = False
                        # while d_ind < len(distances):
                        #     d = distances[d_ind]
                        #     if position_distances[pos_dist_ind][0] <= d < position_distances[pos_dist_ind][1]:
                        #
                        #         # if there is already a player there then need to do something otherwise can just skip and assign the player
                        #         original_position_ind = position_ind
                        #         while position_fill[position_ind] != -1:
                        #
                        #             # move towards the center if positive, won't be anyone there because of how for loop steps through
                        #             if np.average(position_distances[pos_dist_ind]) > 0:
                        #                 position_ind += 1
                        #                 pos_dist_ind += 1
                        #             else:
                        #                 # start with position to the immediate left
                        #                 position_ind_copy = position_ind - 1
                        #
                        #                 # mark how many players need to be moved left
                        #                 while position_fill[position_ind_copy] != -1:
                        #                     position_ind_copy -= 1
                        #
                        #                     # if we've gone all the way to left then just need to move to the right instead
                        #                     if position_ind_copy == 0:
                        #                         flag = True
                        #                         break
                        #                 if flag:
                        #                     position_ind += 1
                        #                     continue
                        #
                        #                 # re arrange going left to right leaving original position_ind to be -1
                        #                 for i in range(position_ind_copy, position_ind):
                        #                     position_fill[i] = position_fill[i + 1]
                        #                     position_fill[i + 1] = -1
                        #
                        #         # fill in the position
                        #         position_fill[position_ind] = d_ind
                        #         d_ind += 1
                        #         position_ind, pos_dist_ind = original_position_ind, original_position_ind
                        #
                        #         # note that position_ind does not ever move left, only right
                        #
                        #     # if no player at the position move position_ind to the right
                        #     else:
                        #         pos_dist_ind += 1
                        #         position_ind += 1
                        #
                        # for pn_ind, pn in enumerate(position_names):
                        #     #if no receiver in the position label it as missing otherwise fill in the proper information
                        #     if position_fill[pn_ind] == -1:
                        #         new_row['player{}_position'.format(pn)] = None
                        #         new_row['player{}_route'.format(pn)] = None
                        #         new_row['player{}_depth'.format(pn)] = None
                        #         new_row['player{}_distance'.format(pn)] = None
                        #         new_row['player{}_intended'.format(pn)] = None
                        #         continue
                        #
                        #     new_row['player{}_position'.format(pn)] = all_receivers_list[position_fill[pn_ind]][0]
                        #     new_row['player{}_route'.format(pn)] = all_receivers_list[position_fill[pn_ind]][1]
                        #     new_row['player{}_depth'.format(pn)] = all_receivers_list[position_fill[pn_ind]][3]
                        #     new_row['player{}_distance'.format(pn)] = distances[position_fill[pn_ind]]
                        #     new_row['player{}_intended'.format(pn)] = all_receivers_list[position_fill[pn_ind]][4]

                if bool(new_row):
                    #if the play was a sack skip it here
                    game_rows.append(new_row)
                else:
                    continue

            out_filename = r'C:\Users\Mitch\Documents\UofM\Fall 2018\NFL\Data\all_plays\all_plays_{}.csv'.format(str(gameId))
            out_df = pd.DataFrame(game_rows)
            out_df.to_csv(out_filename, index=False)

            print('All done with pass {} which is gameId {}'.format(num_doc, str(gameId)))













