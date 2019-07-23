from shapely.geometry import LineString, Point, MultiPoint
from shapely import affinity
import pandas as pd
import os
from string_to_LineString import s2LS

def add_more_points(LS):
    #add more points between pre defined routes to help with neareast neighbor approach
    num_intra_points = 4

    Xs, Ys = LS.coords.xy

    new_Xs, new_Ys = [], []

    for i in range(len(Xs) - 1):
        X1, X2 = Xs[i], Xs[i + 1]
        Y1, Y2 = Ys[i], Ys[i + 1]

        delta_x = X2 - X1
        delta_y = Y2 - Y1

        spacing_x = delta_x/ (num_intra_points + 1)
        spacing_y = delta_y/ (num_intra_points + 1)

        new_Xs.append(X1)
        new_Xs += [X1 + spacing_x*j for j in range(1, num_intra_points + 1)]

        new_Ys.append(Y1)
        new_Ys += [Y1 + spacing_y*j for j in range(1, num_intra_points + 1)]

    new_Xs.append(Xs[-1])
    new_Ys.append(Ys[-1])

    return LineString([(x, y) for x, y in zip(new_Xs, new_Ys)])


#route trees started out complicated but when I got to analysis there was too many combinations so I decided to simply the tree
#note that these are scaled up quite large so that when they are resized they are guaranteed to be shrunk
route_tree_dic = {
    'flat': LineString([(0, 0), (-100, 70)]),
    'slant': LineString([(0, 0), (100, 70)]),
    'curl': LineString([(0, 0), (0, 100), (30, 70)]),
    'comeback': LineString([(0, 0), (0, 100), (-30, 70)]),
    'out': LineString([(0, 0), (0, 100), (-40, 100)]),
    'dig': LineString([(0, 0), (0, 100), (70, 100)]),
    'post': LineString([(0, 0), (0, 100), (60, 200)]),
    'corner': LineString([(0, 0), (0, 100), (-40, 200)]),
    'streak': LineString([(0, 0), (0, 200)]),
    # 'slugo': LineString([(0, 0), (10, 7), (10, 20)]),
    'wheel': LineString([(0, 0), (-100, 70), (-100, 200)]),
    # 'in.and.out': LineString([(0, 0), (10, 7), (10, 12), (0, 20)]),
    # 'out.and.in': LineString([(0, 0), (-10, 7), (-10, 12), (0, 20)]),
    # 'stick.in': LineString([(0, 0), (10, 5), (10, 3)]),
    # 'stick.out': LineString([(0, 0), (-10, 5), (-10, 3)]),
}

for key, value in route_tree_dic.items():
    route_tree_dic[key] = add_more_points(value)

rb_route_tree_dic = {
    'flat': LineString([(0, 0), (-100, 70)]),
    'slant': LineString([(0, 0), (100, 70)]),
    'curl': LineString([(0, 0), (0, 100), (30, 70)]),
    'comeback': LineString([(0, 0), (0, 100), (-30, 70)]),
    'out': LineString([(0, 0), (0, 100), (-40, 100)]),
    'dig': LineString([(0, 0), (0, 100), (70, 100)]),
    'post': LineString([(0, 0), (0, 100), (60, 200)]),
    'corner': LineString([(0, 0), (0, 100), (-40, 200)]),
    # 'in.and.out': LineString([(0, 0), (10, 7), (10, 12), (0, 20)]),
    # 'out.and.in': LineString([(0, 0), (-10, 7), (-10, 12), (0, 20)]),
    # 'stick.in': LineString([(0, 0), (10, 5), (10, 3)]),
    # 'stick.out': LineString([(0, 0), (-10, 5), (-10, 3)]),
    # 'swing': LineString([(0, 0), (-10, 0), (-13, 3), (-13, 10)]),
    # 'v.out': LineString([(0, 0), (-5, 5), (0, 8)]),
    # 'v.in': LineString([(0, 0), (5, 5), (0, 8)]),
    # 'slugo': LineString([(0, 0), (10, 7), (10, 20)]),
    'wheel': LineString([(0, 0), (-100, 70), (-100, 200)]),
}

for key, value in rb_route_tree_dic.items():
    rb_route_tree_dic[key] = add_more_points(value)

def df_to_ls(df):
    '''
    dataframe to dataframe of Linestrings
    :param: wr route only dataframe
    :return: save and return dataframe of linestrings
    '''
    num_rows = df.shape[0]
    current_row = 1
    out_df = pd.DataFrame()
    while current_row < num_rows:
        nflId = df['nflId'].iloc[current_row]
        playId = df['playId'].iloc[current_row]
        position = df['position'].iloc[current_row]
        running_right = df['running_right'].iloc[current_row]
        current_points = []
        while (current_row < num_rows) and (df['nflId'].iloc[current_row] == nflId) and (df['playId'].iloc[current_row] == playId):
            current_points.append((df['x'].iloc[current_row], df['y'].iloc[current_row]))
            current_row += 1

        ls = LineString(current_points)
        out_df = out_df.append(pd.DataFrame({'nflId': [int(nflId)], 'playId': [int(playId)], 'LineString': [ls], 'position': [position], 'running_right': [running_right]}))

    return out_df



def rotate_route(ls, running_right):
    #I wanted all routes to be seen as if they were running up from the left side of the ball
    #need to find where they were orignially run from and transform accordingly
    coords = list(ls.coords)
    centroidx, _ = ls.centroid.xy
    starting_point = coords[0]
    running_right = bool(running_right)
    high_field = starting_point[1] > 53.3/2

    #translate the points so the starting point is at 0, 0
    moved_coords = LineString([(c[0] - starting_point[0], c[1] - starting_point[1]) for c in coords])

    #error checking
    #draw_routes(ls.coords)

    #routes are designed to be from 'left' side of the ball
    if high_field:
        if running_right:
            #do a 90 degree turn
            transformed_coords = affinity.rotate(moved_coords, angle=90, origin=(0, 0))
        else:
            #do a 270 degree turn and mirror
            transformed_coords = affinity.affine_transform(affinity.rotate(moved_coords, angle=270, origin=(0, 0)), (-1, 0, 0, 1, 0, 0))
    else:
        if running_right:
            #do a 90 degree turn and mirror
            transformed_coords = affinity.affine_transform(affinity.rotate(moved_coords, angle=90, origin=(0, 0)), (-1, 0, 0, 1, 0, 0))
            pass
        else:
            #do a 270 degree turn
            transformed_coords = affinity.rotate(moved_coords, angle=270, origin=(0, 0))

    #error checkign
    #draw_routes(transformed_coords)

    return transformed_coords

import matplotlib.pyplot as plt

#was nice for error checking to see if the functions were doing what I wanted
def draw_routes(player_route):
    fig = plt.figure()
    ax = fig.add_subplot(111, xlim=(-10, 130), ylim=(-5, 58.3))
    # field = plt.Rectangle(xy=(0, 0), width=120, height=53.3, fill=False, lw=2)
    # hash_lines = plt.Rectangle(xy=(10, 23.36667), width=100, height=6.6, fill=False, lw=1, linestyle='--')
    # goal_lines = plt.Rectangle(xy=(10, 0), width=100, height=53.3, fill=False, lw=1)
    # ax.add_patch(field)
    # ax.add_patch(hash_lines)
    # ax.add_patch(goal_lines)

    player_x, player_y = player_route.xy
    ax.plot(player_x, player_y, 'b-')
    fig.show()
    plt.close(fig)

def resize_route_wrapper_wrte(wr_route):
    #wrapper to resize all pre defined routes according to wrte route
    transformed_wr_tree = {}
    for name, dic_route in route_tree_dic.items():
        transformed_wr_tree[name] = resize_route(wr_route, dic_route)
        #draw_routes(wr_route, transformed_wr_tree[name], name)
    return transformed_wr_tree

def resize_route_wrapper_rb(rb_route):
    # wrapper to resize all pre defined routes according to rb route
    transformed_rb_tree = {}
    for name, dic_route in rb_route_tree_dic.items():
        transformed_rb_tree[name] = resize_route(rb_route, dic_route)
        #draw_routes(rb_route, transformed_rb_tree[name], name)
    return transformed_rb_tree

def resize_route(wr_route, dic_route):
    #resize the pre defined route based on the largest x or y discrepency and by ratio to maintain proper shape

    #get (minx, maxx, miny, maxy)
    bounding_box_wr = wr_route.bounds
    bounding_box_dic = dic_route.bounds

    mp_dic = MultiPoint(dic_route.coords)

    new_points = []

    #make sure they have the same minx, miny coordinate
    min_x_delta = bounding_box_wr[0] - bounding_box_dic[0]
    min_y_delta = bounding_box_wr[1] - bounding_box_dic[1]

    maxabs_x_delta = abs(bounding_box_wr[2] - (bounding_box_dic[2] + min_x_delta))
    maxabs_y_delta = abs(bounding_box_wr[3] - (bounding_box_dic[3] + min_y_delta))

    #change ratio depending on biggest discrepancy
    if maxabs_x_delta <= maxabs_y_delta:
        ratio = (bounding_box_wr[3] - bounding_box_wr[1])/(bounding_box_dic[3] - bounding_box_dic[1])

    else:
        ratio = (bounding_box_wr[2] - bounding_box_wr[0])/(bounding_box_dic[2] - bounding_box_dic[0]) if (bounding_box_dic[2] - bounding_box_dic[0]) != 0 else (bounding_box_wr[3] - bounding_box_wr[1])/(bounding_box_dic[3] - bounding_box_dic[1])

    #assign new points using projected coordinate position along its own boundary, ratio and new minimum
    for p in mp_dic:
        x, y = p.x, p.y

        new_x = ratio*(x - bounding_box_dic[0]) + bounding_box_wr[0]
        new_y = ratio*(y - bounding_box_dic[1]) + bounding_box_wr[1]

        new_points.append((new_x, new_y))

    # #trying to see if starting the routes at the same point helps ~sorta
    # first_wr_point = Point(new_points[0])
    #
    # diffx, diffy = first_wr_point.x, first_wr_point.y
    #
    # for i, np in enumerate(new_points):
    #     new_points[i] = (np[0] - diffx, np[1] - diffy)

    return LineString(new_points)


def find_closest_route_wrte(wr_route):
    #based on nearest neighbor approach find distance between all points of wr te route and resized pre defined routes
    #whichever pre defined route had the smallest total distance, assign the wr te route that name
    resized_route_dic = resize_route_wrapper_wrte(wr_route)
    route_distance_list = []
    wr_mp = MultiPoint(wr_route.coords)

    #if the end point is less than 4 yards away or the player ends the play behind the line of scrimmage assign them to be
    #blocking or attempting to catch a screen/ bubble route
    if wr_mp[-1].distance(Point((0,0))) <= 4 or wr_mp[-1].y <= 0:
        return 'bubble.block'
    for route_name, route_ls in resized_route_dic.items():
        total_distance_wr = 0
        for P1 in wr_mp:
            total_distance_wr += P1.distance(route_ls)
        average_distance_wr = total_distance_wr/len(wr_mp)
        total_distance_dic = 0
        dic_mp = MultiPoint(route_ls.coords)
        for P2 in dic_mp:
            total_distance_dic += P2.distance(wr_route)
        average_distance_dic = total_distance_dic/len(dic_mp)
        route_distance_list.append((route_name, average_distance_wr + average_distance_dic))
    closest_route_str = min(route_distance_list, key= lambda t: t[1])[0]

    return closest_route_str

def find_closest_route_rb(rb_route):
    resized_route_dic = resize_route_wrapper_rb(rb_route)
    route_distance_list = []
    rb_mp = MultiPoint(rb_route.coords)
    if rb_mp[-1].distance(Point((0,0))) <= 4 or rb_mp[-1].y <= 0:
        return 'bubble.block'
    for route_name, route_ls in resized_route_dic.items():
        total_distance_rb = 0
        for P1 in rb_mp:
            total_distance_rb += P1.distance(route_ls)
        average_distance_rb = total_distance_rb/len(rb_mp)
        total_distance_dic = 0
        dic_mp = MultiPoint(route_ls.coords)
        for P2 in dic_mp:
            total_distance_dic += P2.distance(rb_route)
        average_distance_dic = total_distance_dic/len(dic_mp)
        route_distance_list.append((route_name, average_distance_rb + average_distance_dic))
    closest_route_str = min(route_distance_list, key= lambda t: t[1])[0]

    return closest_route_str

def assign_routes_wrte(wrte_df):
    #wrapper function to assign all routes and use all functions above
    filename = r'C:\Users\Mitch\Documents\UofM\Fall 2018\NFL\Data\wrteRoutes_{}.csv'.format(int(wrte_df['gameId'].iloc[0]))
    if os.path.exists(filename):
        return pd.read_csv(filename)
    ls_df = df_to_ls(wrte_df)
    ls_df['route'] = ls_df.apply(lambda row: find_closest_route_wrte(rotate_route(row['LineString'], row['running_right'])), axis=1)

    ls_df.drop(columns=['running_right'], inplace=True)
    #ls_df.to_csv(filename)

    return ls_df

def assign_routes_rb(rb_df):
    filename = r'C:\Users\Mitch\Documents\UofM\Fall 2018\NFL\Data\rbRoutes_{}.csv'.format(int(rb_df['gameId'].iloc[0]))
    if os.path.exists(filename):
        return pd.read_csv(filename)
    ls_df = df_to_ls(rb_df)
    ls_df['route'] = ls_df.apply(lambda row: find_closest_route_rb(rotate_route(row['LineString'], row['running_right'])), axis=1)

    ls_df.drop(columns=['running_right'], inplace=True)
    #ls_df.to_csv(filename)

    return ls_df

def game_play_id_routes(gameId, playId, running_right, routeInd=None):
    routeID = [routeInd] if routeInd is not None else list(range(5))

    route_df = pd.read_csv(r'C:\Users\Mitch\Documents\UofM\Fall 2018\NFL\Data\Routes_{}.csv'.format(gameId))

    play_df = route_df.loc[route_df['playId'] == playId]

    route_line = play_df.iloc[routeID]

    if isinstance(route_line, pd.DataFrame):
        for i in range(route_line.shape[0]):
            ls = rotate_route(s2LS(route_line['LineString'].iloc[i]), running_right=running_right)
            labeled_route_name = route_line['route'].iloc[i]
            tree = resize_route_wrapper_wrte(ls)
            new_label = find_closest_route_wrte(ls)
            print('')

            for name, route_ls in tree.items():
                fig = plt.figure()
                plt.title('Proposed route: {}, Labeled route: {}'.format(name, labeled_route_name))

                assert isinstance(ls, LineString)
                assert isinstance(route_ls, LineString)

                ax = fig.add_subplot(111, xlim=(min(ls.bounds[0], route_ls.bounds[0]) - 5, max(ls.bounds[2], route_ls.bounds[2]) + 5),
                                     ylim=(min(ls.bounds[1], route_ls.bounds[1]) - 5, max(ls.bounds[3], route_ls.bounds[3]) + 5))

                player_x, player_y = ls.xy
                dic_x, dic_y = route_ls.xy
                ax.plot(player_x, player_y, 'b-')
                ax.plot(dic_x, dic_y, 'r-')
                fig.show()
                plt.close(fig)

    return -1


if __name__ == '__main__':
    #if this file is run as main create Route.csv files for each game
    # from load_data import Data
    #
    # header = r'C:\Users\Mitch\Documents\UofM\Fall 2018\NFL\Data'
    # tracking_files = [f for f in os.listdir(header) if f[0] == 't']
    # for num_doc, tf in enumerate(tracking_files):
    #     gameId = tf[16:-4]
    #     filename = r'C:\Users\Mitch\Documents\UofM\Fall 2018\NFL\Data\Routes_{}.csv'.format(gameId)
    #
    #     if os.path.exists(filename):
    #         continue
    #
    #     data = Data()
    #     data.load_tracking_info(str(gameId))
    #     data.load_wr_te_routes()
    #     data.load_rb_routes()
    #     wrte = data.get_wr_te_routes()
    #     rb = data.get_rb_routes()
    #     df_wrte = wrte[str(gameId)]
    #     df_rb = rb[str(gameId)]
    #     wrte_route_df = assign_routes_wrte(df_wrte)
    #     rb_route_df = assign_routes_rb(df_rb)
    #
    #     combined_df = wrte_route_df.append(rb_route_df)
    #     combined_df = combined_df.sort_values(by='playId', ascending=True)
    #     combined_df.to_csv(filename, index=False)
    #
    #     del data
    #
    #     print('All done with pass {} which is gameId {}'.format(num_doc, str(gameId)))

    game_play_id_routes(2017100108, 4098, False)