from shapely.geometry import LineString

def s2LS(s):
    '''
    takes in a string that has been saved in a .csv and turns it back into a LineString object
    :param s:
    :return: a LineString object
    '''

    #get only values within parenthesis
    indL = s.index('(')
    indR = s.index(')')
    onlyNumbers = s[(indL + 1):indR]
    pairs = onlyNumbers.split(',')
    split_pairs = [p.lstrip().split(' ') for p in pairs]
    split_pairs_float = [list(map(float, sublist)) for sublist in split_pairs]
    ls = LineString(split_pairs_float)
    return ls

