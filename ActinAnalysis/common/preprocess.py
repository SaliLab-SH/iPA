# coding = utf-8


'''
preprocess of the volumn and .mrc from raw data

By angdi  and bing

'''


from matplotlib.cbook import maxdict
import numpy as np
import scipy
import scipy.ndimage
import skimage.feature
import skimage.morphology
import pandas as pd
import matplotlib.pyplot as plt
import mrcfile, tifffile
import xml
import xml.etree.ElementTree as ET
import json
import os, sys
import copy
from common.parser import arg
from scipy.spatial.distance import cdist

def read_dataid():

    idlst = pd.read_csv(os.path.join(arg.root_dir, arg.parameter_dir, 'dataid.csv'))

    idlst = [ idlst['ID'].iloc[i] for i in range(0,idlst.shape[0])  ]

    return idlst


def obtain_dir():
    for maindir, subdir, file_name_list in os.walk(os.path.join(arg.root_dir, arg.data_root_dir), topdown=False):
        locate = np.array(subdir)

    subdir1 = [  os.path.join(arg.root_dir, arg.data_root_dir, subdir) for subdir in locate]
    subdir2 = []
    for dir in subdir1:
        for maindir, subdir, file_name_list in os.walk(dir, topdown=False):
            condition = np.array(subdir)
        for singlecondition in condition:
            subdir2.append( os.path.join(dir, singlecondition) )
    subdir3 = []
    for dir in subdir2:
        for maindir, subdir, file_name_list in os.walk(dir, topdown=False):
            replica = np.array(subdir)
        for singlereplica in replica:
            subdir3.append( os.path.join(dir, singlereplica) )

    return subdir2, subdir3 



def load_files(data_dir, check_ = False):

    mrcname = ''
    mrcname_fullfill = ''
    filamentdata_xmlname = ''
    filamentdata_jsonname_full = ''
    mtdata_xmlname = ''
    MT_jsonname_full = ''
    mito_mrcname = ''
    er_mrcname = ''
    for maindir, subdir, file_name_list in os.walk(data_dir, topdown=False):
        filelist = np.array(file_name_list)

    for name in filelist:
        if '.mrc' in name and 'mito' not in name.lower() and 'er' not in name.lower():
            if 'fullfill' not in name:
                mrcname = os.path.join(data_dir, name)   
            elif 'fullfill' in name:
                mrcname_fullfill = os.path.join(data_dir, name)      
        elif '.xml' in name and ('filament' in name.lower() or 'actin' in name.lower()):
            if 'fill' not in name:
                filamentdata_xmlname = os.path.join(data_dir, name)
        elif '.json' in name and 'filament' in name and 'filled' in name:
                filamentdata_jsonname_full = os.path.join(data_dir, name)
        
        elif '.xml' in name and ('mt' in name.lower() or 'microtube' in name.lower()):
            if 'fill' not in name:
                mtdata_xmlname = os.path.join(data_dir, name)
        elif '.json' in name and 'mt' in name.lower() and 'filled' in name:
                MT_jsonname_full = os.path.join(data_dir, name)
        elif '.mrc' in name and 'mito_filled' in name.lower():
                mito_mrcname = os.path.join(data_dir, name)
        elif '.mrc' in name and 'er' in name.lower():
                er_mrcname = os.path.join(data_dir, name)
        else:
            pass

    if check_:
        if os.path.exists(mrcname): print('isg mrcname:',mrcname) 
        if os.path.exists(mrcname_fullfill): print('isg mrcname processed:', mrcname_fullfill)
        if os.path.exists(filamentdata_xmlname): print('actin xmlname:', filamentdata_xmlname)
        if os.path.exists(filamentdata_jsonname_full): print('actin xml processed:', filamentdata_jsonname_full)
        if os.path.exists(mtdata_xmlname): print('MT name', mtdata_xmlname)
        if os.path.exists(MT_jsonname_full): print('MT processed:', MT_jsonname_full)


    filenamelst = [mrcname, mrcname_fullfill, filamentdata_xmlname, 
                    filamentdata_jsonname_full, 
                    mtdata_xmlname, MT_jsonname_full]

    return filenamelst



def get_filelist(dir, Filelist):
    newDir = dir
    if os.path.isfile(dir):
        Filelist.append(dir)

    elif os.path.isdir(dir):
        for s in os.listdir(dir):
            newDir=os.path.join(dir,s)
            get_filelist(newDir, Filelist)
    return Filelist


def load_parameters(file_idx):
    parpath_ = os.path.join(arg.root_dir, arg.parameter_dir)
    df = pd.read_excel(f'{parpath_}/Mrc_offset.xlsx')
    shift_data=df.values

    shift_xml_xyz = [0,0,0]
    testn = 0
    for line in shift_data:
        if file_idx in line:
            testn = 1
            shift_xml_xyz[0] = line[3]
            shift_xml_xyz[1] = line[2]
            shift_xml_xyz[2] = line[1]

    assert testn == 1

    df2 = pd.read_excel(f'{parpath_}/Voxel_Size.xlsx')
    voxel_data=df2.values
    voxel_size_xyz = [0,0,0]
    for line in voxel_data:
        if file_idx in line:
            voxel_size_xyz[0] = line[1]
            voxel_size_xyz[1] = line[2]
            voxel_size_xyz[2] = line[3]
    return shift_xml_xyz, voxel_size_xyz



def fill_actin_point_gap(filename, resolution = 40):
    '''
    fill pints gap to generate representative points for filaments
    reorganize coordinates zyx to xyz
    
    '''

    xmlname = f'{filename}'
    f = open(xmlname, encoding='utf-8') 
    xml_txt = f.read()
    root = ET.fromstring(xml_txt)

    pagenamelst = ['nodes', 'points', 'segments']
    pageidxlst = []
    for ii, pagename in enumerate(pagenamelst):
        for jj, child in enumerate(root):
            if pagename in str(child.attrib).lower():
                pageidxlst.append(jj)

    assert len(pageidxlst) == 3

    Nodes = root[int(pageidxlst[0])]
    Points = root[int(pageidxlst[1])]
    Segments = root[int(pageidxlst[2])]

    assert 'table' in str(list(Points)[0].tag).lower()
    tablee = list(Points)[0]

    rowws = list(tablee) #tablee.getchildren()

    start_n = 0
    for row_i, roww in enumerate(rowws):
        if 'row' not in str(roww.tag).lower():
            start_n = row_i+1
        else: break


    columnname = []
    for roww in rowws[:start_n + 1]:
        if 'row' in str(roww.tag).lower():
            cells = list(roww)  #roww.getchildren()

            for cell in cells:
                datas = list(cell)  #cell.getchildren()
                for data in datas:
                    columnname.append(data.text)


    Points_pd =pd.DataFrame(columns=columnname)

    # obtain data
    for i, roww in enumerate(rowws[start_n + 1:]):
        cells = roww.getchildren()
        templst = []
        for cell in cells:
            datas = cell.getchildren()
            for data in datas:
                templst.append(data.text)
        Points_pd.loc[i]= templst

    Points_pd_x = Points_pd['X Coord']
    Points_pd_y = Points_pd['Y Coord']
    Points_pd_z = Points_pd['Z Coord']

    Points_pd_x = pd.to_numeric(Points_pd_x) * 10
    Points_pd_y = pd.to_numeric(Points_pd_y) * 10
    Points_pd_z = pd.to_numeric(Points_pd_z) * 10 

    assert 'table' in str(list(Points)[0].tag).lower()

    tablee = Segments.getchildren()[0]
    rowws = tablee.getchildren()


    start_n = 0
    for row_i, roww in enumerate(rowws):
        if 'row' not in str(roww.tag).lower():
            start_n = row_i+1
        else: break


    for roww in rowws[:start_n+1]:

        if 'row' in str(roww.tag).lower():
            cells = roww.getchildren()
            columnname = []
            for cell in cells:
                datas = cell.getchildren()
                for data in datas:
                    columnname.append(data.text)

    segments_pd =pd.DataFrame(columns=columnname)


    # obtain data
    for i, roww in enumerate(rowws[start_n+1:]):
        cells = roww.getchildren()
        templst = []
        for cell in cells:
            datas = cell.getchildren()
            for data in datas:
                templst.append(data.text)
        segments_pd.loc[i]= templst

    points = Points_pd

    coords_lst = []
    for idx in range(len(points)):
        coords_lst.append([float(Points_pd_z.iloc[idx]), float(Points_pd_y.iloc[idx]), float(Points_pd_x.iloc[idx])])  #zyx from csv

    coords = np.array(coords_lst).reshape(-1,3)



    #  
    # input relations 
    filament_idx_lst = []
    filament_coord_lst = []
    segments = segments_pd

    correspond_id = segments['Segment ID']
    correspond_points_id = segments['Point IDs']


    for idx in range(len(segments)):
    # for idx in range(5):
        assert idx ==  int(correspond_id[idx]) ## match 
        coord_idx_lst = [ int(pointidx) for pointidx in correspond_points_id[idx].split(',') ]# str
        temp_coordlst = []
        for pointidx in coord_idx_lst:
            temp_coordlst.append(coords[pointidx])  #   coords[pointidx] unit A  , voxel_coords

        filament_idx_lst.append(coord_idx_lst)
        filament_coord_lst.append(temp_coordlst)




    # fill gap

    
    def get_line(point1, point2):

        vect = np.array([point2[0] - point1[0], point2[1] - point1[1], point2[2] - point1[2]])
        unit_vect = vect / np.linalg.norm(vect) * resolution  # 4 nm
        thresh_len = resolution * 0.3 
        points_sum = np.linalg.norm(vect) // resolution + 1
        pointslst = [ point1 + unit_vect * i for i in range(int(points_sum))]
        while True:
            if np.linalg.norm(pointslst[-1] - point1) >= np.linalg.norm(vect):
                pointslst = pointslst[:-1]
            else:
                break
        if np.linalg.norm(pointslst[-1] - np.array(point2)) < thresh_len:
            pointslst = pointslst[:-1]
        return pointslst


    # extend coord to narrow interval
    filament_coord_extend_lst = [] 

    for filament_coord in filament_coord_lst:
        single_filament_coord_extend_lst = []
        for i in range(1, len(filament_coord)):
            coord1, coord2 = filament_coord[i-1], filament_coord[i]  # array
            temp_coord_extend_lst = get_line(coord1, coord2)

            single_filament_coord_extend_lst.extend(temp_coord_extend_lst)
        single_filament_coord_extend_lst.append(filament_coord[-1])

        filament_coord_extend_lst.append(single_filament_coord_extend_lst)

    # filament_coord_lst
    coords_intotal1 = 0
    for single_filament_coord_extend in filament_coord_lst:
        coords_intotal1 += len(single_filament_coord_extend)
    print(f'coord intotal before: {coords_intotal1}')
    coords_intotal = 0
    for single_filament_coord_extend in filament_coord_extend_lst:
        coords_intotal += len(single_filament_coord_extend)
    print(f'coord intotal after: {coords_intotal}')


    ## save file
    img_columnname = [f'filament_{i}' for i in range(len(filament_coord_extend_lst))]
    points_all_dict = {}

    for idx, filament_vects in enumerate(filament_coord_extend_lst):
        points_all_dict[f'{img_columnname[idx]}'] =[]
        for coord in filament_vects:
            points_all_dict[f'{img_columnname[idx]}'] .append(list(coord))

    return points_all_dict




def get_actin_to_pm_angle(filename):
    '''
    obtain actin angle to plasma membrane from segmentation page
    '''

    xmlname = f'{filename}'
    f = open(xmlname, encoding='utf-8') 
    xml_txt = f.read()
    root = ET.fromstring(xml_txt)


    pagenamelst = ['nodes', 'points', 'segments']
    pageidxlst = []
    for ii, pagename in enumerate(pagenamelst):
        for jj, child in enumerate(root):

            if pagename in str(child.attrib).lower():
                pageidxlst.append(jj)

    assert len(pageidxlst) == 3

    Nodes = root[int(pageidxlst[0])]
    Points = root[int(pageidxlst[1])]
    Segments = root[int(pageidxlst[2])]


    assert 'table' in str(list(Points)[0].tag).lower()

    tablee = Segments.getchildren()[0]
    rowws = tablee.getchildren()

    start_n = 0
    for row_i, roww in enumerate(rowws):
        if 'row' not in str(roww.tag).lower():
            start_n = row_i+1
        else: break


    for roww in rowws[:start_n+1]:
        if 'row' in str(roww.tag).lower():
            cells = roww.getchildren()
            columnname = []
            for cell in cells:
                datas = cell.getchildren()
                for data in datas:
                    columnname.append(data.text)

    segments_pd =pd.DataFrame(columns=columnname)

    # obtain data
    for i, roww in enumerate(rowws[start_n+1:]):
        cells = roww.getchildren()
        templst = []
        for cell in cells:
            datas = cell.getchildren()
            for data in datas:
                templst.append(data.text)
        segments_pd.loc[i]= templst

    # input relations 
    filament_idx_lst = []
    filament_coord_lst = []
    segments = segments_pd

    correspond_id = segments['Segment ID']
    correspond_points_id = segments['Point IDs']
    correspond_angle = segments['OrientationTheta']


    ## save file
    img_columnname = [f'filament_{i}' for i in range(correspond_angle.shape[0])]


    angle_all_dict = {}

    for idx, name in enumerate(img_columnname):
        angle_all_dict[f'{name}'] = correspond_angle.iloc[idx]


    return angle_all_dict



def fill_microtube_point_gap(filename,  resolution = 40):
    '''
    fill pints gap to generate representative points for filaments
    reorganize coordinates zyx to xyz
    
    '''
    xmlname = filename
    f = open(xmlname, encoding='utf-8') 
    xml_txt = f.read()
    root = ET.fromstring(xml_txt)


    pagenamelst = ['nodes', 'points', 'segments']
    pageidxlst = []
    for ii, pagename in enumerate(pagenamelst):
        for jj, child in enumerate(root):
            if pagename in str(child.attrib).lower():
                pageidxlst.append(jj)
    assert len(pageidxlst) == 3

    Nodes = root[int(pageidxlst[0])]
    Points = root[int(pageidxlst[1])]
    Segments = root[int(pageidxlst[2])]


    assert 'table' in str(list(Points)[0].tag).lower()
    tablee = list(Points)[0]

    rowws = list(tablee) #tablee.getchildren()


    start_n = 0
    for row_i, roww in enumerate(rowws):
        if 'row' not in str(roww.tag).lower():
            start_n = row_i+1
        else: break

    columnname = []
    for roww in rowws[:start_n + 1]:
        if 'row' in str(roww.tag).lower():
            cells = list(roww)  #roww.getchildren()

            for cell in cells:
                datas = list(cell)  #cell.getchildren()
                for data in datas:
                    columnname.append(data.text)

    Points_pd =pd.DataFrame(columns=columnname)



    # obtain data
    for i, roww in enumerate(rowws[start_n + 1:]):
        cells = roww.getchildren()
        templst = []
        for cell in cells:
            datas = cell.getchildren()
            for data in datas:
                templst.append(data.text)
        Points_pd.loc[i]= templst

    Points_pd_x = Points_pd['X Coord']
    Points_pd_y = Points_pd['Y Coord']
    Points_pd_z = Points_pd['Z Coord']

    Points_pd_x = pd.to_numeric(Points_pd_x) * 10
    Points_pd_y = pd.to_numeric(Points_pd_y) * 10
    Points_pd_z = pd.to_numeric(Points_pd_z) * 10

    # Segments = root[3]
    assert 'table' in str(list(Points)[0].tag).lower()

    tablee = Segments.getchildren()[0]
    rowws = tablee.getchildren()


    start_n = 0
    for row_i, roww in enumerate(rowws):
        if 'row' not in str(roww.tag).lower():
            start_n = row_i+1
        else: break


    for roww in rowws[:start_n+1]:
        if 'row' in str(roww.tag).lower():
            cells = roww.getchildren()
            columnname = []
            for cell in cells:
                datas = cell.getchildren()
                for data in datas:
                    columnname.append(data.text)

    segments_pd =pd.DataFrame(columns=columnname)
    # obtain data
    for i, roww in enumerate(rowws[start_n+1:]):
        cells = roww.getchildren()
        templst = []
        for cell in cells:
            datas = cell.getchildren()
            for data in datas:
                templst.append(data.text)
        segments_pd.loc[i]= templst



    points = Points_pd

    coords_lst = []
    for idx in range(len(points)):
        coords_lst.append([float(Points_pd_z.iloc[idx]), float(Points_pd_y.iloc[idx]), float(Points_pd_x.iloc[idx])])  #zyx from s

    coords = np.array(coords_lst).reshape(-1,3)

    # input relations 
    filament_idx_lst = []
    filament_coord_lst = []
    segments = segments_pd
    correspond_id = segments['Segment ID']
    correspond_points_id = segments['Point IDs']

    for idx in range(len(segments)):
        assert idx ==  int(correspond_id[idx]) ## match 
        coord_idx_lst = [ int(pointidx) for pointidx in correspond_points_id[idx].split(',') ]
        temp_coordlst = []
        for pointidx in coord_idx_lst:
            temp_coordlst.append(coords[pointidx]) 

        filament_idx_lst.append(coord_idx_lst)
        filament_coord_lst.append(temp_coordlst)

    # fill gap
    def get_line(point1, point2):
        vect = np.array([point2[0] - point1[0], point2[1] - point1[1], point2[2] - point1[2]])
        unit_vect = vect / np.linalg.norm(vect) * resolution 
        thresh_len = resolution * 0.3 
        points_sum = np.linalg.norm(vect) // resolution + 1
        pointslst = [ point1 + unit_vect * i for i in range(int(points_sum))]
        if np.linalg.norm(pointslst[-1] - point2) < thresh_len:
            pointslst = pointslst[:-1]
        return pointslst


    # extend coord to narrow interval
    filament_coord_extend_lst = [] 

    for filament_coord in filament_coord_lst:
        single_filament_coord_extend_lst = []
        for i in range(1, len(filament_coord)):
            coord1, coord2 = filament_coord[i-1], filament_coord[i]
            temp_coord_extend_lst = get_line(coord1, coord2)
            single_filament_coord_extend_lst.extend(temp_coord_extend_lst)
        single_filament_coord_extend_lst.append(filament_coord[-1])
        filament_coord_extend_lst.append(single_filament_coord_extend_lst)

    # filament_coord_lst
    coords_intotal1 = 0
    for single_filament_coord_extend in filament_coord_lst:
        coords_intotal1 += len(single_filament_coord_extend)
    coords_intotal = 0
    for single_filament_coord_extend in filament_coord_extend_lst:
        coords_intotal += len(single_filament_coord_extend)


    ## save file
    img_columnname = [f'MT_{i}' for i in range(len(filament_coord_extend_lst))]
    points_all_dict = {}

    for idx, filament_vects in enumerate(filament_coord_extend_lst):
        points_all_dict[f'{img_columnname[idx]}'] =[]
        for coord in filament_vects:
            points_all_dict[f'{img_columnname[idx]}'] .append(list(coord))

    
    return points_all_dict
