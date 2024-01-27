# coding = utf-8


'''
To calculate filament distance to vesicle membrane.
To calculate filaments points around individual surrounding spaces.

By angdi  and bing

'''

import numpy as np
import scipy.ndimage
import scipy
import math 
from scipy.spatial.distance import cdist
import copy 
import pandas as pd
from torch import threshold

def shift_bias(filament_coords, shift):
    # shift should be minused 
    filament_coords_modified ={}    
    filamentnames = list(filament_coords.keys())
    for idx in range(len(filamentnames)):
        curfilamentcoordslst = filament_coords[f'{filamentnames[idx]}']
        filament_coords_modified[f'{filamentnames[idx]}'] =[]

        for i, coord in enumerate(curfilamentcoordslst):
            coord_new = [coord[0] - shift[0], coord[1] - shift[1], coord[2] - shift[2]]
            # coord_new = [coord_new[0]/voxelsize[0], coord_new[1]/voxelsize[1], coord_new[2]/voxelsize[2]]
            filament_coords_modified[f'{filamentnames[idx]}'].append(coord_new)

    return filament_coords_modified



def coords_to_mem_distance_generator(filament_coords, mask, voxelsize):
    '''
    input: 
        filament coord: dict
        mask: array 
    '''

    filament_coords_modified ={}    
    filamentnames = list(filament_coords.keys())
    for idx in range(len(filamentnames)):
        curfilamentcoordslst = filament_coords[f'{filamentnames[idx]}']
        filament_coords_modified[f'{filamentnames[idx]}'] =[]
        for i, coord_new in enumerate(curfilamentcoordslst):
            coord_new = [coord_new[0]/voxelsize[0], coord_new[1]/voxelsize[1], coord_new[2]/voxelsize[2]]
            filament_coords_modified[f'{filamentnames[idx]}'].append(coord_new)


    coords_lst = []
    for idx in range(len(filamentnames)):
        curfilamentcoordslst = filament_coords_modified[f'{filamentnames[idx]}']
        coords_lst.extend(curfilamentcoordslst)

    coords = np.array(coords_lst).reshape(-1,3)
    # print('line56', coords)
    voxel_coords = coords

    coords_x = [ voxel_coords[i][0] for i in range(len(voxel_coords))]
    coords_y = [ voxel_coords[i][1] for i in range(len(voxel_coords))]
    coords_z = [ voxel_coords[i][2] for i in range(len(voxel_coords))]

    if len(coords_x) > 0:
        assert np.max(coords_x) <= mask.shape[0] and np.min(coords_x)>= 0
        assert np.max(coords_y) <= mask.shape[1] and np.min(coords_y)>= 0
        assert np.max(coords_z) <= mask.shape[2] and np.min(coords_z)>= 0
    else:
        print('no filament')

    img01 = mask
    reverse_img01 = img01 * (-1) + 1
    img01_rev_edt = scipy.ndimage.distance_transform_edt(reverse_img01)
    print(f'11111img01_rev_edt: {img01_rev_edt.shape}')
    avesize = np.average(voxelsize)
    dist = [ img01_rev_edt[int(voxel_coords[i][0])][int(voxel_coords[i][1])][int(voxel_coords[i][2])]  for i in range(len(voxel_coords)) ]
    dist = [dist[i] * avesize/10 for i in range(len(dist))] # nm
    print('111111distance num: ', len(dist))

    return dist




def coords_to_mem_distance_angle_pair_generator(filament_coords, filament_angle_data, mask, voxelsize):
    '''
    input: 
        filament coord: dict
        mask: array 
    actin_data_new, actin_angle_data, vesicle_mem_mask, voxel_size_xyz
    '''

    filament_coords_modified ={}    
    filament_angle_modified = {}
    filamentnames = list(filament_coords.keys())
    for idx in range(len(filamentnames)):
        curfilamentcoordslst = filament_coords[f'{filamentnames[idx]}']
        curfilamentangle = filament_angle_data[f'{filamentnames[idx]}']
        
        filament_coords_modified[f'{filamentnames[idx]}'] =[]
        filament_angle_modified[f'{filamentnames[idx]}'] =[]
        for i, coord_new in enumerate(curfilamentcoordslst):
            coord_new = [coord_new[0]/voxelsize[0], coord_new[1]/voxelsize[1], coord_new[2]/voxelsize[2]]
            filament_coords_modified[f'{filamentnames[idx]}'].append(coord_new)
            filament_angle_modified[f'{filamentnames[idx]}'].append(curfilamentangle)


    coords_lst = []
    angle_lst = []
    for idx in range(len(filamentnames)):
        curfilamentcoordslst = filament_coords_modified[f'{filamentnames[idx]}']
        curfilamentanglelst = filament_angle_modified[f'{filamentnames[idx]}']
        coords_lst.extend(curfilamentcoordslst)
        angle_lst.extend(curfilamentanglelst)

    assert len(coords_lst) == len(angle_lst)

    coords = np.array(coords_lst).reshape(-1,3)
    # print('line56', coords)
    voxel_coords = coords

    coords_x = [ voxel_coords[i][0] for i in range(len(voxel_coords))]
    coords_y = [ voxel_coords[i][1] for i in range(len(voxel_coords))]
    coords_z = [ voxel_coords[i][2] for i in range(len(voxel_coords))]



    if len(coords_x) > 0:
        assert np.max(coords_x) <= mask.shape[0] and np.min(coords_x)>= 0
        assert np.max(coords_y) <= mask.shape[1] and np.min(coords_y)>= 0
        assert np.max(coords_z) <= mask.shape[2] and np.min(coords_z)>= 0
    else:
        print('no filament')

    img01 = mask
    reverse_img01 = img01 * (-1) + 1
    img01_rev_edt = scipy.ndimage.distance_transform_edt(reverse_img01)


    avesize = np.average(voxelsize)
    dist = [ img01_rev_edt[int(voxel_coords[i][0])][int(voxel_coords[i][1])][int(voxel_coords[i][2])]  for i in range(len(voxel_coords)) ]
    
    angle = [ angle_lst[i] for i in range(len(voxel_coords))]
    
    dist = [dist[i] * avesize/10 for i in range(len(dist))] # nm
    # print('distance num', len(dist))

    return dist, angle



def pm_to_actin_angle_angle_pair_generator(filament_coords, filament_angle_data, mask_up, mask_down, voxelsize, threshold, threshold_actin):
    '''
    input: 
        filament coord: dict
        mask: array 
    actin_data_new, actin_angle_data, vesicle_mem_mask, voxel_size_xyz
    '''

    coords_lst = []
    angle_lst = []

    filamentnames = list(filament_coords.keys())

    for idx in range(len(filamentnames)):
        curfilamentcoordslst = filament_coords[f'{filamentnames[idx]}']
        curfilamentangle = filament_angle_data[f'{filamentnames[idx]}']
        
        coords_lst.append(np.array(curfilamentcoordslst[0])/voxelsize[0])
        coords_lst.append(np.array(curfilamentcoordslst[-1])/voxelsize[0])
        angle_lst.append(curfilamentangle)

    assert len(coords_lst) == len(angle_lst)*2

    coords = np.array(coords_lst).reshape(-1,3)
    
    # print('line56', coords)
    voxel_coords = coords

    coords_x = [ voxel_coords[i][0] for i in range(len(voxel_coords))]
    coords_y = [ voxel_coords[i][1] for i in range(len(voxel_coords))]
    coords_z = [ voxel_coords[i][2] for i in range(len(voxel_coords))]

    # mask_down_new = np.ones([mask_down.shape[2], mask_down.shape[1], mask_down.shape[0]])
    # mask_down_new = np.array(mask_down[2])

    if not (np.max(coords_x) <= mask_down.shape[0] and np.min(coords_x)>= 0) or not (np.max(coords_y) <= mask_down.shape[1] and np.min(coords_y)>= 0) or not (np.max(coords_z) <= mask_down.shape[2] and np.min(coords_z)>= 0):
        print('x',np.max(coords_x),np.min(coords_x), mask_down.shape[0])
        print('y',np.max(coords_y),np.min(coords_y), mask_down.shape[1])
        print('z',np.max(coords_z),np.min(coords_z), mask_down.shape[2])

    if len(coords_x) > 0:
        assert np.max(coords_x) <= mask_down.shape[0] and np.min(coords_x)>= 0
        assert np.max(coords_y) <= mask_down.shape[1] and np.min(coords_y)>= 0
        assert np.max(coords_z) <= mask_down.shape[2] and np.min(coords_z)>= 0
    else:
        print('no filament')

    img_up = mask_up
    reverse_img_up = img_up * (-1) + 1
    img_up_rev_edt = scipy.ndimage.distance_transform_edt(reverse_img_up)

    img_down = mask_down
    reverse_img_down = img_down * (-1) + 1
    img_down_rev_edt = scipy.ndimage.distance_transform_edt(reverse_img_down)

    ## angle_mask

    avesize = np.average(voxelsize)
    dist_up0 = [ img_up_rev_edt[int(voxel_coords[i][0])][int(voxel_coords[i][1])][int(voxel_coords[i][2])]  for i in range(len(voxel_coords)) ]
    dist_down0 = [ img_down_rev_edt[int(voxel_coords[i][0])][int(voxel_coords[i][1])][int(voxel_coords[i][2])]  for i in range(len(voxel_coords)) ]

    dist_up = [dist_up0[i] * avesize/10 for i in range(len(dist_up0))] # nm
    dist_down = [dist_down0[i] * avesize/10 for i in range(len(dist_down0))] # nm

    angle_pair_all_down, angle_all_down = [], []
    angle_pair_all_up, angle_all_up = [], []
    idx_lst_up, idx_lst_down = [], []
    dist_record_up, dist_record_down = [], [] # distance of each actin to PM
    
    count = 0
    for i in range(len(dist_up)):
        if int(dist_up[i]) <= threshold or int(dist_down[i]) <= threshold:
            count += 1


    for i in range(len(dist_up)):
        if int(dist_up[i]) <= threshold:
            # dist_min = dist_down[i] if dist_up[i] > dist_down[i] else dist_up[i] 
            idx_lst_up.append(i)
            dist_record_up.append(dist_up[i])

    for i in range(len(dist_down)):
        if int(dist_down[i]) <= threshold:
            # dist_min = dist_down[i] if dist_up[i] > dist_down[i] else dist_up[i] 
            idx_lst_down.append(i)
            dist_record_down.append(dist_down[i])

    for ii in range(len(idx_lst_down)-1):
        if (idx_lst_down[ii] % 2) == 0: 
            if dist_record_down[ii] > dist_record_down[ii]+1: 
                coord_00 = voxel_coords[idx_lst_down[ii]]
                coord_01 = voxel_coords[int(idx_lst_down[ii]+1)]
                angle_actin = angle_lst[int(idx_lst_down[ii]/2)]
            else:
                coord_00 = voxel_coords[int(idx_lst_down[ii]+1)]
                coord_01 = voxel_coords[idx_lst_down[ii]]
                angle_actin = angle_lst[int(idx_lst_down[ii]/2)]
        else: 
            if dist_record_down[ii] > dist_record_down[ii]-1: 
                coord_00 = voxel_coords[idx_lst_down[ii]]
                coord_01 = voxel_coords[idx_lst_down[ii]-1]
                angle_actin = angle_lst[int((idx_lst_down[ii]-1)/2)]
            else:
                coord_01 = voxel_coords[idx_lst_down[ii]]
                coord_00 = voxel_coords[idx_lst_down[ii]-1]
                angle_actin = angle_lst[int((idx_lst_down[ii]-1)/2)]

        for jj in range(len(idx_lst_down)):

            if (idx_lst_down[jj] % 2) == 0:
                if dist_record_down[jj] > dist_record_down[jj]+1: 
                    coord_10 = voxel_coords[idx_lst_down[jj]]
                    coord_11 = voxel_coords[int(idx_lst_down[jj]+1)]
                else:
                    coord_10 = voxel_coords[int(idx_lst_down[jj]+1)]
                    coord_11 = voxel_coords[idx_lst_down[jj]]
            else: 
                if dist_record_down[jj] > dist_record_down[jj]-1: 
                    coord_10 = voxel_coords[idx_lst_down[jj]]
                    coord_11 = voxel_coords[idx_lst_down[jj]-1]
                else:
                    coord_11 = voxel_coords[idx_lst_down[jj]]
                    coord_10 = voxel_coords[idx_lst_down[jj]-1]

            if idx_lst_down[ii] == idx_lst_down[jj]:
                continue
            else:
                
                dis_actin2actin = ((coord_01[0] - coord_11[0])**2 + (coord_01[1] - coord_11[1])**2 + (coord_01[2] - coord_11[2])**2)**0.5

                if dis_actin2actin < threshold_actin:                
                    vect1 = np.array(coord_00) - np.array(coord_01)
                    vect2 = np.array(coord_10) - np.array(coord_11)
                    

                    vector_dot_product = np.dot(vect1, vect2)
                    arccos = np.arccos(vector_dot_product / (np.linalg.norm(vect1) * np.linalg.norm(vect2)))
                    angle = np.degrees(arccos)

                    angle_pair_all_down.append(angle)
                    angle_all_down.append(angle_actin)

    for ii in range(len(idx_lst_up)-1):
        if (idx_lst_up[ii] % 2) == 0: # oushu
            if dist_record_up[ii] > dist_record_up[ii]+1: 
                coord_00 = voxel_coords[idx_lst_up[ii]]
                coord_01 = voxel_coords[int(idx_lst_up[ii]+1)]
                angle_actin = angle_lst[int(idx_lst_up[ii]/2)]
            else:
                coord_00 = voxel_coords[int(idx_lst_up[ii]+1)]
                coord_01 = voxel_coords[idx_lst_up[ii]]
                angle_actin = angle_lst[int(idx_lst_up[ii]/2)]
        else: # jishu 
            if dist_record_up[ii] > dist_record_up[ii]-1: 
                coord_00 = voxel_coords[idx_lst_up[ii]]
                coord_01 = voxel_coords[idx_lst_up[ii]-1]
                angle_actin = angle_lst[int((idx_lst_up[ii]-1)/2)]
            else:
                coord_01 = voxel_coords[idx_lst_up[ii]]
                coord_00 = voxel_coords[idx_lst_up[ii]-1]
                angle_actin = angle_lst[int((idx_lst_up[ii]-1)/2)]

        for jj in range(len(idx_lst_up)):
            if (idx_lst_up[jj] % 2) == 0:
                if dist_record_up[jj] > dist_record_up[jj]+1: 
                    coord_10 = voxel_coords[idx_lst_up[jj]]
                    coord_11 = voxel_coords[int(idx_lst_up[jj]+1)]
                else:
                    coord_10 = voxel_coords[int(idx_lst_up[jj]+1)]
                    coord_11 = voxel_coords[idx_lst_up[jj]]
            else: # jishu 
                if dist_record_up[jj] > dist_record_up[jj]-1: 
                    coord_10 = voxel_coords[idx_lst_up[jj]]
                    coord_11 = voxel_coords[idx_lst_up[jj]-1]
                else:
                    coord_11 = voxel_coords[idx_lst_up[jj]]
                    coord_10 = voxel_coords[idx_lst_up[jj]-1]

            if idx_lst_up[ii] == idx_lst_up[jj]:
                continue
            else:

                dis_actin2actin = ((coord_01[0] - coord_11[0])**2 + (coord_01[1] - coord_11[1])**2 + (coord_01[2] - coord_11[2])**2)**0.5

                if dis_actin2actin < threshold_actin:

                    vect1 = np.array(coord_00) - np.array(coord_01)
                    vect2 = np.array(coord_10) - np.array(coord_11)
                    
                    vector_dot_product = np.dot(vect1, vect2)
                    arccos = np.arccos(vector_dot_product / (np.linalg.norm(vect1) * np.linalg.norm(vect2)))
                    angle = np.degrees(arccos)

                    angle_pair_all_up.append(angle)
                    angle_all_up.append(angle_actin)

    return angle_all_up, angle_pair_all_up, angle_all_down, angle_pair_all_down # pm2actin, actin2actin


def pm_to_actin_number(filament_coords, filament_angle_data, mask_up, mask_down, voxelsize, threshold, threshold_actin):
    '''
    input: 
        filament coord: dict
        mask: array 
    actin_data_new, actin_angle_data, vesicle_mem_mask, voxel_size_xyz
    '''

    coords_lst = []
    angle_lst = []

    filamentnames = list(filament_coords.keys())
    for idx in range(len(filamentnames)):
        curfilamentcoordslst = filament_coords[f'{filamentnames[idx]}']
        curfilamentangle = filament_angle_data[f'{filamentnames[idx]}']
        
        coords_lst.append(np.array(curfilamentcoordslst[0])/voxelsize[0])
        coords_lst.append(np.array(curfilamentcoordslst[-1])/voxelsize[0])
        angle_lst.append(curfilamentangle)

    assert len(coords_lst) == len(angle_lst)*2

    coords = np.array(coords_lst).reshape(-1,3)
    
    voxel_coords = coords

    coords_x = [ voxel_coords[i][0] for i in range(len(voxel_coords))]
    coords_y = [ voxel_coords[i][1] for i in range(len(voxel_coords))]
    coords_z = [ voxel_coords[i][2] for i in range(len(voxel_coords))]

    if not (np.max(coords_x) <= mask_down.shape[0] and np.min(coords_x)>= 0) or not (np.max(coords_y) <= mask_down.shape[1] and np.min(coords_y)>= 0) or not (np.max(coords_z) <= mask_down.shape[2] and np.min(coords_z)>= 0):
        print('x',np.max(coords_x),np.min(coords_x), mask_down.shape[0])
        print('y',np.max(coords_y),np.min(coords_y), mask_down.shape[1])
        print('z',np.max(coords_z),np.min(coords_z), mask_down.shape[2])

    if len(coords_x) > 0:
        assert np.max(coords_x) <= mask_down.shape[0] and np.min(coords_x)>= 0
        assert np.max(coords_y) <= mask_down.shape[1] and np.min(coords_y)>= 0
        assert np.max(coords_z) <= mask_down.shape[2] and np.min(coords_z)>= 0
    else:
        print('no filament')

    img_up = mask_up
    reverse_img_up = img_up * (-1) + 1
    img_up_rev_edt = scipy.ndimage.distance_transform_edt(reverse_img_up)

    img_down = mask_down
    reverse_img_down = img_down * (-1) + 1
    img_down_rev_edt = scipy.ndimage.distance_transform_edt(reverse_img_down)

    avesize = np.average(voxelsize)
    dist_up0 = [ img_up_rev_edt[int(voxel_coords[i][0])][int(voxel_coords[i][1])][int(voxel_coords[i][2])]  for i in range(len(voxel_coords)) ]
    dist_down0 = [ img_down_rev_edt[int(voxel_coords[i][0])][int(voxel_coords[i][1])][int(voxel_coords[i][2])]  for i in range(len(voxel_coords)) ]

    dist_up = [dist_up0[i] * avesize/10 for i in range(len(dist_up0))] # nm
    dist_down = [dist_down0[i] * avesize/10 for i in range(len(dist_down0))] # nm
    
    count = 0
    for i in range(int(len(dist_up)/2)):
        if int(dist_up[int(2*i)]) <= threshold or int(dist_down[int(2*i)]) <= threshold or \
            int(dist_up[int(2*i+1)]) <= threshold or int(dist_down[int(2*i+1)]) <= threshold:
            count += 1
    
    return np.array(count)




def actin_to_actin_distance(filament_coord_extend_dict, imgsize, voxel_size_xyz):
    'input: filaments distance dict'
    edge = 500
    x_bound = [0 + edge, imgsize[0]*voxel_size_xyz[0] - edge]
    y_bound = [0 + edge, imgsize[1]*voxel_size_xyz[1] - edge]
    z_bound = [0 + edge, imgsize[2]*voxel_size_xyz[2] - edge]

    filament_coord_extend_lst = []
    filamentnames = list(filament_coord_extend_dict.keys())

    for idx in range(len(filamentnames)):
        # curfilamentcoordslst = filament_coord_extend_dict[f'{filamentnames[idx]}']
        curfilamentcoordslst = [ np.array(coord) for coord in filament_coord_extend_dict[f'{filamentnames[idx]}'] ]
        filament_coord_extend_lst.append(curfilamentcoordslst)


    filament_coord_lst = filament_coord_extend_lst
    filaments_vect_lst = [] 

    for single_filament_coord in filament_coord_lst:
        # coords1, coord2 = single_filament_coord[0], single_filament_coord[-1]
        vect = np.array(single_filament_coord[-1]) - np.array(single_filament_coord[0])
        unit_vect = vect / np.linalg.norm(vect)
        # unit_vect = vect
        filaments_vect_lst.append(unit_vect)

    filament_img_vect = np.array([0,0,0])
    for vect1 in filaments_vect_lst:
        filament_img_vect = filament_img_vect + vect1



    filament_vect_map_lst = []

    for i, curfilaments_coords in enumerate(filament_coord_extend_lst):
        rest_filament_coord_extend_lst = copy.deepcopy(filament_coord_extend_lst)
        _ = rest_filament_coord_extend_lst.pop(i)
        cur_filament_surr_coords = []
        vect_single_filament_map_lst = []
        vect_single_filament_angle_map_lst = []
        check_normvect = []
        curfilaments_coords_filtered= []
        for coords in curfilaments_coords:

            if (coords[0] >= x_bound[0] and  coords[0] <= x_bound[1]) or \
                (coords[1] >= y_bound[0] and  coords[1] <= y_bound[1]) or \
                (coords[2] >= z_bound[0] and  coords[2] <= z_bound[1]):

                curfilaments_coords_filtered.append(coords)
        curfilaments_coords = curfilaments_coords_filtered
        if len(curfilaments_coords) > 0:

            for j, rest_single_filament_coords in enumerate(rest_filament_coord_extend_lst):

                dists = cdist(curfilaments_coords, rest_single_filament_coords, metric='euclidean') # return array, size [a,b]
                # print(dists)
                filtered_coords = np.where(dists <= edge ) # 500 = 50nm  #100nm
                surrfilament_idx_check_lst = []
                for ii in range(len(filtered_coords[0])): 
                    if filtered_coords[1][ii] not in surrfilament_idx_check_lst:
                        surrfilament_idx_check_lst.append(filtered_coords[1][ii])
                        cur_filament_surr_coords.append(rest_single_filament_coords[filtered_coords[1][ii]])

            # get all surround coords and set vector
            for k, surround_coord in enumerate(cur_filament_surr_coords):
                distss = cdist([surround_coord], curfilaments_coords, metric='euclidean')
                min_d = np.min(distss)
                min_d_idx = np.where(distss == min_d)[1][0]

                cur_fila_coord = curfilaments_coords[min_d_idx] 
                cur_fila_correspond_coord = surround_coord

                if min_d_idx - 1 >= 0:
                    cur_direction_vect = curfilaments_coords[min_d_idx] - curfilaments_coords[min_d_idx-1]
                else:
                    cur_direction_vect = curfilaments_coords[min_d_idx+1] - curfilaments_coords[min_d_idx]
                
                vect_direction = np.dot(cur_direction_vect,filament_img_vect) 
                temp_vect_in_single_vect = cur_fila_correspond_coord - cur_fila_coord


                vect_single_filament_map_lst.append(temp_vect_in_single_vect) # rotate among cur filament coord


        filament_vect_map_lst.append(vect_single_filament_map_lst) # vect to surrounding in a flat

    return filament_vect_map_lst



def actin_to_microtube_distance(filament_coord_extend_dict, MT_coord_extend_dict):
    'input: filaments distance dict'
    
    filament_coord_extend_lst = []
    filamentnames = list(filament_coord_extend_dict.keys())

    for idx in range(len(filamentnames)):
        curfilamentcoordslst = [ np.array(coord) for coord in filament_coord_extend_dict[f'{filamentnames[idx]}'] ]
        filament_coord_extend_lst.append(curfilamentcoordslst)


    MT_coord_extend_lst = [] # [mt1 coords, mt2 coords, ..., mt3 coords ]
    MTnames = list(MT_coord_extend_dict.keys())

    for idx in range(len(MTnames)):
        # curfilamentcoordslst = filament_coord_extend_dict[f'{filamentnames[idx]}']
        curMTcoordslst = [ np.array(coord) for coord in MT_coord_extend_dict[f'{MTnames[idx]}'] ]
        MT_coord_extend_lst.append(curMTcoordslst) 

    



    MT_coord_lst = MT_coord_extend_lst
    MT_vect_lst = [] 

    for single_MT_coord in MT_coord_lst:
        # coords1, coord2 = single_filament_coord[0], single_filament_coord[-1]
        vect = np.array(single_MT_coord[-1]) - np.array(single_MT_coord[0])
        unit_vect = vect / np.linalg.norm(vect)
        # unit_vect = vect
        MT_vect_lst.append(unit_vect)

    MT_img_vect = np.array([0,0,0])
    for vect1 in MT_vect_lst:
        MT_img_vect = MT_img_vect + vect1





    # calculate dist map
    filament_vect_2_MT_map_lst = []

    for i, curMT_coords in enumerate(MT_coord_extend_lst):
        cur_filament_surr_coords = []
        vect_single_filament_map_lst = []
        vect_single_filament_angle_map_lst = []
        check_normvect = []


        for j, single_filament_coords in enumerate(filament_coord_extend_lst):
            dists = cdist(curMT_coords, single_filament_coords, metric='euclidean') # return array, size [a,b]

            filtered_coords = np.where(dists <= 1000 ) #100nm
            surrfilament_idx_check_lst = []
            for ii in range(len(filtered_coords[0])): 
                if filtered_coords[1][ii] not in surrfilament_idx_check_lst:
                    surrfilament_idx_check_lst.append(filtered_coords[1][ii])
                    cur_filament_surr_coords.append(single_filament_coords[filtered_coords[1][ii]])

        # get all surround coords and set vector
        for k, surround_coord in enumerate(cur_filament_surr_coords):
            distss = cdist([surround_coord], curMT_coords, metric='euclidean')
            min_d = np.min(distss)
            min_d_idx = np.where(distss == min_d)[1][0]

            cur_fila_coord = curMT_coords[min_d_idx] 
            cur_fila_correspond_coord = surround_coord

            if min_d_idx > 1:
                cur_direction_vect = curMT_coords[min_d_idx] - curMT_coords[min_d_idx-1]
            else:
                cur_direction_vect = curMT_coords[min_d_idx+1] - curMT_coords[min_d_idx]
            
            # depend on single filament direction to )filaments overall direction
            vect_direction = np.dot(cur_direction_vect, MT_img_vect) 
            temp_vect_in_single_vect = cur_fila_correspond_coord - cur_fila_coord

            vect_single_filament_map_lst.append(temp_vect_in_single_vect) # rotate among cur filament coord


        filament_vect_2_MT_map_lst.append(vect_single_filament_map_lst) # vect to surrounding in a flat
    return filament_vect_2_MT_map_lst



def MTdistance2mem(mt_coords, mask, voxelsize):
    '''
    input: 
        mt coord: dict
        mask: array 
    '''
    filament_coords = mt_coords
    filament_coords_modified ={}    
    filamentnames = list(filament_coords.keys())
    for idx in range(len(filamentnames)):
        curfilamentcoordslst = filament_coords[f'{filamentnames[idx]}']
        filament_coords_modified[f'{filamentnames[idx]}'] =[]
        for i, coord_new in enumerate(curfilamentcoordslst):
            coord_new = [coord_new[0]/voxelsize[0], coord_new[1]/voxelsize[1], coord_new[2]/voxelsize[2]]
            filament_coords_modified[f'{filamentnames[idx]}'].append(coord_new)




    coords_lst = []
    for idx in range(len(filamentnames)):
        curfilamentcoordslst = filament_coords_modified[f'{filamentnames[idx]}']
        coords_lst.extend(curfilamentcoordslst)

    coords = np.array(coords_lst).reshape(-1,3)

    voxel_coords = coords

    coords_x = [ voxel_coords[i][0] for i in range(len(voxel_coords))]
    coords_y = [ voxel_coords[i][1] for i in range(len(voxel_coords))]
    coords_z = [ voxel_coords[i][2] for i in range(len(voxel_coords))]

    
    print('mt coords range', [[np.round(np.min(coords_x)), np.round(np.max(coords_x))],[np.round(np.min(coords_y)), np.round(np.max(coords_y))],[np.round(np.min(coords_z)), np.round(np.max(coords_z))]] )
    print('mask size', mask.shape)
    print('coords in total', len(voxel_coords))


    assert np.max(coords_x) <= mask.shape[0] and np.min(coords_x)>= 0
    assert np.max(coords_y) <= mask.shape[1] and np.min(coords_y)>= 0
    assert np.max(coords_z) <= mask.shape[2] and np.min(coords_z)>= 0


    img01 = mask
    reverse_img01 = img01 * (-1) + 1
    img01_rev_edt = scipy.ndimage.distance_transform_edt(reverse_img01)

    avesize = np.average(voxelsize)
    dist = [ img01_rev_edt[int(voxel_coords[i][0])][int(voxel_coords[i][1])][int(voxel_coords[i][2])]  for i in range(len(voxel_coords)) ]
    dist = [dist[i] * avesize/10 for i in range(len(dist))] # nm
    print('distance num', len(dist))

    return dist


def actinangle2dist(filament_coord_extend_dict, imgsize, voxel_size_xyz):
    'input: filaments distance dict'
    edge = 0.
    x_bound = [0 + edge, imgsize[0]*voxel_size_xyz[0] - edge]
    y_bound = [0 + edge, imgsize[1]*voxel_size_xyz[1] - edge]
    z_bound = [0 + edge, imgsize[2]*voxel_size_xyz[2] - edge]

    filament_coord_extend_lst = []
    filamentnames = list(filament_coord_extend_dict.keys())

    for idx in range(len(filamentnames)):
        # curfilamentcoordslst = filament_coord_extend_dict[f'{filamentnames[idx]}']
        curfilamentcoordslst = [ np.array(coord) for coord in filament_coord_extend_dict[f'{filamentnames[idx]}'] ]
        filament_coord_extend_lst.append(curfilamentcoordslst)

    filament_coord_lst = filament_coord_extend_lst
    filaments_vect_lst = []  # vect for every point in every filament

    for single_filament_coord in filament_coord_lst:
        filament_vectlst = []
        for jj, point in enumerate(single_filament_coord):
            if jj < 1:
                tempvect = single_filament_coord[jj+1] - point
            else:
                tempvect = point- single_filament_coord[jj-1]
            # unit_vect = vect / np.linalg.norm(vect)
            filament_vectlst.append(tempvect)  
        filaments_vect_lst.append(filament_vectlst)  # not unit vect 


    filaments_min_d_lst = []
    filaments_angle_lst = []

    for i, curfilaments_coords in enumerate(filament_coord_extend_lst):
        curfilaments_vect = filaments_vect_lst[i]

        rest_filament_coord_extend_lst = copy.deepcopy(filament_coord_extend_lst)
        rest_filament_vect_extend_lst = copy.deepcopy(filaments_vect_lst)
        _ = rest_filament_coord_extend_lst.pop(i)
        _ = rest_filament_vect_extend_lst.pop(i)

        # curfilaments_coords = [point1, point2, ...]
        filament_dist_lst = []
        filament_angle_lst = []
        for ll, point in enumerate(curfilaments_coords):
            vect1 = curfilaments_vect[ll]
            min_dist = np.inf
            cur_angle = np.inf
            for kk, filaments_coords in enumerate(rest_filament_coord_extend_lst):
                distss = cdist([point], filaments_coords, metric='euclidean')
                min_d = np.min(distss)
                min_d_idx = np.where(distss == min_d)[1][0]

                if min_d < min_dist:
                    min_dist = min_d
                    vect2 = rest_filament_vect_extend_lst[kk][min_d_idx]
                    vector_dot_product = np.dot(vect1, vect2)
                    arccos = np.arccos(vector_dot_product / (np.linalg.norm(vect1) * np.linalg.norm(vect2)))
                    angle = np.degrees(arccos)
                    if angle >= 180:
                        cur_angle = 0
                    elif angle > 90:
                        cur_angle = 180 - angle
                    else:
                        cur_angle = angle

            filament_dist_lst.append(min_dist)
            filament_angle_lst.append(cur_angle)

        filaments_min_d_lst.extend(filament_dist_lst)
        filaments_angle_lst.extend(filament_angle_lst)

    return  filaments_min_d_lst, filaments_angle_lst

