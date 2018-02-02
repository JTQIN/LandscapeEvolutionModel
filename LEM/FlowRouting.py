import numpy as np
import timeit


def slope_lap(Z , L):
    Ni = Z.shape[0]
    Nj = Z.shape[1]

    Z_2 = np.append(Z[0 , :].reshape(1 , Nj) , Z[0 : Ni - 1  , :].reshape(Ni - 1 , Nj) , axis = 0)
    Z_4 = np.append(Z[: , 1].reshape(Ni , 1) , Z[: , 0 : Nj - 1].reshape(Ni , Nj - 1) , axis = 1)
    Z_6 = np.append(Z[: , 1 : Nj].reshape(Ni , Nj - 1) , Z[: , Nj - 2].reshape(Ni , 1) , axis = 1)
    Z_8 = np.append(Z[1 : Ni  , :].reshape(Ni - 1 , Nj) , Z[Ni - 1 , :].reshape(1 , Nj) , axis = 0)

    Z_1 = np.append(Z_2[: , 0].reshape(Ni , 1) , Z_2[: , 0 : Nj - 1].reshape(Ni , Nj - 1) , axis = 1)
    Z_3 = np.append(Z_2[: , 1 : Nj].reshape(Ni , Nj - 1) , Z_2[: , Nj - 1].reshape(Ni , 1) , axis = 1)
    Z_5 = np.copy(Z)
    Z_7 = np.append(Z_8[: , 0].reshape(Ni , 1) , Z_8[: , 0 : Nj - 1].reshape(Ni , Nj - 1) , axis = 1)
    Z_9 = np.append(Z_8[: , 1 : Nj].reshape(Ni , Nj - 1) , Z_8[: , Nj - 1].reshape(Ni , 1) , axis = 1)


    #A = ((Z_1  + Z_3  + Z_7  + Z_9 ) / 4  - (Z_2  + Z_4  + Z_6  + Z_8 ) / 2 + Z_5 ) / L ** 4
    #B = ((Z_1  + Z_3  - Z_7  - Z_9 ) /4 - (Z_2  - Z_8 ) /2) / L ** 3
    #C = ((-Z_1  + Z_3  - Z_7  + Z_9 ) /4 + (Z_4  - Z_6 ) / 2) / L ** 3
    D = ((Z_4  + Z_6 )/2 - Z_5 ) / L ** 2
    E = ((Z_2  + Z_8 ) /2 - Z_5 ) / L ** 2
    F = (-Z_1  + Z_3  + Z_7  - Z_9 ) / (4 * L ** 2)
    G = (-Z_4  + Z_6 ) / (2 * L)
    H = (Z_2  - Z_8 ) / (2 * L)
    #I = Z_5

    del Z_1 , Z_2 , Z_3 ,Z_4 , Z_5 , Z_6 , Z_7 , Z_8 , Z_9

    hx = G
    hxx = 2 * D
    hxy = F
    hy = H
    hyy = 2 * E

    del G , D , F ,H , E

    Lap = hxx + hyy
    Slope = (hx** 2 + hy** 2) ** 0.5

    del hxx , hyy, hx , hy
    
    return Lap , Slope

def fill(Z):
    Ni = Z.shape[0]
    Nj = Z.shape[1]
    border = np.zeros((Ni , Nj)).astype(np.int8)
    border[: , 0] = 1
    border[: , -1] = 1
    border[0 , :] = 1
    border[-1 , :] = 1
    W = np.where(border == 0 , np.max(Z) + 10 , Z)
    eps = 0.000001
    smt_done = 1
    while smt_done == 1:
        smt_done = 0
        proc_ext = np.where((W > Z)& (border==0) , 1 , 0).astype(np.int8)
        list_nb = neighbour_list(W , -1)
        for nb in range(0 , 8):
            case_1 = np.where((proc_ext==1)&(Z >= list_nb[nb] + eps) , 1 , 0).astype(np.int8)
            case_2 = np.where((proc_ext==1)&(case_1==0)&(W > list_nb[nb] + eps) , 1 , 0).astype(np.int8)
            Wnew = np.where(case_1 == 1 , Z , W)
            Wnew = np.where(case_2 == 1 , list_nb[nb] + eps , Wnew)
            if np.sum(np.abs(W - Wnew)) > 0:
                smt_done = 1
            W = np.copy(Wnew)
            list_nb = neighbour_list(W , -1)
    return W
    
    

def neighbour_list(A , border_value):

    Ni = A.shape[0]
    Nj = A.shape[1]
    
    A_7 = np.append(border_value * np.ones((1 , Nj)) , A[0 : Ni - 1  , :].reshape(Ni - 1 , Nj) , axis = 0)
    A_5 = np.append(border_value * np.ones((Ni , 1)) , A[: , 0 : Nj - 1].reshape(Ni , Nj - 1) , axis = 1)
    A_1 = np.append(A[: , 1 : Nj].reshape(Ni , Nj - 1) , border_value * np.ones((Ni , 1)) , axis = 1)
    A_3 = np.append(A[1 : Ni  , :].reshape(Ni - 1 , Nj) , border_value * np.ones((1 , Nj)) , axis = 0)

    A_6 = np.append(border_value * np.ones((Ni , 1)) , A_7[: , 0 : Nj - 1].reshape(Ni , Nj - 1) , axis = 1)
    A_8 = np.append(A_7[: , 1 : Nj].reshape(Ni , Nj - 1) , border_value * np.ones((Ni , 1)) , axis = 1)
    A_4 = np.append(border_value * np.ones((Ni , 1)) , A_3[: , 0 : Nj - 1].reshape(Ni , Nj - 1) , axis = 1)
    A_2 = np.append(A_3[: , 1 : Nj].reshape(Ni , Nj - 1) , border_value * np.ones((Ni , 1)) , axis = 1)

    list_neighbour = []

    list_neighbour.append(A_1)
    list_neighbour.append(A_2)
    list_neighbour.append(A_3)
    list_neighbour.append(A_4)
    list_neighbour.append(A_5)
    list_neighbour.append(A_6)
    list_neighbour.append(A_7)
    list_neighbour.append(A_8)

    return list_neighbour

def neighbour_list_flowdir(A):
    list_neighbour = []

    list_neighbour.append(A_1)
    list_neighbour.append(A_2)
    list_neighbour.append(A_3)
    list_neighbour.append(A_4)
    list_neighbour.append(A_5)
    list_neighbour.append(A_6)
    list_neighbour.append(A_7)
    list_neighbour.append(A_8)
    list_neighbour.append(A_9)

    return list_neighbour

def rs_cal(e0, e1 , e2,dx):
    s1 = np.where(e0 == e1 , 1*10**(-10) , (e0 - e1) / dx)
    s2 = (e1 - e2) / dx
    r = np.where(s1 > 0 , np.arctan(s2/s1) , 1000.0)
    s = np.where(np.minimum(e2 , e1) <= e0 ,   np.sqrt(s1**2+s2**2), -1.0)
    s = np.where((e2 == e1) &  (e2 == e0) ,   -2.0 , s)
    r_new = np.where((r < 0) & (s >0) , 0 , r)
    r_new = np.where((r > np.arctan(dx/dx)) & (s >0) , np.arctan(dx/dx) , r_new)
    s_new = np.where((r < 0)  & (s >0), s1 , s)
    s_new = np.where((r > np.arctan(dx/dx)) & (s >0) , (e0 - e2)/(dx * np.sqrt(2)) , s_new)
    return s_new , r_new

def rg_cal(r, ac , af):
    rg = af * r + ac * np.pi/2.
    return rg

def rg_update(s , r , ac , af , s_max , rg):
    rg = np.where(s>s_max , rg_cal(r , ac , af) , rg)
    s_max = np.maximum(s , s_max)
    return rg , s_max

def Dinf_flowdir(Z , dx):

    # There is free drinage at the boundaries

    Ni = Z.shape[0]
    Nj = Z.shape[1]
    
    Z_border = np.min(Z) -1.
    Z_2 = np.append(Z_border * np.ones((1 , Nj)) , Z[0 : Ni - 1  , :].reshape(Ni - 1 , Nj) , axis = 0)
    Z_4 = np.append(Z_border * np.ones((Ni , 1)) , Z[: , 0 : Nj - 1].reshape(Ni , Nj - 1) , axis = 1)
    Z_6 = np.append(Z[: , 1 : Nj].reshape(Ni , Nj - 1) , Z_border * np.ones((Ni , 1)) , axis = 1)
    Z_8 = np.append(Z[1 : Ni  , :].reshape(Ni - 1 , Nj) , Z_border * np.ones((1 , Nj)) , axis = 0)
    Z_1 = np.append(Z_2[: , 0].reshape(Ni , 1) , Z_2[: , 0 : Nj - 1].reshape(Ni , Nj - 1) , axis = 1)
    Z_3 = np.append(Z_2[: , 1 : Nj].reshape(Ni , Nj - 1) , Z_2[: , Nj - 1].reshape(Ni , 1) , axis = 1)
    Z_5 = np.copy(Z)
    Z_7 = np.append(Z_8[: , 0].reshape(Ni , 1) , Z_8[: , 0 : Nj - 1].reshape(Ni , Nj - 1) , axis = 1)
    Z_9 = np.append(Z_8[: , 1 : Nj].reshape(Ni , Nj - 1) , Z_8[: , Nj - 1].reshape(Ni , 1) , axis = 1)
    
    s1 , r1 = rs_cal(Z_5, Z_6 , Z_3 ,dx)
    rg = rg_cal(r1 , 0.0 , 1.0)
    s_max = np.copy(s1)

    s2 , r2 = rs_cal(Z_5, Z_2 , Z_3 ,dx)
    rg , s_max = rg_update(s2 , r2 , 1.0 , -1.0 , s_max , rg)

    s3 , r3 = rs_cal(Z_5, Z_2 , Z_1 ,dx)
    rg , s_max = rg_update(s3 , r3 , 1.0 , 1.0 , s_max , rg)

    s4 , r4 = rs_cal(Z_5, Z_4 , Z_1 ,dx)
    rg , s_max = rg_update(s4 , r4 , 2.0 , -1.0 , s_max , rg)

    s5 , r5 = rs_cal(Z_5, Z_4 , Z_7 ,dx)
    rg , s_max = rg_update(s5 , r5 , 2.0 , 1.0 , s_max , rg)

    s6 , r6 = rs_cal(Z_5, Z_8 , Z_7 ,dx)
    rg , s_max = rg_update(s6 , r6 , 3.0 , -1.0 , s_max , rg)

    s7 , r7 = rs_cal(Z_5, Z_8 , Z_9 ,dx)
    rg , s_max = rg_update(s7 , r7 , 3.0 , 1.0 , s_max , rg)

    s8 , r8 = rs_cal(Z_5, Z_6 , Z_9 ,dx)
    rg , s_max = rg_update(s8 , r8 , 4.0 , -1.0 , s_max , rg)

    rg = np.where(s_max >= 0 , rg , -1.)

    ## find sinks
    row_array = np.array([np.arange(0 , Ni),]*Nj).transpose()
    col_array = np.array([np.arange(0 , Nj),]*Ni)
    New_row_array = row_array - np.sin(rg)
    New_col_array = col_array + np.cos(rg)
    sink = np.where((rg < 0)|(New_row_array<0)|(New_row_array>=Ni)|(New_col_array<0)|(New_col_array>=Nj) , 1 , 0)
    sink_row = sink.nonzero()[0]
    sink_col = sink.nonzero()[1]
    pit = np.where((rg < 0) , 1 , 0)
    
    return rg , s_max, sink_row , sink_col, pit

def acc_recur(i , j):
    global upArea, done, nei_row, nei_col , p_list
    if done[i , j] == 0:
        done[i , j] = 1
        for nb in range (0 , 8):
            p = p_list[nb][i , j]
            if p > 0:
                ii = nei_row[nb]+i
                jj = nei_col[nb]+j
                acc_recur(ii , jj)
                upArea[i , j] = upArea[i , j] + p * upArea[ii , jj]

def Dinf_flowacc(fdir , sink_row , sink_col):
    global upArea, done, nei_row, nei_col , p_list
    nei_row = np.array([0 , 1 , 1 ,  1 ,  0 , -1 , -1, -1])
    nei_col = np.array([1 , 1 , 0 , -1 , -1 , -1 , 0 , 1])
    nei_facet_1 = np.array([4 ,3 , 2 , 1 , 8 , 7 , 6 , 5])
    nei_facet_2 = np.array([5 ,4 , 3 , 2 , 1 , 8 , 7 , 6])

    fdir_list = neighbour_list(fdir , 1000.)
    upst_list = []
    facet_list = []
    p_list = []
    upst_code = np.zeros_like(fdir).astype(np.uint16)
    for nb in range(0 , 8):
        facet_list.append( np.fix(fdir_list[nb] / (np.pi/4.)) + 1)
        upst_list.append(((facet_list[nb] == nei_facet_1[nb])|(facet_list[nb] == nei_facet_2[nb])).astype(np.int8))
        upst_code = upst_code + upst_list[nb] * (2**nb)
        p_list.append(np.where((facet_list[nb] == nei_facet_1[nb])|(facet_list[nb] == nei_facet_2[nb]) , \
                 np.abs((facet_list[nb] == nei_facet_2[nb]) - (fdir_list[nb] - (facet_list[nb] - 1) * np.pi/4.) / (np.pi/4.)), -1))


    upArea = np.ones_like(fdir)
    upArea[0 , :] = 0
    upArea[-1 , :] = 0
    upArea[: , 0] = 0
    upArea[: , -1] = 0
    
    done = np.zeros_like(fdir).astype(np.int8)
    for sn in range(0 , sink_row.shape[0]):
        row_s = sink_row[sn]
        col_s = sink_col[sn]
        acc_recur(row_s , col_s)
    return upArea
        

    
def flow_routing(Z , dx , option):
    if option == 'Fill':
        Z = fill(Z)
        fdir , slope , sink_row , sink_col = Dinf_flowdir(Z , dx)
        upArea = Dinf_flowacc(fdir, sink_row , sink_col)
    elif option == 'NoFill':
        fdir , slope , sink_row , sink_col , pit = Dinf_flowdir(Z , dx)
        upArea = Dinf_flowacc(fdir, sink_row , sink_col)
    return Z , upArea , fdir , pit

#### Main function
##DEM_raster = 'C:\\Data\\Codes\\WorkingCodes\\LEM\\FlowRout\\In\\Fill_Dem.tif'
##dx = 5.
##Z = arcpy.RasterToNumPyArray(DEM_raster , nodata_to_value=0)
##
##start = timeit.default_timer()
##Z , upArea = flow_routing(Z , dx , 'Fill')
##stop = timeit.default_timer()
##print stop - start



