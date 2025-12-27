import sys
import os

# ================= 1. 基础几何参数 =================
B = 0.03            # 建筑宽 30mm
H = 2.67 * B        # 建筑高
W = B               # 街道宽
Domain_H = 25 * B   # 域高

L_inlet = 33.3 * B
L_outlet = 121.3 * B
Num_Rows = 6

# ================= 2. 网格参数 =================
# 目标：所有壁面（地面、屋顶、侧墙、迎风面、背风面）第一层网格均为 0.4mm
First_Layer_Thickness = 0.0004 

# 网格数
# Z方向 (垂直)
N_z_build = 35      # 0-H
N_z_sky = 65        # H-Top

# Y方向 (横向)
N_y_inner = 14      # 0-0.5B
N_y_outer = 14      # 0.5B-B

# X方向 (流向) 这里把街道拆分成前后两个部分
N_x_inlet = 60
N_x_build = 16      # 建筑顶部
N_x_street_half = 8 # 半个街道 (总共16)
N_x_outlet = 120

# ================= 3. 计算网格simpleGrading参数 =================
def calc_ratio(total_len, n_cells, first_layer):    # (长度, 网格数, 第一层网格长度)
    """网格变化为等比数列"""
    
    if n_cells <= 1: return 1
    
    def get_len(r):
        if abs(r-1) < 1e-9: return n_cells * first_layer
        return first_layer * (pow(r, n_cells) - 1) / (r - 1)

    low, high = 1.00001, 2.0
    r_sol = 1.0
    for _ in range(100):
        mid = (low + high) / 2
        if get_len(mid) < total_len: low = mid
        else: high = mid
    r_sol = (low + high) / 2
    return pow(r_sol, n_cells - 1)  # 返回值为最后一层网格与第一层网格的比

# ================= 4. 预计算各方向simpleGrading参数 =================
R_z_ground = calc_ratio(H, N_z_build, First_Layer_Thickness)
R_z_sky = calc_ratio(Domain_H - H, N_z_sky, First_Layer_Thickness)

_r_temp = calc_ratio(0.5*B, N_y_inner, First_Layer_Thickness)
R_y_wall = 1.0 / _r_temp 

R_x_inlet = 1.0 / calc_ratio(L_inlet, N_x_inlet, First_Layer_Thickness)
R_x_outlet = calc_ratio(L_outlet, N_x_outlet, First_Layer_Thickness)

R_x_street_wake = calc_ratio(0.5*W, N_x_street_half, First_Layer_Thickness)
R_x_street_appr = 1.0 / R_x_street_wake

R_x_build = 1.0 

# ================= 5. 生成坐标点 =================
x_coords = []
# 0. Inlet段起始点
x_coords.append(-L_inlet)

curr_x = 0.0
block_types = ['inlet'] # Block 0 (Inlet -> Build1 Start)

for i in range(Num_Rows):
    # 1. 建筑段
    # 添加建筑起始坐标
    x_coords.append(curr_x)     
    
    # 记录该段为建筑
    block_types.append('build') 
    
    curr_x += B
    # 添加建筑结束坐标
    x_coords.append(curr_x)     
    
    # 2. 街道段
    if i < Num_Rows - 1:
        # 左半街道
        block_types.append('street_wake')
        curr_x += 0.5 * W
        x_coords.append(curr_x) # 街道中点
        
        # 右半街道
        block_types.append('street_appr')
        curr_x += 0.5 * W 


# 3. Outlet段
block_types.append('outlet')
x_coords.append(curr_x + L_outlet)

y_coords = [0, 0.5 * B, B]
z_coords = [0, H, Domain_H]

# ================= 6. 生成 BlockMeshDict =================
def get_v(ix, iy, iz, nx, ny):
    return iz * (nx * ny) + iy * nx + ix

def generate():

    print(r"""/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2012                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
convertToMeters 1;

vertices
(""")
    n_x = len(x_coords)
    n_y = len(y_coords)
    n_z = len(z_coords)

    for z in z_coords:
        for y in y_coords:
            for x in x_coords:
                print(f"    ({x} {y} {z})")
    print(");")

    print("\nblocks\n(")

    # blocks
    for i in range(n_x - 1):
        b_type = block_types[i]
        
        if b_type == 'inlet':
            nx, gx = N_x_inlet, R_x_inlet
        elif b_type == 'outlet':
            nx, gx = N_x_outlet, R_x_outlet
        elif b_type == 'build':
            nx, gx = N_x_build, R_x_build
        elif b_type == 'street_wake':
            nx, gx = N_x_street_half, R_x_street_wake
        elif b_type == 'street_appr':
            nx, gx = N_x_street_half, R_x_street_appr
        
        is_solid = (b_type == 'build')

        for j in range(n_y - 1):
            is_inner_y = (j == 0)
            ny = N_y_inner if is_inner_y else N_y_outer
            gy = R_y_wall if is_inner_y else 1.0 

            for k in range(n_z - 1):
                is_ground = (k == 0)
                nz = N_z_build if is_ground else N_z_sky
                gz = R_z_ground if is_ground else R_z_sky

                if is_solid and is_inner_y and is_ground:
                    continue

                v0 = get_v(i,   j,   k,   n_x, n_y)
                v1 = get_v(i+1, j,   k,   n_x, n_y)
                v2 = get_v(i+1, j+1, k,   n_x, n_y)
                v3 = get_v(i,   j+1, k,   n_x, n_y)
                v4 = get_v(i,   j,   k+1, n_x, n_y)
                v5 = get_v(i+1, j,   k+1, n_x, n_y)
                v6 = get_v(i+1, j+1, k+1, n_x, n_y)
                v7 = get_v(i,   j+1, k+1, n_x, n_y)

                print(f"    hex ({v0} {v1} {v2} {v3} {v4} {v5} {v6} {v7}) "
                      f"({nx} {ny} {nz}) simpleGrading ({gx} {gy} {gz})")

    print(");")
    print("\nedges\n();")
    print("\nboundary\n(")

    # 1. Inlet
    print("    inlet { type patch; faces (")
    for j in range(n_y-1):
        for k in range(n_z-1):
            v0, v3 = get_v(0, j, k, n_x, n_y), get_v(0, j+1, k, n_x, n_y)
            v4, v7 = get_v(0, j, k+1, n_x, n_y), get_v(0, j+1, k+1, n_x, n_y)
            print(f"        ({v0} {v4} {v7} {v3})")
    print("    ); }")

    # 2. Outlet
    print("    outlet { type patch; faces (")
    for j in range(n_y-1):
        for k in range(n_z-1):
            v1, v2 = get_v(n_x-1, j, k, n_x, n_y), get_v(n_x-1, j+1, k, n_x, n_y)
            v5, v6 = get_v(n_x-1, j, k+1, n_x, n_y), get_v(n_x-1, j+1, k+1, n_x, n_y)
            print(f"        ({v1} {v2} {v6} {v5})")
    print("    ); }")

    # 3. Top
    print("    top { type symmetry; faces (")
    for i in range(n_x-1):
        for j in range(n_y-1):
            v4, v5 = get_v(i, j, n_z-1, n_x, n_y), get_v(i+1, j, n_z-1, n_x, n_y)
            v6, v7 = get_v(i+1, j+1, n_z-1, n_x, n_y), get_v(i, j+1, n_z-1, n_x, n_y)
            print(f"        ({v4} {v5} {v6} {v7})")
    print("    ); }")

    # 4. Buildings & Ground
    print("    buildings_and_ground { type wall; faces (")
    # Ground
    for i in range(n_x-1):
        b_type = block_types[i]
        is_solid = (b_type == 'build')
        for j in range(n_y-1):
            if is_solid and (j == 0): continue
            v0, v1 = get_v(i, j, 0, n_x, n_y), get_v(i+1, j, 0, n_x, n_y)
            v2, v3 = get_v(i+1, j+1, 0, n_x, n_y), get_v(i, j+1, 0, n_x, n_y)
            print(f"        ({v0} {v3} {v2} {v1})")
    # Building Faces
    for i in range(n_x-1):
        b_type = block_types[i]
        is_solid = (b_type == 'build')
        if is_solid:
            # Roof
            v4, v5 = get_v(i, 0, 1, n_x, n_y), get_v(i+1, 0, 1, n_x, n_y)
            v6, v7 = get_v(i+1, 1, 1, n_x, n_y), get_v(i, 1, 1, n_x, n_y)
            print(f"        ({v4} {v7} {v6} {v5})")
            # Side Wall
            v0, v3 = get_v(i, 1, 0, n_x, n_y), get_v(i+1, 1, 0, n_x, n_y)
            v4, v7 = get_v(i, 1, 1, n_x, n_y), get_v(i+1, 1, 1, n_x, n_y)
            print(f"        ({v0} {v4} {v7} {v3})")
            # Front/Back Walls
            v1_l, v2_l = get_v(i, 0, 0, n_x, n_y), get_v(i, 1, 0, n_x, n_y)
            v5_l, v6_l = get_v(i, 0, 1, n_x, n_y), get_v(i, 1, 1, n_x, n_y)
            print(f"        ({v1_l} {v2_l} {v6_l} {v5_l})")
            v0_r, v3_r = get_v(i+1, 0, 0, n_x, n_y), get_v(i+1, 1, 0, n_x, n_y)
            v4_r, v7_r = get_v(i+1, 0, 1, n_x, n_y), get_v(i+1, 1, 1, n_x, n_y)
            print(f"        ({v0_r} {v4_r} {v7_r} {v3_r})")
    print("    ); }")

    # 5. Symmetry
    print("    symmetry { type symmetry; faces (")
    # Y=0
    for i in range(n_x-1):
        b_type = block_types[i]
        is_solid = (b_type == 'build')
        for k in range(n_z-1):
            if is_solid and (k == 0): continue
            v0, v3 = get_v(i, 0, k, n_x, n_y), get_v(i+1, 0, k, n_x, n_y)
            v4, v7 = get_v(i, 0, k+1, n_x, n_y), get_v(i+1, 0, k+1, n_x, n_y)
            print(f"        ({v0} {v4} {v7} {v3})")
    # Y=Max
    for i in range(n_x-1):
        for k in range(n_z-1):
            v1, v2 = get_v(i, n_y-1, k, n_x, n_y), get_v(i+1, n_y-1, k, n_x, n_y)
            v5, v6 = get_v(i, n_y-1, k+1, n_x, n_y), get_v(i+1, n_y-1, k+1, n_x, n_y)
            print(f"        ({v1} {v2} {v6} {v5})")
    print("    ); }")

    print(");")
    print("mergePatchPairs();")

if __name__ == "__main__":
    if os.path.exists("system"):
        output_file = os.path.join("system", "blockMeshDict")
    else:
        output_file = "blockMeshDict"
    
    with open(output_file, "w") as f:
        sys.stdout = f
        generate()
    
    sys.stdout = sys.__stdout__
    print(f"网格生成完毕>>>{output_file}")