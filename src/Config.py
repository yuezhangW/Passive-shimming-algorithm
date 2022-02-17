class Config:
    FILES_NAMES = ['6009map1.txt', '9006map1.txt']  # dam文件的文件名
    FILES_INDEX = 0  # 选取的文件序号

    NMR_SIZE = [40, 24]  # NMR匀场片的分布数量（phi方向，z方向）
    NMR_THICKNESS = [0.0001, 0.007]  # NMR匀场片最小厚度和空腔最大厚度 (单位：m)
    NMR_MZ = 2.06  # Mz强度
    NMR_PIECE_AREA = 0.05 * 0.04  # 匀场片的单位面积 (单位：m^2)
    NMR_Z = 1.74  # NMR的z方向长度（单位：m）
    NMR_R = 0.485  # NMR圆柱内半径（单位：m）

    DSV_R = 0.225  # DSV区域球半径（单位：m）

    SHIMMING_KAPPA = 0.001
    SHIMMING_TARGET_INHOMOGENEITY_MIN = 15  # 搜索线性解的最低目标不均匀度
    SHIMMING_TARGET_INHOMOGENEITY_MAX = 23  # 搜索线性解的最高目标不均匀度
    SHIMMING_STEP = 200000  # 搜索步长

    TABLE_AXIS_HOM = "Unhomogeneity(ppm)"
    TABLE_AXIS_X = "Consumption of silicon steels(m)"
    TABLE_AXIS_TIME = "Cost time(second)"
    SELECT_POINTS_NUM = 15

    GA_AXIS_X = "iterations"
    GA_AXIS_Y = "Unhomogeneity(ppm)"
