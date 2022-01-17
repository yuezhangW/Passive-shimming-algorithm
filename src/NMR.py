class NMR:
    def __init__(self, size, thickness, mz, pieceArea, z, r):
        self.size = size  # phi方向匀场片个数和z方向匀场片个数 list格式
        self.thickness = thickness  # 匀场片最小厚度和空腔的最大厚度 list格式
        self.mz = mz
        self.pieceArea = pieceArea  # 匀场片的面积
        self.z = z  # NMR的z方向长度（单位：m）
        self.r = r  # NMR的内半径（单位：m）
