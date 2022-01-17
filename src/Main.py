from src.Config import Config
from src.NMR import NMR
from src.Dataloader import Dataloader
from src.SensitivityCoefficientMatrix import SensitivityCoefficientMatrix
from src.Solver import Solver

from scipy.optimize import linprog
import numpy as np


# 1：计算裸磁场数据的匀场前不均匀度
def test1():
    for i, filename in enumerate(Config.FILES_NAMES):
        magneticFieldData = Dataloader(
            filename=filename,
            r=Config.DSV_R
        )
        print("第" + str(i) + "组数据匀场前不均匀度：" + str(magneticFieldData.getUniformity()))


# 2：Bt为定值时，可匀场的最低不均匀度
def test2():
    for i, filename in enumerate(Config.FILES_NAMES):
        # 配置NMR各项参数
        nmr = NMR(
            size=Config.NMR_SIZE,
            thickness=Config.NMR_THICKNESS,
            mz=Config.NMR_MZ,
            pieceArea=Config.NMR_PIECE_AREA,
            z=Config.NMR_Z,
            r=Config.NMR_R
        )

        # 配置DSV采样点数据
        magneticFieldData = Dataloader(
            filename=filename,
            r=Config.DSV_R
        )

        # 获取灵敏度系数矩阵
        senMat = SensitivityCoefficientMatrix(
            dataloader=magneticFieldData,
            NMR=nmr
        )

        solver = Solver(
            senMat=senMat,
            magneticField=magneticFieldData,
            nmr=nmr
        )

        c, A, b, bounds = solver.getOptConfig(
            method=3,
            epsilon=0,
            kappa=Config.SHIMMING_KAPPA
        )
        resNew = linprog(c, A_ub=A, b_ub=b, bounds=bounds, method="highs")
        hom = resNew.x[-1] / 0.5 / 1e-6
        print("第" + str(i) + "组数据在Bt为定值情况下(Bt=B_avg)可匀场的最低不均匀度为：" + str(hom))


# 3：Bt为变量时，可匀场的最低不均匀度
def test3():
    for i, filename in enumerate(Config.FILES_NAMES):
        # 配置NMR各项参数
        nmr = NMR(
            size=Config.NMR_SIZE,
            thickness=Config.NMR_THICKNESS,
            mz=Config.NMR_MZ,
            pieceArea=Config.NMR_PIECE_AREA,
            z=Config.NMR_Z,
            r=Config.NMR_R
        )

        # 配置DSV采样点数据
        magneticFieldData = Dataloader(
            filename=filename,
            r=Config.DSV_R
        )

        # 获取灵敏度系数矩阵
        senMat = SensitivityCoefficientMatrix(
            dataloader=magneticFieldData,
            NMR=nmr
        )

        solver = Solver(
            senMat=senMat,
            magneticField=magneticFieldData,
            nmr=nmr
        )

        res = solver.getBaseRes(Config.SHIMMING_TARGET_INHOMOGENEITY_MIN, Config.SHIMMING_TARGET_INHOMOGENEITY_MAX, Config.SHIMMING_STEP, Config.SHIMMING_KAPPA, 2)


# 4：非线性规划
def test4():
    for i, filename in enumerate(Config.FILES_NAMES):
        # 配置NMR各项参数
        nmr = NMR(
            size=Config.NMR_SIZE,
            thickness=Config.NMR_THICKNESS,
            mz=Config.NMR_MZ,
            pieceArea=Config.NMR_PIECE_AREA,
            z=Config.NMR_Z,
            r=Config.NMR_R
        )

        # 配置DSV采样点数据
        magneticFieldData = Dataloader(
            filename=filename,
            r=Config.DSV_R
        )

        # 获取灵敏度系数矩阵
        senMat = SensitivityCoefficientMatrix(
            dataloader=magneticFieldData,
            NMR=nmr
        )

        solver = Solver(
            senMat=senMat,
            magneticField=magneticFieldData,
            nmr=nmr
        )

        #res = solver.getNonlinearResult()
        solver.getLeastSquare()



if __name__ == '__main__':
    test1()
    test4()
    # # 配置NMR各项参数
    # nmr = NMR(
    #     size=Config.NMR_SIZE,
    #     thickness=Config.NMR_THICKNESS,
    #     mz=Config.NMR_MZ,
    #     pieceArea=Config.NMR_PIECE_AREA,
    #     z=Config.NMR_Z,
    #     r=Config.NMR_R
    # )
    #
    # # 配置DSV采样点数据
    # magneticFieldData = Dataloader(
    #     filename=Config.FILES_NAMES[Config.FILES_INDEX],
    #     r=Config.DSV_R,
    #     bt=Config.DSV_BT
    # )
    #
    # # 获取灵敏度系数矩阵
    # senMat = SensitivityCoefficientMatrix(
    #     dataloader=magneticFieldData,
    #     NMR=nmr
    # )
    #
    # solver = Solver(
    #     senMat=senMat,
    #     magneticField=magneticFieldData,
    #     nmr=nmr
    # )
    # print(magneticFieldData.getReshapeData().shape)
    # # print(senMat.getSenMat().shape, magneticFieldData.getMagneticFieldMatrix().shape)
    # print(solver.getBaseRes(Config.SHIMMING_TARGET_INHOMOGENEITY_MIN, Config.SHIMMING_TARGET_INHOMOGENEITY_MAX, Config.SHIMMING_STEP, Config.SHIMMING_KAPPA))
