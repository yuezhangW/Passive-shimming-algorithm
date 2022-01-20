from src.Config import Config
from src.NMR import NMR
from src.Dataloader import Dataloader
from src.SensitivityCoefficientMatrix import SensitivityCoefficientMatrix
from src.Solver import Solver

from scipy.optimize import linprog
import os, time
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


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

        senMat.getSenMat()
        startTime = time.time()
        c, A, b, bounds = solver.getOptConfig(
            method=3,
            epsilon=0,
            kappa=Config.SHIMMING_KAPPA
        )
        resNew = linprog(c, A_ub=A, b_ub=b, bounds=bounds, method="highs")
        endTime = time.time()
        hom = resNew.x[-1] / 0.5 / 1e-6
        print("第" + str(i) + "组数据在Bt为定值情况下(Bt=B_avg)可匀场的最低不均匀度为：" + str(hom))
        print("第" + str(i) + "组耗时:"+str(endTime - startTime)+"s")


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
        senMat.getSenMat()

        startTime = time.time()
        res = solver.getBaseRes(Config.SHIMMING_TARGET_INHOMOGENEITY_MIN, Config.SHIMMING_TARGET_INHOMOGENEITY_MAX, Config.SHIMMING_STEP, Config.SHIMMING_KAPPA, 2)
        endTime = time.time()
        print("第" + str(i) + "组耗时:" + str(endTime - startTime) + "s")

# 4：基于L1范数最小二乘正则化算法的无源匀场算法
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

        res = solver.getNonlinearResult(epsilon=20)

# 5：Bt为定值时，计算不均匀度和最小匀场片厚度
def test5():
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

        epsilons = np.linspace(10, 100, 90, endpoint=False)
        dataList = []
        for epsilon in epsilons:
            startTime = time.time()
            c, A, b, bounds = solver.getOptConfig(
                method=1,
                epsilon=epsilon * 0.5 * 1e-6,
                kappa=Config.SHIMMING_KAPPA
            )
            res = linprog(c, A_ub=A, b_ub=b, bounds=bounds, method="highs")
            if res.success:
                # 计算总厚度
                thicknessTotal = np.sum(res.x) * nmr.thickness[0]
                # 计算不均匀度
                magneticFieldDataNow = magneticFieldData.getReshapeData() + senMat.getSenMat().T * nmr.thickness[0] @ res.x
                hom = (np.max(magneticFieldDataNow) - np.min(magneticFieldDataNow)) / magneticFieldData.bt * 1e6
                endTime = time.time()  # 记录程序结束运行时间

                dataList.append([hom, thicknessTotal, endTime - startTime])
        df = pd.DataFrame(dataList, columns=[Config.TABLE_AXIS_HOM, Config.TABLE_AXIS_X, Config.TABLE_AXIS_TIME])
        folderPath = os.path.join(os.path.join(os.path.dirname(os.getcwd()), 'resources'), 'ExcelFiles')
        filePath = os.path.join(folderPath, magneticFieldData.filename + "_hom&x_FTMF.xlsx")
        df.to_excel(filePath)
        print("完成"+magneticFieldData.filename+"的excel记录，单次运行时间均值为：" + str(df[Config.TABLE_AXIS_TIME].mean()) + "s")

        x, y = df[Config.TABLE_AXIS_HOM], df[Config.TABLE_AXIS_X]
        plt.xlabel(Config.TABLE_AXIS_HOM)
        plt.ylabel(Config.TABLE_AXIS_X)
        plt.plot(x, y, 'ob-')
        plt.show()

# 6：Bt为变量时，计算不均匀度和最小匀场片厚度
def test6():
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

        epsilons = np.linspace(10, 100, 90, endpoint=False)
        dataList = []
        for epsilon in epsilons:
            c, A, b, bounds = solver.getOptConfig(
                method=2,
                epsilon=epsilon * 0.5 * 1e-6,
                kappa=Config.SHIMMING_KAPPA
            )
            res = linprog(c, A_ub=A, b_ub=b, bounds=bounds, method="highs")
            if res.success:
                # 计算总厚度
                x = res.x[:-1]
                thicknessTotal = np.sum(x) * nmr.thickness[0]
                # 计算不均匀度
                magneticFieldDataNow = magneticFieldData.getReshapeData() + senMat.getSenMat().T * nmr.thickness[0] @ x
                hom = (np.max(magneticFieldDataNow) - np.min(magneticFieldDataNow)) / magneticFieldData.bt * 1e6
                dataList.append([hom, thicknessTotal])
        df = pd.DataFrame(dataList, columns=[Config.TABLE_AXIS_HOM, Config.TABLE_AXIS_X])
        folderPath = os.path.join(os.path.join(os.path.dirname(os.getcwd()), 'resources'), 'ExcelFiles')
        filePath = os.path.join(folderPath, magneticFieldData.filename + "_hom&x_OTMF.xlsx")
        df.to_excel(filePath)
        print("完成"+magneticFieldData.filename+"的excel记录")

        x, y = df[Config.TABLE_AXIS_HOM], df[Config.TABLE_AXIS_X]
        plt.xlabel(Config.TABLE_AXIS_HOM)
        plt.ylabel(Config.TABLE_AXIS_X)
        plt.plot(x, y, 'ob-')
        plt.show()


# 7：画图
def test7():
    filenames = [["6009map1.txt_hom&x_FTMF.xlsx", "6009map1.txt_hom&x_OTMF.xlsx"],
                 ["9006map1.txt_hom&x_FTMF.xlsx", "9006map1.txt_hom&x_OTMF.xlsx"]]
    imageFilenames = ['6009map1_hom&x_linear.png', '9006map1_hom&x_linear.png']
    for i, filename in enumerate(filenames):
        folderPath = os.path.join(os.path.join(os.path.dirname(os.getcwd()), 'resources'), 'ExcelFiles')
        filePath = os.path.join(folderPath, filename[0])
        dfFTMF = pd.read_excel(filePath)

        xFTMF, yFTMF = dfFTMF[Config.TABLE_AXIS_HOM], dfFTMF[Config.TABLE_AXIS_X]
        xFTMF, yFTMF = xFTMF[np.ceil(np.linspace(0, len(xFTMF)-1, Config.SELECT_POINTS_NUM))], yFTMF[np.ceil(np.linspace(0, len(yFTMF)-1, Config.SELECT_POINTS_NUM))]

        filePath = os.path.join(folderPath, filename[1])
        dfOTMF = pd.read_excel(filePath)

        xOTMF, yOTMF = dfOTMF[Config.TABLE_AXIS_HOM], dfOTMF[Config.TABLE_AXIS_X]
        xOTMF, yOTMF = xOTMF[np.ceil(np.linspace(0, len(xOTMF)-1, Config.SELECT_POINTS_NUM))], yOTMF[
            np.ceil(np.linspace(0, len(yOTMF)-1, Config.SELECT_POINTS_NUM))]

        plt.xlabel(Config.TABLE_AXIS_HOM)
        plt.ylabel(Config.TABLE_AXIS_X)
        plt.plot(xFTMF, yFTMF, 'ob-', xOTMF, yOTMF, 'or-')
        plt.legend(['FTMF(Bt=Bavg)', 'OTMF'])


        folderPath = os.path.join(os.path.join(os.path.dirname(os.getcwd()), 'resources'), 'ImageFiles')
        filePath = os.path.join(folderPath, imageFilenames[i])
        plt.savefig(filePath)
        plt.show()


if __name__ == '__main__':
    test1()
    test3()
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
