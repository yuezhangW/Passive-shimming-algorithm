import os
import time

import numpy as np
import pandas as pd
import cvxpy as cp
import scipy.optimize
from scipy.optimize import linprog
from sko.GA import GA
from sko.PSO import PSO
from sko.DE import DE
from sko.SA import SA
from sko.IA import IA_TSP
from sko.AFSA import AFSA
from matplotlib import pyplot as plt

from src.Config import Config



class Solver:
    def __init__(self, senMat, magneticField, nmr):
        self.senMat = senMat
        self.magneticField = magneticField
        self.nmr = nmr

    """
    步长迭代的方式搜索理论最低下界
    """

    def getBaseRes(self, stepNum, kappa, method):
        targetHomList = np.linspace(0, self.magneticField.getUniformity(), stepNum)

        # 利用二分法进行搜索
        def binarySearch(arr, left, right):
            i = 0
            while left < right:
                i += 1
                mid = int((left + right) / 2)
                c, A, b, bounds = self.getOptConfig(
                    method=method,
                    epsilon=arr[mid] * 0.5 * 1e-6,
                    kappa=kappa
                )
                res = linprog(c, A_ub=A, b_ub=b, bounds=bounds, method="highs")
                print("第" + str(i) + "轮搜索，不均匀度为：" + str(arr[mid]) + "是否存在可行解？" + str(res.success))
                if res.success:
                    right = mid
                else:
                    left = mid + 1
            return left

        idx = binarySearch(targetHomList, 0, len(targetHomList))
        c, A, b, bounds = self.getOptConfig(
            method=method,
            epsilon=targetHomList[idx] * 0.5 * 1e-6,
            kappa=kappa
        )
        res = linprog(c, A_ub=A, b_ub=b, bounds=bounds, method="highs")
        print("最低不均匀度为：" + str(targetHomList[idx]) + str(res.success))
        return res

    """
    返回线性规划包需要用到的各项参数
    """

    def getOptConfig(self, method, epsilon, kappa):
        senMat, magneticFieldMat = self.senMat.getSenMat(), self.magneticField.getReshapeData()
        thickness, thicknessMax = self.nmr.thickness[0], self.nmr.thickness[1]
        bt = self.magneticField.bt

        # bt和T都为定值
        if method == 1:
            c = np.ones((senMat.shape[0]), np.float64)
            A = np.hstack((senMat, -senMat)).T
            b = np.hstack((-magneticFieldMat + (1 + epsilon) * bt, magneticFieldMat - (1 - epsilon) * bt))

            A = A * thickness
            lower = np.zeros((senMat.shape[0]), np.float64)
            upper = np.ones((senMat.shape[0]), np.float64) * thicknessMax / thickness
            bounds = [(lower[i], upper[i]) for i in range(lower.shape[0])]

        # bt为变量T为定值
        elif method == 2:
            c = np.ones((senMat.shape[0] + 1), np.float64)
            c[-1] = 0
            A = np.hstack((
                np.vstack((senMat, -(1 + epsilon) * np.ones((1, senMat.shape[1])))),
                -np.vstack((senMat, -(1 - epsilon) * np.ones((1, senMat.shape[1]))))
            )).T
            b = np.hstack((-magneticFieldMat, magneticFieldMat))

            A[:, :-1] = A[:, :-1] * thickness
            lower = np.zeros((senMat.shape[0] + 1), np.float64)
            lower[-1] = bt * (1 - kappa)
            upper = np.ones((senMat.shape[0] + 1), np.float64) * thicknessMax / thickness
            upper[-1] = bt * (1 + kappa)
            bounds = [(lower[i], upper[i]) for i in range(lower.shape[0])]

        # bt为定值T为变量
        elif method == 3:
            c = np.zeros((senMat.shape[0] + 1), np.float64)
            c[-1] = 1
            A = np.hstack((
                np.vstack((senMat, -bt * np.ones((1, senMat.shape[1])))),
                -np.vstack((senMat, bt * np.ones((1, senMat.shape[1]))))
            )).T
            b = np.hstack((-magneticFieldMat + bt, magneticFieldMat - bt))

            A[:, :-1] = A[:, :-1] * thickness
            lower = np.zeros((senMat.shape[0] + 1), np.float64)
            upper = np.ones((senMat.shape[0] + 1), np.float64) * thicknessMax / thickness
            upper[-1] = 1
            bounds = [(lower[i], upper[i]) for i in range(lower.shape[0])]

        return c, A, b, bounds

    """
    最小二乘法minimize（非线性）
    """

    def getNonlinearResult(self, epsilon, l1, l2):
        def func(args):
            A, Bm, Bt = args
            v = lambda x: np.sum(np.square(A @ x + Bm - Bt)) + l1 * np.sum(np.abs(A @ x + Bm - Bt)) + l2 * np.sum(np.abs(x))
            return v

        def cons(args):
            A, Bm, Bt, epsilon = args
            # ineq:大于等于0
            con = [
                {'type': 'ineq', 'fun': lambda x: -A @ x + (1 + epsilon) * Bt - Bm}
            ] + [
                {'type': 'ineq', 'fun': lambda x: A @ x - (1 - epsilon) * Bt + Bm}
            ]
            return con

        senMat, magneticFieldMat = self.senMat.getSenMat(), self.magneticField.getReshapeData()
        thickness, thicknessMax = self.nmr.thickness[0], self.nmr.thickness[1]
        bt = self.magneticField.bt

        senMat = senMat * thickness
        con = cons((senMat.T,
            magneticFieldMat,
            bt,
            epsilon * 0.5 * 1e-6
        ))

        args = ((
            senMat.T,
            magneticFieldMat,
            bt
        ))

        lower = np.zeros((senMat.shape[0]), np.float64)
        upper = np.ones((senMat.shape[0]), np.float64) * thicknessMax / thickness
        bounds = [(lower[i], upper[i]) for i in range(lower.shape[0])]
        x0 = np.zeros((senMat.shape[0],), np.float64)
        res = scipy.optimize.minimize(func(args), x0, constraints=con, bounds=bounds)
        print("result:" + str(res.success))
        hom = senMat.T @ res.x + magneticFieldMat
        homValue = (np.max(hom) - np.min(hom)) / bt * 1e6
        print("不均匀度：" + str(homValue))

        thicknessTotal = thickness * np.sum(res.x)
        print("匀场片总厚度：" + str(thicknessTotal))

    """
    直接使用启发式算法做整数规划
    """
    def integerProgrammingBySKO(self):
        thickness, thicknessMax = self.nmr.thickness[0], self.nmr.thickness[1]
        A = self.senMat.getSenMat().T * thickness
        Bm = self.magneticField.getReshapeData()
        bt = self.magneticField.bt
        dim = A.shape[1]

        def schaffer(p):
            y = A @ p + Bm
            y = (np.max(y) - np.min(y)) / bt * 1e6
            return y

        lb = [0 for i in range(dim)]
        ub = [(thicknessMax / thickness) for i in range(dim)]
        precision = 1

        ga = GA(
            func=schaffer,
            n_dim=dim,
            size_pop=50,
            max_iter=1000,
            prob_mut=0.001,
            lb=lb,
            ub=ub,
            precision=precision
        )
        best_x, best_y = ga.run()
        print('best_x:', best_x, '\n', 'best_y:', best_y)

        hom = self.senMat.getSenMat().T * thickness @ best_x + self.magneticField.getReshapeData()
        homValue = (np.max(hom) - np.min(hom)) / self.magneticField.bt * 1e6
        print("不均匀度：" + str(homValue))

        thicknessTotal = thickness * np.sum(best_x)
        print("匀场片总厚度：" + str(thicknessTotal))


        Y_history = pd.DataFrame(ga.all_history_Y)
        fig, ax = plt.subplots(2, 1)
        ax[0].plot(Y_history.index, Y_history.values, '.', color='red')
        Y_history.min(axis=1).cummin().plot(kind='line')
        ax[0].set_ylabel(Config.GA_AXIS_Y)
        ax[1].set_xlabel(Config.GA_AXIS_X)
        ax[1].set_ylabel(Config.GA_AXIS_Y)

        folderPath = os.path.join(os.path.join(os.path.dirname(os.getcwd()), 'resources'), 'ImageFiles')
        filePath = os.path.join(folderPath, self.magneticField.filename.split('.')[0]+"_GA_1-Stage.png")
        plt.savefig(filePath)
        plt.show()

    def mixedIntegerAndNonProgramming(self, inteLen):
        nonInteResult = self.getBaseRes(Config.SHIMMING_STEP, Config.SHIMMING_KAPPA, 2)
        nonInteX, bt = nonInteResult.x[:-1], nonInteResult.x[-1]
        thickness, thicknessMax = self.nmr.thickness[0], self.nmr.thickness[1]

        hom = self.senMat.getSenMat().T * thickness @ nonInteX + self.magneticField.getReshapeData()
        homValue = (np.max(hom) - np.min(hom)) / bt * 1e6
        print("全局最优的非整数线性规划解")
        print("不均匀度：" + str(homValue))
        thicknessTotal = thickness * np.sum(nonInteX)
        print("匀场片总厚度：" + str(thicknessTotal))

        print("如果直接对所有非整数解向上取整")
        ceilX = np.ceil(nonInteX)
        hom = self.senMat.getSenMat().T * thickness @ ceilX + self.magneticField.getReshapeData()
        homValue = (np.max(hom) - np.min(hom)) / bt * 1e6
        print("不均匀度：" + str(homValue))
        thicknessTotal = thickness * np.sum(ceilX)
        print("匀场片总厚度：" + str(thicknessTotal))

        print("如果直接对所有非整数解向下取整")
        floorX = np.floor(nonInteX)
        hom = self.senMat.getSenMat().T * thickness @ floorX + self.magneticField.getReshapeData()
        homValue = (np.max(hom) - np.min(hom)) / bt * 1e6
        print("不均匀度：" + str(homValue))
        thicknessTotal = thickness * np.sum(floorX)
        print("匀场片总厚度：" + str(thicknessTotal))

        A = self.senMat.getSenMat().T * thickness
        Bm = self.magneticField.getReshapeData()
        def schaffer(p):
            y = A @ p + Bm
            y = (np.max(y) - np.min(y)) / bt * 1e6
            return y

        def schafferInt(p):
            p = np.round(p)
            y = A @ p + Bm
            y = (np.max(y) - np.min(y)) / bt * 1e6
            return y

        lb = [np.max([0, np.floor(i - inteLen)]) for i in nonInteX]
        ub = [np.min([(thicknessMax / thickness), np.ceil(i + inteLen)]) for i in nonInteX]
        precision = 1
        dim = len(nonInteX)

        algorithms = [
            GA(
                func=schaffer, n_dim=dim, size_pop=50, max_iter=800, prob_mut=0.001, lb=lb, ub=ub, precision=precision
            ),
            DE(
                func=schafferInt, n_dim=dim, size_pop=50, max_iter=800, lb=lb, ub=ub
            ),
            PSO(
                func=schafferInt, n_dim=dim, pop=40, max_iter=150, lb=lb, ub=ub, w=0.8, c1=0.5, c2=0.5
            ),
            SA(
                func=schaffer, x0=[0 for i in range(dim)], lb=lb, ub=ub, L=300, max_stay_counter=150,  precision=precision
            ),
        ]
        names = ['遗传算法', '差分进化算法', '粒子群算法', '模拟退火算法']

        for i, algorithm in enumerate(algorithms):
            startTime = time.time()
            best_x, best_y = algorithm.run()
            print(names[i])

            hom = self.senMat.getSenMat().T * thickness @ best_x + self.magneticField.getReshapeData()
            homValue = (np.max(hom) - np.min(hom)) / bt * 1e6
            print("不均匀度：" + str(homValue))

            thicknessTotal = thickness * np.sum(best_x)
            print("匀场片总厚度：" + str(thicknessTotal))

            endTime = time.time()
            print("耗时：" + str(endTime - startTime))


        # Y_history = pd.DataFrame(ga.all_history_Y)
        # fig, ax = plt.subplots(2, 1)
        # ax[0].plot(Y_history.index, Y_history.values, '.', color='red')
        # Y_history.min(axis=1).cummin().plot(kind='line')
        # plt.show()