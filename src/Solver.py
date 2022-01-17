import numpy as np
import cvxpy as cp
import scipy.optimize
from scipy.optimize import linprog
from src.Config import Config


class Solver:
    def __init__(self, senMat, magneticField, nmr):
        self.senMat = senMat
        self.magneticField = magneticField
        self.nmr = nmr

    """
    步长迭代的方式搜索理论最低下界
    """

    def getBaseRes(self, targetHomMin, targetHomMax, stepNum, kappa, method):
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
                res = linprog(c, A_ub=A, b_ub=b, bounds=bounds,method="highs")
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
            c = np.zeros((senMat.shape[0] + 1), np.float64)
            c[-1] = 1
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
    def getLeastSquare(self):
        senMat, magneticFieldMat = self.senMat.getSenMat(), self.magneticField.getReshapeData()
        thickness, thicknessMax = self.nmr.thickness[0], self.nmr.thickness[1]
        bt = self.magneticField.bt

        A = senMat.T * thickness
        x = cp.Variable(senMat.shape[0])

        objective = cp.Minimize(cp.sum_squares(A @ x + magneticFieldMat - bt))
        epsilon = 100*0.5*1e-6
        constraints = [0 <= x, x <= (thicknessMax / thickness), A @ x - (1 + epsilon) * bt <= -magneticFieldMat, -A @ x + (1 - epsilon) * bt <= magneticFieldMat]
        prob = cp.Problem(objective, constraints)

        # The optimal objective value is returned by `prob.solve()`.
        result = prob.solve()
        print(result)
        # The optimal value for x is stored in `x.value`.
        print(x.value)
        # The optimal Lagrange multiplier for a constraint is stored in
        # `constraint.dual_value`.
        hom = A @ x.value + magneticFieldMat
        homValue = (np.max(hom) - np.min(hom)) / bt / 1e-6
        print(homValue)
        #print(constraints[0].dual_value)


    # 非线性规划（暂未生效）
    def getNonlinearResult(self):
        def func(args):
            A, Bm, Bt = args
            v = lambda x: np.sum(np.square(A @ x + Bm - Bt))
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
            25 * 0.5 * 1e-6
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
        res = scipy.optimize.minimize(func(args), x0, method='SLSQP', constraints=con, bounds=bounds)
        print(res)
        print("不均匀度：")
        hom = senMat.T @ res.x + magneticFieldMat
        print(hom)
        homValue = (np.max(hom) - np.min(hom)) / bt * 1e6
        print(homValue)

