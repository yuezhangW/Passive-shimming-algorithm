import numpy as np


# 灵敏性系数矩阵类，根据输入数据集计算灵敏度系数矩阵
class SensitivityCoefficientMatrix:
    def __init__(self, dataloader, NMR):
        self.dataloader = dataloader
        self.NMR = NMR
        self.data = None

    """
        返回灵敏度系数矩阵
        P:DSV区域采样点
        Q:shimming piece中心点
        """

    def getSenMat(self):
        if self.data is not None:
            return self.data
        qPoints, pPoints, mz = self.getQPoints(), self.getPPoints(), self.NMR.mz
        senMat = np.zeros((qPoints.shape[0], pPoints.shape[0]), np.float64)
        for i, qPoint in enumerate(qPoints):
            for j, pPoint in enumerate(pPoints):
                senMat[i, j] = self.NMR.pieceArea * self.getBzP2P(mz, qPoint, pPoint)
        self.data = senMat
        return senMat

    """
    计算点对点的Bz
    """

    def getBzP2P(self, mz, qPoint, pPoint):
        pAxes = [
            pPoint[2] * np.sin(pPoint[1]) * np.cos(pPoint[0]),
            pPoint[2] * np.sin(pPoint[1]) * np.sin(pPoint[0]),
            pPoint[2] * np.cos(pPoint[1])
        ]

        qAxes = [
            qPoint[2] * np.sin(qPoint[1]) * np.cos(qPoint[0]),
            qPoint[2] * np.sin(qPoint[1]) * np.sin(qPoint[0]),
            qPoint[2] * np.cos(qPoint[1])
        ]

        rSquare = (pAxes[0] - qAxes[0]) ** 2 + (pAxes[1] - qAxes[1]) ** 2
        zSquare = (pAxes[2] - qAxes[2]) ** 2

        bz = (mz / (4 * np.pi)) * (
                (3 * zSquare) / ((rSquare + zSquare) ** 2.5) - 1 / ((rSquare + zSquare) ** 1.5))
        return bz



    """
    使用默认均匀尺寸生成Q的坐标(极坐标)
    phiNums: phi方向piece个数
    zNums: z方向piece个数
    r: 腔体内径(半径 m)
    z: 腔体总长度(m)
    """

    def getQPoints(self):
        phiNums, zNums = self.NMR.size[0], self.NMR.size[1]
        r, z = self.NMR.r, self.NMR.z
        pieceNums = phiNums * zNums
        qPoints = np.zeros((pieceNums, 3), np.float64)

        phiList = np.linspace(0, 2 * np.pi, phiNums, endpoint=False)
        zList = np.linspace(-z / 2, z / 2, zNums, endpoint=False) + (z / zNums / 2)

        for i, phi in enumerate(phiList):
            for j, z in enumerate(zList):
                qPoints[i * len(zList) + j, 0] = phi
                qPoints[i * len(zList) + j, 1] = np.arctan2(r, z)
                qPoints[i * len(zList) + j, 2] = np.sqrt(z ** 2 + r ** 2)

        return qPoints

    """
    使用默认均匀尺寸生成P的坐标(极坐标)
    phiNums: phi方向piece个数
    thetaNums: theta方向piece个数
    r: DSV区域内径(半径 m)
    """

    def getPPoints(self):
        phiNums, thetaNums = self.dataloader.getSize()[0], self.dataloader.getSize()[1]
        r = self.dataloader.r
        sampleNums = phiNums * thetaNums
        pPoints = np.zeros((sampleNums, 3), np.float64)

        phiList = np.linspace(0, 2 * np.pi, phiNums, endpoint=False)
        thetaList = np.linspace(0, np.pi, thetaNums, endpoint=False) + (np.pi / thetaNums / 2)

        for i, phi in enumerate(phiList):
            for j, theta in enumerate(thetaList):
                pPoints[i * len(thetaList) + j, 0] = phi
                pPoints[i * len(thetaList) + j, 1] = theta
                pPoints[i * len(thetaList) + j, 2] = r

        return pPoints
