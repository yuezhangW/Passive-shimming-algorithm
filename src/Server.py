import json

from fastapi import FastAPI
from pydantic import BaseModel
import numpy as np
import uvicorn


class MatrixConfig(BaseModel):
    nmrNumOfPhi: int
    nmrNumOfZ: int
    nmrLengthOfZ: float
    nmrLengthOfR: float
    nmrThicknessMin: float
    nmrThicknessMax: float
    nmrMz: float
    nmrPieceArea: float

    dsvNumOfPhi: int
    dsvNumOfTheta: int
    dsvLengthOfR: float


app = FastAPI()


@app.post("/task/matrix")
async def createMatrix(matrixConfig: MatrixConfig):
    nmrNumOfPhi = matrixConfig.nmrNumOfPhi
    nmrNumOfZ = matrixConfig.nmrNumOfZ
    nmrLengthOfZ = matrixConfig.nmrLengthOfZ
    nmrLengthOfR = matrixConfig.nmrLengthOfR
    nmrThicknessMin = matrixConfig.nmrThicknessMin
    nmrThicknessMax = matrixConfig.nmrThicknessMax
    nmrMz = matrixConfig.nmrMz
    nmrPieceArea = matrixConfig.nmrPieceArea

    dsvNumOfPhi = matrixConfig.dsvNumOfPhi
    dsvNumOfTheta = matrixConfig.dsvNumOfTheta
    dsvLengthOfR = matrixConfig.dsvLengthOfR

    qPoints, pPoints = getQPoints(int(nmrNumOfPhi), int(nmrNumOfZ), np.float(nmrLengthOfR),
                                  np.float(nmrLengthOfZ)), getPPoints(int(dsvNumOfPhi), int(dsvNumOfTheta),
                                                                      np.float(dsvLengthOfR))

    senMat = np.zeros((qPoints.shape[0], pPoints.shape[0]), np.float64)
    print(senMat.shape)
    for i, qPoint in enumerate(qPoints):
        for j, pPoint in enumerate(pPoints):
            senMat[i, j] = nmrPieceArea * getBzP2P(nmrMz, qPoint, pPoint)
    print(senMat.tolist())
    return {"senMat": senMat.tolist()}


"""
计算点对点的Bz
"""

def getBzP2P(mz, qPoint, pPoint):
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

def getQPoints(nmrNumOfPhi, nmrNumOfZ, nmrLengthOfR, nmrLengthOfZ):
    pieceNums = nmrNumOfPhi * nmrNumOfZ
    qPoints = np.zeros((pieceNums, 3), np.float64)

    phiList = np.linspace(0, 2 * np.pi, nmrNumOfPhi, endpoint=False)
    zList = np.linspace(-nmrLengthOfZ / 2, nmrLengthOfZ / 2, nmrNumOfZ, endpoint=False) + (nmrLengthOfZ / nmrNumOfZ / 2)

    for i, phi in enumerate(phiList):
        for j, z in enumerate(zList):
            qPoints[i * len(zList) + j, 0] = phi
            qPoints[i * len(zList) + j, 1] = np.arctan2(nmrLengthOfR, z)
            qPoints[i * len(zList) + j, 2] = np.sqrt(z ** 2 + nmrLengthOfR ** 2)

    return qPoints

"""
使用默认均匀尺寸生成P的坐标(极坐标)
phiNums: phi方向piece个数
thetaNums: theta方向piece个数
r: DSV区域内径(半径 m)
"""

def getPPoints(dsvNumOfPhi, dsvNumOfTheta, dsvLengthOfR):
    sampleNums = dsvNumOfPhi * dsvNumOfTheta
    pPoints = np.zeros((sampleNums, 3), np.float64)

    phiList = np.linspace(0, 2 * np.pi, dsvNumOfPhi, endpoint=False)
    thetaList = np.linspace(0, np.pi, dsvNumOfTheta, endpoint=False) + (np.pi / dsvNumOfTheta / 2)

    for i, phi in enumerate(phiList):
        for j, theta in enumerate(thetaList):
            pPoints[i * len(thetaList) + j, 0] = phi
            pPoints[i * len(thetaList) + j, 1] = theta
            pPoints[i * len(thetaList) + j, 2] = dsvLengthOfR

    return pPoints


if __name__ == '__main__':
    uvicorn.run("Server:app", host="127.0.0.1", port=8080)
