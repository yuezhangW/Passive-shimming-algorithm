import numpy as np
from src.Config import Config
import os


class Dataloader:
    def __init__(self, filename, r):
        self.folderPath = os.path.join(os.path.join(os.path.dirname(os.getcwd()), 'resources'), 'DamFiles')
        self.filePath = os.path.join(self.folderPath, filename)

        # DSV区域半径(m)
        self.r = r

        # 传感器的部分参数(暂无用，仅做保存)
        self.measurementParameters = {
            "Pr": 0,
            "NCY": 0,
            "MCF": 0,
            "MDA": 0,
            "Period": 0,
        }
        # 裸磁场矩阵数据（MHz格式）
        self.dataMat = None
        # meanValue比值矩阵 （T/Mhz）
        self.meanValueMat = None
        # 通过解析Dam文件，给measurementParameters和dataMat赋值
        self.readDamFile()
        # 目标磁场强度，这里用B_avg替代
        self.bt = np.mean(self.getMagneticFieldMatrix())

    # 返回磁场强度矩阵数据 （phi方向，theta方向个数）
    def getSize(self):
        return [len(self.getMagneticFieldMatrix()), len(self.getMagneticFieldMatrix()[0])]

    # 返回磁场强度矩阵的一维形式
    def getReshapeData(self):
        return self.getMagneticFieldMatrix().reshape((self.getSize()[0] * self.getSize()[1],), order='C')  # order=C表示优先行，'F'优先列

    # 返回磁场强度矩阵
    def getMagneticFieldMatrix(self):
        res = np.array(self.dataMat, np.float64)
        for i, item in enumerate(res):
            res[i] = res[i] * self.meanValueMat[i]
        return res

    # 返回该组磁场数据的不均匀度(ppm)
    def getUniformity(self):
        magneticFieldMatrix = self.getMagneticFieldMatrix()
        return ((np.max(magneticFieldMatrix) - np.min(magneticFieldMatrix)) / np.mean(magneticFieldMatrix)) * 1e6

    # 将dam文件读取并存储到内存数组中,并初始化各项参数
    def readDamFile(self):
        try:
            with open(self.filePath, "r") as f:
                # step1: 将所有字符串以二维list形式存储到dataTol中
                dataTol, dataEach, startRecord = [], [], False
                for line in f.readlines():
                    if line.find("Measurement parameters") != -1:
                        startRecord = True
                        if dataEach:
                            dataTol.append(dataEach)
                        dataEach = []
                    if startRecord:
                        dataEach.append(line.strip('\n'))  # 去掉列表中每一个元素的换行符
                if dataEach:
                    dataTol.append(dataEach)

                # step2: 将传感器采集的参数（MeasurementParameters）保存
                sensorParam = dataTol[0][1].split('\t')
                self.measurementParameters['Pr'] = int(sensorParam[0])
                self.measurementParameters['NCY'] = int(sensorParam[1])
                self.measurementParameters['MCF'] = int(sensorParam[2])
                self.measurementParameters['MDA'] = int(sensorParam[3])
                self.measurementParameters['Period'] = int(sensorParam[4])

                # step3: 获取每轮采集数据起始index
                startIndex, meanValueIndex = [], []
                for i, temp in enumerate(dataTol):
                    for j, line in enumerate(temp):
                        if line.find("Valid cycles") != -1:
                            startIndex.append([i, j + 2])
                        if line.find("Mean value") != -1:
                            meanValueIndex.append([i, j + 1])

                # step4: 解析数据
                self.dataMat = [[i for i in range(self.measurementParameters['Pr'])] for j in range(len(dataTol))]
                self.meanValueMat = []
                for i, index in enumerate(startIndex):
                    dataEach = dataTol[index[0]][index[1]:index[1]+self.measurementParameters['Pr']]
                    meanValue = float(dataTol[index[0]][meanValueIndex[i][1] + 1].split('\t')[0]) / float(dataTol[index[0]][meanValueIndex[i][1]].split('\t')[0])
                    self.meanValueMat.append(meanValue)
                    for j, sensorEach in enumerate(dataEach):
                        sensorEach = [float(item) for item in sensorEach.split('\t')]
                        self.dataMat[i][j] = sensorEach[1]
        except IOError as e:
            print(e, "Error:读取文件出错")




if __name__ == '__main__':
    a = Dataloader(Config.FILES_NAMES[1])
    print(len(a.dataMat), len(a.dataMat[0]))
    a.getMagneticFieldMatrix()
    print(a.getUniformity())
    print(a.getSize())

    #
    # file_name = '6009map1.txt'
    # folder_dir = os.path.join(os.getcwd(), 'DamFiles')
    # file_dir = os.path.join(folder_dir, file_name)
    #
    # start, end = 0, 0
    # ans = []
    # fin_data = []
    # with open(file_dir, "r") as f:
    #     for line in f.readlines():
    #         line = line.strip('\n')  # 去掉列表中每一个元素的换行符
    #
    #         if line.find('Measurement parameters') != -1:
    #             start = 1
    #             end = 0
    #             count = 0
    #
    #         if start:
    #             ans.append(line)
    #             count += 1
    #             if count == 47:
    #                 line_hz = ans[3]
    #                 hz = float(line_hz.split('\t')[0])
    #
    #                 line_Bz = ans[4]
    #                 Bz = float(line_Bz.split('\t')[0])
    #
    #                 bare_data_lines = ans[15:15+32]
    #
    #                 for i, bare_data_line in enumerate(bare_data_lines):
    #                     bare_data = float(bare_data_line.split('\t')[1])
    #                     bare_data = bare_data / hz * Bz
    #                     fin_data.append(bare_data)
    #                 #print(bare_data_lines)
    #
    #                 start = 0
    #                 end = 1
    #                 count = 0
    #
    #                 ans = []
    #
    # file_name = 'dam1.txt'
    # file_dir = os.path.join(folder_dir, file_name)
    # with open(file_dir, "w") as f:
    #     for i in fin_data:
    #         f.write(str(i) + " ")
    #
    # print(fin_data, len(fin_data))
