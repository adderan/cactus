from toil.job import Job

class WritePermanentFile(Job):
    def __init__(self, fileID, permanentPath):
        Job.__init__(self)
        self.fileID = fileID
        self.permanentPath = permanentPath
    def run(self, fileStore):
        assert fileStore.globalFileExists(self.fileID)
        fileStore.readGlobalFile(self.fileID, localFilePath=self.permanentPath)
        #assert len(open(self.permanentPath).readlines()) > 0

class WriteToFilestore(Job):
    def __init__(self, fileID, filePath):
        Job.__init__(self)
        self.filePath = filePath
        self.fileID = fileID
    def run(self, fileStore):
        fileStore.updateGlobalFile(self.fileID, self.filePath)
