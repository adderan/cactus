from toil.job import Job

class WritePermanentFile(Job):
    def __init__(self, fileID, permanentPath):
        Job.__init__(self)
        self.fileID = fileID
        self.permanentPath = permanentPath
    def run(self, fileStore):
        fileStore.readGlobalFile(self.fileID, localFilePath=self.permanentPath)
        assert len(open(self.permanentPath).readlines()) > 0
    
