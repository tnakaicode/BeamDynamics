class boundaryCondition(object):
    def __init__(self, *args):
        self.bcType = args[0]
        self.values = args[1]

        # set something like dimension based on dimension of values?
        # set x,y,z and 0 or N?
    def getType(self):
        return self.bcType

    def getValues(self):
        return self.values
