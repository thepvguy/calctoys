from . import common

class FlangeExtensionParams:
    def __init__(self, number_of_steps = 0):
        self.number_of_steps = number_of_steps


class FlangedExtension:
    def __init__(self, params):
        self.params = params

    def getResults(self):
        pass