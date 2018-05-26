import enum

class Configuration(enum.Enum):
    a = 1
    b = 2
    c = 3
    d = 4
    e = 5
    f = 6
    A = 7
    B = 8
    C = 9
    D = 10

class Condition(enum.Enum):
    NEW = 1
    CORRODED = 2

class LoadCaseType(enum.Enum):
    OPERATING = 2
    DESIGN = 1

LoadCaseNumber = [1,2,3,4]

class PitchType(enum.Enum):
    TRIANGLE = 1
    SQUARE = 2
    ROTATED_TRIANGLE = 3
    ROTATED_SQUARE = 4

class AttachmentType(enum.Enum):
    CYLINDER = 1
    HEMI_HEAD = 2

if __name__=="__main__":
    pass