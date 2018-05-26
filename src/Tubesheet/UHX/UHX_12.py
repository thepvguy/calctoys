from ._UHX_common import Configuration as Configuration
from ._UHX_common import AttachmentType as AttachmentType

class UHX_12_calcs:
    def __init__(self, params):
        self.UHX11 = params.UHX11
        self.config = params.config
        self.tubesideAttachmentType = params.tubesideAttachmentType
        self.shellsideAttachmentType = params.shellsideAttachmentType
        self.D_o = params.UHX11.D_o()
        self.D_s = params.D_s
        self.G_s = params.G_s

    def rho_s(self):
        """

        :return: Ratio of shell ID to OTL
        """
        if self.config in [Configuration.a, Configuration.b, Configuration.c]:

            ratio = self.D_s / self.D_o
        elif self.config in [Configuration.d, Configuration.e, Configuration.f]:
            ratio = self.G_s / self.D_o

        else:
            raise ValueError("Invalid configuration <<%r>> specified" % self.config)

        return ratio