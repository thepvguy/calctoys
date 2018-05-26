

class cylinderParams:
    def __init__(self,
                 inside_diameter: float,
                 internal_pressure: float,
                 allowable_stress: float,
                 joint_efficiency_long: float,
                 joint_efficiency_circ: float,
                 thickness: float,
                 internal_corrosion: float,
                 external_corrosion: float):

        self.insideDiameter = inside_diameter
        self.internalPressure = internal_pressure
        self.allowableStress = allowable_stress
        self.jointEfficiencyLong = joint_efficiency_long
        self.jointEfficiencyCirc = joint_efficiency_circ
        self.thickness = thickness
        self.internalCorrosion = internal_corrosion
        self.externalCorrosion = external_corrosion


class UG27Calcs:
    def __init__(self, cylinder_params: cylinderParams):
        self.R = cylinder_params.insideDiameter * 0.5
        self.P = cylinder_params.internalPressure
        self.S = cylinder_params.allowableStress
        self.E_long = cylinder_params.jointEfficiencyLong
        self.E_circ = cylinder_params.jointEfficiencyCirc
        self.t = cylinder_params.thickness
        self.internalCorrosion = cylinder_params.internalCorrosion
        self.externalCorrosion = cylinder_params.externalCorrosion

    def design_thickness(self):
        return max(self.min_thickness_circ(), self.min_thickness_long())

    def min_thickness_circ(self):
        return (self.P * (self.R + self.internalCorrosion) / (
            self.S * self.E_long - 0.6 * self.P)) + self.internalCorrosion + self.externalCorrosion

    def min_thickness_long(self):
        return (self.P * (self.R + self.internalCorrosion) / (
            2 * self.S * self.E_long + 0.4 * self.P)) + self.internalCorrosion + self.externalCorrosion

    def mawp_long(self):
        res = (2 * self.S * self.E_circ * (self.t - self.externalCorrosion - self.internalCorrosion) /
               ((self.R + self.internalCorrosion) - 0.4 * (
                    self.t - self.internalCorrosion - self.externalCorrosion)))
        return res

    def mawp_circ(self):
        res = (self.S * self.E_long * (self.t - self.externalCorrosion - self.internalCorrosion) /
               ((self.R + self.internalCorrosion) + 0.6 * (self.t - self.internalCorrosion - self.externalCorrosion)))
        return res

    def mawp(self):
        return min(self.mawp_circ(), self.mawp_long())


if __name__ == '__main__':

    params = cylinderParams(
        inside_diameter=24.0,
        internal_pressure=100.0,
        allowable_stress=20000.0,
        joint_efficiency_long=1.0,
        joint_efficiency_circ=1.0,
        thickness=1.0,
        internal_corrosion=0.125,
        external_corrosion=0.25
    )

    cyl = UG27Calcs(params)
    print(cyl.design_thickness())
    print(cyl.mawp())
