class UG27params:
    def __init__(self,
                 inside_diameter: float,
                 internal_pressure: float,
                 allowable_stress: float,
                 joint_efficiency_long: float,
                 joint_efficiency_circ: float,
                 thickness: float,
                 ):
        self.insideDiameter = inside_diameter
        self.internalPressure = internal_pressure
        self.allowableStress = allowable_stress
        self.jointEfficiencyLong = joint_efficiency_long
        self.jointEfficiencyCirc = joint_efficiency_circ
        self.thickness = thickness


class UG27Calcs:
    def __init__(self, UG27_params: UG27params):
        self.R = UG27_params.insideDiameter * 0.5
        self.P = UG27_params.internalPressure
        self.S = UG27_params.allowableStress
        self.E_long = UG27_params.jointEfficiencyLong
        self.E_circ = UG27_params.jointEfficiencyCirc
        self.t = UG27_params.thickness

    def design_thickness(self):
        return max(self.min_thk_circ_stress(), self.min_thk_longitudinal_stress())

    def min_thk_circ_stress(self):
        return self.P * self.R / (self.S * self.E_long - 0.6 * self.P)

    def min_thk_longitudinal_stress(self):
        return self.P * self.R / (2 * self.S * self.E_long + 0.4 * self.P)

    def max_pressure_longitudinal_stress(self):
        return 2 * self.S * self.E_circ * self.t / (self.R - 0.4 * self.t)

    def max_pressure_circ_stress(self):
        return self.S * self.E_long * self.t / (self.R + 0.6 * self.t)

    def max_pressure(self):
        return min(self.max_pressure_circ_stress(), self.max_pressure_longitudinal_stress())


if __name__ == '__main__':
    params = UG27params(
        inside_diameter=24.0,
        internal_pressure=100.0,
        allowable_stress=20000.0,
        joint_efficiency_long=1.0,
        joint_efficiency_circ=1.0,
        thickness=1.0
    )

    cyl = UG27Calcs(params)
    print(cyl.design_thickness())
    print(cyl.max_pressure())
