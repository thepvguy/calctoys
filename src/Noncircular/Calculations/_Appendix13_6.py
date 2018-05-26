

class MultiDiameterHole:
    def __init__(self, params):
        self.t = params.t
        self.p = params.p
        self.diameter_thickness_pairs = params.diameter_thickness_pairs

    def isSane(self):
        isSane = True
        floatFuzz = 0.000001
        DTPairThicknessSum = 0

        for pair in self.diameter_thickness_pairs:
            DTPairThicknessSum += pair[1]

        if DTPairThicknessSum - floatFuzz >= self.t:
            isSane = False

        return isSane

    def __b_n(self, dia):
        return self.p - dia

    def X_bar(self):
        """

        :return: Appendix 13-6 equation 6, effective metal area between holes
        """

        X_bar_num = 0
        X_bar_den = 0
        for step in range(len(self.diameter_thickness_pairs)):
            b_n = self.__b_n(self.diameter_thickness_pairs[step][0])
            T_n = self.diameter_thickness_pairs[step][1]
            X_bar_num += b_n * T_n * self.__effectiveStepThickness(step)
            X_bar_den += b_n * T_n

        return X_bar_num / X_bar_den

    def __effectiveStepThickness(self, index):
        effectiveThickness = self.diameter_thickness_pairs[index][1] * 0.5

        if index + 1 != len(self.diameter_thickness_pairs):
            for i in range(len(self.diameter_thickness_pairs), index + 1):
                effectiveThickness = effectiveThickness + self.diameter_thickness_pairs[i][1]

        return effectiveThickness

    def I(self):
        """

        :return: Appendix 13-6, effective moment of inertia of metal between multidiameter holes
        """
        X_bar = self.X_bar()
        firstTerm = 0
        otherSum = 0
        if len(self.diameter_thickness_pairs) > 1:
            for step in range(len(self.diameter_thickness_pairs) - 1, 1):
                b_n = self.__b_n(self.diameter_thickness_pairs[step][0])
                T_n = self.diameter_thickness_pairs[step][1]
                firstTerm += b_n * T_n ** 3
                otherSum += b_n * T_n * (self.__effectiveStepThickness(step) - X_bar) ** 2
        firstTerm += self.__b_n(self.diameter_thickness_pairs[-1][0]) * self.diameter_thickness_pairs[-1][1] ** 3
        otherSum += self.__b_n(self.diameter_thickness_pairs[-1][0]) * self.diameter_thickness_pairs[-1][1] *(X_bar - 0.5 * self.diameter_thickness_pairs[-1][1]) ** 2

        return (1/12) * firstTerm + otherSum

    def c(self):
        X_bar = self.X_bar()
        return max(X_bar, self.t - X_bar)

    def D_E_bending(self):
        return self.p - (6.0 * self.I()) / ((self.t ** 2) *(self.c()))

    def D_E_membrane(self):
        D_E = 0
        for pair in self.diameter_thickness_pairs:
            D_E += pair[0] * pair[1]
        return D_E / self.t

    def e_b(self):
        return (self.p - self.D_E_bending()) / self.p

    def e_m(self):
        return (self.p - self.D_E_membrane()) /self.p


