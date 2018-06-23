import math
import abc

class Load(abc.ABC):
    def __init__(self):
        pass

    @abc.abstractmethod
    def get_loads(self):
        pass


class LateralLoad(Load):
    def __init__(self, magnitude, distance_from_datum, angle,  name="Load"):
        super().__init__()
        self.__magnitude = magnitude
        self.__elevation = distance_from_datum
        self.__angle = angle
        self.__name = name

    def get_loads(self):
        return [self]

    @property
    def magnitude(self):
        return self.__magnitude

    @property
    def elevation(self):
        return self.__elevation

    @property
    def angle(self):
        return self.__angle

    @property
    def moment_x(self, at_elevation=None):
        if not at_elevation:
            at_elevation = self.elevation
        return self.force_x * at_elevation

    @property
    def moment_y(self, at_elevation=None):
        if not at_elevation:
            at_elevation = self.elevation

        return self.force_y * at_elevation

    @property
    def force_x(self):
        return self.magnitude * math.cos(math.radians(self.angle))

    @property
    def force_y(self):
        return self.magnitude * math.sin(math.radians(self.angle))


class VerticalLoad(Load):
    def __init__(self, magnitude, distance_from_datum, distance_from_cl, angle, name="Load"):
        super().__init__()
        self.__magnitude = magnitude
        self.__elevation = distance_from_datum
        self.__moment_arm = distance_from_cl
        self.__angle = angle
        self.__name = name

    def get_loads(self):
        return [self]

    @property
    def magnitude(self):
        return self.__magnitude

    @property
    def elevation(self):
        return self.__elevation

    @property
    def angle(self):
        return self.__angle

    @property
    def moment_x(self):
        return self.force_x * self.__moment_arm

    @property
    def moment_y(self):
        return self.force_y * self.__moment_arm

    @property
    def force_x(self):
        return self.magnitude * math.cos(math.radians(self.angle))

    @property
    def force_y(self):
        return self.magnitude * math.sin(math.radians(self.angle))


class PlatformLoad(Load):
    def __init__(self,
                 start_angle,
                 end_angle,
                 vessel_od,
                 platform_clearance,
                 platform_width,
                 railing_unit_weight,
                 grating_unit_weight,
                 distance_from_datum,
                 name="Platform Load"
                 ):
        super().__init__()
        self.__start_angle = start_angle
        self.__end_angle = end_angle
        self.__vessel_od = vessel_od
        self.__platform_clearance = platform_clearance
        self.__platform_width = platform_width
        self.__railing_unit_weight = railing_unit_weight
        self.__grating_unit_weight = grating_unit_weight
        self.__distance_from_datum = distance_from_datum
        self.__name = name

    @property
    def __r_i(self):
        """

        :return: Inside Radius
        """
        return 0.5 * self.__vessel_od + self.__platform_clearance

    @property
    def __r_o(self):
        """

        :return: Outside Radius
        """
        return self.__r_i + self.__platform_width

    @property
    def __alpha_radians(self):
        """

        :return: half the included angle of the platform
        """
        return math.radians(abs(self.__end_angle - self.__start_angle)) * 0.5

    def grating_area(self):
        return self.__alpha_radians * self.__platform_width * (2 * self.__r_o - self.__platform_width)

    def outer_railing_length(self):
        return 2 * self.__alpha_radians * self.__r_o

    def inner_railing_length(self):
        # 2 * self.__alpha_radians * self.__r_i
        return 0

    def radial_railing_length(self):
        return self.__platform_width * 2

    def grating_centroid(self):
        return self.__r_o - \
               self.__r_o * \
               (
                       1 -
                       ((2 * math.sin(self.__alpha_radians))/(3 * self.__alpha_radians)) *
                        (1 - (self.__platform_width/self.__r_o) + (1/(2 - (self.__platform_width/self.__r_o))))
               )

    def outer_railing_centroid(self):
        return self.__r_o * math.sin(self.__alpha_radians) / self.__alpha_radians

    def inner_railing_centroid(self):
        # return self.__r_i * math.sin(self.__alpha_radians) / self.__alpha_radians
        return 0

    def radial_railing_centroid(self):
        return 0.5 * (self.__r_o + self.__r_i) * math.cos(self.__alpha_radians)

    def grating_weight(self):
        return self.__grating_unit_weight * self.grating_area()

    def outer_railing_weight(self):
        return self.outer_railing_length() * self.__railing_unit_weight

    def inner_railing_weight(self):
        return self.inner_railing_length() * self.__railing_unit_weight

    def radial_railing_weight(self):
        return self.radial_railing_length() * self.__railing_unit_weight

    def total_weight(self):
        return sum([self.grating_weight(),
                    self.outer_railing_weight(),
                    self.inner_railing_weight(),
                    self.radial_railing_weight()]
                   )

    def effective_centroid(self):
        return (self.grating_weight() * self.grating_centroid() +
                self.outer_railing_weight() * self.outer_railing_centroid() +
                self.inner_railing_weight() * self.inner_railing_centroid() +
                self.radial_railing_weight() * self.radial_railing_centroid()) / self.total_weight()

    def load_angle(self):
        return 0.5 * (self.__end_angle - self.__start_angle) + self.__start_angle

    def get_loads(self):
        return [VerticalLoad(self.total_weight(), self.__distance_from_datum, self.effective_centroid(), self.load_angle(), self.__name)]


class LoadingModel:
    def __init__(
            self,
            support_length,
            support_attachment_elevation,
            vessel_tangent_distance_below_attachment_point
    ):
        self.__support_length = support_length
        self.__support_attachment_elevation = support_attachment_elevation
        self.__vessel_tangent_distance_below_attachment_point = vessel_tangent_distance_below_attachment_point
        self.__loads = []

    def add_load_source(self, load_source):
        if not load_source:
            return

        for load in load_source.get_loads():
            self.__loads.append(load)
        sorted(self.__loads, key=lambda x: x.elevation)

    def total_force_x(self):
        return sum([s.force_x for s in self.__loads])

    def total_force_y(self):
        return sum([s.force_y for s in self.__loads])

    def total_moment_x(self):
        return sum([x.moment_x for x in self.__loads])

    def total_moment_y(self):
        return sum([y.moment_y for y in self.__loads])

    def base_shear(self):
        return max(abs(self.total_force_x()), abs(self.total_force_y()))

    def direction(self):
        return math.degrees(math.atan2(self.total_force_y(), self.total_force_x()))

    def total_moment(self):
        return ((self.total_moment_x() ** 2) + (self.total_moment_y() ** 2)) ** 0.5


if __name__ == "__main__":
    load_model = LoadingModel(support_length=72, support_attachment_elevation=5, vessel_tangent_distance_below_attachment_point=2)

    platform_load = PlatformLoad(
        start_angle=270,
        end_angle=360,
        vessel_od=73,
        platform_clearance=4,
        platform_width=48,
        railing_unit_weight=(18/12),
        grating_unit_weight=(13.9/144),
        distance_from_datum=110
    )
    load_model.add_load_source(platform_load)
    '''
    load_model.add_load_source(LateralLoad(magnitude=2500.0, distance_from_datum=725.0, angle=0.0, name="1"))
    load_model.add_load_source(LateralLoad(magnitude=2500.0, distance_from_datum=612.50, angle=0.0, name="2"))
    load_model.add_load_source(VerticalLoad(magnitude=250.0, distance_from_datum=500.0, angle=180.0, name="3", distance_from_cl=500))
    load_model.add_load_source(VerticalLoad(magnitude=250.0, distance_from_datum=500.0, angle=180.0, name="4", distance_from_cl=500))
    '''

    print("\nPlatform:")
    print("grating_area is : {}".format(platform_load.grating_area()))
    print("outer_railing_length is : {}".format(platform_load.outer_railing_length()))
    print("inner_railing_length is : {}".format(platform_load.inner_railing_length()))
    print("radial_railing_length is : {}".format(platform_load.radial_railing_length()))
    print("grating_centroid is : {}".format(platform_load.grating_centroid()))
    print("outer_railing_centroid is : {}".format(platform_load.outer_railing_centroid()))
    print("inner_railing_centroid is : {}".format(platform_load.inner_railing_centroid()))
    print("radial_railing_centroid is : {}".format(platform_load.radial_railing_centroid()))
    print("grating_weight is : {}".format(platform_load.grating_weight()))
    print("outer_railing_weight is : {}".format(platform_load.outer_railing_weight()))
    print("inner_railing_weight is : {}".format(platform_load.inner_railing_weight()))
    print("radial_railing_weight is : {}".format(platform_load.radial_railing_weight()))
    print("total_weight is : {}".format(platform_load.total_weight()))
    print("effective_centroid is : {}".format(platform_load.effective_centroid()))
    print("load_angle is : {}".format(platform_load.load_angle()))

    print("\nLoad model:")
    print("total_force_x is : {0}".format(load_model.total_force_x()))
    print("total_force_y is : {0}".format(load_model.total_force_y()))
    print("base_shear is : {0}".format(load_model.base_shear()))
    print("direction is : {0}".format(load_model.direction()))
    print("total_moment_x is : {0}".format(load_model.total_moment_x()))
    print("total_moment_y is : {0}".format(load_model.total_moment_y()))
    print("total_moment is : {0}".format(load_model.total_moment()))
