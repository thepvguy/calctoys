if __name__ == "__main__":
    from Flange.Clamped.Appendix24 import *
    import os

    print("cwd:")
    print(os.getcwd())


    def pretty(d, indent=0):
        for key, value in d.items():
            print('\t' * indent + str(key))
            if isinstance(value, dict):
                pretty(value, indent + 1)
            else:
                print('\t' * (indent + 1) + str(value))


    condition = DesignCondition(
        is_operating=True,
        temperature=200,
        ambient_temperature=75,
        pressure=3000,
        bolting_is_controlled=False,
        has_retainer=False
    )

    bolt_material = Material(1)
    bolt = Bolting(
        material=bolt_material,
        diameter=1.75,
        root_area=1.98,
    )

    hub_material = Material(2)
    left_hub = Hub(
        sketch=Hub.Sketch.c,
        material=hub_material,
        outside_diameter=(18 + (12.75 + 2.75) * 2),
        inner_diameter=18,
        hub_cross_section_corner_radius=0.25,
        small_end_thickness=12.75,
        length_of_small_end=15,
        taper_length=2.75,
        neck_length=15,
        neck_thickness_at_shoulder=12.75,
        shoulder_thickness=7.321,
        shoulder_height=2.75,
        transition_angle=10,
        shoulder_angle=10,
        friction_angle=5
    )

    clamp_material = Material(3)
    clamp = Clamp(
        sketch=Clamp.Sketch.a,
        material=clamp_material,
        bolt_circle_radius=32.25,
        clamp_inside_diameter=43.75,
        clamp_width=28,
        effective_clamp_thickness=7.625,
        effective_clamp_gap=14,
        corner_radius=0.25,
        distance_from_bolt_circle_to_clamp_od=3.7,
        effective_lip_length=2.75,
        lug_height=15,
        lug_width=28
    )

    gasket = Gasket(
        outer_diameter=20,
        inner_diameter=18,
        gasket_factor=0,
        seating_stress=0
    )

    half_assy = HubAndClamp(
        design_data=condition,
        bolting=bolt,
        gasket=gasket,
        hub=left_hub,
        clamp=clamp
    )

    pretty(half_assy.get_result(), 0)
    print("x")
