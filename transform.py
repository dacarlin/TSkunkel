import math
from autoprotocol import UserError
from modules.utils import *


def transform(protocol, params):

    # general parameters
    constructs = params['constructs']
    num_constructs = len(constructs)
    plates = list(set([construct.container for construct in constructs]))
    if len(plates) != 1:
        raise UserError('You can only transform aliquots from one common container.')

    # **** need to be able to check if plate is sealed to add run-chaining ****
    mm_mult = 1.3

    transformation_plate = protocol.ref("transformation_plate", None, "96-pcr", discard=True)
    protocol.incubate(transformation_plate, "cold_20", "10:minute")

    transformation_wells = transformation_plate.wells_from(0, num_constructs)
    for i in range(num_constructs):
        protocol.provision("rs16pbjc4r7vvz", transformation_wells[i], "50:microliter")

    for i, well in enumerate(constructs):
        protocol.transfer(well, transformation_wells[i], "2.0:microliter",
                          dispense_speed="10:microliter/second",
                          mix_after=False,
                          new_group=det_new_group(i))
        if well.name:
            transformation_wells[i].set_name(well.name)
        else:
            transformation_wells[i].set_name('construct_%s' % (i+1))

    # NEED to confirm second de-seal is working OR move to cover/uncover 96-flat
    protocol.seal(transformation_plate)
    protocol.incubate(transformation_plate, "cold_4", "20:minute", shaking=False, co2=0)
    protocol.unseal(transformation_plate)
    protocol.dispense_full_plate( transformation_plate, 'soc', '50:microliter' )
    protocol.seal(transformation_plate)
    protocol.incubate(transformation_plate, "warm_37", "10:minute", shaking=True)
    protocol.unseal(transformation_plate)

    # spread on agar plates

    # kan "ki17rs7j799zc2"
    # amp "ki17sbb845ssx9"
    # specto "ki17sbb9r7jf98"
    # cm "ki17urn3gg8tmj"
    # "noAB" "ki17reefwqq3sq"

    agar_plates = []
    agar_wells = WellGroup([])
    for well in range(0, len(transformation_wells), 6):
        agar_name = "agar-%s_%s" % (len(agar_plates), printdatetime(time=False))
        agar_plate = ref_kit_container(protocol, agar_name, "6-flat", "ki17rs7j799zc2", discard=False, store='cold_4')
        agar_plates.append(agar_plate)
        for i, w in enumerate(transformation_wells[well:well + 6]):
            protocol.spread(w, agar_plate.well(i), "100:microliter")
            agar_wells.append(agar_plate.well(i).set_name(w.name))

    for agar_p in agar_plates:
        protocol.incubate( agar_p, 'warm_37', '12:hour' )
        protocol.image_plate( agar_p, mode='top', dataref=agar_p.name )

    # return agar plates to end protocol
    return agar_plates

if __name__ == '__main__':
    from autoprotocol.harness import run
    run(transform, 'Transform')
