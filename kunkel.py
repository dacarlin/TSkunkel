import math
from datetime import datetime

from transcriptic_utils import *
from autoprotocol import UserError
from autoprotocol.container import WellGroup
from autoprotocol.pipette_tools import aspirate_source, depth

def the_date():
    return datetime.now().strftime('%Y-%m-%d')

def kunkel_full(protocol, params):
    growth_media = params["construct_setup"]['growth_media']
    ssDNA = params["construct_setup"]['ssDNA']
    mutant_constructs = []

    # make mutant objects for accessibility
    construct_collect = {}
    for csv_row in params["construct_setup"]['mutant_upload']:
        if csv_row["mutant_label"] not in construct_collect.keys():
            construct_collect[csv_row["mutant_label"]] = []
            construct_collect[csv_row["mutant_label"]].append(
                {
                "sequence": csv_row["sequence"],
                "purification": csv_row["purification"],
                "scale": csv_row["scale"],
                "oligo_label": csv_row["oligo_label"]

            })
        else:
            construct_collect[csv_row["mutant_label"]].append(
                {
                "sequence": csv_row["sequence"],
                "purification": csv_row["purification"],
                "scale": csv_row["scale"],
                "oligo_label": csv_row["oligo_label"]

            }
                )

    oligo_collect = {}
    for row in params["construct_setup"]["mutant_upload"]:
        if (row["sequence"] not in oligo_collect.keys() and row["oligo_label"] in protocol.refs.keys()):
            raise RuntimeError("You cannot specify two different "
                   "oligos to be synthesized with the "
                   "same name %s" % row['oligo_label'])
        elif row["sequence"] not in oligo_collect.keys():
            oligo_collect[row["sequence"]] = {
                "sequence": row["sequence"],
                "purification": row["purification"],
                "scale": row["scale"],
                "destination": protocol.ref(row["oligo_label"], None, "micro-2.0", storage="cold_4").well(0)
            }

    for mut in construct_collect.keys():
        mut_oligos = [o for o in construct_collect[mut]]
        mutant = Mutant(mut)
        for oligo in mut_oligos:
            mutant.add_oligos(oligo_collect[oligo["sequence"]]["destination"])
        mutant_constructs.append(mutant)


    oligos_to_synthesize = []
    for o in oligo_collect.keys():
        scale_default(len(oligo_collect[o]["sequence"]), oligo_collect[o]["scale"], oligo_collect[o]["destination"].container.name)
        oligos_to_synthesize.append(oligo_collect[o])
    protocol.oligosynthesize(oligos_to_synthesize)

    assemble_params = {
        'ssDNA': ssDNA,
        'constructs': [{ 'mutant_name': mu.name, 'oligos': mu.oligos } for mu in mutant_constructs],
        'mutant_objs': mutant_constructs
    }

    # print assemble_params

    # this looks like
    # {
    #   'constructs':
    #     {
    #       'name': 'mutant_1',
    #       'oligos': [ Well objects ] # used when diluting
    #     },
    #   'mutant_objs': [ Mutant objects ] # used when?
    # }

    def assemble(protocol, params):

        # provision atp for entire protocol
        atp_reagents = {'atp': {"resource_id": 'rs16pccshb6cb4',
                        'reagent_ratio': 0.1},
                        'water': {"resource_id": 'rs17gmh5wafm5p', 'reagent_ratio': 0.9}}

        reagents = {'pnkbuffer': {"resource_id": 'rs16pc9rd5sg5d', "reagent_ratio": 3},
                    'water': {"resource_id": 'rs17gmh5wafm5p', "reagent_ratio": 18},
                    'pnk': {"resource_id": 'rs16pc9rd5hsf6', "reagent_ratio": 1}}

        pol_reagents = {"buffer": {"resource_id": 'rs17sh5rzz79ct', "reagent_ratio": 0.6},
                    "t4ligase": {"resource_id": 'rs16pc8krr6ag7', "reagent_ratio": 0.4},
                    "t7polymerase": {"resource_id": 'rs16pca2urcz74', "reagent_ratio": 0.4},
                    "dntp": {"resource_id": 'rs16pcb542c5rd', "reagent_ratio": 0.4}
                    }

        # general parameters
        ssDNA = params['ssDNA']
        constructs = [construct['oligos'] for construct in params['constructs']]
        num_constructs = len(constructs)
        flattened = [val for oligo in constructs for val in oligo]
        oligos = list({v: v for v in flattened}.values())
        num_oligos = len(oligos)
        mm_mult = 1.3
        num_rxts_plus = 3

        # refs
        water = protocol.ref( 'water', cont_type='96-deep', discard=True )
        atp = protocol.ref("atp_10mM", cont_type='micro-1.5', discard=True).well(0)
        atp_rxts = (num_oligos + num_constructs)

        # provisioning
        protocol.dispense_full_plate( water, 'water', '1000:microliter' )
        provision_reagents(atp_reagents, atp, atp_rxts, mm_mult, num_rxts_plus=6)
        kinase_mix = []
        for i in range(int(math.ceil(num_oligos/60.0))):
            kinase_mix.append(protocol.ref("kinase_mix-{}".format(i + 1), None, "micro-1.5", discard=True).well(0))
        provision_reagents(reagents, kinase_mix, num_oligos, mm_mult, num_rxts_plus)
        protocol.transfer(atp, kinase_mix, "%s:microliter" % ((num_oligos + num_rxts_plus) * 1 * mm_mult), new_group=True)

        #   kinase oligos
        kinase_oligo_plate = protocol.ref( "kinase_oligo_plate_{}".format( the_date() ), None, "96-pcr", discard=True)
        wells_to_kinase = kinase_oligo_plate.wells_from( 0, num_oligos )
        protocol.transfer( kinase_mix, wells_to_kinase, "23:microliter", blowout_buffer=True, pre_buffer='15:microliter', one_tip=True, one_source=True )

        for i, oligo in enumerate( oligos ):
            protocol.transfer( oligo, wells_to_kinase[i], "7:microliter",
                                  mix_after=False, new_group=det_new_group(i),
                                  aspirate_source=aspirate_source(depth=depth("ll_following",
                                                                  lld="pressure",
                                                                  distance="0.0:meter")),
                                  blowout_buffer=True,
                                  pre_buffer=11 )

        protocol.seal(kinase_oligo_plate)
        protocol.incubate( kinase_oligo_plate, 'warm_37', '65:minute' )

        # dilute oligos
        protocol.unseal(kinase_oligo_plate)

        diluted_oligo_plate = protocol.ref("dilute_oligo_plate", None, "96-flat", discard=True)
        p.dispense_full_plate( diluted_oligo_plate, 'water', '200:microliter' )
        diluted_oligo_wells = diluted_oligo_plate.wells_from(0, num_constructs)

        for i, m in enumerate(constructs):
            for kin_oligo in m:
                index = list(oligos).index(kin_oligo)
                protocol.transfer(kinase_oligo_plate.well(index), diluted_oligo_plate.well(i),"2:microliter",mix_after=False,mix_vol="5:microliter")
            diluted_oligo_plate.well(i).set_name(params['constructs'][i]['mutant_name'])

        protocol.cover(diluted_oligo_plate)
        protocol.spin(diluted_oligo_plate, "250:g", "2:minute")

        mm_mult_ssDNA = 1.5

        mix_plate = protocol.ref("mix_plate", None, "96-pcr", discard=True)
        ssDNA_mix = mix_plate.well(0)
        protocol.transfer(ssDNA,
                          ssDNA_mix,
                          "%s:microliter" % ((num_constructs + num_rxts_plus) * 2.0 * mm_mult_ssDNA),
                          blowout_buffer=True,
                          **transfer_kwargs( num_constructs + 1 ) )
        protocol.provision('rs17sh5rzz79ct', ssDNA_mix, "%s:microliter" % ((num_constructs + num_rxts_plus) * 0.2 * mm_mult_ssDNA))

        # anneal oligos
        protocol.uncover(diluted_oligo_plate)
        annealing_plate = protocol.ref("annealing_oligo_plate", None, "384-pcr", storage="cold_20")
        anneal_wells = annealing_plate.wells_from(0, num_constructs)
        protocol.transfer(ssDNA_mix,
                          anneal_wells.wells,
                          "2.2:microliter",
                          dispense_speed="50:microliter/second",
                          blowout_buffer=True,
                          **transfer_kwargs(7, True, True))

        for dil_oligo, reaction in zip(diluted_oligo_wells.wells, anneal_wells.wells):
            protocol.transfer(dil_oligo,
                              reaction,
                              "2:microliter",
                              aspirate_source=aspirate_source(depth("ll_bottom", distance=".001:meter")),
                              mix_after=True,
                              mix_vol="2:microliter",
                              flowrate="50:microliter/second",
                              repetitions=2,
                              blowout_buffer=True,
                              new_group=det_new_group(i),
                              **transfer_kwargs(5))
            reaction.set_name(dil_oligo.name)

        protocol.seal(annealing_plate)
        protocol.spin(annealing_plate, "250:g", "2:minute")
        protocol.thermocycle(annealing_plate, [{
            "cycles": 1,
            "steps": thermocycle_ramp("95:celsius", "25:celsius", "60:minute", "4:minute")
            }],
            volume="5:microliter",
            dataref=None,
            dyes=None)

        # polymerize
        protocol.unseal(annealing_plate)
        polymerize_MM = mix_plate.well(12)

        provision_reagents(pol_reagents, polymerize_MM, num_constructs, mm_mult, num_rxts_plus)
        protocol.transfer(atp, polymerize_MM, "%s:microliter" % ((num_constructs + num_rxts_plus) * 0.4 * mm_mult), new_group=True)

        for reaction in anneal_wells.wells:
            protocol.transfer(polymerize_MM, reaction, "2.2:microliter",
                              mix_after=False,
                              blowout_buffer=True,
                              pre_buffer=0)

            if 'mutant_objs' in params.keys():
                mut = next(m for m in params['mutant_objs'] if m.name == reaction.name)
                mut.anneal_well = reaction

        protocol.seal(annealing_plate)
        protocol.incubate(annealing_plate, "ambient", "1.5:hour")

        # pass plate back for unsealing
        return annealing_plate

    annealing_plate = assemble(protocol, assemble_params)
    protocol.unseal(annealing_plate)

    transform_params = {
        'growth_media': growth_media,
        'constructs': [ mu.anneal_well for mu in mutant_constructs ],
        'mutant_objs': mutant_constructs
    }

    ### begin transform

    def transform(protocol, params):
        # general parameters
        constructs = params['constructs']
        num_constructs = len(constructs)
        plates = list(set([construct.container for construct in constructs]))
        if len(plates) != 1:
            raise UserError('You can only transform aliquots from one common container.')

        mm_mult = 1.3

        transformation_plate = protocol.ref("transformation_plate", None, "96-pcr", discard=True)
        protocol.incubate( transformation_plate, "cold_20", "10:minute" )

        transformation_wells = transformation_plate.wells_from(0, num_constructs)
        for i in range(num_constructs):
            protocol.provision("rs16pbjc4r7vvz", transformation_wells[i], "10:microliter")

        for i, well in enumerate(constructs):
            protocol.transfer(well, transformation_wells[i], "5.0:microliter",
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
        protocol.dispense_full_plate( transformation_plate, 'soc', '90:microliter' )
        protocol.seal(transformation_plate)
        protocol.incubate(transformation_plate, "warm_37", "10:minute", shaking=True )
        protocol.unseal(transformation_plate)

        agar_plates = []
        agar_wells = WellGroup([])

        for well in range(0, len(transformation_wells), 6):
            agar_plate = ref_kit_container(protocol,
                                           "agar-{}_{}".format( len( agar_plates ), the_date() ),
                                           "6-flat",
                                           return_agar_plates(6)[params['growth_media']],
                                           discard=False, store='cold_4')
            agar_plates.append(agar_plate)
            for i, w in enumerate( transformation_wells[well:well + 6] ) :
                protocol.spread( w, agar_plate.well(i), "100:microliter" )
                agar_wells.append( agar_plate.well(i).set_name(w.name) )

        for agar_p in agar_plates:
            protocol.incubate( agar_p, 'warm_37', '12:hour' )
            protocol.image_plate( agar_p, mode='top', dataref=agar_p.name )
            protocol.cover( agar_plate )

if __name__ == '__main__':
    from autoprotocol.harness import run
    run(kunkel_full, 'Kunkel')
