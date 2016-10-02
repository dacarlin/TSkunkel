# Kunkel mutagenesis
# Originally by Transcriptic, Inc., being cleaned up by Alex Carlin
# input is designed oligos as text file and prepped aliquot of ssDNA
# output is agar plates with Kunkel transformants
#
# limitations of this script
# - maximum number of mutants at once is 96
# -
#

# if we want to implement plate oligos, then
# - minimum number of oligos in single order is 24

import math
from autoprotocol.pipette_tools import aspirate_source, depth
from autoprotocol import UserError
from modules.utils import *
from autoprotocol.protocol import Protocol
from autoprotocol.instruction import Instruction
from collections import OrderedDict

# shitty functions we need to get rid of
def provision_reagents(reagents, dest, num_rxts, mm_mult=1.3, num_rxts_plus=3.0):
    for reagent in reagents.values():
        protocol.provision(reagent['resource_id'], dest, "%s:microliter" % ((num_rxts + num_rxts_plus) * reagent['reagent_ratio'] * mm_mult))

def transfer_kwargs(pre_buffer, one_tip=False, one_source=False):
    kwargs = {"one_tip": one_tip,
          "one_source": one_source,
          "pre_buffer": "%s:microliter" % pre_buffer,
          "blowout_buffer": True}
    return(kwargs)


import pandas

# Kunkel protocol
def my_kunkel( p, run_params ):

    # inputs
    ssDNA = run_params[ "construct_setup" ][ 'ssDNA' ]
    antibiotic = run_params[ 'construct_setup' ][ 'growth_media' ]
    df = pandas.DataFrame( run_params[ 'construct_setup' ][ 'mutant_upload' ] )

    # calcuate indexes for later
    df[ 'mutant_index' ] = pandas.factorize( df['mutant_label'] )[ 0 ]
    df[ 'agar_plate' ] = df.mutant_index.map( lambda x: x // 6 )
    df[ 'agar_well' ] = df.mutant_index.map( lambda x: x % 6 )
    #print df

    # oligo order
    my_order = []
    for i in df.sequence.unique():
        my_params = {
            'sequence': i, 'destination': p.ref( i, None, 'micro-2.0', storage='ambient' ),
            'scale': '25nm', 'purification': 'standard' }
        my_order.append( my_params )

    p.oligosynthesize( my_order )

    # resuspend oligos

    # assembly

    num_oligos = len( my_order )
    num_constructs = len( df )
    n1 = lambda x: '{}:microliter'.format( x * ( len( df ) + 1 ) )
    my_list = range( len( df ) )

    # provision kinase mix

    kinase_reagents = [
        ( 'pnkbuffer', 'rs16pc9rd5sg5d', 3 ),
        ( 'water', 'rs17gmh5wafm5p', 18 ),
        ( 'pnk', 'rs16pc9rd5hsf6', 1 ),
        ( 'atp', 'rs16pccshb6cb4', 0.1 ) # FIX ME: this is 100 mM ATP!!
    ]

    kinase_mx = p.ref( 'kinase_mx', None, 'micro-2.0', discard=True )
    for __, reagent, ratio in kinase_reagents:
        p.provision( reagent, kinase_mx.well(0), n1( ratio ) )

    # kinase

    kinase_pt = p.ref( 'kinase_pt', None, '96-pcr', discard=True )
    my_oligos = [ i['destination'] for i in my_order ]
    p.transfer( kinase_mx.well(0), kinase_pt.wells( my_list ), '23:microliter', pre_buffer='15:microliter' )
    for i, o in enumerate( my_oligos ):
        p.transfer( o.well(0), kinase_pt.well( i ), '7:microliter', pre_buffer='10:microliter' )
    p.seal( kinase_pt )
    p.incubate( kinase_pt, 'warm_37', '1:hour' )

    # dilute
    dilute_pt = p.ref( 'dilute_pt', None, '96-flat', discard=True )
    p.dispense_full_plate( dilute_pt, 'water', '200:microliter' )
    for idx, row in df.iterrows():
        p.transfer( p.refs[ row[ 'sequence' ] ].container.well( 0 ), dilute_pt.well( row[ 'mutant_index' ] ), '2:microliter' )
    p.cover( dilute_pt )
    p.spin( dilute_pt, '250:g', '2:minute' )

    # optional: start with plasmid DNA
    # Nb.BsrDI nicks bottom strand of pET29
    start_with_dsDNA = False
    if start_with_dsDNA:
        p.provision( 'nicking endonuclease Nb.BsrDI' )
        p.provision( 'exonuclease' )
        p.thermocycle( 'mix', '65:celsius', '1:hour' )

    #
    #
    # protocol.cover(diluted_oligo_plate)
    # protocol.spin(diluted_oligo_plate, "250:g", "2:minute")
    #
    # # make ssDNA_mastermix
    # mm_mult_ssDNA = 1.5
    #
    # mix_plate = protocol.ref("mix_plate", None, "96-pcr", discard=True)
    # ssDNA_mix = mix_plate.well(0)
    # protocol.transfer(ssDNA,
    #                   ssDNA_mix,
    #                   "%s:microliter" % ((num_constructs + num_rxts_plus) * 2.0 * mm_mult_ssDNA),
    #                   **transfer_kwargs((num_constructs + 1) * 1))
    # protocol.provision('rs17sh5rzz79ct', ssDNA_mix, "%s:microliter" % ((num_constructs + num_rxts_plus) * 0.2 * mm_mult_ssDNA))
    #
    # # anneal oligos
    # protocol.uncover(diluted_oligo_plate)
    # annealing_plate = protocol.ref("annealing_oligo_plate", None, "384-pcr", storage="cold_20")
    # anneal_wells = annealing_plate.wells_from(0, num_constructs)
    # protocol.transfer(ssDNA_mix,
    #                   anneal_wells.wells,
    #                   "2.2:microliter",
    #                   dispense_speed="50:microliter/second",
    #                   **transfer_kwargs(7, True, True))
    #
    # for dil_oligo, reaction in zip(diluted_oligo_wells.wells, anneal_wells.wells):
    #     protocol.transfer(dil_oligo,
    #                       reaction,
    #                       "2:microliter",
    #                       aspirate_source=aspirate_source(depth("ll_bottom", distance=".001:meter")),
    #                       mix_after=True,
    #                       mix_vol="2:microliter",
    #                       flowrate="50:microliter/second",
    #                       repetitions=2,
    #                       new_group=det_new_group(i),
    #                       **transfer_kwargs(5))
    #     reaction.set_name(dil_oligo.name)
    #
    # protocol.seal(annealing_plate)
    # protocol.spin(annealing_plate, "250:g", "2:minute")
    # protocol.thermocycle(annealing_plate, [{
    #     "cycles": 1,
    #     "steps": thermocycle_ramp("95:celsius", "25:celsius", "60:minute", "4:minute")
    #     }],
    #     volume="5:microliter",
    #     dataref=None,
    #     dyes=None)
    #
    # # polymerize
    # protocol.unseal(annealing_plate)
    # polymerize_MM = mix_plate.well(12)
    # reagents = {"buffer": {"resource_id": 'rs17sh5rzz79ct', "reagent_ratio": 0.6},
    #             "t4ligase": {"resource_id": 'rs16pc8krr6ag7', "reagent_ratio": 0.4},
    #             "t7polymerase": {"resource_id": 'rs16pca2urcz74', "reagent_ratio": 0.4},
    #             "dntp": {"resource_id": 'rs16pcb542c5rd', "reagent_ratio": 0.4}
    #             }
    #
    # provision_reagents(reagents, polymerize_MM, num_constructs, mm_mult, num_rxts_plus)
    # protocol.transfer(atp, polymerize_MM, "%s:microliter" % ((num_constructs + num_rxts_plus) * 0.4 * mm_mult), new_group=True)
    #
    # for reaction in anneal_wells.wells:
    #     protocol.transfer(polymerize_MM,
    #                       reaction,
    #                       "2.2:microliter",
    #                       mix_after=False,
    #                       **transfer_kwargs(10))
    #     if 'mutant_objs' in params.keys():
    #         mut = next(m for m in params['mutant_objs'] if m.name == reaction.name)
    #         mut.anneal_well = reaction
    #
    # protocol.seal(annealing_plate)
    # protocol.incubate(annealing_plate, "ambient", "1.5:hour")
    # protocol.unseal(annealing_plate)
    #
    # now the transformations
    # transformation_plate = protocol.ref("transformation_plate", None, "96-pcr", discard=True)
    # protocol.incubate(transformation_plate, "cold_20", "10:minute")
    # protocol.stamp( anneal_pt, transformation_plate, '3:microliter' )
    # protocol.incubate( transformation_plate, 'cold_20', '10:minute' )
    # protocol.dispense( transformation_plate, 'soc', '200:microliter' )
    # protocol.incubate( transformation_plate, '1:hour', 'warm_37', shaking=True )
    #
    # agar_plates = []
    # for well in range(0, len(transformation_wells), 6):
    #     agar_plate = ref_kit_container(protocol,
    #                                    "agar-%s_%s" % (len(agar_plates), printdatetime(time=False)),
    #                                    "6-flat",
    #                                    return_agar_plates(6)[params['growth_media']],
    #                                    discard=False, store='cold_4')
    #     agar_plates.append(agar_plate)
    #     for i, w in enumerate(transformation_wells[well:well + 6]):
    #         protocol.spread(w, agar_plate.well(i), "100:microliter")
    #         agar_wells.append(agar_plate.well(i).set_name(w.name))
    #
    # for agar_p in agar_plates:
    #     protocol.incubate( agar_p, 'warm_37', '12:hour' )
    #     protocol.image_plate( agar_p, mode='top', dataref=agar_p.name )
    #     protocol.cover( agar_plate )

if __name__ == '__main__':
    from autoprotocol.harness import run
    run( my_kunkel, 'unkel')
