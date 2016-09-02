from kunkel_recipies import kinase_mix_resus, nick_mix, pol_mix, provision_helper, ramp
from autoprotocol.container import Container
from autoprotocol.protocol import Ref
import pandas

def single_kunkel( p, params ):

    # before running this protocol
    # transform your plasmid into CJ236
    # grow a large culture (1 L)
    # purify plasmid DNA from the culture

    dsdna = True # set to False if you have already prepped ssDNA

    # format of input CSV of order:
    #
    #   mutant_label        arbitrary string
    #   sequence            sequence of mutageneic oligo

    # df = pandas.read_csv( params[ 'order_csv' ] )
    df = pandas.read_csv( 'example2.csv' )
    df = pandas.concat( [ df, df, df ] ) # makes 96 mutants for debugging
    df.index = range( len( df ) )

    ML = list( df.index ) # we will need this list of int a lot :)
    n = len( ML )

    # order the oligos into a plate, they come dry, resuspend in kinase mix
    oligo_pt = p.ref( 'oligo_pt', None, '96-pcr', storage='cold_20' )
    p.oligosynthesize( [ dict( destination='oligo_pt/{}'.format(idx), sequence=x.sequence, scale='10nm' ) for idx, x in df.iterrows() ] )

    # kinase
    kin_mix = p.ref( 'kin_mix', None, 'micro-2.0', discard=True )
    for reagent, vol, my_id in kinase_mix_resus:
        p.provision( my_id, kin_mix.well(0), '{}:microliter'.format( vol * ( n + 2 ) ) )
    p.transfer( kin_mix.well( 0 ), oligo_pt.wells( ML ), '10:microliter', mix_after=True )
    p.seal( oligo_pt )
    p.spin( oligo_pt, '2000:g', '30:second' )
    p.incubate( oligo_pt, 'warm_37', '60:minute' )
    # may have to remove some liquid here to bring down to  1 uL
    p.unseal( oligo_pt )

    # dilute
    p.dispense_full_plate( oligo_pt, 'water', '95:microliter' )

    # reagents and ssDNA master mix for annealing
    nk_mix = p.ref( 'nick_mix', None, 'micro-2.0', discard=True )
    if dsdna:
        for reagent, vol, my_id in nick_mix:
            p.provision( my_id, nk_mix.well( 0 ), '{}:microliter'.format( vol * ( n + 1 ) ) )
    else:
        p.provision( 'rs17sh5rzz79ct', nk_mix.well( 0 ), '{}:microliter'.format( .2 * ( n + 1 ) ) )
    my_ssdna = p.ref( 'my_ssdna', 'ct18bf7kgp82fs', 'micro-1.5', storage='cold_20' )
    p.transfer( my_ssdna.well( 0 ), nk_mix.well( 0 ), '{}:microliter'.format( 2 * ( n + 1 ) ) )
    if dsdna:
        p.incubate( nk_mix, 'warm_37', '1:hour' )

    # anneal
    anneal_pt = p.ref( 'anneal_pt', None, '384-pcr', storage='cold_20' )
    p.distribute( nk_mix.well( 0 ), anneal_pt.wells( ML ), '2:microliter' )
    p.stamp( oligo_pt, anneal_pt.well(0), '2:microliter' )
    p.seal( anneal_pt )
    p.thermocycle( anneal_pt, [ { "cycles": 1, "steps": ramp } ] )

    # polymerize
    pmix = provision_helper( p, pol_mix, n )
    p.unseal( anneal_pt )
    p.transfer( pmix.well(0), anneal_pt.wells( ML ), '2:microliter' )
    p.seal( anneal_pt )
    p.incubate( anneal_pt, 'ambient', '90:minute' )

    # transform, recover, and plate
    p.incubate( anneal_pt, 'cold_4', '10:minute' )
    p.unseal( anneal_pt )
    p.provision( 'rs16pbj944fnny', anneal_pt.wells( ML ), '10:microliter' )
    p.seal( anneal_pt )
    p.incubate( anneal_pt, 'cold_4', '10:minute' )
    grow_pt = p.ref( 'grow_pt', None, '96-deep', discard=True )
    p.dispense_full_plate( grow_pt, 'soc', '100:microliter' )
    p.unseal( anneal_pt )
    p.stamp( anneal_pt.well(0), grow_pt, '15:microliter' )
    p.cover( grow_pt )
    p.incubate( grow_pt, 'warm_37', '90:minute', shaking=True )

    # some trickery still around agar plates at Transcriptic
    my_plates = []
    for i in list( set( map( lambda x: x/6, ML ) ) ): # used to get indexes for the 6-well plates
        name = 'agar_{}'.format( i )
        my_ref = Ref( name, { 'reserve': 'ki17rs7j799zc2', 'store': { 'where': 'cold_4' } }, Container( None, p.container_type( '6-flat' ) ) )
        p.refs[ name ] = my_ref
        my_plates.append( my_ref.container )

    # sample-centric spreading (can't think of a simple way to do it another way)
    for m in ML:
        p.spread( grow_pt.well( m ), my_plates[ m / 6 ].well( m % 6 ), '55:microliter' )

    # incubate and image
    for agar_pt in my_plates:
        p.incubate( agar_pt, 'warm_37', '10:hour' )

    for idx, agar_pt in enumerate( my_plates ):
        #p.uncover( agar_pt )
        p.image_plate( agar_pt, mode='top', dataref='agar_{}'.format( idx ) )
        #p.cover( agar_pt )

    #grow_pt2 = p.ref( 'grow_pt2', None, '96-deep', discard=True )


if __name__ == '__main__':
    from autoprotocol.harness import run
    run( single_kunkel, 'single_kunkel' )
