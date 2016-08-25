from kunkel_recipies import kinase_mix, pol_mix, provision_helper, ramp
from autoprotocol.container import Container
from autoprotocol.protocol import Ref
import pandas

def single_kunkel( p, params ):

    # format of input CSV of order:
    #
    #   mutant_label        arbitrary string
    #   sequence            sequence of mutageneic oligo
    #   ssDNA               Transcriptic ID of ssDNA to use

    # df = pandas.read_csv( params[ 'order_csv' ] )
    df = pandas.read_csv( 'example.csv', index_col=0 )

    df = pandas.concat( [ df, df, df ] ) # makes 96 mutants for debugging
    df.index = range( len( df ) )        # for debugging only!

    # order the oligos
    oligo_pt = p.ref( 'oligo_pt', None, '96-pcr', storage='cold_20' )
    p.oligosynthesize( [ dict( destination='oligo_pt/{}'.format(idx), sequence=x.sequence, scale='10nm' ) for idx, x in df.iterrows() ] )
    # oligos come normalized, wet in 96-well plate
    # see http://www.idtdna.com/pages/products/dna-rna/96-and-384-well-plates
    # may need p.unseal( oligo_pt )

    # kinase
    ML = list( df.index ) # we will need this list of int a lot :)
    kin_mix = provision_helper( p, kinase_mix, len( df ) )
    kin_pt = p.ref( 'kin_pt', None, '96-pcr', discard=True )
    p.transfer( kin_mix.well( 0 ), kin_pt.wells( ML ), '23:microliter' )
    #p.transfer( oligo_pt.wells( ML ), kin_pt.wells( ML ), '7:microliter' )
    p.stamp( oligo_pt, kin_pt, '7:microliter' )
    p.seal( kin_pt )
    p.incubate( kin_pt, 'warm_37', '60:minute' )
    p.unseal( kin_pt )

    # dilute
    dil_pt = p.ref( 'dil_pt', None, '96-flat', discard=True )
    p.dispense_full_plate( dil_pt, 'water', '200:microliter' )
    p.stamp( kin_pt, dil_pt, '2:microliter' )
    p.cover( dil_pt )
    p.spin( dil_pt, '2000:g', '30:second' )
    p.uncover( dil_pt )

    # anneal
    my_ssdna = p.ref( 'my_ssdna', 'ct18bf7kgp82fs', 'micro-1.5', discard=True )
    mix = p.ref( 'mix', None, 'micro-1.5', storage='cold_20' )
    anneal_pt = p.ref( 'anneal_pt', None, '384-pcr', storage='cold_20' )
    p.transfer( mix.well(0), anneal_pt.wells( ML ), '2:microliter' )
    p.transfer( my_ssdna.well(0), anneal_pt.wells( ML ), '2.2:microliter' )
    p.seal( anneal_pt )
    p.thermocycle( anneal_pt, [ { "cycles": 1, "steps": ramp } ] )

    # polymerize
    pmix = provision_helper( p, pol_mix, len( df ) )
    p.unseal( anneal_pt )
    p.transfer( pmix.well(0), anneal_pt.wells_from( 0, len(df) ), '2:microliter', one_tip=False )
    p.seal( anneal_pt )
    p.incubate( anneal_pt, 'ambient', '90:minute' )

    # transform
    p.incubate( anneal_pt, 'cold_4', '20:minute' )
    p.provision( 'rs16pbj944fnny', anneal_pt.wells( ML ), '10:microliter' )
    p.cover( anneal_pt )
    p.incubate( anneal_pt, 'cold_4', '20:minute' )
    grow_pt = p.ref( 'grow_pt', None, '96-deep', discard=True )
    p.dispense_full_plate( grow_pt, 'tb_50ug_ml_kan', '100:microliter' )
    p.stamp( anneal_pt, grow_pt, '15:microliter' )
    p.cover( grow_pt )
    p.incubate( grow_pt, 'warm_37', '10:hour', shaking=True )

    agar_pt_names = []
    kit_item = Container( None, p.container_type( '6-flat' ) )
    for i in range( len(df)/6 ):
        name = 'agar_{}'.format( i )
        my_ref = Ref( name, { 'reserve': 'ki17rs7j799zc2', 'store': { 'where': 'cold_4' } }, kit_item )
        p.refs[ name ] = my_ref
        agar_pt_names.append( name )
    for i in range( len( df ) ):
        p.spread( grow_pt.well( i ), p.refs[ 'agar_{}'.format( i/6 ) ].container.well( i%6 ), '55:microliter' )

    agar_pts = [ p.refs[ n ] for n in agar_pt_names ]
    for agar_pt in agar_pts:
        #print agar_pt.name
        #print dir( agar_pt )
        p.cover( agar_pt )
        p.incubate( agar_pt, 'warm_37', '10:hour' )
        p.uncover( agar_pt )
        p.image_plate( agar_pt, mode='top', dataref='jeff' )
        p.cover( agar_pt )

if __name__ == '__main__':
    from autoprotocol.harness import run
    run( single_kunkel, 'single_kunkel' )
