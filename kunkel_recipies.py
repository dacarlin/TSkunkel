from uuid import uuid4

kinase_mix = [
    ( 'pnk_buffer', 3, 'rs16pc9rd5sg5d' ),
    ( 'atp', 1, 'rs16pccshb6cb4' ),
    ( 'pnk', 1, 'rs16pc9rd5hsf6' ),
    ( 'water', 18, 'rs17gmh5wafm5p' ),
]

pol_mix = [
    ( 'lig_buffer', 0.6, 'rs17sh5rzz79ct' ),
    ( 'dntps', 0.4, 'rs16pcb542c5rd' ),
    ( 'atp', 0.4, 'rs16pccshb6cb4' ),
    ( 'pol', 0.4, 'rs16pca2urcz74' ),
    ( 'ligase', 0.4, 'rs16pc8krr6ag7' ),
]

def provision_helper( p, recipe, n ):
    n += 1 # make N+1
    my_cont = p.ref( str(uuid4()), None, 'micro-2.0', discard=True )
    for reagent, vol, id in recipe:
        p.provision( id, my_cont.well(0), '{}:microliter'.format( vol * n ) )
    return my_cont

ramp = [ dict( temperature='{}:celsius'.format( 95-i ), duration='1:minute' ) for i in range( 70 ) ]
