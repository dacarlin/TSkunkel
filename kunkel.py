from autoprotocol.protocol import Protocol
from autoprotocol.instruction import Instruction
from assemble import assemble
from transform import transform
from autoprotocol.container import WellGroup
from mutant import Mutant
from collections import OrderedDict
from autoprotocol import UserError
from modules.utils import printdatetime

def scale_default(length, scale, label):
    ok = True
    if scale == '25nm':
        ok = True if (length >= 15 and length <= 60) else False
    elif scale == '100nm':
        ok = True if (length >= 10 and length <= 90) else False
    elif scale == '250nm':
        ok = True if (length >= 5 and length <= 100) else False
    elif scale == '1um':
        ok = True if (length >= 5 and length <= 100) else False
    else:
        ok = False
    if not ok:
        raise UserError("The specified oligo '%s' is %s base pairs long."
                           " This sequence length is invalid for the scale "
                           "of synthesis chosen (%s)." % (label, length, scale))


def kunkel_full(protocol, params):
    growth_media = params["construct_setup"]['growth_media']
    #num_colonies = params["construct_setup"]['num_colonies']
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
        'constructs': [{
            'mutant_name': mu.name,
            'oligos': mu.oligos} for mu in mutant_constructs],
        'mutant_objs': mutant_constructs
    }

    annealing_plate = assemble(protocol, assemble_params)
    protocol.unseal(annealing_plate)

    transform_params = {
        #'num_colonies': num_colonies,
        'growth_media': growth_media,
        'constructs': [mu.anneal_well for mu in mutant_constructs],
        'mutant_objs': mutant_constructs

    }

    # get agar plates back from transform protocol
    agar_plates = transform( protocol, transform_params )

    for agar_plate in agar_plates:
        protocol.cover( agar_plate )

if __name__ == '__main__':
    from autoprotocol.harness import run
    run(kunkel_full, 'Kunkel')
