import json 
import screed 
import pandas 

ssDNA_id = 'ct19384uh8abx3' # actually water 
growth_media = 'lb-broth-50ug-ml-kan'
mutants = [] 

df = pandas.read_csv( 'raw_collect.txt', sep='\s+' )
err = []

for index, series in df.iterrows():
  try:
    seq = [ record.sequence for record in screed.open( 'my_oligos/{}.fasta'.format( series.name ) ) ][0]
    my_oligo = dict( oligo_label=series.name, mutant_label=series.name, sequence=seq, scale='25nm',purification='standard' ) 
    mutants.append( my_oligo ) 
  except Exception as e:
    err += [ series ] 

refs=dict(ssDNA_source={'id':ssDNA_id, 'type': 'micro-2.0', 'store': 'cold_20' })
parameters=dict(construct_setup=dict(growth_media=growth_media,ssDNA='ssDNA_source/0',mutant_upload=mutants))
result = dict(refs=refs, parameters=parameters)
print json.dumps( result, indent=2 ) 

with open( 'error_log', 'w' ) as fn:
  fn.write( str( err ) ) 
