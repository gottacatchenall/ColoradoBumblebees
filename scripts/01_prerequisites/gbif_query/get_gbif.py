from pygbif import occurrences
from pygbif import species

import pandas as pd

print("Warning: this script requires you enter your GBIF credentials. Never upload those credentials to Github. It is recommended that you use an environment variable for this purpose.")

plants = [x for x in pd.read_csv("plants.csv")["Plants"]]
bees = [x for x in pd.read_csv("bees.csv")["Bee"]]

# drop outgroup species
plants = plants[0:len(plants)-1]
bees = bees[0:len(bees)-1]

splist = bees + plants

bbox = 'POLYGON((-110.5 34,-103.5 34,-103.5 44,-110.5 44,-110.5 34))'

keys = []

for x in splist:
    print(x)
    s = species.name_suggest(x)[0]['key']
    keys += [s]

for k in keys:
    occurrences.search(
                    taxonKey = k, 
                    geometry=bbox,
                    hasCoordinate=True,
                )

query = { 
    "type": "and",
    "predicates": [
            {  
                "type": "in",
                "key": "TAXON_KEY",
                "values": keys
            },
            {
                "type": "within",
                "geometry": 'POLYGON((-110.5 34,-103.5 34,-103.5 44,-110.5 44,-110.5 34))'
            },
            {
                "type": "equals",
                "key": "HAS_COORDINATE",
                "value": True
            },
            {
                "type": "in",
                "key": "BASIS_OF_RECORD",
                "values": ["HUMAN_OBSERVATION", "PRESERVED_SPECIMEN"]
            }
    ]
}

occurrences.download(query,
                    user='<YOUR-USERNAME-HERE>',
                    email="<YOUR-EMAIL-HERE>", 
                    pwd='<YOUR-INFO-HERE>')




plants = [x for x in pd.read_csv("plants.csv")["Plants"]]





