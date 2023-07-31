The `data` directory is organized as follows.

```
data
├── embargo
│   └── interactions
│       ├── raw                   # Required to start
│       │   ├─── pikespeak.csv
│       │   ├─── gothic.csv
│       │   └─── elkmeadows.csv
│       └── clean
│           ├─── pikespeak.csv
│           ├─── gothic.csv
│           └─── elkmeadows.csv
├── public
│   ├── occurence                 # Required to start
│   │   └─── occurrence.csv 
│   ├── phylogeny                 
│   │   ├── raw_sequences         # Required to start
│   │   │   ├─── bee_sequences.csv
│   │   │   └─── plant_sequences.csv
│   │   └── aligned_sequences     # Required to start
│   │       ├─── bee_alignment.nex
│   │       └─── plant_alignment.nex
│   ├── species_lists             # Required to start
│   │  ├─── bees.csv
│   │  └─── plants.csv
│   └── chelsa            
│       ├─── raw                  # Required to start
│       └─── pca
├── artifacts
│   ├── sdms                
│   ├── species_representations                
│   ├── interaction_prediction    
│   ├── projected_overlap                
│   └── interaction_richness
└── misc
    └── geojsons
        ├─── co               
        │    ├─── counties.json                 
        │    └─── states.json 
        └─── usa
             ├─── counties.json
             └─── states.json