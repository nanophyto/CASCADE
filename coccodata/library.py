import pandas as pd
import numpy as np
import sys
from collections import namedtuple
from yaml import load, Loader
import yaml
import math
import glob, os

class library:

    def __init__(self, phases, family, sizes):
        self.phases_path = phases
        self.family_path = family
        self.sizes_path = sizes
        self.Study = namedtuple("Study", ['id', 'mean', 'sd', 'method', 'species', 'genera', 'family', 'phase', 'alternate_phase'])


        def def_grouping():

            with open(self.phases_path, 'r') as f:
                phases = load(f, Loader=Loader)

            with open(self.family_path, 'r') as f:
                families = load(f, Loader=Loader)

            d = pd.DataFrame.from_dict(phases, orient='index')
            d = d.rename_axis("species").reset_index()
            d['genera'] = d['species'].str.split(" ").str[0]

            inverse = {}
            for k,v in families.items():
                for x in v:
                    inverse.setdefault(x, []).append(k)

            df = pd.DataFrame.from_dict(inverse, orient='index')
            df = df.rename_axis("genera").reset_index()
            df = df.rename(columns={0: "family"})
            d = pd.merge(d, df, on='genera', how="outer")
            d = d.where(pd.notnull(d), None)

            library = list(d.itertuples(name='species', index=False))

            return(d)

        groups = def_grouping()

        species = groups.species
        species = {x for x in species if x is not None}
        species_list = list(species)
        species_list.sort()
        self.species_list = species_list

        def import_data(path):

            all_files = glob.glob(os.path.join(path, "*.csv"))

            d = pd.concat((pd.read_csv(f) for f in all_files), ignore_index=True)
            d = d.fillna(0)

            return(d)

        d = import_data(self.sizes_path)

        self.list_of_studies = d['reference'].unique()

        def size_and_method(d, study, species):

            d = d[(d['species'] == species) & (d['reference']==study)]

            keys_to_extract = ['mean', 'sd', 'method']
            extracted_data = {}

            for key in keys_to_extract:
                try:
                    extracted_data[key] = d.get(key, None).item()
                except:
                    extracted_data[key] = None
            return extracted_data['mean'], extracted_data['sd'], extracted_data['method']

        def classification(groups, species):

            groups = groups[groups['species'] == species]
            keys_to_extract = ['genera', 'family', 'phase', 'alternate_phase']
            extracted_data = {}

            for key in keys_to_extract:
                try:
                    extracted_data[key] = groups.get(key, None).item()
                except:
                    extracted_data[key] = None

            return extracted_data['genera'], extracted_data['family'], extracted_data['phase'], extracted_data['alternate_phase']

        def fill_namedtuple(groups, species_list):

            studies = []
            for id in self.list_of_studies:
                for species in species_list:
                    genera, family, phase, alternate_phase = classification(groups, species)
                    studies.append(self.Study(id, *size_and_method(d, id, species), species, genera, family, phase, alternate_phase))

            return(studies)

        self.library = fill_namedtuple(groups, species_list)
    
    def return_ntpl(self):
        return(self.library)

    def return_species_list(self):
        return(self.species_list)
    
    def export_yml(self, path):
        spp_list = []
        #sort species list alphabetically:
        for i in range(len(self.species_list)):

            name = self.species_list[i]
            species_library =  [t for t in self.library  if t.species == name]

            sizes = {study.id:study._asdict() for study in species_library }
            for id in self.list_of_studies:
                del sizes[id]['id']
                del sizes[id]['species']
                del sizes[id]['genera']
                del sizes[id]['family']
                del sizes[id]['phase']
                del sizes[id]['alternate_phase']

            #d = asdict(dc, dict_factory=lambda x: {k: v for (k, v) in x if v is not None})


            species =  {name: {
                    'genera': species_library[0].genera,
                    'family': species_library[0].family,
                'phase': species_library[0].phase,
                    'alternate_phase': species_library[0].alternate_phase,
                    'size' : sizes
                }}
            spp_list.append(species)

        with open(path, 'w') as outfile:
            yaml.dump(spp_list, outfile, default_flow_style=False)

        print("exported yml to: " + str(path))



# def find_undefined_spp(library):
#     for i in range(len(library)):
#         ntpl = library[i]
#         print(ntpl.mean)

#         if (all(ntpl.mean) is None):
#                 print(str(ntpl.species))


# find_undefined_spp(library)

print("fin")