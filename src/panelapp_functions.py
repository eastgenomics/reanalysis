#!usr/bin/env python


import pandas as pd

from panelapp import Panelapp
from panelapp import api
from panelapp import queries


""" Functions over all PA panels """


def write_all_attribute_values():
    """ Iterate over every current panel version in PA and identify all
    possible values of the 'mode_of_inheritance', 'penetrance' and
    'mode_of_pathogenicity' gene/region/STR attributes. """

    panels = queries.get_all_panels()

    confs = []  # confidence level
    mois = []  # mode of inheritance
    mops = []  # mode of pathogenicity
    pens = []  # penetrance

    str_names = []  # entity name
    str_seqs = []  # repeated sequence
    str_normals = []  # normal no. repeats
    str_paths = []  # pathogenic no. repeats

    region_names = []  # entity name
    region_haplos = []  # haploinsufficiency score
    region_triplos = []  # triplosensitivity score
    region_overlaps = []  # required overlap percentage
    region_types = []  # variant types

    for id, panel_object in panels.items():

        panel = panel_object.get_data()

        # 4 properties are common to multiple entity types
        for entity_type in [panel['genes'], panel['strs'], panel['regions']]:
            for entity in entity_type:

                if entity['confidence_level'] not in confs:
                    confs.append(entity['confidence_level'])

                if entity['mode_of_inheritance'] not in mois:
                    mois.append(entity['mode_of_inheritance'])

                if entity['penetrance'] not in pens:
                    pens.append(entity['penetrance'])

                try:
                    if entity['mode_of_pathogenicity'] not in mops:
                        mops.append(entity['mode_of_pathogenicity'])

                except KeyError:  # STRs don't have a mode of pathogenicity
                    pass

        # Properties specific to STRs
        for entity in panel['strs']:

            if entity['entity_name'] not in str_names:
                str_names.append(entity['entity_name'])

            if entity['repeated_sequence'] not in str_seqs:
                str_seqs.append(entity['repeated_sequence'])

            if entity['normal_repeats'] not in str_normals:
                str_normals.append(entity['normal_repeats'])

            if entity['pathogenic_repeats'] not in str_paths:
                str_paths.append(entity['pathogenic_repeats'])

        # Properties specific to regions
        for entity in panel['regions']:

            if entity['entity_name'] not in region_names:
                region_names.append(entity['entity_name'])

            if entity['haploinsufficiency_score'] not in region_haplos:
                region_haplos.append(entity['haploinsufficiency_score'])

            if entity['triplosensitivity_score'] not in region_triplos:
                region_triplos.append(entity['triplosensitivity_score'])

            if entity['required_overlap_percentage'] not in region_overlaps:
                region_overlaps.append(entity['required_overlap_percentage'])

            if entity['type_of_variants'] not in region_types:
                region_types.append(entity['type_of_variants'])

    all_lists = [
        confs,
        mops,
        pens,
        str_names,
        str_seqs,
        str_normals,
        str_paths,
        region_names,
        region_haplos,
        region_triplos,
        region_overlaps,
        region_types]

    for list in all_lists:
        try:
            list.sort()

        except TypeError:
            pass

    # Dump results out to text
    with open('search_all_panels.txt', 'w') as writer:
        writer.write('Possible confidence values:\n\n{}\n'.format(str(confs)))
        writer.write('\nPossible MOI values:\n\n{}\n'.format(str(mois)))
        writer.write('\nPossible MOP values:\n\n{}\n'.format(str(mops)))
        writer.write('\nPossible penetrance values:\n\n{}\n'.format(str(pens)))

        writer.write('\nSTRs\n')
        writer.write('\nPossible entity names:\n\n{}\n'.format(str(str_names)))
        writer.write('\nPossible sequences:\n\n{}\n'.format(str(str_seqs)))
        writer.write('\nPossible normal counts:\n\n{}\n'.format(str(str_normals)))
        writer.write('\nPossible pathogenic counts:\n\n{}\n'.format(str(str_paths)))

        writer.write('\nRegions\n')
        writer.write('\nPossible entity names:\n\n{}\n'.format(str(region_names)))
        writer.write('\nPossible haplo scores:\n\n{}\n'.format(str(region_haplos)))
        writer.write('\nPossible triple scores:\n\n{}\n'.format(str(region_triplos)))
        writer.write('\nPossible overlap percents:\n\n{}\n'.format(str(region_overlaps)))
        writer.write('\nPossible variant types:\n\n{}\n'.format(str(region_types)))


def write_pa_objects_list():
    """ Write out a list of all panel IDs and their associated panel
    objects. 'panels' is a dict where each key is a panel ID and each
    value is a panel object. """

    panels = queries.get_all_panels()

    with open('pa_panel_objects_list.txt', 'w') as writer:
        for key, value in panels.items():
            writer.write('{}: {}\n'.format(key, value))


def write_panel_entity_counts():
    """ Write out the number of genes, regions and STRs covered by each
    PanelApp panel. """

    panels = queries.get_all_panels()

    smallest_panel = ''
    smallest_panel_size = 10000

    with open('feature_counts.txt', 'w') as writer:
        writer.write('Number of genes, regions and STRs in each PA panel\n')

    for id, panel_object in panels.items():

        panel = panel_object.get_data()

        non_empty_features = 0  # range 0-3 for genes, regions & strs
        features_sum = 0  # sum of all genes, regions & strs values

        with open('feature_counts.txt', 'a') as writer:

            writer.write('\nPanel {}: {}\n'.format(panel['id'], panel['name']))

            for key, value in panel['stats'].items():
                writer.write('{}: {}\n'.format(key, value))

                features_sum += value

                if value > 0:
                    non_empty_features += 1

            # find smallest panel where genes/regions/strs are all non-empty
            if (non_empty_features == 3) and\
                (features_sum < smallest_panel_size):

                smallest_panel = '{}: {}'.format(panel['id'], panel['name'])
                smallest_panel_size = features_sum

    print('The smallest panel where genes, regions and STRs are all \
        non-empty is: {}'.format(smallest_panel))


def write_grs_panels_list():
    """ Write out a list of PA panels which contain at least one gene,
    region and STR. """

    panels = queries.get_all_panels()

    grs_panels = []

    for id, panel_object in panels.items():

        panel = panel_object.get_data()

        if (panel['genes'] != []) and \
            (panel['regions'] != []) and \
            (panel['strs'] != []):

            info = {
                'id' : panel['id'],
                'name' : panel['name'],
                'version' : panel['version'],
                'gene_count' : len(panel['genes']),
                'region_count' : len(panel['regions']),
                'str_count' : len(panel['strs'])
                }

            grs_panels.append(info)

    with open('pa_grs_panels.txt', 'w') as writer:
        for panel in grs_panels:
            for key, value in panel.items():
                writer.write('{}: {}\n'.format(key, value))


def print_all_strs():
    """ Identify the number of PA panels which cover at least one STR.
    For each STR covered by any panel, print out its name and the gene
    it affects. """

    panels = queries.get_all_panels()

    all_strs = []
    panels_with_strs = 0
    total_strs = 0

    for id, panel_object in panels.items():
        panel = panel_object.get_data()

        if len(panel['strs']) > 0:
            panels_with_strs += 1

            for repeat in panel['strs']:
                total_strs += 1

                info = (
                        repeat['entity_name'],
                        repeat['gene_data']['gene_symbol'])

                if info not in all_strs:
                    all_strs.append(info)

    print('{} panels have >0 STRs, panels cover {} total STRs'.format(
        panels_with_strs,
        total_strs))

    for item in all_strs:
        print(item)


""" Single-panel functions """


def get_panelapp_panel(id, version = None):
    """ Retrieve panel object representing specified version of a
    PanelApp panel ('Panelapp.Panel' doesn't always work, because some
    older panel versions don't contain the hgnc_symbol element)

    args:
        id [str/int]: the panel's PanelApp ID

        version [str/float] (optional): version of the panel to use - if
            not supplied, retrieves current version

    returns:
        panel: dict of PanelApp data on specified version of panel
    """

    if version:
        path = ["panels", str(id)]
        param = {"version": str(version)}

        url = api.build_url(path, param)
        panel = api.get_panelapp_response(url)

    else:
        panel = Panelapp.Panel(str(id))
        panel.get_data()

    return panel


def write_pa_panel(id, version = None):
    """ Dump out the contents of a PanelApp panel to a text file

    args:
        id [str/int]: PanelApp ID for a panel
        version [str/float] (optional): specific panel version to get
    """

    panel = get_panelapp_panel(id, version)

    filename = 'pa_dump_panel_{}_version_{}.txt'.format(id, panel['version'])

    with open(filename, 'w') as writer:
        writer.write(str(panel))


def get_all_green_entities(hgnc_df, id, version = None):
    """ Get details of all 'green' genes, regions and STRs (confidence
    level 3) from a PanelApp panel. Also returns panel version because
    this is otherwise not known for current panel versions (since
    function is called with 'version=None').

    args:
        hgnc_df: pandas df containing dump of HGNC site
        id [str/int]: PanelApp ID for a panel
        version [str/float] (optional): panel version, defaults to current

    returns:
        panel_version [str]: version of panel in PA
        panel_entities [list of dicts]: green genes, regions and STRs
    """

    panel = get_panelapp_panel(id, version)

    panel_version = panel['version']
    panel_entities = []

    for type in ['genes', 'regions', 'strs']:
        for entity in panel['{}'.format(type)]:
            if entity['confidence_level'] == '3':

                info = get_entity_dict(hgnc_df, entity)
                panel_entities.append(info)

    return panel_version, panel_entities


def get_green_gene_hgncs(hgnc_df, panel_object):
    """ Get a list of the green genes in a PanelApp panel

    args:
        hgnc_df: pandas df containing dump of HGNC site
        panel_object: PanelApp panel object

    returns:
        panel_genes: list of dicts, one for each green panel gene
            {entity_type,
            entity_name,
            associated_gene,
            chromosome,
            grch37_coordinates}
    """

    panel_genes = []

    for gene in panel_object['genes']:
        if gene['confidence_level'] == '3':

            try:
                gene_id = gene['gene_data']['hgnc_id']

            except KeyError:

                gene_symbol = gene['gene_symbol']

                gene_id = get_entity_hgnc(hgnc_df, gene_symbol)

            if gene_id not in panel_genes:
                panel_genes.append(gene_id)

    return panel_genes


def compare_green_genes(id_1, id_2, version_1 = None, version_2 = None):
    """ Compare the lists of green genes in two PanelApp panel objects
    (can be two versions of the same panel)

    args:
        pa_panel_1: dict representing a PanelApp panel
        pa_panel_2: dict representing a PanelApp panel

    returns:
        bool for 'the two gene lists are identical'
    """

    panel_1 = get_panelapp_panel(id_1, version_1)
    panel_2 = get_panelapp_panel(id_2, version_2)

    gene_lists = []

    for panel in [panel_1, panel_2]:

        genes = get_green_gene_hgncs(panel)
        genes.sort()
        gene_lists.append(genes)

    if gene_lists[0] == gene_lists[1]:
        return True

    else:
        return False


""" Single-entity functions """


def get_entity_dict(hgnc_df, entity):
    """ Get a dict of details about a gene/region/STR from a PA panel

    args:
        hgnc_df: pandas df containing dump of HGNC site
        entity [dict]: single gene/region/STR from a PA panel object

    returns:
        info [dict]: details about that entity
            {type,
            name,
            chromosome,
            grch37_coordinates,
            associated_gene}
    """

    # Initialise output dict

    info = {'type' : entity['entity_type']}

    # Get human-readable entity name (different for regions)

    if entity['entity_type'] == 'region':
        info['name'] = entity['verbose_name']

    else:
        info['name'] = entity['entity_name']

    # Get entity's genomic coordinates (different for genes)

    if entity['entity_type'] == 'gene':

        try:
            location = entity['gene_data']['ensembl_genes']['GRch37']['82'][\
                'location']

            split_location = location.split(':')
            split_coordinates = split_location[1].split('-')

            info['chromosome'] = split_location[0]

            info['grch37_coordinates'] = [
                int(split_coordinates[0]),
                int(split_coordinates[1])]

        except KeyError:

            info['chromosome'] = None
            info['grch37_coordinates'] = None

    else:
        info['chromosome'] = entity['chromosome']
        info['grch37_coordinates'] = entity['grch37_coordinates']

    # Get HGNC ID for entity's associated gene (if there is one)

    if entity['gene_data'] != None:

        gene_symbol = entity['gene_data']['gene_symbol']

        info['associated_gene'] = get_entity_hgnc(hgnc_df, gene_symbol)

    else:
        info['associated_gene'] = None

    return info


def get_entity_hgnc(hgnc_df, gene_symbol):
    """ Get the HGNC IF for a supplied gene symbol, if one exists.

    args:
        hgnc_df: pandas df containing dump of HGNC site
        gene_symbol [str]: any gene symbol

    returns:
        hgnc_id [str], or None: HGNC ID of gene associated with entity
    """

    hgnc_id = None

    # If a row exists where this gene is official gene symbol,
    try:

        # Get that row's index and associated HGNC ID
        target_index = hgnc_df.index[hgnc_df['symbol'] == gene_symbol]
        hgnc_id = hgnc_df.loc[target_index[0], 'hgnc_id']

    # If this gene isn't an official symbol,
    except IndexError:

        i = 0

        # Look through 'previous official gene symbols' column
        for value in hgnc_df['prev_symbol']:

            # If gene appears in this column, get that row's HGNC ID
            if gene_symbol in str(value):

                hgnc_id = hgnc_df.iloc[i].loc['hgnc_id']
                break

            i += 1

    return hgnc_id


def main():

    panel_id = 90
    version = 1.5

    # hgnc_dump = 'data/hgnc_dump_210727.txt'
    # hgnc_df = pd.read_csv(hgnc_dump, sep='\t')

    # write_all_attribute_values()
    # write_pa_objects_list()
    # write_grs_panels_list()
    # write_panel_entity_counts()
    # print_all_strs()

    # get_panelapp_panel(id, version)
    # write_pa_panel(id, version)
    # get_all_green_entities(hgnc_df, id, version)
    # get_green_gene_hgncs(hgnc_df, panel_object)
    # compare_green_genes(id_1, id_2, version_1, version_2)

    # get_entity_dict(hgnc_df, entity)
    # get_entity_hgnc(hgnc_df, gene_symbol)


if __name__ == '__main__':
    main()
