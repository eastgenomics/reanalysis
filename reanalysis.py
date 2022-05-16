#!usr/bin/env python

"""
PART 1. How does the new pipeline compare to what currently exists?

Get case data from OpenCGA:
-Which panels were used (ID, name, version)
-Which phenotypes were reported

Get panel data from PanelApp:
-need genes, regions and STRs used in both original and current versions
-all of these entities can, but don't always, have an associated gene
-might be easier to use genomic coordinates for each entity?
-this would standardise dealing with them all
-it would also make constructing a bed file easier

Get VCF of case variants - Liv is sourcing from GEl

Filtering on the VCF:
-Apply the bed file
-Start with QUAL score == 20?
-If any poor quality variants are present they may be tagged as such
-Need to apply a bed file - construct via BioMart, include 3-5bp padding
either side to ensure splice sites are included

Prioritising filtered variants:
-VEP has an API
-Filter on...
    gnomAD AF (make sure to use genome, not exome)
    known clinical significance
    missense predictors? only makes sense if the variant is missense
-Exomiser - needs phenotype data, can get from primary findings

Compare output with original findings:
-How many variants were originally identified in each tier?
-How has that changed after reanalysis?

Returning output:
-Variants which are known to be pathogenic in the context of the phen

------------------------------------------------------------------------

PART 2. How does the new pipeline compare to another proposed
alternative (TierUp)?

"""


import subprocess
import pandas as pd

# from bs4 import BeautifulSoup

from panelapp import Panelapp
from panelapp import api

from pyopencga.opencga_config import ClientConfiguration
from pyopencga.opencga_client import OpencgaClient


def session_setup(user, host):
    """ Define username and host server to use for OpenCGA login. Log in
    generates password prompt in terminal. Session token then generated and
    used to define OpenCGA client instance for the program.

    args:
        user [string]: username for OpenCGA login
        host [string]: host server URL

    returns:
        oc: OpenCGA client instance based on user credentials and token
    """

    config_dict = {'rest': {'host': host}}
    config = ClientConfiguration(config_dict)

    oc = OpencgaClient(config)
    oc.login(user)

    token = oc.token
    oc = OpencgaClient(configuration=config, token=token)

    return oc


def get_hgnc_id(hgnc_df, gene_symbol):
    """ Get the HGNC ID for a supplied gene symbol, if one exists.

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


""" PanelApp or PanelApp panel object functions """


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
        panel = Panelapp.Panel(str(id)).get_data()

    return panel


def get_all_green_entities(hgnc_df, id, version = None):
    """ Get details of all 'green' genes, regions and STRs (i.e.
    confidence level 3) from a PanelApp panel.

    args:
        hgnc_df: pandas df containing dump of HGNC site
        id [str/int]: PanelApp ID for a panel
        version [str/float] (optional): panel version, defaults to current

    returns:
        panel_version [str]: version of panel in PA
        panel_entities [list of dicts]: green genes, regions and STRs
    """

    panel = get_panelapp_panel(id, version)

    panel_version = panel['version']  # in case version not supplied as an arg
    panel_entities = []

    for type in ['genes', 'regions', 'strs']:
        for entity in panel['{}'.format(type)]:
            if entity['confidence_level'] == '3':

                info = get_entity_dict(hgnc_df, entity)
                panel_entities.append(info)

    return panel_version, panel_entities


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

        info['associated_gene'] = get_hgnc_id(hgnc_df, gene_symbol)

    else:
        info['associated_gene'] = None

    return info


""" Study-level functions """


def write_all_case_ids(oc, study_id):
    """ Write out the IDs of all cases within a study

    args:
        oc: OpenCGA client
        study_id [str]: id specifying project and study
    """

    cases = oc.clinical.search(study = study_id, include = 'id').get_results()

    case_ids = [case['id'] for case in cases]
    case_ids.sort()

    with open('all_case_ids.txt', 'w') as writer:
        for id in case_ids:
            writer.write(id + '\n')


def get_single_case_data(oc, study_id, case_id):
    """ Retrieve all data from specified case in specified study

    args:
        oc: OpenCGA client
        study_id [str]: id specifying project and study
        case_id [str]: id specifying case within that study

    returns:
        case: dict containing all case data
    """

    case = oc.clinical.info(
        study = study_id,
        clinical_analysis = case_id
        ).get_results()[0]

    return case


""" Case-level functions """


def get_case_panels(case, hgnc_dump):
    """ Retrieve details about each panel used in the specified case.
    Includes panel ID and name; version of original panel; list of all
    genes, regions and STRs with confidence level 3 in the original
    version; and the same for the current version.

    args:
        case [dict]: all data from a single case in an OpenCGA study

    returns:
        panels_list [list]: list of panel dicts-
            {'panel_id',
            'panel_name',
            'original_version',
            'original_entities',
            'current_version',
            'current_entities'}
    """

    panels_list = []

    for panel in case['panels']:

        id = panel['source']['id']
        name = panel['name']

        original_version, original_entities = get_all_green_entities(
            hgnc_dump,
            id,
            panel['source']['version'])

        current_version,  current_entities = get_all_green_entities(
            hgnc_dump,
            id)

        sorted_original = sorted(original_entities, key = lambda x: x['name'])
        sorted_current = sorted(current_entities, key = lambda x: x['name'])

        panel_info = {
            'panel_id' : id,
            'panel_name' : name,
            'original_version' : original_version,
            'original_entities' : sorted_original,
            'current_version' : current_version,
            'current_entities' : sorted_current
            }

        panels_list.append(panel_info)

    return panels_list


def get_all_case_intpns(case):
    """ Return a list of all interpretation dicts in a case

    args:
        case [dict]: all data from a single case in an OpenCGA study

    returns:
        all_intpns: list of dicts, each dict is one case interpretation
    """

    all_intpns = []

    all_intpns.append(case['interpretation'])

    for interpretation in case['secondaryInterpretations']:
        all_intpns.append(interpretation)

    return all_intpns


def get_case_phenotypes(case):
    """ Return a list of the phenotypes associated with a case. All
    phenotypes associated with the case are identified by iterating over
    all evidences in all variants in all interpretations.

    args:
        case [dict]: all data from a single case in an OpenCGA study

    returns:
        phenotypes [list]: each element is a str of one phenotype
    """

    phenotypes = []
    interpretations = get_all_case_intpns(case)

    for intpn in interpretations:
        for variant in intpn['primaryFindings']:
            for evidence in variant['evidences']:
                for phenotype in evidence['phenotypes']:

                    if phenotype['id'] not in phenotypes:
                        phenotypes.append(phenotype['id'])

    return phenotypes


""" bed file functions """


def get_hgnc_string(entity_list):
    """

    args:
        entity_list [list]:

    returns:
        hgnc_string [str]:
    """

    hgnc_string = ''

    for entity in entity_list:
        if entity['associated_gene']:

            hgnc = entity['associated_gene'][6:]

            if entity != entity_list[-1]:
                hgnc_string += '{},'.format(hgnc)

            else:
                hgnc_string += '{}'.format(hgnc)

    return hgnc_string


def get_biomart_output(case_id, panel_list):
    """ Using the output of get_case_panels(), create bed files for the
    entities (genes and STRs) covered in the original and current
    versions of the panels in the case.

    Bed files are created by querying Ensembl BioMart (GRCh37) with
    lists of HGNC IDs.

    Can't currently include regions in bed files, because they appear to
    have neither an associated gene nor a GRCh37 location...

    args:
        case_id: ID of the OpenCGA case (for file naming purposes)

        panel_list: output of get_case_panels() - each element is a dict
            containing details on one panel (id, name, original &
            current versions, original & current entities)

    returns:
        list of bed file filename strings
    """

    # given an HGNC ID, how does BioMart turn that into a bed file region?
    # does the region include all exons and introns? what about padding?

    # get lists of all entities in original and current versions of panels

    original = []
    current = []

    for panel in panel_list:

        for entity in panel['original_entities']:
            if entity not in original:
                original.append(entity)

        for entity in panel['current_entities']:
            if entity not in current:
                current.append(entity)

    # define destinations for output bed files

    arg_dicts = [
        {'input' : original,
        'output' : 'data/case_{}_original.bed'.format(case_id)},

        {'input' : current,
        'output' : 'data/case_{}_current.bed'.format(case_id)}
        ]

    # query BioMart API using subprocess for each list of entities

    for dict in arg_dicts:

        hgnc_string = get_hgnc_string(dict['input'])
        filename = dict['output']

        query = (
            'https://grch37.ensembl.org/biomart/martservice?query=<?xml '
            'version="1.0" encoding="utf-8"?><!DOCTYPE Query>'
            '<Query count="" datasetConfigVersion="0.6" formatter="TSV" '
            'header="0" uniqueRows="0" virtualSchemaName="default">'
            '<Dataset interface="default" name="hsapiens_gene_ensembl">'
            '<Filter name="hgnc_id" value="{}"/>'
            '<Attribute name="chromosome_name"/>'
            '<Attribute name="start_position"/>'
            '<Attribute name="end_position"/>'
            '</Dataset></Query>'.format(hgnc_string))

        biomart_call = [
            'wget',
            '-O',
            '{}'.format(filename),
            '{}'.format(query)
            ]

        subprocess.run(biomart_call)

    return [dict['output'] for dict in arg_dicts]


def main():

    user = 'jmiles'  # terminal will prompt for password
    host = 'https://uat.eglh.app.zettagenomics.com/opencga'

    project = 'reanalysis'
    study = 'rd37'
    study_id = 'emee-glh@{}:{}'.format(project, study)

    # case_id = 'SAP-3904-1'  # randomly chosen example case
    case_id = 'OPA-9965-2'  # has multiple panels with regions and STRs

    hgnc_dump = 'data/hgnc_dump_210727.txt'


    ## Set up an OpenCGA session
    oc = session_setup(user, host)

    ## Retrieve all data from a specific case
    case = get_single_case_data(oc, study_id, case_id)

    ## Get info for all panels in a case including original/current genes
    hgnc_df = pd.read_csv(hgnc_dump, sep='\t')
    panel_list = get_case_panels(case, hgnc_df)

    ## Construct bed files for original and current gene lists
    file_list = get_biomart_output(case_id, panel_list)

    ## Get a list of phenotype IDs for the case
    # phenotypes = get_case_phenotypes(case)

    ## Also need a VCF for each case - Liv is trying to get hold of these

    ## Write out current and original panels and entities
    # with open('case_{}_panel_details.txt'.format(case_id), 'w') as writer:
    #     writer.write(str(panel_list))


if __name__ == '__main__':
    main()
