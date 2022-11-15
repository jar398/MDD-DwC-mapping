#!/usr/bin/env python3

import os, csv, argparse, sys, re
from util import find_index, del_list_indexes, generate_hash, generate_taxonID

MDD_DwC_mapping = {
    'Authority_author': 'scientificNameAuthor',   # scientificNameAuthor,namePublishedInYear and authorityParentheses will be combined to create a new column scientificNameAuthorship
    'Authority_fullcitation': 'namePublishedIn',
    'Authority_link': 'authoritySpeciesLink',
    'Authority_sp_author': 'scientificNameAuthor',
    'Authority_sp_fullcitation': 'namePublishedIn',
    'Authority_sp_link': 'authoritySpeciesLink',
    'Authority_sp_year': 'namePublishedInYear',
    'Authority_year': 'namePublishedInYear',
    'authorityParentheses': 'authorityParentheses',
    'authoritySpeciesAuthor': 'scientificNameAuthor',
    'authoritySpeciesCitation': 'namePublishedIn',
    'authoritySpeciesLink': 'authoritySpeciesLink',
    'authoritySpeciesYear': 'namePublishedInYear',
    'biogeographicRealm':'bioGeographicRealm',
    'CMW_sciName': 'scientificNameCMW',
    'countryDistribution': 'countryDistribution',
    'diffSinceCMW': 'diffSinceCMW',
    'diffSinceMSW3': 'diffSinceMSW3',
    'domestic': 'domestic',
    'domestic?': 'domestic', 
    'extinct': 'extinct', 
    'extinct?': 'extinct',
    'Family': 'family',
    'family': 'family', 
    'flagged': 'flagged',
    'flagged?': 'flagged',
    'genus': 'genus', 
    'Genus': 'genus',
    'genusTransfersinceMSW3?': 'genusTransferSinceMSW3',
    'Geo_distribution': 'bioGeographicRealm',
    'Holotype_voucher': 'holotypeVoucher',
    'holotypeVoucher': 'holotypeVoucher',
    'id': 'taxonID',
    'ID_number': 'taxonID',
    'IfNew_category': 'ifNewCategory',
    'IfNew_described_SciName': 'ifNewDescribedScientificName',
    'IfNew_evidenceAuthors': 'ifNewEvidenceAuthors',
    'IfNew_evidenceCitation': 'ifNewEvidenceCitation',
    'IfNew_evidenceLink': 'ifNewEvidenceLink',
    'IfNew_evidenceYear': 'ifNewEvidenceYear',
    'IfNew_GeoRegion': 'ifNewGeoRegion',
    'IfNew_nameAuthors': 'ifNewNameAuthors',
    'IfNew_nameCitation': 'ifNewNameCitation',
    'IfNew_nameLink': 'ifNewNameLink', 
    'IfNew_nameYear': 'ifNewNameYear',
    'IfNew_valid_SciName': 'ifNewValidScientificName',
    'IfTransfer_evidenceCitation': 'ifTransferEvidenceCitation',
    'IfTransfer_oldIUCN_ID': 'ifTransferOldIucnID',
    'IfTransfer_oldSciName': 'ifTransferOldScientificName',
    'infraclass': 'infraclass', 
    'infraorder': 'infraorder',
    'iucnStatus': 'iucnStatus',
    'juniorSynonym?': 'juniorSynonym',
    'magnorder': 'magnorder',
    'mainCommonName': 'vernacularName',
    'majorSubtype': 'superorder',
    'MajorSubtype': 'superorder',
    'MajorType': 'infraclass',
    'majorType': 'infraclass',
    'MDDv1': 'MDDv1',
    'MSW3_matchtype': 'MSW3MatchType',
    'MSW3_SciName': 'MSW3ScientificName',
    'MSW3_sciName': 'MSW3ScientificName',
    'newSppSinceMSW3?': 'newSpeciesSinceMSW3',
    'nominalNames': 'nominalNames',      # this will be removed to create new rows for synonyms, and a new column 'acceptedNameUsageID' will be added 
    'Order': 'order',
    'order': 'order',
    'originalNameCombination': 'originalNameCombination',
    'otherCommonNames': 'otherCommonNames',
    'parvorder': 'parvorder',
    'phylosort': 'phylosort',
    'sciName': 'canonicalName',
    'SciName': 'canonicalName',
    'specific_epithet': 'specificEpithet',
    'specificEpithet': 'specificEpithet',
    'subclass': 'subclass',
    'subfamily': 'subfamily',
    'Subfamily': 'subfamily', 
    'subgenus': 'subgenus',
    'suborder': 'suborder',
    'superfamily': 'superfamily',
    'superorder': 'superorder',
    'synonymOf?': 'synonymOf',
    'TaxonomyNotes': 'taxonRemarks',
    'taxonomyNotes': 'taxonRemarks',
    'TaxonomyNotes_Citation': 'taxonRemarksCitation',
    'taxonomyNotesCitation': 'taxonRemarksCitation',
    'Tribe': 'tribe',
    'tribe': 'tribe',
    'typeLocality': 'typeLocality',
    'typeLocalityLatitude': 'typeLocalityLatitude',
    'typeLocalityLongitude': 'typeLocalityLongitude',
    '\ufeffMajorType': 'majorType',
}

MDD_new_columns = [
    'parentNameUsageID',
    'acceptedNameUsageID',
    'scientificName',
    'scientificNameAuthorship',
    'genericName',
    'taxonRank',
    'taxonomicStatus',
    'nomenclaturalStatus',
    'scientificNameHashKey',
]

taxonID_entries, synonym_entries, subgenus_entries, genus_entries, tribe_entries, subfamily_entries, family_entries, superfamily_entries, parvorder_entries, infraorder_entries, \
suborder_entries, order_entries, superorder_entries, magnorder_entries, infraclass_entries, subclass_entries = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

taxonomic_rank_based_global_entries = {
    'subclass': subclass_entries, 
    'infraclass': infraclass_entries,
    'magnorder': magnorder_entries,
    'superorder': superorder_entries,
    'order': order_entries,
    'suborder': suborder_entries,
    'infraorder': infraorder_entries,
    'parvorder': parvorder_entries,
    'superfamily': superfamily_entries,
    'family': family_entries,
    'subfamily': subfamily_entries,
    'tribe': tribe_entries,
    'genus': genus_entries,
    'subgenus': subgenus_entries
}

def map_MDD_to_DwC(infile, outfile):
  reader = csv.reader(infile)
  in_header = next(reader)
  DwC_header = [MDD_DwC_mapping[col] for col in in_header if col in MDD_DwC_mapping]

  id_index = find_index(DwC_header, 'taxonID')
  author_index = find_index(DwC_header, 'scientificNameAuthor')
  year_index = find_index(DwC_header, 'namePublishedInYear')
  paranthesis_index = find_index(DwC_header, 'authorityParentheses')
  cn_index = find_index(DwC_header, 'canonicalName')
  ep_index = find_index(DwC_header, 'specificEpithet')
  if ep_index == None: ep_index = len(DwC_header)
  genus_index = find_index(DwC_header, 'genus')
  subgenus_index = find_index(DwC_header, 'subgenus')
  nn_index = find_index(DwC_header, 'nominalNames')

  subclass_index = find_index(DwC_header, 'subclass')
  infraclass_index = find_index(DwC_header, 'infraclass')
  magnorder_index = find_index(DwC_header, 'magnorder')
  superorder_index = find_index(DwC_header, 'superorder')
  order_index = find_index(DwC_header, 'order')
  suborder_index = find_index(DwC_header, 'suborder')
  infraorder_index = find_index(DwC_header, 'infraorder')
  parvorder_index = find_index(DwC_header, 'parvorder')
  superfamily_index = find_index(DwC_header, 'superfamily')
  family_index = find_index(DwC_header, 'family')
  subfamily_index = find_index(DwC_header, 'subfamily')
  tribe_index = find_index(DwC_header, 'tribe')

  out_header = [col for col in DwC_header]
  for col in MDD_new_columns:
    out_header.append(col)

  out_pnu_index = find_index(out_header, 'parentNameUsageID')
  out_anu_index = find_index(out_header, 'acceptedNameUsageID')
  out_sn_index = find_index(out_header, 'scientificName')
  out_sna_index = find_index(out_header, 'scientificNameAuthorship')
  out_gn_index = find_index(out_header, 'genericName')
  out_tr_index = find_index(out_header, 'taxonRank')
  out_ts_index = find_index(out_header, 'taxonomicStatus')
  out_ns_index = find_index(out_header, 'nomenclaturalStatus')
  out_snhk_index = find_index(out_header, 'scientificNameHashKey')

  out_rows =[]
  file_taxonID_list = []
  syn_species_dict = {}
  global synonym_entries, subgenus_entries, genus_entries, tribe_entries, subfamily_entries, family_entries, superfamily_entries, parvorder_entries, \
  infraorder_entries, suborder_entries, order_entries, superorder_entries, magnorder_entries, infraclass_entries, subclass_entries, taxonID_entries

  def add_taxon_entry_to_output_list(canonicalName, taxonomicRank, parentNameUsageID, index):
    is_taxon_added_to_output_list = [taxon_entry for taxon_entry in out_rows if taxon_entry[out_tr_index]==taxonomicRank and taxon_entry[cn_index]==canonicalName]
    if is_taxon_added_to_output_list:
      parent_id = is_taxon_added_to_output_list[0][id_index]
    else:
      taxon_entry = create_higher_taxa_entry(canonicalName, taxonomicRank, parentNameUsageID)
      output_taxon_row = [None for col in out_header]
      indexes_to_be_filled = [id_index, cn_index, out_tr_index, out_ts_index, index, out_pnu_index]
      values = [taxon_entry['taxonID'], taxon_entry['canonicalName'], taxon_entry['taxonomicRank'], taxon_entry['taxonomicStatus'], \
              taxon_entry[taxonomicRank], taxon_entry['parentNameUsageID']]
      for indx, value in zip(indexes_to_be_filled, values):
        output_taxon_row[indx] = value
      out_rows.append(output_taxon_row)
      parent_id = taxon_entry['taxonID']
    return parent_id

  for row in reader:
    if len(row) != len(DwC_header):
      print(f'Unexpected number of columns in a row: expected {len(DwC_header)} vs actual {len(row)}', file=sys.stderr)
      print(f'{row}', file=sys.stderr)
      break

    # deals with empty or duplicate taxonID
    row_id = row[id_index]
    if row_id is None or row_id in file_taxonID_list:
      row_id, taxonID_entries = generate_taxonID(taxonID_entries)
    file_taxonID_list.append(row_id)

    if ep_index == len(row):
      row.append(row[cn_index].split('_')[1])

    out_row = row + [None for i in range(len(MDD_new_columns))]
    out_row[id_index] = row_id
    columns_to_capitalize = [i for i in [subclass_index, infraclass_index, magnorder_index, superorder_index, order_index, suborder_index, infraorder_index, \
                  parvorder_index, superfamily_index, family_index, subfamily_index, tribe_index] if i is not None]

    for index in columns_to_capitalize:
      if out_row[index] is not None and out_row[index] != 'NA':
        out_row[index] = out_row[index].capitalize()
    out_row[cn_index] = ' '.join(out_row[cn_index].split('_'))
    if row[genus_index] is not None: out_row[out_gn_index] = row[genus_index] #generic name is same as genus name for accepted name usage
    out_row[out_tr_index] = 'species'
    out_row[out_ts_index] = 'accepted'

    canonical_name = ' '.join(row[cn_index].split('_'))
    author = row[author_index]
    year = row[year_index]
    sn_authorship = f'({author}, {year})' if paranthesis_index is not None and int(row[paranthesis_index]) == 1 else f'{author}, {year}'
    scientific_name = f'{canonical_name} {sn_authorship}'
    out_row[out_sn_index] = scientific_name
    out_row[out_sna_index] = sn_authorship

    epithet = row[ep_index]
    if epithet[-2:] == 'us': 
      epithet = epithet[:-2]+'x'
    elif epithet[-1:] == 'a':
      epithet = epithet[:-1]+'x'
    sn_hash_key = generate_hash([epithet, author, year])
    out_row[out_snhk_index] = sn_hash_key

    parent_id = None

    taxa = [
        [subclass_index, 'subclass'],
        [infraclass_index, 'infraclass'],
        [magnorder_index, 'magnorder'],
        [superorder_index, 'superorder'],
        [order_index, 'order'],
        [suborder_index, 'suborder'],
        [infraorder_index, 'infraorder'],
        [parvorder_index, 'parvorder'],
        [superfamily_index, 'superfamily'],
        [family_index, 'family'],
        [subfamily_index, 'subfamily'],
        [tribe_index, 'tribe'],
        [genus_index, 'genus'],
        [subgenus_index, 'subgenus'],
      ]

    for taxa_details in taxa:
      if taxa_details[0] is not None and row[taxa_details[0]] is not None and row[taxa_details[0]] != 'NA':
        parent_id = add_taxon_entry_to_output_list(row[taxa_details[0]].capitalize(), taxa_details[1], parent_id, taxa_details[0])

    out_row[out_pnu_index] = parent_id
    out_rows.append(out_row)

    
    if nn_index is not None and row[nn_index] is not None:
      syn_list = row[nn_index].split('|')
      unique_syn_list = [] # to check multiple entries of the same synonym for a species

      for syn in syn_list:
        if syn.strip() == '': continue
        syn_epithet, syn_authorship = syn.split(' ', 1)
        if syn_epithet == row[ep_index]: continue

        if re.search('\[(.*?)[\]\)]$', syn_authorship):
          syn_authorship, nomenclatural_status = syn_authorship.split('[')
          if syn_authorship is None: 
            syn_authorship = ''
          else:
            syn_authorship = syn_authorship[:-1]
          nomenclatural_status = nomenclatural_status[:-1]
        elif re.search('\[', syn_authorship):
          assert False, f'nominal name - {syn} - formating error'
        else:
          nomenclatural_status = ''

        if f'{syn_epithet} {syn_authorship}' in unique_syn_list:
          print(f'row with duplicate synonyms: {out_row[out_sn_index]}')
          continue
        else:
          unique_syn_list.append(f'{syn_epithet} {syn_authorship}')
        
        scientific_name = f'{row[genus_index]} {syn_epithet} {syn_authorship}'
        canonical_name = f'{row[genus_index]} {syn_epithet}'

        if scientific_name in syn_species_dict:
          syn_species_dict[scientific_name].append(out_row[cn_index])
          continue
        else:
          syn_species_dict[scientific_name] = [out_row[cn_index]]

        #Find if a synonym is junior or senior synonym
        find_syn_year = re.match(r'.*([1-3][0-9]{3})', syn_authorship)
        if find_syn_year is not None and row[year_index] is not None:
          syn_year = find_syn_year.group(1)
          if int(syn_year) < int(row[year_index]):
            ts_status = 'senior synonym'
          elif int(syn_year) > int(row[year_index]):
            ts_status = 'junior synonym'
          else:
            ts_status = 'synonym'
        else:
          ts_status = 'synonym'

        exist_synonym = is_synonym_row_available(scientific_name, row_id)
        syn_row = [None for col in out_header]
        syn_indexes = [id_index, out_pnu_index, out_anu_index, out_sn_index, out_sna_index, cn_index, out_gn_index, ep_index, out_tr_index, out_ts_index, out_ns_index, genus_index]
        if not exist_synonym:
          syn_taxon_id, taxonID_entries = generate_taxonID(taxonID_entries)
        else:
          syn_taxon_id = exist_synonym['taxonID']
        syn_values = [syn_taxon_id, out_row[out_pnu_index], row_id, scientific_name, syn_authorship, canonical_name, row[genus_index], syn_epithet, \
                'species', ts_status, nomenclatural_status, row[genus_index]]
        for index, value in zip(syn_indexes, syn_values):
          syn_row[index] = value

        out_rows.append(syn_row)

  print("===============================================")
  for k, v in syn_species_dict.items():
    if len(v)>1:
      print(f'{k}  -  {", ".join(v)}')
  writer = csv.writer(outfile)
  indexes_to_delete = {author_index, year_index, paranthesis_index, nn_index}
  out_header = del_list_indexes(out_header, indexes_to_delete)
  writer.writerow(out_header)
  for row in out_rows:
    trimmed_row = del_list_indexes(row, indexes_to_delete)
    writer.writerow(trimmed_row)

def is_taxon_available(name, taxonomic_rank):
  gobal_entries_to_search = taxonomic_rank_based_global_entries[taxonomic_rank]
  if len(gobal_entries_to_search) == 0: return False

  for entry in gobal_entries_to_search:
    if name == entry['canonicalName']:
      return entry

  return False

def is_synonym_row_available(synonym_scientific_name, anu_id):
  if len(synonym_entries)==0: return False

  for synonym_row in synonym_entries:
    if synonym_scientific_name==synonym_row['scientificName'] and anu_id==synonym_row['acceptedNameUsageID']:
      return synonym_row

  return False

def create_higher_taxa_entry(canonicalName, taxonomicRank, parentNameUsageID):
  global genus_entries, tribe_entries, subfamily_entries, family_entries, superfamily_entries, parvorder_entries, subgenus_entries, \
  infraorder_entries, suborder_entries, order_entries, superorder_entries, magnorder_entries, infraclass_entries, subclass_entries, taxonID_entries

  exist_taxon =  is_taxon_available(canonicalName, taxonomicRank)
  if exist_taxon:
    taxon_id = exist_taxon['taxonID']
  else:
    taxon_id, taxonID_entries = generate_taxonID(taxonID_entries)

  taxon_entry = {'taxonID': taxon_id, 'canonicalName': canonicalName, 'taxonomicRank': taxonomicRank, 'taxonomicStatus': 'accepted',\
              taxonomicRank: canonicalName, 'parentNameUsageID': parentNameUsageID}
  taxonomic_rank_based_global_entries[taxonomicRank].append(taxon_entry)

  return taxon_entry

def managed_entries_across_MDD(outpath):
  global synonym_entries, subgenus_entries, genus_entries, tribe_entries, subfamily_entries, family_entries, superfamily_entries, parvorder_entries, \
  infraorder_entries, suborder_entries, order_entries, superorder_entries, magnorder_entries, infraclass_entries, subclass_entries

  for file_ in os.listdir(outpath):
    if not file_.startswith('.')  and os.path.isfile(os.path.join(outpath, file_)):
      infile = os.path.join(outpath, file_)
      with open(infile) as f:
        print("Reading %s" % infile, file=sys.stderr)
        reader = csv.DictReader(f)
        for row in reader:
          if 'taxonRank' not in row:
            break
          else:
            if row['taxonRank']=='species' and row['taxonomicStatus'] in ['junior synonym', 'senior synonym', 'synonym']:
              if any([True for dict_ in synonym_entries if dict_['scientificName']==row['scientificName']]):
                continue
              else:
                synonym_entries.append(row)
            elif row['taxonRank']=='subgenus':
              if all([True for dict_ in subgenus_entries if row['canonicalName'] != dict_['canonicalName']]):
                subgenus_entries.append(row)
              else:
                continue
            elif row['taxonRank']=='genus':
              if all([True for dict_ in genus_entries if row['canonicalName'] != dict_['canonicalName']]):
                genus_entries.append(row)
              else:
                continue
            elif row['taxonRank']=='tribe':
              if all([True for dict_ in tribe_entries if row['canonicalName'] != dict_['canonicalName']]):
                tribe_entries.append(row)
              else:
                continue
            elif row['taxonRank']=='subfamily':
              if all([True for dict_ in subfamily_entries if row['canonicalName'] != dict_['canonicalName']]):
                subfamily_entries.append(row)
              else:
                continue
            elif row['taxonRank']=='family':
              if all([True for dict_ in family_entries if row['canonicalName'] != dict_['canonicalName']]):
                family_entries.append(row)
              else:
                continue
            elif row['taxonRank']=='superfamily':
              if all([True for dict_ in superfamily_entries if row['canonicalName'] != dict_['canonicalName']]):
                superfamily_entries.append(row)
              else:
                continue
            elif row['taxonRank']=='parvorder':
              if all([True for dict_ in parvorder_entries if row['canonicalName'] != dict_['canonicalName']]):
                parvorder_entries.append(row)
              else:
                continue
            elif row['taxonRank']=='infraorder':
              if all([True for dict_ in infraorder_entries if row['canonicalName'] != dict_['canonicalName']]):
                infraorder_entries.append(row)
              else:
                continue
            elif row['taxonRank']=='suborder':
              if all([True for dict_ in suborder_entries if row['canonicalName'] != dict_['canonicalName']]):
                suborder_entries.append(row)
              else:
                continue
            elif row['taxonRank']=='order':
              if all([True for dict_ in order_entries if row['canonicalName'] != dict_['canonicalName']]):
                order_entries.append(row)
              else:
                continue
            elif row['taxonRank']=='superorder':
              if all([True for dict_ in superorder_entries if row['canonicalName'] != dict_['canonicalName']]):
                superorder_entries.append(row)
              else:
                continue
            elif row['taxonRank']=='magnorder':
              if all([True for dict_ in magnorder_entries if row['canonicalName'] != dict_['canonicalName']]):
                magnorder_entries.append(row)
              else:
                continue
            elif row['taxonRank']=='infraclass':
              if all([True for dict_ in infraclass_entries if row['canonicalName'] != dict_['canonicalName']]):
                infraclass_entries.append(row)
              else:
                continue
            elif row['taxonRank']=='subclass':
              if all([True for dict_ in subclass_entries if row['canonicalName'] != dict_['canonicalName']]):
                subclass_entries.append(row)
              else:
                continue
            else:
              continue

def taxonIDs_across_MDD(inpath, outpath):
  taxonIDs = set()

  for file_ in os.listdir(inpath):
    if not file_.startswith('.')  and os.path.isfile(os.path.join(inpath, file_)):
      infile = os.path.join(inpath, file_)
      with open(infile) as f:
        print(". Reading %s" % infile, file=sys.stderr)
        reader = csv.DictReader(f)
        for row in reader:
          id_ = [row[key] for key in row if key in ['id', 'taxonID', 'ID_number']]
          if '' in id_: id_.remove('')
          if id_:
            taxonIDs.add(int(id_[0]))

  for file_ in os.listdir(outpath):
    if not file_.startswith('.')  and os.path.isfile(os.path.join(outpath, file_)):
      infile = os.path.join(outpath, file_)
      with open(infile) as f:
        print(". Reading %s" % infile, file=sys.stderr)
        reader = csv.DictReader(f)
        for row in reader:
          if row['taxonID']:
            taxonIDs.add(int(row['taxonID']))
          else:
            print("Missing taxonID: %s" % row['canonicalName'], file=sys.stderr)
  return list(taxonIDs)
  
if __name__=='__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--input', default='MDD_v1_6495species_JMamm.csv', \
            help='name of the input file to be mapped to DwC format. The file must be in CSV format')
  parser.add_argument('--output', help='name of the output file')
  args=parser.parse_args()

  infile = args.input
  outfile = args.output
  if outfile is None: outfile = f'{args.input.rsplit(".", 1)[0]}_inDwC.csv'

  JAR = True
  if JAR:
    qual_infile = infile
    qual_outfile = outfile
    input_dir_path = os.path.dirname(infile)
    output_dir_path = os.path.dirname(outfile)
  else:
    input_dir_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../data')
    output_dir_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../MDD_DwC_versions')
    qual_infile = os.path.join(input_dir_path, infile)
    qual_outfile = os.path.join(output_dir_path, outfile)

  managed_entries_across_MDD(output_dir_path)
  taxonID_entries = taxonIDs_across_MDD(input_dir_path, output_dir_path)

  with open(qual_infile) as f1, open(qual_outfile, 'w') as f2:
    map_MDD_to_DwC(f1, f2)
