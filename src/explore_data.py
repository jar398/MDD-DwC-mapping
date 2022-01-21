#!/usr/bin/env python3

import csv, argparse, sys, re
from util import find_index, del_list_indexes, generate_hash, generate_taxonID

MDD_vers= [
		'MDD_v1_6495species_JMamm.csv', 
		'MDD_v1.1_6526species.csv', 
		'MDD_v1.2_6485species.csv', 
		'MDD_v1.3_6513species.csv', 
		'MDD_v1.31_6513species.csv', 
		'MDD_v1.4_6533species.csv', 
		'MDD_v1.5_6554species.csv',
		'MDD_v1.6_6557species.csv', 
		# 'MDD_v1.7_6567species.csv',
		# 'MDD_v1_6495species_JMamm_inDwC.csv', 
		# 'MDD_v1.1_6526species_inDwC.csv', 
		# 'MDD_v1.2_6485species_inDwC.csv', 
		# 'MDD_v1.3_6513species_inDwC.csv', 
		# 'MDD_v1.31_6513species_inDwC.csv', 
		# 'MDD_v1.4_6533species_inDwC.csv', 
		# 'MDD_v1.5_6554species_inDwC.csv',
		# 'MDD_v1.6_6557species_inDwC.csv', 
		'MDD_v1.7_6567species_inDwC.csv',
]

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
		'subgenus': 'subgenus', # update to include genus name as per DwC?
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
	count = 0

	def add_entry_to_output(canonicalName, taxonomicRank, parentNameUsageID, index):
		added_to_output_list = [out_row for out_row in out_rows if out_row[out_tr_index]==taxonomicRank and out_row[cn_index]==canonicalName]
		if added_to_output_list:
			parent_id = added_to_output_list[0][id_index]
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
		sn_authorship = f'({author}, {year})' if paranthesis_index is not None and row[paranthesis_index] == 1 else f'{author}, {year}'
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
				parent_id = add_entry_to_output(row[taxa_details[0]].capitalize(), taxa_details[1], parent_id, taxa_details[0])

		# # Create new row for subclass if it doesn't already exist
		# if subclass_index is not None and row[subclass_index] is not None and row[subclass_index] != 'NA':
		# 	parent_id = add_entry_to_output(row[subclass_index], 'subclass', parent_id, subclass_index)
		# 	subclass_in_out_rows = [out_row for out_row in out_rows if out_row[out_tr_index]=='subclass' and out_row[cn_index]==row[subclass_index]]
		# 	if subclass_in_out_rows:
		# 		parent_id = subclass_in_out_rows[0][id_index]
		# 	else:
		# 		subclass_entry = create_higher_taxa_entry(row[subclass_index], 'subclass', parent_id)
		# 		out_subclass_row = [None for col in out_header]
		# 		indexes_to_be_filled = [id_index, cn_index, out_tr_index, out_ts_index, subclass_index, out_pnu_index]
		# 		values = [subclass_entry['taxonID'], subclass_entry['canonicalName'], subclass_entry['taxonomicRank'], subclass_entry['taxonomicStatus'], \
		# 					subclass_entry['subclass'], subclass_entry['parentNameUsageID']]
		# 		for index, value in zip(indexes_to_be_filled, values):
		# 			out_subclass_row[index] = value
		# 		out_rows.append(out_subclass_row)
		# 		parent_id = subclass_entry['taxonID']

		# # create new row for infraclass if it doesn't already exist
		# if infraclass_index is not None and row[infraclass_index] is not None and row[infraclass_index] != 'NA':
		# 	infraclass_in_out_rows = [out_row for out_row in out_rows if out_row[out_tr_index]=='infraclass' and out_row[cn_index]==row[infraclass_index]]
		# 	if infraclass_in_out_rows:
		# 		parent_id = infraclass_in_out_rows[0][id_index]
		# 	else:
		# 		infraclass_entry = create_higher_taxa_entry(row[infraclass_index], 'infraclass', parent_id)
		# 		out_infraclass_row = [None for col in out_header]
		# 		indexes_to_be_filled = [id_index, cn_index, out_tr_index, out_ts_index, infraclass_index, out_pnu_index]
		# 		values = [infraclass_entry['taxonID'], infraclass_entry['canonicalName'], infraclass_entry['taxonomicRank'], infraclass_entry['taxonomicStatus'], \
		# 					infraclass_entry['infraclass'], infraclass_entry['parentNameUsageID']]
		# 		for index, value in zip(indexes_to_be_filled, values):
		# 			out_infraclass_row[index] = value
		# 		out_rows.append(out_infraclass_row)
		# 		parent_id = infraclass_entry['taxonID']

		# # create new row for infraclass if it doesn't already exist
		# if magnorder_index is not None and row[magnorder_index] is not None and row[magnorder_index] != 'NA':
		# 	magnorder_in_out_rows = [out_row for out_row in out_rows if out_row[out_tr_index]=='magnorder' and out_row[cn_index]==row[magnorder_index]]
		# 	if magnorder_in_out_rows:
		# 		parent_id = magnorder_in_out_rows[0][id_index]
		# 	else:
		# 		magnorder_entry = create_higher_taxa_entry(row[magnorder_index], 'magnorder', parent_id)
		# 		out_magnorder_row = [None for col in out_header]
		# 		indexes_to_be_filled = [id_index, cn_index, out_tr_index, out_ts_index, magnorder_index, out_pnu_index]
		# 		values = [magnorder_entry['taxonID'], magnorder_entry['canonicalName'], magnorder_entry['taxonomicRank'], magnorder_entry['taxonomicStatus'], \
		# 					magnorder_entry['magnorder'], magnorder_entry['parentNameUsageID']]
		# 		for index, value in zip(indexes_to_be_filled, values):
		# 			out_magnorder_row[index] = value
		# 		out_rows.append(out_magnorder_row)
		# 		parent_id = magnorder_entry['taxonID']

		# # create new row for infraclass if it doesn't already exist
		# if superorder_index is not None and row[superorder_index] is not None and row[superorder_index] != 'NA':
		# 	superorder_in_out_rows = [out_row for out_row in out_rows if out_row[out_tr_index]=='superorder' and out_row[cn_index]==row[superorder_index]]
		# 	if superorder_in_out_rows:
		# 		parent_id = superorder_in_out_rows[0][id_index]
		# 	else:
		# 		superorder_entry = create_higher_taxa_entry(row[superorder_index], 'superorder', parent_id)
		# 		out_superorder_row = [None for col in out_header]
		# 		indexes_to_be_filled = [id_index, cn_index, out_tr_index, out_ts_index, superorder_index, out_pnu_index]
		# 		values = [superorder_entry['taxonID'], superorder_entry['canonicalName'], superorder_entry['taxonomicRank'], superorder_entry['taxonomicStatus'], \
		# 					superorder_entry['superorder'], superorder_entry['parentNameUsageID']]
		# 		for index, value in zip(indexes_to_be_filled, values):
		# 			out_superorder_row[index] = value
		# 		out_rows.append(out_superorder_row)
		# 		parent_id = superorder_entry['taxonID']

		# # create new row for infraclass if it doesn't already exist
		# if order_index is not None and row[order_index] is not None and row[order_index] != 'NA':
		# 	order_in_out_rows = [out_row for out_row in out_rows if out_row[out_tr_index]=='order' and out_row[cn_index]==row[order_index]]
		# 	if order_in_out_rows:
		# 		parent_id = order_in_out_rows[0][id_index]
		# 	else:
		# 		order_entry = create_higher_taxa_entry(row[order_index], 'order', parent_id)
		# 		out_order_row = [None for col in out_header]
		# 		indexes_to_be_filled = [id_index, cn_index, out_tr_index, out_ts_index, order_index, out_pnu_index]
		# 		values = [order_entry['taxonID'], order_entry['canonicalName'], order_entry['taxonomicRank'], order_entry['taxonomicStatus'], \
		# 					order_entry['order'], order_entry['parentNameUsageID']]
		# 		for index, value in zip(indexes_to_be_filled, values):
		# 			out_order_row[index] = value
		# 		out_rows.append(out_order_row)
		# 		parent_id = order_entry['taxonID']

		# # create new row for infraclass if it doesn't already exist
		# if suborder_index is not None and row[suborder_index] is not None and row[suborder_index] != 'NA':
		# 	suborder_in_out_rows = [out_row for out_row in out_rows if out_row[out_tr_index]=='suborder' and out_row[cn_index]==row[suborder_index]]
		# 	if suborder_in_out_rows:
		# 		parent_id = suborder_in_out_rows[0][id_index]
		# 	else:
		# 		suborder_entry = create_higher_taxa_entry(row[suborder_index], 'suborder', parent_id)
		# 		out_suborder_row = [None for col in out_header]
		# 		indexes_to_be_filled = [id_index, cn_index, out_tr_index, out_ts_index, suborder_index, out_pnu_index]
		# 		values = [suborder_entry['taxonID'], suborder_entry['canonicalName'], suborder_entry['taxonomicRank'], suborder_entry['taxonomicStatus'], \
		# 					suborder_entry['suborder'], suborder_entry['parentNameUsageID']]
		# 		for index, value in zip(indexes_to_be_filled, values):
		# 			out_suborder_row[index] = value
		# 		out_rows.append(out_suborder_row)
		# 		parent_id = suborder_entry['taxonID']

		# # create new row for infraclass if it doesn't already exist
		# if infraorder_index is not None and row[infraorder_index] is not None and row[infraorder_index] != 'NA':
		# 	infraorder_in_out_rows = [out_row for out_row in out_rows if out_row[out_tr_index]=='infraorder' and out_row[cn_index]==row[infraorder_index]]
		# 	if infraorder_in_out_rows:
		# 		parent_id = infraorder_in_out_rows[0][id_index]
		# 	else:
		# 		infraorder_entry = create_higher_taxa_entry(row[infraorder_index], 'infraorder', parent_id)
		# 		out_infraorder_row = [None for col in out_header]
		# 		indexes_to_be_filled = [id_index, cn_index, out_tr_index, out_ts_index, infraorder_index, out_pnu_index]
		# 		values = [infraorder_entry['taxonID'], infraorder_entry['canonicalName'], infraorder_entry['taxonomicRank'], infraorder_entry['taxonomicStatus'], \
		# 					infraorder_entry['infraorder'], infraorder_entry['parentNameUsageID']]
		# 		for index, value in zip(indexes_to_be_filled, values):
		# 			out_infraorder_row[index] = value
		# 		out_rows.append(out_infraorder_row)
		# 		parent_id = infraorder_entry['taxonID']

		# # create new row for infraclass if it doesn't already exist
		# if parvorder_index is not None and row[parvorder_index] is not None and row[parvorder_index] != 'NA':
		# 	parvorder_in_out_rows = [out_row for out_row in out_rows if out_row[out_tr_index]=='parvorder' and out_row[cn_index]==row[parvorder_index]]
		# 	if parvorder_in_out_rows:
		# 		parent_id = parvorder_in_out_rows[0][id_index]
		# 	else:
		# 		parvorder_entry = create_higher_taxa_entry(row[parvorder_index], 'parvorder', parent_id)
		# 		out_parvorder_row = [None for col in out_header]
		# 		indexes_to_be_filled = [id_index, cn_index, out_tr_index, out_ts_index, parvorder_index, out_pnu_index]
		# 		values = [parvorder_entry['taxonID'], parvorder_entry['canonicalName'], parvorder_entry['taxonomicRank'], parvorder_entry['taxonomicStatus'], \
		# 					parvorder_entry['parvorder'], parvorder_entry['parentNameUsageID']]
		# 		for index, value in zip(indexes_to_be_filled, values):
		# 			out_parvorder_row[index] = value
		# 		out_rows.append(out_parvorder_row)
		# 		parent_id = parvorder_entry['taxonID']

		# # create new row for infraclass if it doesn't already exist
		# if superfamily_index is not None and row[superfamily_index] is not None and row[superfamily_index] != 'NA':
		# 	superfamily_in_out_rows = [out_row for out_row in out_rows if out_row[out_tr_index]=='superfamily' and out_row[cn_index]==row[superfamily_index]]
		# 	if superfamily_in_out_rows:
		# 		parent_id = superfamily_in_out_rows[0][id_index]
		# 	else:
		# 		superfamily_entry = create_higher_taxa_entry(row[superfamily_index], 'superfamily', parent_id)
		# 		out_superfamily_row = [None for col in out_header]
		# 		indexes_to_be_filled = [id_index, cn_index, out_tr_index, out_ts_index, superfamily_index, out_pnu_index]
		# 		values = [superfamily_entry['taxonID'], superfamily_entry['canonicalName'], superfamily_entry['taxonomicRank'], superfamily_entry['taxonomicStatus'], \
		# 					superfamily_entry['superfamily'], superfamily_entry['parentNameUsageID']]
		# 		for index, value in zip(indexes_to_be_filled, values):
		# 			out_superfamily_row[index] = value
		# 		out_rows.append(out_superfamily_row)
		# 		parent_id = superfamily_entry['taxonID']

		# # create new row for infraclass if it doesn't already exist
		# if family_index is not None and row[family_index] is not None and row[family_index] != 'NA':
		# 	family_in_out_rows = [out_row for out_row in out_rows if out_row[out_tr_index]=='family' and out_row[cn_index]==row[family_index]]
		# 	if family_in_out_rows:
		# 		parent_id = family_in_out_rows[0][id_index]
		# 	else:
		# 		family_entry = create_higher_taxa_entry(row[family_index], 'family', parent_id)
		# 		out_family_row = [None for col in out_header]
		# 		indexes_to_be_filled = [id_index, cn_index, out_tr_index, out_ts_index, family_index, out_pnu_index]
		# 		values = [family_entry['taxonID'], family_entry['canonicalName'], family_entry['taxonomicRank'], family_entry['taxonomicStatus'], \
		# 					family_entry['family'], family_entry['parentNameUsageID']]
		# 		for index, value in zip(indexes_to_be_filled, values):
		# 			out_family_row[index] = value
		# 		out_rows.append(out_family_row)
		# 		parent_id = family_entry['taxonID']

		# # create new row for infraclass if it doesn't already exist
		# if subfamily_index is not None and row[subfamily_index] is not None and row[subfamily_index] != 'NA':
		# 	subfamily_in_out_rows = [out_row for out_row in out_rows if out_row[out_tr_index]=='subfamily' and out_row[cn_index]==row[subfamily_index]]
		# 	if subfamily_in_out_rows:
		# 		parent_id = subfamily_in_out_rows[0][id_index]
		# 	else:
		# 		subfamily_entry = create_higher_taxa_entry(row[subfamily_index], 'subfamily', parent_id)
		# 		out_subfamily_row = [None for col in out_header]
		# 		indexes_to_be_filled = [id_index, cn_index, out_tr_index, out_ts_index, subfamily_index, out_pnu_index]
		# 		values = [subfamily_entry['taxonID'], subfamily_entry['canonicalName'], subfamily_entry['taxonomicRank'], subfamily_entry['taxonomicStatus'], \
		# 					subfamily_entry['subfamily'], subfamily_entry['parentNameUsageID']]
		# 		for index, value in zip(indexes_to_be_filled, values):
		# 			out_subfamily_row[index] = value
		# 		out_rows.append(out_subfamily_row)
		# 		parent_id = subfamily_entry['taxonID']

		# # create new row for infraclass if it doesn't already exist
		# if tribe_index is not None and row[tribe_index] is not None and row[tribe_index] != 'NA':
		# 	tribe_in_out_rows = [out_row for out_row in out_rows if out_row[out_tr_index]=='tribe' and out_row[cn_index]==row[tribe_index]]
		# 	if tribe_in_out_rows:
		# 		parent_id = tribe_in_out_rows[0][id_index]
		# 	else:
		# 		tribe_entry = create_higher_taxa_entry(row[tribe_index], 'tribe', parent_id)
		# 		out_tribe_row = [None for col in out_header]
		# 		indexes_to_be_filled = [id_index, cn_index, out_tr_index, out_ts_index, tribe_index, out_pnu_index]
		# 		values = [tribe_entry['taxonID'], tribe_entry['canonicalName'], tribe_entry['taxonomicRank'], tribe_entry['taxonomicStatus'], \
		# 					tribe_entry['tribe'], tribe_entry['parentNameUsageID']]
		# 		for index, value in zip(indexes_to_be_filled, values):
		# 			out_tribe_row[index] = value
		# 		out_rows.append(out_tribe_row)
		# 		parent_id = tribe_entry['taxonID']

		# # create new row for infraclass if it doesn't already exist
		# if genus_index is not None and row[genus_index] is not None and row[genus_index] != 'NA':
		# 	genus_in_out_rows = [out_row for out_row in out_rows if out_row[out_tr_index]=='genus' and out_row[cn_index]==row[genus_index]]
		# 	if genus_in_out_rows:
		# 		parent_id = genus_in_out_rows[0][id_index]
		# 	else:
		# 		genus_entry = create_higher_taxa_entry(row[genus_index], 'genus', parent_id)
		# 		out_genus_row = [None for col in out_header]
		# 		indexes_to_be_filled = [id_index, cn_index, out_tr_index, out_ts_index, genus_index, out_pnu_index]
		# 		values = [genus_entry['taxonID'], genus_entry['canonicalName'], genus_entry['taxonomicRank'], genus_entry['taxonomicStatus'], \
		# 					genus_entry['genus'], genus_entry['parentNameUsageID']]
		# 		for index, value in zip(indexes_to_be_filled, values):
		# 			out_genus_row[index] = value
		# 		out_rows.append(out_genus_row)
		# 		parent_id = genus_entry['taxonID']

		# # create new row for infraclass if it doesn't already exist
		# if subgenus_index is not None and row[subgenus_index] is not None and row[subgenus_index] != 'NA':
		# 	subgenus_in_out_rows = [out_row for out_row in out_rows if out_row[out_tr_index]=='subgenus' and out_row[cn_index]==row[subgenus_index]]
		# 	if subgenus_in_out_rows:
		# 		parent_id = subgenus_in_out_rows[0][id_index]
		# 	else:
		# 		subgenus_entry = create_higher_taxa_entry(row[subgenus_index], 'subgenus', parent_id)
		# 		out_subgenus_row = [None for col in out_header]
		# 		indexes_to_be_filled = [id_index, cn_index, out_tr_index, out_ts_index, subgenus_index, out_pnu_index]
		# 		values = [subgenus_entry['taxonID'], subgenus_entry['canonicalName'], subgenus_entry['taxonomicRank'], subgenus_entry['taxonomicStatus'], \
		# 					subgenus_entry['subgenus'], subgenus_entry['parentNameUsageID']]
		# 		for index, value in zip(indexes_to_be_filled, values):
		# 			out_subgenus_row[index] = value
		# 		out_rows.append(out_subgenus_row)
		# 		parent_id = subgenus_entry['taxonID']

		out_row[out_pnu_index] = parent_id


		# Create new row for genus and subgenus, and add parentNameUsageID to the outgoing row
		# exist_genus = is_genus_row_available(row[genus_index])
		# genus_row_indexes_to_be_filled = [id_index, cn_index, out_tr_index, out_ts_index, genus_index, out_gn_index]
		# if subgenus_index is None or row[subgenus_index] is None or row[subgenus_index]=='NA':
		# 	# check if genus already exists in any of the MDD versions, and if so returns a dict, else False
		# 	out_genus_row = [None for col in out_header]
		# 	if not exist_genus:
		# 		taxon_id, taxonID_entries = generate_taxonID(taxonID_entries)
		# 		out_row[out_pnu_index] = taxon_id
		# 		values = [taxon_id, row[genus_index], 'genus', 'accepted', row[genus_index], row[genus_index]]
		# 		for index, value in zip(genus_row_indexes_to_be_filled, values):
		# 			out_genus_row[index] = value
		# 		out_rows.append(out_genus_row)
		# 		genus_entries.append({'taxonID': taxon_id, 'managedID': taxon_id, 'canonicalName': row[genus_index], 'taxonRank': 'genus', 'taxonomicStatus': 'accepted', \
		# 										'genus': row[genus_index], 'genericName': row[genus_index]})
		# 	else:
		# 		genus_row_added_to_out_rows = [out_row for out_row in out_rows if out_row[out_tr_index]=='genus' and out_row[cn_index]==row[genus_index]]
		# 		if not genus_row_added_to_out_rows:
		# 			values = [exist_genus['taxonID'], row[genus_index], 'genus', 'accepted', row[genus_index], row[genus_index]]
		# 			for index, value in zip(genus_row_indexes_to_be_filled, values):
		# 				out_genus_row[index] = value
		# 			out_rows.append(out_genus_row)
		# 		out_row[out_pnu_index] = exist_genus['taxonID']
		# else:
		# 	exist_subgenus = is_subgenus_row_available(row[subgenus_index])
		# 	out_genus_row = [None for col in out_header]
		# 	out_subgenus_row = [None for col in out_header]
		# 	out_subgenus_row_indexes_to_be_filled = [id_index, cn_index, out_tr_index, out_ts_index, genus_index, subgenus_index, out_gn_index, out_pnu_index]
		# 	if exist_subgenus:
		# 		assert exist_subgenus and exist_subgenus
		# 		genus_row_added_to_out_rows = [out_row for out_row in out_rows if out_row[out_tr_index]=='genus' and out_row[cn_index]==row[genus_index]]
		# 		subgenus_row_added_to_out_rows = [out_row for out_row in out_rows if out_row[out_tr_index]=='subgenus' and out_row[cn_index]==row[subgenus_index]]
		# 		assert len(genus_row_added_to_out_rows)<2
		# 		assert len(subgenus_row_added_to_out_rows)<2

		# 		if not genus_row_added_to_out_rows and not subgenus_row_added_to_out_rows:
		# 			genus_values = [exist_genus['taxonID'], row[genus_index], 'genus', 'accepted', row[genus_index], row[genus_index]]
		# 			for index, value in zip(genus_row_indexes_to_be_filled, genus_values):
		# 				out_genus_row[index] = value
		# 			subgenus_values = [exist_subgenus['taxonID'], row[subgenus_index], 'subgenus', 'accepted', row[genus_index], row[subgenus_index], \
		# 								row[genus_index], exist_subgenus['parentNameUsageID']]
		# 			for index, value in zip(genus_row_indexes_to_be_filled, subgenus_values):
		# 				out_subgenus_row[index] = value
		# 			out_rows.extend([out_genus_row, out_subgenus_row])
		# 		elif genus_row_added_to_out_rows and not subgenus_row_added_to_out_rows:
		# 			subgenus_values = [exist_subgenus['taxonID'], row[subgenus_index], 'subgenus', 'accepted', row[genus_index], row[subgenus_index], \
		# 								row[genus_index], exist_subgenus['parentNameUsageID']]
		# 			for index, value in zip(genus_row_indexes_to_be_filled, subgenus_values):
		# 				out_subgenus_row[index] = value
		# 			out_rows.append(out_subgenus_row)
		# 		out_row[out_pnu_index] = exist_subgenus['taxonID']

		# 	elif exist_genus and not exist_subgenus:
		# 		genus_row_added_to_out_rows = [out_row for out_row in out_rows if out_row[out_tr_index]=='genus' and out_row[cn_index]==row[genus_index]]
		# 		assert len(genus_row_added_to_out_rows)<2
		# 		if not genus_row_added_to_out_rows:
		# 			genus_values = [exist_genus['taxonID'], row[genus_index], 'genus', 'accepted', row[genus_index], row[genus_index]]
		# 			for index, value in zip(genus_row_indexes_to_be_filled, genus_values):
		# 				out_genus_row[index] = value
		# 			out_rows.append(out_genus_row)
				
		# 		taxon_id, taxonID_entries = generate_taxonID(taxonID_entries)
		# 		out_row[out_pnu_index] = taxon_id
		# 		subgenus_values = [taxon_id, row[subgenus_index], 'subgenus', 'accepted', row[genus_index], row[subgenus_index], row[genus_index], exist_genus['taxonID']]
		# 		for index, value in zip(out_subgenus_row_indexes_to_be_filled, subgenus_values):
		# 			out_subgenus_row[index] = value
		# 		out_rows.append(out_subgenus_row)
		# 		subgenus_entries.append({'taxonID': taxon_id, 'canonicalName': row[subgenus_index], 'taxonRank': 'subgenus', 'taxonomicStatus': 'accepted', \
		# 										'genus': exist_genus['canonicalName'], 'subgenus': row[subgenus_index], 'genericName': row[genus_index], \
		# 										'parentNameUsageID': exist_genus['taxonID']})

		# 	else:
		# 		genus_taxon_id, taxonID_entries = generate_taxonID(taxonID_entries)
		# 		subgenus_taxon_id, taxonID_entries = generate_taxonID(taxonID_entries)
		# 		genus_values = [genus_taxon_id, row[genus_index], 'genus', 'accepted', row[genus_index], row[genus_index]]
		# 		for index, value in zip(genus_row_indexes_to_be_filled, genus_values):
		# 			out_genus_row[index] = value
		# 		subgenus_values = [subgenus_taxon_id, row[subgenus_index], 'subgenus', 'accepted', row[genus_index], row[subgenus_index], row[genus_index], genus_taxon_id]
		# 		for index, value in zip(out_subgenus_row_indexes_to_be_filled, subgenus_values):
		# 			out_subgenus_row[index] = value
		# 		out_row[out_pnu_index] = subgenus_taxon_id
		# 		out_rows.append(out_genus_row)
		# 		out_rows.append(out_subgenus_row)
		# 		genus_entries.append({'taxonID': genus_taxon_id, 'canonicalName': row[genus_index], 'taxonRank': 'genus', 'taxonomicStatus': 'accepted', \
		# 										'genus': row[genus_index], 'genericName': row[genus_index]})
		# 		subgenus_entries.append({'taxonID': subgenus_taxon_id, 'canonicalName': row[subgenus_index], 'taxonRank': 'subgenus', 'taxonomicStatus': 'accepted', \
		# 										'genus': row[genus_index], 'subgenus': row[subgenus_index], 'genericName': row[genus_index], \
		# 										'parentNameUsageID': genus_taxon_id})
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

				# find if a synonym is repeated in the file, perhaps linked to more than one species
				# if [row for row in out_rows if row[out_ts_index] in ['junior synonym', 'senior synonym', 'synonym'] and row[out_sna_index]==syn_authorship and \
				# 		row[ep_index]==syn_epithet]:
				# 	# print(f'{syn_epithet} {syn_authorship}  -  {out_row[cn_index]}')
				# 	continue

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

	gobal_entries_to_search = taxonomic_rank_based_global_entries[taxonomic_rank]
	if len(gobal_entries_to_search) == 0: return False

	for entry in gobal_entries_to_search:
		if name == entry['canonicalName']:
			return entry

	return False

def create_higher_taxa_entry(canonicalName, taxonomicRank, parentNameUsageID):

	global genus_entries, tribe_entries, subfamily_entries, family_entries, superfamily_entries, parvorder_entries, subgenus_entries, \
	infraorder_entries, suborder_entries, order_entries, superorder_entries, magnorder_entries, infraclass_entries, subclass_entries, taxonID_entries

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
	'subgenus': subgenus_entries,
	}

	exist_taxon =  is_taxon_available(canonicalName, taxonomicRank)
	if exist_taxon:
		taxon_id = exist_taxon['taxonID']
	else:
		taxon_id, taxonID_entries = generate_taxonID(taxonID_entries)

	taxon_entry = {'taxonID': taxon_id, 'canonicalName': canonicalName, 'taxonomicRank': taxonomicRank, 'taxonomicStatus': 'accepted',\
							taxonomicRank: canonicalName, 'parentNameUsageID': parentNameUsageID}
	taxonomic_rank_based_global_entries[taxonomicRank].append(taxon_entry)

	return taxon_entry

# def is_subclass_entry_available(subclass_name):
# 	if len(subgenus_entries)==0: return False

# 	for entry in subclass_entries:
# 		if subclass_name==entry['canonicalName']:
# 			return entry

# 	return False 


# def is_genus_row_available(genus_name):
# 	if len(genus_entries)==0: return False

# 	for genus_row in genus_entries:
# 		if genus_name == genus_row['canonicalName']:
# 			return genus_row

# 	return False

# def is_subgenus_row_available(subgenus_name):
# 	if len(subgenus_entries)==0: return False

# 	for subgenus_row in subgenus_entries:
# 		if subgenus_name == subgenus_row['canonicalName']:
# 			return subgenus_row

# 	return False

def is_synonym_row_available(synonym_scientific_name, anu_id):
	if len(synonym_entries)==0: return False

	for synonym_row in synonym_entries:
		if synonym_scientific_name==synonym_row['scientificName'] and anu_id==synonym_row['acceptedNameUsageID']:
			return synonym_row

	return False

def managed_entries_across_MDD():
	global synonym_entries, subgenus_entries, genus_entries, tribe_entries, subfamily_entries, family_entries, superfamily_entries, parvorder_entries, \
	infraorder_entries, suborder_entries, order_entries, superorder_entries, magnorder_entries, infraclass_entries, subclass_entries

	for file_ in MDD_vers:
		with open(file_) as f:
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
					elif row['taxonRank']=='infravorder':
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


def taxonIDs_across_MDD():
	taxonIDs = set()

	for file_ in MDD_vers:
		with open(file_) as f:
			reader = csv.DictReader(f)
			for row in reader:
				taxonIDs.add([row[key] for key in row if key in ['id', 'taxonID', 'ID_number']][0])

	return list(taxonIDs)

if __name__=='__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--input', default='MDD_v1_6495species_JMamm.csv', \
						help='name of the input file to be mapped to DwC format. The file must be in CSV format')
	parser.add_argument('--output', help='name of the output file')
	args=parser.parse_args()

	outfile = args.output
	if outfile is None: outfile =f'{args.input.rsplit(".", 1)[0]}_inDwC.csv'

	managed_entries_across_MDD()
	taxonID_entries = taxonIDs_across_MDD()

	with open(args.input) as f1, open(outfile, 'w') as f2:
		map_MDD_to_DwC(f1, f2)

	# for file_ in MDD_vers:
	# 	with open(file_) as f:
	# 		reader = csv.DictReader(f)
	# 		tax_ids = []
	# 		for row in reader:
	# 			tax_ids.append([row[key] for key in row if key in ['id', 'taxonID', 'ID_number']][0])

	# 		if len(list(set(tax_ids))) == len(tax_ids):
	# 			print(f'No duplicate key in {file_}')
	# 		else:
	# 			unique_ids = set()
	# 			for id_ in tax_ids:
	# 				if id_ is not None and id_ not in unique_ids:
	# 					unique_ids.add(id_)
	# 				else:
	# 					print(f'duplicate id in {file_} : {id_}')






