#!/usr/bin/env python3

import csv, argparse, sys, re
from util import find_index, del_list_indexes, generate_hash, generate_taxonID

MDD_vers= [
		'MDD_v1_6495species_JMamm_inDwC.csv', 
		'MDD_v1.1_6526species_inDwC.csv', 
		'MDD_v1.2_6485species_inDwC.csv', 
		'MDD_v1.3_6513species_inDwC.csv', 
		'MDD_v1.31_6513species_inDwC.csv', 
		'MDD_v1.4_6533species_inDwC.csv', 
		'MDD_v1.5_6554species_inDwC.csv',
		'MDD_v1.6_6557species_inDwC.csv', 
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
		'majorSubtype': 'majorSubType',
		'MajorSubtype': 'majorSubType',
		'MajorType': 'majorType',
		'majorType': 'majorType',
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

subgenus_entries, genus_entries, taxonID_entries = [], [], []

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
	if ep_index == None: ep_index = -1
	genus_index = find_index(DwC_header, 'genus')
	subgenus_index = find_index(DwC_header, 'subgenus')
	nn_index = find_index(DwC_header, 'nominalNames')

#Columns to be capitalized
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
	global subgenus_entries, genus_entries, taxonID_entries

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

		if ep_index == -1:
			row[ep_index] = row[cn_index].split('_')[1]

		out_row = row + [None for i in range(len(MDD_new_columns))]
		out_row[id_index] = row_id
		columns_to_capitalize = [i for i in [order_index, suborder_index, infraorder_index, parvorder_index, superfamily_index, family_index, subfamily_index, tribe_index] if i is not None]

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

		# Create new row for genus and subgenus, and add parentNameUsageID to the outgoing row
		exist_genus = is_genus_row_available(row[genus_index])
		genus_row_indexes_to_be_filled = [id_index, cn_index, out_tr_index, out_ts_index, genus_index, out_gn_index]
		if subgenus_index is None or row[subgenus_index] is None or row[subgenus_index]=='NA':
			# check if genus already exists in any of the MDD versions, and if so returns a dict, else False
			out_genus_row = [None for col in out_header]
			if not exist_genus:
				taxon_id, taxonID_entries = generate_taxonID(taxonID_entries)
				out_row[out_pnu_index] = taxon_id
				values = [taxon_id, row[genus_index], 'genus', 'accepted', row[genus_index], row[genus_index]]
				for index, value in zip(genus_row_indexes_to_be_filled, values):
					out_genus_row[index] = value
				out_rows.append(out_genus_row)
				genus_entries.append({'taxonID': taxon_id, 'canonicalName': row[genus_index], 'taxonRank': 'genus', 'taxonomicStatus': 'accepted', \
												'genus': row[genus_index], 'genericName': row[genus_index]})
			else:
				genus_row_added_to_out_rows = [out_row for out_row in out_rows if out_row[out_tr_index]=='genus' and out_row[cn_index]==row[genus_index]]
				if not genus_row_added_to_out_rows:
					values = [exist_genus['taxonID'], row[genus_index], 'genus', 'accepted', row[genus_index], row[genus_index]]
					for index, value in zip(genus_row_indexes_to_be_filled, values):
						out_genus_row[index] = value
					out_rows.append(out_genus_row)
				out_row[out_pnu_index] = exist_genus['taxonID']
		else:
			exist_subgenus = is_subgenus_row_available(row[subgenus_index])
			out_genus_row = [None for col in out_header]
			out_subgenus_row = [None for col in out_header]
			out_subgenus_row_indexes_to_be_filled = [id_index, cn_index, out_tr_index, out_ts_index, genus_index, subgenus_index, out_gn_index, out_pnu_index]
			if exist_subgenus:
				assert exist_subgenus and exist_subgenus
				genus_row_added_to_out_rows = [out_row for out_row in out_rows if out_row[out_tr_index]=='genus' and out_row[cn_index]==row[genus_index]]
				subgenus_row_added_to_out_rows = [out_row for out_row in out_rows if out_row[out_tr_index]=='subgenus' and out_row[cn_index]==row[subgenus_index]]
				assert len(genus_row_added_to_out_rows)<2
				assert len(subgenus_row_added_to_out_rows)<2

				if not genus_row_added_to_out_rows and not subgenus_row_added_to_out_rows:
					genus_values = [exist_genus['taxonID'], row[genus_index], 'genus', 'accepted', row[genus_index], row[genus_index]]
					for index, value in zip(genus_row_indexes_to_be_filled, genus_values):
						out_genus_row[index] = value
					subgenus_values = [exist_subgenus['taxonID'], row[subgenus_index], 'subgenus', 'accepted', row[genus_index], row[subgenus_index], \
										row[genus_index], exist_subgenus['parentNameUsageID']]
					for index, value in zip(genus_row_indexes_to_be_filled, subgenus_values):
						out_subgenus_row[index] = value
					out_rows.extend([out_genus_row, out_subgenus_row])
				elif genus_row_added_to_out_rows and not subgenus_row_added_to_out_rows:
					subgenus_values = [exist_subgenus['taxonID'], row[subgenus_index], 'subgenus', 'accepted', row[genus_index], row[subgenus_index], \
										row[genus_index], exist_subgenus['parentNameUsageID']]
					for index, value in zip(genus_row_indexes_to_be_filled, subgenus_values):
						out_subgenus_row[index] = value
					out_rows.append(out_subgenus_row)
				out_row[out_pnu_index] = exist_subgenus['taxonID']

			elif exist_genus and not exist_subgenus:
				genus_row_added_to_out_rows = [out_row for out_row in out_rows if out_row[out_tr_index]=='genus' and out_row[cn_index]==row[genus_index]]
				assert len(genus_row_added_to_out_rows)<2
				if not genus_row_added_to_out_rows:
					genus_values = [exist_genus['taxonID'], row[genus_index], 'genus', 'accepted', row[genus_index], row[genus_index]]
					for index, value in zip(genus_row_indexes_to_be_filled, genus_values):
						out_genus_row[index] = value
					out_rows.append(out_genus_row)
				
				taxon_id, taxonID_entries = generate_taxonID(taxonID_entries)
				out_row[out_pnu_index] = taxon_id
				subgenus_values = [taxon_id, row[subgenus_index], 'subgenus', 'accepted', row[genus_index], row[subgenus_index], row[genus_index], exist_genus['taxonID']]
				for index, value in zip(out_subgenus_row_indexes_to_be_filled, subgenus_values):
					out_subgenus_row[index] = value
				out_rows.append(out_subgenus_row)
				subgenus_entries.append({'taxonID': taxon_id, 'canonicalName': row[subgenus_index], 'taxonRank': 'subgenus', 'taxonomicStatus': 'accepted', \
												'genus': exist_genus['canonicalName'], 'subgenus': row[subgenus_index], 'genericName': row[genus_index], \
												'parentNameUsageID': exist_genus['taxonID']})

			else:
				genus_taxon_id, taxonID_entries = generate_taxonID(taxonID_entries)
				subgenus_taxon_id, taxonID_entries = generate_taxonID(taxonID_entries)
				genus_values = [genus_taxon_id, row[genus_index], 'genus', 'accepted', row[genus_index], row[genus_index]]
				for index, value in zip(genus_row_indexes_to_be_filled, genus_values):
					out_genus_row[index] = value
				subgenus_values = [subgenus_taxon_id, row[subgenus_index], 'subgenus', 'accepted', row[genus_index], row[subgenus_index], row[genus_index], genus_taxon_id]
				for index, value in zip(out_subgenus_row_indexes_to_be_filled, subgenus_values):
					out_subgenus_row[index] = value
				out_row[out_pnu_index] = subgenus_taxon_id
				out_rows.append(out_genus_row)
				out_rows.append(out_subgenus_row)
				genus_entries.append({'taxonID': genus_taxon_id, 'canonicalName': row[genus_index], 'taxonRank': 'genus', 'taxonomicStatus': 'accepted', \
												'genus': row[genus_index], 'genericName': row[genus_index]})
				subgenus_entries.append({'taxonID': subgenus_taxon_id, 'canonicalName': row[subgenus_index], 'taxonRank': 'subgenus', 'taxonomicStatus': 'accepted', \
												'genus': row[genus_index], 'subgenus': row[subgenus_index], 'genericName': row[genus_index], \
												'parentNameUsageID': genus_taxon_id})
		out_rows.append(out_row)
		
		if nn_index is not None and row[nn_index] is not None:
			synonyms_list = row[nn_index].split('|')
			for syn in synonyms_list:
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
				

				scientific_name = f'{row[genus_index]} {syn_epithet} {syn_authorship}'
				canonical_name = f'{row[genus_index]} {syn_epithet}'

				syn_row = [None for col in out_header]
				syn_taxon_id, taxonID_entries = generate_taxonID(taxonID_entries)
				syn_indexes = [id_index, out_pnu_index, out_anu_index, out_sn_index, out_sna_index, cn_index, out_gn_index, ep_index, out_tr_index, out_ts_index, out_ns_index, genus_index]
				syn_values = [syn_taxon_id, out_row[out_pnu_index], row_id, scientific_name, syn_authorship, canonical_name, row[genus_index], syn_epithet, \
								'species', 'synonym', nomenclatural_status, row[genus_index]]
				for index, value in zip(syn_indexes, syn_values):
					syn_row[index] = value
				out_rows.append(syn_row)

	writer = csv.writer(outfile)
	indexes_to_delete = {author_index, year_index, paranthesis_index, nn_index}
	out_header = del_list_indexes(out_header, indexes_to_delete)
	writer.writerow(out_header)
	for row in out_rows:
		trimmed_row = del_list_indexes(row, indexes_to_delete)
		writer.writerow(trimmed_row)

def is_genus_row_available(genus_name):
	if len(genus_entries)==0: return False

	for genus_row in genus_entries:
		if genus_name == genus_row['canonicalName']:
			return genus_row

	return False

def is_subgenus_row_available(subgenus_name):
	if len(subgenus_entries)==0: return False

	for subgenus_row in subgenus_entries:
		if subgenus_name == subgenus_row['canonicalName']:
			return subgenus_row

	return False

def genus_subgenus_entries_across_MDD():
	genus_entries = []	
	subgenus_entries = []
	for file_ in MDD_vers:
		with open(file_) as f:
			reader = csv.DictReader(f)
			for row in reader:
				if 'taxonRank' not in row:
					break
				else:
					if row['taxonRank']=='genus' and any([True for dict_ in genus_entries if row['canonicalName'] == dict_['canonicalName']]):
						genus_entries.append(row)
					if row['taxonRank']=='subgenus' and any([True for dict_ in subgenus_entries if row['canonicalName'] == dict_['canonicalName']]):
						subgenus_entries.append(row)

	return genus_entries, subgenus_entries

def taxonIDs_across_MDD():
	taxonIDs = set()

	for file_ in MDD_vers:
		with open(file_) as f:
			reader = csv.DictReader(f)
			for row in reader:
				taxonIDs.add([row[key] for key in row if key in ['id', 'taxonID', 'ID_number']][0])

	return list(taxonIDs)

if __name__=='__main__':
	# parser = argparse.ArgumentParser()
	# parser.add_argument('--input', default='MDD_v1_6495species_JMamm.csv', \
	# 					help='name of the input file to be mapped to DwC format. The file must be in CSV format')
	# parser.add_argument('--output', help='name of the output file')
	# args=parser.parse_args()

	# outfile = args.output
	# if outfile is None: outfile =f'{args.input.rsplit(".", 1)[0]}_inDwC.csv'

	# genus_entries, subgenus_entries = genus_subgenus_entries_across_MDD()
	# taxonID_entries = taxonIDs_across_MDD()

	# with open(args.input) as f1, open(outfile, 'w') as f2:
	# 	map_MDD_to_DwC(f1, f2)

	for file_ in MDD_vers:
		with open(file_) as f:
			reader = csv.DictReader(f)
			tax_ids = []
			for row in reader:
				tax_ids.append([row[key] for key in row if key in ['id', 'taxonID', 'ID_number']][0])

			if len(list(set(tax_ids))) == len(tax_ids):
				print(f'No duplicate key in {file_}')
			else:
				unique_ids = set()
				for id_ in tax_ids:
					if id_ is not None and id_ not in unique_ids:
						unique_ids.add(id_)
					else:
						print(f'duplicate id in {file_} : {id_}')






