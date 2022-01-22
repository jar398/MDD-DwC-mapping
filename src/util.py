#!/usr/bin/env python3

import csv, hashlib

taxonID_count = 1006599

def find_index(header, column):
	if column in header:
		return header.index(column)
	else:
		return None

def del_list_indexes(list_, indexs):
	indexs = [i for i in indexs if i]
	if not indexs:
		return list_
	assert max(indexs) < len(list_), 'index out of range'
	return [j for i, j in enumerate(list_) if i not in indexs]

def generate_hash(list_):
  return hashlib.sha1(" ".join(list_).encode('utf-8')).hexdigest()

def generate_taxonID(existing_ids):
	global taxonID_count
	while True:
		if taxonID_count not in existing_ids:
			existing_ids.append(taxonID_count)
			return taxonID_count, existing_ids
		taxonID_count += 1

def has_duplicate_taxonID(file_):
		
		with open(file_) as f:
			reader = csv.DictReader(f)
			tax_ids = []

			for row in reader:
				id_ = [row[key] for key in row if key in ['id', 'taxonID', 'ID_number']]
				if '' in id_: id_.remove('')
				if id_:
				tax_ids.append(int(id_[0]))

			if len(list(set(tax_ids))) == len(tax_ids):
				print(f'No duplicate taxonID in {file_}')
			else:
				unique_ids = set()
				for id_ in tax_ids:
					if id_ is not None and id_ not in unique_ids:
						unique_ids.add(id_)
					else:
						print(f'duplicate id in {file_} : {id_}')

def are_taxonID_managed_between_versions(file_version_1, file_version_2):
 # check if the two entries from the two files have same id, then they must have same details; if two entries from the two files represent the same taxon, 
 # then they must have the same id
 pass





