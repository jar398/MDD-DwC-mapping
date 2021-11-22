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
