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

# Allocate a fresh taxonID.  Modified by JAR
def generate_taxonID(all_taxonIDs):
  global taxonID_count
  (max_id, id_set) = all_taxonIDs
  new_id = max_id + 1
  # id_set.add(new_id) - not needed 
  all_taxonIDs[0] = new_id      # Update max
  return (new_id, all_taxonIDs)

