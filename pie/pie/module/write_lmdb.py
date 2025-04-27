import os
import lmdb
import pickle


def write_lmdb(lmdb_list: list, chrom: list, output: str) -> None:
    lmdb_dict = dict()
    for i, j in zip(chrom, lmdb_list):
        lmdb_dict[i] = j    

    env = lmdb.open(os.path.abspath(output + ".lmdb"), map_size = 1e10)

    with env.begin(write=True) as txn:
        for c, pos_dict in lmdb_dict.items():
            for pos, features in pos_dict.items():
                # fill 0 to nine digits 00000100
                key = f"{c}:{pos:09d}".encode()
                value = pickle.dumps(features)
                txn.put(key, value)
        
