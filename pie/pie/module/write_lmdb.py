import os
import lmdb
import pickle


def write_lmdb(lmdb_list: list, chrom: list, output: str, bed_target_list=None) -> None:
    """Write event-keyed results, independently scoped by chromosome or BED row."""
    if bed_target_list is None:
        scopes = [("chromosome", str(c)) for c in chrom]
    else:
        scopes = [("region", str(c), int(start), int(end))
                  for c, start, end in bed_target_list]
    if len(scopes) != len(lmdb_list):
        raise ValueError("LMDB results and chromosome/region scopes must have equal lengths")
    if len(set(scopes)) != len(scopes):
        raise ValueError("LMDB chromosome/region scopes must be unique")

    schema = {
        "schema_version": 2,
        "event_key_fields": ("pos1", "variant_hash"),
        "position_coordinates": "1-based VCF POS",
        "region_coordinates": "0-based half-open BED",
        "chromosome_key": "chromosome:{chrom}:{pos1:09d}:{variant_hash}",
        "region_key": "region:{chrom}:{start0:09d}:{end0:09d}:{pos1:09d}:{variant_hash}",
        "scopes": scopes,
        "records_per_scope": [len(events) for events in lmdb_list],
    }
    env = lmdb.open(os.path.abspath(output + ".lmdb"), map_size=10**10)
    try:
        main_db = env.open_db(None)
        with env.begin(write=True) as txn:
            # Replace this output atomically, including any keys from older runs.
            txn.drop(main_db, delete=False)
            txn.put(b"__schema__", pickle.dumps(schema))
            for scope, events in zip(scopes, lmdb_list):
                for event_key, features in events.items():
                    pos, variant_hash = event_key
                    pos, variant_hash = int(pos), str(variant_hash)
                    if scope[0] == "chromosome":
                        key = f"chromosome:{scope[1]}:{pos:09d}:{variant_hash}"
                    else:
                        key = (f"region:{scope[1]}:{scope[2]:09d}:{scope[3]:09d}:"
                               f"{pos:09d}:{variant_hash}")
                    value = dict(features)
                    value["event_key"] = (pos, variant_hash)
                    value["scope"] = scope
                    txn.put(key.encode("utf-8"), pickle.dumps(value))
    finally:
        env.close()


