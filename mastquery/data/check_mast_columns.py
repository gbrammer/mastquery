import yaml
import numpy as np
from mastquery import jwst
import previous_defaults

filters = jwst.CALIB_FILTERS
filters += jwst.FULL_SUBARRAY
filters += jwst.FINE_GUIDE

filters += [{"paramName": "expstart", "values": [{"min": 59730.4, "max": 59735.4}]}]

res = jwst.query_all_jwst(filters=filters, fix=False, columns="*")

inst_query = {}

for inst in ["NIS", "NRC", "MIR", "NRS"]:
    print("Instrument: ", inst)
    inst_query[inst] = jwst.query_jwst(
        instrument=inst,
        columns="*",
        filters=filters,
        extra={"format": "json", "pagesize": 10},
    )

########
# FGS queries
print("Instrument: FGS")
fgs = {}

fgs["FGS"] = jwst.query_jwst(
    instrument="FGS",
    columns="*",
    filters=filters,
)

print("Instrument: GS")
fgs["GS"] = jwst.query_jwst(
    instrument="GS",
    columns="*",
    filters=[
        {"paramName": "date_obs_mjd", "values": [{"min": 59730.4, "max": 59731.4}]},
        {"paramName": "program", "values": [1074]}
    ],
)

all_columns = []

for inst in inst_query:
    all_columns += inst_query[inst].colnames

missing_from_all = ~np.isin(previous_defaults.ALL_COLUMNS, all_columns)
print(
    "ALL_COLUMNS not in current result: ",
    np.array(previous_defaults.ALL_COLUMNS)[missing_from_all],
)

# generated columns
pop_columns = ["dt", "prog_pi", "filter-pupil", "inst-mode"]

for c in np.unique(all_columns):
    if c.endswith("_mjd"):
        ck = c.replace("_mjd", "")
        if ck in all_columns:
            pop_columns.append(ck)

# In all instruments
common_cols = np.unique(all_columns)
for inst in inst_query:
    common_cols = common_cols[np.isin(common_cols, inst_query[inst].colnames)]

inst_columns = {inst: inst_query[inst].colnames for inst in inst_query}

for inst in inst_query:
    cols_ = inst_columns[inst]
    keep = ~np.isin(cols_, pop_columns)
    inst_columns[inst] = np.array(cols_)[keep].tolist()

for inst in fgs:
    cols_ = np.array(fgs[inst].colnames)
    keep = ~np.isin(cols_, pop_columns)
    inst_columns[inst] = cols_[keep].tolist()
    
    previous_defaults.DEFAULT_COLUMNS[inst] = inst_columns[inst]

for inst in inst_query:
    inst_cols = inst_query[inst].colnames

    missing_from_inst = ~np.isin(previous_defaults.DEFAULT_COLUMNS[inst], inst_cols)

    print(
        f"DEFAULT_COLUMNS[{inst}] not in current result: ",
        np.array(previous_defaults.DEFAULT_COLUMNS[inst])[missing_from_inst],
    )

for inst in inst_query:
    print(
        f"{inst} {len(inst_query[inst]):>5} {len(inst_query[inst].colnames):>3}"
        + f"  {len(previous_defaults.DEFAULT_COLUMNS[inst]):>3}"
    )

column_collections = {
    "common": {
        "all": common_cols.tolist(),
        "default": np.array(common_cols)[
            np.isin(common_cols, previous_defaults.ALL_COLUMNS)
        ].tolist(),
    }
}

for inst in inst_columns:
    column_collections[inst] = {"all": inst_columns[inst]}
    keep = np.isin(inst_columns[inst], previous_defaults.DEFAULT_COLUMNS[inst])
    column_collections[inst]["default"] = np.array(inst_columns[inst])[keep].tolist()

column_collections["FGS"] = {
    "all": fgs["FGS"].colnames,
    "default": fgs["FGS"].colnames
}

with open("mast_columns.yaml", "w") as fp:
    yaml.dump(column_collections, stream=fp)
