from pathlib import Path
import os
import csv
import gzip
import json
import math
from collections import defaultdict

import numpy as np
from scipy.stats import chi2, norm

ROOT = Path(os.environ["P5_DATA_PARENT"]) / "paper5_evolution_spatial_bridge"
P3 = Path(os.environ["P5_DATA_PARENT"]) / "paper3_GSE300977_raw"

SCORE = Path(
    "P5_STAGE1A_FULL_TMD_EVOLUTION_SCORE.tsv"
)

MAP = (
    P3 /
    "pseudo_transfer_preregister/bridge/"
    "GSE300977_REAL_PSEUDO_PRIMARY_PLUS75_UNITS_v1.tsv"
)

DENS = (
    P3 /
    "tmd_spatial_primary_reveal_v1/output/"
    "PRIMARY_PLUS75_UNIT_DENSITIES_v1.tsv"
)

TMD = (
    ROOT /
    "03_TMD_MAPPING_FEASIBILITY/"
    "P3_TO_P5_PRIMARY_PLUS75_HUMAN_TMD_COORDINATES_FROZEN.tsv"
)

ALIGN_CACHE = Path(
    "ENSEMBL_FULL_ALIGNMENT_CACHE"
)

for p in (
    SCORE,
    MAP,
    DENS,
    TMD,
    ALIGN_CACHE,
):
    if not p.exists():
        raise FileNotFoundError(p)

# ============================================================
# THIS IS THE FIRST FORMAL P3 PHENOTYPE REVEAL FOR P5.
#
# Frozen before execution:
# - evolutionary predictor
# - POST30..90 primary window
# - 4505-TMD join universe
# - P3 D formula
# - location/library mapping
# - primary/adjusted/within/pseudo models
# - STOP rule
# ============================================================

LOCATIONS = [
    "ER",
    "mitochondria",
    "cytosol",
]

PAIRS = {
    "ER": [
        ("B01", "B02"),
        ("B11", "B12"),
    ],

    "mitochondria": [
        ("B09", "B10"),
        ("B07", "B08"),
    ],

    "cytosol": [
        ("B05", "B06"),
        ("B03", "B04"),
    ],
}

KD = {
    "I": 4.5,
    "V": 4.2,
    "L": 3.8,
    "F": 2.8,
    "C": 2.5,
    "M": 1.9,
    "A": 1.8,
    "G": -0.4,
    "T": -0.7,
    "S": -0.8,
    "W": -0.9,
    "Y": -1.3,
    "P": -1.6,
    "H": -3.2,
    "E": -3.5,
    "Q": -3.5,
    "D": -3.5,
    "N": -3.5,
    "K": -3.9,
    "R": -4.5,
}


def split_values(x):
    x = str(x).strip()

    if not x:
        return []

    return [
        y.strip()
        for y in x.split(",")
        if y.strip()
    ]


def fnum(x):
    try:
        v = float(x)

        if math.isfinite(v):
            return v

    except Exception:
        pass

    return None


def ungap(x):
    return (
        str(x)
        .replace("-", "")
        .replace(".", "")
        .replace(" ", "")
        .replace("\n", "")
        .replace("\r", "")
        .upper()
    )


def load_gz(path):
    with gzip.open(
        path,
        "rt",
        encoding="utf-8"
    ) as fh:
        return json.load(fh)


def z_reference(values):
    arr = np.asarray(
        [
            x for x in values
            if x is not None
            and math.isfinite(x)
        ],
        dtype=float
    )

    if len(arr) < 2:
        raise RuntimeError(
            "cannot standardize fewer than 2 values"
        )

    mu = float(
        np.mean(arr)
    )

    sd = float(
        np.std(
            arr,
            ddof=1
        )
    )

    if not math.isfinite(sd) or sd <= 0:
        raise RuntimeError(
            "standardization SD <= 0"
        )

    return mu, sd


def z_apply(x, mu, sd):
    if x is None:
        return None

    return (
        float(x) - mu
    ) / sd


def holm_adjust(pvals):
    n = len(pvals)

    order = sorted(
        range(n),
        key=lambda i:
            pvals[i]
    )

    adj = [None] * n
    running = 0.0

    for rank, i in enumerate(
        order,
        start=1
    ):

        value = (
            n - rank + 1
        ) * pvals[i]

        running = max(
            running,
            value
        )

        adj[i] = min(
            running,
            1.0
        )

    return adj


# ============================================================
# Cluster-robust OLS
# CR1 finite-sample correction.
# ============================================================

def fit_cluster_ols(
    X,
    y,
    clusters,
    names,
):

    X = np.asarray(
        X,
        dtype=float
    )

    y = np.asarray(
        y,
        dtype=float
    )

    clusters = np.asarray(
        clusters,
        dtype=object
    )

    n,k = X.shape

    if n <= k:
        raise RuntimeError(
            f"OLS n<=k: {n}<={k}"
        )

    xtx_inv = np.linalg.pinv(
        X.T @ X
    )

    beta = (
        xtx_inv
        @ X.T
        @ y
    )

    resid = (
        y
        -
        X @ beta
    )

    unique_clusters = np.unique(
        clusters
    )

    G = len(
        unique_clusters
    )

    if G <= 1:
        raise RuntimeError(
            "cluster count <=1"
        )

    meat = np.zeros(
        (k,k),
        dtype=float
    )

    for g in unique_clusters:

        idx = (
            clusters == g
        )

        Xg = X[idx]
        ug = resid[idx]

        score = (
            Xg.T
            @ ug
        )

        meat += np.outer(
            score,
            score
        )

    correction = (
        G / (G - 1)
    ) * (
        (n - 1)
        /
        (n - k)
    )

    cov = (
        correction
        *
        xtx_inv
        @ meat
        @ xtx_inv
    )

    se = np.sqrt(
        np.maximum(
            np.diag(cov),
            0
        )
    )

    return {
        "beta":
            beta,

        "cov":
            cov,

        "se":
            se,

        "names":
            list(names),

        "n":
            n,

        "k":
            k,

        "clusters":
            G,
    }


def wald_test(
    fit,
    R,
    q=None,
):

    beta = fit[
        "beta"
    ]

    cov = fit[
        "cov"
    ]

    R = np.asarray(
        R,
        dtype=float
    )

    if R.ndim == 1:
        R = R.reshape(
            1,
            -1
        )

    if q is None:
        q = np.zeros(
            R.shape[0]
        )

    else:
        q = np.asarray(
            q,
            dtype=float
        )

    diff = (
        R @ beta
        - q
    )

    V = (
        R
        @ cov
        @ R.T
    )

    W = float(
        diff.T
        @ np.linalg.pinv(V)
        @ diff
    )

    df = int(
        np.linalg.matrix_rank(
            R
        )
    )

    p = float(
        chi2.sf(
            W,
            df
        )
    )

    return W, df, p


def coefficient_stats(
    fit,
    indices,
):

    rows = []

    for idx in indices:

        beta = float(
            fit[
                "beta"
            ][idx]
        )

        se = float(
            fit[
                "se"
            ][idx]
        )

        if se > 0:

            z = beta / se

            p = float(
                2
                *
                norm.sf(
                    abs(z)
                )
            )

        else:

            z = float("nan")
            p = float("nan")

        rows.append({
            "parameter":
                fit[
                    "names"
                ][idx],

            "beta":
                beta,

            "SE":
                se,

            "CI_2.5":
                beta
                - 1.959963984540054
                * se,

            "CI_97.5":
                beta
                + 1.959963984540054
                * se,

            "z":
                z,

            "p":
                p,
        })

    valid_ps = [
        x["p"]
        for x in rows
    ]

    if all(
        math.isfinite(x)
        for x in valid_ps
    ):

        adj = holm_adjust(
            valid_ps
        )

        for row,padj in zip(
            rows,
            adj
        ):
            row[
                "Holm_p"
            ] = padj

    else:

        for row in rows:
            row[
                "Holm_p"
            ] = float("nan")

    return rows


# ============================================================
# 1. Frozen P5 predictor table
# ============================================================

score_rows = {}

with open(
    SCORE,
    encoding="utf-8"
) as fh:

    r = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in r:

        if (
            row[
                "is_internal"
            ] != "1"
            or
            row[
                "primary_score_available"
            ] != "1"
        ):

            continue

        key = (
            row[
                "human_ensembl_gene"
            ],
            row[
                "human_ensembl_transcript"
            ],
            int(
                row[
                    "tmd_index"
                ]
            ),
        )

        if key in score_rows:
            raise RuntimeError(
                f"duplicate P5 score key {key}"
            )

        score_rows[
            key
        ] = row


if len(
    score_rows
) != 6511:

    raise RuntimeError(
        f"expected 6511 scored internal TMDs, "
        f"got {len(score_rows)}"
    )


# ============================================================
# 2. Frozen P3 paired-unit metadata
# ============================================================

key_to_units = defaultdict(
    set
)

unit_meta = {}

with open(
    MAP,
    encoding="utf-8"
) as fh:

    r = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in r:

        txs = split_values(
            row[
                "target_transcripts"
            ]
        )

        genes = split_values(
            row[
                "gene_ids"
            ]
        )

        tids = split_values(
            row[
                "tmd_indices"
            ]
        )

        unit = row[
            "paired_real_unit_id"
        ]

        unit_meta[
            unit
        ] = row

        if (
            len(txs) == 1
            and
            len(genes) == 1
            and
            len(tids) == 1
        ):

            key = (
                genes[0],
                txs[0],
                int(
                    tids[0]
                )
            )

            key_to_units[
                key
            ].add(
                unit
            )


# ============================================================
# 3. FIRST REVEAL:
#    Read unit-level P3 density values.
# ============================================================

density = {}

with open(
    DENS,
    encoding="utf-8"
) as fh:

    r = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in r:

        unit = row[
            "paired_real_unit_id"
        ]

        if unit in density:
            raise RuntimeError(
                f"duplicate density unit: {unit}"
            )

        density[
            unit
        ] = row


# ============================================================
# 4. Build exactly prelocked 4505 bridge TMDs
# ============================================================

bridge = []

for key,row in sorted(
    score_rows.items()
):

    units = key_to_units.get(
        key,
        set()
    )

    if len(units) != 1:
        continue

    unit = next(
        iter(
            units
        )
    )

    if unit not in density:
        continue

    bridge.append({
        "key":
            key,

        "gene":
            key[0],

        "transcript":
            key[1],

        "tmd_index":
            key[2],

        "paired_real_unit_id":
            unit,

        "score":
            fnum(
                row[
                    "POST30_90_4D_EXCESS_CONSERVATION_Z_MEDIAN"
                ]
            ),

        "aa_identity":
            fnum(
                row[
                    "median_window_aa_identity"
                ]
            ),

        "GC3":
            fnum(
                row[
                    "human_window_GC3"
                ]
            ),

        "CpG":
            fnum(
                row[
                    "human_window_CpG_density"
                ]
            ),

        "density":
            density[
                unit
            ],
    })


if len(
    bridge
) != 4505:

    raise RuntimeError(
        f"expected frozen bridge universe 4505, "
        f"got {len(bridge)}"
    )


if any(
    x[
        "score"
    ] is None
    for x in bridge
):

    raise RuntimeError(
        "missing frozen primary score in bridge universe"
    )


# ============================================================
# 5. TMD length + Kyte-Doolittle hydropathy
#    Outcome-independent frozen covariates.
# ============================================================

coord = {}

with open(
    TMD,
    encoding="utf-8"
) as fh:

    r = csv.DictReader(
        fh,
        delimiter="\t"
    )

    bridge_keys = {
        x[
            "key"
        ]
        for x in bridge
    }

    for row in r:

        key = (
            row[
                "human_ensembl_gene"
            ],
            row[
                "human_ensembl_transcript"
            ],
            int(
                row[
                    "tmd_index"
                ]
            ),
        )

        if key not in bridge_keys:
            continue

        coord[
            key
        ] = {
            "start0":
                int(
                    row[
                        "human_tmd_aa_start0"
                    ]
                ),

            "end0":
                int(
                    row[
                        "human_tmd_aa_end0_exclusive"
                    ]
                ),

            "length":
                int(
                    row[
                        "human_tmd_length_aa"
                    ]
                ),
        }


if len(
    coord
) != 4505:

    raise RuntimeError(
        f"TMD coordinate coverage {len(coord)}/4505"
    )


protein_by_gene = {}

for gene in sorted({
    x[
        "gene"
    ]
    for x in bridge
}):

    obj = load_gz(
        ALIGN_CACHE /
        f"{gene}.protein_alignment.json.gz"
    )

    entries = obj.get(
        "data",
        []
    )

    if len(entries) != 1:
        raise RuntimeError(
            f"bad alignment schema {gene}"
        )

    proteins = set()

    for h in entries[
        0
    ].get(
        "homologies",
        []
    ):

        if h.get(
            "type"
        ) != "ortholog_one2one":
            continue

        src = h.get(
            "source",
            {}
        )

        aln = src.get(
            "align_seq",
            ""
        )

        if aln:
            proteins.add(
                ungap(
                    aln
                )
            )

    if len(
        proteins
    ) != 1:

        raise RuntimeError(
            f"non-unique frozen human protein {gene}"
        )

    protein_by_gene[
        gene
    ] = next(
        iter(
            proteins
        )
    )


for x in bridge:

    c = coord[
        x[
            "key"
        ]
    ]

    x[
        "TMD_length"
    ] = float(
        c[
            "length"
        ]
    )

    protein = protein_by_gene[
        x[
            "gene"
        ]
    ]

    seq = protein[
        c[
            "start0"
        ]:
        c[
            "end0"
        ]
    ]

    if len(seq) != c[
        "length"
    ]:

        raise RuntimeError(
            f"TMD slice length mismatch {x['key']}"
        )

    if all(
        aa in KD
        for aa in seq
    ):

        x[
            "hydropathy"
        ] = float(
            np.mean(
                [
                    KD[
                        aa
                    ]
                    for aa in seq
                ]
            )
        )

    else:

        x[
            "hydropathy"
        ] = None


# ============================================================
# 6. Freeze predictor/covariate standardization references
#    BEFORE considering outcome validity.
# ============================================================

score_mu,score_sd = z_reference(
    [
        x[
            "score"
        ]
        for x in bridge
    ]
)

for x in bridge:

    x[
        "score_z"
    ] = z_apply(
        x[
            "score"
        ],
        score_mu,
        score_sd
    )


COVS = [
    "aa_identity",
    "GC3",
    "CpG",
    "tmd_index",
    "TMD_length",
    "hydropathy",
]


cov_reference = {}

for cov in COVS:

    values = []

    for x in bridge:

        if cov == "tmd_index":
            v = float(
                x[
                    "tmd_index"
                ]
            )
        else:
            v = x[
                cov
            ]

        if (
            v is not None
            and
            math.isfinite(
                float(v)
            )
        ):

            values.append(
                float(v)
            )

    mu,sd = z_reference(
        values
    )

    cov_reference[
        cov
    ] = (
        mu,
        sd
    )

    for x in bridge:

        if cov == "tmd_index":
            v = float(
                x[
                    "tmd_index"
                ]
            )
        else:
            v = x[
                cov
            ]

        x[
            cov
            +
            "_z"
        ] = (
            None
            if v is None
            else
            z_apply(
                float(v),
                mu,
                sd
            )
        )


# ============================================================
# 7. Construct P3 outcomes from frozen formula
# ============================================================

def dval(
    row,
    lib,
    side,
):
    return fnum(
        row[
            f"{lib}_{side}_density"
        ]
    )


for x in bridge:

    row = x[
        "density"
    ]

    x[
        "outcomes"
    ] = {}

    for loc in LOCATIONS:

        rep_results = []

        for pd_lib,in_lib in PAIRS[
            loc
        ]:

            R_pd = dval(
                row,
                pd_lib,
                "real"
            )

            R_in = dval(
                row,
                in_lib,
                "real"
            )

            P_pd = dval(
                row,
                pd_lib,
                "pseudo"
            )

            P_in = dval(
                row,
                in_lib,
                "pseudo"
            )

            vals = [
                R_pd,
                R_in,
                P_pd,
                P_in,
            ]

            valid = (
                all(
                    v is not None
                    and v > 0
                    for v in vals
                )
            )

            if not valid:

                rep_results.append(
                    None
                )

                continue

            E_real = math.log2(
                R_pd
                /
                R_in
            )

            E_pseudo = math.log2(
                P_pd
                /
                P_in
            )

            D = (
                E_real
                -
                E_pseudo
            )

            rep_results.append({
                "E_real":
                    E_real,

                "E_pseudo":
                    E_pseudo,

                "D":
                    D,
            })


        if all(
            z is not None
            for z in rep_results
        ):

            r1,r2 = rep_results

            x[
                "outcomes"
            ][
                loc
            ] = {
                "D_rep1":
                    r1[
                        "D"
                    ],

                "D_rep2":
                    r2[
                        "D"
                    ],

                "D_mean":
                    (
                        r1[
                            "D"
                        ]
                        +
                        r2[
                            "D"
                        ]
                    ) / 2,

                "E_real_mean":
                    (
                        r1[
                            "E_real"
                        ]
                        +
                        r2[
                            "E_real"
                        ]
                    ) / 2,

                "E_pseudo_mean":
                    (
                        r1[
                            "E_pseudo"
                        ]
                        +
                        r2[
                            "E_pseudo"
                        ]
                    ) / 2,
            }

        else:

            x[
                "outcomes"
            ][
                loc
            ] = None


# ============================================================
# 8. Stacked primary dataset
# ============================================================

stacked = []

for x in bridge:

    for loc in LOCATIONS:

        outcome = x[
            "outcomes"
        ][
            loc
        ]

        if outcome is None:
            continue

        stacked.append({
            "gene":
                x[
                    "gene"
                ],

            "transcript":
                x[
                    "transcript"
                ],

            "tmd_index":
                x[
                    "tmd_index"
                ],

            "unit":
                x[
                    "paired_real_unit_id"
                ],

            "location":
                loc,

            "score":
                x[
                    "score"
                ],

            "score_z":
                x[
                    "score_z"
                ],

            "D_mean":
                outcome[
                    "D_mean"
                ],

            "D_rep1":
                outcome[
                    "D_rep1"
                ],

            "D_rep2":
                outcome[
                    "D_rep2"
                ],

            "E_real_mean":
                outcome[
                    "E_real_mean"
                ],

            "E_pseudo_mean":
                outcome[
                    "E_pseudo_mean"
                ],

            **{
                cov:
                    x[
                        cov
                    ]
                    if cov
                    !=
                    "tmd_index"
                    else
                    float(
                        x[
                            "tmd_index"
                        ]
                    )

                for cov in COVS
            },

            **{
                cov
                +
                "_z":
                    x[
                        cov
                        +
                        "_z"
                    ]

                for cov in COVS
            },
        })


# ============================================================
# Design builders
# ============================================================

def unadjusted_design(
    rows,
    outcome_name
):
    names = [
        "intercept_ER",
        "intercept_mitochondria",
        "intercept_cytosol",
        "score_ER",
        "score_mitochondria",
        "score_cytosol",
    ]

    X = []
    y = []
    cl = []

    for r in rows:

        loc_i = LOCATIONS.index(
            r[
                "location"
            ]
        )

        z = [0.0] * 6

        z[
            loc_i
        ] = 1.0

        z[
            3
            +
            loc_i
        ] = r[
            "score_z"
        ]

        X.append(
            z
        )

        y.append(
            r[
                outcome_name
            ]
        )

        cl.append(
            r[
                "gene"
            ]
        )

    return (
        np.asarray(
            X
        ),
        np.asarray(
            y
        ),
        np.asarray(
            cl,
            dtype=object
        ),
        names
    )


def adjusted_design(
    rows
):

    clean = []

    for r in rows:

        if all(
            r[
                cov
                +
                "_z"
            ]
            is not None

            for cov in COVS
        ):

            clean.append(
                r
            )

    block_names = [
        "intercept",
        "score",
    ] + COVS

    names = []

    for loc in LOCATIONS:

        for b in block_names:

            names.append(
                f"{b}_{loc}"
            )

    B = len(
        block_names
    )

    X = []
    y = []
    cl = []

    for r in clean:

        loc_i = LOCATIONS.index(
            r[
                "location"
            ]
        )

        z = [
            0.0
        ] * (
            len(
                LOCATIONS
            )
            *
            B
        )

        off = (
            loc_i
            *
            B
        )

        z[
            off
        ] = 1.0

        z[
            off + 1
        ] = r[
            "score_z"
        ]

        for j,cov in enumerate(
            COVS,
            start=2
        ):

            z[
                off + j
            ] = r[
                cov
                +
                "_z"
            ]

        X.append(
            z
        )

        y.append(
            r[
                "D_mean"
            ]
        )

        cl.append(
            r[
                "gene"
            ]
        )

    return (
        clean,
        np.asarray(
            X
        ),
        np.asarray(
            y
        ),
        np.asarray(
            cl,
            dtype=object
        ),
        names,
        B
    )


# ============================================================
# 9. PRIMARY
# ============================================================

X,y,cl,names = unadjusted_design(
    stacked,
    "D_mean"
)

primary = fit_cluster_ols(
    X,
    y,
    cl,
    names
)

primary_slope_idx = [
    3,
    4,
    5,
]

R = np.zeros(
    (
        3,
        len(
            names
        )
    )
)

for j,idx in enumerate(
    primary_slope_idx
):

    R[
        j,
        idx
    ] = 1.0


primary_W,primary_df,primary_p = (
    wald_test(
        primary,
        R
    )
)


# Compartment specificity:
# ER = mito = cytosol
Rc = np.zeros(
    (
        2,
        len(
            names
        )
    )
)

Rc[
    0,
    3
] = 1

Rc[
    0,
    4
] = -1

Rc[
    1,
    3
] = 1

Rc[
    1,
    5
] = -1

comp_W,comp_df,comp_p = (
    wald_test(
        primary,
        Rc
    )
)


primary_coef = coefficient_stats(
    primary,
    primary_slope_idx
)


# ============================================================
# 10. ADJUSTED
# ============================================================

(
    adjusted_rows,
    X,
    y,
    cl,
    adj_names,
    B
) = adjusted_design(
    stacked
)

adjusted = fit_cluster_ols(
    X,
    y,
    cl,
    adj_names
)

adj_score_idx = [
    loc_i * B + 1
    for loc_i in range(
        3
    )
]

Ra = np.zeros(
    (
        3,
        len(
            adj_names
        )
    )
)

for j,idx in enumerate(
    adj_score_idx
):

    Ra[
        j,
        idx
    ] = 1.0


adj_W,adj_df,adj_p = (
    wald_test(
        adjusted,
        Ra
    )
)

adjusted_coef = coefficient_stats(
    adjusted,
    adj_score_idx
)


# ============================================================
# 11. WITHIN-PROTEIN
# Primary within test is unadjusted,
# exactly matching the frozen primary-within specification.
# ============================================================

within_rows = []

for loc in LOCATIONS:

    locrows = [
        r
        for r in stacked
        if r[
            "location"
        ] == loc
    ]

    by_gene = defaultdict(
        list
    )

    for r in locrows:

        by_gene[
            r[
                "gene"
            ]
        ].append(
            r
        )

    for gene,rows in by_gene.items():

        if len(
            rows
        ) < 2:
            continue

        mean_y = float(
            np.mean(
                [
                    r[
                        "D_mean"
                    ]
                    for r in rows
                ]
            )
        )

        mean_s = float(
            np.mean(
                [
                    r[
                        "score_z"
                    ]
                    for r in rows
                ]
            )
        )

        for r in rows:

            within_rows.append({
                **r,

                "within_y":
                    r[
                        "D_mean"
                    ]
                    -
                    mean_y,

                "within_score":
                    r[
                        "score_z"
                    ]
                    -
                    mean_s,
            })


within_names = [
    "within_score_ER",
    "within_score_mitochondria",
    "within_score_cytosol",
]

X = []
y = []
cl = []

for r in within_rows:

    z = [
        0.0,
        0.0,
        0.0,
    ]

    loc_i = LOCATIONS.index(
        r[
            "location"
        ]
    )

    z[
        loc_i
    ] = r[
        "within_score"
    ]

    X.append(
        z
    )

    y.append(
        r[
            "within_y"
        ]
    )

    cl.append(
        r[
            "gene"
        ]
    )


within_fit = fit_cluster_ols(
    X,
    y,
    cl,
    within_names
)

Rw = np.eye(
    3
)

within_W,within_df,within_p = (
    wald_test(
        within_fit,
        Rw
    )
)

within_coef = coefficient_stats(
    within_fit,
    [
        0,
        1,
        2
    ]
)


# ============================================================
# 12. Replicate sensitivities
# ============================================================

rep_results = {}

for outcome_name in [
    "D_rep1",
    "D_rep2",
]:

    X,y,cl,names_rep = unadjusted_design(
        stacked,
        outcome_name
    )

    fit = fit_cluster_ols(
        X,
        y,
        cl,
        names_rep
    )

    Rr = np.zeros(
        (
            3,
            len(
                names_rep
            )
        )
    )

    for j,idx in enumerate(
        [
            3,
            4,
            5
        ]
    ):

        Rr[
            j,
            idx
        ] = 1

    W,df,p = wald_test(
        fit,
        Rr
    )

    rep_results[
        outcome_name
    ] = {
        "fit":
            fit,

        "W":
            W,

        "df":
            df,

        "p":
            p,

        "coef":
            coefficient_stats(
                fit,
                [
                    3,
                    4,
                    5
                ]
            ),
    }


# ============================================================
# 13. Pseudo negative control
# ============================================================

X,y,cl,names_pseudo = unadjusted_design(
    stacked,
    "E_pseudo_mean"
)

pseudo_fit = fit_cluster_ols(
    X,
    y,
    cl,
    names_pseudo
)

Rp = np.zeros(
    (
        3,
        len(
            names_pseudo
        )
    )
)

for j,idx in enumerate(
    [
        3,
        4,
        5
    ]
):

    Rp[
        j,
        idx
    ] = 1

pseudo_W,pseudo_df,pseudo_p = (
    wald_test(
        pseudo_fit,
        Rp
    )
)

pseudo_coef = coefficient_stats(
    pseudo_fit,
    [
        3,
        4,
        5
    ]
)


# ============================================================
# 14. Real component companion
# ============================================================

X,y,cl,names_real = unadjusted_design(
    stacked,
    "E_real_mean"
)

real_fit = fit_cluster_ols(
    X,
    y,
    cl,
    names_real
)

Rr = np.zeros(
    (
        3,
        len(
            names_real
        )
    )
)

for j,idx in enumerate(
    [
        3,
        4,
        5
    ]
):

    Rr[
        j,
        idx
    ] = 1

real_W,real_df,real_p = (
    wald_test(
        real_fit,
        Rr
    )
)

real_coef = coefficient_stats(
    real_fit,
    [
        3,
        4,
        5
    ]
)


# ============================================================
# 15. Per-location N
# ============================================================

location_n = {}

location_genes = {}

for loc in LOCATIONS:

    rr = [
        x
        for x in stacked
        if x[
            "location"
        ] == loc
    ]

    location_n[
        loc
    ] = len(
        rr
    )

    location_genes[
        loc
    ] = len({
        x[
            "gene"
        ]
        for x in rr
    })


within_location_n = {}

within_location_genes = {}

for loc in LOCATIONS:

    rr = [
        x
        for x in within_rows
        if x[
            "location"
        ] == loc
    ]

    within_location_n[
        loc
    ] = len(
        rr
    )

    within_location_genes[
        loc
    ] = len({
        x[
            "gene"
        ]
        for x in rr
    })


# ============================================================
# 16. Frozen adjudication
# ============================================================

no_primary_ci = all(
    not (
        r[
            "CI_2.5"
        ] > 0
        or
        r[
            "CI_97.5"
        ] < 0
    )

    for r in primary_coef
)


STOP = (
    primary_p >= 0.05
    and
    no_primary_ci
    and
    within_p >= 0.05
)


STRONG = (
    primary_p < 0.05
    and
    adj_p < 0.05
    and
    within_p < 0.05
)


if STOP:

    verdict = (
        "STOP_PRIMARY_HYPOTHESIS"
    )

elif STRONG:

    verdict = (
        "STRONG_BRIDGE_SUPPORT"
    )

else:

    verdict = (
        "INTERMEDIATE_BRIDGE_RESULT"
    )


# ============================================================
# 17. Output model table
# ============================================================

model_rows = []


def add_model_rows(
    label,
    coef_rows,
):

    for i,row in enumerate(
        coef_rows
    ):

        loc = LOCATIONS[
            i
        ]

        model_rows.append({
            "model":
                label,

            "location":
                loc,

            "beta_per_1SD_evolution_score":
                row[
                    "beta"
                ],

            "cluster_robust_SE":
                row[
                    "SE"
                ],

            "CI_2.5":
                row[
                    "CI_2.5"
                ],

            "CI_97.5":
                row[
                    "CI_97.5"
                ],

            "p":
                row[
                    "p"
                ],

            "Holm_p":
                row[
                    "Holm_p"
                ],
        })


add_model_rows(
    "PRIMARY_D_MEAN",
    primary_coef
)

add_model_rows(
    "ADJUSTED_D_MEAN",
    adjusted_coef
)

add_model_rows(
    "WITHIN_PROTEIN_D_MEAN",
    within_coef
)

add_model_rows(
    "REPLICATE1_D",
    rep_results[
        "D_rep1"
    ][
        "coef"
    ]
)

add_model_rows(
    "REPLICATE2_D",
    rep_results[
        "D_rep2"
    ][
        "coef"
    ]
)

add_model_rows(
    "PSEUDO_CONTROL_E_PSEUDO_MEAN",
    pseudo_coef
)

add_model_rows(
    "REAL_COMPONENT_E_REAL_MEAN",
    real_coef
)


with open(
    "P5_STAGE1B_FIRST_BRIDGE_MODEL_RESULTS.tsv",
    "w",
    newline="",
    encoding="utf-8"
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=list(
            model_rows[
                0
            ].keys()
        ),
        delimiter="\t",
        lineterminator="\n"
    )

    w.writeheader()
    w.writerows(
        model_rows
    )


# ============================================================
# 18. Reveal dataset
# ============================================================

bridge_output = []

for r in stacked:

    bridge_output.append({
        "human_ensembl_gene":
            r[
                "gene"
            ],

        "human_ensembl_transcript":
            r[
                "transcript"
            ],

        "tmd_index":
            r[
                "tmd_index"
            ],

        "paired_real_unit_id":
            r[
                "unit"
            ],

        "location":
            r[
                "location"
            ],

        "evolution_score":
            r[
                "score"
            ],

        "evolution_score_z":
            r[
                "score_z"
            ],

        "D_rep1":
            r[
                "D_rep1"
            ],

        "D_rep2":
            r[
                "D_rep2"
            ],

        "D_mean":
            r[
                "D_mean"
            ],

        "E_real_mean":
            r[
                "E_real_mean"
            ],

        "E_pseudo_mean":
            r[
                "E_pseudo_mean"
            ],

        "median_window_aa_identity":
            r[
                "aa_identity"
            ],

        "human_window_GC3":
            r[
                "GC3"
            ],

        "human_window_CpG_density":
            r[
                "CpG"
            ],

        "TMD_length_aa":
            r[
                "TMD_length"
            ],

        "TMD_mean_Kyte_Doolittle_hydropathy":
            r[
                "hydropathy"
            ],
    })


bridge_output.sort(
    key=lambda x:(
        x[
            "human_ensembl_gene"
        ],
        int(
            x[
                "tmd_index"
            ]
        ),
        x[
            "location"
        ],
    )
)


with open(
    "P5_STAGE1B_FIRST_BRIDGE_REVEALED_DATA.tsv",
    "w",
    newline="",
    encoding="utf-8"
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=list(
            bridge_output[
                0
            ].keys()
        ),
        delimiter="\t",
        lineterminator="\n"
    )

    w.writeheader()
    w.writerows(
        bridge_output
    )


# ============================================================
# 19. Human-readable first reveal summary
# ============================================================

def fmt(x):
    if x is None:
        return "NA"

    if isinstance(
        x,
        str
    ):
        return x

    return f"{x:.8g}"


lines = [
    "P5 STAGE1B FIRST FROZEN BRIDGE REVEAL",
    "=" * 76,
    "",
    "FIRST_FORMAL_P3_PHENOTYPE_REVEAL = TRUE",
    "",
    f"PRELOCKED_BRIDGE_TMD = {len(bridge)}",
    f"VALID_STACKED_TMD_LOCATION_ROWS = {len(stacked)}",
    "",
]

for loc in LOCATIONS:

    lines.append(
        f"{loc}_VALID_TMD = "
        f"{location_n[loc]}"
    )

    lines.append(
        f"{loc}_VALID_GENES = "
        f"{location_genes[loc]}"
    )


lines += [
    "",
    "PRIMARY MODEL — D_mean",
    "-" * 76,
    f"N = {primary['n']}",
    f"GENE_CLUSTERS = {primary['clusters']}",
    f"JOINT_WALD = {fmt(primary_W)}",
    f"JOINT_DF = {primary_df}",
    f"JOINT_P = {fmt(primary_p)}",
    "",
]


for loc,row in zip(
    LOCATIONS,
    primary_coef
):

    lines.append(
        f"{loc}: "
        f"beta={fmt(row['beta'])} "
        f"SE={fmt(row['SE'])} "
        f"CI=[{fmt(row['CI_2.5'])},"
        f"{fmt(row['CI_97.5'])}] "
        f"p={fmt(row['p'])} "
        f"Holm_p={fmt(row['Holm_p'])}"
    )


lines += [
    "",
    "COMPARTMENT SPECIFICITY",
    "-" * 76,
    f"WALD = {fmt(comp_W)}",
    f"DF = {comp_df}",
    f"P = {fmt(comp_p)}",
    "",
    "ADJUSTED MODEL",
    "-" * 76,
    f"N = {adjusted['n']}",
    f"GENE_CLUSTERS = {adjusted['clusters']}",
    f"JOINT_WALD = {fmt(adj_W)}",
    f"JOINT_DF = {adj_df}",
    f"JOINT_P = {fmt(adj_p)}",
    "",
]


for loc,row in zip(
    LOCATIONS,
    adjusted_coef
):

    lines.append(
        f"{loc}: "
        f"beta={fmt(row['beta'])} "
        f"SE={fmt(row['SE'])} "
        f"CI=[{fmt(row['CI_2.5'])},"
        f"{fmt(row['CI_97.5'])}] "
        f"p={fmt(row['p'])} "
        f"Holm_p={fmt(row['Holm_p'])}"
    )


lines += [
    "",
    "WITHIN-PROTEIN MODEL",
    "-" * 76,
    f"N = {within_fit['n']}",
    f"GENE_CLUSTERS = {within_fit['clusters']}",
    f"JOINT_WALD = {fmt(within_W)}",
    f"JOINT_DF = {within_df}",
    f"JOINT_P = {fmt(within_p)}",
]


for loc in LOCATIONS:

    lines.append(
        f"{loc}_WITHIN_TMD = "
        f"{within_location_n[loc]}"
    )

    lines.append(
        f"{loc}_WITHIN_GENES = "
        f"{within_location_genes[loc]}"
    )


lines.append("")


for loc,row in zip(
    LOCATIONS,
    within_coef
):

    lines.append(
        f"{loc}: "
        f"beta={fmt(row['beta'])} "
        f"SE={fmt(row['SE'])} "
        f"CI=[{fmt(row['CI_2.5'])},"
        f"{fmt(row['CI_97.5'])}] "
        f"p={fmt(row['p'])} "
        f"Holm_p={fmt(row['Holm_p'])}"
    )


lines += [
    "",
    "REPLICATE SENSITIVITY",
    "-" * 76,
    f"REP1_JOINT_P = "
    f"{fmt(rep_results['D_rep1']['p'])}",
    f"REP2_JOINT_P = "
    f"{fmt(rep_results['D_rep2']['p'])}",
    "",
]


for outcome_name in [
    "D_rep1",
    "D_rep2",
]:

    lines.append(
        outcome_name + ":"
    )

    for loc,row in zip(
        LOCATIONS,
        rep_results[
            outcome_name
        ][
            "coef"
        ]
    ):

        lines.append(
            f"  {loc}: "
            f"beta={fmt(row['beta'])} "
            f"p={fmt(row['p'])}"
        )


lines += [
    "",
    "PSEUDO NEGATIVE CONTROL",
    "-" * 76,
    f"JOINT_WALD = {fmt(pseudo_W)}",
    f"JOINT_DF = {pseudo_df}",
    f"JOINT_P = {fmt(pseudo_p)}",
    "",
]


for loc,row in zip(
    LOCATIONS,
    pseudo_coef
):

    lines.append(
        f"{loc}: "
        f"beta={fmt(row['beta'])} "
        f"SE={fmt(row['SE'])} "
        f"p={fmt(row['p'])}"
    )


lines += [
    "",
    "REAL COMPONENT COMPANION",
    "-" * 76,
    f"JOINT_WALD = {fmt(real_W)}",
    f"JOINT_DF = {real_df}",
    f"JOINT_P = {fmt(real_p)}",
    "",
]


for loc,row in zip(
    LOCATIONS,
    real_coef
):

    lines.append(
        f"{loc}: "
        f"beta={fmt(row['beta'])} "
        f"SE={fmt(row['SE'])} "
        f"p={fmt(row['p'])}"
    )


lines += [
    "",
    "FROZEN ADJUDICATION",
    "=" * 76,
    f"PRIMARY_JOINT_P = {fmt(primary_p)}",
    f"ADJUSTED_JOINT_P = {fmt(adj_p)}",
    f"WITHIN_PROTEIN_JOINT_P = {fmt(within_p)}",
    f"PSEUDO_CONTROL_JOINT_P = {fmt(pseudo_p)}",
    f"COMPARTMENT_SPECIFICITY_P = {fmt(comp_p)}",
    "",
    f"P5_PRIMARY_HYPOTHESIS = {verdict}",
    "",
    "PRIMARY_WINDOW_RETUNING_AUTHORIZED = FALSE",
    "PREDICTOR_REDEFINITION_AUTHORIZED = FALSE",
    "POST_HOC_COMPARTMENT_SELECTION_AUTHORIZED = FALSE",
    "CAUSAL_CLAIM_AUTHORIZED = FALSE",
]


Path(
    "P5_STAGE1B_FIRST_FROZEN_BRIDGE_REVEAL_SUMMARY.txt"
).write_text(
    "\n".join(
        lines
    ) + "\n",
    encoding="utf-8"
)


print(
    "\n".join(
        lines
    ),
    flush=True
)
