from pathlib import Path
import os
from collections import defaultdict
import gzip
import hashlib

import numpy as np
import pandas as pd
import statsmodels.api as sm


SOURCE_PARENT = Path(os.environ["P5_DATA_PARENT"])

ROOT = (
    SOURCE_PARENT
    / "paper5_evolution_spatial_bridge"
    / "05_GSE297497_MPT_BRIDGE"
)

MPT_DIR = (
    ROOT
    / "13_MPT_OUTCOME_REVEAL"
)

MPT_TABLE = (
    MPT_DIR
    / "P5_V2_TMD_OUTCOME_TABLE.tsv"
)

REFGENE = (
    ROOT
    / "09_FORMAL_REFERENCE"
    / "P5_V2_CURRENT_UCSC_hg38_refGene.txt.gz"
)

RUNROOT = (
    ROOT
    / "14_OSTA_SPECIFICITY"
    / "runs"
)

SPEC = (
    ROOT
    / "14_OSTA_SPECIFICITY"
    / "P5_V2_OSTA_SPECIFICITY_FREEZE.yaml"
)

PREREVEAL_SPEC = (
    ROOT
    / "07_PREREVEAL_FREEZE"
    / "P5_V2_PREREVEAL_MPT_MODEL_SPEC.yaml"
)

OCCDIR = Path("osta_occupancy_npz")
OCCDIR.mkdir(exist_ok=True)

EXPECTED_PREREVEAL_SHA = (
    "181ed0bef6bc81e4750f421c1e0b3f79"
    "a7bfd5fc55432d60664a191d7854b912"
)

PROFILES = {
    "OST48_rep1":
        ("SRR33621390", "SRR33621389"),

    "OST48_rep2":
        ("SRR33621388", "SRR33621387"),

    "RPN2_rep1":
        ("SRR33621386", "SRR33621385"),
}

ALL_RUNS = sorted(
    {
        r
        for pair in PROFILES.values()
        for r in pair
    }
)

CANONICAL = {
    *(f"chr{i}" for i in range(1,23)),
    "chrX",
    "chrY",
    "chrM",
}


def sha256_file(path):
    h = hashlib.sha256()

    with open(path, "rb") as f:
        for chunk in iter(
            lambda: f.read(1024*1024),
            b""
        ):
            h.update(chunk)

    return h.hexdigest()


def qc_value(run, key):

    p = (
        RUNROOT
        / run
        / "final"
        / f"{run}.QC.txt"
    )

    d = {}

    for line in p.read_text().splitlines():

        if " = " in line:

            k, v = line.split(
                " = ",
                1
            )

            d[k.strip()] = v.strip()

    return d[key]


def cluster_ols(
    df,
    ycol
):

    d = (
        df[
            [
                ycol,
                "score_z",
                "human_ensembl_gene",
            ]
        ]
        .replace(
            [np.inf, -np.inf],
            np.nan
        )
        .dropna()
        .copy()
    )

    X = sm.add_constant(
        d[["score_z"]].astype(float),
        has_constant="add"
    )

    fit = sm.OLS(
        d[ycol].astype(float),
        X
    ).fit(
        cov_type="cluster",
        cov_kwds={
            "groups":
                d[
                    "human_ensembl_gene"
                ],
            "use_correction":
                True,
        },
        use_t=True,
    )

    ci = fit.conf_int().loc[
        "score_z"
    ]

    return {
        "N":
            len(d),

        "genes":
            d[
                "human_ensembl_gene"
            ].nunique(),

        "beta":
            float(
                fit.params[
                    "score_z"
                ]
            ),

        "se":
            float(
                fit.bse[
                    "score_z"
                ]
            ),

        "ci_low":
            float(
                ci.iloc[0]
            ),

        "ci_high":
            float(
                ci.iloc[1]
            ),

        "p":
            float(
                fit.pvalues[
                    "score_z"
                ]
            ),
    }


print("=" * 80)
print("P5-V2 OST-A SPECIFICITY REVEAL")
print("=" * 80)


# ============================================================
# 1. LOCKS
# ============================================================

pre_sha = sha256_file(
    PREREVEAL_SPEC
)

print(
    "PREREVEAL_SPEC_SHA256 =",
    pre_sha
)

assert (
    pre_sha
    == EXPECTED_PREREVEAL_SHA
)

assert SPEC.exists()

print(
    "SPECIFICITY_FREEZE_PRESENT = TRUE"
)


# ============================================================
# 2. REQUIRE ALL 6 OST-A RUNS
# ============================================================

for run in ALL_RUNS:

    assert (
        RUNROOT
        / run
        / ".DONE"
    ).exists()

    assert (
        RUNROOT
        / run
        / "final"
        / f"{run}.asite_plus15.tsv"
    ).is_file()

print(
    "FORMAL_OSTA_RUNS_PRESENT = 6"
)


# ============================================================
# 3. LOAD MPT TMD TABLE
# ============================================================

m = pd.read_csv(
    MPT_TABLE,
    sep="\t"
)

assert len(m) == 5538

if "author_tx" not in m.columns:

    m["author_tx"] = (
        m[
            "author_transcript_names"
        ]
        .astype(str)
        .str.strip()
    )

targets = sorted(
    m["author_tx"].unique()
)

assert len(targets) == 1400

target_set = set(targets)

print(
    "PRELOCKED_TMD =",
    len(m)
)

print(
    "TARGET_TRANSCRIPTS =",
    len(targets)
)


# ============================================================
# 4. REFSEQ TRANSCRIPT MODELS
# ============================================================

ref_rows = defaultdict(list)

with gzip.open(
    REFGENE,
    "rt"
) as fh:

    for line in fh:

        f = line.rstrip(
            "\n"
        ).split("\t")

        if len(f) < 11:
            continue

        tx = f[1]

        if (
            tx not in target_set
            or f[2] not in CANONICAL
        ):
            continue

        ref_rows[tx].append(f)


models = {}

for tx in targets:

    rr = ref_rows[tx]

    assert len(rr) == 1

    f = rr[0]

    chrom = f[2]
    strand = f[3]

    cds_start = int(f[6])
    cds_end = int(f[7])

    exon_starts = [
        int(x)
        for x in f[9]
        .rstrip(",")
        .split(",")
        if x
    ]

    exon_ends = [
        int(x)
        for x in f[10]
        .rstrip(",")
        .split(",")
        if x
    ]

    genomic_exons = list(
        zip(
            exon_starts,
            exon_ends,
        )
    )

    ordered = (
        genomic_exons
        if strand == "+"
        else list(
            reversed(
                genomic_exons
            )
        )
    )

    exon_maps = []

    t_cursor = 0

    cds_intervals = []

    for gs, ge in ordered:

        exon_maps.append(
            {
                "gstart": gs,
                "gend": ge,
                "tstart": t_cursor,
            }
        )

        is0 = max(
            gs,
            cds_start
        )

        ie0 = min(
            ge,
            cds_end
        )

        if is0 < ie0:

            if strand == "+":

                ts = (
                    t_cursor
                    + is0 - gs
                )

                te = (
                    t_cursor
                    + ie0 - gs
                )

            else:

                ts = (
                    t_cursor
                    + ge - ie0
                )

                te = (
                    t_cursor
                    + ge - is0
                )

            cds_intervals.append(
                (ts, te)
            )

        t_cursor += (
            ge - gs
        )

    cds_t_start = min(
        x[0]
        for x in cds_intervals
    )

    cds_t_end = max(
        x[1]
        for x in cds_intervals
    )

    protein_len = (
        cds_t_end
        - cds_t_start
    ) // 3

    # RefSeq CDS contains stop codon:
    # previous MPT audit showed +3 nt in all 1400.
    if (
        cds_t_end
        - cds_t_start
    ) % 3 == 0:

        # Remove terminal stop codon
        protein_len -= 1

    models[tx] = {
        "chrom":
            chrom,

        "strand":
            strand,

        "tx_len":
            t_cursor,

        "cds_t_start":
            cds_t_start,

        "protein_len":
            protein_len,

        "exons":
            exon_maps,
    }


# ============================================================
# 5. OST-A OCCUPANCY
# ============================================================

def build_occupancy(run):

    cache = (
        OCCDIR
        / f"{run}.target_codon_rpm.npz"
    )

    if cache.exists():

        z = np.load(
            cache,
            allow_pickle=False
        )

        out = {
            k: z[k]
            for k in z.files
        }

        print(
            run,
            "CACHE = REUSED"
        )

        return out


    path = (
        RUNROOT
        / run
        / "final"
        / f"{run}.asite_plus15.tsv"
    )

    total = int(
        qc_value(
            run,
            "ASITE_ASSIGNED_READS"
        )
    )

    print(
        run,
        "ASITE =",
        total
    )

    df = pd.read_csv(
        path,
        sep="\t",
        usecols=[
            "chrom",
            "strand",
            "asite_0based",
        ]
    )

    grouped = {}

    for (
        chrom,
        strand
    ), g in df.groupby(
        [
            "chrom",
            "strand"
        ],
        sort=False
    ):

        arr = g[
            "asite_0based"
        ].to_numpy(
            dtype=np.int64
        )

        arr.sort()

        grouped[
            (
                str(chrom),
                str(strand)
            )
        ] = arr

    del df

    result = {}

    mapped = 0

    for n, tx in enumerate(
        targets,
        start=1
    ):

        model = models[tx]

        arr = grouped.get(
            (
                model["chrom"],
                model["strand"]
            )
        )

        chunks = []

        if arr is not None:

            for ex in model["exons"]:

                gs = ex["gstart"]
                ge = ex["gend"]
                ts = ex["tstart"]

                lo = np.searchsorted(
                    arr,
                    gs,
                    side="left"
                )

                hi = np.searchsorted(
                    arr,
                    ge,
                    side="left"
                )

                if hi <= lo:
                    continue

                pos = arr[lo:hi]

                if (
                    model["strand"]
                    == "+"
                ):

                    tc = (
                        ts
                        + pos
                        - gs
                    )

                else:

                    tc = (
                        ts
                        + ge
                        - 1
                        - pos
                    )

                chunks.append(
                    tc.astype(
                        np.int64
                    )
                )

        if chunks:

            coords = np.concatenate(
                chunks
            )

            mapped += len(coords)

            starts = np.maximum(
                coords - 50,
                0
            )

            ends = np.minimum(
                coords + 51,
                model["tx_len"]
            )

            diff = np.zeros(
                model["tx_len"] + 1,
                dtype=np.int64
            )

            np.add.at(
                diff,
                starts,
                1
            )

            np.add.at(
                diff,
                ends,
                -1
            )

            cov = np.cumsum(
                diff[:-1]
            )

        else:

            cov = np.zeros(
                model["tx_len"],
                dtype=np.int64
            )


        cs = model[
            "cds_t_start"
        ]

        ncod = model[
            "protein_len"
        ]

        coding = cov[
            cs:
            cs + 3 * ncod
        ]

        if len(coding) != (
            3 * ncod
        ):
            raise RuntimeError(
                f"{run} {tx}: "
                "coding length error"
            )

        rpm = (
            coding
            .reshape(
                ncod,
                3
            )
            .mean(axis=1)
            * 1_000_000.0
            / total
        ).astype(
            np.float32
        )

        result[tx] = rpm

        if n % 200 == 0:

            print(
                run,
                "TARGETS_DONE =",
                n
            )

    print(
        run,
        "TARGET_MAPPING_EVENTS =",
        mapped
    )

    np.savez_compressed(
        cache,
        **result
    )

    return result


occ = {}

for run in ALL_RUNS:

    occ[run] = build_occupancy(
        run
    )

print(
    "ALL_6_OSTA_OCCUPANCIES_READY = TRUE"
)


# ============================================================
# 6. CONSTRUCT +30...+90 OST-A SCORES
# ============================================================

profile_cols = []

for pname, (
    input_run,
    ip_run
) in PROFILES.items():

    col = (
        pname
        + "_TMD_score"
    )

    profile_cols.append(col)

    vals = []

    for row in m.itertuples(
        index=False
    ):

        tx = row.author_tx

        inp = occ[
            input_run
        ][tx]

        ip = occ[
            ip_run
        ][tx]

        A = float(
            np.mean(inp)
        )

        end = int(
            row.author_end_1based
        )

        plen = len(inp)

        if (
            A > 0
            and end + 90 <= plen
        ):

            s = (
                end + 29
            )

            e = (
                end + 90
            )

            iw = inp[s:e]
            pw = ip[s:e]

            assert len(iw) == 61
            assert len(pw) == 61

            loge = np.log2(
                (
                    pw + A
                )
                /
                (
                    iw + A
                )
            )

            vals.append(
                float(
                    np.mean(loge)
                )
            )

        else:

            vals.append(
                np.nan
            )

    m[col] = vals


m["OSTA_profiles_valid"] = (
    m[
        profile_cols
    ]
    .notna()
    .sum(axis=1)
)

m["OSTA_TMD_score"] = np.where(
    m[
        "OSTA_profiles_valid"
    ] == 3,
    m[
        profile_cols
    ].mean(axis=1),
    np.nan
)


print(
    "OSTA_VALID_TMD =",
    int(
        m[
            "OSTA_TMD_score"
        ]
        .notna()
        .sum()
    )
)


# ============================================================
# 7. SHARED MPT / OST-A UNIVERSE
# ============================================================

shared = m[
    m[
        "primary_MPT_TMD_score"
    ].notna()
    &
    m[
        "OSTA_TMD_score"
    ].notna()
].copy()

print(
    "SHARED_MPT_OSTA_TMD =",
    len(shared)
)

print(
    "SHARED_MPT_OSTA_GENES =",
    shared[
        "human_ensembl_gene"
    ].nunique()
)


shared[
    "MPT_minus_OSTA_score"
] = (
    shared[
        "primary_MPT_TMD_score"
    ]
    -
    shared[
        "OSTA_TMD_score"
    ]
)


# ============================================================
# 8. FORMAL SPECIFICITY TESTS
# ============================================================

mpt = cluster_ols(
    shared,
    "primary_MPT_TMD_score"
)

osta = cluster_ols(
    shared,
    "OSTA_TMD_score"
)

diff = cluster_ols(
    shared,
    "MPT_minus_OSTA_score"
)


# ============================================================
# 9. PROFILE DESCRIPTIVES
# ============================================================

profile_results = {}

for col in profile_cols:

    profile_results[col] = (
        cluster_ols(
            shared,
            col
        )
    )


# ============================================================
# 10. INTERPRETATION
# ============================================================

osta_weaker = (
    abs(
        osta["beta"]
    )
    <
    abs(
        mpt["beta"]
    )
)

difference_sig = (
    diff["p"] < 0.05
)

mpt_same_as_original = (
    np.sign(
        mpt["beta"]
    )
    ==
    np.sign(
        -0.0212987
    )
)

if (
    osta_weaker
    and difference_sig
):

    interpretation = (
        "MPT_SPECIFICITY_SUPPORTED"
    )

elif osta_weaker:

    interpretation = (
        "MPT_SPECIFICITY_DIRECTIONALLY_SUGGESTIVE"
    )

else:

    interpretation = (
        "MPT_SPECIFICITY_NOT_SUPPORTED"
    )


# ============================================================
# 11. OUTPUT
# ============================================================

def fmt(r):

    return (
        "N={N} genes={genes} "
        "beta={beta:.6g} "
        "SE={se:.6g} "
        "95%CI=[{ci_low:.6g},{ci_high:.6g}] "
        "p={p:.6g}"
    ).format(**r)


print()
print("=" * 80)
print("P5-V2 OST-A SPECIFICITY REVEAL")
print("=" * 80)

print(
    "MPT_SHARED:",
    fmt(mpt)
)

print(
    "OSTA_SHARED:",
    fmt(osta)
)

print(
    "MPT_MINUS_OSTA:",
    fmt(diff)
)

print()
print("OST-A PROFILES:")

for col in profile_cols:

    print(
        col + ":",
        fmt(
            profile_results[col]
        )
    )

print()

print(
    "OSTA_ABS_BETA_WEAKER_THAN_MPT =",
    str(
        osta_weaker
    ).upper()
)

print(
    "MPT_MINUS_OSTA_DIFFERENCE_P_LT_0_05 =",
    str(
        difference_sig
    ).upper()
)

print(
    "MPT_SHARED_DIRECTION_MATCHES_ORIGINAL =",
    str(
        mpt_same_as_original
    ).upper()
)

print(
    "OSTA_SPECIFICITY_INTERPRETATION =",
    interpretation
)

print()

print(
    "FORMAL_P5_V2_PRIMARY_ADJUDICATION_REMAINS = "
    "INTERMEDIATE_BRIDGE_RESULT"
)

print(
    "OSTA_MAY_RESCUE_PRIMARY = FALSE"
)

print(
    "CAUSAL_CLAIM_AUTHORIZED = FALSE"
)


m.to_csv(
    "P5_V2_OSTA_TMD_OUTCOME_TABLE.tsv",
    sep="\t",
    index=False
)

pd.DataFrame(
    [
        {
            "analysis":
                "MPT_SHARED",
            **mpt
        },
        {
            "analysis":
                "OSTA_SHARED",
            **osta
        },
        {
            "analysis":
                "MPT_MINUS_OSTA",
            **diff
        },
    ]
).to_csv(
    "P5_V2_OSTA_SPECIFICITY_CORE_RESULTS.tsv",
    sep="\t",
    index=False
)

pd.DataFrame(
    [
        {
            "analysis":
                col,
            **profile_results[col]
        }
        for col in profile_cols
    ]
).to_csv(
    "P5_V2_OSTA_PROFILE_RESULTS.tsv",
    sep="\t",
    index=False
)

print()
print(
    "OSTA_SPECIFICITY_OUTPUTS_WRITTEN = TRUE"
)
