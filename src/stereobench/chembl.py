"""ChEMBL lookups for the three things StereoAtlas does not carry.

StereoAtlas gives compound ids and a mean activity. It gives no per-member
SMILES, no `assay_chembl_id`, and no activity qualifier. All three are needed
before a candidate row can be called a pair:

- **SMILES** so `stereo.classify_pair` can confirm enantiomerism exactly,
  rather than inferring it from descriptor strings that miscall meso compounds.
- **assay id** so both members come from one assay. Same target and same
  activity type is not enough: cross-lab IC50 variance routinely exceeds the
  stereochemical effect, so a pair spanning two papers can show a difference
  that has nothing to do with chirality.
- **qualifier** so censored measurements (`>`, `<`) are excluded. A `>10000 nM`
  reading is a lower bound, not a value, and these are the likeliest source of
  the implausible fold-differences in the archive.

Read-only, public, unauthenticated. Nothing identifying is sent.
"""

from __future__ import annotations

import json
import math
import ssl
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass

import certifi

BASE = "https://www.ebi.ac.uk/chembl/api/data"

# Generic client name. Deliberately carries no contact details: this is a
# third-party service and the user's identity is not its business.
USER_AGENT = "stereobench/0.1 (research benchmark; python-urllib)"

# ChEMBL is a public good running on someone else's budget. Be a good citizen.
REQUEST_DELAY_S = 0.25
MAX_RETRIES = 3

# Qualifiers that denote an actual measurement rather than a bound.
UNCENSORED_RELATIONS = ("=", None, "")


# A framework Python on macOS has no system CA bundle wired into OpenSSL, so
# verification fails with CERTIFICATE_VERIFY_FAILED. Point at certifi's bundle
# explicitly rather than the usual "fix" of disabling verification.
_SSL_CONTEXT = ssl.create_default_context(cafile=certifi.where())


def _get(path: str, params: dict) -> dict:
    url = f"{BASE}/{path}?{urllib.parse.urlencode(params)}"
    last_error: Exception | None = None
    for attempt in range(MAX_RETRIES):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(  # noqa: S310
                request, timeout=60, context=_SSL_CONTEXT
            ) as response:
                payload = json.loads(response.read().decode())
            time.sleep(REQUEST_DELAY_S)
            return payload
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            last_error = exc
            time.sleep(2.0 * (attempt + 1))
    raise RuntimeError(f"ChEMBL request failed after {MAX_RETRIES} tries: {url}") from last_error


def fetch_smiles(chembl_ids: list[str], *, batch_size: int = 40) -> dict[str, str]:
    """Isomeric SMILES for a list of compounds, batched.

    ChEMBL's canonical_smiles field is isomeric -- stereo tags are present --
    which is what `classify_pair` needs. Compounds with no structure on record
    are simply absent from the result.
    """
    out: dict[str, str] = {}
    for start in range(0, len(chembl_ids), batch_size):
        batch = chembl_ids[start : start + batch_size]
        payload = _get(
            "molecule.json",
            {
                "molecule_chembl_id__in": ",".join(batch),
                "limit": len(batch),
                "only": "molecule_chembl_id,molecule_structures",
            },
        )
        for molecule in payload.get("molecules", []):
            structures = molecule.get("molecule_structures") or {}
            smiles = structures.get("canonical_smiles")
            if smiles:
                out[molecule["molecule_chembl_id"]] = smiles
    return out


@dataclass(frozen=True)
class Measurement:
    """One activity record, kept only if it is a real number from one assay."""

    assay_chembl_id: str
    standard_type: str
    standard_value: float
    standard_units: str
    standard_relation: str | None
    document_chembl_id: str | None

    @property
    def is_censored(self) -> bool:
        return self.standard_relation not in UNCENSORED_RELATIONS


def fetch_activities(
    molecule_chembl_id: str, target_chembl_id: str, standard_type: str
) -> list[Measurement]:
    """Every nM measurement of one compound against one target, by assay."""
    payload = _get(
        "activity.json",
        {
            "molecule_chembl_id": molecule_chembl_id,
            "target_chembl_id": target_chembl_id,
            "standard_type": standard_type,
            "limit": 1000,
            "only": (
                "assay_chembl_id,standard_type,standard_value,standard_units,"
                "standard_relation,document_chembl_id"
            ),
        },
    )
    out = []
    for activity in payload.get("activities", []):
        value, units = activity.get("standard_value"), activity.get("standard_units")
        if value is None or units != "nM":
            continue
        out.append(
            Measurement(
                assay_chembl_id=activity.get("assay_chembl_id"),
                standard_type=activity.get("standard_type"),
                standard_value=float(value),
                standard_units=units,
                standard_relation=activity.get("standard_relation"),
                document_chembl_id=activity.get("document_chembl_id"),
            )
        )
    return out


@dataclass(frozen=True)
class SharedAssay:
    """The chosen assay, plus what the choice was made from.

    `n_shared_assays` and `n_measurements` exist so the selection can be
    audited after the fact. A pair resolved from a single candidate assay was
    not selected at all; one resolved from six was, and the reader deserves to
    know which they are looking at.
    """

    eutomer: "Measurement"
    distomer: "Measurement"
    n_shared_assays: int
    n_measurements: int


def find_shared_assay(
    eutomer: list[Measurement], distomer: list[Measurement]
) -> SharedAssay | None:
    """One assay measuring both members, uncensored on both sides.

    This is the filter that makes a pair worth anything. Where an assay reports
    a member more than once, the measurements are taken at their geometric mean
    -- these are log-scale quantities, so an arithmetic mean would be biased
    upward by the weakest replicate.

    **Where several assays qualify, the best-supported one wins**, meaning the
    one backed by the most measurements across both members. An earlier version
    took the assay with the LARGEST separation, on the reasoning that a pilot
    wants clear cases. That reasoning is wrong in a way worth spelling out: the
    hypothesis under test is that enantiomers separate, so choosing the assay
    that separates them most puts a thumb on the scale in precisely the
    direction of the expected result. The effect is not hypothetical -- one
    pair in the set reads 0.87 or 4.66 depending only on which assay is picked,
    a 6000-fold difference in the reported effect arising from the selection
    rule rather than from chemistry.

    Replicate count is a defensible criterion because it is independent of the
    outcome: it says how well measured the pair is, not how much it differs.
    Ties break on the assay id, so the choice is deterministic across runs.
    """

    def by_assay(measurements: list[Measurement]) -> dict[str, list[Measurement]]:
        grouped: dict[str, list[Measurement]] = {}
        for m in measurements:
            if m.is_censored or not m.assay_chembl_id or m.standard_value <= 0:
                continue
            grouped.setdefault(m.assay_chembl_id, []).append(m)
        return grouped

    eu_by_assay, di_by_assay = by_assay(eutomer), by_assay(distomer)
    shared = set(eu_by_assay) & set(di_by_assay)
    if not shared:
        return None

    def collapse(measurements: list[Measurement]) -> Measurement:
        if len(measurements) == 1:
            return measurements[0]

        geometric = math.exp(
            sum(math.log(m.standard_value) for m in measurements) / len(measurements)
        )
        first = measurements[0]
        return Measurement(
            assay_chembl_id=first.assay_chembl_id,
            standard_type=first.standard_type,
            standard_value=geometric,
            standard_units=first.standard_units,
            standard_relation="=",
            document_chembl_id=first.document_chembl_id,
        )

    def support(assay_id: str) -> int:
        return len(eu_by_assay[assay_id]) + len(di_by_assay[assay_id])

    # Outcome-independent criterion, then a deterministic tiebreak. Sorting on
    # the id rather than relying on set order matters: set iteration order
    # varies between runs, so a tie would otherwise resolve differently each
    # time and the pair set would not reproduce.
    best = max(sorted(shared), key=support)

    return SharedAssay(
        eutomer=collapse(eu_by_assay[best]),
        distomer=collapse(di_by_assay[best]),
        n_shared_assays=len(shared),
        n_measurements=support(best),
    )


def fetch_target_uniprot(target_chembl_id: str) -> str | None:
    """The single-protein UniProt accession behind a ChEMBL target, if there is one.

    Returns None for targets that are not a single protein -- complexes, cell
    lines, protein families. Those cannot be modelled as one chain and have no
    place in this benchmark.
    """
    payload = _get(
        "target.json",
        {"target_chembl_id": target_chembl_id, "limit": 1, "only": "target_components,target_type"},
    )
    targets = payload.get("targets", [])
    if not targets or targets[0].get("target_type") != "SINGLE PROTEIN":
        return None
    components = targets[0].get("target_components") or []
    if len(components) != 1:
        return None
    return components[0].get("accession")


def fetch_uniprot_sequence(accession: str) -> str | None:
    """Canonical sequence from UniProt. Separate service, same read-only terms."""
    url = f"https://rest.uniprot.org/uniprotkb/{urllib.parse.quote(accession)}.fasta"
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(  # noqa: S310
            request, timeout=60, context=_SSL_CONTEXT
        ) as response:
            text = response.read().decode()
    except (urllib.error.URLError, TimeoutError):
        return None
    time.sleep(REQUEST_DELAY_S)
    lines = [line for line in text.splitlines() if line and not line.startswith(">")]
    return "".join(lines) or None
