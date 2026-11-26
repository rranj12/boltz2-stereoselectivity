"""Tests for the generated run script and the YAMLs it drives.

These exist because the failures they catch are discovered late and cost GPU
hours: a run script that drops the manifest's seeds, or YAMLs whose arms do
not share an MSA, both produce a run that completes and means nothing.
"""

from __future__ import annotations

import pandas as pd
import pytest

from build_run import build_ccd_yaml, ligand_code, write_run_script


@pytest.fixture
def manifest():
    return pd.DataFrame(
        [
            {"run_name": "pair1__a_native__r0", "boltz_seed": 1000,
             "yaml": "/runs/inputs/a.yaml", "msa": "/data/msa/T1.csv"},
            {"run_name": "pair1__a_native__r1", "boltz_seed": 1001,
             "yaml": "/runs/inputs/b.yaml", "msa": "/data/msa/T1.csv"},
        ]
    )


def test_every_input_gets_its_own_recorded_seed(tmp_path, manifest):
    """The regression: seeds lived in the manifest and nowhere else."""
    script = tmp_path / "run.sh"
    write_run_script(manifest, script, out_root=tmp_path / "out", uses_cached_msa=True)
    text = script.read_text()
    assert "--seed 1000" in text
    assert "--seed 1001" in text


def test_seeds_are_attached_to_the_right_input(tmp_path, manifest):
    script = tmp_path / "run.sh"
    write_run_script(manifest, script, out_root=tmp_path / "out", uses_cached_msa=True)
    first = script.read_text().split("pair1__a_native__r1")[0]
    assert "--seed 1000" in first and "--seed 1001" not in first


def test_cached_msas_mean_no_server_flag(tmp_path, manifest):
    script = tmp_path / "run.sh"
    write_run_script(manifest, script, out_root=tmp_path / "out", uses_cached_msa=True)
    assert "--use_msa_server" not in script.read_text()


def test_without_cached_msas_the_server_is_used(tmp_path, manifest):
    script = tmp_path / "run.sh"
    write_run_script(manifest, script, out_root=tmp_path / "out", uses_cached_msa=False)
    assert "--use_msa_server" in script.read_text()


def test_script_is_executable_and_fails_fast(tmp_path, manifest):
    script = tmp_path / "run.sh"
    write_run_script(manifest, script, out_root=tmp_path / "out", uses_cached_msa=True)
    assert script.stat().st_mode & 0o111
    assert "set -euo pipefail" in script.read_text()


def test_script_skips_completed_inputs(tmp_path, manifest):
    """A killed GPU session must be restartable with the same command."""
    script = tmp_path / "run.sh"
    write_run_script(manifest, script, out_root=tmp_path / "out", uses_cached_msa=True)
    assert "already done" in script.read_text()


def test_ligand_codes_are_five_characters():
    """Boltz resolves a ccd code against its own dictionary first."""
    assert ligand_code(0) == "SB000"
    assert len(ligand_code(999)) == 5


def test_ccd_yaml_carries_the_msa_when_given_one():
    yaml = build_ccd_yaml("MKT", "SB001", msa_path="/data/msa/T1.csv")
    assert "msa: /data/msa/T1.csv" in yaml


def test_ccd_yaml_omits_msa_when_not_given():
    assert "msa:" not in build_ccd_yaml("MKT", "SB001")


def test_ccd_yaml_requests_affinity_for_the_ligand():
    yaml = build_ccd_yaml("MKT", "SB001")
    assert "ccd: SB001" in yaml
    assert "binder: B" in yaml
